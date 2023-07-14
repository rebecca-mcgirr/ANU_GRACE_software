  subroutine get_sca_data(interp_flag,calling_prog,which_file,SCA_file)

! subroutine to read the L1B sca data (using existing subroutines), then infill any data gaps. It returns 
! (via mod_subs/accelerom_mod sca_obs array) an array of star camera observations with no gaps
!
! P. Tregoning
! 4 November 2019

  use accred_mod  ! defines allocatable arrays for star camera obs (sca_obs) and the number of obs (STARn)


  implicit none

  logical      , intent(in)  :: interp_flag     ! T: interpolate values for SCA1B gaps, F: don't interpolate
  integer*4    , intent(in)  :: which_file      ! 1: GRACE A/C; 2: GRACE B/D. Determines into which array the vals are stored
  character(*) , intent(in ) :: SCA_file        ! name of SCA file
  character(*) , intent(in)  :: calling_prog    ! name of program calling this subroutine

! local variables
  integer*4,parameter        :: lu_sca=60       ! unit number for SCA file
  integer*4                  :: n_sca           ! number of SCA obs as listed in the SCA1B header
  integer*4                  :: n_sca_full      ! number of obs in SCA data span (including epochs with no data)
  integer*4                  :: mission         ! determined when reading SCA1B header file

! variables for temporary read of SCA data
  real(kind=8),allocatable   :: sca_obs_tmp(:,:),sca_obs_tmp2(:,:)

! counters etc
  integer*4     :: irow,iepoch,i,start_gap,end_gap,npoints
  character*150 :: message
  
! open the star camera file 
  call level1B_open(lu_sca,calling_prog,SCA_file)

! read the header and allocate the temporary array
  mission = -1
  call sca_read_hdr(lu_sca,calling_prog,SCA_file,mission,n_sca,sca_step )
  allocate(sca_obs_tmp(n_sca,6))

! sca_obs contains:
! gracesec, 4xquaternions  star_camera_flag
  call sca_read_data(lu_sca,calling_prog,SCA_file,mission,n_sca, sca_obs_tmp)
  sca_step = sca_obs_tmp(2,1)-sca_obs_tmp(1,1)
  do iepoch=3,n_sca
     sca_step = min(sca_step, sca_obs_tmp(iepoch,1)-sca_obs_tmp(iepoch-1,1))
  enddo


! update the number of star camera observations to now include missing data
  n_sca_full = nint( (sca_obs_tmp(n_sca,1)-sca_obs_tmp(1,1))/sca_step) + 1

! PT180620: now, fill out the SCA data, leaving zero entries where there are data gaps (to be filled in by interpolation later)
  allocate(sca_obs_tmp2(n_sca_full,6))
  sca_obs_tmp2 = 0.d0
  sca_obs_tmp2(1,:) = sca_obs_tmp(1,:)
  do i=2,n_sca
    irow = nint( (sca_obs_tmp(i,1)-sca_obs_tmp(1,1))/sca_step) + 1
    sca_obs_tmp2( irow,:) = sca_obs_tmp(i,:)
!print*,'get_sca_data i, irow',i,irow,sca_obs_tmp2( irow,1:3)
  enddo

  write(message,'(a,i10)')" Number of star camera observations (including gap epochs): ",n_sca_full
  call status_update('STATUS',calling_prog,'get_sca_data',' ',message,0)

  if(interp_flag)then
! now, fill in any missing values using a 2nd-order quadratic interpolation
    n_sca = n_sca_full
    do i=1,n_sca-1
      if(nint(sca_obs_tmp2(i,1)) == 0)then      ! it is a missing epoch
        if(i == 1)then
          sca_obs_tmp2(i,:) = sca_obs_tmp2(i+1,:)    ! just make it all the same as epoch 2
        else
          ! using 800 points either side of the gap, fit a quadratic to the data and use modelled values to infill missing obs
          ! first, find the end of the data gap
          start_gap = i
          end_gap = i
          do while (end_gap < n_sca .and. nint(sca_obs_tmp2(end_gap,1)) == 0)
            end_gap = end_gap + 1
          enddo
          ! now adjust by 800 epochs before/after
          npoints = 100
          if (start_gap > npoints)then
            start_gap = start_gap - npoints
          else
            start_gap = 1
          endif
          if(int((n_sca - end_gap)/sca_step) > npoints)then
            end_gap = end_gap + npoints
          else
            end_gap = n_sca
          endif
          ! fit the quadratic model to the span of data
          call infill_TS(.false.,calling_prog,"quadratic ","SCA",end_gap-start_gap+1,5 &
                         ,sca_obs_tmp2(start_gap:end_gap,:),sca_step)
        endif
      endif
    enddo

    call status_update('STATUS',calling_prog,'get_sca_data',' ',"Have infilled missing star camera epochs",0)
  endif

! copy to the appropriate sca array defined in mod_subs/accred_mod.f90
  if(which_file == 1)then
    allocate(sca_obs(n_sca_full,6))
    sca_obs = sca_obs_tmp2
  else if (which_file == 2) then
    allocate(sca_obs_2(n_sca_full,6))
    sca_obs_2 = sca_obs_tmp2
  endif

! deallocate arrays
  deallocate(sca_obs_tmp)
  deallocate(sca_obs_tmp2)

! close the SCA file
  close(lu_sca)

  return

  end subroutine get_sca_data