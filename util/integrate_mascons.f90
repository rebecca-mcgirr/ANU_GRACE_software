  program integrate_mascons

! program to integrate the change in volume for a provided list of mascons. It will output the volume change and surface area
! for a given set of epochs. It is designed to interface with util/grab_mascons (which will create the list of mascons in a
! requested region) to enable the integration of, for example, the Amazon, Caspian Sea, Greenland etc ....
!
! P. Tregoning
! 17 January 2020

  implicit none

  character*150             :: mascon_list,mascon_file,fit_list,soln_list
  character*150             :: message

! unit numbers
  integer, parameter :: lumsc_list=10,lusoln_list=11,lu_addnorm=12

! mascon parameters
  real(kind=8),allocatable    :: msc_to_integrate(:,:)
               
! parameters to read the EWH solutions
  integer*4                 :: n_files,n_msc
  real(kind=8), allocatable :: msc_vals(:,:,:)
  character*150             :: fit_file

! counters
  integer*4  :: ifile,imsc,n_msc_to_int

! other stuff
  real(kind=8)  :: sum_vol,sum_area,epoch
  character*100 :: line
  integer*4     :: trimlen
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! decode command line
! mascon list file
  call getarg(1,mascon_list)
  if(mascon_list(1:1) == "")then
    call status_update('FATAL','UTIL','integrate_mascons',' ' &
             ,"Runstring: integrate_mascons GRN_msc.list soln_batch5z_zeroapr_iter2_min20.list",0)
  else
    open(lumsc_list,file=mascon_list(1:trimlen(mascon_list)),status='old')
  endif

! list of solution files to read
  call getarg(2,soln_list)
  open(lusoln_list,file=soln_list(1:trimlen(soln_list)),status='old')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read in the numbers of primary mascons to be integrated
  ! the number of mascons to integrate
  read(lumsc_list,*)n_msc_to_int

  ! allocate arrays
  allocate(msc_to_integrate(n_msc_to_int,2))

  do imsc = 1,n_msc_to_int

    read(lumsc_list,*)msc_to_integrate(imsc,1:2)  ! primary mascon number and the area of the primary mascon
  enddo
  close(lumsc_list)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read in all the solution EWH values. This duplicates code in compute_ts_WRMS ....
! read the number of files (which is the number of epochs)
  read(lusoln_list,*)n_files

! read the first file name
  read(lusoln_list,'(a)')fit_file

! open it
  open(lu_addnorm,file=fit_file,status='old')

! read the header through to the first line of a mascon. Subroutine reads right through to find how
! many mascons there are, then rewinds, then leaves the file at the line of the first mascon.
  call read_fitfile_to_msc(lu_addnorm,n_msc,epoch,1)
  write(message,'(a,i7,a)')"There are ",n_msc," mascons in the solution file: "
  call status_update('STATUS','UTIL','integrate_mascons',fit_file,message,0)
  
! allocate the mascon array
  allocate(msc_vals(n_msc,n_files,3))   ! we store the EWH value, uncertainty and epoch of each estimate

! set it all back to the beginning
  close(lu_addnorm)
  rewind(lusoln_list)
  read(lusoln_list,'(a)')line
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read in all the mascon EWH values from all the fit files
  call status_update('STATUS','UTIL','integrate_mascons'," ","Reading all mascon solution values",0)

  do ifile=1,n_files
    read(lusoln_list,'(a)')fit_file
    open(lu_addnorm,file=fit_file,status='old')
    call read_fitfile_to_msc(lu_addnorm,n_msc,epoch,2)
    msc_vals(:,ifile,3) = epoch
    call read_fitfile_msc_vals(lu_addnorm,n_msc,msc_vals(:,ifile,:))
    close(lu_addnorm)
  enddo
  close(lusoln_list)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ok, so we have our list of mascons to integrate and we have all the solution EWH values. Add em up
  do ifile = 1, n_files
    sum_vol  = 0.d0
    sum_area = 0.d0

    do imsc = 1,n_msc_to_int
      sum_vol = sum_vol   + msc_vals(nint(msc_to_integrate(imsc,1)),ifile,1)*msc_to_integrate(imsc,2)
      sum_area = sum_area + msc_to_integrate(imsc,2)

!print*,ifile,nint(msc_to_integrate(imsc,1)),msc_vals(nint(msc_to_integrate(imsc,1)),ifile,1) &
!             ,msc_to_integrate(imsc,2) &
!             ,msc_vals(nint(msc_to_integrate(imsc,1)),ifile,1)*msc_to_integrate(imsc,2)/361.9e12,sum_vol/361.9e12


    enddo

! output
    write(*,'(f15.6,e17.8,f10.4, f20.4,a,f20.4,a)')msc_vals(nint(msc_to_integrate(1,1)),ifile,3),sum_vol,sum_vol/sum_area &
          ,sum_vol*1.e-9,' GTonnes',-1.d0*sum_vol/361.9e12,' m GSL'

  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call status_update('STATUS','UTIL','integrate_mascons',' ',"End of integrate_mascons",0)

  end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

