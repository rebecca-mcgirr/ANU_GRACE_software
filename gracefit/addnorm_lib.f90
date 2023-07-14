!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! subroutines related to the ADDNORM program
!
! P. Tregoning
! 9 July 2014
!
!  use_tidal_msc     :  routine to make a list of tidal amplitudes to be estimated
!  read_msc_ts_WRMS  : routine to read a variety of values, then apply a mascon-specific constraint
!  read_addnorm_cmd  : read the new (211206) addnorm command file and pass info back to main program
!  read_msc_ts_WRMS  : reads in solution-based EWH values for mascons
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_addnorm_cmd(lucmd )

! subroutine to read all the required information from the addnorm command file and
! return it to the main program
!
! P. Tregoning
! 6 December 2021
!
! PT220608: added output_correlation and use_correl_in_reg


  use addnorm_mod  ! provides all the variable declarations that are to be read from the addnorm command file

  implicit none

! passed variables
  integer*4   :: LUCMD             ! unit number of command file

! local variables
  real(kind=8):: tmp_R8,tmp_tik
  integer*4   :: i,lcmd     ! what is this?
  character*200 :: message,tmpchar,line


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RM201020: initialise logicals to false
  msc_solve     = .false.
  IC_solve      = .false.
  tmsc_solve    = .false.
  output_snorm  = .false.
  output_SVD    = .false.
  output_ascii  = .false.
  output_correl = .false.
  tikhonov      = .false.
  use_correl_in_reg = .false.
  norm_netcdf   = .false.
    
! PT211206: read the command file using the gracefit command_read subroutine
  call command_storeCommand(LUCMD,'est_IC',tmp_R8,1,1,'Estimate ICs?',0)
  if(nint(tmp_R8) == 1)IC_solve = .true.

  call command_storeCommand(LUCMD,'est_msc',tmp_R8,1,1,'Estimate mascons?',0)
  if(nint(tmp_R8) == 1)msc_solve = .true.

  call command_storeCommand(LUCMD,'est_tides',tmp_R8,1,1,'Estimate tidal mscs?',0)
  if(nint(tmp_R8) == 1)tmsc_solve = .true.

  call command_storeCommand(LUCMD,'output_snorm',tmp_R8,1,1,'Output stacked normal equations?',0)
  if(nint(tmp_R8) == 1)output_snorm = .true.

  call command_storeCommand(LUCMD,'output_SVD',tmp_R8,1,1,'Output SVD info?',0)
  if(nint(tmp_R8) == 1)output_SVD = .true.

  call command_storeCommand(LUCMD,'output_ascii_norm',tmp_R8,1,1,'Output ascii normal equations?',0)
  if(nint(tmp_R8) == 1)output_ascii = .true.

  call command_storeCommand(LUCMD,'output_correl',tmp_R8,1,1,'Output correlation matrix?',0)
  if(nint(tmp_R8) == 1)output_correl = .true.

  call command_storeCommand(LUCMD,'norm_netcdf',tmp_R8,1,1,'netcdf format for normal equations input?',0)
  if(nint(tmp_R8) == 1)norm_netcdf = .true.

  !!!!!! Output information to the screen
  if(IC_solve)  then
    call status_update('STATUS','addnorm','read_addnorm_cmd',' ','Estimate ICs                    : Y',0)
  else
    call status_update('STATUS','addnorm','read_addnorm_cmd',' ','Estimate ICs                    : N',0)
  endif
  if(msc_solve) then
    call status_update('STATUS','addnorm','read_addnorm_cmd',' ','Estimate mascons                : Y',0)
  else
    call status_update('STATUS','addnorm','read_addnorm_cmd',' ','Estimate mascons                : N',0)
  endif
  if(tmsc_solve)then
    call status_update('STATUS','addnorm','read_addnorm_cmd',' ','Estimate tidal msc              : Y',0)
  else
    call status_update('STATUS','addnorm','read_addnorm_cmd',' ','Estimate tidal msc              : N',0)
  endif

  if(output_snorm)then
    call status_update('STATUS','addnorm','read_addnorm_cmd',' ','Output stacked normal equations : Y',0)
  else
    call status_update('STATUS','addnorm','read_addnorm_cmd',' ','Output stacked normal equations : N',0)
  endif
  if(output_SVD)then
    call status_update('STATUS','addnorm','read_addnorm_cmd',' ','Output SVD information          : Y',0)
  else
    call status_update('STATUS','addnorm','read_addnorm_cmd',' ','Output SVD information          : N',0)
  endif
  if(output_ascii)then
    call status_update('STATUS','addnorm','read_addnorm_cmd',' ','Output ascii normal equations   : Y',0)
  else
    call status_update('STATUS','addnorm','read_addnorm_cmd',' ','Output ascii normal equations   : N',0)
  endif
  if(output_correl)then
    call status_update('STATUS','addnorm','read_addnorm_cmd',' ','Output correlation matrix       : Y',0)
  else
    call status_update('STATUS','addnorm','read_addnorm_cmd',' ','Output correlation matrix       : N',0)
  endif
  if(norm_netcdf)then
    call status_update('STATUS','addnorm','read_addnorm_cmd',' ','NetCDF format normal equations  : Y',0)
  else
    call status_update('STATUS','addnorm','read_addnorm_cmd',' ','NetCDF format normal equations  : N',0)
  endif


  !!!!!!!!!!!!!!!!!!!!!!!!
  !!   regularisation  !!!
  !!!!!!!!!!!!!!!!!!!!!!!!
  !name of the regularisation file
  tmpchar = " "
  call get_keyword(LUCMD,'reg_file',tmpchar,lcmd,1)
  lcmd = 1
  call first_word(tmpchar,msc_constraint_file,lcmd)

  ! scale factor for regularisation file
  call command_storeCommand(LUCMD,'scale_reg_file',msc_scl_factor,1,1,'Mascon regularisation scale factor',0)

  ! tikhonov regularisation
  call command_storeCommand(LUCMD,'tikhonov_reg',tmp_tik,1,2,'Tikhonov regularisation',0)
  if(nint(tmp_tik) == 1)then
    tikhonov = .true.
    msc_constraint_file = 'tikhonov'
    call command_storeCommand(LUCMD,'tikhonov_reg',lambda,2,2,'Tikhonov regularisation',0)
    write(message,'(a,f15.4)')'Tikhonov regularisation with lambda = ',lambda
    call status_update('STATUS','addnorm','read_addnorm_cmd',' ',message,0)
  else    
    write(message,'(a,f7.2)')'Mascon regularisation file with scale factor = ',msc_scl_factor
    call status_update('STATUS','addnorm','read_addnorm_cmd',msc_constraint_file,message,0)
  endif

  ! use adaptive regularisation in building regularisation matrix?
  call command_storeCommand(LUCMD,'use_adaptive_reg',tmp_R8,1,1,'Use adaptive reg. to build regularisation matrix?',0)
  if(nint(tmp_R8) == 1)then
    use_adaptive_reg = .true.
    call command_storeCommand(LUCMD,'use_adaptive_reg',adapt_sigma,2,2,'Adaptive regularisation sigma',0)
  endif

  if(use_adaptive_reg)then
    call status_update('STATUS','addnorm','read_addnorm_cmd',' ','Use adaptive reg to build reg   : Y',0)
  else
    call status_update('STATUS','addnorm','read_addnorm_cmd',' ','Use adaptive reg to build reg   : N',0)
  endif
   
  ! use mascon correlations in building regularisation matrix?
  call command_storeCommand(LUCMD,'use_correl_in_reg',tmp_R8,1,1,'Use mascon correlations to build regularisation matrix?',0)
  if(nint(tmp_R8) == 1)use_correl_in_reg = .true.
  if(use_correl_in_reg)then
    call status_update('STATUS','addnorm','read_addnorm_cmd',' ','Use mascon correl to build reg  : Y',0)
  else
    call status_update('STATUS','addnorm','read_addnorm_cmd',' ','Use mascon correl to build reg  : N',0)
  endif
   
  !!!!!!!!!!!!!!!!!!!!!!!!
  !!  msc constraints  !!!
  !!!!!!!!!!!!!!!!!!!!!!!!
  call command_storeCommand(lucmd,'apr_mascons',msc_param_const,1,1,'Mascon Constraint',0)
  write(message,'(a,f17.3,a)')'Mascon parameter constraint = ',msc_param_const,' (metres)'
  call status_update('STATUS','addnorm','read_addnorm_cmd',' ',message,0)


  !!!!!!!!!!!!!!!!!!!!!!!!
  !!  conserv. of mass !!!
  !!!!!!!!!!!!!!!!!!!!!!!!
  ! conservation of mass uncertainty
  call command_storeCommand(LUCMD,'cons_mass_sigma',cons_mass_sigma,1,1,'Conservation of mass sigma',0)
  ! conservation of mass offset (for cases where the a priori model does not conserve mass
  call command_storeCommand(LUCMD,'cons_mass_offset',mass_offset,1,1,'conservation of mass offset',0)
  write(message,'(a,2f17.2)')'Conservation of mass sigma and offset (kg) = ',cons_mass_sigma,mass_offset
  call status_update('STATUS','addnorm','read_addnorm_cmd',' ',message,0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the normal equation file information
  call command_storeCommand(LUCMD,'num_files',tmp_R8,1,1,'number of normal equation files',0)
  nfiles = nint(tmp_R8)
  allocate(normeq_files(nfiles))

  ! now, we need to find the line "num_files" in the command file, then read the normal equations files
  rewind(lucmd)
  message = " "
  do while (message(1:9) /= "num_files")
    read(lucmd,'(a)')message
    do i=1,20
      if(message(i:i+8) == "num_files")message(1:9) = "num_files"
    enddo
  enddo

  ! ok, so we are now at the line of the first normal equations filename
  ! PT220124: make it ignore files commented out with a "#"
 i = 0
  do while (i < nfiles)
     read(lucmd,'(a)')line
     if(line(1:1) /= "#")then
       i = i+1
       read(line,'(a)')normeq_files(i)
     endif
  enddo
  write(message,'(a,i5,a)')'Have ',nfiles,' normal equation files to be read'
  call status_update('STATUS','addnorm','read_addnorm_cmd',' ',message,0)
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   



  return
  end subroutine read_addnorm_cmd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_msc_ts_WRMS(luWRMS,msc_WRMS_file,constraint_model,soln_epoch,nparam_t,nmascon_t,msc_const)

! read in a file of mascon time series metrics, then apply mascon-specific constraints
! 
! P. Tregoning/J. Pfeffer
! 22 January 2020
!
! MODS:
! PT200225: changed to read cosine/sine amplitudes separately rather than the annual amplitude (allows to reconstruct the signal)
! PT200225: modified to read different header information
! PT220608: changed dimensioning of msc_const to be for all parameters, not just mascons

  implicit none

! passed variables
  integer*4     :: luWRMS                            ! unit number for mascon WRMS file
  character*150 :: msc_WRMS_file                     ! file containing metrics of mascon time series
  integer*4     :: nparam_t                          ! dimension of msc_const
  integer*4     :: nmascon_t                         ! number of mascons expected in file
  real(kind=8)  :: msc_const(nparam_t,nparam_t)      ! mascon constraint matrix
  character*20  :: constraint_model                  ! which combination of metrics to use to assign the constraint
  real(kind=8)  :: soln_epoch                        ! decimal year of the solution epoch

! variables to read in mascon metric information
  real(kind=8),allocatable :: msc_crds(:,:)          ! mascon longitude/latitude
  real(kind=8),allocatable :: msc_annual(:,:)          ! cosine/sine amplitudes of annual signal for each mascon
  real(kind=8),allocatable :: msc_wrms(:,:)          ! the wrms values (**list them)
  real(kind=8),allocatable :: msc_monthly(:,:)       ! monthly low-frequency signal values from mascon time series
  real(kind=8),allocatable :: epochs(:)

  integer*4                :: nmonths                ! number of monthly values of low-freq signal
  integer*4                :: nmsc                   ! number of mascons in the WRMS file
  integer*4                :: nrms                   ! number of RMS values in file (between annual ampl and low-freq values)
  integer*4                :: ntot                   ! total number of columns
  character*100            :: line

! counters
  integer*4                :: imsc,imsc2,iRMS,imonth,ioerr
  character*100            :: message

! variables to interpolate the low-frequency msc metrics to the solution epoch
  real(kind=8)             :: dt
  integer*4                :: i1,i2
  real(kind=8),allocatable :: lowfrq(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!
! read header info. First line is a comment, giving column headers. Second line contains number of mascons.
  open(luWRMS,file=msc_WRMS_file,status='old',iostat=ioerr)
  if(ioerr /= 0)call status_update('FATAL','ADDNORM','addnorm',msc_WRMS_file,"Error opening file",0)

  read(luWRMS,'(a)')line
  read(luWRMS,*)nmsc,nmonths,ntot
  
! the "ntot" value in the file includes lat/lon/density/cosine_ampl/sine_ampl
  nrms = ntot - 5
  if(nmsc /= nmascon_t)then
    write(message,'(a,i7,a,i7,a)')"Mismatch in number of mascons expected (",nmascon_t,") and found (",nmsc,")"
    call status_update('FATAL','UTIL','read_msc_ts_WRMS',msc_WRMS_file,message,0)
  else
    allocate(msc_crds(nmsc,3))          ! lat, lon, density
    allocate(msc_annual(nmsc,2))        ! cosine/sine amplitudes
    allocate(msc_WRMS(nmsc,nrms))
    allocate(msc_monthly(nmsc,nmonths))
    allocate(epochs(nmonths))
    allocate(lowfrq(nmsc))
  endif

! read the epochs contained in the file. They are in decimal years
  read(luWRMS,*)(epochs(imonth),imonth=1,nmonths)

! ok, now loop through mascons and read in all the information
  do imsc=1,nmsc
    read(luWRMS,*)msc_crds(imsc,:),msc_annual(imsc,:),(msc_wrms(imsc,iRMS),iRMS=1,nrms) &
                ,(msc_monthly(imsc,imonth),imonth=1,nmonths)
  enddo

! now, assign a mascon constraint to each mascon, depending on what model was requested
  do imsc=1,nmsc
    if(constraint_model(1:4) == "WRMS")then
      ! we just want to use the total WRMS of the time series. This is the first RMS value, fifth column in the file.
      if(imsc == 1)call status_update('STATUS','ADDNORM','read_msc_ts_WRMS',msc_WRMS_file &
                   ,'Applying WRMS of time series as msc constraint',0)
      msc_const(imsc,imsc) = msc_const(imsc,imsc)+1.d0/(msc_wrms(imsc,1))**2

    else if (constraint_model(1:8) == "landWRMS")then
      ! we just want to use the total WRMS of the time series but only for continental mascons. set ocean mascons to 0.1 m
      if(imsc == 1)call status_update('STATUS','ADDNORM','read_msc_ts_WRMS',msc_WRMS_file &
                   ,'Applying WRMS of time series as msc constraint to land, 0.1m to ocean',0)
      if(msc_crds(imsc,3) < 1010.d0)then
        msc_const(imsc,imsc) = msc_const(imsc,imsc)+1.d0/(msc_wrms(imsc,1))**2
      else
        msc_const(imsc,imsc) = msc_const(imsc,imsc)+1.d0/(0.1d0)**2    ! set a 0.1 m constraint on oceans if we use land-only WRMS
      endif

    else if (constraint_model(1:6) == "LOWFRQ")then
      if(imsc == 1)then
        call status_update('STATUS','ADDNORM','read_msc_ts_WRMS',msc_WRMS_file &
                   ,'Applying low-frequency signal of time series as msc constraint',0)

        ! find the two monthly values that straddle the required epoch
        imonth = 1
        do while (epochs(imonth) < soln_epoch .and. imonth <= nmonths)
          if(epochs(imonth) < soln_epoch)imonth = imonth+1
        enddo
        ! interpolate between i2 and i1
        if(imonth == 1)then
          write(message,'(a,2(f15.4,a))')"Requested epoch (",soln_epoch,") before first epoch in file (",epochs(1),")"
          call status_update('STATUS','ADDNORM','read_msc_ts_WRMS',constraint_model,message,0)
          lowfrq(:) = msc_monthly(:,1)
        else if (imonth > nmonths)then
          write(message,'(a,2(f15.4,a))')"Requested epoch (",soln_epoch,") after last epoch in file (",epochs(nmonths),")"
          call status_update('STATUS','ADDNORM','read_msc_ts_WRMS',constraint_model,message,0)
          lowfrq(:) = msc_monthly(:,nmonths)
        else
          i2 = imonth
          i1 = imonth - 1
          dt = epochs(i2)-epochs(i1)
          write(message,'(3(a,f15.4))')"Requested epoch (",soln_epoch,") between epochs",epochs(i1)," and",epochs(i2)
          call status_update('STATUS','ADDNORM','read_msc_ts_WRMS',constraint_model,message,0)
          do imsc2=1,nmsc 
            lowfrq(imsc2) = msc_monthly(imsc2,i1) + (soln_epoch-epochs(i1)) * (msc_monthly(imsc2,i2)-msc_monthly(imsc2,i1))/dt
          enddo
        endif
      endif

      ! just the low-frequency part
      if(constraint_model(1:11) == "LOWFRQ_ONLY")then
        msc_const(imsc,imsc) = msc_const(imsc,imsc)+1.d0/lowfrq(imsc)**2

      elseif(constraint_model(1:13) == "LOWFRQ+ANNUAL")then
        msc_const(imsc,imsc) = msc_const(imsc,imsc)+1.d0/( dsqrt(msc_annual(imsc,1)**2+msc_annual(imsc,2)**2) &
                                                           +2.d0*dabs(lowfrq(imsc)) )**2

      elseif(constraint_model(1:13) == "LOWFRQ+ANN+HF")then
        msc_const(imsc,imsc) = msc_const(imsc,imsc) &
                              +1.d0/( dsqrt(msc_annual(imsc,1)**2+msc_annual(imsc,2)**2) &
                              +dabs(msc_wrms(imsc,4))+4.d0*dabs(lowfrq(imsc)) )**2  ! annual_ampl+hf_rms+lowfreq_signal

      elseif(constraint_model(1:17) == "LOWFRQ_land_ANNHF")then
        if(msc_crds(imsc,3) < 1010.d0)then
          msc_const(imsc,imsc) = msc_const(imsc,imsc) &
                              +1.d0/( dsqrt(msc_annual(imsc,1)**2+msc_annual(imsc,2)**2)+dabs(msc_wrms(imsc,4)) &
                                    +4.d0*dabs(lowfrq(imsc)) )**2  ! annual_ampl+hf_rms+lowfreq_signal
        else
          msc_const(imsc,imsc) = 1.d0/0.1d0**2
        endif
      endif

    endif

  enddo

  end subroutine read_msc_ts_WRMS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine use_tidal_msc(nmascon_t,max_msc_tide,nmsc_tidal_t,mcon_tides,tidal_lookup,prm_tidal,tide_list_epoch)

  ! makes a complete list of all tidal amplitudes per mascon to be estimated. The list
  ! gets added to each time that a new set of normal equations are read in if there are
  ! new entries in the normal equations that are not in the list.
  !
  ! P. Tregoning
  ! 9 July 2014

  implicit none

  integer*4   , intent(in)     ::  nmascon_t                                    ! number of mascons
  integer*4   , intent(in)     ::  max_msc_tide                                 ! max number of tidal constituents that can be estimated at each mascon
  integer*4   , intent(inout)  ::  nmsc_tidal_t                                 ! total number of tidal amplitudes to date that will be estimated
  ! tidal amplitudes estimated in this set of normal equations
  integer*4   , intent(in)     ::  mcon_tides(nmascon_t)
  ! lookup table 
  integer*4   , intent(inout)  ::  tidal_lookup(nmascon_t,max_msc_tide,2)       ! lookup table for which tidal amplitudes are estimable (and how many times)
  character*30, intent(inout)  ::  prm_tidal(nmascon_t*max_msc_tide*2)          ! parameter label for the tidal amplitudes
  integer*4   ,intent(out)     ::  tide_list_epoch(nmascon_t*max_msc_tide*2,2)  ! an actual list of the mascons and tidal amplitudes estimated for this particular set of normal equations     
  integer*4, allocatable       ::  tidal_pointer(:,:)

  ! local variables
  integer*4      ::  imsc,i,j,k
  logical        ::  bitmap
  character*150  ::  message
  character*2    ::  tide_label(5)
  integer*4      ::  tide_count

  ! define the tidal constituents
  tide_label(1) = "M2"
  tide_label(2) = "O1"
  tide_label(3) = "S2"
  tide_label(4) = "K1"
  tide_label(5) = "K2"

  allocate(tidal_pointer(nmascon_t*max_msc_tide,2))

  ! initialise variables the first time we enter the subroutine
  if (nmsc_tidal_t == 0 ) then
     tidal_lookup = 0
     tide_count = 0
  endif

  ! loop through the bit-mapped mascon values and store as appropriate
  do imsc=1,nmascon_t
     if(mcon_tides(imsc) > 0) then   ! we have some estimated tidal amplitudes for this mascon
        do j=1,max_msc_tide
           if(bitmap(mcon_tides(imsc),j)) then
              tidal_lookup(imsc,j,1) = tidal_lookup(imsc,j,1) + 1
              tide_count = tide_count + 1
              tide_list_epoch(tide_count,1) = imsc
              tide_list_epoch(tide_count,2) = j
           endif
        enddo
     endif
  enddo

  ! now loop through and find how many non-zero elements in tidal_lookup. This tells us how many tidal amplitudes are to be estimated thus far.
  nmsc_tidal_t = 0
  do imsc=1,nmascon_t
     do j=1,max_msc_tide
        if(tidal_lookup(imsc,j,1) > 0)then
           nmsc_tidal_t = nmsc_tidal_t + 1
           tidal_pointer(nmsc_tidal_t,1) = imsc   ! the mascon number for this tidal amplitude parameter
           tidal_pointer(nmsc_tidal_t,2) = j      ! the tidal constituent for this tidal amplitude parameter
           ! PT160216: changed alignment of print statement to match format of gracefit files (chagned i6 to i5, added space after (m)  )
           write(prm_tidal((nmsc_tidal_t-1)*2+1),'(i5,a,i4,1x,a2,a)')nmascon_t+(nmsc_tidal_t-1)*2+1,". TMC" &
                ,imsc,tide_label(j),"sin (m) "
           write(prm_tidal((nmsc_tidal_t-1)*2+2),'(i5,a,i4,1x,a2,a)')nmascon_t+(nmsc_tidal_t-1)*2+2,". TMC"  &
                ,imsc,tide_label(j),"cos (m) "
        endif
     enddo
  enddo


  ! now transfer the row number for each constituent back into the second leaf of the tidal_lookup array
  tidal_lookup(:,:,2) = 0
  do i=1,nmsc_tidal_t
     tidal_lookup(tidal_pointer(i,1),tidal_pointer(i,2),2) = i   ! this assigns the row number in the total of the summed tidal amplitudes
  enddo

  !DEBUG
  !  print*,'DEBUG  rows for tidal amplitudes for MC4535 = ',tidal_lookup(4535,:,2),' #',tidal_lookup(4535,:,1)
  !  print*,'DEBUG  rows for tidal amplitudes for MC4536 = ',tidal_lookup(4536,:,2),' #',tidal_lookup(4536,:,1)

  !DEBUG
  !  do i=1,nmsc_tidal_t*2
  !    print*, "DEBUG  TMC label is:  ",prm_tidal(i)
  !  enddo
  return
end subroutine use_tidal_msc





