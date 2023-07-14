  program apriori_mascons

! program to generate a file (in vcv format) of mascon values. The input information will be
! a set of parameters for a quadratic+annual model which is a best-fit to a time series of
! mass change estimates (in terms of EWH) for each mascon in the V002 mascon file.
!
! P. Tregoning
! 25 September 2019
!
! MODS
! PT200225: modified to read a second input file type: from the output of get_rms_and_lowpass_values.py. This file contains
!           information on the annual amplitude signal plus a low-frequency component from a prior assessment of time series
!           for each mascon
! PT200830: output zero values for ocean mascons
! PT211019: use a threshold value - from command line - to determine whether to output mascon updated values or zero.

  use mascon_mod          ! defines all the mascon arrays

  implicit none

! command line arguments
  character*150 :: mascon_file      ! file containing mascon geometry
  character*150 :: msc_model_file   ! file containing the quadratic+annual model parameters for the corresponding mascons
  character*150 :: user_mascon_file ! a possible different mascon file onto which we want to calculate the EWH
  integer*4     :: date(3),doy      ! epoch onto which to calculate the mascon EWH values (yr,mo,dd)
  character*150 :: usr_msc_file     ! possible different mascon file onto which to calculate the EWH values
  character*150 :: output_file      ! hardwired file of form "msc_apriori_YYYY_MM.vcv"
  integer*4     :: model            ! type of model required 
  character*4   :: zero_ocean       ! RM221101: flag to zero ocean mascons
  real(kind=8)  :: msc_threshold    ! dabs(mascon values) must exceed this value to be output, otherwise zero.

! mascon EWH model variables
  real(kind=8),allocatable :: msc_model(:,:)       ! array to hold the mascon EWH model (total_prim, 5)
  real(kind=8),allocatable :: apriori_msc(:)       ! array to store the calculated a priori EWH mascon values
  real(kind=8),allocatable :: msc_model_apriori(:) ! array to hold the a priori values as computed from the low-freq file
  logical                  :: quad_model           ! flag of what type of input model file was provided

! mascon file and user mascon geometry variables
  integer*4                :: tern_number
  real(kind=8),allocatable :: mcon_prim_usr(:,:)
  integer*4                :: total_prim_usr
  logical                  :: same_msc_geom
  real(kind=8)             :: lat_spacing
  integer*4                :: model_msc_number

! variables for ensuring conservation of mass
  real(kind=8)  :: land_area,land_vol,ocean_area,ocean_vol
  real(kind=8)  :: vol_error,ocean_vol_orig,vol_error_orig,ocean_vol_corr,mass_error,mass_error_orig
  real(kind=8),allocatable :: mass_corr(:)

! variables related to using a file containing information on annual amplitude and low-frequency signals for the a priori model
  character*10  :: file_type        ! indicator whether it is a simple quadratic model file or one of RMS values and low-freq model
  integer*4     :: nmsc_model       ! the number of mascons in the low-freq model file
  integer*4     :: n_LF_epochs      ! number of epochs in the low-frequency model
  integer*4     :: n_columns        ! number of columns per line in the low-frequency file
!  real(kind=8),allocatable :: LF_epochs(:)       ! low-frequency epochs (in decimal years)
!  real(kind=8),allocatable :: LF_coords(:,:)     ! lat/lon/density of mascons in the low-frequency file (not needed in this program?)
!  real(kind=8),allocatable :: LF_annual(:,2)     ! cosine/sine amplitudes for each mascon
!  real(

! local variables
  real(kind=8)  :: dec_yr           ! decimal year of required epoch
  integer*4     :: ioerr
  integer*4     :: imsc
  real(kind=8)  :: pi,w,dt,lon,lat
  character*110 :: arg,message
  integer*4     :: trimlen
  character*100 :: line

! unit numbers
  integer*4, parameter :: lumsc=10, lumsc_model=11,lumsc_usr=12,luout=13
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialise variables and declare some constants
  date = 0
  pi = 4.d0*datan(1.d0)
  w = 2.d0*pi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the command line
  call getarg(1,mascon_file)
  if(mascon_file(1:1) == " ")then
    print*,"apriori_mascons mascons_stage4_V002 msc_model.all 2016 01 15 mascons_stage4_V002b bit-mapped_model zero_ocean [msc_threshold (in m)]"
    print*," where bit-mapped models are: "
    print*,"       1: low-frequency component"
    print*,"       2: annual amplitude"
    print*,"       4: WRMS of total time series"
    print*,"       8: WRMS of land time series, 0.1m ocean"
    print*,"      16: WRMS of high-frequecy component"
    print*,"      32: annual amplitude if > 0.2 m"
    print*," "
    print*," default mascon threshold value is 0.0 (i.e. every mascon gets an updated value)"
    stop
  endif

  ! file containing msc EWH WRMS and/or a priori model
  call getarg(2,msc_model_file)  ! mascon file of modelled values
  open(lumsc_model,file=msc_model_file,status='old',iostat=ioerr)
  if(ioerr /= 0)then
    call status_update('FATAL','UTIL','apriori_mascons',msc_model_file &
                     ,"Error opening mascon EWH/apriori model file. Does it exist?",ioerr)
  else
    call status_update('STATUS','UTIL','apriori_mascons',msc_model_file,"Opened mascon EWH/apriori model file. ",ioerr)
  endif

  ! required epoch
  call getarg(3,arg)           ! year
  read(arg,*)date(1)
  call getarg(4,arg)           ! month
  read(arg,*)date(2)
  call getarg(5,arg)           ! day
  read(arg,*)date(3)

  call getarg(6,usr_msc_file)  ! mascon file of geometry onto which we calculate the values
  if(usr_msc_file(1:trimlen(usr_msc_file)) == mascon_file(1:trimlen(mascon_file)))then
    call status_update('STATUS','UTIL','apriori_mascons',mascon_file,"Will calculate mascons on same coordinates in file.",0)
    usr_msc_file = mascon_file
    same_msc_geom = .true.
  else
    call status_update('STATUS','UTIL','apriori_mascons',usr_msc_file,"Will calculate mascons on coordinates in different file.",0)
    same_msc_geom = .false.
  endif

  ! bit-mapped model required
  call getarg(7,arg)
  read(arg,*)model

  !RM221101: zero ocean flag
  call getarg(8,zero_ocean)
  if(zero_ocean == "Y")then
    call status_update('STATUS','UTIL','apriori_mascons',msc_model_file,"Will set all ocean mascons to zero.",0)
  endif
    

  ! mascon threshold value for updating a priori values
  arg = " "
  msc_threshold = 0.d0
  call getarg(9,arg)
  if(arg(1:1) == " ")then
    call status_update('STATUS','UTIL','apriori_mascons',' ',"No mascon threshold input. All mascon a priori values will update",0)
  else
    read(arg,*)msc_threshold
    write(message,'(a,f10.4,a)')"Only mascons with abs(value) > ",msc_threshold &
                               ," will be updated. All other mascons will have zero value"
    call status_update('STATUS','UTIL','apriori_mascons',' ',message,0)
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! check whether the mascon geometry file for the output is the same as that of the msc_model file. If not, read the output mascon geometry
! and store it somewhere.
  if(.not. same_msc_geom)then
    call read_msc_hdr(lumsc_usr,usr_msc_file,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                                ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
    call allocate_mascon_arrays
    ! PT180220: using trimlen here doesn't work - the subroutine is expecting a C*150, so pass the whole thing through
    !  call read_mascon_file(lumsc_in,mascon_file(1:trimlen(mascon_file)))
    call read_mascon_file(lumsc_usr,mascon_file)

    ! now allocate a separate set of mascon arrays in which to store these values. We just need the primary mascon coords (lat/lon)
    allocate(mcon_prim_usr(max_prim,4))    ! primary mascon info (lat/lon/rad/area/hgt/density/#sec/#tern/#first secondary/tidal flag/%land/geoid)
    mcon_prim_usr(:,1) = mcon_prim(:,1)*180.d0/pi    ! save the latitudes
    mcon_prim_usr(:,2) = mcon_prim(:,2)*180.d0/pi    ! save the longitudes
    mcon_prim_usr(:,3) = mcon_prim(:,4)    ! save the area of the primary mascon
    mcon_prim_usr(:,4) = mcon_prim(:,6)    ! save the density of the primary mascon
    total_prim_usr = total_prim
    write(message,'(a,i8,a)')"Mascons will be computed using mascon file with ",total_prim_usr," mascons"
    call status_update('STATUS','UTIL','apriori_mascons',usr_msc_file,message,0)

    ! deallocate the arrays
    call deallocate_mascon_arrays
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the mascon file that relates to the mascon EWH model
  call read_msc_hdr(lumsc,mascon_file,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                                ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
  call allocate_mascon_arrays
! PT180220: using trimlen here doesn't work - the subroutine is expecting a C*150, so pass the whole thing through
!  call read_mascon_file(lumsc_in,mascon_file(1:trimlen(mascon_file)))
  call read_mascon_file(lumsc,mascon_file)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set up the mcon_prim_usr array in the event that the same mascon geometry was used for the EWH model and the output. We
! can then use the array "mcon_prim_usr" for either case.
  if(same_msc_geom)then
    total_prim_usr = total_prim
    allocate(mcon_prim_usr(total_prim_usr,4))
    mcon_prim_usr(:,1) = mcon_prim(:,1)    ! save the latitudes
    mcon_prim_usr(:,2) = mcon_prim(:,2)    ! save the longitudes
    mcon_prim_usr(:,3) = mcon_prim(:,4)    ! save the area of the primary mascon
    mcon_prim_usr(:,4) = mcon_prim(:,6)    ! save the density of the primary mascon
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! determine the solution epoch
  call ymd_to_doy(date,doy)
  dec_yr = dble(date(1)) + dble(doy)/365.d0
  dt = dec_yr - 2003.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! first, determine from the first line whether it is a quadratic model file or the output from get_rms_and_lowpass_values.py
  read(lumsc_model,'(a16)')line
  backspace(lumsc_model)
  if(line(1:16) == "# n_msc n_epochs")then
    allocate(msc_model_apriori(total_prim))
    call status_update('STATUS','UTIL','apriori_mascons',msc_model_file,"reading lowpass and RMS file",0)
    close(lumsc_model)
    call read_msc_model('APRIORI_MASCONS',lumsc_model,msc_model_file,'APR',     model,dec_yr,total_prim,msc_model_apriori)
    quad_model = .false.
    msc_model_apriori = msc_model_apriori*1.d3  ! convert to mm to be compatible with the quad_model files

  else
    ! read in the mascon EWH parameter model. File has format:
    ! longitude  latitude  offset  rate  ampl_sine  ampl_cos  accel cubic
    allocate(msc_model(total_prim,6))
    call status_update('STATUS','UTIL','apriori_mascons',msc_model_file,"reading basic linear/quadratic model file",0)
    do imsc=1,total_prim
      read(lumsc_model,*)lon,lat,msc_model(imsc,:)
    enddo

    quad_model = .true.   ! set this to invoke computation of the model lower down in the program
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write a VCV file header
  write(output_file,'(a,i4,a,i2.2,a)')"msc_apriori_",date(1),"_",date(2),".vcv"
  open(luout,file=output_file,status='unknown')
  write(luout,'(a,i4,i3.2,a,a)')"V2 APRIORI  VCV  created for ",date(1:2)," from file: ",msc_model_file(1:trimlen(msc_model_file))
! PT200117: add the number of mascons into the second line of the file
  write(luout,'(i6,a)')total_prim_usr," mascons in file"
  write(luout,'(a,a)')"Output primary mascon geometry file: ",usr_msc_file(1:trimlen(usr_msc_file))
  write(luout,'(a)')"An incomplete VCV file, containing only mascon values to be read"

  write(luout,'(a)')"SOLUTION A PRIORI AND VECTOR:"
  write(luout,'(a)')"PARAMETER                     A PRIORI             VECTOR            SIGMA"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  now assign the apriori model on each primary mascon - either using the same or user-entered mascon geometry
  ! if output mascon geometry is the same as msc_model_file geometry
  if(same_msc_geom)then
    allocate(apriori_msc(total_prim))
    do imsc=1,total_prim
      if(quad_model)then
         apriori_msc(imsc) =  msc_model(imsc,1) + msc_model(imsc,2)*dt + msc_model(imsc,3)*dsin(w*dt)  &
                         + msc_model(imsc,4)*dcos(w*dt) + msc_model(imsc,5)*dt**2 + msc_model(imsc,6)*dt**3
      else
        if(mcon_prim_usr(imsc,4) > 1000.d0 .and. zero_ocean == "Y")then
          apriori_msc(imsc) = 0.d0
        else
          apriori_msc(imsc) = msc_model_apriori(imsc)
        endif
      endif
    enddo

  ! if user has specified a different output mascon geometry
  else
    allocate(apriori_msc(total_prim_usr))
    ! set the ternary latitudinal spacing
    lat_spacing = 10.d0/60.d0     ! default is 10' ternary masons in latitude

    ! establish the array that contains how many ternary latitude bands and how many ternary mascons per band
    call tern_lat_bands_ell(lat_spacing)

    do imsc=1,total_prim_usr

      ! need to find out which primary mascon in the msc_model_file geometry this user primary mascon resides in

      ! find which ternary it is in
      call calc_which_ternary(.false.,mcon_prim_usr(imsc,1),mcon_prim_usr(imsc,2),lat_spacing,tern_number)
      model_msc_number = mcon_tern_ptr(tern_number,1)

      ! now calculate the apriori EWH using the appropriate primary mascon model values
      if(quad_model)then
        apriori_msc(imsc) =  msc_model(model_msc_number,1) + msc_model(model_msc_number,2)*dt &                     ! offset and rate
                         + msc_model(model_msc_number,3)*dsin(w*dt) + msc_model(model_msc_number,4)*dcos(w*dt) &    ! annual variation
                         + msc_model(model_msc_number,5)*dt**2     &                                                ! quadratic term
                         + msc_model(model_msc_number,6)*dt**3                                                      ! cubic term 
      else
! PT200830: if an ocean mascon, output a zero value
        if(mcon_prim_usr(imsc,4) > 1000.d0 .and. zero_ocean == "Y")then
          apriori_msc(imsc) = 0.d0
        else
          apriori_msc(imsc) = msc_model_apriori(model_msc_number)
        endif
      endif
    enddo

  endif 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT211019: set to zero any mascons whose value is below the msc_threshold
  do imsc=1,total_prim_usr
    if( (dabs(apriori_msc(imsc))*1.d-3 - msc_threshold)  < 0.001d0)then
! print*,msc_threshold,'set mascon',imsc,' of value ',apriori_msc(imsc), ' to zero',(dabs(apriori_msc(imsc))-msc_threshold)*1.d-3
     apriori_msc(imsc) = 0.d0   ! set 1 mm to avoid roundoff errors
!    else
!print*,'update mascon',imsc,' to value ',apriori_msc(imsc),msc_threshold,(dabs(apriori_msc(imsc)) - msc_threshold)*1.d-3
    endif
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! need to make sure that the a priori model is mass conserving before we write it out to the vcv file ...

! calculate the sum of the volume changes across all mascons. That is, sum(a priori EWH x mascon area) should be zero
! PT211104: this should be done on mass, not on volume
  land_area = 0.d0
  land_vol = 0.d0
  ocean_area = 0.d0
  ocean_vol = 0.d0
  do imsc=1,total_prim_usr
    if(mcon_prim_usr(imsc,4) < 1010.d0)then   ! it is a land primary mascon
      land_area = land_area + mcon_prim_usr(imsc,3)
      land_vol  = land_vol  + mcon_prim_usr(imsc,3)*apriori_msc(imsc)*1.d-3   ! volume in m^3
    else
      ocean_area = ocean_area + mcon_prim_usr(imsc,3)
      ocean_vol  = ocean_vol  + mcon_prim_usr(imsc,3)*apriori_msc(imsc)*1.d-3   ! volume in m^3
    endif
  enddo

! was mass conserved?
  vol_error = land_vol + ocean_vol
  vol_error_orig = vol_error
  mass_error = land_vol*1000.d0 + ocean_vol*1029.d0
  mass_error_orig = mass_error
  write(message,'(a,e15.6,a,f8.4,a)')"Original mass conservation error of ",mass_error*1.d0/(land_area+ocean_area) &
                                     ," kg across Earth, or " &
                                     ,1.d3*mass_error/(1029.d0*ocean_area)," mm GSL"
  call status_update('STATUS','UTIL','apriori_mascons',' ',message,0)

! take the residual volume and assign a correction across all ocean mascons, proportional to the area of each ocean mascon
  allocate(mass_corr(total_prim_usr))
  mass_corr = 0.d0

  do imsc=1,total_prim_usr
    if(mcon_prim_usr(imsc,4) > 1010.d0)then
      mass_corr(imsc) = -1.d0*mass_error/(1029.d0*ocean_area) ! * mcon_prim_usr(imsc,3)/ocean_area
    endif   
    ! write out in vcv format
    write(luout,'(i5,a4,i5.5,a,12x,f17.7,f17.7,f17.7)')imsc,". MC",imsc," (m)",0.d0 &
                ,apriori_msc(imsc)*1.d-3 + mass_corr(imsc),1.d0
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! check that it is now mass conserving
  ocean_vol_orig = ocean_vol
  land_area = 0.d0
  land_vol = 0.d0
  ocean_area = 0.d0
  ocean_vol = 0.d0
  ocean_vol_corr = 0.d0
  do imsc=1,total_prim_usr
    if(mcon_prim_usr(imsc,4) < 1010.d0)then   ! it is a land primary mascon
      land_area = land_area + mcon_prim_usr(imsc,3)
      land_vol  = land_vol  + mcon_prim_usr(imsc,3)*apriori_msc(imsc)*1.d-3   ! volume in m^3
    else
      ocean_area = ocean_area + mcon_prim_usr(imsc,3)
      ocean_vol  = ocean_vol  + mcon_prim_usr(imsc,3)*(apriori_msc(imsc)*1.d-3+mass_corr(imsc))  ! volume in m^3
      ocean_vol_corr = ocean_vol_corr + mcon_prim_usr(imsc,3)*mass_corr(imsc)
!print*,'ocean mass correction:',imsc,mcon_prim_usr(imsc,3)*mass_corr(imsc),mcon_prim_usr(imsc,3),mass_corr(imsc)
    endif
  enddo
! was mass conserved?
  vol_error = land_vol + ocean_vol
  mass_error = land_vol*1000.d0 + ocean_vol*1029.d0

  write(message,'(a,e15.8,a,e15.8,a)')"Adjusted mass conservation error of ",mass_error*1.d0/(land_area+ocean_area) &
                                     ," kg across Earth, or " &
                                     ,1.d3*mass_error/(1029.d0*ocean_area)," mm GSL"
  call status_update('STATUS','UTIL','apriori_mascons',' ',message,0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call status_update('STATUS','UTIL','apriori_mascons',output_file,"Written out apriori  file.",0)

  end
