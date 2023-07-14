  program tidecdftoasc

! program to read the new generation of tide grids - in netcdf format - and output in ascii format ready to be plotted somehow.
!
! P. Tregoning
! 14 September 2016

  use tide_netcdf_mod
  use mascon_mod

  implicit none

  character*150 :: tide_netcdf_file         ! input tide grid file in netcdf format
  character*150 :: ocean_mascon_file        ! ocean mascon file related to the tide netcdf grid
  character*1   :: flag                     ! M for spatial map of one epoch, T for time series at a single point
  integer*4     :: tide_epoch               ! epoch of the file to be read
  character*100 :: arg
  integer*4     :: n_mascons                ! number of ternary mascons in the tide netcdf file

! lat/lon for time series generation
  real(kind=8)        :: lat,lon
  integer*4     :: itern

! variables for getting tide beneath a GRACE groundtrack
  character*150 :: grdtrack_file            ! kbrr file containing epochs and satellite coordinates

! local variables
  real(kind=8)  :: pi,rad_fact
  integer*4     :: num_ocean_tern           ! the pointer to the ocean ternary 
  character*100 :: message,line
  integer*4     :: iepoch                   ! epoch loop counter (goes from 1 to 145)
  logical       :: debug

  pi = 4.0*atan(1.0)
  rad_fact = pi/180.d0

! get command line arguments
  call getarg(1,tide_netcdf_file)
  if(tide_netcdf_file(1:1) == " ")then
    print*,'Runstring: tidecdftoasc tide_netcdf_file ocean_mascon_file M/T epoch or lat/long '
    stop
  endif

! open the input tide grid file (netcdf format)
  tidedata = T_open(tide_netcdf_file)
! get header information from the file
  print *, T_get_info(tidedata, "User")
  print *, T_get_info(tidedata, "Mascons_version")
  print *, T_get_info(tidedata, "Producer")
! get the number of ternary mascons in the tide netcdf file
  n_mascons = T_get_nmascons(tidedata)
  print*,'There are ',n_mascons,' ternary mascons in the tide grid file'

  allocate(ternary_ht1(n_mascons))

! get the name of the ocean mascon file which relates to the tide netcdf file
  call getarg(2,ocean_mascon_file)
  open(22,file=ocean_mascon_file,status='old')
  ! read array sizes from the header
  call read_msc_hdr(22,ocean_mascon_file,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code, &
       max_prim,max_sec,max_tern,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
!  read(22,'(a100)') line
!  write(*,*) line
!  read(line(14:100),*) max_prim,max_sec,max_tern
!  max_tern=1500000
  write(*,*) "array sizes: ",max_prim,max_sec,max_tern

!  allocate the array sizes for the mascons
  call allocate_mascon_arrays

! rewind mascon file because read_ocean_mascons reads the header too  
  rewind(22) 

  ! read in the ocean mascons
  call read_ocean_mascons(22,ocean_mascon_file)


! get the flag of whether we want a spatial map of values or a time series at a single point
  call getarg(3,flag)
  if(flag == "M")then  ! we want a map
    ! next argument is the epoch of the tide grid file that we want to extract (starts at 1, last epoch is 145)
    call getarg(4,arg)
    read(arg,*)tide_epoch
    write(message,'(a,i7)')'plot a map of epoch ',tide_epoch
    call status_update('STATUS','UTIL','tidecdftoasc',ocean_mascon_file,message,0)
  else if (flag == "T")then  ! we want a time series at a point
    tide_epoch = 0
    call getarg(4,arg)
    read(arg,*)lat
    call getarg(5,arg)
    read(arg,*)lon
    write(message,'(a,2f10.4)')'will generate a time series of tide heights at location:',lat,lon
    call status_update('STATUS','UTIL','tidecdftoasc',ocean_mascon_file,message,0)
  else if (flag == "G")then  ! we want a time series of tide heights beneath the groundtrack of a GRACE orbit
    tide_epoch = 0
    call getarg(4,grdtrack_file)
  else
    write(message,'(a)')'second argument must be either M (map), T (time series) or G (groundtrack)'
    call status_update('FATAL','UTIL','tidecdftoasc',ocean_mascon_file,message,0)
    stop
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! global map of ocean heights for the requested epoch
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(flag == "M")then
    call T_read(tidedata,tide_epoch,ternary_ht1)
    do itern=1,max_tern   !n_mascons(tidedata)
      if(mcon_ocean_tern_ptr(itern,2) /= 0) then  ! it is an ocean ternary and should be in the tide netcdf file
        num_ocean_tern = mcon_ocean_tern_ptr(itern,2)
        write(*,'(2f10.4,f8.3,2i9,f10.4,a)')mcon_ocean_tern(num_ocean_tern,1:2)/rad_fact,ternary_ht1(num_ocean_tern)/100. &
                                    ,tide_epoch,itern,mcon_ocean_tern(num_ocean_tern,6)," ocean_ternary"
      endif
    enddo
    stop "Program Finished"
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! time series at a point
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(flag == "T")then
!   find the relevant ocean ternary number for this location. 
!   First calculate which ternary it is in
    ! work out how many ternary latitude bands
    call tern_lat_bands_ell(10.d0/60.d0)
    debug = .false.
    call calc_which_ternary(debug,lat,lon,10.d0/60.d0,itern)
!   Now, which ocean ternary is this?
    num_ocean_tern = mcon_ocean_tern_ptr(itern,2)

! get the time series of tide heights for this mascon
    call T_get_TS(tidedata,num_ocean_tern, ternary_ht1(1:145))  ! there are 145 epochs in the netcdf file, 10 mins apart
    do iepoch = 1,145    
      write(*,'(2f10.4,f8.3,2i9,f10.4,a)')mcon_ocean_tern(num_ocean_tern,1:2)/rad_fact,ternary_ht1(num_ocean_tern)/100. &
                                    ,iepoch,itern,mcon_ocean_tern(num_ocean_tern,6)," ocean_ternary TS"
    enddo
    stop "Extracted time series at a point. Program Finished"
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! time series beneath a groundtrack
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  if(flag == "G")then
!!   we want a tide height at a different location for each epoch
!    open(24,file=grdtrack_file,status='old')
!!   skip over the header lines of the input plt*.kb file
!    do i=1,7
!      read(24,'(a)')line
!    enddo
!
!! now, loop over all epochs in the kbrr file, get the GRACE A lat/long and then get the tide height for this epoch and lat/lon
!   ioerr = 0
!  do while (ioerr == 0)
!    read(24,*)
!  enddo !!! THIS ISN"T FINISHED !!!
!  
!! get tide height from the netcdf file
!    ternary_ht1(i) = T_get_msc_time(tidedata,num_ocean_tern, mjd,secs)
!
!  endif



  end





