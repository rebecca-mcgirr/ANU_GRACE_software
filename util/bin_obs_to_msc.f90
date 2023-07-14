  program bin_obs_to_msc

! Given a list of lat and lon and a mascon file calculate number of satellite passes over each mascon, based on ~/gt/util/which_mascon
! To be used in conjunction with ~/gt/util/xyz_to_llr which converts cartesian coordinates to spherical lat, lon, radius
!
! R. McGirr
! 8 February 2021

  use mascon_mod

  implicit none

  character*50  :: arg
  character*150 :: coord_file
  character*150 :: mascon_file
  integer*4     :: tern_number
  integer*4     :: nobs   ! reader from header of input file?
  integer*4     :: mission
  
  real(kind=8), allocatable      :: lat(:) ! latitudes read from input file
  real(kind=8), allocatable      :: lon(:) ! longitudes read from input file
  real(kind=8), allocatable      :: bin_prim(:,:) ! will contain the lat lon of each prim and the number of sat passes
  integer*4, allocatable      :: ll_to_prim(:) ! will contain the prim location of each observation

  real(kind=8) :: pi
  real(kind=8) :: lat_spacing
  real(kind=8) :: percent_passed
  character*150 :: message
  character*100 :: junk
  integer*4     :: i,imsc,step
  integer*4     :: prim_passed
! define pi
  pi = 4.d0*atan(1.d0)

! Get the command line information
  call getarg(1,arg)
  if(arg(1:1) == " ")then
    print*,"Runstring: bin_obs_to_msc groundtrack.latlon nobs mascon_file"
    stop
  endif
  
  read(arg,*)coord_file
  call getarg(2,arg)
  read(arg,*)nobs
  call getarg(3,mascon_file)
  call getarg(4,arg)
  read(arg,*)mission
  
! allocate arrays
  allocate(lat(nobs))
  allocate(lon(nobs))
  allocate(ll_to_prim(nobs))

! Read in coord_file here
  open(15,file=coord_file,status='old')
  read(15,*)junk
  do i=1,nobs
    read(15,*)lat(i),lon(i)
  enddo
  close(15)

! Read in mascon file 
  call read_msc_hdr(20,mascon_file,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                                ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
  call allocate_mascon_arrays
  call read_mascon_file(20,mascon_file) 
  !print*,"DEBUG max_prim: ",max_prim

! allocate arrays
  allocate(bin_prim(max_prim,3))

! set the ternary latitudinal spacing
  lat_spacing = 10.d0/60.d0     ! default is 10' ternary masons in latitude

! establish the array that contains how many ternary latitude bands and how many ternary mascons per band
  call tern_lat_bands_ell(lat_spacing)

! loop through each coordinate for each observation, find corresponding ternary which points to its primary
  if (mission == 0) then
    step = 1
  else if (mission == 1)then
    step = 5
  else
    STOP "MISSION NUMBER NOT RECOGNISED"
  endif

  do i = 1,nobs,step
  ! find which ternary it is in
    call calc_which_ternary(.false.,lat(i),lon(i),lat_spacing,tern_number)

    if(mcon_tern_ptr(tern_number,1) == 0)then
      write(message,'(a,2f11.4)')"Error: there is no ternary mascon in the input mascon file at location:",lat(i),lon(i)
      call status_update('FATAL','UTIL','bin_obs_to_msc',' ',message,0)
    endif
  ! convert to primary number
    ll_to_prim(i) = mcon_tern_ptr(tern_number,1)
  enddo

! loop through each mascon, and count number of times that mascon occurs in ll_to_prim
  bin_prim(:,1) = mcon_prim(:,1)*180.d0/pi ! copy lat
  bin_prim(:,2) = mcon_prim(:,2)*180.d0/pi ! copy lon
  bin_prim(:,3) = 0.d0 ! initialise bin count to 0 for all mascons
  !print*,"NOBS: ",nobs
  do i = 1,nobs,step
    bin_prim(ll_to_prim(i),3) = bin_prim(ll_to_prim(i),3) + 1.d0
    !if(bin_prim(ll_to_prim(i),3) > 0)print*,bin_prim(ll_to_prim(i),:),ll_to_prim(i)
  enddo

! calculate percentage of mascons passed
  prim_passed = 0
  do i = 1,max_prim
    if (bin_prim(i,3) > 0) then
      prim_passed = prim_passed + 1
    endif
  enddo
  percent_passed = float(prim_passed) / float(max_prim) * 100.d0

  write(message,'(a,i6,a,i6,a,f8.4,a)')"Number of mascons observed: ",prim_passed," of ",max_prim,", or ",percent_passed &
                                         ,"% of mascons observed"
  call status_update('STATUS','UTIL','bin_obs_to_msc',' ',message,0)

! write out 
  open(26,file='bin_obs_to_msc.out',status='unknown')
  write(26,'(a)')"# prim_passed, max_prim, perecent_passed"
  write(26,'(i6,i6,f9.4,a)')prim_passed,max_prim,percent_passed,"%"
  write(26,'(a)')"# lat, lon, no. passes"
  do i = 1,max_prim
    write(26,*)bin_prim(i,:)
  enddo
  close(26)
  end
