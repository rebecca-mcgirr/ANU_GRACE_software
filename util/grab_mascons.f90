  program grab_mascons

! program to identify all the ternary mascons within a particular, user-defined polygon region, extract them from their
! current primary mascons and put them all into a newly created primary mascon. This will be useful for creating primary
! mascons of particular drainage basins.
!
! Many of the subroutines and some of the logic replicates what is done in util/reshape_mascons
!
! P. Tregoning
! 18 December 2016
!
! MODS:
! PT170209: allow over-ride of density of mascons within and without through command line control
! PT170420: allow a region to be specified by coords of a rectangle rather than by a polygon (to enable me to split up the Shelf region)
! PT190822: make option to just output the primary mascon numbers contained with the region
! 
! PT200129: major change in logic. Just check whether each ternary is within the region. If so, use it. Don't worry about checking
!           primary mascon CoM anymore - it gives the wrong results sometimes. Forced this to happen by setting in_or_out to "1" in all cases.

  use mascon_mod     ! provides the declarations of the variables needed to read in and store all the mascons

  implicit none

! command line arguments
  character    :: mascon_file*150               ! mascon input file
  character    :: msc_file_out*150              ! mascon output file
  character    :: polygon_file*150              ! file containing the polygon of the required region
  character    :: region*15                     ! region to re-define (PineIs, Thwaites, Totten, Helheim etc)
  real(kind=8) :: density_in,density_out        ! over-riding controls on densities required for ternarys within/outside region
  character    :: tern_type*4                   ! flag whether we want to grab land or ocean ternarys

  integer*4    :: ioerr
  integer*4    :: imsc,itern

! parameters related to defining the geographical region
  real(kind=8), allocatable :: region_poly(:,:)  ! vertices of the polygon region
  integer*4                 :: npoints           ! number of polygon vertices
  integer*4                 :: in_or_out         ! 1 if point in polygon, 0 if on the edge of the polygon, -1 if out of polygon
  integer*4                 :: nmsc_region       ! number of primary mascons found within the requested region
  integer*4, allocatable    :: region_prims(:)   ! list of primary mascon numbers found within the requested region
  real(kind=8)              :: region_area       ! sum of all ternary mascon areas found within the requested region
  logical                   :: cross_0E          ! logical to know whether to translate the region polygon so that it no longer crosses 0E

! parameters for counting changes in mascon numbers etc
  integer*4  :: new_tern_per_prim ! new maximum number of ternary mascons per primary
  integer*4  :: new_tern_per_sec  ! new maximum number of ternary mascons per secondary
  integer*4  :: new_max_tern_per_prim ! maximum number of ternary mascons per primary (and, hence, per secondary) in the reshaped file

! variables to hold the reshaped mascons
  real(kind=8),allocatable  :: msc_reshaped_prims(:,:)            ! the reshaped primary mascons
  real(kind=8),allocatable  :: msc_tern_reshaped(:,:,:)     ! the reshaped ternary mascons
  integer*4                 :: ntern_reshape                ! number of ternary mascons to be reshaped 

! variables related to outputting mascons using write_mascon_record
  real(kind=8) :: msc_crds(3)       ! temp storage of mascon lat/lon/rad
  integer*4    :: colours(2)        ! colour of secondary and primary mascon
  integer*4    :: pointers(2)       ! tern2prim and tern2sec pointers
  integer*4    :: nmsc_current      ! (re)counting of the primary mascons as we output them
  integer*4    :: nsec,ntern        ! temporary value for the secondary mascon number and the ternary mascon number
  real(kind=8) :: col_range(2)      ! variable used to calculate new colour for primary mascons
  integer*4    :: prim_colour       ! new colour for primary mascon after reshaping
  character*6  :: sec_flag          ! single, temporary value assigned to the type of secondary mascon
  integer*4    :: tern_number       ! temporary value of the actual ternary mascon number
  character*12 :: hashcode          ! temporary variable used to calculate the new character code for the output mascon file


! unit numbers
  integer*4, parameter      :: lumsc_in  = 10       ! unit number of  input mascon file
  integer*4, parameter      :: lumsc_out = 11       ! unit number of output mascon file
  integer*4, parameter      :: lupoly    = 12       ! unit number of file containing polygons around continents etc

! counters and logicals
  integer*4 :: i,j,k
  integer*4 :: nfound                            ! count how many of the requested primary mascons have been found in the ternary file
  logical   :: use_primary     

! other stuff
  integer*4 :: indx, trimlen
  character :: line*200, message*250,arg*100
  logical   :: debug,debug_rms
  real(kind=8)         :: pi

! variables for writing date and time into header of output file
  character(8)  :: date
  character(10) :: time
  character(5)  :: timezone

! PT161218: variables to densify the points along the polygon that defines the region
  integer*4     :: npoints_dense          ! number of points required along the polygon so that they are no further apart than 25 km
  real(kind=8)  :: region_poly_dense(100000,2)   ! array to hold the dense polygon coordinates
  integer*4     :: region_prims_dense(100000)    ! array to hold the primary mascon numbers that reside within the dense polygon
  real(kind=8)  :: gcirc                  ! great circle angle (in radians) between two coordinates

! variables to grab the ternarys from existing masons and move them
  logical, allocatable :: remove_tern(:,:)
  integer*4            :: all_in
  real(kind=8)         :: dlon
  integer*4            :: n_reshaped,n_to_reshape
  integer*4            :: n_other_density               ! number of ternarys in an original mascon that don't meet criteria to be grabbed
  integer*4            :: new_n_prim,new_n_sec          ! new numbers of primary and secondary mascons ???
  integer*4            :: points_needed                 ! number of points to densify between two polygon vertices
  real(kind=8),allocatable :: mcon_tern_inside(:,:,:)   ! array to hold the ternary mascons that are found to be inside the polygon
  real(kind=8),allocatable :: mcon_tern_tmp(:,:,:)
  character*6          :: sec_flag_tmp                  ! temporary variable to construct the secondary mascon flag based on the primary mascon

! variables to copy some of the original mcon_prim information
  real(kind=8), allocatable :: mcon_prim_orig(:)        ! keep a copy of the original numbers of ternary mascons in each primary

! some values
  real(kind=8),parameter :: earthrad = 6378137.45d0     ! radius of the Earth
  logical      :: found
  real(kind=8) :: req_area                              ! required area for a primary mascon ???
  real(kind=8) :: rms_limit
  real(kind=8) :: tern_spacing                          ! spacing in latitude (in decimal degrees) of ternary mascons
  integer*4    :: n_tern_tmp
  integer*4    :: ntern_inside                          ! total number of ternary mascons found within the polygon region
  integer*4    :: ntern_inside_tmp
  integer*4    :: not_in_polygon
  logical      :: only_prim_nums                        ! logical to output just the primary mascon numbers contained within the region

! PT170420: coords of a rectangle to describe a region
  real(kind=8) :: minlat,minlon,maxlat,maxlon

  pi = 4.d0*atan(1.d0)
  tern_spacing = 10.d0/60.d0      ! we use 10' ternary mascons

  new_max_tern_per_prim = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  get the command line and cmdfile information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call getarg(1,mascon_file)
  if(mascon_file(1:1) == " ")then
    call status_update('FATAL','UTIL','grab_mascons',' ' &
              ,'Runstring: grab_mascons mascon_file_in mascon_file_out polygon_file region Land/Ocea rms_limit',0)
  else
    call getarg(2,msc_file_out)
    region = " "
    call getarg(3,polygon_file)

! PT170420: check whether reading the name of a polygon file or a flag that we want to break up the continental shelf mascon
    if(polygon_file(1:5) /= "Shelf")then
      call getarg(4,region)
      call getarg(5,tern_type)
      call getarg(6,arg)
      read(arg,*)rms_limit

! PT190822: use a negative rms_limit as a flag to output only the primary mascons contained within the region
      if(rms_limit < 0.d0)then
        only_prim_nums = .true.
        write(message,'(a,a,a)')'Will output only primary mascon numbers contained in region "',region,'"'
        call status_update('STATUS','UTIL','grab_mascons',mascon_file,message,0)
      else
        only_prim_nums = .false.
        write(message,'(a,a,a)')'Will create primary mascon file for region "',region,'"'
        call status_update('STATUS','UTIL','grab_mascons',mascon_file,message,0)
      endif
    else
      call getarg(4,arg)
      read(arg,*)minlat
      call getarg(5,arg)
      read(arg,*)minlon
      call getarg(6,arg)
      read(arg,*)maxlat
      call getarg(7,arg)
      read(arg,*)maxlon
      call getarg(8,tern_type)
      call getarg(9,arg)
      read(arg,*)rms_limit

      write(message,'(a,4f8.1,a,f6.0)')'Will create primary mascon file for continental shelf region "' &
        ,minlat,minlon,maxlat,maxlon,'". Including all ternarys with density: ',density_in
      call status_update('STATUS','UTIL','grab_mascons',mascon_file,message,0)
    endif
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  open the input mascon file, read the header, allocate arrays and read in the mascons
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call read_msc_hdr(lumsc_in,mascon_file,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                                ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
  call allocate_mascon_arrays
  call read_mascon_file(lumsc_in,mascon_file)

! allocate a variable and copy the original number of ternary mascons in the primary mascons
  allocate(mcon_prim_orig(max_prim))
  mcon_prim_orig(:) = mcon_prim(:,8)
  allocate(mcon_tern_tmp(1,1000000,nvar_tern))     ! we don't know how large the polygon area is yet, so don't know how many ternarys there might be.
  mcon_tern_tmp = 0.d0

! set up the ternary latitude bands
  call tern_lat_bands_ell(10.d0/60.d0)


! open the output file
  open(lumsc_out,file=msc_file_out,status='unknown')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! allocate some arrays
  allocate(region_prims(total_prim))    ! make it as big as it could possibly need to be.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! define coordinates of a polygon that defines the requested region
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(polygon_file(1:5) /= "Shelf" .and. polygon_file /= "SthnO" .and. polygon_file(1:5) /= "Antar" .and. polygon_file(1:5) /= "Arcti")then

    open(lupoly,file=polygon_file,status='old')

! read the file to find the line for the region of interest
    line = " "
    npoints = -999
    cross_0E = .false.
    do while (line(9:19) /= region)
      read(lupoly,'(a)',iostat=ioerr,end=1000)line
      if(ioerr == 0 .and. line(9:19) == region) then
        read(line(1:6),*)npoints
        allocate(region_poly(npoints,2))
! PT170608: use a value from the command line instead
!        read(line(20:40),*)rms_limit
        do i=1,npoints
          read(lupoly,*)region_poly(i,:)
          if(region_poly(i,2) < 0.d0 )cross_0E = .true.
        enddo
        write(message,'(a,a,a,i6,a)')"Polygon around region ",region," has ",npoints," vertices"
        call status_update('STATUS','UTIL','grab_mascons','mascon.polygons',message,0)
print*,"cross_0E: ",cross_0E

! PT161118: check whether there are any -ive longitudes in the polygon vertices. If so, the polygon
!           crosses 0E. In this case, shift it 90degrees to the east so that the in_polygon algorithm
!           will work.
        if(cross_0E)region_poly(:,2) = region_poly(:,2) + 90.d0

      endif
    enddo

1000 continue


    if(npoints == -999)then
      write(message,'(a,a,a)')"End of reading file and region ",region," not found"
      call status_update('FATAL','UTIL','grab_mascons',polygon_file,message,0)
    endif

  elseif(polygon_file(1:5) == "Shelf")then
    !!!!!!!! shelf region to be defined by a rectangle !!!!
    print*,'Fill in the region_poly variable for the Shelf rectangle region'
    npoints = 5
    allocate(region_poly(npoints,2))
    rms_limit = 1000.d0
    cross_0E = .false.
    region_poly(1,1) = minlat
    region_poly(1,2) = minlon
    region_poly(2,1) = maxlat
    region_poly(2,2) = minlon
    region_poly(3,1) = maxlat
    region_poly(3,2) = maxlon
    region_poly(4,1) = minlat
    region_poly(4,2) = maxlon
    region_poly(5,1) = minlat
    region_poly(5,2) = minlon
    do i=1,npoints
      if(region_poly(i,2) < 0.d0 )cross_0E = .true.
    enddo
    if(cross_0E)region_poly(:,2) = region_poly(:,2) + 90.d0
  else
      write(message,'(a,a,a)')"No polygon for region ",region,". Dealt with explicitly in the code."
      call status_update('STATUS','UTIL','grab_mascons','mascon.polygons',message,0)
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ensure that the vertices of the polygon are no further apart than what might be a typical
! high-resolution mascon size. I'm going to make it 25 km, which should be smaller than we ever need.
! By making it very short, we should ensure that all primary mascons partly in/partly out of the 
! region are captured by the code below and included in the list of mascons that donate ternary
! mascons to the new region.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(region(1:5) /= "SthnO" .and. region(1:5) /= "Antar" .and. region(1:5) /= "Arcti")then
    npoints_dense = 0

! PT170118: convert the polygon region coords to radians
    region_poly(:,1:2) = region_poly(:,1:2)*pi/180.d0

    do i=2,npoints

! add the polygon vertex into the dense polygon coord list
      npoints_dense = npoints_dense + 1
      region_poly_dense(npoints_dense,1:2) = region_poly(i-1,1:2)

! calculate the great circle distance between this point and the next one
      dlon = region_poly(i,2) - region_poly(i-1,2)
      gcirc = dacos( dcos(region_poly(i-1,1))*dcos(region_poly(i,1))  &
             +dsin(region_poly(i-1,1))*dsin(region_poly(i,1))*dcos(dlon) )


! add more points in between if necessary
      if(gcirc > 10.d3/earthrad) then       ! it is more than 10 km along the great circle between the two points

      ! we need to divide it up into a more dense set of points along the line. 
        points_needed = int(gcirc/(10.d3/earthrad))

      ! now just linearly interpolate the coords from the start to the end to infill with this many points
        do j=1,points_needed
          npoints_dense = npoints_dense + 1
          region_poly_dense(npoints_dense,1) = region_poly(i-1,1)+j*(region_poly(i,1)-region_poly(i-1,1))/(dble(1+points_needed))
          region_poly_dense(npoints_dense,2) = region_poly(i-1,2) + j*dlon/(dble(1+points_needed))
        enddo
      else
        npoints_dense = npoints_dense + 1
        region_poly_dense(npoints_dense,1:2) = region_poly(i,1:2)
      endif
    enddo

! add the last polygon vertex into the dense polygon coord list
    npoints_dense = npoints_dense + 1
    region_poly_dense(npoints_dense,1:2) = region_poly(npoints,1:2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now check the segment from the last back to the first vertex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dlon = region_poly(1,2) - region_poly(npoints,2)
    gcirc = dacos( dcos(region_poly(npoints,1))*dcos(region_poly(1,1))  &
             +dsin(region_poly(npoints,1))*dsin(region_poly(1,1))*dcos(dlon) )
! add more points in between if necessary
    if(gcirc > 10.d3/earthrad) then       ! it is more than 25 km along the great circle between the two points
      ! we need to divide it up into a more dense set of points along the line. 
      points_needed = int(gcirc/(10.d3/earthrad))

      ! now just linearly interpolate the coords from the start to the end to infill with this many points
      do j=1,points_needed
        npoints_dense = npoints_dense + 1
        region_poly_dense(npoints_dense,1) = region_poly(npoints,1) &
                                           + j*(region_poly(1,1)-region_poly(npoints,1))/(dble(1+points_needed))
        region_poly_dense(npoints_dense,2) = region_poly(npoints,2) + j*dlon/(dble(1+points_needed))
      enddo
! add the first polygon vertex into the dense polygon coord list
      npoints_dense = npoints_dense + 1
      region_poly_dense(npoints_dense,1:2) = region_poly(1,1:2)

    else
      npoints_dense = npoints_dense + 1
      region_poly_dense(npoints_dense,1:2) = region_poly(1,1:2)
    endif

! PT170118: convert region_poly_dense from radians back to decimal degrees
    region_poly_dense(:,1:2) = region_poly_dense(:,1:2)*180.d0/pi
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! loop over all polygon coordinates of the region of interest and identify in which primary mascon
! they reside. This will form part of the list of primary mascons from which to extract the ternary
! mascons for the new region.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  n_to_reshape = 0
  ntern_reshape = 0
  do i = 1,npoints_dense
    debug = .false.
    call calc_which_ternary(debug,region_poly_dense(i,1),region_poly_dense(i,2),tern_spacing,tern_number)
    imsc = mcon_tern_ptr(tern_number,1)

    ! now we have identified the primary in which this ternary resides. Check whether it is already in the list.
    found = .false.
    do j=1,n_to_reshape
      if(region_prims_dense(j) == imsc)found = .true.
    enddo

    if(.not.found)then
      n_to_reshape = n_to_reshape + 1
      region_prims_dense(n_to_reshape) = imsc
      ntern_reshape = ntern_reshape + nint(mcon_prim(imsc,8))

!        if(only_prim_nums)then
!          if(mcon_prim(imsc,6) < 1010.d0 .and. ( tern_type(1:4) == "Land" .or. tern_type(1:3) == "All") )then
!            print*,"Primary mascon",imsc," is in region. Area =",mcon_prim(imsc,4),mcon_prim(imsc,6)
!          else if (mcon_prim(imsc,6) > 1010.d0 .and. (tern_type(1:4) == "Ocea" .or. tern_type(1:3) == "All") )then
!            print*,"Primary mascon",imsc," is in region. Area =",mcon_prim(imsc,4),mcon_prim(imsc,6)
!          endif
!        endif

    endif
  enddo

  write(message,'(a,i6,a)')"Found ",n_to_reshape," primary mascons containing vertices of the region polygon"
  if(npoints_dense >  0)call status_update('STATUS','UTIL','grab_mascons',' ',message,0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! loop over all primary mascons and identify the ones within the region of interest that might not
! cross the polygon line (ie are contained wholly within the region
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  n_to_reshape  = 0
!  ntern_reshape = 0
!  region_area = 0.d0
!  new_max_tern_per_prim = 0

  all_in = 0
  do imsc = 1,total_prim

! is this mascon in the region of interest?
    if(region(1:5) /= "SthnO" .and. region(1:5) /= "Antar" .and. region(1:5) /= "Arcti")then
      call msc_in_region(cross_0E,mcon_prim(imsc,1)*180.d0/pi,mcon_prim(imsc,2)*180.d0/pi,region_poly_dense(:,1) &
           ,region_poly_dense(:,2),npoints_dense,in_or_out)
    elseif(region(1:5) == "SthnO")then
      if(mcon_prim(imsc,1)*180.d0/pi < -60.d0 .and. mcon_prim(imsc,6) > 1010.d0 .or. prim_flags(imsc) == "POcean" )then
        in_or_out = 1
      else
        in_or_out = 0
      endif
    elseif(region(1:5) == "Antar")then
      if(mcon_prim(imsc,1)*180.d0/pi < -60.d0 .and. mcon_prim(imsc,6) < 1010.d0 )then
        in_or_out = 1
      else
        in_or_out = 0
      endif
    elseif(region(1:5) == "Arcti")then
      if(mcon_prim(imsc,1)*180.d0/pi > 65.d0 .or. ( mcon_prim(imsc,1)*180.d0/pi > 56.d0 .and. (mcon_prim(imsc,2)*180.d0/pi > 295.d0 &
                         .or. mcon_prim(imsc,2)*180.d0/pi < 17.d0)).and. mcon_prim(imsc,6) > 1010.d0)then
        in_or_out = 1
      !elseif(mcon_prim(imsc,1)*180.d0/pi > 50.d0 .and. prim_flags(imsc) == "PGrndO")then
      !  print*,"Found prim PGrndO"
      !  in_or_out = 1
      else
        in_or_out = 0
      endif    
    endif

! PT200129: just force the primary mascon to be considered
    !in_or_out = 1

    if(in_or_out == 1)then

! Check whether it is already in the list.
      found = .false.
      do j=1,n_to_reshape
        if(region_prims_dense(j) == imsc)then
          found = .true.
        endif
      enddo

      if(.not.found)then
        all_in = all_in + 1
        n_to_reshape  = n_to_reshape + 1
        region_prims_dense(n_to_reshape) = imsc
        ntern_reshape = ntern_reshape + nint(mcon_prim(imsc,8))
        prim_flags(imsc) = "P-----"
!PT190822: output the primary mascon number
!        if(only_prim_nums)then
!          if(mcon_prim(imsc,6) < 1010.d0 .and. ( tern_type(1:4) == "Land" .or. tern_type(1:3) == "All") )then
!            print*,"Primary mascon",imsc," is in region. Area =",mcon_prim(imsc,4),mcon_prim(imsc,6)
!          else if (mcon_prim(imsc,6) > 1010.d0 .and. (tern_type(1:4) == "Ocea" .or. tern_type(1:3) == "All") )then
!            print*,"Primary mascon",imsc," is in region. Area =",mcon_prim(imsc,4),mcon_prim(imsc,6)
!          endif
!        endif
      endif
    endif
  enddo

  if(only_prim_nums)then
    do j=1,n_to_reshape
      imsc = region_prims_dense(j)
          if(mcon_prim(imsc,6) < 1010.d0 .and. ( tern_type(1:4) == "Land" .or. tern_type(1:3) == "All") )then
            print*,"Land primary mascon",imsc," is in region. Area =",mcon_prim(imsc,4),mcon_prim(imsc,6),j &
                  ,mcon_prim(imsc,1:2)*180.d0/pi,mcon_prim(imsc,6)
          else if (mcon_prim(imsc,6) > 1010.d0 .and. (tern_type(1:4) == "Ocea" .or. tern_type(1:3) == "All") )then
            print*,"Ocean primary mascon",imsc," is in region. Area =",mcon_prim(imsc,4),mcon_prim(imsc,6),j &
                  ,mcon_prim(imsc,1:2)*180.d0/pi,mcon_prim(imsc,6)

          endif
    enddo

    print*,"All done. Grab_mascons can stop now."
    stop
  endif

! allocate an array which will be used to store whether to remove or keep a ternary in its original primary mascon
  allocate(remove_tern(n_to_reshape,max_tern_per_prim))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                              !!
!!          separate ternary mascons as being within or outside the requested polygon           !!
!!                                                                                              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  remove_tern(:,:) = .false.
  ntern_inside = 0
  req_area = 0.d0
  n_other_density = 0
  do i=1,n_to_reshape
    imsc = region_prims_dense(i)
    n_tern_tmp = 0
    do itern=1,nint(mcon_prim(imsc,8))
      if(region(1:5) /= "SthnO" .and. region(1:5) /= "Antar" .and. region(1:5) /= "Arcti")then
        call msc_in_region(cross_0E,mcon_tern(imsc,itern,1)*180.d0/pi,mcon_tern(imsc,itern,2)*180.d0/pi,region_poly_dense(:,1) &
                           ,region_poly_dense(:,2),npoints_dense,in_or_out)

! PT170607: special case for the southern ocean
      elseif(region == "SthnO")then
        if( (mcon_tern(imsc,itern,1)*180.d0/pi < -60.d0) .and. prim_flags(imsc) /= "PShelf" &
                                                         .and. mcon_tern(imsc,itern,6)> 1010.d0 )then
!print*,'want ternary:',mcon_tern(imsc,itern,1)*180.d0/pi,mcon_tern(imsc,itern,2)*180.d0/pi
          in_or_out = 1
        else
          in_or_out = 0
        endif
      elseif(region == "Antarctica")then
        if( (mcon_tern(imsc,itern,1)*180.d0/pi < -60.d0) .and. prim_flags(imsc) /= "PShelf" &
                                                         .and. mcon_tern(imsc,itern,6) < 1010.d0 )then
!print*,'want ternary:',mcon_tern(imsc,itern,1)*180.d0/pi,mcon_tern(imsc,itern,2)*180.d0/pi
          in_or_out = 1
        else
          in_or_out = 0
        endif
      elseif(region == "Arctic")then
        if( mcon_tern(imsc,itern,1)*180.d0/pi > 65.d0 .or. (mcon_tern(imsc,itern,1)*180.d0/pi > 56.d0 &
                       .and. (mcon_tern(imsc,itern,2)*180.d0/pi > 295.d0 .or. mcon_tern(imsc,itern,2)*180.d0/pi < 17.d0)) &
                       .and. prim_flags(imsc) /= "PShelf" .and. mcon_tern(imsc,itern,6) > 1010.d0 )then

!print*,'want ternary:',mcon_tern(imsc,itern,1)*180.d0/pi,mcon_tern(imsc,itern,2)*180.d0/pi
          in_or_out = 1
          !print*,(mcon_tern(imsc,itern,2)*180.d0/pi)
        !elseif(mcon_tern(imsc,itern,1)*180.d0/pi > 50.d0)then 
        !  print*,"Found tern PGrndO"
        !  in_or_out = 1
        else
          in_or_out = 0
        endif
      endif
      if(in_or_out == 1 )then        ! it is inside the polygon 
!if(imsc == 2933)print*,'imsc,itern,crds:',imsc,itern,mcon_tern(imsc,itern,1)*180.d0/pi,mcon_tern(imsc,itern,2)*180.d0/pi &
!               ," ",tern_type,mcon_tern(imsc,itern,6)
! PT170607: include only ternarys with an appropriate density
        if((tern_type == "Land" .and. mcon_tern(imsc,itern,6) < 1010.) .or.  &
                       (tern_type == "Ocea" .and. mcon_tern(imsc,itern,6) > 1010.) .or. &
                        tern_type(1:3) == "All")then
          remove_tern(i,itern) = .true.
          ntern_inside = ntern_inside + 1
          !PT170120: add the ternary information to a temporary array for the ternarys inside the region
          mcon_tern_tmp(1,ntern_inside,:) = mcon_tern(imsc,itern,:)
          ! PT170119: add the ternary area to the total for the polygon region
          req_area = req_area + mcon_tern(imsc,itern,4)
        else 
          n_tern_tmp = n_tern_tmp + 1
          mcon_tern(imsc,n_tern_tmp,:) = mcon_tern(imsc,itern,:)
          ntern_reshape = ntern_reshape - 1
        endif
          
      else
! PT170119: we now shuffle the ternarys external to the polygon to make the new external primary mascon
!print*,'imsc,itern not what we want:',imsc,itern,mcon_tern(imsc,itern,1)*180.d0/pi,mcon_tern(imsc,itern,2)*180.d0/pi &
!             ,mcon_tern(imsc,itern,6)
        n_tern_tmp = n_tern_tmp + 1
        mcon_tern(imsc,n_tern_tmp,:) = mcon_tern(imsc,itern,:)
        ntern_reshape = ntern_reshape - 1
      endif
    enddo
    if(nint(mcon_prim(imsc,8))-n_tern_tmp > 0)then
      write(message,'(a,i7,a,i7,2a,i7,a)')'Primary mascon',imsc,' has ',n_tern_tmp,' ternary mascons outside the region' &
          ,' (or wrong density) and ',nint(mcon_prim(imsc,8))-n_tern_tmp,' ternary mascons within the region with correct density'
      call status_update('STATUS','UTIL','grab_mascons',' ',message,0)
    endif

! PT170119: update the number of ternarys remaining in the external primary
    mcon_prim(imsc,8) = n_tern_tmp
    mcon_sec(imsc,1,8) = mcon_prim(imsc,8)

! RM200114: count the original primary mascons that are either not required density or partly outside the region
    if( n_tern_tmp > 0)n_other_density = n_other_density + 1

  enddo
!stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                                          !! 
!!     E S T A B L I S H     A     S I N G L E    M A S C O N     I N S I D E    T H E    P O L Y G O N     !!
!!                                                                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! allocate an array to hold all the ternary mascons inside the polygon
  allocate(mcon_tern_inside(1,ntern_inside,nvar_tern))
  write(message,'(a,i8,a)')"Single primary mascon inside the polygon has",ntern_inside,' ternary mascons'
  call status_update('STATUS','UTIL','grab_mascons',msc_file_out,message,0)

! update the max_tern_per_prim variable if necessary
  if(ntern_inside > max_tern_per_prim)max_tern_per_prim = ntern_inside

! fill out the mcon_tern_inside array which is properly dimensioned for the number of ternary mascons
  ntern_inside_tmp = 0
  do itern = 1,ntern_inside
    mcon_tern_inside(1,itern,:) = mcon_tern_tmp(1,itern,:)  
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                          O U T P U T    H E A D E R    I N F O R M A T I O N
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write out the first header line, containing the code and the numbers of mascons etc
! RM200114: the new number of primary mascons is the total minus only those with the right density and completely inside the region, plus the new mascon
  new_n_prim = total_prim + 1 - (n_to_reshape - n_other_density)  ! we have added one new mascon .... but may need to remove any that were wholly within the required region
  new_n_sec  = new_n_prim

! add a new line to the header to represent what has been done in running this program
  call date_and_time(date,time,timezone)
  write(msc_hdr_lines(n_hdr_lines+1),'(a,a,a,i6,a,e10.4,a,1x,a8,1x,a10,1x,a5)') &
       "# grab_mascons   , region: ",region," created and contains",ntern_reshape &
       ," ternary mascons. Area: ~",req_area/1.d6," km^2. Timetag: ",date,time,timezone
! generate new unique code for the output program
  call hash_header(msc_hdr_lines(2:n_hdr_lines+1),n_hdr_lines+1,hashcode)

! write the new code and mascon numbers to the first line of the header
  write(msc_hdr_lines(1)(1:81),'(a,i6,5i8)')hashcode(1:8),new_n_prim,new_n_sec,max_tern,max_sec_per_prim &
                                           ,max_tern_per_prim,max_tern_per_prim

! and write out all header lines to the output file
  do i=1,n_hdr_lines+1
    write(lumsc_out,'(a)')msc_hdr_lines(i)
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                              !! 
!!             O U T P U T     M A S C O N S    O U T S I D E    T H E    P O L Y G O N         !!
!!                                                                                              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nmsc_current = 0
  call status_update('STATUS','UTIL','grab_mascons',msc_file_out,"Write out mascons outside of the polygon",0)
  do imsc = 1,total_prim 
    if(nint(mcon_prim(imsc,8)) > 0)then

! PT170119: we should really recompute the CoM of the primary mascons, but I haven't done it yet !!!

      nmsc_current = nmsc_current + 1
!     primary record
      nsec = 1
      msc_crds(1:2) = mcon_prim(imsc,1:2)*180.d0/pi
      msc_crds(3)   = mcon_prim(imsc,3)
      pointers = nmsc_current
      colours(1:2) = 0
! PT190306: Reset the primary mascon flag using the ternary flag, for any primary mascon that was partially outside the polygon
!          [Note: we should recompute the CoM of the primary mascon here as well .... but that is for another day! ]
      if(prim_flags(imsc)(2:4) == "---")then
        prim_flags(imsc)(2:6) = tern_flag(nint(mcon_tern(imsc,1,7)))(2:6)
      endif
      call write_mascon_record(lumsc_out,prim_flags(imsc),nmsc_current,nsec,nint(mcon_prim(imsc,8)),msc_crds     &
                              ,mcon_prim(imsc,4),mcon_prim(imsc,5),mcon_prim(imsc,12),mcon_prim(imsc,6),mcon_prim(imsc,11) &
                              ,mcon_prim(imsc,10),mcon_region(imsc),pointers,colours)
!     secondary record for this primary. We still can only have one secondary per primary at this point.
      pointers(1) = nmsc_current
      pointers(2) = 0
! PT161121: we don't want this, we want the secondary mascons to be sequential along with the primary mascons
!      nsec = mcon_sec(imsc,1,7)           ! this is the unique secondary mascon number as read from the input mascon file
!      colours = sec_colour(nsec)

! PT170426: change secondary flag to be exactly the primary flag 
      sec_flag_tmp = "S"//prim_flags(imsc)(2:6)
      call write_mascon_record(lumsc_out,sec_flag_tmp,nmsc_current,nmsc_current,nint(mcon_prim(imsc,8)),msc_crds &
                              ,mcon_prim(imsc,4),mcon_prim(imsc,5),mcon_prim(imsc,12),mcon_prim(imsc,6),mcon_prim(imsc,11) &
                              ,mcon_prim(imsc,10),mcon_region(imsc),pointers,colours)
!     ternary records for this secondary
      do itern = 1,nint(mcon_sec(imsc,1,8))
        pointers(1) = nmsc_current
        pointers(2) = nsec
        ntern = nint(mcon_tern(imsc,itern,7))
        msc_crds(1:2) = mcon_tern(imsc,itern,1:2)*180.d0/pi
        msc_crds(3)   = mcon_tern(imsc,itern,3)
        call write_mascon_record(lumsc_out,tern_flag(ntern),ntern,nmsc_current,nint(mcon_prim(imsc,8)),msc_crds &
                              ,mcon_tern(imsc,itern,4),mcon_tern(imsc,itern,5),mcon_tern(imsc,itern,8),mcon_tern(imsc,itern,6) &
                              ,mcon_prim(imsc,11),mcon_prim(imsc,10),mcon_region(imsc),pointers,tern_colours(ntern,:))
      enddo
    endif
  enddo
  write(message,'(a,i6,a)')"Written out ",nmsc_current,' primary mascons that are outside the polygon'
  call status_update('STATUS','UTIL','grab_mascons',msc_file_out,message,0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                             !!!
!!         O U T P U T    M A S C O N    W I T H I N    T H E    P O L Y G O N                 !!!
!!                                                                                             !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  primary record
! calculate the centre of mass of the ternarys within the polygon
  mcon_prim(1,8) = ntern_inside    ! set the number of ternarys so that the CoM can be calculated
  call MASCON_calc_CoM(ntern_inside ,nvar_prim,nvar_tern,mcon_tern_inside(1,:,:),mcon_prim(1,:) )
  write(message,'(a,2f10.3)')'Polygon CoM (lat/lon): ',mcon_prim(1,1:2)*180.0/pi
  call status_update('STATUS','UTIL','grab_mascons',' ',message,0)

! fill in the mcon_prim information
  mcon_prim(1,3) = sum(mcon_tern_inside(1,:,3))/ntern_inside      ! average radius  of all the ternarys
  mcon_prim(1,4) = sum(mcon_tern_inside(1,:,4))                   ! sum of area     of all the ternarys
  mcon_prim(1,5) = sum(mcon_tern_inside(1,:,5))/ntern_inside      ! average hgt     of all the ternarys
  mcon_prim(1,6) = sum(mcon_tern_inside(1,:,6))/ntern_inside      ! average density of all the ternarys
  mcon_prim(1,7) = 1               ! # secondary mascons
  mcon_prim(1,8) = ntern_inside    ! # ternary mascons
  mcon_prim(1,9) = new_n_prim      ! # of first secondary mascon (= primary mascon number if there is only one secondary in each primary)
  mcon_prim(1,10) = 0              ! turn off the tidal amplitude estimates
  mcon_prim(1,11) = 1.           ! hardwire to 100% land for now
  mcon_prim(1,12) = sum(mcon_tern_inside(1,:,8))/ntern_inside      ! average density of all the ternarys
  nsec = 1
  msc_crds(1:2) = mcon_prim(1,1:2)*180.d0/pi
  msc_crds(3)   = mcon_prim(1,3)
  pointers = nmsc_current + 1
  colours(1:2) = 0
! define a primary mascon type, based on the density of the primary mascon
  if(region(1:4) == "Pine" )then
    prim_flags(1) = "PPineI"
    mcon_prim(1,10) = 100.d0             
  elseif(region(1:5) == "Thwai" )then
    prim_flags(1) = "PThwai"
    mcon_prim(1,10) = 100.d0             
  elseif(region(1:5) == "SthnO" )then
    prim_flags(1) = "PSthnO"
    mcon_prim(1,10) = 0.d0             
    colours(1:2) = 20
  elseif(region(1:5) == "Arcti" )then
    prim_flags(1) = "PArcti"
    mcon_prim(1,10) = 0.d0             
    colours(1:2) = 30
  elseif(region(1:5) == "Antar" )then
    prim_flags(1) = "PAntar"
    mcon_prim(1,10) = 0.d0             
    colours(1:2) = 12
  elseif(region(1:5) == "Ocean")then
    prim_flags(1) = "POcean"
    mcon_prim(1,10) = 0.d0             
    colours(1:2) = 50
  elseif(region(1:5) == "GrndO")then
    prim_flags(1) = "PGrndO"
    mcon_prim(1,10) = 0.d0             
    colours(1:2) = 70
  elseif(region(1:5) == "A400k")then
    prim_flags(1) = "PA400k"
    mcon_prim(1,10) = 0.d0             
    colours(1:2) = 80
  else
    prim_flags(1) = "P"//region(1:5)
    mcon_prim(1,10) = 100.d0                     
  endif
  call write_mascon_record(lumsc_out,prim_flags(1),nmsc_current+1,nsec,nint(mcon_prim(1,8)) &
                              ,msc_crds,mcon_prim(1,4),mcon_prim(1,5),mcon_prim(1,12),mcon_prim(1,6) &
                              ,mcon_prim(1,11)*100.d0 &
                              ,mcon_prim(1,10),region,pointers,colours)
!!   secondary records for this primary
  pointers(1) = nmsc_current + 1
  pointers(2) = 0
  nsec = nmsc_current + 1           
! define a secondary mascon type, based on the density of the primary mascon
  sec_flag = "S"//prim_flags(1)(2:6)
  call write_mascon_record(lumsc_out,sec_flag,nmsc_current+1,nmsc_current+1,nint(mcon_prim(1,8)) &
                              ,msc_crds,mcon_prim(1,4),mcon_prim(1,5),mcon_prim(1,12),mcon_prim(1,6) &
                              ,mcon_prim(1,11) &
                              ,mcon_prim(1,10),region,pointers,colours)
!     ternary records for this secondary
  do itern = 1,nint(mcon_prim(1,8))
    pointers(1) = nmsc_current + 1
    pointers(2) = nsec
    msc_crds(1:2) = mcon_tern_inside(1,itern,1:2)*180.d0/pi
    msc_crds(3)   = mcon_tern_inside(1,itern,3)
    tern_number = nint(mcon_tern_inside(1,itern,7))
    tern_flag(1) = "T"//prim_flags(1)(2:6)

! set different colours for different specific regions
    if(region(1:4) == "Pine")then
      colours(1:2) = 1750
    else if (region(1:5) == "Thwait")then
      colours(1:2) = 1725
    else if (region(1:5) == "SthnO")then
      colours(1:2) = 12
    else if (region(1:5) == "Arcti")then
      colours(1:2) = 12
    else if (region(1:5) == "Antar")then
      colours(1:2) = 1800
    else
      colours(1:2) = 1600
    endif

    call write_mascon_record(lumsc_out,tern_flag(1),tern_number,nmsc_current+1,nint(mcon_prim(1,8)) &
                              ,msc_crds &
                              ,mcon_tern_inside(1,itern,4),mcon_tern_inside(1,itern,5),mcon_tern_inside(1,itern,8)    &
                              ,mcon_tern_inside(1,itern,6),mcon_prim(1,11) &
                              ,mcon_prim(1,10),region,pointers,colours)
  enddo
  write(message,'(a,a)')"Written out single mascon within polygon for region: ",region
  call status_update('STATUS','UTIL','grab_mascons',msc_file_out,message,0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end




