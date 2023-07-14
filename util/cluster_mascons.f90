  program cluster_mascons

! program to take a set of primary mascons that define a particular region (e.g. a continent, an island or a drainage basin) and 
! reshape them to be roughly equal area and equal shape. This is done by reassigning the ternary mascons of the oversize primary
! mascons to the undersize primary mascons, until they are all roughly the same.
!
! Many of the subroutines and some of the logic replicates what is done in util/configure_mascons
!
! P. Tregoning
! 5 November 2014 as "reshape_mascons". 
! 11 February 2020 as "cluster_mascons"
!
! MODS:
! PT161108: changed to read the new mascon file format. The selection of mascons to reshape is now done within this program
!           rather than in the calling script.
! PT170123: add mascon codes to leave untouched to the command line
! PT170205: add the capability to enter a mascon number to reshape, rather than a region
! PT/HMcQ181123: fixed bug in numbering of secondary mascons
! 
! PT200211: replace the voronoi cell algorithm with a clustering algorithm that relies on knowing the distance from each ternary to the
!           1000 closest ternarys (as computed by program util/dist_flag)

  use mascon_mod     ! provides the declarations of the variables needed to read in and store all the mascons

  implicit none

! command line arguments
  character    :: mascon_file*150               ! mascon input file
  character    :: msc_file_out*150              ! mascon input file
  character    :: region*15                     ! region to adjust (all/Australasia/Europe/.... etc)
  real(kind=8) :: req_area                      ! required area of reshaped primary mascons (will dictate the size of the output mascons)
  integer*4    :: ioerr
  integer*4    :: n_untouched                        ! number of mascon codes to leave untouched by this program
  character*5,allocatable :: msc_codes_untouched(:)  ! codes to leave untouched by this program

! output primary and ternary mascon information
  integer*4                 :: nmsc_prim_out     ! number of primary mascons in the recomputed primary mascon configuration. 
  real(kind=8), allocatable :: msc_prim_out(:,:) ! number of adjusted primary mascon, lat/lon/hgt/area/density/tide_adjust?
  real(kind=8)              :: sum_area          ! total area of the input primary mascons within the requested region, calcuated from the input ternary mascons
  integer*4                 :: imsc

! parameters related to defining geographical regions
  real(kind=8), allocatable :: region_poly(:,:)  ! vertices of the polygon region
  integer*4                 :: npoints           ! number of polygon vertices
  real(kind=8)              :: rms_limit         ! rms lower limit of rms of centroid coord changes, below which we stop the iteration
  integer*4                 :: in_or_out         ! 1 if point in polygon, 0 if on the edge of the polygon, -1 if out of polygon
  integer*4                 :: nmsc_region       ! number of primary mascons found within the requested region
  integer*4, allocatable    :: region_prims(:)   ! list of primary mascon numbers found within the requested region
  integer*4, allocatable    :: region_prims_reshape(:)! list of primary mascon numbers after reshaping
  character*6, allocatable  :: prim_flags_reshape(:)! 6-char flags for type of mascon (PDeep/Pshelf etc) 
  real(kind=8)              :: region_area       ! sum of all ternary mascon areas found within the requested region
  logical                   :: cross_0E          ! logical to know whether to translate the region polygon so that it no longer crosses 0E
  logical                   :: region_density_good ! flag to say whether a primary mascon has the requisite density for a requested region
  logical                   :: leave_alone         ! flag to indicate whether this mascon should remain untouched

! parameters for counting changes in mascon numbers etc
  integer*4  :: new_n_prim        ! the new total number of primary mascons
  integer*4  :: new_n_sec         ! the same value as new_n_prim
  integer*4  :: new_tern_per_prim ! new maximum number of ternary mascons per primary
  integer*4  :: new_tern_per_sec  ! new maximum number of ternary mascons per secondary
  integer*4  :: n_reshaped        ! number of primary mascons after reshaping in the region
  integer*4  :: n_to_reshape      ! number of primary mascons within the region that will be reshaped
  integer*4  :: new_max_tern_per_prim ! maximum number of ternary mascons per primary (and, hence, per secondary) in the reshaped file
  integer*4  :: n_empty_prim      ! number of empty primary mascons after reshape
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
  character*6  :: sec_flag_tmp      ! temporary secondary mascon code based on primary mascon code
  integer*4    :: tide_flag         ! set to zero for now (ie no tidal amplitude modelling)

! parameters concerned with transferring mascons
  integer*4    :: i_big
  real(kind=8) :: tmp_big
  real(kind=8), allocatable :: dist_pr_2_ter(:,:) ! distances between primary CoM and all ternarys within the primary. Second col is pointer to the ternary #
  integer*4    :: itern

! PT170205: variables required to reshape just a single primary mascon
  integer*4    :: req_mascon        ! number of mascon to reshpae, or zero otherwise
 
! unit numbers
  integer*4, parameter      :: lumsc_in  = 10       ! unit number of  input mascon file
  integer*4, parameter      :: lumsc_out = 11       ! unit number of output mascon file
  integer*4, parameter      :: lupoly    = 12       ! unit number of file containing polygons around continents etc

! counters and logicals
  integer*4 :: i,j,k,counter
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

! PT170424: logical variable to indicate that we want to reshape the Shelf mascons
  logical       :: reshape_shelf

! RM190207: logical variable to indicate that primaries will be reshaped to equal ternary mascons
  logical       :: prim_to_tern

! RM190207: variables to hold ternary information for prim_to_tern
  real(kind=8),allocatable :: tern_crds(:,:) ! array of all the ternary mascons and attributes
  integer*4,   allocatable :: n_terns(:)     ! number of ternary mascons in each reshaped primary mascon   

! PT200211: variables related to the new clustering algorithm



  pi = 4.d0*atan(1.d0)
  new_max_tern_per_prim = 0
  n_empty_prim = 0
  reshape_shelf = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  get the command line and cmdfile information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call getarg(1,mascon_file)
  if(mascon_file(1:1) == " ")then
    call status_update('WARNING','UTIL','cluster_mascons',' ', &
                   'Runstring: cluster_mascons mascon_file mascon_file_out region max_area/tern rms_limit [n_untouched + msc codes]',0)
    call status_update('FATAL','UTIL','cluster_mascons',' ', &
                        'use "MC" as region (followed by mascon number e.g. MC0010) to reshape a single mascon',0)
  else
    call getarg(2,msc_file_out)
    region = " "
    call getarg(3,region)
    call getarg(4,arg)
      if(arg(1:4) == "tern")then
        prim_to_tern = .true.
        req_area = 2.25e8
      else
        prim_to_tern = .false.
        read(arg,*)req_area
      endif
    call getarg(5,arg)
    read(arg,*)rms_limit
! PT170205: check whether it is a region or a mascon number
    if(region(1:2) == "MC" )then
      read(region(3:7),*)req_mascon
      write(message,'(a,i4,a,e8.2,a)')'Will re-cluster primary mascon number "',req_mascon,'" to size ',req_area, " km^2"
! PT170424: check whether it is the continental shelf region
    else if (region(1:2)  == "SH")then
      read(region(3:6),*)req_mascon
      write(message,'(a,i4,a,e8.2,a)')'Will re-cluster shelf mascon number "',req_mascon,'" to size ',req_area, " km^2"
    else
      req_mascon = 0
      write(message,'(a,a,a,e8.2,a)')'Will re-cluster primary mascon file in region "',region,'" to size ',req_area, " km^2"
    endif
    call status_update('STATUS','UTIL','reshape_mascons',mascon_file,message,0)
    call getarg(6,arg)
    if(arg(1:1) /= " ")then
      read(arg,*)n_untouched   ! the number of mascon codes to be left alone by this program
      if(n_untouched > 0)then
        ! allocate the character array
        allocate(msc_codes_untouched(n_untouched))
        ! read them in from the command line
        message = " "
        msc_codes_untouched = "     "
        do i=1,n_untouched
          call getarg(6+i,msc_codes_untouched(i))
          write( message(52+(i-1)*6:51+(i-1)*6+5),'(a5)')msc_codes_untouched(i)
        enddo
      endif
      write(message(1:51),'(a,i5,a)')'Will not adjust ',n_untouched,' mascons with specific codes: '
      call status_update('STATUS','UTIL','reshape_mascons',' ',message,0)
    endif
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  open the input ternary mascon file, read the header, allocate arrays and read in the mascons
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call read_msc_hdr(lumsc_in,mascon_file,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                                ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
  call allocate_mascon_arrays
! PT180307: remove the trimlen in the name of the mascon file
  call read_mascon_file(lumsc_in,mascon_file)

! open the output file
  open(lumsc_out,file=msc_file_out,status='unknown')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! allocate some arrays
  ! PT170607: first, deallocate the ocean arrays, since we are short on memory and we don't use them
  !           in this program
  deallocate(mcon_ocean_prim)
  deallocate(mcon_ocean_tern)
  deallocate(mcon_ocean_tern_ptr)
  deallocate(ocean_eftid_part)
  deallocate(mcon_tide_name)
  deallocate(mcon_tide_period)
  deallocate(msc_tide)

  ! now allocate the primary mascons for the region
  allocate(region_prims(total_prim))    ! make it as big as it could possibly need to be.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if(req_mascon == 0 .and. region(1:5) /="SthnO")then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! define coordinates of a polygon that defines the requested region
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    open(lupoly,file='mascon.polygons',status='old')

! read the file to find the line for the region of interest
    line = " "
    npoints = -999
    cross_0E = .false.
    do while (line(9:19) /= region)
      read(lupoly,'(a)',iostat=ioerr,end=1000)line
      if(ioerr == 0 .and. line(9:19) == region) then
        read(line(1:6),*)npoints
        allocate(region_poly(npoints,2))
! PT170608: use the value from the command line instead
!        read(line(20:40),*)rms_limit
        do i=1,npoints
          read(lupoly,*)region_poly(i,:)
          if(region_poly(i,2) < 0.d0 )cross_0E = .true.
        enddo
        write(message,'(a,a,a,i6,a)')"Polygon around region ",region," has ",npoints," vertices"
        call status_update('STATUS','UTIL','cluster_mascons','mascon.polygons',message,0)

! PT161118: check whether there are any -ive longitudes in the polygon vertices. If so, the polygon
!           crosses 0E. In this case, shift it 90degrees to the east so that the in_polygon algorithm
!           will work.
        if(cross_0E)region_poly(:,2) = region_poly(:,2) + 90.d0

      endif
    enddo

1000 continue

    if(npoints == -999)then
      write(message,'(a,a,a)')"End of reading file and region ",region," not found"
      call status_update('FATAL','UTIL','cluster_mascons','mascon.polygons',message,0)
    endif

  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! loop over all primary mascons and identify the ones within the region of interest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  n_to_reshape  = 0
  ntern_reshape = 0
  region_area = 0.d0
  new_max_tern_per_prim = 0

  do imsc = 1,total_prim

! PT170123: determine whether the mascon code indicates that the mascon is one that we should not modify
    leave_alone = .false.
    do i=1,n_untouched
!print*,'imsc, prim_flags(imsc),i,msc_codes_untouched(i)',imsc, prim_flags(imsc),i,msc_codes_untouched(i)
      if(prim_flags(imsc)(2:6) == msc_codes_untouched(i))then
        leave_alone = .true.
        write(message,'(a,i6,1x,a)')'leave this mascon alone:',imsc,prim_flags(imsc)
        call status_update('STATUS','UTIL','cluster_mascons','',message,0)
      endif
    enddo
! determine whether the mascon density is appropriate for the requested region
    region_density_good = .false.
    if(region(1:5) == "Pacif" .or. region(1:5) == "Ross_" .or.region(1:5) == "Amund" .or. region(1:6) == "Arctic" &
             .or.prim_flags(imsc) == "PShelf" .or. region(1:5) == "SthnO".or. region(1:5) == "N_Atl" .or. region(1:5) == "Atlan" &
             .or. region(1:5) == "S_Atl" .or. region(1:6) == "Indian")then
      if (nint(mcon_prim(imsc,6)) > 1000)then
        region_density_good = .true.
      endif
    else
      if (nint(mcon_prim(imsc,6))  == 1000)then
        region_density_good = .true.
      endif
    endif

! PT170205: force some settings if we just want to reshape a particular mascon
    if(req_mascon /= 0 .and. imsc /= req_mascon)then
      leave_alone = .true.
    endif
 
! is this mascon in the region of interest? Does it contain land?
    if(prim_flags(imsc) /= "PShelf" .and. region(1:5) /= "Antar" .and. region(1:6) /= "Arctic" .and. region(1:5) /= "SthnO" &
                       .and. region(1:5) /= "N_Atl" .and. region(1:6) /= "S_Atl" .and. region(1:5) /= "Atlan"  &
                       .and. region(1:6) /= "Indian"  &
                       .and. region_density_good .and. .not. leave_alone )then   
      call msc_in_region(cross_0E,mcon_prim(imsc,1)*180.d0/pi,mcon_prim(imsc,2)*180.d0/pi,region_poly(:,1) &
         ,region_poly(:,2),npoints,in_or_out)
      if(in_or_out == 1)then
        n_to_reshape  = n_to_reshape  + 1
        region_prims(n_to_reshape ) = imsc
        ntern_reshape = ntern_reshape + nint(mcon_prim(imsc,8))
        region_area = region_area + sum(mcon_tern(imsc,1:nint(mcon_prim(imsc,8)),4))
        prim_flags(imsc) = "P-----"
      else
! mascon is not in the region. see whether the number of ternarys per primary is greater than the current greatest
        if(nint(mcon_prim(imsc,8)) > new_max_tern_per_prim)new_max_tern_per_prim = mcon_prim(imsc,8)
      endif
!PT170607: need to break the "ocean" into north and south of the equator - otherwise we hit memory errors
    else if (prim_flags(imsc) == "PDeep " .and. (region(1:5) == "Ocean" .or.region(1:6) == "Indian") &
            .and..not. leave_alone ) then       ! all the ocean ones
        n_to_reshape  = n_to_reshape  + 1
        region_prims(n_to_reshape ) = imsc
        ntern_reshape = ntern_reshape + nint(mcon_prim(imsc,8))
        region_area = region_area + sum(mcon_tern(imsc,1:nint(mcon_prim(imsc,8)),4))
        prim_flags(imsc) = "P-----"
    else if (region(1:5)=="SthnO".and.mcon_prim(imsc,1) < -1.d0*pi/3.0 .and. mcon_prim(imsc,6) > 1010.d0  &
                 .and. .not. leave_alone .and. prim_flags(imsc) /= "PA400k"  )then      ! all the Southern Ocean ones
        n_to_reshape  = n_to_reshape  + 1
        region_prims(n_to_reshape ) = imsc
        ntern_reshape = ntern_reshape + nint(mcon_prim(imsc,8))
        region_area = region_area + sum(mcon_tern(imsc,1:nint(mcon_prim(imsc,8)),4))
        prim_flags(imsc) = "P-----"
    else if (prim_flags(imsc) ==  "PShelf" .and. region(1:2) == "SH" .and. .not. leave_alone ) then       ! all the shelf ones
        n_to_reshape  = n_to_reshape  + 1
        region_prims(n_to_reshape ) = imsc
        ntern_reshape = ntern_reshape + nint(mcon_prim(imsc,8))
        region_area = region_area + sum(mcon_tern(imsc,1:nint(mcon_prim(imsc,8)),4))
        prim_flags(imsc) = "P-----"
        reshape_shelf = .true.
    else if (region(1:5) == "Antar" .and. mcon_prim(imsc,1) < -pi/3.d0 .and. prim_flags(imsc) == "PLand " &
               .and. .not. leave_alone )then  ! all mascons south of S59
        n_to_reshape  = n_to_reshape  + 1
        region_prims(n_to_reshape ) = imsc
        ntern_reshape = ntern_reshape + nint(mcon_prim(imsc,8))
        region_area = region_area + sum(mcon_tern(imsc,1:nint(mcon_prim(imsc,8)),4))
        prim_flags(imsc) = "P-----"
! PT180320: decrease the latitude down to below Greenland (just under 60N)
    else if (region(1:6) == "Arctic" & 
                 .and. (mcon_prim(imsc,1) > 65.d0*pi/180.d0 .or. (mcon_prim(imsc,1) > 58.d0*pi/180.d0 &
                         .and. mcon_prim(imsc,2) > 295.d0*pi/180.d0 .and. mcon_prim(imsc,2) < 360.d0*pi/180.d0 ) ) &
                 .and. prim_flags(imsc) /= "PLand " .and. prim_flags(imsc) /= "PGrndO")then  ! all ocean mascons north of N65 and north of 58N around Greenland
        n_to_reshape  = n_to_reshape  + 1
        region_prims(n_to_reshape ) = imsc
        ntern_reshape = ntern_reshape + nint(mcon_prim(imsc,8))
        region_area = region_area + sum(mcon_tern(imsc,1:nint(mcon_prim(imsc,8)),4))
        prim_flags(imsc) = "P-----"
    else 
! check the max_tern_per_prim
      if(nint(mcon_prim(imsc,8)) > new_max_tern_per_prim)new_max_tern_per_prim = mcon_prim(imsc,8)
    endif
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT170205: add some code to set up the above values for the case when we want to reshape one 
!           mascon in particular
  if(req_mascon /= 0)then
    n_to_reshape = 1
    region_prims(1) = req_mascon
    ntern_reshape = nint(mcon_prim(req_mascon,8))
    region_area = sum(mcon_tern(req_mascon,1:nint(mcon_prim(req_mascon,8)),4))
    prim_flags(req_mascon) = "P-----"
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RM190212: from the area of the region, work out how many new primary mascons we need and allocate arrays
  if(prim_to_tern)then
    n_reshaped = ntern_reshape
  else
    n_reshaped = nint(region_area/req_area)
  endif
! PT180507: set this to one if it is zero but a specific mascon was requested
  if(n_reshaped == 0 .and. req_mascon /= 0)n_reshaped = 1

  if(n_reshaped < 1) then
    write(message,'(a)')"No ternary mascons found inside requested area. Nothing to do. Copy input file to output"
    call status_update('STATUS','UTIL','cluster_mascons',' ',message,0)
    write(message,'(a,a,a,a)')"cp ",mascon_file(1:trimlen(mascon_file))," ",msc_file_out(1:trimlen(msc_file_out))
    call system(message)
    stop
  else
    write(message,'(a,e10.4,a,i6,a,a,i6,a)')"Region area (",region_area/1.d3," km^2) will be broken into",n_reshaped &
        ,' primary mascons.'," There were originally",n_to_reshape,' primary mascons in the region.'
    call status_update('STATUS','UTIL','cluster_mascons',' ',message,0)
  endif
! allocate the arrays
  allocate(msc_reshaped_prims(n_reshaped,nvar_prim))
! PT161123: can we reduce the second dimension to ntern_reshape/4 so that we won't exceed memorylimits?
  allocate(msc_tern_reshaped(n_reshaped,nint(ntern_reshape/1.0),nvar_tern))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! RM190204: allocate arrays using no of reshaped mascons (causing loss of ternaries?)
  allocate(region_prims_reshape(n_reshaped))
  allocate(prim_flags_reshape(n_reshaped))

  do i=1,n_reshaped
    if(i <= total_prim) then
      region_prims_reshape(i) = region_prims(i)
      prim_flags_reshape(i) = prim_flags(i)
    else
      region_prims_reshape(i) = 0
      prim_flags_reshape(i) = ""
    endif
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! for the primary mascons within the requested region, do something with them .....
  debug = .false.
  debug_rms = .true.

! RM190207: determine if we want to make all primaries equal to existing ternaries or use the vornoi subroutine
  if(prim_to_tern)then
    write(message,'(a)')"Reshaping primaries to mimic each ternary in region"
    call status_update('STATUS','UTIL','cluster_mascons',' ',message,0)
! RM190212: allocate arrays    
    allocate(tern_crds(ntern_reshape,nvar_tern))
    allocate(n_terns(n_reshaped))
! PT161118: shift all longitudes 90 deg east if the mascons span the 0 deg longitude line    
    if(cross_0E)then
      write(message,'(a)')"Shifting longitudes by 90 deg for region "
      call status_update('STATUS','UTIL','cluster_mascons',' ',message,0)
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RM190212: transfer all the ternary mascons into one array
    counter = 0
    do i=1,n_to_reshape
      do j=1,nint(mcon_prim(region_prims_reshape(i),8))
        counter = counter+1
        tern_crds(counter,1:8) = mcon_tern(region_prims_reshape(i),j,1:8)
! PT161118: shift all longitudes 90 deg east if the mascons span the 0 deg longitude line
        if(cross_0E)then
          tern_crds(counter,2) = tern_crds(counter,2) + pi/2.d0
          if(tern_crds(counter,2) > 2.d0*pi)tern_crds(counter,2) = tern_crds(counter,2) - 2.d0*pi
        endif
      enddo
    enddo
! RM190212: copy ternary parameters (lat,lon,rad,area,height,density and geoid) to new primaries
    n_terns=0
    do imsc=1,n_reshaped
      msc_reshaped_prims(imsc,1:6) = tern_crds(imsc,1:6)
      msc_reshaped_prims(imsc,12) = tern_crds(imsc,8)
      ! assign ternary to its primary equivalent
      n_terns(imsc) = n_terns(imsc) + 1
      tern_crds(imsc,9) = imsc
      do k=1,8
        msc_tern_reshaped(imsc,n_terns(imsc),k) = tern_crds(imsc,k)
      enddo
! RM190212: determine number of ternaries in each primary (should = 1), calculate %land from height
      msc_reshaped_prims(imsc,8) = n_terns(imsc)
      if(msc_reshaped_prims(imsc,5) > 0.d0)then
        msc_reshaped_prims(imsc,11) = 1.d0
      endif
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT161118: shift all longitudes 90 deg west if the mascons span the 0 deg longitude line
    if(cross_0E)then
      write(message,'(a)')"Shifting longitudes back by 90 deg for region "
      call status_update('STATUS','UTIL','cluster_mascons',' ',message,0)
! first the primary mascons
      do i=1,n_reshaped
        msc_reshaped_prims(i,2) = msc_reshaped_prims(i,2) - pi/2.d0
        if(msc_reshaped_prims(i,2) < 0.d0)msc_reshaped_prims(i,2) = msc_reshaped_prims(i,2) + 2.d0*pi
      enddo
! now the ternary mascons
      do i=1,n_reshaped
        msc_tern_reshaped(i,:,2) = msc_tern_reshaped(i,:,2) - pi/2.d0
        do j=1,nint(msc_reshaped_prims(i,8))
          if(msc_tern_reshaped(i,j,2) < 0.d0)msc_tern_reshaped(i,j,2) = msc_tern_reshaped(i,j,2) +2.d0*pi
        enddo
      enddo
    endif

! PT200211: create the new clustered primaries with specified area
  else
    call cluster_primaries(.false.,region,rms_limit,cross_0E,n_to_reshape,region_prims_reshape,msc_reshaped_prims,n_reshaped &
                        ,msc_tern_reshaped,nint(ntern_reshape/1.0),debug,debug_rms)
  endif

  call srand(n_reshaped)
  if(region(1:3) == "Oce" .or.region(1:5) == "Pacif"  .or. region(1:5) == "Amund" .or.region(1:5) == "Ross_"  &
                          .or.region(1:6) == "Arctic" .or. region(1:5) == "SthnO" .or. region(1:5) == "N_Atl" &
                          .or.region(1:5) == "S_Atl"  .or. region(1:5) == "Atlan" .or. region(1:6) == "Indian")then
    col_range(1) = 0.
    col_range(2) = 450.
  else if (region(1:5) == "Shelf" .or. region(1:2) == "SH" )then
    col_range(1) = 500.
    col_range(2) = 750.
  else
    col_range(1) = 750.
    col_range(2) = 1500.
  endif

  do imsc = 1,n_reshaped
! assign a colour to each ternary in a primary
! PT190320: if we reshape a single mascon that is ocean then we need to identify it by density and assign the colour
! RM190527: added option to colour GIA mascons
    if(region(1:2) == "MC" .and. msc_reshaped_prims(imsc,6) > 3000.d0)then
      prim_colour = nint(500. + rand(0) * (750.-500.))
    else if(region(1:2) == "MC" .and. msc_reshaped_prims(imsc,6) > 1010.d0)then
      prim_colour = nint(0. + rand(0) * (450.-0.))
    else
      prim_colour = nint(col_range(1) + rand(0) * (col_range(2)-col_range(1)))
    endif

    do itern=1,nint(msc_reshaped_prims(imsc,8))
!print*,'imsc,itern,msc_tern_reshaped(imsc,itern,7)',imsc,itern,msc_tern_reshaped(imsc,itern,7),msc_reshaped_prims(imsc,8) &
!        ,msc_reshaped_prims(imsc,1:2)*180./pi
      tern_number = nint(msc_tern_reshaped(imsc,itern,7))
!      if(tern_number > 0)then
        tern_colours(tern_number,:) = prim_colour
!      else
! PT161119: put a BAD fix in here. I don't know why msc_reshaped_prims(imsc,8) is so big in some cases, but it 
!           is causing tern_numbers to look beyond where ternarys have been defined in the array. Set it back to 
!           a reasonable number
!        msc_reshaped_prims(imsc,8) = itern - 1
!      endif

    enddo

! assign the new region name to each primary/secondary/ternary
!    mcon_region(imsc) = region

! by the way, we need to assign a mascon number to any extra primary mascons that were created
    if(imsc > n_to_reshape)region_prims_reshape(imsc) = imsc

! and, finally, check whether the new_max_tern_per_prim needs updating
    if(nint(msc_reshaped_prims(imsc,8)) > new_max_tern_per_prim)new_max_tern_per_prim = nint(msc_reshaped_prims(imsc,8))

! RM190212: check if any of the reshaped primaries contain < 1 ternary
    if(nint(msc_reshaped_prims(imsc,8)) < 1)n_empty_prim = n_empty_prim + 1
  enddo

! RM190212: warn user if any empty primary mascons exist
  if(n_empty_prim > 0)then
    write(message,'(a,i2,a)')"Found ",n_empty_prim," empty reshaped primary mascon(s)"
    call status_update('WARNING','UTIL','cluster_mascons',' ',message,0)
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                          O U T P U T    A L L    M A S C O N S
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write out the first header line, containing the code and the numbers of mascons etc
! PT190319: subtract off the number of empty reshaped primary mascons. We don't write these ones out.
  new_n_prim = total_prim - n_to_reshape + n_reshaped - n_empty_prim
  new_n_sec  = new_n_prim

!!!!!!!!!!!!!!!!!!!!!
! write the header !!
!!!!!!!!!!!!!!!!!!!!! 
! add a new line to the header to represent what has been done in running this program
  call date_and_time(date,time,timezone)
  write(msc_hdr_lines(n_hdr_lines+1),'(a,a,a,i6,a,e12.6,a,1x,a8,1x,a10,1x,a5)') &
       "# cluster_mascons, region: ",region," has ",n_reshaped &
       ," primary mascons of size ~",req_area/1.d0," m^2. Timetag: ",date,time,timezone
! generate new unique code for the output program
  call hash_header(msc_hdr_lines(2:n_hdr_lines+1),n_hdr_lines+1,hashcode)

! write the new code and mascon numbers to the first line of the header
  write(msc_hdr_lines(1)(1:81),'(a,i6,5i8)')hashcode(1:8),new_n_prim,new_n_sec,max_tern,max_sec_per_prim &
                                           ,new_max_tern_per_prim,new_max_tern_per_prim

! and write out all header lines to the output file
  do i=1,n_hdr_lines+1
    write(lumsc_out,'(a)')msc_hdr_lines(i)
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write out all the mascons. 
! First, the ones untouched
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(message,'(a,i6,a)')"Writing out ",total_prim-n_to_reshape," untouched primary mascons"
  call status_update('STATUS','UTIL','cluster_mascons',msc_file_out,message,0)
  nmsc_current = 0
  do imsc = 1,total_prim 
    if(prim_flags(imsc) /= "P-----")then
      nmsc_current = nmsc_current + 1
!     primary record
      nsec = 1
      msc_crds(1:2) = mcon_prim(imsc,1:2)*180.d0/pi
      msc_crds(3)   = mcon_prim(imsc,3)
      pointers = nmsc_current
      colours(1:2) = 0
      tide_flag    = 0  ! PT170606: for now, don't turn on tidal amplitude modelling
      call write_mascon_record(lumsc_out,prim_flags(imsc),nmsc_current,nsec,nint(mcon_prim(imsc,8)),msc_crds     &
                              ,mcon_prim(imsc,4),mcon_prim(imsc,5),mcon_prim(imsc,12),mcon_prim(imsc,6) &
                              ,mcon_prim(imsc,11),mcon_prim(imsc,10),mcon_region(imsc),pointers,colours)
!     secondary record for this primary. We still can only have one secondary per primary at this point.
      pointers(1) = nmsc_current
      pointers(2) = 0
! PT161121: we don't want this, we want the secondary mascons to be sequential along with the primary mascons
!      nsec = mcon_sec(imsc,1,7)           ! this is the unique secondary mascon number as read from the input mascon file
!      colours = sec_colour(nsec)
! PT170426: change secondary flag to be exactly the primary flag 
      sec_flag_tmp = "S"//prim_flags(imsc)(2:6)
      call write_mascon_record(lumsc_out,sec_flag_tmp,nmsc_current,nsec,nint(mcon_prim(imsc,8)),msc_crds &
                              ,mcon_prim(imsc,4),mcon_prim(imsc,5),mcon_prim(imsc,12),mcon_prim(imsc,6) &
                              ,mcon_prim(imsc,11),mcon_prim(imsc,10),mcon_region(imsc),pointers,colours)
!     ternary records for this secondary
      do itern = 1,nint(mcon_sec(imsc,1,8))
        pointers(1) = nmsc_current
        pointers(2) = nmsc_current
        ntern = nint(mcon_tern(imsc,itern,7))
        msc_crds(1:2) = mcon_tern(imsc,itern,1:2)*180.d0/pi
        msc_crds(3)   = mcon_tern(imsc,itern,3)
        call write_mascon_record(lumsc_out,tern_flag(ntern),ntern,nsec,nint(mcon_prim(imsc,8)),msc_crds &
                              ,mcon_tern(imsc,itern,4),mcon_tern(imsc,itern,5),mcon_tern(imsc,itern,8),mcon_tern(imsc,itern,6) &
                              ,mcon_prim(imsc,11),mcon_prim(imsc,10),mcon_region(imsc), pointers,tern_colours(ntern,:))
      enddo
    endif
  enddo
  write(message,'(a,i6,a)')"Written out",total_prim-n_to_reshape,' untouched mascons to output'
  call status_update('STATUS','UTIL','cluster_mascons',msc_file_out,message,0)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! now the reshaped ones for the region
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(message,'(a,i6,a)')"Writing out ",n_reshaped,"  reshaped primary mascons"
  call status_update('STATUS','UTIL','cluster_mascons',msc_file_out,message,0)
  do imsc = 1,n_reshaped

! PT190319: only write out the new reshaped primaries if there are actually ternary mascons contained within them
    if(msc_tern_reshaped(imsc,1,7) > 0.d0)then
      nmsc_current = nmsc_current + 1

!print*,'imsc,msc_reshaped_prims(i,:)',imsc,msc_reshaped_prims(imsc,:)
!!   primary record
      nsec = 1
      msc_crds(1:2) = msc_reshaped_prims(imsc,1:2)*180.d0/pi
      msc_crds(3)   = msc_reshaped_prims(imsc,3)
      pointers = nmsc_current 
      colours(1:2) = 0
! define a primary mascon type, based on the density of the primary mascon
      if(msc_reshaped_prims(imsc,6) == 1000.d0)then
        if(region(1:4) == "Pine" )then
          prim_flags_reshape(imsc) = "PPineI"
        else if (region(1:5) == "Thwai")then
          prim_flags_reshape(imsc) = "PThwai"
        else if (region(1:5) == "Third")then
          prim_flags_reshape(imsc) = "PThird"
!! PT190308: try using the ternary flag to name the primary flag for the new region
        else if (tern_flag(nint(msc_tern_reshaped(imsc,1,7)))(1:5) /= "TDeep" .and. &
                 tern_flag(nint(msc_tern_reshaped(imsc,1,7)))(1:5) /= "TLand" )then
          prim_flags_reshape(imsc)(1:6) = "P"//tern_flag(nint(msc_tern_reshaped(imsc,1,7)))(2:6)
        else
          prim_flags_reshape(imsc) = "PLand "
        endif
        msc_reshaped_prims(imsc,10) = 0             
      else
        if (reshape_shelf .or. region(1:5) == "Shelf") then
          prim_flags_reshape(imsc) = "PShelf"
          msc_reshaped_prims(imsc,10) = 0             ! for now, don't turn on the tidal amplitude estimates  
          region = "Shelf"
        else if (region(1:5) == "Arcti")then
          prim_flags_reshape(imsc) = "PArcti"
          msc_reshaped_prims(imsc,10) = 0             ! for now, don't turn on the tidal amplitude estimates  
          region = "Arctic"
        else if (region(1:5) == "SthnO")then
          prim_flags_reshape(imsc) = "PSthnO"
          msc_reshaped_prims(imsc,10) = 0             ! for now, don't turn on the tidal amplitude estimates  
! PT190309: allow Casp, Amazo, Murra etc to retain their codes
        else if (tern_flag(nint(msc_tern_reshaped(imsc,1,7)))(1:5) /= "TDeep" .and. &
                 tern_flag(nint(msc_tern_reshaped(imsc,1,7)))(1:5) /= "TLand" )then
          prim_flags_reshape(imsc)(2:6) = tern_flag(nint(msc_tern_reshaped(imsc,1,7)))(2:6)
          msc_reshaped_prims(imsc,10) = 0             ! for now, don't turn on the tidal amplitude estimates  
!          region = "Southern_ocean "
        else if (region(1:5) /= "Shelf" )then
          prim_flags_reshape(imsc) = "PDeep "
          msc_reshaped_prims(imsc,10) = 0  
        endif
      endif
      call write_mascon_record(lumsc_out,prim_flags_reshape(imsc),nmsc_current,nsec & 
                              ,nint(msc_reshaped_prims(imsc,8)) &
                              ,msc_crds,msc_reshaped_prims(imsc,4),msc_reshaped_prims(imsc,5),msc_reshaped_prims(imsc,12) &
                              ,msc_reshaped_prims(imsc,6),msc_reshaped_prims(imsc,11)*100.d0 &
                              ,msc_reshaped_prims(imsc,10),region,pointers,colours)
!!   secondary records for this primary
      pointers(1) = nmsc_current
      pointers(2) = 0
      nsec = nmsc_current           
      colours = 0
! define a secondary mascon type, based on the density of the primary mascon
      sec_flag = "S"//prim_flags_reshape(imsc)(2:6)
!      if(msc_reshaped_prims(imsc,6) == 1000.d0)then
!        sec_flag = "SLand "
!      else
!        sec_flag = "SDeep "
!      endif
      call write_mascon_record(lumsc_out,sec_flag,nmsc_current,nsec,nint(msc_reshaped_prims(imsc,8)) &
                              ,msc_crds,msc_reshaped_prims(imsc,4),msc_reshaped_prims(imsc,5),msc_reshaped_prims(imsc,12) &
                              ,msc_reshaped_prims(imsc,6),msc_reshaped_prims(imsc,11)*100.d0 &
                              ,msc_reshaped_prims(imsc,10),region,pointers,colours)
!     ternary records for this secondary
      do itern = 1,nint(msc_reshaped_prims(imsc,8))
        pointers(1) = nmsc_current 
        pointers(2) = nmsc_current
        msc_crds(1:2) = msc_tern_reshaped(imsc,itern,1:2)*180.d0/pi
        msc_crds(3)   = msc_tern_reshaped(imsc,itern,3)
        tern_number = nint(msc_tern_reshaped(imsc,itern,7))
! PT170203: ensure that the correct ternary code is written out for the Shelf, PineI and Thwai regions
        if(region(1:4) == "Pine" )then
          tern_flag(tern_number) = "TPineI"
        else if (region(1:5) == "Thwai")then
          tern_flag(tern_number) = "TThwai"
        else if (region(1:5) == "Third")then
          tern_flag(tern_number) = "TThird"
        else if (region(1:5) == "Arcti")then
          tern_flag(tern_number) = "TArcti"
        else if (region(1:5) == "SthnO")then
          tern_flag(tern_number) = "TSthnO"
        else if (reshape_shelf)then
          tern_flag(tern_number) = "TShelf"
        endif
          call write_mascon_record(lumsc_out,tern_flag(tern_number),tern_number,nsec &
                              ,nint(msc_reshaped_prims(imsc,8)),msc_crds &
                              ,msc_tern_reshaped(imsc,itern,4),msc_tern_reshaped(imsc,itern,5),msc_tern_reshaped(imsc,itern,8)   &
                              ,msc_tern_reshaped(imsc,itern,6),msc_reshaped_prims(imsc,11) &
                              ,msc_reshaped_prims(imsc,10),region,pointers,tern_colours(tern_number,:))
      enddo
    endif
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end




