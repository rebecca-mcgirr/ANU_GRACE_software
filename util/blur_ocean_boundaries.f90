  program blur_ocean_boundaries

! program to identify all the primary mascons that have a CoM near the edges of ocean segments created defined by the user and 
! merge them so that they can be reshped and region boundaries become less obvious in the mascon geometries
!
! This is a modified version of grab_mascons that takes a start point and end point to define a line along a longitude or
! latitude and only grabs ocean mascons
!
! R. McGirr
! 6 April 2020

  use mascon_mod     ! provides the declarations of the variables needed to read in and store all the mascons

  implicit none

! command line arguments
  character    :: mascon_file*150               ! mascon input file
  character    :: msc_file_out*150              ! mascon output file
  real(kind=8) :: dist
  real(kind=8) :: lon1, lat1, lon2, lat2
  integer*4    :: imsc,itern
  character    :: region*15                     ! region to re-define (PineIs, Thwaites, Totten, Helheim etc)
! parameters related to defining the geographical region
  integer*4                 :: in_or_out         ! 1 if point in polygon, 0 if on the edge of the polygon, -1 if out of polygon
  integer*4, allocatable    :: region_prims(:)   ! list of primary mascon numbers found within the requested region

! variables to hold the reshaped mascons
  integer*4                 :: ntern_reshape                ! number of ternary mascons to be reshaped 

! variables related to outputting mascons using write_mascon_record
  real(kind=8) :: msc_crds(3)       ! temp storage of mascon lat/lon/rad
  integer*4    :: colours(2)        ! colour of secondary and primary mascon
  integer*4    :: pointers(2)       ! tern2prim and tern2sec pointers
  integer*4    :: nmsc_current      ! (re)counting of the primary mascons as we output them
  integer*4    :: nsec,ntern        ! temporary value for the secondary mascon number and the ternary mascon number
  character*6  :: sec_flag          ! single, temporary value assigned to the type of secondary mascon
  integer*4    :: tern_number       ! temporary value of the actual ternary mascon number
  character*12 :: hashcode          ! temporary variable used to calculate the new character code for the output mascon file
  real(kind=8) :: col_range(2)      ! variable used to calculate new colour for primary mascons
  integer*4    :: prim_colour       ! new colour for primary mascon after reshaping

! unit numbers
  integer*4, parameter      :: lumsc_in  = 10       ! unit number of  input mascon file
  integer*4, parameter      :: lumsc_out = 11       ! unit number of output mascon file

! counters and logicals
  integer*4 :: i,j,k

! other stuff
  character     :: line*200,message*250,arg*100
  real(kind=8)  :: pi

! variables for writing date and time into header of output file
  character(8)  :: date
  character(10) :: time
  character(5)  :: timezone

! PT161218: variables to densify the points along the polygon that defines the region
  integer*4     :: region_prims_dense(100000)    ! array to hold the primary mascon numbers that reside within the dense polygon

! variables to grab the ternarys from existing masons and move them
  logical, allocatable :: remove_tern(:,:)
  integer*4            :: all_in
  integer*4            :: n_to_reshape
  integer*4            :: n_other_density               ! number of ternarys in an original mascon that don't meet criteria to be grabbed
  integer*4            :: new_n_prim,new_n_sec          ! new numbers of primary and secondary mascons ???
  real(kind=8),allocatable :: mcon_tern_inside(:,:,:)   ! array to hold the ternary mascons that are found to be inside the polygon
  real(kind=8),allocatable :: mcon_tern_tmp(:,:,:)
  character*6          :: sec_flag_tmp                  ! temporary variable to construct the secondary mascon flag based on the primary mascon

! variables to copy some of the original mcon_prim information
  real(kind=8), allocatable :: mcon_prim_orig(:)        ! keep a copy of the original numbers of ternary mascons in each primary

! some values
  logical      :: found
  real(kind=8) :: req_area                              ! required area for a primary mascon ???
  real(kind=8) :: info_flag
  real(kind=8) :: tern_spacing                          ! spacing in latitude (in decimal degrees) of ternary mascons
  real(kind=8) :: lon, minlon, maxlon
  real(kind=8) :: lat, minlat, maxlat
  integer*4    :: n_tern_tmp
  integer*4    :: ntern_inside                          ! total number of ternary mascons found within the polygon region
  integer*4    :: ntern_inside_tmp
  logical      :: only_prim_nums                        ! logical to output just the primary mascon numbers contained within the region
  logical      :: line_of_lon, line_of_lat

  pi = 4.d0*atan(1.d0)
  tern_spacing = 10.d0/60.d0      ! we use 10' ternary mascons
  line_of_lon = .false.
  line_of_lat = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  get the command line and cmdfile information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call getarg(1,mascon_file)
  if(mascon_file(1:1) == " ")then
    call status_update('FATAL','UTIL','blur_ocean_boundaries',' ' &
              ,'Runstring: blur_ocean_boundaries mascon_file_in mascon_file_out dist_degrees lat1 lon1 lat2 lon2 info_flag',0)
  else
    call getarg(2,msc_file_out)
    call getarg(3,arg)
    read(arg,*)dist
    call getarg(4,arg)
    read(arg,*)lat1
    call getarg(5,arg)
    read(arg,*)lon1
    call getarg(6,arg)
    read(arg,*)lat2
    call getarg(7,arg)
    read(arg,*)lon2
    call getarg(8,arg)
    read(arg,*)info_flag

! PT190822: use a negative info_flag as a flag to output only the primary mascons contained within the region
    if(info_flag < 0.d0)then
      only_prim_nums = .true.
      write(message,'(a,f6.4,a)')'Will output only primary mascon numbers within ',dist,' degrees of ocean segment boundaries'
      call status_update('STATUS','UTIL','blur_ocean_boundaries',mascon_file,message,0)
    else
      only_prim_nums = .false.
      write(message,'(a,f6.4,a)')'Will create primary mascon file for mascons within ',dist,' degrees of ocean segment boundaries'
      call status_update('STATUS','UTIL','blur_ocean_boundaries',mascon_file,message,0)
    endif

    if(lat1 == lat2 .and. lon1 /= lon2)then
      line_of_lat = .true.
      lat = lat1
      minlon = min(lon1, lon2)
      maxlon = max(lon1, lon2)
      write(message,'(a,f8.4,a,f8.4,a,f8.4,a)')'Ocean segment defined along ',lat,' degrees latitude from ',minlon,' to ',maxlon,' degrees longitude'
      call status_update('STATUS','UTIL','blur_ocean_boundaries',mascon_file,message,0)
    else if(lon1 == lon2 .and. lat1 /= lat2)then
      line_of_lon = .true.
      lon = lon1
      minlat = min(lat1, lat2)
      maxlat = max(lat1, lat2)
      write(message,'(a,f8.4,a,f8.4,a,f8.4,a)')'Ocean segment defined along ',lon,' degrees longitude from ',minlat,' to ',maxlat,' degrees latitude'
      call status_update('STATUS','UTIL','blur_ocean_boundaries',mascon_file,message,0)
    else if(lat1 > lat2 .or. lon1 > lon2)then
      call status_update('FATAL','UTIL','blur_ocean_boundaries',' ' &
              ,'The first pair of coords must be less than the second pair (i.e. 0 60 0 120 or -60 0 60 0)',0)
    else
      call status_update('FATAL','UTIL','blur_ocean_boundaries',' ' &
              ,'Can only define a straight segment along a single line of longitude/latitude',0)
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
! loop over all primary mascons and identify the ones within the region of interest that might not
! cross the polygon line (ie are contained wholly within the region
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  region = "OceanBlur"
  n_to_reshape = 0
  ntern_reshape = 0
  all_in = 0
  do imsc = 1,total_prim
    if( line_of_lat .and. (mcon_prim(imsc,1)*180.d0/pi >= lat-dist .and. mcon_prim(imsc,1)*180.d0/pi <= lat+dist &
                 .and. mcon_prim(imsc,2)*180.d0/pi >= minlon .and. mcon_prim(imsc,2)*180.d0/pi <= maxlon) &
                 .and. mcon_prim(imsc,6) > 1010.d0 )then
      in_or_out = 1
    else if( line_of_lon .and. (mcon_prim(imsc,2)*180.d0/pi >= lon-dist .and.  mcon_prim(imsc,2)*180.d0/pi <= lon+dist &
                 .and. mcon_prim(imsc,1)*180.d0/pi >= minlat .and. mcon_prim(imsc,1)*180.d0/pi <= maxlat) &
                 .and. mcon_prim(imsc,6) > 1010.d0 )then
      in_or_out = 1
    else
      in_or_out = 0
    endif
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
      endif
    endif
  enddo

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
      if(mcon_tern(imsc,itern,6) > 1010.)then
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
    enddo
    if(nint(mcon_prim(imsc,8))-n_tern_tmp > 0)then
      write(message,'(a,i7,a,i7,2a,i7,a)')'Primary mascon',imsc,' has ',n_tern_tmp,' ternary mascons outside the region' &
          ,' (or wrong density) and ',nint(mcon_prim(imsc,8))-n_tern_tmp,' ternary mascons within the region with correct density'
      call status_update('STATUS','UTIL','blur_ocean_boundaries',' ',message,0)
    endif

! PT170119: update the number of ternarys remaining in the external primary
    mcon_prim(imsc,8) = n_tern_tmp
    mcon_sec(imsc,1,8) = mcon_prim(imsc,8)

! RM200114: count the original primary mascons that are either not required density or partly outside the region
    if( n_tern_tmp > 0)n_other_density = n_other_density + 1

  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                                          !! 
!!     E S T A B L I S H     A     S I N G L E    M A S C O N     I N S I D E    T H E    P O L Y G O N     !!
!!                                                                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! allocate an array to hold all the ternary mascons inside the polygon
  allocate(mcon_tern_inside(1,ntern_inside,nvar_tern))
  write(message,'(a,i8,a)')"Single primary mascon inside the polygon has",ntern_inside,' ternary mascons'
  call status_update('STATUS','UTIL','blur_ocean_boundaries',msc_file_out,message,0)

! update the max_tern_per_prim variable if necessary
  if(ntern_inside > max_tern_per_prim)max_tern_per_prim = ntern_inside

! fill out the mcon_tern_inside array which is properly dimensioned for the number of ternary mascons
  ntern_inside_tmp = 0
  do itern = 1,ntern_inside
    mcon_tern_inside(1,itern,:) = mcon_tern_tmp(1,itern,:)  
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(info_flag < 0.d0)then
  write(message,'(a)')"All done. Blur_ocean_boundaries can stop now."
  call status_update('STATUS','UTIL','blur_ocean_mascons',msc_file_out,message,0)
  stop
endif

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
       "# blur_oceans    , region: ",region," created and contains",ntern_reshape &
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
  call status_update('STATUS','UTIL','blur_ocean_boundaries',msc_file_out,"Write out mascons outside of the polygon",0)
  do imsc = 1,total_prim 
    if(mcon_prim(imsc,6) > 1010.d0)then
      col_range(1) = 0.
      col_range(2) = 450.
    else if(mcon_prim(imsc,6) < 1010.d0)then
      col_range(1) = 750.
      col_range(2) = 1500.
    else
      col_range(1) = 500.
      col_range(2) = 700.
    endif
    prim_colour = nint(col_range(1) + rand(0) * (col_range(2)-col_range(1)))
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
                              ,mcon_prim(imsc,11),mcon_prim(imsc,10),mcon_region(imsc),pointers,prim_colour)
      enddo
    endif
  enddo
  write(message,'(a,i6,a)')"Written out ",nmsc_current,' primary mascons that are outside the polygon'
  call status_update('STATUS','UTIL','blur_ocean_boundaries',msc_file_out,message,0)
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
  call status_update('STATUS','UTIL','blur_ocean_boundaries',' ',message,0)

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
  colours(1:2) = 1600
  prim_flags(1) = "P"//region(1:5)
  mcon_prim(1,10) = 100.d0                     
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
    call write_mascon_record(lumsc_out,tern_flag(1),tern_number,nmsc_current+1,nint(mcon_prim(1,8)) &
                              ,msc_crds &
                              ,mcon_tern_inside(1,itern,4),mcon_tern_inside(1,itern,5),mcon_tern_inside(1,itern,8)    &
                              ,mcon_tern_inside(1,itern,6),mcon_prim(1,11) &
                              ,mcon_prim(1,10),region,pointers,colours)
  enddo
  write(message,'(a,a)')"Written out single mascon within polygon for region: ",region
  call status_update('STATUS','UTIL','blur_ocean_boundaries',msc_file_out,message,0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end





