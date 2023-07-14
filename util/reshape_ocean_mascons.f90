  program reshape_ocean_mascons

! program to a single primary ocean mascon and break it up into smaller primary mascons, so that tidal amplitudes can
! be estimated on specific mascons. 
!
! There are only a few places in the world's oceans where we think we are going to need to estimate corrections to the
! ocean tide model(s): Arctic and Antarctic oceans, plus perhaps some other regions (to be identified).
!
! This program has been written specifically to extract the ternary mascons from a single primary ocean mascon for Antarctica
! or the Arctic. Other regions can be added if required.
!
! P. Tregoning
! 17 May 2019

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

! PT190522: extra variables for ocean mascon reshaping
  real(kind=8) :: nonregion_area


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
    call status_update('WARNING','UTIL','reshape_ocean_mascons',' ', &
       'Runstring: reshape_ocean_mascons mascon_file mascon_file_out region max_area/tern rms_limit [n_untouched + msc codes]',0)
    call status_update('FATAL','UTIL','reshape_ocean_mascons',' ', &
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
      read(region(3:6),*)req_mascon
      write(message,'(a,i4,a,e8.2,a)')'Will reshape primary mascon number "',req_mascon,'" to size ',req_area, " km^2"
! PT170424: check whether it is the continental shelf region
    else if (region(1:2)  == "SH")then
      read(region(3:6),*)req_mascon
      write(message,'(a,i4,a,e8.2,a)')'Will reshape shelf mascon number "',req_mascon,'" to size ',req_area, " km^2"
    else
      req_mascon = 0
      write(message,'(a,a,a,e8.2,a)')'Will reshape primary mascon file in region "',region,'" to size ',req_area, " km^2"
    endif
    call status_update('STATUS','UTIL','reshape_ocean_mascons',mascon_file,message,0)
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
      call status_update('STATUS','UTIL','reshape_ocean_mascons',' ',message,0)
    endif
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  open the input ternary mascon file, read the header, allocate arrays and read in the mascons
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call read_msc_hdr(lumsc_in,mascon_file,msc_ocn_hdr_lines,max_ocn_hdr_lines,n_ocn_hdr_lines,msc_ocn_hdr_code &
                                ,max_ocean_prim,max_ocean_sec,max_ocean_tern &
                                ,max_ocean_tern_per_prim,max_ocean_tern_per_sec,max_ocean_sec_per_prim)
! PT190516: because we have not read in the full complement of ternarys, we don't have a good value for the non-ocean mascon values. Set them here.
  max_tern = 1485118
  max_prim = 10000
  call allocate_mascon_arrays
! PT180307: remove the trimlen in the name of the mascon file
  call read_ocean_mascons(lumsc_in,mascon_file)

! open the output file
  open(lumsc_out,file=msc_file_out,status='unknown')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! allocate some arrays
  ! now allocate the primary mascons for the region
  total_prim = max_ocean_prim
  allocate(region_prims(max_ocean_prim))    ! make it as big as it could possibly need to be.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! because we start with an ocean mascon file that contains only one primary mascon, we actually
! need only loop over all the ternarys within it to find out whether each ternary resides within
! the requested region ....
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! we first need to know how many ternarys there will be in the region requested. Read through them all
  ! and count them
  n_to_reshape = 1
  ntern_reshape = 0
  region_area = 0.d0
print*,'mcon_ocean_prim(1,1):',mcon_ocean_prim(1,1)," ",region(1:5),mcon_ocean_tern(100,1),mcon_ocean_tern(100,1)*180.d0/pi
  do itern=1,mcon_ocean_prim(1,1)

    ! Arctic ocean
    if(region(1:6) == "Arctic" .and. mcon_ocean_tern(itern,1)*180.d0/pi > 58.d0)then         ! it is an Arctic ocean ternary north of 58N
      ntern_reshape = ntern_reshape + 1
      region_area = region_area + mcon_ocean_tern(itern,4)
    else if(region(1:4) == "Anta" .and. mcon_ocean_tern(itern,1)*180.d0/pi < -60.d0)then   ! it is an Antarctic ocean ternary south of 60S
      ntern_reshape = ntern_reshape + 1
      region_area = region_area + mcon_ocean_tern(itern,4)
    else
    ! PT190517: add other regions here as required

    endif
  enddo
  write(message,'(a,a,a,i7)')'Number of ocean ternarys found in region ',region(1:6),": ",ntern_reshape
  call status_update('STATUS','UTIL','reshape_ocean_mascons','',message,0)

  ! now, dimension the required arrays for the new primary ocean mascons
  n_reshaped = nint(region_area/req_area)
  write(message,'(a,e10.4,a,i6,a)')"Region area (",region_area/1.d3," km^2) will be broken into",n_reshaped &
        ,' primary oean mascons.'
  call status_update('STATUS','UTIL','reshape_ocean_mascons',' ',message,0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! allocate the arrays, now that we know how many new ocean primaries there will be
  allocate(msc_reshaped_prims(n_reshaped,nvar_prim))
  allocate(msc_tern_reshaped(n_reshaped,nint(ntern_reshape/1.0),nvar_tern))

! RM190204: allocate arrays using no of reshaped mascons
  allocate(region_prims_reshape(n_reshaped))
  allocate(prim_flags_reshape(n_reshaped))

  do i=1,n_reshaped
      region_prims_reshape(i) = i+1  ! the new ocean primary  mascon numbers will start from 2
      prim_flags_reshape(i) = ""
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now, loop through all the ternarys again and siphon off the attributes of the ones in the right area
  region_area = 0.d0
  nonregion_area = 0.d0
  ntern_reshape = 0
  do itern=1,mcon_ocean_prim(1,1)

    ! Arctic ocean
    if(region(1:6) == "Arctic" .and. mcon_ocean_tern(itern,1)*180.d0/pi > 58.d0)then         ! it is an Arctic ocean ternary north of 58N
      ntern_reshape = ntern_reshape + 1
      msc_tern_reshaped(1,ntern_reshape,1:7) = mcon_ocean_tern(itern,1:7)
      tern_flag_ocean(itern) = "P-----"
      region_area = region_area + mcon_ocean_tern(itern,4)
    else if(region(1:4) == "Anta" .and. mcon_ocean_tern(itern,1)*180.d0/pi < -60.d0)then   ! it is an Antarctic ocean ternary south of 60S
      ntern_reshape = ntern_reshape + 1
      msc_tern_reshaped(1,ntern_reshape,1:7) = mcon_ocean_tern(itern,1:7)
      tern_flag_ocean(itern) = "P-----"
      region_area = region_area + mcon_ocean_tern(itern,4)

    else
    ! add up the ocean ternary area outside the region
      nonregion_area = nonregion_area + mcon_ocean_tern(itern,4)
    endif
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! for the primary mascons within the requested region, do something with them .....
  debug = .false.
  debug_rms = .true.
  n_to_reshape = 1

print*,region,rms_limit,cross_0E,n_to_reshape,n_reshaped,ntern_reshape,nint(ntern_reshape/1.0),debug,debug_rms
  call voronoi_reshape(.true.,region,rms_limit,cross_0E,n_to_reshape,region_prims_reshape,msc_reshaped_prims,n_reshaped &
                        ,msc_tern_reshaped,nint(ntern_reshape/1.0),debug,debug_rms)

  call srand(n_reshaped)
  col_range(1) = 0.
  col_range(2) = 450.

  do imsc = 1,n_reshaped
    ! assign a colour to each ternary in a primary
    prim_colour = nint(col_range(1) + rand(0) * (col_range(2)-col_range(1)))
    do itern=1,nint(msc_reshaped_prims(imsc,8))
      tern_number = nint(msc_tern_reshaped(imsc,itern,7))
      tern_colours(tern_number,:) = prim_colour
    enddo

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
    call status_update('WARNING','UTIL','reshape_mascons',' ',message,0)
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
  write(msc_ocn_hdr_lines(n_ocn_hdr_lines+1),'(a,a,a,i6,a,e12.6,a,1x,a8,1x,a10,1x,a5)') &
       "# reshape_ocean_mascons, region: ",region," has ",n_reshaped &
       ," primary mascons of size ~",req_area/1.d0," m^2. Timetag: ",date,time,timezone
! generate new unique code for the output program
  call hash_header(msc_ocn_hdr_lines(2:n_ocn_hdr_lines+1),n_ocn_hdr_lines+1,hashcode)

! write the new code and mascon numbers to the first line of the header
  write(msc_hdr_lines(1)(1:81),'(a,i6,5i8)')hashcode(1:8),new_n_prim,new_n_sec,max_ocean_tern,max_sec_per_prim &
                                           ,new_max_tern_per_prim,new_max_tern_per_prim

! and write out all header lines to the output file
  do i=1,n_ocn_hdr_lines+1
    write(lumsc_out,'(a)')msc_ocn_hdr_lines(i)
    write(*,'(i6,a)')i,msc_ocn_hdr_lines(i)
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write out all the mascons. 
! First, create a primary containing all the untouched ternarys
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(message,'(a,i6,a)')"Writing out ",n_to_reshape," untouched ternary ocean mascons"
  call status_update('STATUS','UTIL','reshape_ocean_mascons',msc_file_out,message,0)
  nmsc_current = 0
  ! first, the primary mascon line. We don't have some of the information, so just make it up here ....
  nmsc_current = 1
  nsec = 1
  msc_crds(1) = 90.d0
  msc_crds(1) =  0.d0
  msc_crds(3)   = 6378000.d0
  pointers = nmsc_current
  colours(1:2) = 0
  tide_flag    = 0  ! PT170606: for now, don't turn on tidal amplitude modelling
  ! primary mascon line
  call write_mascon_record(lumsc_out,prim_flags_ocean(1),nmsc_current,nsec,mcon_ocean_prim(1,1)-ntern_reshape,msc_crds     &
                           ,nonregion_area,0.d0,0.d0,1029.d0   &
                           ,0.d0,tide_flag,mcon_region(imsc),pointers,colours)
  ! secondary mascon line
  sec_flag_tmp = "S"//prim_flags_ocean(1)(2:6)
  call write_mascon_record(lumsc_out,sec_flag_tmp,nmsc_current,nsec,mcon_ocean_prim(1,1)-ntern_reshape,msc_crds &
                              ,nonregion_area,0.d0,0.d0,1029.d0   &
                              ,0.d0,tide_flag,mcon_region(imsc),pointers,colours)
  ! all the ternary records
  do itern = 1,mcon_ocean_prim(1,1) 
    if(prim_flags_ocean(imsc) /= "P-----")then
      pointers(1) = nmsc_current
      pointers(2) = nmsc_current
      ntern = nint(mcon_ocean_tern(itern,7))
      msc_crds(1:2) = mcon_ocean_tern(itern,1:2)*180.d0/pi
      msc_crds(3)   = mcon_ocean_tern(itern,3)
      call write_mascon_record(lumsc_out,tern_flag_ocean(itern),itern,nsec,mcon_ocean_prim(1,1),msc_crds &
                       ,mcon_ocean_tern(itern,4),mcon_ocean_tern(itern,5),mcon_ocean_tern(itern,8) &
                       ,mcon_ocean_tern(itern,6),mcon_ocean_prim(imsc,11),mcon_ocean_prim(imsc,10),mcon_region(imsc) &
                       , pointers,tern_colours(ntern,:))
    endif
  enddo
  write(message,'(a,i6,a)')"Written out",1,' primary mascon of untouched ocean ternary mascons'
  call status_update('STATUS','UTIL','reshape_ocean_mascons',msc_file_out,message,0)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! now the reshaped ones for the region
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(message,'(a,i6,a)')"Writing out ",n_reshaped,"  reshaped primary mascons"
  call status_update('STATUS','UTIL','reshape_ocean_mascons',msc_file_out,message,0)
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
          prim_flags_reshape(imsc)(2:6) = tern_flag(nint(msc_tern_reshaped(imsc,1,7)))(2:6)
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




