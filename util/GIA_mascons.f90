program GIA_mascons

! Generates a GIA mascon file that contains one primary mascon in region of interest, output file can then be used  
! as input in reshape_mascons. The height of the GIA ternary mascons are defined at the water-rock interface whose
! depth are extracted from BEDMAP and GEBCO models.
! 
! The GIA_mascons program performs the following tasks:
! 1. Opens and reads input stage4 mascon file
! 2. Ternary mascons that exist in the specified region are duplicated to from the GIA ternary mascons
!    i.   The radius of the ocean and land ternaries are adjusted for either:
!         a) The geoid-ellipsoid seperation and ice-thickness/altitude for land/ocean ternaries, or
!         b) A user defined depth in km below the ellipsoid surface for all ternaries
!    ii.  Density of all ternary mascons are set to 3300 kg/m^3
!    iii. Number of total ternaries is increased
!    iv.  All other ternary mascon characteristics are unchanged 
! 3. A single GIA primary mascon is created which contains all GIA ternary mascons, radius is set to ...
! 4. The new GIA primary mascon is coloured
! 5. Mascons are written to an intermediate GIA mascon outfile, ready to be used as input to reshape_mascons
!
! R. McGirr
! 5 April 2019
!
! TO DO:
! 1. Read in ice thickness values for non-zero ternaries, subtract from GIA land ternaries
! 2. Adjust GIA mascon radius under ocean for altitude 
! 3. Figure out appropriate radius to set primary mascon

use mascon_mod     ! provides the declarations of the variables needed to read in and store all the mascons

implicit none
! command line arguments
  character*150             :: mascon_file           ! input mascon file
  character*150             :: output_file           ! output GIA mascon file
  character*15              :: region                ! GIA mascon region
  real(kind=8)              :: depth                 ! depth of GIA mascons
  logical                   :: ice                   ! GIA mascons at rock-water interface

! mascon characteristics and counters
  integer*4                 :: imsc,itern,i,j,k      ! mascon counters
  integer*4                 :: ngia,nice             ! number of GIA mascons and mascons with ice cover
  integer*4                 :: nocean,nlandice       ! number of ocean and land ice GIA mascons
  integer*4                 :: nmsc,nsec,ntern       ! number of prim, sec and tern mascons
  integer*4                 :: extra_prim            ! temporary random ternary number
  real(kind=8)              :: density = 3300.d0     ! desnity of GIA mascons
  real(kind=8)              :: rad_e                 ! temp radius of the ellipsoidal surface
  real(kind=8)              :: rad_g                 ! temp radius of the GIA surface
  real(kind=8)              :: region_area           ! area of the GIA mascon
  character*150             :: ice_file              ! input ice thickness file
  character*150             :: hdr_lines(2)          ! vector to store hdr lines of ice file

! parameters used to define GIA region
  integer*4,   allocatable  :: region_prims(:)       ! list of primary mascons found within the requested region
  integer*4,   allocatable  :: region_prims_gia(:)   ! list of primary mascons after reshaping
  integer*4,   allocatable  :: tern_gia_colours(:,:) ! colour of corresponding sec and prim mascon for each GIA tern  
  real(kind=8),allocatable  :: ternary_ice(:,:)    ! contains ice thickness and associated ternary number
  real(kind=8),allocatable  :: msc_gia_prims(:,:)    ! gia primary mascons
  real(kind=8),allocatable  :: msc_gia_terns(:,:,:)  ! gia ternary mascons
  real(kind=8),allocatable  :: tern_crds(:,:)        ! array of all the ternary mascons and attributes in GIA region
  character*6, allocatable  :: prim_flags_gia(:)     ! 6 character flags for type of GIA prim
  character*6, allocatable  :: tern_flags_gia(:)     ! 6 character flags for type of GIA tern
  logical,     allocatable  :: rand_numbers(:)       ! flags to say whether a particular rand number has been used yet

! parameters for counting changes in mascon numbers
  integer*4                 :: new_total_prim        ! new total number of primary mascons
  integer*4                 :: new_total_sec         ! new total number of secondary mascons
  integer*4                 :: new_total_tern        ! new total number of ternary mascons
  integer*4                 :: new_max_tern_per_prim ! new max number of ternarys per primary
  integer*4                 :: n_reshaped            ! number of primary GIA mascons after reshaping
  integer*4                 :: ntern_reshape         ! number of ternary GIA mascons to be reshaped

! parameters used to output mascons
  integer*4                 :: colours(2)            ! colour of secondary and primary mascon
  integer*4                 :: pointers(2)           ! tern to prim and tern to sec pointers
  integer*4                 :: prim_colour           ! colour for gia primary
  integer*4                 :: tern_number           ! temp storage of the actual ternary mascon number
  integer*4                 :: tide_flag             ! set to zero for now (ie no tidal amplitude modelling
  real(kind=8)              :: msc_crds(3)           ! temp storage of mascon lat/lon/rad
  real(kind=8)              :: col_range(2)          ! variable used to calculate new colour for primary mascons
  character*6               :: sec_flag              ! single, temporary value assigned to the type pf secondary mascon
  character*12              :: hashcode              ! character code of output mascon file

! constants for calculating mascon radius and area
  real(kind=8), parameter   :: pi = 4.d0*atan(1.d0)  ! approximation of pi
  real(kind=8), parameter   :: Ra = 6378.137d3       ! semi-major axis WGS84 ellipsoid
  real(kind=8), parameter   :: Rb = 6356.7523142d3   ! semi-minor axis WGS84 ellipsoid

! unit numbers for input and output mascon files
  integer*4,    parameter   :: lu_ice    = 9         ! unit number for ice thickness file
  integer*4,    parameter   :: lumsc_in  = 10        ! unit number of input mascon file
  integer*4,    parameter   :: lumsc_out = 11        ! unit number of output mascon file

! variables for writing date and time into header of output file
  character(8)              :: date
  character(10)             :: time
  character(5)              :: timezone

! variables used for debugging, errors and writing progress to console
  integer*4                 :: ioerr
  character*50              :: arg
  character*150             :: message
  logical                   :: debug = .true.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! decipher the command line arguments
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call getarg(1,mascon_file)
  if ( mascon_file == " " ) then
    call status_update('WARNING','UTIL','GIA_mascons',' ','Runstring: GIA_mascons mascon_file output_file region depth',0)
    call status_update('WARNING','UTIL','GIA_mascons',' ','E.g.: mascons_stage4 mascons_stage4_GIA Antarctica ice',0)
    call status_update('FATAL','UTIL','GIA_mascons',' ','Either specify ice (ice/rock interface) or depth (km below ellipsoid surface)',0)
  endif
  call getarg(2,output_file)
  call getarg(3,region)
  if (region(1:5) == "Antar") then
    write(message,'(a,e12.6,a)')'Will create Antarctic GIA mascons of density: ',density,' kg/m^3'
    call status_update('STATUS','UTIL','GIA_mascons',' ',message,0)
  elseif (region(1:5) == "Green") then
    call status_update('FATAL','UTIL','GIA_mascons',' ','Region not coded, exiting',0)
  elseif (region(1:5) == "Laure") then
    call status_update('FATAL','UTIL','GIA_mascons',' ','Region not coded, exiting',0)
  elseif (region(1:5) == "Fenno") then
    call status_update('FATAL','UTIL','GIA_mascons',' ','Region not coded, exiting',0)
  else
    call status_update('FATAL','UTIL','GIA_mascons',' ','Region options: Antarctica, Greenland, Laurentide, Fennoscandia',0)
  endif
  call getarg(4,arg)
    if (arg(1:3) == "ice") then
      ice = .true.
      write(message,'(a)')'GIA mascon surface will be adjusted to the depth of the ice/rock interface'
      call status_update('STATUS','UTIL','GIA_mascons',' ',message,0)
    else
      read(arg,*)depth
      write(message,'(a,e12.6,a)')'GIA mascon surface will be adjusted to ',depth,' km below the ellipsoid surface'
      call status_update('STATUS','UTIL','GIA_mascons',' ',message,0)
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! open and read ice thickness file 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  if (ice .eqv. .true.) then
    ice_file = "ternary_ice_point"
    open(lu_ice,file=ice_file,status='old')
! read header lines including number of ice ternaries
    read(lu_ice,'(a)')hdr_lines(1)
    read(hdr_lines(1)(8:13),*)nice 
    read(lu_ice,'(a)')hdr_lines(2) 
! allocate array that holds ice thickness and tern number
    allocate(ternary_ice(nice,2))
! read in ice thicknesses corresponding to each ternary
    do imsc = 1, nice
      read(lu_ice,*)ternary_ice(imsc,1:2)
!print*,ternary_ice(imsc,:)
    enddo
    close(lu_ice)
    if(debug)print*,'DEBUG: nice:',nice
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! open and read input file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! read mascon file header 
  call read_msc_hdr(lumsc_in,mascon_file,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec, &
                    max_tern,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
! allocate mascon arrays
  call allocate_mascon_arrays
! read mascon file
  call read_mascon_file(lumsc_in,mascon_file)
! open the output file
  open(lumsc_out,file=output_file,status='unknown')

! allocate some arrays
  allocate(region_prims(total_prim))    ! make it as big as it could possibly need to be

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Loop over all ternary mascons in each primary and identify the ones within the region of interest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ntern = 0
  nmsc = 0
  region_area = 0.d0
  region_prims = 0

! first, loop through all primarys and find any that contain a ternary in the region of interest
  do imsc = 1, total_prim
    if( any(mcon_tern(imsc,:,1) < -60.d0*pi/180.d0 )) then
      nmsc = nmsc + 1
      region_prims(nmsc) = imsc
    endif
  enddo

! now, just loop through the primarys that contain ternarys in the region of interest and find specific ternarys
  do imsc = 1, nmsc
    do itern = 1, nint(mcon_prim(region_prims(imsc),8))
      if ( mcon_tern(region_prims(imsc),itern,1) < -60.d0*pi/180.d0 ) then ! test program with small area
        ntern = ntern + 1
        region_area = region_area + mcon_tern(region_prims(imsc),itern,4) ! need to do area adjustment
      endif
    enddo 
  enddo

  if(debug)print*,'DEBUG: nmsc:',nmsc,'ntern:',ntern,' region_area:',region_area

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do something with the mascons in the GIA region
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

! create one GIA mascon
  ngia = 1

! allocate the arrays
  allocate(msc_gia_terns(ngia,ntern,nvar_tern))
  allocate(msc_gia_prims(ngia,nvar_prim))
  allocate(region_prims_gia(ngia))
  allocate(prim_flags_gia(ngia))
  allocate(tern_flags_gia(ntern))
  allocate(rand_numbers(ntern))
  
  region_prims_gia(ngia) = 0
  prim_flags_gia(ngia) = ""
  msc_gia_prims(ngia,:) = 0.d0

! transfer all the ternary mascons into one array
  nocean = 0
  nlandice = 0
  itern = 0

  do i = 1, nmsc
    do j = 1, nint(mcon_prim(region_prims(i),8))
      if( mcon_tern(region_prims(i),j,1) < -60.d0*pi/180.d0 ) then
      	itern = itern + 1
! lat, lon, altitude and geoid-ellipsoid separation remain the same, density set to upper mantle, ntern increases
      	msc_gia_terns(ngia,itern,1:2) = mcon_tern(region_prims(i),j,1:2)       ! lat, lon in radians
      	msc_gia_terns(ngia,itern,5) = mcon_tern(region_prims(i),j,5)           ! altitude in m
      	msc_gia_terns(ngia,itern,6) = density                                  ! density of upper mantle
      	msc_gia_terns(ngia,itern,7) = mcon_tern(region_prims(i),j,7) + ntern   ! no. of ternarys
      	msc_gia_terns(ngia,itern,8) = mcon_tern(region_prims(i),j,8)           ! geoid-ellipsoid separation

! adjust radius depending on option specified by user
        if(ice .eqv. .true.) then
! radius of all mascons are adjusted for the geoid/ellipsoid separation
          msc_gia_terns(ngia,itern,3) = mcon_tern(region_prims(i),j,3) + mcon_tern(region_prims(i),j,8)
! radius of all ocean mascons are further adjusted for altitude
      	  if( mcon_tern(region_prims(i),j,6) > 1010.d0 ) then
            msc_gia_terns(ngia,itern,3) = msc_gia_terns(ngia,itern,3) + mcon_tern(region_prims(i),j,5)
            nocean = nocean + 1
! radius of land mascons are further adjusted for their associated ice height
          else
            do k = 1, nice
              if( mcon_tern(region_prims(i),j,7) == ternary_ice(k,1) )then
                msc_gia_terns(ngia,itern,3) = msc_gia_terns(ngia,itern,3) - ternary_ice(k,2)
                nlandice = nlandice + 1
              endif                                                                  
            enddo
          endif
        else
! radius adjusted by user-defined depth below the ellipsoid surface
          msc_gia_terns(ngia,itern,3) = mcon_tern(region_prims(i),j,3) - depth*1.e3
        endif

! area is adjusted by the ratio of the ellipsoidal radius to the GIA radius squared
        rad_e = mcon_tern(region_prims(i),j,3)
        rad_g = msc_gia_terns(ngia,itern,3)
      	msc_gia_terns(ngia,itern,4) = mcon_tern(region_prims(i),j,4) * (rad_g/rad_e)**2
      endif
    enddo
  enddo
  if(debug)print*,'DEBUG: total land ternaries with ice cover, without ice cover, ocean',nlandice, &
                   ntern-nocean-nlandice,nocean
  if(debug)print*,'DEBUG: msc_gia_terns:',msc_gia_terns(ngia,ntern,1:8)

! fill in primary mascon info
  msc_gia_prims(ngia,8) = ntern   ! number of ternarys in primary
  msc_gia_prims(ngia,6) = density ! desnity of GIA mascon
  msc_gia_prims(ngia,11) = -1.d0   ! percent land in GIA mascons

! choose a random ternary mascon to use for the initial CoM
  call srand(ntern)
  rand_numbers = .false.
  call get_random_value(rand_numbers,ntern,extra_prim)

! assign the coord information
  msc_gia_prims(ngia,1:3) = msc_gia_terns(ngia,extra_prim,1:3)

! calculate mascon centre of mass
  call MASCON_calc_COM(ntern,nvar_prim,nvar_tern,msc_gia_terns(ngia,:,:),msc_gia_prims(ngia,:))

! assign mascon number to gia primary mascon that was created
  if(ngia > nmsc)region_prims_gia(ngia) = ngia 

! check whether the max ternarys per primarys need updating
  if(nint(msc_gia_prims(ngia,8)) > new_max_tern_per_prim)new_max_tern_per_prim = nint(msc_gia_prims(ngia,8))
  
  if(debug)print*,'DEBUG: msc_gia_prims:',msc_gia_prims(ngia,1:12)

! calculate new total primary and secondary mascons  
  new_total_prim = total_prim + ngia
  new_total_sec = new_total_prim
  new_total_tern = max_tern + ntern

  if(debug)print*,'DEBUG: new_total_prim:',new_total_prim,'new_total_tern:',new_total_tern

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Colour GIA mascons
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(tern_gia_colours(new_total_tern,2))

  call srand(ngia)
! colour all gia_prims, don't discriminate between ocean or continent
  col_range(1) = 500.
  col_range(2) = 750.
! assign a colour to each individual ternary within the same primary
  prim_colour = nint(col_range(1) + rand(0) * (col_range(2)-col_range(1))) ! find random colour for each primary
  do itern = 1, nint(msc_gia_prims(ngia,8)) ! loop through each ternary belonging to current primary
    tern_number = nint(msc_gia_terns(ngia,itern,7))
    tern_gia_colours(tern_number,:) = prim_colour
  enddo

! check whether the max ternarys per primarys need updating
  if(nint(msc_gia_prims(ngia,8)) > new_max_tern_per_prim)new_max_tern_per_prim = nint(msc_gia_prims(ngia,8))

! assign mascon numbers to gia primary mascons that were created
  if(ngia > nmsc)region_prims_gia(ngia) = ngia
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Output all mascons
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! add a new line to the header to represent what has been done in running this program
  call date_and_time(date,time,timezone)
  write(msc_hdr_lines(n_hdr_lines+1),'(a,a,a,i6,a,e12.6,a,1x,a8,1x,a10,1x,a5)')"# GIA_mascons    , region: ", &
        region," has ",ngia," primary GIA mascon of size ~",region_area/1.d0," m^2. Timetag: ",date,time,timezone
! generate new unique code for the output program
  call hash_header(msc_hdr_lines(2:n_hdr_lines+1),n_hdr_lines+1,hashcode)

! write the new code and mascon numbers to the first line of the header
  write(msc_hdr_lines(1)(1:81),'(a,i6,5i8)')hashcode(1:8),new_total_prim,new_total_sec,new_total_tern,max_sec_per_prim, &
        new_max_tern_per_prim,new_max_tern_per_prim

! and write out all header lines to the output file
  do i=1,n_hdr_lines+1
    write(lumsc_out,'(a)')msc_hdr_lines(i)
  enddo
  
! write out all the mascons, first the untouched
  write(message,'(a,i6,a)')"Writing out ",total_prim," untouched primary mascons"
  call status_update('STATUS','UTIL','GIA_mascons',output_file,message,0)
  nmsc = 0
  do imsc = 1,total_prim
    if(prim_flags(imsc) /= "P-----")then
      nmsc = nmsc + 1
! write out primary record
      nsec = 1
      msc_crds(1:2) = mcon_prim(imsc,1:2)*180.d0/pi
      msc_crds(3) = mcon_prim(imsc,3)
      pointers = nmsc
      colours(1:2) = 0
      tide_flag = 0
      call write_mascon_record(lumsc_out,prim_flags(imsc),nmsc,nsec,nint(mcon_prim(imsc,8)),msc_crds, &
           mcon_prim(imsc,4),mcon_prim(imsc,5),mcon_prim(imsc,12),mcon_prim(imsc,6),mcon_prim(imsc,11), &
           mcon_prim(imsc,10),mcon_region(imsc),pointers,colours)
      pointers(1) = nmsc
      pointers(2) = 0
! write out secondary record
      sec_flag = "S"//prim_flags(imsc)(2:6)
      call write_mascon_record(lumsc_out,sec_flag,nmsc,nsec,nint(mcon_prim(imsc,8)),msc_crds, &
           mcon_prim(imsc,4),mcon_prim(imsc,5),mcon_prim(imsc,12),mcon_prim(imsc,6),mcon_prim(imsc,11), &
           mcon_prim(imsc,10),mcon_region(imsc),pointers,colours)
! write ternary records for this secondary
      do itern = 1,nint(mcon_sec(imsc,1,8))
        pointers(1) = nmsc
        pointers(2) = nmsc
        ntern = nint(mcon_tern(imsc,itern,7))
        msc_crds(1:2) = mcon_tern(imsc,itern,1:2)* 180.d0/pi
        msc_crds(3) = mcon_tern(imsc,itern,3)
        call write_mascon_record(lumsc_out,tern_flag(ntern),ntern,nsec,nint(mcon_prim(imsc,8)),msc_crds, &
             mcon_tern(imsc,itern,4),mcon_tern(imsc,itern,5),mcon_tern(imsc,itern,8),mcon_tern(imsc,itern,6), &
             mcon_prim(imsc,11),mcon_prim(imsc,10),mcon_region(imsc),pointers,tern_colours(ntern,:))
      enddo
    endif
  enddo
  write(message,'(a,i6,a)')"Written out ",total_prim," untouched mascons to output"
  call status_update('STATUS','UTIL','GIA_mascons',output_file,message,0)
  
! now write out the GIA mascons
  write(message,'(a,i6,a)')"Writing out",1," GIA primary mascon"
  call status_update('STATUS','UTIL','GIA_mascons',output_file,message,0)
  nmsc = nmsc + 1
! write out primary record
  nsec = 1
  msc_crds(1:2) = msc_gia_prims(ngia,1:2)*180.d0/pi
  msc_crds(3) = msc_gia_prims(ngia,3)
  pointers = nmsc
  colours(1:2) = 0
  prim_flags_gia(ngia) = "PGIA  "
  msc_gia_prims(ngia,10) = 0
  call write_mascon_record(lumsc_out,prim_flags_gia(ngia),nmsc,nsec,nint(msc_gia_prims(ngia,8)), &
       msc_crds,msc_gia_prims(ngia,4),msc_gia_prims(ngia,5),msc_gia_prims(ngia,12),msc_gia_prims(ngia,6), &
       msc_gia_prims(ngia,11),msc_gia_prims(ngia,10),region,pointers,colours)
! write out secondary record
  pointers(1) = nmsc
  pointers(2) = 0
  nsec = nmsc
  colours = 0
  sec_flag = "S"//prim_flags_gia(ngia)(2:6)
  call write_mascon_record(lumsc_out,sec_flag,nmsc,nsec,nint(msc_gia_prims(ngia,8)),msc_crds, &
       msc_gia_prims(ngia,4),msc_gia_prims(ngia,5),msc_gia_prims(ngia,12),msc_gia_prims(ngia,6), &
       msc_gia_prims(ngia,11),msc_gia_prims(ngia,10),region,pointers,colours)
! write ternary records for this secondary
  do itern = 1, nint(msc_gia_prims(ngia,8))
    pointers(1) = nmsc
    pointers(2) = nsec
    msc_crds(1:2) = msc_gia_terns(ngia,itern,1:2)*180.d0/pi
    msc_crds(3) = msc_gia_terns(ngia,itern,3)
    ntern = nint(msc_gia_terns(ngia,itern,7))
    tern_flags_GIA(itern) = "TGIA  "
    call write_mascon_record(lumsc_out,tern_flags_GIA(itern),ntern,nsec,nint(msc_gia_prims(ngia,8)),msc_crds, &
         msc_gia_terns(ngia,itern,4),msc_gia_terns(ngia,itern,5),msc_gia_terns(ngia,itern,8), &
         msc_gia_terns(ngia,itern,6),msc_gia_prims(ngia,11),msc_gia_prims(ngia,10),region,pointers, &
         tern_gia_colours(ntern,:))
  enddo
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End of program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
  write(message,'(a)')'Normal end of GIA_mascons'
  call status_update('STATUS','UTIL','GIA_mascons',' ',message,0)

end