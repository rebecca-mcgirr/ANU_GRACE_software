  program separate_mascons

! program to re-distribute ternary mascons so that all primary mascons are either all ocean or all land.
! It goes like this:
!
! 1. We start with a mascon file that contains a generic configuration of primary/ternary mascons
! 2. For the ternary mascons we know
!      2.1  the height/depth of the ternary mascon (i.e. is it land or water)
!      2.2  the primary mascon to which it belongs currently
! 3. Compute whether each primary mascon is predominantly land or water (is this done already?)
! 4. Identify which primary mascons are a mix of land and water. For each of the ternary mascons in those mixed mascons:
!      4.1  calculate the distance to the surrounding primary mascon coordinates
!      4.2  reassign the minority ternary mascon to the nearest primary mascon of the same type
! 5. Once all primary mascons contain only the same type of ternary mascons we are done
! 6. Recompute the properties of the adjusted primary mascons to get new primary mascon coords and surface area.
!
!  The output file can then be used as input into the voronoi_cell program to reconfigure the shape/area of the primary mascons. This
!  process ensures that our primary mascons do not cross coastlines, since we don't allow a mix of land and ocean ternary mascons.
!
!  Adapted from program configure_mascons (essentially updated to use the new mascon file)
!
!  P. Tregoning  (with help from Herb McQueen for the algorithm concept)
!  28 October 2016
!
! PT210902: rather than merging the minority of a mascon affected by the coastline and keeping the - now small - majority,
!           merge the ternaries of the small majority part as well. That is, eliminate the mixed mascons
!           altogether and merge them into nearby mascons of similar type.

  use mascon_mod

  implicit none

! input/output files
  character :: msc_input*150       ! name of input  primary mascon file containing primary mascon coords
  character :: msc_output*150      ! name of output primary mascon file containing the recomputed CoM of the adjusted primary mascons
  character :: seed_file*100       ! file of additional primary mascons to soak up the difficult mixed regions of the world (Iceland, NZ, islands etc)
  integer*4 :: min_terns           ! minimum number of ternarys per primary required for primary mascon to survive

! maximum number of surrounding mascons
  integer*4,parameter :: max_surround_msc=1000

! primary mascon coords
  real(kind=8),allocatable :: msc_prim_in(:,:)       ! variable to store the decimal input  primary mascon coordinates. lat/lon/radius/area/height/density/% land/#ternaries/bit-mapped tide
  real(kind=8),allocatable :: msc_prim_out(:,:)      ! variable to store the decimal output primary mascon coordinates. lat/lon/radius/area/height/density/% land/#ternaries/bit-mapped tide 
  integer*4 :: nland                                 ! number of ternary mascons within a primary mascon that are land
  integer*4 :: current_total_prim                    ! variable to track the number of primary mascons that still exist (goes down if all ternary mascons are continental shelf)
! number of secondary mascons
  integer*4 :: nsec

! ternary mascon coords
  real(kind=8),allocatable :: msc_tern_in(:,:,:)  ! variable to store the decimal input  ternary mascon coordinates. Lat, lon, rad, area, density, depth
  real(kind=8),allocatable :: msc_tern_out(:,:,:) ! variable to store the decimal output ternary mascon coordinates. Lat, lon, rad, area, density, depth
  logical   :: transfer_tern                      ! indicates whether a ternary mascon can be transferred to a surrounding primary mascon.
  real(kind=8) :: tern_to_primary(max_surround_msc)             ! distances from ternary mascon to surrounding primary mascons

! PT161104: variable to increase the max number of ternary mascons above that of the input file value
  integer*4, parameter :: msc_tern_mult=10
  integer*4            :: new_max_tern_per_prim

! output variables
  real(kind=8) :: density             ! density of water for this primary mascon (based on whether ocean or land)
  integer*4 :: tides                  ! indicator to graceorb as to whether or not to estimate tidal constituents (only ON for shallow water primary mascons)
  character :: prim_flag*5            ! flag to indicate type of primary mascon [Deep, Land, Shal(low) ]
  integer*4 :: tmp_sec_colour,prim_colour ! colour values for plotting the mascons
  integer*4 :: new_prim_total         ! count of the total number of primary mascons remaining after the land/ocean reconfiguration

! variables for mascons surrounding a mixed mascon
  integer*4    :: imsc                  ! do loop index
  real(kind=8) :: msc_surrounding(max_surround_msc,4 )! primary mascon coords of those surrounding a mixed mascon. 2nd dimesion variables are: lat,lon,primary msc #, % land
  real(kind=8) :: tmpnewlon,dlat,dlon   ! varaibles for determining the surrounding mascons
  integer*4    :: n_surround            ! number of surrounding mascons within a given distance (determined in MASCON_nearest_msc)
  logical             :: found
  logical,allocatable :: prim_in_type(:)       ! logical variable to be true if predominantly land, false if predominantly ocean
  logical             :: surround_type(max_surround_msc)     ! logical for surrounding mascons. True if land, false if ocean
  integer*4           :: give_to               ! temporary variable being the primary mascon number to receive a transferred ternary mascon
  real(kind=8)        :: min_dist              ! minimum distance between a mascon and its surrounding mascons.
  real(kind=8)        :: tmp_dist
  integer*4           :: ntern_given           ! number of ternary mascons of a single primary mascon that are transferred to different primary mascons
  integer*4           :: ntern_ocean           ! count the number of ocean ternary mascons in a land primary mascon
  integer*4           :: ntern_land            ! count the number of land ternary mascons in an ocean primary mascon
  real(kind=8)        :: max_distance          ! maximum distance away for which ternary mascons can be transferred into this mascon

! variables for continental shelf mascons
  integer*4            :: n_shelf                ! number of continental shelf ternary mascons
  real(kind=8)         :: shelf_depth            ! maximum depth of a mascon deemed to be on the continental shelf
  integer*4, parameter :: max_shelf = 100000     ! 80000 is enough for a continental shelf depth of -200 m
  real(kind=8)         :: msc_shelf(max_shelf,8) ! storage of continental shelf mascon information
  integer*4            :: tmp_shelf

! variables for inland seas
  integer*4            :: n_casp,n_eyre          ! number of ternary mascons in the Caspian Sea and Lake Eyre
  integer*4,parameter  :: max_seas = 10000        ! maximum number of ternary mascons in each inland sea
  real(kind=8)         :: msc_casp(max_seas,8)   ! storage for Caspian sea ternary mascons
  real(kind=8)         :: msc_eyre(max_seas,8)   ! storage for Lake Eyre ternary mascons
  integer*4            :: tmp_other 

! file unit numbers
  integer*4 :: lumsc_in, lumsc_out, lutern_in, lutern_out

! variables associated with seed primary mascons
  integer*4            :: luseed               ! unit number of seed file
  integer*4, parameter :: max_seed = 1000      ! max_number of seed mascons
  integer*4            :: nseed                ! actual number of seed mascons in file
  real(kind=8)         :: msc_seed(max_seed,6) ! coords and density of seed mascons (colatdeg,colatmin,londeg,lonmin,density,transfer_dist)

! variables associated with identifying the perimeter of primary mascons
  integer*4 :: nvertex
  integer*4 :: vertex_prim(500)                ! numbers of ternary mascons that make up the perimeter of the primary mascon (clockwise order, apparently)

! variables to write out the mascon information
  real(kind=8) :: msc_crds(3)
  character*6  :: msc_code
  integer*4    :: pointers(2),colours(2)
  character*15 :: region

  real(kind=8), allocatable :: msc_primary_type(:,:)   ! percentage of primary mascon that is ternary land mascons

! counters
  integer*4 :: i,j,k
  integer*4 :: imsc_prim       ! number of input primary mascon when reading the input ternary mascon file
  integer*4 :: ntern           ! number of input ternary mascons per input primary mascon
  integer*4 :: nmixed          ! number primary mascons that are a mix of land and water
  integer*4 :: itern

! other stuff
  integer*4 :: ioerr
  character :: message*250
  character :: arg*100
  real(kind=8) :: tmprad,tmparea,tmpdensity,tmpdepth,tmpdummy
  integer*4    :: trimlen
  real(kind=8) :: pi

! debug variables
  real(kind=8) :: percent_land
  logical   :: debug




 debug = .false.
 pi = 4.d0*datan(1.d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      decode command line
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call getarg(1,msc_input)
  if (msc_input(1:1) == " ")then
    write(message,*)'Runstring: separate_mascons msc_infile msc_outfile shelf_depth min_terns'
    call status_update('FATAL','UTIL','msc_ocean_vs_water',' ',message,0)
  endif
  call getarg(2,msc_output)
  call getarg(3,arg)
  read(arg,*)shelf_depth
  call getarg(4,arg)
  read(arg,*)min_terns
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      open files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  lumsc_in   = 10
  lumsc_out  = 11

! output mascon file
  open(lumsc_out,file=msc_output,status='unknown',iostat=ioerr)
  if (ioerr /= 0)then
    call status_update('FATAL','UTIL','separate_mascons',msc_output,'Error opening output mascon file',0)
  endif

  call status_update('STATUS','UTIL','separate_mascons',msc_output,'Opened output mascon file',0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      read the mascon information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call read_msc_hdr(lumsc_in,msc_input,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                                ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
  call allocate_mascon_arrays
  call read_mascon_file(lumsc_in,msc_input)
! PT161111: set up a local variable that tracks the number of total primary mascons currently. It starts as total_prim but can decrease if all ternary 
!           mascons are continental shelf (and, therefore, are kicked out of the primary mascon)
  current_total_prim = total_prim
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      allocate the mascon arrays for the output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT170531: this is now defined in mod_subs/mascon_mods.f90
!  nvar_prim = 10
  allocate(msc_prim_in(total_prim,nvar_prim))          ! 2nd dimension variables are: lat/lon/rad/area/hgt/density/#sec/#tern/#first secondary/tidal flag/%land/geoid
  allocate(msc_prim_out(total_prim,nvar_prim))         ! 2nd dimension variables are: lat/lon/rad/area/hgt/density/#sec/#tern/#first secondary/tidal flag/%land/geoid
  allocate( msc_tern_in(total_prim*msc_tern_mult,max_tern_per_prim*msc_tern_mult,nvar_tern))  ! allow up to max_tern ternary mascons per primary mascon. 3rd dimension variables are lat/lon/radius/area/depth/density/geoid.
  allocate(msc_tern_out(total_prim*msc_tern_mult,max_tern_per_prim*msc_tern_mult,nvar_tern))  ! allow up to max_tern ternary mascons per primary mascon. 3rd dimension variables are lat/lon/radius/area/depth/density/geoid.
  allocate(prim_in_type(total_prim))             ! will be set to true if predominantly land, false if ocean
!  allocate(max_dists(total_prim))                ! set to 700 km for input primary mascons, value in seed file for seed primary mascons
  max_distance = 1000000.    ! set to 2,000 km
! allocate the array of primary (PDeep/PShelf/PCasp/PEyre/PLand) and ternary mascon flags (TDeep/TShelf/TCasp/TEyre/TLand)
!  allocate(prim_flags(max_prim+10))  ! increase it by 10 beyond the actual number of primary mascons in the input file (just to be sure ...)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   transfer the input primary mascon information to msc_prim_in array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT161028: note that the order of msc_prim_in from 7-9 is NOT the same as for mcon_prim ......
  msc_prim_in(1:total_prim,1:6) = mcon_prim(1:total_prim,1:6)    ! lat/lon/radius/area/height/density
  msc_prim_in(1:total_prim,7) = mcon_prim(1:total_prim,11)        ! % land
  msc_prim_in(1:total_prim,8) = mcon_prim(1:total_prim,8)         ! number of ternarys
  msc_prim_in(1:total_prim,9) = mcon_prim(1:total_prim,10)        ! bit-mapped tide flag
! PT170531: add the geoid/ellipsoid separation
  msc_prim_in(1:total_prim,10) = mcon_prim(1:total_prim,12)        ! geoid/ellipsoid separation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      sort through the ternary mascon information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  n_shelf = 0
  n_casp  = 0
  n_eyre  = 0
  call status_update('STATUS','UTIL','separate_mascons',' ','Split off continental shelf ternarys',0)
  do imsc = 1, total_prim
!    if(mod(imsc,500) == 0)then
!      write(message,'(a,i6,a,i6,a,i7,a)')'Reading mascon ',imsc,' of',total_prim,' input primary mascons.' &
!            ,n_shelf,' continental shelf mascons to date.' 
!      call status_update('STATUS','UTIL','separate_mascons',' ',message,0)
!    endif

    tmp_shelf = 0.0
    tmp_other = 0.0
    do itern=1,nint(mcon_prim(imsc,8))

! strip off the continental shelf ternary mascons
! PT161105: this part is now more sophisticated. We separate here the "normal" ternary mascons (ocean or land)
!           from the "Shelf", Caspian Sea, Lake Eyre and ...
      if(mcon_tern(imsc,itern,5) < 0.0 )then

! PT161105: a box around the Caspian sea
        if( 1 .eq.2 .and. mcon_tern(imsc,itern,1) > 36.d0*pi/180. .and. mcon_tern(imsc,itern,1) < 54.d0*pi/180.d0 &
            .and.  mcon_tern(imsc,itern,2) > 44.d0*pi/180.d0 .and.  mcon_tern(imsc,itern,2) < 56.d0*pi/180.d0)then
!          tmp_other = tmp_other + 1
!!print*,'Caspian Sea ternary: ',imsc,itern,mcon_tern(imsc,itern,1:2)*180.d0/pi
!          n_casp  = n_casp  + 1
!          msc_casp(n_casp ,1) = mcon_tern(imsc,itern,1) ! lat
!          msc_casp(n_casp ,2) = mcon_tern(imsc,itern,2) ! lon
!          msc_casp(n_casp ,3) = mcon_tern(imsc,itern,3) ! radius
!          msc_casp(n_casp ,4) = mcon_tern(imsc,itern,4) ! area
!          msc_casp(n_casp ,5) = mcon_tern(imsc,itern,5) ! depth/height 
!          msc_casp(n_casp ,6) = 1000.                   ! density. Over-ride the original to make this fresh water
!          msc_casp(n_casp ,7) = mcon_tern(imsc,itern,7) ! ternary number
!          msc_casp(n_casp ,8) = mcon_tern(imsc,itern,8) ! geoid
!          tern_flag(nint(mcon_tern(imsc,itern,7))) = "TCasp "
          msc_tern_in(imsc,itern-tmp_other,1) = mcon_tern(imsc,itern,1) ! lat
          msc_tern_in(imsc,itern-tmp_other,2) = mcon_tern(imsc,itern,2) ! lon
          msc_tern_in(imsc,itern-tmp_other,3) = mcon_tern(imsc,itern,3) ! radius
          msc_tern_in(imsc,itern-tmp_other,4) = mcon_tern(imsc,itern,4) ! area
          msc_tern_in(imsc,itern-tmp_other,5) = mcon_tern(imsc,itern,5) ! depth/height
          msc_tern_in(imsc,itern-tmp_other,6) = 1000.                   ! density. Over-ride the original to make this land again.
          msc_tern_in(imsc,itern-tmp_other,7) = mcon_tern(imsc,itern,7) ! ternary number
          msc_tern_in(imsc,itern-tmp_other,8) = mcon_tern(imsc,itern,8) ! geoid
          tern_flag(nint(mcon_tern(imsc,itern,7))) = "TLand "
! a few points east of the Caspian Sea that should go back into normal mascons. There are some points out at ~80E ...
        else if (mcon_tern(imsc,itern,1) >  36.d0*pi/180. .and. mcon_tern(imsc,itern,1) <  54.d0*pi/180.d0 &
            .and.  mcon_tern(imsc,itern,2) > 56.d0*pi/180.d0 .and.  mcon_tern(imsc,itern,2) <  100.5d0*pi/180.d0)then
!print*,'East of Caspian Sea ternary: ',imsc,itern,mcon_tern(imsc,itern,1:2)*180.d0/pi
          msc_tern_in(imsc,itern-tmp_other,1) = mcon_tern(imsc,itern,1) ! lat
          msc_tern_in(imsc,itern-tmp_other,2) = mcon_tern(imsc,itern,2) ! lon
          msc_tern_in(imsc,itern-tmp_other,3) = mcon_tern(imsc,itern,3) ! radius
          msc_tern_in(imsc,itern-tmp_other,4) = mcon_tern(imsc,itern,4) ! area
          msc_tern_in(imsc,itern-tmp_other,5) = mcon_tern(imsc,itern,5) ! depth/height
          msc_tern_in(imsc,itern-tmp_other,6) = 1000.                   ! density. Over-ride the original to make this land again.
          msc_tern_in(imsc,itern-tmp_other,7) = mcon_tern(imsc,itern,7) ! ternary number
          msc_tern_in(imsc,itern-tmp_other,8) = mcon_tern(imsc,itern,8) ! geoid
          tern_flag(nint(mcon_tern(imsc,itern,7))) = "TLand "
! Lake Eyre
! PT210902: don't make Lake Eyre a separate for the moment - it is too small to estimate on its own
        else if (mcon_tern(imsc,itern,1) > -30.d0*pi/180. .and. mcon_tern(imsc,itern,1) < -27.d0*pi/180.d0 &
            .and.  mcon_tern(imsc,itern,2) > 136.d0*pi/180.d0 .and.  mcon_tern(imsc,itern,2) < 138.5d0*pi/180.d0)then
!          tmp_other = tmp_other + 1
!!print*,'Lake Eyre ternary: ',imsc,itern,mcon_tern(imsc,itern,1:2)*180.d0/pi,mcon_tern(imsc,itern,5)
!          n_eyre  = n_eyre  + 1
!          msc_eyre(n_eyre ,1) = mcon_tern(imsc,itern,1) ! lat
!          msc_eyre(n_eyre ,2) = mcon_tern(imsc,itern,2) ! lon
!          msc_eyre(n_eyre ,3) = mcon_tern(imsc,itern,3) ! radius
!          msc_eyre(n_eyre ,4) = mcon_tern(imsc,itern,4) ! area
!          msc_eyre(n_eyre ,5) = mcon_tern(imsc,itern,5) ! depth/height 
!          msc_eyre(n_eyre ,6) = 1000.                   ! density. Over-ride the original to make this fresh water
!          msc_eyre(n_eyre ,7) = mcon_tern(imsc,itern,7) ! ternary number
!          msc_eyre(n_eyre ,8) = mcon_tern(imsc,itern,8) ! geoid
!          tern_flag(nint(mcon_tern(imsc,itern,7))) = "TEyre "
          msc_tern_in(imsc,itern-tmp_other,1) = mcon_tern(imsc,itern,1) ! lat
          msc_tern_in(imsc,itern-tmp_other,2) = mcon_tern(imsc,itern,2) ! lon
          msc_tern_in(imsc,itern-tmp_other,3) = mcon_tern(imsc,itern,3) ! radius
          msc_tern_in(imsc,itern-tmp_other,4) = mcon_tern(imsc,itern,4) ! area
          msc_tern_in(imsc,itern-tmp_other,5) = mcon_tern(imsc,itern,5) ! depth/height
          msc_tern_in(imsc,itern-tmp_other,6) = 1000.                   ! density. Over-ride the original to make this land again.
          msc_tern_in(imsc,itern-tmp_other,7) = mcon_tern(imsc,itern,7) ! ternary number
          msc_tern_in(imsc,itern-tmp_other,8) = mcon_tern(imsc,itern,8) ! geoid
          tern_flag(nint(mcon_tern(imsc,itern,7))) = "TLand "
! a few points east of Lake Eyre that should go back into normal mascons
! PT210902: just set all of Lake Eyre back to being land and part of the general mascon pattern (rather than its own mascon)
        else if (mcon_tern(imsc,itern,1) > -30.d0*pi/180. .and. mcon_tern(imsc,itern,1) < -27.d0*pi/180.d0 &
            .and.  mcon_tern(imsc,itern,2) > 138.d0*pi/180.d0 .and.  mcon_tern(imsc,itern,2) < 142.5d0*pi/180.d0)then
!print*,'East of Lake Eyre ternary: ',imsc,itern,mcon_tern(imsc,itern,1:2)*180.d0/pi
          msc_tern_in(imsc,itern-tmp_other,1) = mcon_tern(imsc,itern,1) ! lat
          msc_tern_in(imsc,itern-tmp_other,2) = mcon_tern(imsc,itern,2) ! lon
          msc_tern_in(imsc,itern-tmp_other,3) = mcon_tern(imsc,itern,3) ! radius
          msc_tern_in(imsc,itern-tmp_other,4) = mcon_tern(imsc,itern,4) ! area
          msc_tern_in(imsc,itern-tmp_other,5) = mcon_tern(imsc,itern,5) ! depth/height
          msc_tern_in(imsc,itern-tmp_other,6) = 1000.                   ! density. Over-ride the original to make this land again.
          msc_tern_in(imsc,itern-tmp_other,7) = mcon_tern(imsc,itern,7) ! ternary number
          msc_tern_in(imsc,itern-tmp_other,8) = mcon_tern(imsc,itern,8) ! geoid
          tern_flag(nint(mcon_tern(imsc,itern,7))) = "TLand "
! points below sea level in Egypt that should go back into normal mascons
        else if (mcon_tern(imsc,itern,1) >  28.d0*pi/180. .and. mcon_tern(imsc,itern,1) < 30.66d0*pi/180.d0 &
            .and.  mcon_tern(imsc,itern,2) > 20.6d0*pi/180.d0 .and.  mcon_tern(imsc,itern,2) < 32.15d0*pi/180.d0)then
!print*,'Egyptian ternary: ',imsc,itern,mcon_tern(imsc,itern,1:2)*180.d0/pi
          msc_tern_in(imsc,itern-tmp_other,1) = mcon_tern(imsc,itern,1) ! lat
          msc_tern_in(imsc,itern-tmp_other,2) = mcon_tern(imsc,itern,2) ! lon
          msc_tern_in(imsc,itern-tmp_other,3) = mcon_tern(imsc,itern,3) ! radius
          msc_tern_in(imsc,itern-tmp_other,4) = mcon_tern(imsc,itern,4) ! area
          msc_tern_in(imsc,itern-tmp_other,5) = mcon_tern(imsc,itern,5) ! depth/height
          msc_tern_in(imsc,itern-tmp_other,6) = 1000.                   ! density. Over-ride the original to make this land again.
          msc_tern_in(imsc,itern-tmp_other,7) = mcon_tern(imsc,itern,7) ! ternary number
          msc_tern_in(imsc,itern-tmp_other,8) = mcon_tern(imsc,itern,8) ! geoid
          tern_flag(nint(mcon_tern(imsc,itern,7))) = "TLand "
! points below sea level in California
        else if (mcon_tern(imsc,itern,1) > 32.d0*pi/180. .and. mcon_tern(imsc,itern,1) < 34.d0*pi/180.d0 &
            .and.  mcon_tern(imsc,itern,2) > 243.5d0*pi/180.d0 .and.  mcon_tern(imsc,itern,2) < 245.0d0*pi/180.d0)then
!print*,'A Californian ternary: ',imsc,itern,mcon_tern(imsc,itern,1:2)*180.d0/pi
          msc_tern_in(imsc,itern-tmp_other,1) = mcon_tern(imsc,itern,1) ! lat
          msc_tern_in(imsc,itern-tmp_other,2) = mcon_tern(imsc,itern,2) ! lon
          msc_tern_in(imsc,itern-tmp_other,3) = mcon_tern(imsc,itern,3) ! radius
          msc_tern_in(imsc,itern-tmp_other,4) = mcon_tern(imsc,itern,4) ! area
          msc_tern_in(imsc,itern-tmp_other,5) = mcon_tern(imsc,itern,5) ! depth/height
          msc_tern_in(imsc,itern-tmp_other,6) = 1000.                   ! density. Over-ride the original to make this land again.
          msc_tern_in(imsc,itern-tmp_other,7) = mcon_tern(imsc,itern,7) ! ternary number
          msc_tern_in(imsc,itern-tmp_other,8) = mcon_tern(imsc,itern,8) ! geoid
          tern_flag(nint(mcon_tern(imsc,itern,7))) = "TLand "
! PT200128: a dozen or so points in NW africa below sea level
        else if (mcon_tern(imsc,itern,1) > 30.d0*pi/180. .and. mcon_tern(imsc,itern,1) < 35.d0*pi/180.d0 &
            .and.  mcon_tern(imsc,itern,2) > 5.d0*pi/180.d0 .and.  mcon_tern(imsc,itern,2) < 10.0d0*pi/180.d0)then
!print*,'A Californian ternary: ',imsc,itern,mcon_tern(imsc,itern,1:2)*180.d0/pi
          msc_tern_in(imsc,itern-tmp_other,1) = mcon_tern(imsc,itern,1) ! lat
          msc_tern_in(imsc,itern-tmp_other,2) = mcon_tern(imsc,itern,2) ! lon
          msc_tern_in(imsc,itern-tmp_other,3) = mcon_tern(imsc,itern,3) ! radius
          msc_tern_in(imsc,itern-tmp_other,4) = mcon_tern(imsc,itern,4) ! area
          msc_tern_in(imsc,itern-tmp_other,5) = mcon_tern(imsc,itern,5) ! depth/height
          msc_tern_in(imsc,itern-tmp_other,6) = 1000.                   ! density. Over-ride the original to make this land again.
          msc_tern_in(imsc,itern-tmp_other,7) = mcon_tern(imsc,itern,7) ! ternary number
          msc_tern_in(imsc,itern-tmp_other,8) = mcon_tern(imsc,itern,8) ! geoid
          tern_flag(nint(mcon_tern(imsc,itern,7))) = "TLand "

        else if( mcon_tern(imsc,itern,5) > shelf_depth)then
          tmp_other = tmp_other + 1
          tmp_shelf = tmp_shelf + 1.0
          n_shelf = n_shelf + 1
          if(n_shelf > max_tern)then
            write(message,'(a,f6.1,a,i7,a,i7,a)')"Number of continental shelf ternarys shallower than ",-shelf_depth,"m (" &
                                 ,n_shelf,") exceeds max number of ternarys allowed (",max_tern,"). Increase max_term dimension"
            call status_update('FATAL','UTIL','separate_mascons',' ',message,0)
          endif
          msc_shelf(n_shelf,1) = mcon_tern(imsc,itern,1) ! lat
          msc_shelf(n_shelf,2) = mcon_tern(imsc,itern,2) ! lon
          msc_shelf(n_shelf,3) = mcon_tern(imsc,itern,3) ! radius
          msc_shelf(n_shelf,4) = mcon_tern(imsc,itern,4) ! area
          msc_shelf(n_shelf,5) = mcon_tern(imsc,itern,5) ! depth/height 
          msc_shelf(n_shelf,6) = mcon_tern(imsc,itern,6) ! density
          msc_shelf(n_shelf,7) = mcon_tern(imsc,itern,7) ! ternary number
          msc_shelf(n_shelf,8) = mcon_tern(imsc,itern,8) ! geoid
          tern_flag(nint(mcon_tern(imsc,itern,7))) = "TShelf"
        else
          msc_tern_in(imsc,itern-tmp_other,1) = mcon_tern(imsc,itern,1) ! lat
          msc_tern_in(imsc,itern-tmp_other,2) = mcon_tern(imsc,itern,2) ! lon
          msc_tern_in(imsc,itern-tmp_other,3) = mcon_tern(imsc,itern,3) ! radius
          msc_tern_in(imsc,itern-tmp_other,4) = mcon_tern(imsc,itern,4) ! area
          msc_tern_in(imsc,itern-tmp_other,5) = mcon_tern(imsc,itern,5) ! depth/height
          msc_tern_in(imsc,itern-tmp_other,6) = mcon_tern(imsc,itern,6) ! density
          msc_tern_in(imsc,itern-tmp_other,7) = mcon_tern(imsc,itern,7) ! ternary number
          msc_tern_in(imsc,itern-tmp_other,8) = mcon_tern(imsc,itern,8) ! geoid
          tern_flag(nint(mcon_tern(imsc,itern,7))) = "TDeep "
        endif
      else  ! it is a land ternary
        msc_tern_in(imsc,itern-tmp_other,1) = mcon_tern(imsc,itern,1) ! lat
        msc_tern_in(imsc,itern-tmp_other,2) = mcon_tern(imsc,itern,2) ! lon
        msc_tern_in(imsc,itern-tmp_other,3) = mcon_tern(imsc,itern,3) ! radius
        msc_tern_in(imsc,itern-tmp_other,4) = mcon_tern(imsc,itern,4) ! area
        msc_tern_in(imsc,itern-tmp_other,5) = mcon_tern(imsc,itern,5) ! depth/height
        msc_tern_in(imsc,itern-tmp_other,6) = mcon_tern(imsc,itern,6) ! density
        msc_tern_in(imsc,itern-tmp_other,7) = mcon_tern(imsc,itern,7) ! ternary number
        msc_tern_in(imsc,itern-tmp_other,8) = mcon_tern(imsc,itern,8) ! geoid
        tern_flag(nint(mcon_tern(imsc,itern,7))) = "TLand "
      endif
    enddo
! determine the number of (non-continental shelf/Casp/Eyre) ternary mascons for this primary
    msc_prim_in(imsc,8) = mcon_prim(imsc,8) - tmp_other

! calculate and store the percentage of land of this primary mascon
    debug = .false.
    if(msc_prim_in(imsc,8) > 0)then
      call MASCON_percent_land(debug,nvar_tern,nint(msc_prim_in(imsc,8)),max_tern_per_prim*msc_tern_mult,msc_tern_in(imsc,:,:) &
                              ,msc_prim_in(imsc,7))
      if(debug)print*,'imsc,percent land',imsc, msc_prim_in(imsc,7)
    else
      write(message,'(a,i6,2f8.2,a)')'Mascon ',imsc,msc_prim_in(imsc,1:2)*180./pi &
             ,' has no ternary mascons left. They were all continental shelf/casp/eyre mascons'
      call status_update('STATUS','UTIL','separate_mascons',' ',message,0)
      msc_prim_in(imsc,8) = 0
! PT161111: remove one from the number of current total primary mascons
      current_total_prim = current_total_prim - 1
    endif
    if(msc_prim_in(imsc,7) > 0.5d0 )prim_in_type(imsc) = .true.   ! sets to true if > 50% land ternaries
  enddo
  debug = .false.

! output mascon information
  call status_update('STATUS','UTIL','separate_mascons',' ','have read in the ternary mascons',0)
  write(message,'(i8,a,f8.2,a)')n_shelf,' ternary mascons found in the continental shelf zone (shallower than ',shelf_depth,' m)'
  call status_update('STATUS','UTIL','separate_mascons',' ',message,0)
  write(message,'(i8,a)')n_casp,' ternary mascons found in the Caspian Sea (below sea level)'
  call status_update('STATUS','UTIL','separate_mascons',' ',message,0)
  write(message,'(i8,a)')n_eyre,' ternary mascons found in Lake Eyre (below sea level)'
  call status_update('STATUS','UTIL','separate_mascons',' ',message,0)


! transfer the values of the input primary and ternary mascons to the output arrays
  msc_tern_out = msc_tern_in
  msc_prim_out = msc_prim_in
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       calculate the percentage land of each original primary mascon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do imsc=1, total_prim    
    if(msc_prim_in(imsc,8) > 0)then
      call MASCON_percent_land(debug,nvar_tern,nint(msc_prim_in(imsc,8)),max_tern_per_prim*msc_tern_mult,msc_tern_in(imsc,:,:) &
               ,msc_prim_in(imsc,7))
    endif
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      loop over primary mascons with mixed land/water
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nmixed = 0
!  do imsc=1,total_prim         ! loop over only the original input primary mascons 
  do imsc = 4509,4509
    if(msc_prim_in(imsc,1)*180./pi < 100.d0)cycle ! won't separate any mascons

! PT141030: if changed, update the arrays of the original and the adjusted primary and ternary mascons after each loop and recompute the percentage of land for each
    do j=1,total_prim
      if(msc_prim_in(j,8) /= msc_prim_out(j,8))then   
        msc_prim_in(j,:) = msc_prim_out(j,:)
        msc_tern_in(j,:,:) = msc_tern_out(j,:,:)
        debug = .false.
        if(msc_prim_in(j,8) > 0)then
          call MASCON_percent_land(debug,nvar_tern,nint(msc_prim_in(j,8)),max_tern_per_prim*msc_tern_mult,msc_tern_in(j,:,:) &
                ,msc_prim_in(j,7))
        endif
      endif
    enddo


    debug = .false.
    if(debug)print*,'imsc,msc_prim_in(imsc,7)',imsc,msc_prim_in(imsc,1:2)*180./pi,msc_prim_in(imsc,7)


    if(msc_prim_in(imsc,7) > 0.d0 .and. msc_prim_in(imsc,7) < 1.d0)then
      write(message,'(a,i6,a,f6.2,a)')"Mixed primary mascon:",imsc,'. Percent land is',msc_prim_in(imsc,7)*100.d0,'%'
      call status_update('STATUS','UTIL','separate_mascons',' ',message,0)
      nmixed = nmixed + 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! We now have the coords of the surrounding primary mascons. Loop over all misfit ternary 
! mascons to work out to which of the surrounding mascons the misfits should be moved.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

!!! *********************************************************************************************************************
!!! PT210902: change the logic here. We want to eject all ocean parts of a mixed mascon AND all land parts so
!!!           that the mascon no longer exists at all. Change the "if land then ... else if ocean then ..." to 
!!!           just doing both cases. This should empty the mixed primary mascon off all ternaries.
!!! *********************************************************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                    !
!         eject ocean mascons        !
!                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     debug = .false.
     call MASCON_nearest_mscs(debug,imsc,max_surround_msc,total_prim,msc_tern_mult,max_tern_per_prim,nvar_tern &
                             ,nvar_prim,msc_prim_out,msc_tern_out,max_distance,msc_surrounding &
                              ,1029.d0,n_surround)
     if(msc_prim_in(imsc,7) > 0.2d0  )then          ! Primary is < 90% ocean. Therefore, we want to eject ocean ternary mascons
        ntern_given = 0
        ntern_ocean = 0
        do itern = 1, nint(msc_prim_in(imsc,8)  ) 
          if(debug)print*,'top of tern loop for giving ocean ternaries',imsc,itern,' of',nint(msc_prim_in(imsc,8)  )

!          if(msc_tern_in(imsc,itern,5) < 0.d0) then  ! height value is < 0, so it is an ocean ternary. 
          if(tern_flag(nint(msc_tern_in(imsc,itern,7))) == "TDeep " ) then  ! ternary flag is for ocean, so it is an ocean ternary. 
            ntern_ocean = ntern_ocean+1
            
! Is there a neighbouring mascon with the appropriate mix to take the unwanted ternary mascons?
            transfer_tern = .false.
            min_dist = 1.d7  ! set default distances ternary to surrounding primary mascons to large values.
            give_to  = 0
            do j=1,n_surround
              debug = .false.
              if(msc_surrounding(j,4) < 0.5d0 .and.msc_surrounding(j,3) > 0 ) then   ! it is predominantly water. Accept the transfer of the ternary water mascon.
                if(debug)print*,'mascon ',j,' might take unwanted ternary ocean mascons  ',msc_surrounding(j,:)

! calculate the distance from the ternary to the coords of this primary mascon
                call MASCON_msc_to_msc(msc_tern_in(imsc,itern,1),msc_tern_in(imsc,itern,2),msc_surrounding(j,1) &
                               ,msc_surrounding(j,2), tmp_dist)

!               is this surrounding mascon closer than previous ones?
                if(tmp_dist < min_dist)then
                  min_dist = tmp_dist
                  give_to = msc_surrounding(j,3)
                  transfer_tern = .true.
                  if(debug)print*,'mascon ',j,' set to take unwanted ternary ocean mascons ',msc_surrounding(j,:)
                endif

              endif  ! end of loop 
            enddo  ! end of loop over surrounding primary mascons
            if(transfer_tern)then
              ntern_given = ntern_given + 1
              debug = .false.
              if(debug)print*,'Primary ',give_to,' taking ocean ternary',itern,msc_tern_in(imsc,itern,1:2)
              call MASCON_transfer_ternary(debug,imsc,give_to,total_prim*msc_tern_mult,max_tern_per_prim*msc_tern_mult &
                                      ,nvar_prim,nvar_tern,itern &
                                      ,msc_tern_in(imsc,:,:),msc_tern_out(imsc,:,:),msc_tern_out(give_to,:,:) &
                                      ,msc_prim_out(imsc,:),msc_prim_out(give_to,:) ,nint(msc_prim_in(imsc,8)))
            endif

          endif  ! end of test whether ternary is an ocean mascon
        enddo  ! end of ternary loop 

        write(message,'(a,i4,a,i4,a,i5)')'Transferred ',ntern_given,' ocean ternary mascons of',nint(msc_prim_in(imsc,8)  ) &
               ,' ternary mascons out of mixed primary mascon ',imsc
        call status_update('STATUS','UTIL','separate_mascons',' ',message,0)
! check that the area is now 100% land
        debug = .false.
        if(msc_prim_out(imsc,8) > 0)then
          call MASCON_percent_land(debug,nvar_tern,nint(msc_prim_out(imsc,8)) &
                                          ,max_tern_per_prim*msc_tern_mult,msc_tern_out(imsc,:,:),percent_land)
        endif
! PT141110: we must update the percent land of the output primary mascon
        msc_prim_out(imsc,7) = percent_land
        if(percent_land < 1.d0)then
          write(message,'(a,i6,a,f6.2,a,2f10.4)')' Land primary mascon (',imsc,') still mixed (' &
                   ,(1.d0-percent_land)*100.d0,'% ocean)' ,msc_prim_out(imsc,1:2)*180./pi
          call status_update('WARNING','UTIL','separate_mascons',' ',message,0)
        else
          write(message,'(a,i6,a,f6.2,a,2f10.4,i5,a)')'Mixed primary mascon (',imsc,') now ' &
               ,percent_land*100.d0,'% land',msc_prim_out(imsc,1:2)*180./pi,nint(msc_prim_out(imsc,8)  ),' ternaries left.'
!          call status_update('STATUS','UTIL','separate_mascons',' ',message,0)          
        endif
      endif  ! end of if statement checking whether to eject ocean ternaries

      ! PT210902: recalculate the centre of mass of the primary
      !call MASCON_calc_CoM(max_tern_per_prim*msc_tern_mult ,nvar_prim,nvar_tern,msc_tern_out(imsc,:,:),msc_prim_out(imsc,:) )

      ! PT210902: need to update msc_prim_in for this mascon to be what is left after removing ocean ternaries
      msc_prim_in(imsc,:)   = msc_prim_out(imsc,:)
      msc_tern_in(imsc,:,:) = msc_tern_out(imsc,:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                    !
!        eject  land mascons         !
!                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (msc_prim_in(imsc,7) < 0.8d0 .or. nint(msc_prim_in(imsc,8)) < min_terns) then  ! the mixed mascon is < 90% land or too small. Need to eject the land ternary mascons  
!        if(nint(msc_prim_in(imsc,8)) < min_terns) then
!          print*,'**** mascon',imsc,' has only', nint(msc_prim_in(imsc,8)),' ternaries ****'
!        endif

! PT210903: I don't want the land primaries within Iceland to transferred to Greenland primaries
! for 300 km mascons_stage2 it is mascon 222
! for 200 km mascons_stage2_200km it is 479 and 557
        if(total_prim == 4582 .and. (imsc == 221 .or. imsc ==222 .or. imsc==223 .or. imsc==275 .or. imsc==276 &
          .or. imsc==41 ) )then      ! it is likely that this is mascons_stage2_300km
          cycle
        elseif(total_prim == 10314 .and. ((imsc>=478.and.imsc<=480).or.(imsc>=557.and.imsc<=559)) )then ! it is likey that this is mascons_stage2_200km
          cycle
        endif


        
        debug = .false.
        ntern_given = 0
        ntern_land = 0

        ! recalculate the closest mascons, ignoring any that may now be empty of ternaries
        call MASCON_nearest_mscs(debug,imsc,max_surround_msc,total_prim,msc_tern_mult,max_tern_per_prim,nvar_tern &
                              ,nvar_prim,msc_prim_out,msc_tern_out,max_distance,msc_surrounding &
                              ,1000.d0,n_surround)


        do itern = 1, nint(msc_prim_in(imsc,8)  )
          if(debug)print*,'ternary record:',imsc,itern,msc_tern_in(imsc,itern,:)
          if(msc_tern_in(imsc,itern,5) >= 0.d0) then  ! it is a ternary with at least some land. 
            ntern_land = ntern_land+1
! Is there a neighbouring mascon with the appropriate mix to take the unwanted ternary mascons?
            transfer_tern = .false.
            min_dist = 1.d10  ! set default distances ternary to surrounding primary mascons to large values.
            give_to  = 0
            do j=1,n_surround
              if(msc_surrounding(j,4) >= 0.d0) then   ! it has land. Accept the transfer of the ternary land mascon.
                transfer_tern = .true.

! calculate the distance from the ternary to the coords of this primary mascon
                call MASCON_msc_to_msc(msc_tern_in(imsc,itern,1),msc_tern_in(imsc,itern,2),msc_surrounding(j,1) &
                               ,msc_surrounding(j,2), tmp_dist)
                debug=.false.
                if(debug)print*,'imsc,j,tmp_dist',imsc,j,tmp_dist
                if(tmp_dist < min_dist)then
                  min_dist = tmp_dist
                  give_to = msc_surrounding(j,3)
                  if(debug)print*,'surround mascon ',j,' set to take unwanted ternary land  mascons ',msc_surrounding(j,:) &
                                 ,min_dist
                endif

              else
!print*,'surrounding mascon',j,msc_surrounding(j,:)
              endif
            enddo  ! end of loop over surrounding primary mascons
            if(transfer_tern)then
              ntern_given = ntern_given + 1
              if(imsc == 89)debug = .true.
              if(debug)print*,'Primary ',give_to,' taking land  ternary',itern,msc_tern_in(imsc,itern,1:2) &
                              ,(msc_prim_in(give_to,7)),nint(msc_prim_in(give_to,8))
              debug = .false.
              call MASCON_transfer_ternary(debug,imsc,give_to,total_prim*msc_tern_mult,max_tern_per_prim*msc_tern_mult &
                                      ,nvar_prim,nvar_tern &
                                      ,itern,msc_tern_in(imsc,:,:) &
                                      ,msc_tern_out(imsc,:,:),msc_tern_out(give_to,:,:),msc_prim_out(imsc,:) &
                                      ,msc_prim_out(give_to,:),nint(msc_prim_in(imsc,8)))
            else
              write(message,'(a,i7,a,i7)')"no surrounding mascons with enough land to transfer ternary",itern &
                                      ," from mascon",imsc
              call status_update('STATUS','UTIL','separate_mascons',' ',message,0)
            endif
          endif  ! end of test whether ternary is a land mascon
        enddo  ! end of ternary loop 


        if(ntern_given > 0)then
          write(message,'(a,i4,a,i5)')'Transferred ',ntern_given,'  land ternary mascons out of mixed primary mascon ',imsc
        else
          write(message,'(a,i7)')'Could not transfer any land ternary mascons out of ocean primary mascon ',imsc
        endif
        call status_update('STATUS','UTIL','separate_mascons',' ',message,0)

! recalculate the percentage of land
        debug = .false.
        if(msc_prim_out(imsc,8) > 0)call MASCON_percent_land(debug,nvar_tern,nint(msc_prim_out(imsc,8)) &
                                            ,max_tern_per_prim*msc_tern_mult,msc_tern_out(imsc,:,:),percent_land)
! PT141110: we must update the percent land of the output primary mascon
        msc_prim_out(imsc,7) = percent_land
        if (nint(msc_prim_out(imsc,8)) == 0 )then
          write(message,'(a,i6)')'All ternaries transferred out of mixed primary mascon: ',imsc
          call status_update('STATUS','UTIL','separate_mascons',' ',message,0)    
          msc_prim_out(imsc,7) = 0.d0     
        else
          write(message,'(a,i6,a,f8.2,a,2f10.4)')'Mixed primary mascon (',imsc,') changed to ',percent_land*100.d0,'% land' &
                 ,msc_prim_out(imsc,1:2)*180./pi
          call status_update('STATUS','UTIL','separate_mascons',' ',message,0)          
        endif
      endif  ! end of test on whether to eject land ternaries  


    else
      write(message,'(a,i6,a,f6.2,a,i6,a,2f10.4)')'Remaining  primary mascon (',imsc,') unmixed at ',percent_land*100.d0 &
                 ,"% land with",nint(msc_prim_out(imsc,8))," ternaries.",msc_prim_out(imsc,1:2)*180./pi
          call status_update('STATUS','UTIL','separate_mascons',' ',message,0)            
    endif  ! end of mixed mascon if statement

  enddo    ! end of input primary mascon loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                               !
!       Calculate new centre of mass of new Primary mascons             !
!                                                                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call status_update('STATUS','UTIL','separate_mascons',' ','Calculating new centroids of primary mascons',0)
! PT161104: also work out the maximum number of ternary mascons per primary while going through this loop
  new_max_tern_per_prim = 0
  do imsc=1,total_prim
! PT141111: now also calculates the area of the new primary mascons
    call MASCON_calc_CoM(max_tern_per_prim*msc_tern_mult ,nvar_prim,nvar_tern,msc_tern_out(imsc,:,:),msc_prim_out(imsc,:) )
    if(nint(msc_prim_out(imsc,8)) > new_max_tern_per_prim)then
      new_max_tern_per_prim = nint(msc_prim_out(imsc,8))
    endif
!print*,'CoMs imsc,msc_prim_out(21,8)',imsc,msc_prim_out(21,8)
  enddo 
! PT161105: check whether the number of shelf ternary mascons is greater than this number
  if(n_shelf > new_max_tern_per_prim)new_max_tern_per_prim = n_shelf
  if(n_casp  > new_max_tern_per_prim)new_max_tern_per_prim = n_casp 
  max_tern_per_sec = new_max_tern_per_prim 

 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!                                    OUTPUT
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call status_update('STATUS','UTIL','separate_mascons',msc_output,"Outputting ocean/land separated mascon file",0)
  new_prim_total = 0

! update the numbers in the first header line
  if(n_shelf > 0)then
! PT170601: we need to add 4 here, not sure why it isn't three!!!
    max_prim = current_total_prim + 4   ! add one each for Shelf, Capsian and Lake Eyre mascons
    max_sec = max_prim
  else
    max_prim = current_total_prim + 3   ! add one each for Caspian and Lake Eyre mascons
    max_sec = max_prim
  endif
  write(msc_hdr_lines(1),'(a1,a,6i8)')"#",msc_hdr_code,max_prim,max_sec,max_tern &
                                ,max_sec_per_prim,new_max_tern_per_prim,max_tern_per_sec
! write out the header of the file. Not done properly for the moment .....
  write(msc_hdr_lines(n_hdr_lines+1),'(2a,f8.3,a)')"# modified by separate_mascons to remove mixed ocean/land primary mascons. "&
        ," Continental shelf depth: ",shelf_depth
  do j=1,n_hdr_lines+1
    write(lumsc_out,'(a)')msc_hdr_lines(j)
  enddo

! set up our random numbers for the colours of the primary mascons
  call srand(total_prim)

! now the mascons
  new_prim_total = 0
  do imsc=1,total_prim
!print*,'imsc,msc_prim_out(21,8)',imsc,msc_prim_out(21,8),msc_prim_out(imsc,8)

     if(mod(imsc,500) == 0)then
       write(message,'(a,i6,2f8.2)')'   Writing mascon ',imsc, msc_prim_out(imsc,1:2)*180.d0/pi
       call status_update('STATUS','UTIL','separate_mascons',' ',message,0)
     endif
    debug = .false.
! PT141116: don't call this if there are no ternary mascons
    if(msc_prim_out(imsc,8) > 0)then
      call MASCON_percent_land(debug,nvar_tern,nint(msc_prim_out(imsc,8)),max_tern_per_prim*msc_tern_mult,msc_tern_out(imsc,:,:) &
                ,msc_prim_out(imsc,7))
    else
      msc_prim_out(imsc,7) = 0
    endif
! assign a density based on whether it is land or water
    if(msc_prim_out(imsc,7) < 0.5) then
      msc_prim_out(imsc,6) = 1029.d0
      msc_prim_out(imsc,9) = 0
      prim_flags(imsc) = "PDeep "
!print*,'should be water',imsc
    else
!print*,'should be land',imsc
      msc_prim_out(imsc,6) = 1000.d0
! based on mean depth, decide whether to turn on the estimation of tidal constituents. 
      prim_flags(imsc) = "PLand "
      msc_prim_out(imsc,9) = 0
    endif
! PT161103: only write out the primary/secondary/ternary if there are more than zero ternary mascons left in this primary
    if(msc_prim_out(imsc,8) > 0)then
      new_prim_total = new_prim_total + 1
      nsec = 1

! PT161103: assign a colour to this primary mascon
      if(msc_prim_out(imsc,5) < 0.)then
        prim_colour = rand(0)*450.
        prim_colour = rand(0)*450.
      else if(msc_prim_out(imsc,5) >= 0.)then
        prim_colour = 720 + rand(0)*780.
        prim_colour = 720 + rand(0)*780.
      endif

! write out the primary mascon
      msc_code = prim_flags(imsc)
      msc_crds(1:2) = msc_prim_out(imsc,1:2)*180.d0/pi
      msc_crds(3) = msc_prim_out(imsc,3)
      pointers = 0
      colours = 0
      region = mcon_region(imsc)
      call write_mascon_record(lumsc_out,msc_code,new_prim_total,nsec,nint(msc_prim_out(imsc,8)),msc_crds &
            ,msc_prim_out(imsc,4),msc_prim_out(imsc,5),msc_prim_out(imsc,10),msc_prim_out(imsc,6),msc_prim_out(imsc,7)*100.d0 &
            ,msc_prim_out(imsc,9) , region,pointers,colours)
!if(imsc >10 .and. imsc < 30)print*,"wrote out ",msc_code,imsc,new_prim_total,msc_prim_out(imsc,6),msc_prim_out(imsc,7)*100.d0,nint(msc_prim_out(imsc,8))

! PT161030: also write out a secondary mascon line
      pointers(1) = new_prim_total
      pointers(2) = 0
      colours = prim_colour
      msc_code = "S"//prim_flags(imsc)(2:6)
      call write_mascon_record(lumsc_out,msc_code,new_prim_total,nsec,nint(msc_prim_in(imsc,8)),msc_crds &
            ,msc_prim_out(imsc,4),msc_prim_out(imsc,5),msc_prim_out(imsc,10),msc_prim_out(imsc,6),msc_prim_out(imsc,7)*100.d0 &
            ,msc_prim_out(imsc,9) , region,pointers,colours)

! and write out the ternary lines
      do itern=1,nint(msc_prim_out(imsc,8))
! PT161205: do we need to check the density here before assigning the TDeep or TLand etc flag?
        if(msc_tern_out(imsc,itern,6) > 1000.d0)then
          msc_code = "TDeep "   !//prim_flags(imsc)(2:6)
        else
          msc_code = "TLand "   !//prim_flags(imsc)(2:6)
        endif
        msc_crds(1:2) = msc_tern_out(imsc,itern,1:2)*180.d0/pi
        msc_crds(3) = msc_tern_out(imsc,itern,3)
        pointers = new_prim_total
        colours(1) = prim_colour
        colours(2) = prim_colour
        call write_mascon_record(lumsc_out,msc_code,nint(msc_tern_out(imsc,itern,7)),nsec,nint(msc_prim_in(imsc,8)),msc_crds &
            ,msc_tern_out(imsc,itern,4),msc_tern_out(imsc,itern,5),msc_tern_out(imsc,itern,8),msc_tern_out(imsc,itern,6) &
            ,msc_prim_out(imsc,7)*100.d0,msc_prim_out(imsc,9), region,pointers,colours)
      enddo
    else
print*,'not outputting original mascon',imsc
    endif

  enddo  ! end of loop over writing out primary mascons.
  write(message,'(a,i7,a)')"Written to file ",new_prim_total," separated primary mascons"
  call status_update('STATUS','UTIL','separate_mascons',msc_output,message,0)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT161105: output the Caspian Sea ternary mascons as one primary
  if(n_casp > 0)then
    call MASCON_calc_CoM(max_seas,nvar_prim,nvar_tern,msc_casp(:,:),msc_prim_out(1,:) )  ! called here just to compute the centroid of the Caspian Sea ternary mascons
    write(message,'(a,i7,a)')"Writing to file ",n_casp ," Caspian Sea ternary mascons"
    call status_update('STATUS','UTIL','separate_mascons',msc_output,message,0)
    prim_flag = "Casp "
    prim_colour = 200
    tmp_sec_colour = prim_colour

! primary Caspian
    msc_code = "P"//prim_flag
    msc_crds(1:2) = msc_prim_out(1,1:2)*180.d0/pi
    msc_crds(3) = msc_prim_out(1,3)
    pointers = new_prim_total+1
    colours = 200
    region = "Caspian        "
    density = 1000.d0
    msc_prim_out(1,9) = 0     ! don't estimate tides in the Caspian Sea (not yet, at least)
    call write_mascon_record(lumsc_out,msc_code,new_prim_total+1,nsec,n_casp,msc_crds &
            ,msc_prim_in(1,4),msc_prim_in(1,5),msc_prim_in(1,10),density,msc_prim_in(1,7)*100.d0,msc_prim_out(1,9) &
            , region,pointers,colours)

! secondary Caspian
      pointers(1) = new_prim_total+1
      pointers(2) = 0
      colours = prim_colour
      msc_code = "S"//prim_flag
      call write_mascon_record(lumsc_out,msc_code,new_prim_total+1,nsec,n_casp,msc_crds &
            ,msc_prim_in(1,4),msc_prim_in(1,5),msc_prim_in(1,10),density,msc_prim_in(1,7)*100.d0,msc_prim_out(1,9) &
            , region,pointers,colours)


! ternary Caspian
      do itern=1,n_casp 
        msc_code = "T"//prim_flag
        msc_crds(1:2) = msc_casp(itern,1:2)*180.d0/pi
        msc_crds(3) = msc_casp(itern,3)
        pointers = new_prim_total + 1
        colours(1) = prim_colour
        colours(2) = tmp_sec_colour
        call write_mascon_record(lumsc_out,msc_code,nint(msc_casp(itern,7)),nsec,nint(msc_prim_in(1,8)),msc_crds &
            ,msc_casp(itern,4),msc_casp(itern,5),msc_casp(itern,8),density,msc_prim_in(1,7)*100.d0 &
            ,msc_prim_out(1,9), region,pointers,colours)
      enddo
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT161105: output the Lake Eyre ternary mascons as one primary
    if(n_eyre > 0)then
      call MASCON_calc_CoM(max_seas,nvar_prim,nvar_tern,msc_eyre(:,:),msc_prim_out(1,:) )  ! called here just to compute the centroid of the Caspian Sea ternary mascons
      write(message,'(a,i7,a)')"Writing to file ",n_eyre ," Lake Eyre ternary mascons"
      call status_update('STATUS','UTIL','separate_mascons',msc_output,message,0)
      prim_flag = "Eyre "
      prim_colour = 200
      tmp_sec_colour = prim_colour
! primary Lake Eyre
      msc_code = "P"//prim_flag
      msc_crds(1:2) = msc_prim_out(1,1:2)*180.d0/pi
      msc_crds(3) = msc_prim_out(1,3)
      pointers = new_prim_total+2
      colours = 200
    region = "Lake_Eyre      "
      density = 1000.d0
      msc_prim_out(1,9) = 0     ! don't estimate tides in Lake Eyre 
      call write_mascon_record(lumsc_out,msc_code,new_prim_total+2,nsec,n_eyre,msc_crds &
            ,msc_prim_in(1,4),msc_prim_in(1,5),msc_prim_in(1,10),density,msc_prim_in(1,7)*100.d0,msc_prim_out(1,9) &
            , region,pointers,colours)

! secondary Lake Eyre
      pointers(1) = new_prim_total+2
      pointers(2) = 0
      colours = prim_colour
      msc_code = "S"//prim_flag
      call write_mascon_record(lumsc_out,msc_code,new_prim_total+2,nsec,n_eyre,msc_crds &
            ,msc_prim_in(1,4),msc_prim_in(1,5),msc_prim_in(1,10),density,msc_prim_in(1,7)*100.d0,msc_prim_out(1,9) &
            , region,pointers,colours)

! ternary Lake Eyre
      do itern=1,n_eyre
        msc_code = "T"//prim_flag
        msc_crds(1:2) = msc_eyre(itern,1:2)*180.d0/pi
        msc_crds(3) = msc_eyre(itern,3)
        pointers = new_prim_total + 2
        colours(1) = prim_colour
        colours(2) = tmp_sec_colour
        call write_mascon_record(lumsc_out,msc_code,nint(msc_eyre(itern,7)),nsec,nint(msc_prim_in(1,8)),msc_crds &
            ,msc_eyre(itern,4),msc_eyre(itern,5),msc_eyre(itern,8),density,msc_prim_in(1,7)*100.d0 &
            ,msc_prim_out(1,9), region,pointers,colours)
      enddo
    endif

! PT141117: finally, output the continental shelf ternary mascons
  if(n_shelf > 0) then
    call MASCON_calc_CoM(max_shelf,nvar_prim,nvar_tern,msc_shelf(:,:),msc_prim_out(1,:) )  ! called here just to compute the average depth of the continental shelf ternary mascons
    write(message,'(a,i7,a)')"Writing to file ",n_shelf," continental shelf ternary mascons"
    call status_update('STATUS','UTIL','separate_mascons',msc_output,message,0)
    prim_flag = "Shelf"
    nsec = 1
! define the colour for the shelf mascons
    prim_colour = 520 + rand(0)*170.
    prim_colour = 520 + rand(0)*170.
    tmp_sec_colour = prim_colour

    density = 1029.d0
! primary Shelf mascons
    msc_code = "PShelf"      
    msc_crds(1:2) = msc_prim_out(1,1:2)*180.d0/pi
    msc_crds(3) = msc_prim_out(1,3)
    pointers = new_prim_total+3
    colours = 200
    region = "Shelf          "
    msc_prim_out(1,9) = 31     ! estimate tides on the continental shelves and shallow seas 
    call write_mascon_record(lumsc_out,msc_code,new_prim_total+3,nsec,n_shelf,msc_crds &
            ,msc_prim_in(1,4),msc_prim_in(1,5),msc_prim_in(1,10),density,msc_prim_in(1,7)*100.d0,msc_prim_out(1,9) &
            , region,pointers,colours)

! secondary Shelf
    pointers(1) = new_prim_total+3
    pointers(2) = 0
    colours = prim_colour
    msc_code = "S"//prim_flag
    call write_mascon_record(lumsc_out,msc_code,new_prim_total+3,nsec,n_shelf,msc_crds &
            ,msc_prim_in(1,4),msc_prim_in(1,5),msc_prim_in(1,10),density,msc_prim_in(1,7)*100.d0,msc_prim_out(1,9) &
            , region,pointers,colours)

! ternary Shelf mascons
    do itern=1,n_shelf
      msc_code = "T"//prim_flag
      msc_crds(1:2) = msc_shelf(itern,1:2)*180.d0/pi
      msc_crds(3) = msc_shelf(itern,3)
      pointers = new_prim_total + 3
      colours(1) = prim_colour
      colours(2) = tmp_sec_colour
      call write_mascon_record(lumsc_out,msc_code,nint(msc_shelf(itern,7)),nsec,nint(msc_prim_in(1,8)),msc_crds &
            ,msc_shelf(itern,4),msc_shelf(itern,5),msc_shelf(itern,8),msc_shelf(itern,6),msc_prim_in(1,7)*100.d0 &
            ,msc_prim_out(1,9), region,pointers,colours)

    enddo
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  call status_update('STATUS','UTIL','separate_mascons',' ',"End of separate_mascons",0)



  end

