! allocate_mascon_arrays
! deallocate_mascon_arrays
! read_msc_hdr
! ! tern_lat_bands
! calc_which_ternary
! read_mascon_file
! read_ocean_mascons
! write_mascon_record
! calc_mascon_xyz
! msc_in_region
! get_random_value
! tern_lat_bands_ell
! adjust_mascon_surface
! read_msc_model

  subroutine allocate_mascon_arrays

! subroutine to allocate the mascon arrays. It uses different memory to dimension them dynamically
! than to dimension them in the mascon_mod.f90 file, so I'm doing it this way. It also speeds up
! the compilation time for the programs .....
!
! P. Tregoning
! 9 September 2016
!
! MODS

  use mascon_mod

  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                           !!
!!   A L L O C A T E    T H E    M A S C O N    A R R A Y S  !!
!!                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! primary mascon variables
! PT161028: added tidal amplitude flag as an extra column to mcon_prim
! PT170530: added geoid as an extra column (increased to 12)
  allocate(mcon_prim(max_prim,nvar_prim))    ! primary mascon info (lat/lon/rad/area/hgt/density/#sec/#tern/#first secondary/tidal flag/%land/geoid)
  allocate(mcon_prim_EWH(max_prim))        ! EWH a priori value for each mascon
  allocate(mcon_region(max_prim))     ! region character code for each primary mascon
  allocate(prim_flags(max_prim))      ! 6-char flags for type of mascon (PDeep/PShelf etc)

  ! array used in graceorb for all the mix of prim/sec/tern mascons to use in the orbit integration
  allocate(mcon_EWH_vector(max_tern))

! secondary
! PT170530: added geoid as an extra column (increased to 9)
  allocate(mcon_sec(max_prim,max_sec_per_prim,nvar_prim))    ! secondary info ( lat/lon/rad/area/hgt/density/sec_number/#terns in this secondary/geoid )
  allocate(mcon_sec_EWH(max_sec))                   ! EWH a priori value for each mascon
  allocate(mcon_sec_ptr(max_sec))                    ! pointer to primary mascon in which this secondary resides
  allocate(sec_colour(max_sec))                      ! colour of primary mascon in which this secondary resides
  allocate(sec_flags(max_sec))                       ! 6-char flags for type of mascon (PDeep/PShelf etc)

! ternary
! PT161028: added tidal amplitude flag as an extra column to mcon_tern
! PT170530: added geoid as an extra column (increased to 9)
  allocate(mcon_tern(max_prim,max_tern_per_prim,nvar_tern))  ! ternary info ( lat/lon/rad/area/hgt/density/tern_number/geoid ) 
  allocate(mcon_tern_EWH(max_tern))                  ! EWH a priori value for each ternary mascon
  allocate(mcon_tern_ptr(max_tern,3))                ! pointers to primary/secondary/seq. ternary within primary/
  allocate(tern_colours(max_tern,2))                 ! colour of primary and secondary mascon in which this ternary resides
  allocate(tern_flag(max_tern))                      ! 6-char flags for type of ternary mascon (PDeep/PShelf etc)

! the arrays used to compute the accelerations from the vector list of mascons per epoch
  allocate(mcon_xyz(3,800 * max_prim))    ! XYZ Efixed coordinates (m)
  allocate(mcon_area(800 * max_prim))     ! mascon areas (sq. m.)
  allocate(mcon_num(max_prim))            ! mass con counter for gravity & partials
  allocate(mcon_rho(800* max_prim))       ! density of the mascons in the list 

!!!!!!!!!!!!!!!!!!!!
!   OCEAN MASCONS  !
!!!!!!!!!!!!!!!!!!!!
! primary ocean mascons
  allocate(mcon_ocean_prim(max_ocean_prim,2))         ! primary ocean info (lat/lon/rad/area/hgt/#tern/bit-mapped tide
  allocate(prim_flags_ocean(max_ocean_prim))          ! 6-char flags for type of ocean mascon (PDeep/PShelf etc)
! ternary
! PT170530: added geoid as an extra column (increased to 9)
  allocate(mcon_ocean_tern(max_ocean_tern,9))       ! ocean ternary info (lat/lon/rad/area/hgt/density/tern2prim pointer/the ternary number/geoid
! PT190522: mcon_ocean_tern_ptr needs to be allocated according to the total ternarys, NOT the total ocean ternarys
  allocate(mcon_ocean_tern_ptr(max_tern,2))         ! pointer to ocean primary mascons in which the ocean ternary resides (needs to be dimensioned for
  allocate(tern_flag_ocean(max_ocean_tern))         ! 6-char flags for type of ternary mascon (PDeep/PShelf etc)
                                                    !    the maximum number of ALL ternarys, not just ocean ternarys

! set the pointer values to zero
  mcon_ocean_tern_ptr = 0

!!!!!!!!!!!!!!!!!!!!
! tidal amplitudes !
!!!!!!!!!!!!!!!!!!!!
! tidal amplitude partials
  allocate(ocean_eftid_part(3,max_msc_tides,2,max_ocean_prim))    ! partials for tidal amplitude parameters
  allocate(mcon_tide_name(max_msc_tides))                         ! 2-char descriptor of tidal constituents (M2, O1, S2, K1, K2). These are defined in graceorb.f90
  allocate(mcon_tide_period(max_msc_tides))                       ! period of the tidal constituents
  allocate(msc_tide(max_msc_tides,2,max_ocean_prim))              ! a priori values for the amplitudes of each tidal constituent (sine and cosine components)

  return
  end  subroutine allocate_mascon_arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  subroutine deallocate_mascon_arrays

! subroutine to deallocate the mascon arrays.
!
! P. Tregoning
! 5 December 2018
!
! MODS

  use mascon_mod

  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                           !!
!!   DE-A L L O C A T E    T H E    M A S C O N    A R R A Y S  !!
!!                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! primary mascon variables
! PT161028: added tidal amplitude flag as an extra column to mcon_prim
! PT170530: added geoid as an extra column (increased to 12)
  deallocate(mcon_prim)    ! primary mascon info (lat/lon/rad/area/hgt/density/#sec/#tern/#first secondary/tidal flag/%land/geoid)
  deallocate(mcon_prim_EWH)        ! EWH a priori value for each mascon
  deallocate(mcon_region)     ! region character code for each primary mascon
  deallocate(prim_flags)      ! 6-char flags for type of mascon (PDeep/PShelf etc)

! secondary
! PT170530: added geoid as an extra column (increased to 9)
  deallocate(mcon_sec)    ! secondary info ( lat/lon/rad/area/hgt/density/sec_number/#terns in this secondary/geoid )
  deallocate(mcon_sec_ptr)                    ! pointer to primary mascon in which this secondary resides
  deallocate(sec_colour)                      ! colour of primary mascon in which this secondary resides
  deallocate(sec_flags)                       ! 6-char flags for type of mascon (PDeep/PShelf etc)

! ternary
! PT161028: added tidal amplitude flag as an extra column to mcon_tern
! PT170530: added geoid as an extra column (increased to 9)
  deallocate(mcon_tern)  ! ternary info ( lat/lon/rad/area/hgt/density/tern_number/geoid ) 
  deallocate(mcon_tern_ptr)                ! pointers to primary/secondary/seq. ternary within primary/
  deallocate(tern_colours)                 ! colour of primary and secondary mascon in which this ternary resides
  deallocate(tern_flag)                      ! 6-char flags for type of ternary mascon (PDeep/PShelf etc)

! the arrays used to compute the accelerations from the vector list of mascons per epoch
  deallocate(mcon_xyz)    ! XYZ Efixed coordinates (m)
  deallocate(mcon_area)     ! mascon areas (sq. m.)
  deallocate(mcon_num)            ! mass con counter for gravity & partials
  deallocate(mcon_rho)       ! density of the mascons in the list 

! the mascon EWH vectors
  deallocate(mcon_EWH_vector)
  deallocate(mcon_sec_EWH)
  deallocate(mcon_tern_EWH)
  
!!!!!!!!!!!!!!!!!!!!
!   OCEAN MASCONS  !
!!!!!!!!!!!!!!!!!!!!
! primary ocean mascons
  deallocate(mcon_ocean_prim)         ! primary ocean info (lat/lon/rad/area/hgt/#tern/bit-mapped tide
  deallocate(prim_flags_ocean)        ! 6-char flags for type of ocean mascon (PDeep/PShelf etc)
! ternary
! PT170530: added geoid as an extra column (increased to 9)
  deallocate(mcon_ocean_tern)             ! ocean ternary info (lat/lon/rad/area/hgt/density/tern2prim pointer/the ternary number/geoid
  deallocate(mcon_ocean_tern_ptr)         ! pointer to ocean primary mascons in which the ocean ternary resides (needs to be dimensioned for
  deallocate(tern_flag_ocean)             ! 6-char flags for type of ternary mascon (PDeep/PShelf etc)
                                                    !    the maximum number of ALL ternarys, not just ocean ternarys
!!!!!!!!!!!!!!!!!!!!
! tidal amplitudes !
!!!!!!!!!!!!!!!!!!!!
! tidal amplitude partials
  deallocate(ocean_eftid_part)    ! partials for tidal amplitude parameters
  deallocate(mcon_tide_name)                         ! 2-char descriptor of tidal constituents (M2, O1, S2, K1, K2). These are defined in graceorb.f90
  deallocate(mcon_tide_period)                       ! period of the tidal constituents
  deallocate(msc_tide)              ! a priori values for the amplitudes of each tidal constituent (sine and cosine components)

  return
  end subroutine deallocate_mascon_arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_msc_hdr(lumsc,msc_file,hdr_lines,max_hdr_lines,n_hdr_lines,msc_file_code,nprim,nsec,ntern &
                                ,tern_per_prim,tern_per_sec,sec_per_prim)

! subroutine to open the mascon file and read in the header information.
! It needs to be able to work for either temporal gravity mascons or
! ocean mascons
!
! P. Tregoning
! 4 November 2016

  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! passed variables
  integer,       intent(in)    :: max_hdr_lines              ! dimensioning of hdr_lines array     
  character*150, intent(inout) :: hdr_lines(max_hdr_lines)   ! vector of header lines in the mascon file
  integer,       intent(in)    :: lumsc                      ! unit number of mascon file to be opened and read
  character*150, intent(in)    :: msc_file                   ! name of mascon file
  integer*4,     intent(inout)   :: n_hdr_lines                ! actual number of header lines found in the mascon file   
  character*8,   intent(inout)   :: msc_file_code              ! 8-character file code for mascon file
  integer*4,     intent(inout)   :: nprim,nsec,ntern           ! number of ternary, secondary, primary mascons in the file
  integer*4,     intent(inout)   :: tern_per_prim,tern_per_sec,sec_per_prim ! max number of tern/sec per other type of mascon

! local variables
  integer*4      :: ioerr,i
  character*150  :: message
  logical        :: hdr_line
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! open the mascon file
  open (unit=lumsc, file=msc_file, status='old',iostat=ioerr)
  if(ioerr /= 0)then
    call status_update('FATAL','LIB','read_msc_hdr',msc_file,'Error opening mascon file. Does it exist?',0)
  else
    call status_update('STATUS','LIB','read_msc_hdr',msc_file,'Reading mascon information',0)
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the first header line. It contains the mascon file code and the 
! information on the number of primary/secondary/ternary mascons etc
  n_hdr_lines = 1
  read(lumsc,'(a)')hdr_lines(n_hdr_lines)
  ! extract out the mascon file code and the values on the line
  msc_file_code = hdr_lines(1)(1:8)
  read(hdr_lines(1)(9:100),*)nprim ,nsec,ntern,sec_per_prim,tern_per_prim,tern_per_sec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read, store and count the reamining header lines
  hdr_line = .true.
  do while (hdr_line)
    n_hdr_lines = n_hdr_lines + 1
    read(lumsc,'(a)')hdr_lines(n_hdr_lines)
    if(hdr_lines(n_hdr_lines)(1:1) /= "#")then
      hdr_line = .false.
      hdr_lines(n_hdr_lines) = " "
      n_hdr_lines = n_hdr_lines - 1
      backspace(lumsc)
    else
      if(n_hdr_lines < 5)call status_update('STATUS','LIB','read_msc_hdr',' ',hdr_lines(n_hdr_lines),0)
    endif
  enddo

! PT190906: also write out the last four header lines
  if(n_hdr_lines > 7)then
    call status_update('STATUS','LIB','read_msc_hdr',' '," .... The last four header lines ....",0)
    do i=n_hdr_lines-4,n_hdr_lines
      call status_update('STATUS','LIB','read_msc_hdr',' ',hdr_lines(i),0)
    enddo
  endif

  write(message,'(a,i5,a)')"Read ",n_hdr_lines," header lines from input mascon file"
  call status_update('STATUS','LIB','read_msc_hdr',msc_file,message,0)

  return
  end subroutine read_msc_hdr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  subroutine tern_lat_bands(lat_spacing)
!
!! subroutine to establish how many ternary latituide bands and how many ternary mascons per band
!!
!! P. Tregoning
!! 3 September 2016
!
!  use mascon_mod
!
!  implicit none
!
!  real(kind=8) :: lat_spacing     ! separation of ternary mascons (in decimal degrees)
!  integer*4    :: points_per_degree
!  real(kind=8) :: colat
!  integer*4    :: i,num_long,total_nterns
!  real(kind=8) :: del_long
!  real(kind=8) :: pi,rad_fact
!  character*100 :: message
!
!  pi = 4.d0*atan(1.d0)
!  rad_fact = pi/180.d0
!
!
!! work out the number of latitudinal bands
!  n_tern_lat_bands = int(180.e0/lat_spacing) + 1
!  points_per_degree = nint(1.e0/lat_spacing)
!
!  write(message,'(a,f20.17,a,i8)')"  Ternary latitude spacing (degrees):",lat_spacing,' Number of bands: ',n_tern_lat_bands
!  call status_update ('STATUS','LIB','tern_lat_bands',' ',message,0)
!  
!! allocate the latitude bands array
!  allocate(nterns_per_band(n_tern_lat_bands))
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! loop through the latitiude bands and work out how manny ternary  !!
!! mascons in each
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  total_nterns = 0
!  do i = 1, n_tern_lat_bands
!    colat = dble(i - 1)/points_per_degree
!    if ( ( i .eq. 1 ) .or. ( i .eq. n_tern_lat_bands ) ) then
!      nterns_per_band(i) = 1
!    else
!      del_long = 180./(float(n_tern_lat_bands - 1) * sin(colat * rad_fact))
!      nterns_per_band(i) = idint(360./del_long)
!    endif
!    total_nterns = total_nterns + nterns_per_band(i)
!  enddo
!
!  write(message,'(a,i10,a)')'  There are a total of ',total_nterns,' ternary mascons'
!  call status_update ('STATUS','LIB','tern_lat_bands',' ',message,0)
!
!  return
!
!  end subroutine tern_lat_bands
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tern_lat_bands_ell(lat_spacing)

! subroutine to establish how many ternary latituide bands and how many ternary mascons per band
! mascons equally spaced in bands on WGS84 ellipsoid
!
! P. Tregoning
! 3 September 2016
! Mod by H.McQueen  1 Dec 2016


  use mascon_mod

  implicit none

  real(kind=8), parameter         :: Ra = 6378.137e3, inflat = 298.257223563
  real(kind=8) :: lat_spacing     ! separation of ternary mascons (in decimal degrees)
  integer*4    :: points_per_degree
  real(kind=8) :: colat
  integer*4    :: i,num_long,total_nterns
  real(kind=8) :: pi,rad_fact
  character*100 :: message
  real(kind=8) :: bwid,cos_lat,sin_lat,scirc,Rb,Re,Rl

  pi = 4.d0*atan(1.d0)
  rad_fact = pi/180.d0
  Rb=Ra*(1-1/inflat)
  Re=(2*Ra+Rb)/3
  bwid=Re*lat_spacing*rad_fact

! work out the number of latitudinal bands
  n_tern_lat_bands = int(180.d0/lat_spacing) + 1
  points_per_degree = nint(1.d0/lat_spacing)

  write(message,'(a,f20.17,a,i8)')"  Ternary latitude spacing (degrees):",lat_spacing,' Number of bands: ',n_tern_lat_bands
  call status_update ('STATUS','LIB','tern_lat_bands',' ',message,0)

! allocate the latitude bands array
  allocate(nterns_per_band(n_tern_lat_bands))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! loop through the latitiude bands and work out how manny ternary  !!
! mascons in each
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  total_nterns = 0
  do i = 1, n_tern_lat_bands
    colat=(i-1)*lat_spacing
    cos_lat=dcos(colat*rad_fact)
    sin_lat=dsin(colat*rad_fact)
    Rl=(Ra**2*sin_lat)**2+(Rb**2*cos_lat)**2
    Rl=sqrt(Rl/((Ra*sin_lat)**2+(Rb*cos_lat)**2))
    scirc=2*pi*Rl*dsin(colat*rad_fact)
    if ( ( i .eq. 1 ) .or. ( i .eq. n_tern_lat_bands ) ) then
        nterns_per_band(i) = 1
    else
        nterns_per_band(i) = nint(scirc/bwid)
    endif
    total_nterns = total_nterns + nterns_per_band(i)
  enddo
  write(message,'(a,i10,a)')'  There are a total of ',total_nterns,' ternary mascons'
  call status_update ('STATUS','LIB','tern_lat_bands',' ',message,0)

!! verify correct lat band structure
!  open(97,file="junk",status="unknown")
!  do i=1,n_tern_lat_bands
!     write(97,'(2i5)') i,nterns_per_band(i)
!  enddo
!  close(97)
     
  return

  end subroutine tern_lat_bands_ell
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calc_which_ternary(debug,lat,lon,tern_spacing,tern_number)

! subroutine to work out the ternary number, given a lat/lon
!
! P. Tregoning: 2 September 2016
!
! MODS
! HM191111: mascons mislocating by 1 ternary in lat & lon - fixed
    
  use mascon_mod

  implicit none

! passed variables
  logical     , intent(in)  :: debug            ! flag to output debug 
  real(kind=8), intent(in)  :: lat, lon         ! coords of the requested point
  real(kind=8), intent(in)  :: tern_spacing     ! latitudinal width (in dec. degrees) of the ternary mascons
  integer*4,    intent(out) :: tern_number

! other variables
  real(kind=8) :: long_width
  integer*4    :: lat_band,num_long
  character    :: message*250
  real(kind=8) :: pi

  pi = 4.d0*datan(1.d0)

  if(debug)print*,'DEBUG: calc_which_ternary lat, lon, tern_spacing',lat,lon,tern_spacing

  ! Trap to ensure that latitude (not co-latitude) is passed in for satellite location
  if(lat > 90.d0)then
    write(message,'(a,2f8.4,a)')"Methinks that co-latitude has been provided. It must be latitude (lat,lon: " &
                                ,lat,lon,")"
    call status_update('FATAL','LIB','calc_which_ternary',' ',message,0)
  endif

! from the latitude work out which ternary latitude band the point is in
! rounding correction - ensure that N-pole point is in band 1
   
  lat_band = nint((90.d0-lat)/tern_spacing) +1
  if(lat_band == 0)lat_band = 1
  
  if(debug)print*,"DEBUG: colat,spacing,real(colat/spacing)" &
           ,(90.d0-lat),tern_spacing,(90.d0-lat)/tern_spacing,nint((90.d0-lat)/tern_spacing)+1

! from the longitude work out which cell the point is in along the meridian
! wrapping correction - first half cell belongs to last ternary in the band
  
  long_width = 360.d0/nterns_per_band(lat_band)
  num_long = nint(lon/long_width)
  if(num_long.eq.0) num_long=nterns_per_band(lat_band)

! work out the ternary mascon number
  
  tern_number = sum(nterns_per_band(1:lat_band-1)) + num_long 

  if(debug)then
    print*,'DEBUG: latitude band:',lat_band,' of',n_tern_lat_bands,' bands.'
    print*,'DEBUG: ternary is number ',num_long,' mascon on the latitude band'
    print*,'DEBUG: ternary number is',tern_number
  endif


  return
  end subroutine calc_which_ternary
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_mascon_file(luin,mascon_file)

! subroutine to read the primary, secondary and ternary mascon information from 
! the new, single file that contains all information. This replaces the code
! that was in graceorb.f90, which read several different files
!
! P. Tregoning
! 5 September 2016

! MODS
! PT161028: added the bit-mapped tidal amplitude information to the mcon_prim array
! PT161104: mascon file now opened in subroutine read_mascon_header
! PT170530: increased dimensions of mascon vectors to include extra "geoid" column
  use  mascon_mod

  implicit none

  integer*4, intent(in)     :: luin              ! unit number of mascon file (=21)
  character*(*),intent(in)  :: mascon_file

! primary mascon variables
  integer*4    :: num_mcon                       ! running count of how many primary mascons in mascon_file
  character*6  :: tmp_code                       ! primary mascon code (PLand, PDeep, PShelf)
  real(kind=8) :: tmp_prim_info(9)               ! array to read lat/lon/rad/area/hgt/geoid/density/%land/bit-mapped-tide from mascon_file
  integer*4    :: tmp_nsec,tmp_ntern             ! temp read of number of secondary and ternary mascons in a primary mascon
  integer*4    :: tmp_prim                       ! temp read of a primary mascon number
  integer*4    :: current_nsec                   ! running count of how many secondary mascons have been encountered in mascon_file
  character*15 :: tmp_region                     ! name of region assigned to the primary

! secondary mascon variables
  integer*4    :: tmp_sec                        ! temp read of a secondary mascon number
  character*6  :: tmp_sec_code                   ! secondary mascon code (SLand, SDeep, SShelf)
  integer*4    :: tmp_sec_ntern                  ! temp read of number of secondary and ternary mascons in a primary mascon
  real(kind=8) :: tmp_sec_info(7)                ! array to read lat/lon/rad/area/hgt/density/geoid from mascon_file
  integer*4    :: sec2prim                       ! temp read of pointer to primary mascon in which the secondary resides
  integer*4    :: tmp_sec_colour                 ! temp read of colour for secondary mascon

! ternary mascon variables
  integer*4    :: num_tern                       ! running count of how many ternary mascons in mascon_file
  integer*4    :: tmp_tern                       ! temp read of a ternary mascon number
  character*6  :: tmp_tern_code                  ! ternary mascon code (SLand, SDeep, SShelf)
  real(kind=8) :: tmp_tern_info(7)               ! array to read lat/lon/rad/area/hgt/density/geoid from mascon_file
  integer*4    :: tern2prim                      ! temp read of pointer to primary mascon in which the ternary resides
  integer*4    :: tern2sec                       ! temp read of pointer to secondary mascon in which the ternary resides
  integer*4    :: tern_colour(2)                 ! prim/sec colours assigned to the ternary mascon (probably not used in graceorb!)
  integer*4    :: tot_tern                       ! total number of ternary mascons

! local variables
  integer*4                 :: ioerr,isec,itern
  character*300             :: message
  real(kind=8)              :: rad_fact

  rad_fact = (4.0*atan(1.d0)/180.d0)       ! conversion factor of pi/180


! the mascon file has lines as follows (with secondary embedded in primary, and ternary embedded in secondary):
!#00000    4581    4581 1485118       1   69390   69390                                                                                                
!# setup by initialize_mascons  20161103 235759.263 +1100                                                                                              
!# modified by classify_mascons  20161103 235900.834 +1100                                                                                             
!# modified by separate_mascons to remove mixed ocean/land primary mascons.  Continental shelf depth: -150.000                                         
!      1  PDeep       1    227   89.0588  185.7386  6356752.3   87393425632.  -4228.0 1029.   0.0     0       region         
!      1  SDeep     227          89.0588  185.7386  6356752.3   87393425632.  -4228.0 1029.     1   317       region         
!      1  TDeep    90.0000    0.0000  6356752.3     269748154.  -4228.0 1029.     1     1   317     0            region
!      2  TDeep    89.8333   60.0000  6356752.5     359663635.  -4219.0 1029.     1     1   317     0            region

! PT170530: new version with "geoid" column (in metres) between height/depth and density
!      1  Pdeep       1    227   90.0000    0.0000  6356752.3   77604230192.  -3521.1     14.8 1029.   0.0    31          Ocean       
!      1  Sdeep     227          90.0000    0.0000  6356752.3   77604230192.  -3521.1     14.8 1029.     1   580          Ocean       
!      1  Tdeep    90.0000    0.0000  6356752.3     268539828.  -4213.4     14.9 1029.     1     1  1021     2            Ocean       
!      2  Tdeep    89.8333   60.0000  6356752.5     358052557.  -4210.5     15.1 1029.     1     1  1021     3            Ocean       


! loop through the file until the end of file is reached. Store information as we go ....
  num_mcon = 0
  tot_tern = 0
  current_nsec = 1

  do while (ioerr == 0)
    tmp_region = " "
    read(luin,*,iostat=ioerr,end=1000)tmp_prim,tmp_code,tmp_nsec,tmp_ntern,tmp_prim_info,tmp_region
    if(ioerr == 0 .and. tmp_code(1:1) == "P") then
!     this should be a valid primary mascon line. Check that it doesn't exceed the max number of primary mascons
      if(tmp_prim > max_prim)then
        write(message,'(a,i7,a,i7)')"Input file contains too many primary mascons (",tmp_prim &
                                    ,'). Max dimension in mascon_mod.f90:',max_prim      
        call status_update('FATAL','LIB','read_mascon_file',mascon_file,message,0)
      else
!print*,tmp_prim,tmp_code,tmp_nsec,tmp_ntern,tmp_prim_info,tmp_region
      endif
      num_mcon = num_mcon + 1
      num_tern = 0                               ! reset to zero the counter of ternary mascons in this primary mascon

      if ( num_mcon .gt. max_prim ) then
        write(message,'(a,i6,a,i6,a)')"Too many mascons (",num_mcon,"). Max number: ",max_prim,"."
        call status_update('FATAL','LIB','read_mascon_file',' ',message,0)
      endif
!     store the primary mascon information
      mcon_prim(num_mcon,1) = tmp_prim_info(1) * rad_fact    ! latitude of primary mascon
      mcon_prim(num_mcon,2) = tmp_prim_info(2) * rad_fact    ! longitude of primary mascon
      mcon_prim(num_mcon,3) = tmp_prim_info(3)     ! mean radius of primary mascon
      mcon_prim(num_mcon,4) = tmp_prim_info(4)     ! area of primary mascon (as summed up from secondary mascons)
      mcon_prim(num_mcon,5) = tmp_prim_info(5)     ! mean height above sea level
      mcon_prim(num_mcon,6) = tmp_prim_info(7)     ! density of mascon (1000 for continental, 1029 for oceanic
!     pointer information for primary-to-secondary
      mcon_prim(num_mcon,7) = tmp_nsec             ! number of secondary mascons contained within this primary mascon
      mcon_prim(num_mcon,8) = tmp_ntern            ! number of ternary mascons contained within this primary mascon
      mcon_prim(num_mcon,9) = current_nsec         ! (sequential) number of the first secondary mascon in this primary mascon
! bit-mapped tidal information
      mcon_prim(num_mcon,10) = tmp_prim_info(9)    ! bit-mapped tidal amplitude code for primary
! percentage land for the primary mascon
      mcon_prim(num_mcon,11) = tmp_prim_info(8)    ! percentage land of the primary
! PT170530: geoid/ellipsoid separation
      mcon_prim(num_mcon,12) = tmp_prim_info(6)    ! geoid/ellipsoid separation
! 6-char flag for the type of mascon
      prim_flags(num_mcon) = tmp_code
! character region description
      mcon_region(num_mcon)  = tmp_region          ! character*15 region description
!      print*,'imsc,tmp_region',num_mcon," ",mcon_region(num_mcon)
!     increment the number of secondary mascons
      current_nsec = current_nsec + tmp_nsec
    else
print*,'tmp_prim,tmp_code,tmp_nsec,tmp_ntern,tmp_prim_info,tmp_code',tmp_prim,tmp_code,tmp_nsec,tmp_ntern,tmp_prim_info,tmp_code
      write(message,*)"Error somewhere in primary mascon line: " &
        ,tmp_prim,tmp_code,tmp_nsec,tmp_ntern,tmp_prim_info,tmp_code
      call status_update('FATAL','LIB','read_mascon_file',' ',message,ioerr)
    endif

! reda the secondary and ternary mascons that comprise this primary mascon. The next line will be a line with info for
! a secondary, followed by one line for each ternary within the secondary. The pattern repeats until tmp_nsec secondary
! mascon lines - and their corresponding ternary mascons - have been read. At that point, the next line should be again
! a primary mascon line.
    do isec = 1,int(mcon_prim(num_mcon,7))
!     read the secondary mascon line information
      read(luin,*,iostat=ioerr,end=1001)tmp_sec,tmp_sec_code,tmp_sec_ntern,tmp_sec_info(1:nvar_sec-1),sec2prim,tmp_sec_colour
!print*,"the sec line:",isec,tmp_sec,tmp_sec_code,tmp_sec_ntern,tmp_sec_info(1:6),sec2prim,tmp_sec_colour
      if(ioerr == 0) then
!       this should be a valid secondary mascon line
        if(tmp_sec > max_sec)then
          write(message,'(a,i7,a,i7)')"Input file contains too many secondary mascons (",tmp_sec &
                                    ,'). Max dimension in mascon_mod.f90:',max_sec      
          call status_update('FATAL','LIB','read_mascon_file',mascon_file,message,0)
        else
!print*,tmp_sec,tmp_sec_code,tmp_sec_ntern,tmp_sec_info(1:7),sec2prim,tmp_sec_colour
        endif
        mcon_sec(num_mcon,isec,1) = tmp_sec_info(1) * rad_fact ! secondary mascon latitude
        mcon_sec(num_mcon,isec,2) = tmp_sec_info(2) * rad_fact ! secondary mascon longitude
        mcon_sec(num_mcon,isec,3) = tmp_sec_info(3)  ! secondary mascon mean radius 
        mcon_sec(num_mcon,isec,4) = tmp_sec_info(4)  ! secondary mascon area (summed from ternary mascons in secondary)
        mcon_sec(num_mcon,isec,5) = tmp_sec_info(5)  ! secondary mascon mean depth/height above sea level
        mcon_sec(num_mcon,isec,9) = tmp_sec_info(6)  ! secondary mascon mean geoid/ellipsoid separation above sea level    
        mcon_sec(num_mcon,isec,6) = tmp_sec_info(7)  ! secondary mascon density
!       pointer information for secondary-to-primary
        mcon_sec(num_mcon,isec,7) = tmp_sec          ! secondary mascon number
        mcon_sec(num_mcon,isec,8) = tmp_sec_ntern    ! number of ternary mascons in this secondary mascon
        mcon_sec_ptr(tmp_sec) = sec2prim             ! pointer to primary mascon in which this secondary resides
        sec_colour(sec2prim) = tmp_sec_colour        ! colour of this secondary mascon
        sec_flags(tmp_sec)   = tmp_sec_code          ! 6-char code for the secondary mascon (SDeep/SLand etc)
      else
        print*,tmp_sec,tmp_sec_code,tmp_sec_ntern,tmp_sec_info(1:7),sec2prim,tmp_sec_colour
        write(message,'(a,i7)')"Error reading  mascon line for secondary mascon",tmp_sec
        call status_update('FATAL','LIB','read_mascon_file',mascon_file,message,0)
      endif

! now read the ternary mascon information for each ternary contained within this secondary
      do itern=1,int(mcon_sec(num_mcon,isec,8))
        read(luin,*,iostat=ioerr,end=1002)tmp_tern,tmp_tern_code,tmp_tern_info(1:nvar_tern-2) &
           ,tern2prim,tern2sec,tern_colour(1),tern_colour(2)
        if(ioerr == 0)then
          if(tmp_tern > max_tern)then
            write(message,'(a,i7,a,i7)')"Input file contains too many ternary mascons (",tmp_tern &
                                      ,'). Max dimension in mascon_mod.f90:',max_tern      
            call status_update('FATAL','LIB','read_mascon_file',mascon_file,message,0)
          endif
          num_tern = num_tern + 1
          mcon_tern(num_mcon,num_tern,1) = tmp_tern_info(1) * rad_fact ! ternary mascon latitude
          mcon_tern(num_mcon,num_tern,2) = tmp_tern_info(2) * rad_fact ! ternary mascon longitude
          mcon_tern(num_mcon,num_tern,3) = tmp_tern_info(3)  ! ternary mascon radius
          mcon_tern(num_mcon,num_tern,4) = tmp_tern_info(4)  ! ternary mascon area
          mcon_tern(num_mcon,num_tern,5) = tmp_tern_info(5)  ! ternary mascon depth/height above sea level
          mcon_tern(num_mcon,num_tern,6) = tmp_tern_info(7)  ! ternary mascon density
!         pointer information for ternary-to-primary and ternary-to-secondary
          mcon_tern(num_mcon,num_tern,7) = tmp_tern          ! ternary mascon number
! PT170531: add the geoid/ellipsoid separation
          mcon_tern(num_mcon,num_tern,8) = tmp_tern_info(6)  ! ternary mascon geoid/ellipsoid separation
! PT170605: 9th variable is used in voronoi_reshape, so set to zero here
          mcon_tern(num_mcon,num_tern,9) = 0.0               ! used as a pointer in voronoi_reshape

          mcon_tern_ptr(tmp_tern,1) = tern2prim              ! pointer to   primary mascon in which this ternary resides
          mcon_tern_ptr(tmp_tern,2) = tern2sec               ! pointer to secondary mascon in which this ternary resides
          mcon_tern_ptr(tmp_tern,3) = num_tern               ! pointer to the sequential ternary number of this ternary within its primary
          tern_colours(tmp_tern,1)  = tern_colour(1)         ! colour of primary mascon in which the ternary resides
          tern_colours(tmp_tern,2)  = tern_colour(2)         ! colour of secondary mascon in which the ternary resides
          tern_flag(tmp_tern)       = tmp_tern_code          ! 6-char flag for ternary mascon (TDeep/TLand etc)
        else
          print*,tmp_tern,tmp_tern_code,tmp_tern_info(1:nvar_tern-2) &
           ,tern2prim,tern2sec,tern_colour(1),tern_colour(2)
          write(message,'(a,i7)')"Error reading  mascon line for ternary mascon",tmp_tern
          call status_update('FATAL','LIB','read_mascon_file',mascon_file,message,0)
        endif
      enddo  ! end of loop over reading ternary mascons contained within this particular secondary

    enddo  ! end of loop over reading secondary mascon containd within this particular primary

!   increment the total number of ternary mascons
    tot_tern = tot_tern + mcon_prim(num_mcon,8)

  enddo  ! end of generic loop reading through the file (i.e. the loop over "primary mascons" without knowing how many there are in the file)

1000 continue
! store the total number of secondary and primary mascons in the file
  total_sec  = tmp_sec
  total_prim = num_mcon

  write(message,'(a,i7,a)')"Found ",total_prim,"   primary mascons in file "
  call status_update('STATUS','LIB','read_mascon_file',mascon_file(50:150),message,0)
  write(message,'(a,i7,a)')"Found ",total_sec ," secondary mascons in file "
  call status_update('STATUS','LIB','read_mascon_file',mascon_file(50:150),message,0)
  write(message,'(a,i7,a)')"Found ",tot_tern,"   ternary mascons in file "
  call status_update('STATUS','LIB','read_mascon_file',mascon_file(50:150),message,0)

!print*,'subroutine:  mcon_prim(223,:)',mcon_prim(223,:)
!print*,'subroutine: tern_number, mcon_tern_ptr(tern_number,:)',223, mcon_tern_ptr(200:223,:)

! close the mascon_file
  close(luin)
  return

1001 continue
  print*,'must have had a problem reading a secondary mascon line'
  stop

1002 continue
  print*,'must have had a problem reading a ternary mascon line'
  stop


  end subroutine read_mascon_file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_ocean_mascons(luin,ocean_mascon_file)

! subroutine to read in the ocean mascon file. This file contains only primary and
! ternary mascons - we think we don't need the secondary mascons. The plan is to:
!
! 1. compute the tidal accelerations on all ternary mascons at all times. We can
!    put this into an OMP loop to make it efficient. We can eliminate the dist_flag/spotlight
!    computations if we always compute on ternary mascons, so we don't need secondary mascons
! 2. Integrate the partials of the ternary mascons so that tidal amplitudes can be
!    estimated on the primary mascons
!
! P. Tregoning
! 12 September 2016
!
! MODS
! PT170530: added geoid/ellipsoid separation to mascon arrays
! PT170602: removed unnecessary check on ternary mascon number

  use mascon_mod

  implicit none

  integer*4, intent(in)    :: luin
  character*150,intent(in) :: ocean_mascon_file

! primary mascon variables
  integer*4    :: num_ocean_mcon                 ! running count of how many primary mascons in ocean mascon_file
  character*6  :: tmp_code                       ! primary mascon code (PLand, PDeep, PShelf)
  real(kind=8) :: tmp_prim_info(nvar_prim)       ! array to read lat/lon/rad/area/hgt/geoid/density/%land/bit-mapped-tide from mascon_file
  integer*4    :: tmp_nsec,tmp_ntern             ! temp read of number of secondary and ternary mascons in a primary mascon
  integer*4    :: tmp_prim                       ! temp read of a primary mascon number

! ternary mascon variables
  integer*4    :: tmp_tern                       ! temp read of a ternary mascon number
  character*6  :: tmp_tern_code                  ! ternary mascon code (SLand, SDeep, SShelf)
  real(kind=8) :: tmp_tern_info(nvar_tern)       ! array to read lat/lon/rad/area/hgt/geoid/density from mascon_file
  integer*4    :: tern2prim                      ! temp read of pointer to primary mascon in which the ternary resides
  integer*4    :: tern_colour(2)                 ! prim/sec colours assigned to the ternary mascon (probably not used in graceorb!)
  integer*4    :: tot_ocean_tern                 ! total number of ternary mascons
  integer*4    :: itern


! local variables
  real(kind=8)  :: rad_fact
  integer*4     :: ioerr
  character*250 :: message

  rad_fact = (4.0*atan(1.d0)/180.d0)       ! conversion factor of pi/180


!  open (unit=luin, file=ocean_mascon_file, status='old',iostat=ioerr)
  if(ioerr /= 0)then
    call status_update('FATAL','GRACEORB','read_ocean_mascon_file',ocean_mascon_file &
                       ,'Error opening mascon file. Does it exist?',0)
  else
    call status_update('STATUS','GRACEORB','read_ocean_mascons',ocean_mascon_file,'Reading ocean mascon information',0)
  endif

! loop through the file until the end of file is reached. Store information as we go ....
  total_ocean_tern = 0
  total_ocean_prim_ampl = 0

  do while (ioerr == 0)
    read(luin,*,iostat=ioerr,end=1000)tmp_prim,tmp_code,tmp_nsec,tmp_ntern,tmp_prim_info(1:9)
    if(ioerr == 0 .and. tmp_code(1:1) == "P") then
!     this should be a valid primary mascon line. tmp_ntern is how many ternary ocean mascons in this primary ocean mascon.
      num_ocean_mcon = tmp_prim
      mcon_ocean_prim(tmp_prim,1) = tmp_ntern         ! number of ternary mascons in this ocean mascon
! store the bit-mapped tide information
      mcon_ocean_prim(tmp_prim,2) = nint(tmp_prim_info(9))  ! bit-mapped flag for the tidal amplitudes
      prim_flags_ocean(tmp_prim) = tmp_code
      if(mcon_ocean_prim(tmp_prim,2) /= 0)then
        total_ocean_prim_ampl = total_ocean_prim_ampl + 1
        write(message,'(a,i6,a,i6)')" Ocean mascon",tmp_prim," bit-mapped tide flag:",mcon_ocean_prim(tmp_prim,2)
        call status_update('STATUS','LIB','read_ocean_mascons',' ',message,0)
      endif
    else
      write(message,*)"Error somewhere in primary mascon line: " &
        ,tmp_prim,tmp_code,tmp_nsec,tmp_ntern,tmp_prim_info,tmp_code
      call status_update('FATAL','GRACEORB','read_ocean_mascons',' ',message,0)
    endif

! skip the secondary line (NEED TO DO THIS BETTER LATER ON !!!!)
    read(luin,'(a)')message

! read in the ternary ocean mascons
    do itern=1,tmp_ntern
      read(luin,*,iostat=ioerr,end=1002)tmp_tern,tmp_tern_code,tmp_tern_info(1:7),tern2prim
      if(ioerr == 0)then
! PT170602: I don't think this matters - comment it out
!       if(tmp_tern > max_ocean_tern)then
!          write(message,'(a,i7,a,i7)')"Input file contains too many ocean ternary mascons (",tmp_tern &
!                                    ,'). Max dimension in mascon_mod.f90:',max_ocean_tern      
!          call status_update('WARNING','GRACEORB','read_ocean_mascons',ocean_mascon_file,message,0)
!        endif
        total_ocean_tern = total_ocean_tern + 1
        mcon_ocean_tern(total_ocean_tern,1) = tmp_tern_info(1) * rad_fact ! ternary mascon latitude
        mcon_ocean_tern(total_ocean_tern,2) = tmp_tern_info(2) * rad_fact ! ternary mascon longitude
        mcon_ocean_tern(total_ocean_tern,3) = tmp_tern_info(3)            ! ternary mascon radius
        mcon_ocean_tern(total_ocean_tern,4) = tmp_tern_info(4)            ! ternary mascon area
        mcon_ocean_tern(total_ocean_tern,5) = tmp_tern_info(5)            ! ternary mascon depth/height above sea level
        mcon_ocean_tern(total_ocean_tern,9) = tmp_tern_info(6)            ! ternary mascon geoid/ellipsoid separation
        mcon_ocean_tern(total_ocean_tern,6) = tmp_tern_info(7)            ! ternary mascon density
        mcon_ocean_tern(total_ocean_tern,7) = tern2prim                   ! pointer from the sequential ocean ternary number to the ternary number
        mcon_ocean_tern(total_ocean_tern,8) = tmp_tern                    ! the ternary number in the global grid of ternarys
!debug
!print*,'ocean ternary number, ternary number, tern2prim',total_ocean_tern,tmp_tern,tern2prim
!       pointer information for ternary-to-primary and ternary-to-secondary
        mcon_ocean_tern_ptr(tmp_tern,1)   = tern2prim         ! pointer to   primary mascon in which this ternary resides
        mcon_ocean_tern_ptr(tmp_tern,2)   = total_ocean_tern  ! row is ternary mascon number, value is sequential ocean mascon number.
                                                              ! we use this to get tide height at satellite location for debug output
! DEBUG
!print*,'imsc,itern,tmp_tern,mcon_ocean_tern_ptr(tmp_tern,2)',num_ocean_mcon,itern,tmp_tern,mcon_ocean_tern_ptr(tmp_tern,2)

      else
        write(message,'(a,i7)')"Error reading  mascon line for ocean ternary mascon",tmp_tern
        call status_update('FATAL','GRACEORB','read_ocean_mascons',ocean_mascon_file,message,0)
      endif
    enddo  ! end of loop over reading ternary mascons contained within this particular secondary
  enddo

1000 continue
! store the total number of secondary and primary mascons in the file
  total_ocean_prim = num_ocean_mcon

  write(message,'(a,i7,a)')"Found ",total_ocean_prim,"   primary ocean mascons in file "
  call status_update('STATUS','GRACEORB','read_ocean_mascons',ocean_mascon_file(50:150),message,0)
  write(message,'(a,i7,a)')"Found ",total_ocean_tern,"   ternary ocean mascons in file "
  call status_update('STATUS','GRACEORB','read_ocean_mascons',ocean_mascon_file(50:150),message,0)
  write(message,'(a,i6)')" Number of ocean tide/ampl mascons: ",total_ocean_prim_ampl

  return


1002 continue
  print*,'must have had a problem reading a ternary ocean mascon line'
  stop



  end subroutine read_ocean_mascons
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_mascon_record(luout,msc_code,msc_num,nsec,ntern,msc_crds,msc_area,msc_depth,msc_geoid,msc_density &
                                ,percent_land,tide_flag,region,pointers,colours)

! subroutine to write out a primary, secondary or ternary mascon record (ascii)
! to a mascon file. This will standardise the file format across all the mascon
! programs
!
! P. Tregoning
! 9 November 2016

  implicit none

  integer*4,   intent(in) :: luout                ! unit number of output file to which the mascon line will be written
  character*6, intent(in) :: msc_code             ! "PDeep ", "SLand ", "TShelf", "TCasp ", "PEyre ", etc
  integer*4,   intent(in) :: msc_num              ! the number of the actual mascon to be written out
  integer*4,   intent(in) :: nsec,ntern           ! the number of secondary and/or ternary mascons in this actual mascon
  real(kind=8),intent(in) :: msc_crds(3)          ! the mascon coords (lat, lon, radius)
  real(kind=8),intent(in) :: msc_area,msc_depth,msc_geoid,msc_density,percent_land  ! mascon attribute values
  integer*4,   intent(in) :: tide_flag            ! bit-mapped flag for tidal amplitude estimation
  character*15,intent(in) :: region               ! name of region to which mascon relates
  integer*4,   intent(in) :: pointers(2)          ! tern2prim, tern2sec, sec2prim pointers
  integer*4,   intent(in) :: colours(2)           ! up to two colours to assign to the mascon

! local variable
  logical  :: debug

debug = .false.
if(debug )then
  print*,'Arguments passed in:'
  print*,'luout',luout
  print*,' msc_code ',msc_code
  print*,' msc_num',msc_num
  print*,' nsec ',nsec
  print*,' ntern',ntern
  print*,' msc_crds',msc_crds
  print*,' msc_area',msc_area
  print*,' msc_depth',msc_depth
  print*,' msc_geoid',msc_geoid
  print*,' msc_density',msc_density
  print*,' percent_land',percent_land
  print*,' tide_flag',tide_flag
  print*,' region ',region
  print*,' pointers',pointers
  print*,' colours ',colours
endif



! primary mascon line
! PT/HMcQ: increased to i8 for number of mascons in primary and secondary lines
  if(msc_code(1:1) == "P")then
    write(luout,'(i7,2x,a6,2i8,2f10.4,f11.1,f18.0,2f9.1,f6.0,f6.1,i6,2x,a15)')  &
           msc_num,msc_code,nsec,ntern,msc_crds,msc_area,msc_depth,msc_geoid,msc_density,percent_land,tide_flag,region
 
! secondary mascon
  else if(msc_code(1:1) == "S")then
    write(luout,'(i7,2x,a6,i8,8x,2f10.4,f11.1,f18.0,2f9.1,f6.0,2i6,2x,a15)')   &
           msc_num,msc_code,ntern,msc_crds,msc_area,msc_depth,msc_geoid,msc_density,pointers(1),colours(1),region

! ternary mascon
  elseif(msc_code(1:1) == "T")then
    write(luout,'(i7,2x,a6,2f10.4,f11.1,f15.0,2f9.1,f6.0,4i6,9x,a15)')  &
           msc_num,msc_code,msc_crds,msc_area,msc_depth,msc_geoid,msc_density,pointers,colours,region

  endif

 
  end subroutine write_mascon_record
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calc_mascon_xyz(use_ocean_mascons)

! subroutine to convert the input mascon lat/lon/rad into cartesian XYZ
!
! P. Tregoning
! 7 September 2016
!
! MODS
! PT201031: adjust the spherical radii of the ocean ternary mascons for the ellipsoidal height and geoid/ellipsoid separation so 
!           that the ocean tide height changes occur on the geoid, not on a sphere.

  use mascon_mod

  implicit none

  character*1, intent(in) :: use_ocean_mascons
  integer*4               :: iprim,isec,itern
  real(kind=8)            :: tmp_xyz(3)

  do iprim=1,total_prim
! primary mascons
    tmp_xyz(1) = mcon_prim(iprim,3)*cos(mcon_prim(iprim,1))*cos(mcon_prim(iprim,2))
    tmp_xyz(2) = mcon_prim(iprim,3)*cos(mcon_prim(iprim,1))*sin(mcon_prim(iprim,2))
    tmp_xyz(3) = mcon_prim(iprim,3)*sin(mcon_prim(iprim,1))
!   and replace the values of mcon_prim with the computed XYZ values
    mcon_prim(iprim,1:3) = tmp_xyz(1:3)
! secondary mascons
    do isec = 1,int(mcon_prim(iprim,7))
      tmp_xyz(1) = mcon_sec(iprim,isec,3)*cos(mcon_sec(iprim,isec,1))*cos(mcon_sec(iprim,isec,2))
      tmp_xyz(2) = mcon_sec(iprim,isec,3)*cos(mcon_sec(iprim,isec,1))*sin(mcon_sec(iprim,isec,2))
      tmp_xyz(3) = mcon_sec(iprim,isec,3)*sin(mcon_sec(iprim,isec,1))
!     and replace the values of mcon_prim with the computed XYZ values
      mcon_sec(iprim,isec,1:3) = tmp_xyz(1:3)
    enddo

! ternary mascons
    do itern=1,int(mcon_prim(iprim,8))
      tmp_xyz(1) = mcon_tern(iprim,itern,3)*cos(mcon_tern(iprim,itern,1))*cos(mcon_tern(iprim,itern,2))
      tmp_xyz(2) = mcon_tern(iprim,itern,3)*cos(mcon_tern(iprim,itern,1))*sin(mcon_tern(iprim,itern,2))
      tmp_xyz(3) = mcon_tern(iprim,itern,3)*sin(mcon_tern(iprim,itern,1))
!     and replace the values of mcon_prim with the computed XYZ values
      mcon_tern(iprim,itern,1:3) = tmp_xyz(1:3)
    enddo

  enddo
  call status_update('STATUS','GRACEORB','calc_mascon_xyz',' '," Converted all mascons lat/lon/rad into e-fixed XYZ",0)

! check whether there are ternary ocean mascons to convert as well
  if(use_ocean_mascons == "G")then    ! "G" means we are representing the tides on a grid (not as spherical harmonics)
! ternary mascons
    do itern=1,total_ocean_tern
        !PT201031: correct the geocentric radii to the ellopsoid for geoid/ellipsoid separation. This will put the ocean ternarys on the geoid.
        mcon_ocean_tern(itern,3) = mcon_ocean_tern(itern,3) + mcon_ocean_tern(itern,9)

        tmp_xyz(1) = mcon_ocean_tern(itern,3)*cos(mcon_ocean_tern(itern,1))*cos(mcon_ocean_tern(itern,2))
        tmp_xyz(2) = mcon_ocean_tern(itern,3)*cos(mcon_ocean_tern(itern,1))*sin(mcon_ocean_tern(itern,2))
        tmp_xyz(3) = mcon_ocean_tern(itern,3)*sin(mcon_ocean_tern(itern,1))
!       and replace the values of mcon_prim with the computed XYZ values
        mcon_ocean_tern(itern,1:3) = tmp_xyz(1:3)
    enddo
    call status_update('STATUS','GRACEORB','calc_mascon_xyz',''&
                      ," Converted all ternary ocean mascons lat/lon/rad into e-fixed XYZ",0)
  endif



  return
  end subroutine calc_mascon_xyz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine msc_in_region(shift_polygon,lat,lon,poly_lat,poly_lon,npoints,in_or_out)     

! subroutine to determine whether a primary mascon resides within a particular region
! as defined by a polygon of lat/lon coordinates
!
! P. Tregoning
! 11 November 2016

  implicit none

  logical,       intent(in)    :: shift_polygon                       ! indicates whether the polygon was shifted east by 90 deg so that it doesn't cross 0E
  real(kind=8),  intent(in)    :: lat,lon                             ! coordinates of point
  integer*4,     intent(in)    :: npoints                             ! number of vertices of the polygon (first point is not repeated)
  real(kind=8),  intent(inout) :: poly_lat(npoints),poly_lon(npoints) ! coords of vertices of polygon describing region   
  integer*4,     intent(out)   :: in_or_out                           ! -1 if point inside the polygon, 0 if on an edge, 1 otherwise 

! local variables
  integer*4  :: i,j
  real(kind=8),allocatable   :: tmplat(:),tmplon(:)
  real(kind=8) :: temp_lon
  LOGICAL MX,MY,NX,NY                                               

  allocate(tmplat(npoints))
  allocate(tmplon(npoints))

! PT161118: do we need to shift the point's longitude by 90 deg?
  temp_lon = lon
  if(shift_polygon)temp_lon = temp_lon + 90.d0
  if(temp_lon > 360.d0)temp_lon = temp_lon - 360.d0

  do i=1,npoints
    tmplat(i) = poly_lat(i)-lat
    tmplon(i) = poly_lon(i)-temp_lon
  enddo

  in_or_out = -1
  do i=1,npoints
    j=1+mod(i,npoints)

      MX=tmplat(I).GE.0.0                                                    
      NX=tmplat(J).GE.0.0                                                    
      MY=tmplon(I).GE.0.0                                                    
      NY=tmplon(J).GE.0.0                                                    
      IF(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX)) GO TO 2       
      IF(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))) GO TO 3  
      in_or_out=-in_or_out                                                      
      GO TO 2                                                           
3     IF((tmplon(I)*tmplat(J)-tmplat(I)*tmplon(J))/(tmplat(J)-tmplat(I))) 2,4,5                       
4     in_or_out=0                                                           
      RETURN                                                            
5     in_or_out=-in_or_out                                                      


2  enddo

  deallocate(tmplat)
  deallocate(tmplon)
  return
  end subroutine msc_in_region
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine hash_header(headrec,nhr,hashcode)

! generates a hashcode from the supplied header lines
!
!    headrec: array of header records to be hashed
!    nhr: number of header records
!    hashcode: 7 digit CRC checksum prefixed with "#"
!     - trailing blanks in header lines are ignored
!     - uses linux cksum function implementing ISO/IEC 8802-3:1989
!     - ten digit decimal hashcode compressed to base36
!
! H. McQueen
! 15 November 2016

   implicit none

    integer i,nhr,lun
    integer (kind=8) n
    integer m,rem,base
    character*12 hashcode
    character(150), dimension(1:nhr)    :: headrec
    character*1, dimension(1:12)   ::chr

    lun=99
    open(lun,file="headcheck",status="unknown")
    write(lun,'(a120)') (trim(headrec(i)),i=1,nhr)
    close(lun)

    call system("cksum headcheck|awk '{print $1}' >checksum")
    open(lun,file="checksum",status="old")
    read(lun,'(a12)') hashcode
    close(lun)

! convert hashcode to base36

    read(hashcode,*) n
    base=36
    i=0
    do while(n.gt.0)
        m=n/base
        rem=n-m*base
        n=m
        i=i+1
        if(rem.le.9) then
            chr(i)=char(rem+48)
        else
            chr(i)=char(rem+55)
            endif
        enddo
    n=i+1
    do i=1,n-1
        hashcode(i:i)=chr(n-i)
        enddo
    do i=n,12
        hashcode(i:i)=" "
        enddo


    hashcode="#"//hashcode

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_random_value(rand_numbers,nvals,rand_number)

! subroutine to get a random value and check that it hasn't been selected previously

  implicit none

! passed variables
  integer*4, intent(in)     :: nvals
  logical  , intent(inout)  :: rand_numbers(nvals)
  integer*4, intent(out)    :: rand_number

! local variables
  real(kind=8) :: tmp_rand
  logical   :: found 
  integer*4 :: i
  found = .false.
  call srand(nvals)

  do while (.not. found)
    tmp_rand = rand(0)
    rand_number = int(tmp_rand*nvals)+1
!print*,'rand(0),nvals,rand(0)*nvals',tmp_rand,nvals,tmp_rand*nvals,rand_number
    if(.not. rand_numbers(rand_number))then
      rand_numbers(rand_number) = .true.
      found = .true.
    endif

! debug
!    do i=1,nvals
!      print*,'i, T/F: ',i,rand_numbers(i), rand_number,rand_number*nvals
!    enddo

  enddo

  end subroutine get_random_value
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine adjust_mascon_surface(calling_prog,mascon_surface)

! subroutine to adjust the radii of all mascons (primary, secondary, ternary)
! to put them on either the ellipsoid (default), geoid or topographic surface
! The information on geoid/ellipsoid separation is contained in the mascon file,
! along with the topography/bathymetry information. 
!
! For "geoid", we simply add the geoid/ellipsoid separation to the mascon 
! radii (which is what is in the file, placing the mascons on the ellipsoid).
!
! For "topography" we put ocean (and floating ice shelf) mascons on the geoid
! but also add the topographic height above msl to the continental mascons.
!
! P. Tregoning
! 23 October 2018

  use mascon_mod

  implicit none

  character(*) ,intent(in) :: calling_prog
  character*10 ,intent(in) :: mascon_surface           ! input code for which surface the mascons should sit on
  integer*4     :: iprim,isec,itern                    ! counters
  character*250 :: message

! variables to calculate topographic corrections
  integer*4     :: nsec_topo,nprim_topo
  real(kind=8)  :: topo_sec,topo_prim
  real(kind=8),allocatable :: total_prim_adj(:)

! do nothing if the ellipsoidal coordinates are to be used
  if (mascon_surface(1:9) == "ellipsoid")then
    call status_update("STATUS",calling_prog,"adjust_mascon_surface",' ',"Mascons located on the ellipsoid",0)
    return
  else
    allocate(total_prim_adj(max_prim))
    total_prim_adj = 0.d0
  endif


! mascons on the geoid
  if(mascon_surface(1:5) == "geoid" .or. mascon_surface(1:10) == "topography")then
    if(mascon_surface(1:5) == "geoid") &
       call status_update("STATUS",calling_prog,"adjust_mascon_surface",' ',"Mascons located on the geoid",0)
    ! loop through all primary/secondary/ternary mascons and add the geoid/ellipsoid separation to the mascon radius
    do iprim = 1,max_prim
      ! RM190521: only adjust mascons if they aren't GIA mascons
      if (mcon_prim(iprim,6) < 3000.d0) then
        mcon_prim(iprim,3) = mcon_prim(iprim,3) + mcon_prim(iprim,12)
        total_prim_adj(iprim) = total_prim_adj(iprim) + mcon_prim(iprim,12)              ! store the geoid primary mascon adjustment
        do isec = 1,nint(mcon_prim(iprim,7))
          mcon_sec(iprim,isec,3) = mcon_sec(iprim,isec,3) + mcon_sec(iprim,isec,9)
          do itern = 1,nint(mcon_sec(iprim,isec,8))
            mcon_tern(iprim,itern,3) = mcon_tern(iprim,itern,3) + mcon_tern(iprim,itern,8)
          enddo
        enddo
      endif
    enddo
  endif
  write(message,'(a,2(f10.3,a))')"Maximum geoid radial changes to a primary mascon: " &
                            ,minval(total_prim_adj),",",maxval(total_prim_adj)," m."
  if(mascon_surface(1:5) == "geoid")call status_update('STATUS',calling_prog,'adjust_mascon_surface',' ',message,0)

! mascons on the surface topography
  if(mascon_surface(1:10) == "topography")then
    call status_update("STATUS",calling_prog,"adjust_mascon_surface",' ',"Mascons located on the surface topography",0)

! here we need to do nothing for ocean mascons, since the "geoid" approximates mean sea level (close enough anyway).
! For land mascons (including ice-covered land and marine-grounded ice) we need to average the ternary topo values
! (but only those that are "land") and zero for the ocean ones (in the case of a primary with a mix of land and 
! ocean ternarys) to derive the topo correction for the primary mascons. We cannot use the topo values in the 
! mascon file because it is an average of bathymetric depths and topographic heights.

    do iprim = 1,max_prim
      topo_prim  = 0.d0            ! initialise sum of primary topography to zero
      nprim_topo = 0               ! initialise number of land ternarys in this primary to be zero

      do isec =1,nint(mcon_prim(iprim,7))
        topo_sec  = 0.d0            ! initialise sum of secondary topography to zero
        nsec_topo  = 0              ! initialise number of land ternarys in this secondary to be zero

        do itern = 1,nint(mcon_sec(iprim,isec,8))
          if(mcon_tern(iprim,itern,6) < 1010.d0 ) then   ! it is not an ocean or GIA ternary. Only add topo correction for land
            nsec_topo = nsec_topo + 1
            nprim_topo = nprim_topo + 1
            topo_sec  = topo_sec  + mcon_tern(iprim,itern,5)
            topo_prim = topo_prim + mcon_tern(iprim,itern,5)
          endif
        enddo

        ! add the mean topographic correction for this secondary
        if(nsec_topo > 0)mcon_sec(iprim,isec,3) = mcon_sec(iprim,isec,3) + topo_sec/dble(nsec_topo)
      enddo

      ! add the mean topographic correction for this primary
      if(nprim_topo > 0)then
        mcon_prim(iprim,3) = mcon_prim(iprim,3) + topo_prim/dble(nprim_topo)
        total_prim_adj(iprim) = total_prim_adj(iprim) + topo_prim/dble(nprim_topo)
!print*,'mascon ',iprim,' nprim_topo, topo_prim',nprim_topo, topo_prim/dble(nprim_topo)
      endif
    enddo

  endif

! output information concerning the maximum radial changes to the primary mascons
  write(message,'(a,2(f10.3,a))')"Maximum radial changes to a primary mascon: " &
                              ,minval(total_prim_adj),",",maxval(total_prim_adj)," m."
  call status_update('STATUS',calling_prog,'adjust_mascon_surface',' ',message,0)

! c'est tout! 
  deallocate(total_prim_adj)
  return

  end subroutine adjust_mascon_surface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calc_tern_distances(ntern,ternary_dists)

! subroutine to calculate the distances between one ternary and all others
! and then to return them in ascending order of distance
!
! P. Tregoning
! 10 February 2020

  use mascon_mod

  implicit none

! passed variables
  integer*4     :: ntern                              !  number of ternary distances to return for each ternary
  real(kind=8)  :: ternary_dists(max_tern,ntern,2)  !  for each ternary, the ternary number/dist of neighbouring ternarys (in ascending distance)

! local variables
  integer*4     :: itern,itern2,itern3,iprim
  real(kind=8)  :: dists(max_tern),sorted_dists(max_tern)
  real(kind=8)  :: ternary_order(max_tern,2)
  real(kind=8),allocatable  :: all_terns(:,:)
  character*100 :: message

! allocate a single array containing information of all ternarys
  allocate(all_terns(max_tern,4))


! open an output file to write out the information as we go along
  open(33,file="ternary_distances.dat",status='unknown')
  write(33,*)" 1000  closest ternarys and the separation distance"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first, transfer all the ternarys into the single array
  do iprim=1,max_prim
    do itern = 1,nint(mcon_prim(iprim,8))
      all_terns(nint(mcon_tern(iprim,itern,7)),1:3) = mcon_tern(iprim,itern,1:3)
      all_terns(nint(mcon_tern(iprim,itern,7)),4) = mcon_tern(iprim,itern,6)
    enddo
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now, loop over each ternary and calculate the distance to each other ternary.
  do itern = 1,max_tern
    if(mod(itern,100)==0)then
      write(message,'(a,i7)')"Calculating distances to ternary: ",itern
      call status_update('STATUS','UTIL','calc_tern_distances',' ',message,0)
    endif
    dists = 0.d0

    do itern2 = 1,max_tern
      ternary_order(itern2,1) = dble(itern2)
      if(all_terns(itern,4) == all_terns(itern2,4))then
        dists(itern2) = dsqrt( (all_terns(itern,1)-all_terns(itern2,1))**2 + (all_terns(itern,2)-all_terns(itern2,2))**2 &
                         + (all_terns(itern,3)-all_terns(itern2,3))**2)
      else
        dists(itern2) = 4.d5 
      endif

      ! set the distance to 400 km if it was > 400 km. The first 1000 closest ternarys extend to only 332 km distance, so
      ! we can speed up the code by making all values > 400 km the same value.
      if(dists(itern2) > 4.d5)dists(itern2) = 4.d5

    enddo
    dists(itern) = 4.d5     ! set the distance to very large for the ternary with itself, so that it doesn't get used

    ! sort in order of smallest to largest
    !call bubble_sort_R8(max_tern,ternary_order(:,1),dists,ternary_order(:,2),sorted_dists)
    call rapid_sort(.true.,4.d5,max_tern,ternary_order(:,1),dists,ternary_order(:,2),sorted_dists)

    ! transfer the required number of sorted distances to return variable
    do itern3=1,ntern
      ternary_dists(itern,itern3,1) = dble(ternary_order(itern3,2))
      ternary_dists(itern,itern3,2) = sorted_dists(itern3)
!print*,itern,itern3,ternary_dists(itern,itern3,:)
    enddo

! for now, write it out from here
!    write(33,'(i7,1000i7,1000f15.1)')itern,nint(ternary_dists(itern,1:1000,1)),ternary_dists(itern,1:1000,2)
    write(33,*)itern,nint(ternary_dists(itern,1:1000,1)),ternary_dists(itern,1:1000,2)

  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  return
  end subroutine calc_tern_distances
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_msc_model(calling_prog,lumsc_model,msc_model_file,apr_or_RMS,model,soln_epoch,nmascon_t,msc_values)

! read in a file of mascon time series metrics, then apply mascon-specific constraints
! 
! P. Tregoning/J. Pfeffer
! 22 January 2020
!
! MODS:
! PT200225: changed to read cosine/sine amplitudes separately rather than the annual amplitude (allows to reconstruct the signal)
! PT200225: modified to read different header information

  implicit none

! passed variables
  character*15 ,intent(in) :: calling_prog           ! name of main program
  integer*4    ,intent(in) :: lumsc_model            ! unit number for mascon WRMS file
  character*150,intent(in) :: msc_model_file         ! file containing metrics of mascon time series
  character*3  ,intent(in) :: apr_or_RMS             ! flag whether to return a priori mascon values or regularisation constraints
  integer*4    ,intent(in) :: model                  ! bit-mapped combination of metrics to use to assign the mascon value
  real(kind=8) ,intent(in) :: soln_epoch             ! decimal year of the solution epoch
  integer*4    ,intent(in) :: nmascon_t              ! number of mascons expected in file
  real(kind=8) ,intent(out):: msc_values(nmascon_t)  ! mascon values (either apriori, or regularisation sigma)

! variables to read in mascon metric information
  real(kind=8),allocatable :: msc_crds(:,:)          ! mascon longitude/latitude
  real(kind=8),allocatable :: msc_annual(:,:)        ! cosine/sine amplitudes of annual signal for each mascon
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
  character*200            :: message

! variables to interpolate the low-frequency msc metrics to the solution epoch
  real(kind=8)             :: dt
  integer*4                :: i1,i2
  real(kind=8),allocatable :: lowfrq(:)
  real(kind=8)             :: pi

! external functions
  logical                  :: bitmap

  pi = 4.d0*datan(1.d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the file. First line is a comment, giving column headers. Second line contains number of mascons.
  open(lumsc_model,file=msc_model_file,status='old',iostat=ioerr)
  if(ioerr /= 0)call status_update('FATAL','ADDNORM','addnorm',msc_model_file,"Error opening file",0)

  read(lumsc_model,'(a)')line
  read(lumsc_model,*)nmsc,nmonths,ntot

! the "ntot" value in the file includes lat/lon/density/cosine_ampl/sine_ampl
  nrms = ntot - nmonths - 5
  if(nmsc /= nmascon_t)then
    write(message,'(a,i7,a,i7,a)')"Mismatch in number of mascons expected (",nmascon_t,") and found (",nmsc,")"
    call status_update('FATAL','UTIL','read_msc_ts_WRMS',msc_model_file,message,0)
  else
    allocate(msc_crds(nmsc,3))          ! lat, lon, density
    allocate(msc_annual(nmsc,2))        ! cosine/sine amplitudes
    allocate(msc_WRMS(nmsc,nrms))
    allocate(msc_monthly(nmsc,nmonths))
    allocate(epochs(nmonths))
    allocate(lowfrq(nmsc))
  endif

! read the epochs contained in the file. They are in decimal years
  read(lumsc_model,'(a)')line
  read(lumsc_model,*)(epochs(imonth),imonth=1,nmonths)

! ok, now loop through mascons and read in all the information
  read(lumsc_model,'(a)')line
  do imsc=1,nmsc
    read(lumsc_model,*)msc_crds(imsc,:),msc_annual(imsc,:),(msc_wrms(imsc,iRMS),iRMS=1,nrms) &
                ,(msc_monthly(imsc,imonth),imonth=1,nmonths)

  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now, assign a mascon value to each mascon, depending on what model was requested and whether apriori or regularisation
! constraint is what is required
!
  msc_values = 0.d0
  do imsc=1,nmsc

    ! low-frequency component of model. Can be used for both a priori values and regularisation constraints.
    if(bitmap(model,1) )then  
      if(imsc == 1)then

        ! find the two monthly values that straddle the required epoch
        imonth = 1
        !PT220117: fix bug here when extrapolating beyond nmonths
        !do while (epochs(imonth) < soln_epoch .and. imonth <= nmonths)
        !  if(epochs(imonth) < soln_epoch)imonth = imonth+1
        !enddo
        do while ( imonth <= nmonths)
          if(epochs(imonth) < soln_epoch)then
            imonth = imonth+1
          else
            exit ! this will stop it looping
          endif
        enddo
        ! interpolate between i2 and i1
        if(imonth == 1)then
          write(message,'(a,2(f15.4,a))')"Requested epoch (",soln_epoch,") before first epoch in file (",epochs(1) &
                                        ,"). Using model of the first epoch."
          call status_update('STATUS',calling_prog,'read_msc_model',msc_model_file,message,0)
          lowfrq(:) = msc_monthly(:,1)
        else if (imonth > nmonths)then
          write(message,'(a,2(f15.4,a))')"Requested epoch (",soln_epoch,") after last epoch in file (",epochs(nmonths) &
                                        ,"). Using model of the last epoch."
          call status_update('STATUS',calling_prog,'read_msc_model',msc_model_file,message,0)
          lowfrq(:) = msc_monthly(:,nmonths)
        else
          i2 = imonth
          i1 = imonth - 1
          dt = epochs(i2)-epochs(i1)
          write(message,'(3(a,f15.4))')"Requested epoch (",soln_epoch,") between epochs",epochs(i1)," and",epochs(i2)
          call status_update('STATUS',calling_prog,'read_msc_model',msc_model_file,message,0)
          do imsc2=1,nmsc 
            lowfrq(imsc2) = msc_monthly(imsc2,i1) + (soln_epoch-epochs(i1)) * (msc_monthly(imsc2,i2)-msc_monthly(imsc2,i1))/dt
          enddo
        endif
        if(apr_or_RMS(1:3) == "APR")then
          call status_update('STATUS',calling_prog,'read_msc_model',msc_model_file &
                   ,'Including low-frequency signal of time series in a priori model',0)
        else
          call status_update('STATUS',calling_prog,'read_msc_model',msc_model_file &
                   ,'Including low-frequency signal of time series in msc constraint',0)
        endif
      endif
      msc_values(imsc) = msc_values(imsc) + lowfrq(imsc)

    endif

    ! annual variation. Used differently for a priori value vs regularisation constraint. 
    if(bitmap(model,2) )then  
      if(apr_or_RMS(1:3) == "APR")then
          if(imsc == 1)call status_update('STATUS',calling_prog,'read_msc_model',msc_model_file &
                   ,'Including annual signal of time series in a priori model',0)
        msc_values(imsc) = msc_values(imsc) + msc_annual(imsc,1)*dcos(2.d0*pi*soln_epoch) &
                                            + msc_annual(imsc,2)*dsin(2.d0*pi*soln_epoch)
      else
          if(imsc == 1)call status_update('STATUS',calling_prog,'read_msc_model',msc_model_file &
                   ,'Including annual amplitude of time series in msc constraint',0)
        msc_values(imsc) = msc_values(imsc) + dsqrt(msc_annual(imsc,1)**2+msc_annual(imsc,2)**2)
      endif
    endif

    ! wrms of total time series. Only likely to be used for regularisation constraint.
    if(bitmap(model,3) )then  
      if(imsc == 1)call status_update('STATUS',calling_prog,'read_msc_model',msc_model_file &
                   ,'Including WRMS of time series in msc constraint',0)
        msc_values(imsc) = msc_values(imsc) + msc_wrms(imsc,1)**2
    endif

    ! wrms of total time series for land, 0.1m for ocean. Only likely to be used for regularisation constraint.
    if(bitmap(model,4) )then  
      if(imsc == 1)call status_update('STATUS',calling_prog,'read_msc_model',msc_model_file &
                   ,'Including WRMS of time series in msc constraint to land, 0.1m for ocean',0)
      if(msc_crds(imsc,3) < 1010.d0)then
        msc_values(imsc) = msc_values(imsc) + msc_wrms(imsc,1)**2
      else
        msc_values(imsc) = msc_values(imsc) + msc_wrms(imsc,1) + 0.1d0**2    ! set a 0.1 m constraint on oceans if we use land-only WRMS
      endif      
    endif

    ! wrms of high frequency component. Only likely to be used for regularisation constraint.
    if(bitmap(model,5) )then  
      if(imsc == 1)call status_update('STATUS',calling_prog,'read_msc_model',msc_model_file &
                   ,'Including high frequency WRMS of time series in msc constraint',0)
        msc_values(imsc) = msc_values(imsc) + msc_wrms(imsc,4)
    endif
    
    ! annual variation > 0.2 m amplitude 
    if(bitmap(model,6) )then  
if(imsc == 1)print*,"we are in here for model=",model
      if(apr_or_RMS(1:3) == "APR")then
        if(imsc == 1)call status_update('STATUS',calling_prog,'read_msc_model',msc_model_file &
                   ,'Including annual signal > 0.2m of time series in a priori model',0)
        if(dsqrt(msc_annual(imsc,1)**2+msc_annual(imsc,2)**2) > 0.15d0)then
          msc_values(imsc) = msc_values(imsc) + msc_annual(imsc,1)*dcos(2.d0*pi*soln_epoch) &
                                              + msc_annual(imsc,2)*dsin(2.d0*pi*soln_epoch)
        endif
      else
        if(imsc == 1)call status_update('STATUS',calling_prog,'read_msc_model',msc_model_file &
                   ,'Including annual amplitude > 0.2m of time series in msc constraint',0)
        if(dsqrt(msc_annual(imsc,1)**2+msc_annual(imsc,2)**2) > 0.2d0)then
          msc_values(imsc) = msc_values(imsc) + dsqrt(msc_annual(imsc,1)**2+msc_annual(imsc,2)**2)
        endif
      endif
    endif

  enddo

print*,'model was',model

  return
  end subroutine read_msc_model
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





