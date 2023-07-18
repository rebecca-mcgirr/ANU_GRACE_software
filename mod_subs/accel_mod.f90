! Module needed for accelerometer data and rotation matrix
! needed for the partials computation
  module accel_mod

  use mascon_mod

     real(kind=8), dimension(3,3) :: Qrow
     real(kind=8), dimension(0:3) ::  quat
     real(kind=8), dimension(3,1) :: accrom 
     real(kind=8), dimension(3) :: accobs 

! Elements for Mass-con force and partials.  Initially we use just one, but allow for more

! MODIFIED: March 29 2012 - T Purcell to introduce three layers of mascon
! Layer 1 - The primary mascons represent the coarsest resolution for far-field
!           calculations (angular distance > 20 degrees). The mascon cell size
!           resolution is initially set to 330 km squares (3 degrees latitude).
! Layer 2 - The secondary mascons represent an intermediate resolution to be
!           used in the intermediate field (5 deg < angular distance < 20 deg).
!           Secondary mascon cell size is set to 55 km squares (30 minutes lat).
! Layer 3 - The ternary mascons represent the finest resolution to be used in
!           the near field (angular distance < 5 degrees). Ternary mascon cell
!           size is set to 18 km squares (10 arc minutes of latitude). 

    integer*4, parameter:: max_mcon1 = 8000    ! Maximum number of primary mass cons allowed
    integer*4, parameter:: max_mcon2 = 50      ! Maximum number of secondary mass cons per primary cell
    integer*4, parameter:: max_mcon3 = 350     ! Maximum number of ternary mass cons per primary cell
    integer*4, parameter:: ppd_tern = 6        ! Max number of ternary cells per degree
    integer*4, parameter:: max_tern_accelmod = 15000  ! Max number of ternary cells
    integer*4, parameter:: max_tern_lat = 1 + (180 * ppd_tern) ! Max number of ternary latitude steps

    integer :: num_lat_sat, num_lon_sat, num_cell_sat, istep
    real(kind=8) :: r, rad_fact, sat_rad, sat_lat, sat_lon, del_lat, del_lon

    integer*4  num_mcon, tmp_mcon       ! Number of mass cons
! PT140606: stuff
    real(kind=8) bs(3)                        ! XYZ accelerometer biases for current epoch
    real(kind=8) mcon_fract_ocean(max_tern_accelmod)   ! fraction of mascon that is ocean (0 = all land, 1 = all ocean)
    real(kind=8) :: sum_land_area_secondary(max_mcon1*max_mcon2),sum_ocean_area_secondary(max_mcon1*max_mcon2)

! primary mascon variables
    integer*4 mcon_rad1(max_mcon1)     ! distance of primary masons from the centre of the Earth (m)
    real(kind=8) mcon_area1(max_mcon1)       ! Area of primary mass cons (sq. m)
    integer*4 mcon_rho1(max_mcon1)     ! density of primary mass cons (kg/m^3)
    real(kind=8) mcon_dh(max_mcon1)          ! Thickness of mass cons (m)
    real(kind=8) mcon_colat1(max_mcon1)      ! co-Latitudes  of mass cons (rads)
    real(kind=8) mcon_lng1(max_mcon1)        ! Longitudes of mass cons (rads)
    integer*4 mcon_tide1(max_mcon1)    ! bit-mapped flag for whether to estimate tidal corrections (1: M2, 2: O1, 4: S2, 8: K1, 16: K2)
    real(kind=8) mindist_msc2sat(max_mcon1)  ! minimum distance of each mascon to the satellite throughout the orbit integration
    real(kind=8) mcon_fract_ocean1(max_mcon1)! fraction of area of primary mascon that is ocean
! secondary mascon variables
    integer*4 num_tot_sec                       ! total number of secondary mascons
    integer*4 num_mcon2(max_mcon1)              ! Number of secondary mascons in each primary mass con cell
    integer*4 mcon_rad2(max_mcon2, max_mcon1)   ! distance of secondary mascons from the centre of the Earth (m)
    integer*4 mcon_rho2(max_mcon2,max_mcon1)    ! density of secondary mass cons (kg/m^3)
    real(kind=8) mcon_area2(max_mcon2, max_mcon1)     ! Area of mass cons (sq. m)
    real(kind=8) mcon_colat2(max_mcon2, max_mcon1)    ! co-Latitudes of mass cons (rads)
    real(kind=8) mcon_lng2(max_mcon2, max_mcon1)      ! Longitudes of mass cons (rads)
    real(kind=8) mcon_fract_ocean2(max_mcon1*max_mcon2)    ! fraction of area of each secondary mascon that is ocean 

! ternary mascons
    integer*4 num_tot_tern                      ! total number of secondary mascons
    integer*4 num_mcon3(max_mcon1)              ! Number of ternary mascons in each primary mass con cell
    integer*4 mcon_rad3(max_mcon3, max_mcon1)   ! distance of ternary of mascons from the centre of the Earth (m)
    integer*4 mcon_rho3(max_mcon3,max_mcon1)    ! density of ternary mass cons (kg/m^3)
    real(kind=8) mcon_area3(max_mcon3, max_mcon1)     ! Area of mass cons (sq. m)
    real(kind=8) mcon_colat3(max_mcon3, max_mcon1)    ! co-Latitudes  of mass cons (rads)
    real(kind=8) mcon_lng3(max_mcon3, max_mcon1)      ! Longitudes of mass cons (rads)


    real(kind=8) mcon_colat(max_mcon3*max_mcon1)      ! co-Latitudes  of mass cons (rads)
    real(kind=8) mcon_lon(max_mcon3*max_mcon1)        ! Longitudes of mass cons (rads)
    ! an array for the tide heights on coastal ternary mascons ("coastal ternary" is an ocean ternary in a secondary of a mix of ocean/continent)
    real(kind=8) coastal_ternary(30000,2,290)  ! there are around 30,000 of them and 145 tide epochs/day, so this approximate dimension should be big enough)

! This array will allow us to determine which mascon cells should be modeled
! using primary, secondary, or ternary mascon resolution. It will be read
! from a pre-processed file.

    integer*4 pts_per_deg          ! number of ternary latitude steps per degree - gives resolution of ternary cells
    integer*4 num_tern_long(max_tern_lat) ! no. ternary cells for each lat
    integer*4 num_tern_lat         ! no. ternary latitude steps
    integer*4 mcon_count(max_tern_lat)    ! running total count of ternary cells up to current latitude
    integer*4 mcon_index(max_tern_accelmod) ! mass con index giving the primary cell that each ternary cell lies in
    integer*4 dist_flag(800*max_mcon1)   ! mass con counter for gravity & partials

! These arrays will allow us to convert satellite position into a corresponding
! ternary mascon cell, and calculate which primary cell the ternary cell lies
! in. These quantities are read from a pre-processed file. 

    integer*4 :: nd ! number of parameters for partials (plus 1 for pos, vel)


! PT140605: define new arrays, being the pointer of a primary/secondary/ternary mascon to its nearest secondary mascon
!           We need this for applying the ocean tide heights via the mascons
    integer*4 :: mcon_prim2sec(max_mcon1)
    integer*4 :: mcon_sec2sec(max_mcon2,max_mcon1)
    integer*4 :: mcon_tern2sec(max_mcon3,max_mcon1)
    integer*4 :: mcon_tide_ptr(max_mcon3*max_mcon1)

! Now we define the arrays of deformational components of gravitational
! acceleration.
    integer*4, parameter:: max_def_colat = 1200 ! Max. no. colatitudes @ which deformation is tabulated
    integer*4, parameter:: max_def_rad = 50     ! Max. no. altitudes @ which deformation is tabulated

    integer*4 nrad_def, ntheta_def
    real(kind=8) drad_def, dtheta_def, rad0_def, theta0_def

    real(kind=8) def_up(max_def_rad, max_def_colat)   ! deformation terms (radial)
    real(kind=8) def_tang(max_def_rad, max_def_colat) ! deformation terms (tangential)

    real(kind=8) tang_def_vector(3,8 * max_mcon1)     ! XYZ Efixed coordinates (m) for the tangential component of the deformation-induced change in gravity
    real(kind=8) up_def_vector(3,8 * max_mcon1)       ! XYZ Efixed coordinates (m) for the radial component of the deformation-induced change in gravity



  end module accel_mod

