  module mascon_mod

! variables to define and declare the ternary mascon configuration
    integer*4             :: n_tern_lat_bands         ! number of ternary latitude bands
    integer*4,allocatable :: nterns_per_band(:)       ! number of ternary mascons in each latitude band

! variables for the header lines of the mascon file
    integer*4,parameter   :: max_hdr_lines=200              ! maximum number of header lines in the mascon file
    integer*4             :: n_hdr_lines                    ! actual number of header lines in the mascon file
    character*150         :: msc_hdr_lines(max_hdr_lines)   ! array of header lines from/for the mascon file
    character*8           :: msc_hdr_code                   ! unique code identifier for each individual mascon file

! primary mascon array declarations
! PT161104: replace these parameter statements with values read from the header line of the mascon file
!    integer*4,parameter   :: max_prim=8000                  ! maximum number of primary mascons permitted
!    integer*4,parameter   :: max_sec_per_prim=1000          ! maximum number of secondary mascons per primary
!    integer*4,parameter   :: max_tern_per_prim=40000        ! maximum number of ternary mascons per primary
!    integer*4,parameter   :: max_sec = 200000               ! maximum number of secondary mascons
!    integer*4,parameter      :: max_tern=2332800            ! the exact number of 10' ternary mascons (1484568). For a regular 10' grid there are 2332800
    integer*4   :: max_prim              ! maximum number of primary mascons permitted
    integer*4   :: max_sec_per_prim      ! maximum number of secondary mascons per primary
    integer*4   :: max_tern_per_prim     ! maximum number of ternary mascons per primary
    integer*4   :: max_tern_per_sec      ! maximum number of ternary mascons per secondary  (not sure if we need this ... !)
    integer*4   :: max_sec               ! maximum number of secondary mascons
    integer*4   :: max_tern              ! the exact number of 10' ternary mascons (1484568). For a regular 10' grid there are 2332800

! variables that define the number attributes (i.e. columns) of primary, secondary and ternary arrays (lat.lon/density etc)
    integer*4,parameter   :: nvar_prim= 12
    integer*4,parameter   :: nvar_sec =  8
    integer*4,parameter   :: nvar_tern = 9
 
    real(kind=8),allocatable          :: mcon_prim(:,:)     ! primary mascon info (lat/lon/rad/area/hgt/density/#sec/#tern/#first secondary/tidal flag/%land/geoid)
    real(kind=8),allocatable          :: mcon_prim_EWH(:)   ! EWH a priori value for each mascon
    integer*4                         :: total_prim         ! total number of primary mascons (same as num_mcon)
    character*6,allocatable           :: prim_flags(:)      ! an array of primary mascon flags (PDeep/Pland/PCasp/PEyre etc)
    character*15,allocatable          :: mcon_region(:)     ! region character code for each primary mascon
    integer*4,   allocatable          :: mcon_flag(:)       ! flag to leave as primary or break into secondary or ternary
    logical                           :: use_flag_file      ! logical to use or ignore distances between mascons.
    real(kind=8),allocatable          :: mcon_EWH_vector(:) ! EWH variable used in graceorb for a priori values for each mascon (prim/sec/tern) used in orbit integration

! secondary mascon array declarations
    real(kind=8),allocatable          :: mcon_sec(:,:,:)    ! secondary info ( lat/lon/rad/area/hgt/density/sec_number/#terns in this secondary/geoid )
    real(kind=8),allocatable          :: mcon_sec_EWH(:)    ! EWH a priori value for each secondary mascon
    integer*4  ,allocatable           :: mcon_sec_ptr(:)    ! pointer to primary mascon in which this secondary resides
    integer*4                         :: total_sec          ! total number of secondary mascons
    real(kind=8),allocatable          :: sec_colour(:)      ! colour of secondary mascon
    character*6,allocatable           :: sec_flags(:)       ! an array of secondary mascon flags (SDeep/Sland/SCasp/SEyre etc)

! ternary mascon array declarations
    real(kind=8),allocatable :: mcon_tern(:,:,:)            ! ternary info ( lat/lon/rad/area/hgt/density/tern_number/geoid )
    real(kind=8),allocatable :: mcon_tern_EWH(:)            ! EWH a priori value for each ternary mascon
    integer*4,allocatable    :: mcon_tern_ptr(:,:)          ! pointer to primary and secondary mascons in which the ternary resides/seq. ternary number within primary
    real(kind=8),parameter   :: ternary_lat_spacing=10.d0   ! we currently use a 10' ternary latitude spacing
    character*6,allocatable  :: tern_flag(:)                ! an array of ternary mascon flags (TDeep/TLand/TShelf/TCasp etc)
    integer*4,allocatable    :: tern_colours(:,:)           ! colours of primary and secondary mascons in which the ternary mascon resides

! PT140606: the following are used in a list of mascons (mix of prim/second/tern) to be used at a particular epoch.
    real(kind=8),allocatable :: mcon_xyz(:,:)               ! XYZ Efixed coordinates (m)
    real(kind=8),allocatable :: mcon_area(:)                ! mascon areas (sq. m.)
    integer*4,allocatable    :: mcon_num(:)                 ! mass con counter for gravity & partials
    integer*4,allocatable    :: mcon_rho(:)                 ! density of the mascons in the list 

! PT140613: partials that relate mascons to acceleration, tidal amplitudes at mascons to acc, and mascons to position
    real(kind=8),allocatable :: mcon_efacc_part(:,:)     ! partials relating mascons to accelerations (originally computed in ICpart_calc)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                               !!
!! O C E A N    M A S C O N    V A R I A B L E S !!
!!                                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! we will set up the ocean mascons as primary, secondary and ternary but, at this stage, I think we will only be using
! two levels (primary and ternary). The secondary ocean mascons might simply disappear from this file if it turns out 
! that we don't need them.

! variables for the header lines of the ocean mascon file
    integer*4,parameter   :: max_ocn_hdr_lines=30                 ! maximum number of header lines in the mascon file
    integer*4             :: n_ocn_hdr_lines                      ! actual number of header lines in the mascon file
    character*150         :: msc_ocn_hdr_lines(max_ocn_hdr_lines) ! array of header lines from/for the mascon file
    character*8           :: msc_ocn_hdr_code                     ! unique code identifier for each individual mascon file

    integer*4   :: max_ocean_sec_per_prim      ! maximum number of secondary mascons per primary
    integer*4   :: max_ocean_tern_per_prim     ! maximum number of ternary mascons per primary
    integer*4   :: max_ocean_tern_per_sec      ! maximum number of ternary mascons per secondary  (not sure if we need this ... !)
    integer*4   :: max_ocean_sec               ! maximum number of secondary mascons
    integer*4   :: max_ocean_tern              ! the exact number of 10' ternary mascons (1484568). For a regular 10' grid there are 2332800

! primary ocean mascons
! PT161104: replce this parameter statement with the value read from the header line of the mascon file
!  integer*4, parameter     :: max_ocean_prim=8000
  integer*4                :: max_ocean_prim               ! maximum number of ocean primary mascons
  integer*4, allocatable   :: mcon_ocean_prim(:,:)         ! primary ocean info (#tern/bit-mapped tide)
  integer*4                :: total_ocean_prim             ! total actual number of primary ocean mascons
  character*6,allocatable  :: prim_flags_ocean (:)         ! an array of primary ocean mascon flags (PDeep/Pland/PCasp/PEyre etc)

! ternary ocean mascons
  real(kind=8),allocatable :: mcon_ocean_tern(:,:)       ! ocean ternary info (lat/lon/rad/area/hgt/density/geoid/s???
  integer*4,allocatable    :: mcon_ocean_tern_ptr(:,:)   ! pointer to ocean primary  mascons in which the ocean ternary resides
  integer*4                :: total_ocean_tern           ! total number of ocean ternary mascons
  character*6,allocatable  :: tern_flag_ocean(:)         ! an array of ternary ocean mascon flags (TDeep/TShelf etc)

! ocean mascon tidal amplitude partials: ocean_eftid_part(3,max_msc_tides,2,max_ocean_prim)
  real(kind=8),allocatable :: ocean_eftid_part(:,:,:,:)  ! array to hold partial derivatives of tidal amplitude partials 



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! M A S C O N    T I D A L    A M P L I T U D E S
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT140612: variables to describe the tidal amplitudes that we may choose to estimate at certain ocean mascons
  integer*4, parameter     :: max_msc_tides = 5       ! maximum number of tidal constituents that can be estimated
  character*2,allocatable  :: mcon_tide_name(:)       ! 2-char descriptor of tidal constituents (M2, O1, S2, K1, K2). These are defined in graceorb.f90
  real(kind=8),allocatable :: mcon_tide_period(:)     ! period of the tidal constituents
  real(kind=8),allocatable :: msc_tide(:,:,:)         ! a priori values for the amplitudes of each tidal constituent (sine and cosine components)
  integer*4                :: num_mcon_tide           ! total number of mascon tidal constituents that are to be estimated  
  integer*4                :: total_ocean_prim_ampl   ! number of ocean primarys for which we want to estimate tidal amplitudes
  integer*4,   allocatable :: ocean_ternary_no_ampl(:)! ternary numbers of all ternarys for which we want no modelling or  estimation of tidal amplitudes
  integer*4,   allocatable :: ocean_ternary_ampl(:)   ! ternary numbers of all ternarys for which we DO want modelling and estimation of tidal amplitudes
  integer*4                :: n_tern_ampl,n_tern_no_ampl ! the number of ternary with and without tidal amplitudes involved.


  end module mascon_mod

  
