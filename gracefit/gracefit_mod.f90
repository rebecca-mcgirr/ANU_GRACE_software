!------------------------------------------------------------------------------
! ANU/GRACE, Softwares for GRACE data
!------------------------------------------------------------------------------
!
! MODULE: gracefit_mod
!
!> @author
!> Other people, ANU
!> Sebastien Allgeyer, ANU
!
! DESCRIPTION: 
!> Module containing variable declaration need across different subroutine of the program.
!
! REVISION HISTORY:
! 1900-01-01 - Initial Version
! 2016-08-31 - Converted to module file for easiness. 
! 2019-07-22 - add kb_cos_fft variables  PT/RmcG
! 2020-10-14 - add use_obs array, logical that states whether to use KBR/LRI range/range rate/range accel observatins at each epoch
!------------------------------------------------------------------------------
module gracefit_mod

    use orbits_rw

! PT180625: declare the "mission" variable here so that it is ubiquitous in the code
    integer*4      :: mission    ! GRACE: 0, GRACE FO: 1, GRACE II: 2

! Variables dealing with satellite numbers and names
    integer*4      :: SAT_1      ! Number attributed to satellite 1  (A/C/E?)
    integer*4      :: SAT_2      ! Number attributed to satellite 2  (B/D/F?)
    character(1)   :: satnam(2)  ! Names of satellites

    type(orbits_data) , dimension(4) :: orb_files

!!!!!!!!!!!!!! PARAMETERS !!!!!!!!!!!!!!!!!!!!

    integer*4, parameter :: iRANGE  = 1
    integer*4, parameter :: iRRATE  = 2
    integer*4, parameter :: iRACC   = 3

    double precision, parameter :: SEC_IN_DAY = 86400.0

    integer*4, parameter :: BUTTERWORTH_OFFSET = 14

    double precision,parameter :: epoch_interval = 5.d0  ! Epoch interval. Hardcoded to match what is found KBR1B files

    double precision,parameter :: tdtoff = 51.184d0  !@# No idea

    double precision, parameter :: PI = 4.d0*datan(1.d0)

!   Unit conversions

    double precision, parameter :: m_km = 1.d-3
    double precision, parameter :: m_mm = 1.0d3
    double precision, parameter :: m_um = 1.0d6
    double precision, parameter :: m_nm = 1.0d9

    double precision, parameter :: rad_deg = 180.d0/PI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Maximum number for certain parameters (only those absolutely necessary)
    integer*4, parameter :: maxgprm = 6       ! Maximum number of gps parameters (constant number based on GNV1B files)
    integer*4, parameter :: max_SB = 6        ! Maximum number scale and bias parameters
    integer*4, parameter :: maxrvecprm = 10   ! Maximum number of IC's per satellite (x y z for pos vel, 1/rev sin/cos, 2/rev sin/cos)
    integer*4, parameter :: maxsat = 2        ! Maximum number of satellites
    integer*4, parameter :: max_mcon_tides = 5! Maximum number of tidal constituent amplitudes that can be estimated per mascon

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Number of satellites being considered
    integer*4 :: nsat_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Empirically decomposed kbrr components to use
     double precision :: use_kbrr_decomp(8)
     double precision :: use_kb_nr(2)
     double precision :: use_kb_cos_fft(3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Apply length-dependent mascon constraints
     double precision :: use_msc_len_const

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Apply mascon conservation of mass constraints
     double precision :: use_msc_conserve_mass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Apply accelerometer scale dependencies
     double precision :: use_accel_scale_const
     double precision :: sigma_acc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Number of observations considered
    integer*4 :: ngobs_t
    integer*4 :: nkbr_t
    integer*4 :: nkbrr_t
    integer*4 :: nkbra_t
    integer*4 :: nCoM_t
    integer*4 :: ncond_t
    integer*4 :: naobs_t
    integer*4 :: nsatobs_t
    integer*4 :: nkobs_t
    integer*4 :: nobs_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Number of parameters used
    integer*4 :: ngprm_t
    integer*4 :: nScale_t
    integer*4 :: nBias_t
    integer*4 :: n_onepr_t
    integer*4 :: ntwopr_t
    integer*4 :: nantoff_t
    integer*4 :: nmascons_t
    integer*4 :: norbprm_t
    integer*4 :: nsatprm_t
    integer*4 :: nrvecprm_t
    integer*4 :: nrvec_t
    integer*4 :: nparam_t
    integer*4 :: est_msc_tides        ! estimate (1) or don't estimate (0) the tidal mascon amplitudes
    integer*4 :: nmsc_t               ! number of mascons in a GTORB file
    integer*4 :: nmsc_tid_constit_t   ! total number of tidal amplitudes to be estimated (sum of all constituents on all "active" mascons)
    integer*4 :: nrpy_t               ! roll/pitch/yaw parameters (3 per satellite)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Indices for observations
    integer*4 :: igobs
    integer*4 :: ikbr
    integer*4 :: ikbrr
    integer*4 :: ikbra
    integer*4 :: iCoM
    integer*4 :: icond
    integer*4 :: iacobs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! array to indicate whether to use KBR/LRI observations or not
! array is dimensioned in gracefit/gracesim as (nepochs_t,6)
! PT220518: increased to 8 obs to include pos and vel
    logical,allocatable :: use_obs(:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Indices for parameters
    integer*4 :: iICparm
    integer*4 :: iScale
    integer*4 :: iBias
    integer*4 :: ionepr
    integer*4 :: itwopr
    integer*4 :: iantoff
    integer*4 :: imascons
    integer*4 :: imsctide
    integer*4 :: irpy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Observation numbers

    integer*4, parameter :: minGPSobs = 0  ! Minimum number of GPS observations that need to be used (same for both satellites)
    integer*4, parameter :: minKBobs = 0   ! Minimum number of KB observations that are needed

    integer*4        :: nGPSobs_t(maxsat) ! Number of GPS observations to be used (set in header_GTORB, updated in input_readGNV1B)
    integer*4        :: nKBobs_t          ! Number of Kband observations to be used  (set in header_GTORB, updated in input_readKBR1B)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Parameter names and units
!   character(30) :: prmnam((6+3+3+3+3)*2 + 4600 + 1000*5*2)  ! ICs + rpy + mascons + mascon tidal amplitudes
! RM190307: increased length of prmnam
! PT190905: now use declaration of prmnam from soln_mod
!   character(30) :: prmnam((6+3+3+3+3)*2 + 46000 + 1000*5*2)  ! ICs + rpy + mascons + mascon tidal amplitudes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Observation names and units
    character(32) :: obsnam(6*2+3+6*2+2*2+3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Epoch related variables

    integer*4 :: gnv_interval             ! Interval (in number of epochs) at which GPS observations are considered
    integer*4 :: nepochs_t                ! Number of epochs looked at
    integer*4 :: start_neq,end_neq        ! start/end epochs for stacking the normal equations
    double precision :: spanhrs           ! Number of hours that will be looked at
    double precision :: starting_epoch    ! First epoch in grace seconds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! an array to indicate what sort of model to apply to accelerometer obs (in
! order to have them run flat across the day)
  real(kind=8) :: acc_obs_model(3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Debug variables

    integer*4 :: debug_iepoch             ! the iepoch variable passed through an include file for debug purposes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! L1B file names
  character(80) :: gn_fnam(2)          ! GNV1B files to be used
  character(80) :: kb_fnam             ! KBR1B file to be used
  character(80) :: lri_fnam            ! LRI1B file to be used
  character(80) :: acc_fnam(2)         ! ACC1B files to be used
  character(80) :: thr_fnam(2)         ! THR1B files to be used

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! an array to contain the GNV1B epochs (in gracesecs)
  integer*4, allocatable :: gvec_epochs(:,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! a generic character variable that has the name of the calling program
  character*10  :: calling_prog
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! a variable to determine how much of the gracefit process should be run
! 0: all of it; 1: just to calculating prefit residuals; 2: just to stacking normal equations
  integer*4  :: gracefit_step

! PT211206: logical to indicate whether to use outlier.mask
  logical    :: use_outlier_mask
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module gracefit_mod
