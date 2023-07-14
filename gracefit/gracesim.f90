!********************************************************************************************************************************
!  File: gracesim.f90
!
! -------------------------------------------------------------------------------------------------------
!  ->->->  This is GRACEFIT, but changed to read two sets of GTORB files: 
!          1. a "truth" set, from which we generate GPS and KBRR 'observations'
!          2. a "theoretical" set, from which we try to recover the "truth"
!  P. Tregoning
!  15 October 2013.
!
!  I have made modifications to the GRACEFIT code only insofar as necessary to replace the reading and use of the 
!  GNV1B and KBR1B files with the "truth" set of GTORB files. Essentially, I 
!
!  a) read the "truth" GTORB files before the theoretical GTORB 
!  b) store away the rvec vector (position and velocity)
!  c) read the theoretical GTORB as usual
!  d) turn off the read of the GNV1B and KBR1B
!  e) calculate the range rate from the "truth" positions/velocities
!  f) move the "truth" rvec values into gvec (where the GNV1B values are stored in gracefit)
!  g) store the range rate derived from the "truth" GTORB in the KBrange array 
!
! ---------------------------------------------------------------------------------------------------------
!
!
!  Purpose: To fit (in a least squares sense) tabular points from one or more grace
!           orbit files (GNV1B), using partials and a priori
!           tabular points from an GTORB integrated t-file.
!
!  Author:  Thomas Greenspan
!           (this is the restructured version of the program of the same name by
!           Simon McClusky, worked on and modified by Paul Tregoning)
!
! Contact info:   S. McClusky simon.mcclusky@anu.edu.au
!                 P. Tregoning paul.tregoning@anu.edu.au
!
! Input files:  Command file
!               GTORB file(s) (from graceorb)
!               GNV1B file(s)
!               ACC1B file(s)
!               THR1B file(s)
!               KBR1B file
!               ut1, pole, and nutation files
! PT130620: The sun position is taken from file grace/tables/JPLEPH and read using the JPL
!           subroutines in jplephem.f (which, incidentally, both names and opens the file !!) !@#
!
! Output files:  .rms file
!                .svs file
!                .vcv file
!                .fit file
!                .resid file
!                .corr file
!                plot .kb file
!                plot .A/B file(s)
!
! MODS:
! PT181119: add the reading of an actual KBR1B file, from which to extract the list of epochs when KBR observations are actually 
!           available (rather than just creating them for every epoch)
!
!********************************************************************************************************************************

  Program GRACESIM

! PT170609: include the mod file to define the record length of the GTORB files
    use gtorb_mod
! PT171129: update to remove the use of gracefit.h90
  USE gracefit_mod
  use accred_mod      ! declares the variables for the accelerometer observations
  use soln_mod        ! PT190905: adds the declaration of the arrays for the solution (apriori, adjust, soln)


    implicit none

! PT181120: we need this include so that we know the unit number of the KBR1B file
  include 'input.h90'

!********************  Variable declarations ********************

! Paramters
    integer*4, parameter :: NUM_EXPECTED_ARGS = 14  ! Number of expected arguments.

! Status variables
    character(120) :: version         ! Version of program being run
    character(256) :: message         ! Message to be written out to standard output

! Command line related variables
    character(80) :: cmdnam,rmsnam                         ! Files to be read from the command-line
!    character(80) :: gt_fnam(2)                            ! GTORB file names to be used
    character(80) :: truth_fnam(2)                         ! GTORB file names for the "truth" GTORB files
    character(80) :: arg                                   ! Dummy variable to be read in iop from the command-line
    integer*4     :: iop                                   ! Plot resid indicator
    character     :: date_str*10,c_yr*4,c_month*2,c_day*2  ! Variables for reading the date from the command-line

! Variables storing data from command file and GTORB file(s)
    double precision, allocatable :: apr_prm(:)               ! A priori values for parameters
    double precision, allocatable :: apr_const(:)             ! A priori constant values
    double precision, allocatable :: apr_tide_ampl_const(:,:) ! A priori tidal amplitude constant values
    double precision, allocatable :: apr_wght(:)              ! A priori data weights
    double precision :: kbrr_prefit_tol                       ! Orbit misfit and kbrr prefit tolerances

! Data read in or calculated from GTORB file(s)
  integer*4,parameter             :: maxsats = 2
    integer*4                     :: rec_len(2)                    ! record length of each GTORB file
    double precision, allocatable :: rvec(:,:,:)                   ! Vector of positions, velocities and partials
    double precision, allocatable :: truthvec(:,:,:)               ! Vector of positions, velocities and partials for "truth" GTORB
    double precision, allocatable :: sciframe_quat(:,:,:)          ! Quaternions for each satellite for every epoch
    double precision, allocatable :: srf2trf_rotmat(:,:,:,:)       ! SRF to TRF rotation matrix for GRACE A and B
    double precision, allocatable :: srf2trf_deriv_rotmat(:,:,:,:) ! Differentiated SRF to TRF rotation matrix for GRACE A and B
    double precision :: satics_t(maxrvecprm+1,2)                   ! A priori values for IC's
    double precision :: apr_ScaleBias(max_SB,maxsats)               ! A priori values for scale and bias
    double precision :: GPSantoff(7,maxsats)                        ! Quaternion and offsets for antenna 
! PT170213: add mascon filenames and codes
    character*40  :: combined_mascon_file,ocean_mascon_file        ! name of temporal and ocean mascon files used in GRACEORB
    character*8   :: msc_hdr_code,msc_ocn_hdr_code,msc_reg_code    ! mascon file codes of temporal and ocean mascon files used in GRACEORB
    integer*4     :: total_ocean_prim

! Variables used for rpy calculations and antenna offsets
    double precision, allocatable :: rpy(:,:,:)          ! Roll, Pitch, Yaw of satellites for each epoch
    double precision, allocatable :: AOC(:,:,:)          ! Antenna offset corrections
    double precision, allocatable :: GPS_antoff_trf(:,:) ! Antenna offset in TRF
    double precision, allocatable :: GPSant_adj(:,:)     ! Antenna offsets to be read in from command line

! Data read in from ACC1B file(s)
    double precision, allocatable :: acc(:,:,:)    ! Acceleration data

! Data read in or calculated from THR1B file(s)
    double precision, allocatable :: thr(:,:,:)          ! Thruster data  !@# Change index order (does it even need to be here?)
    logical, allocatable          :: thrusters_off(:,:)  ! Data on whether thrusts are affecting accelerations measurements

! Data read in from GNV1B files
    double precision, allocatable :: gvec(:,:,:)      ! Positions and velocities for each satellite at every epoch
    integer*4       , allocatable :: num_gps_used(:)  ! number of GPS observations that will be used (after decimating the 5sec data)
    logical         , ALLOCATABLE :: use_gvec(:,:)          ! logical to indicate whether there are GNV1B obs to match the GTORB obs
    integer*4                     :: igvec             ! counter through the GNV1B records
! PT190528: pass this variable back from input_readGNV1B
    integer*4                     :: n_gnv             ! number of GNV1B epochs read from the input GNV1B file

! Data read in from KBR1B file
    double precision, allocatable :: kbrange(:,:)    ! Range value (biased), rate, and acceleration
    double precision, allocatable :: kblt(:,:)       ! Light time corrections for range, rate and acceleration
    double precision, allocatable :: kbant(:,:)      ! Antenna phase center corrections for range, rate and acceleration
    double precision, allocatable :: kbrange_L1B(:,:)    ! Range value (biased), rate, and acceleration
    double precision, allocatable :: kblt_L1B(:,:)       ! Light time corrections for range, rate and acceleration
    double precision, allocatable :: kbant_L1B(:,:)      ! Antenna phase center corrections for range, rate and acceleration
!   Variables that might be used in the future but not needed in current version program
!    double precision, allocatable :: kbion_cor(:)  ! Ionospheric range correction for Ka frequencies
!    double precision, allocatable :: kbsnr(:,:,:)  ! SNR K and Ka bands for both satellites
    integer*4                     :: kbrr_misfit_num ! Number of kbrr misfits
  ! PT170703: flag to override the filtering of kbr data
  LOGICAL                       :: override_kbr_filt

! PT181120: variables used when reading the actual KBR1B file
    integer*4 :: n_kbr,num_epochs
    logical   :: use_kbr_epochs        ! logical flag as to whether to use only the epochs for which KBR obs exist

! Counter variables
!    integer*4 :: iepoch    ! Counter that runs through epochs (numbered)
    integer*4 :: isatt      ! Counter that runs through satellites
    integer*4 :: i,j       ! Counter variable

! Variables used for calculations
    double precision, allocatable :: pre_omc(:,:)   ! Prefit Observed Minus Computed values to be used in LS algorithm
    DOUBLE PRECISION, ALLOCATABLE :: pre_omc_kb_unfilt(:,:)   ! KBRR unfiltered prefit  OMC, backup copy
    DOUBLE PRECISION, ALLOCATABLE :: post_omc_kb_unfilt(:,:)  ! KBRR unfiltered postfit OMC, backup copy
    real(kind=8 ), allocatable :: normeq(:,:)    ! The left side of the LS algorithm. P_t*W*P
! PT221122: AtWb now declared in soln_mod
!    double precision, allocatable :: AtWb(:)        ! The right side of the LS algorithm. P_t*W*OMC
    double precision, allocatable :: adjust(:)      ! Solution to normal equations
    double precision, allocatable :: post_omc(:,:)  ! Postfit Observed Minus Computed values
    double precision, allocatable :: beta_angle(:)  ! beta angle at each GPS epoch
    DOUBLE PRECISION, ALLOCATABLE :: bvec(:)

! Variables used for application of shadow condition
    integer*4, allocatable :: shadow_counter(:)       ! Variable that counts the number of epochs the satellites are in shadow for
    logical, allocatable   :: apply_shadow_cond(:,:)  ! Data on whether the shadow condition can be applied
    double precision, allocatable :: shadow(:,:)      ! Variable used to indicate whether satellite is in shadow

! PT130813: variables used for inversion using just mascons and emf-decomposed kbrr omc
    double precision, allocatable :: filt_kbrr_omc(:)      ! sum of the decomposed elements of the kbrr omc
    double precision, allocatable :: filt_kbra_omc(:)      ! sum of the decomposed elements of the kbra omc
    double precision, allocatable :: kb_1pr(:)

! PT130826: variables to apply length-dependent mascon constraints
    double precision, allocatable :: msc_const(:,:)
    character                     :: msc_const_file*100
    integer*4 :: trimlen        ! Function that returns the length of a string
    integer*4 :: ioerr

! PT140527: mascon ocean/land identifier
    logical,          allocatable :: mcon_ocean(:)
    integer*4                     :: mcon_tides(8000,2)    ! we have to declare the dimensions of this array before calling header_readGTORB ....
    real(kind=8),     allocatable :: msc_tide_amp(:,:,:)

! PT131202: variables for outputting prefit residuals before stacking normal equations
    double precision :: satlon, satlat, sumsq(14)

! PT130902: variables for openMP stacking of normal equations
! PT201006: Seb apparently increased normeq_tmp and AtWb_tmp to three-dimensional arrays
    real(kind=8 ), allocatable :: normeq_tmp(:,:,:)    ! The left side of the LS algorithm. P_t*W*P     tmp for openMP implementation
    real(kind=8 ), allocatable :: AtWb_tmp(:,:)        ! The right side of the LS algorithm. P_t*W*OMC  tmp for openMP implementation
    integer*4 count_epoch_loop

! PT201006: additional variables from Seb's mods to gracefit, to use less memory and speed up the stacking of normal equations
!  DOUBLE PRECISION, ALLOCATABLE , target:: part(:,:,:)    ! Partials of observables with respect to the parameters
!  double precision, dimension(:,:), pointer :: ptr_part
  DOUBLE PRECISION, ALLOCATABLE :: normeq_tmp2(:,:)    ! The left side of the LS algorithm. P_t*W*P     tmp for openMP implementation
  DOUBLE PRECISION, ALLOCATABLE :: AtWb_tmp2(:)    
  integer :: ompthread, ompmax
  INTEGER, EXTERNAL :: omp_get_thread_num, omp_get_max_threads
  integer*4 :: k

! PT140617: last epoch that we want to process (hard-wired for now, but should set through the command file)
    integer*4 :: last_epoch

! PT131104: variables for reading binary, direct access GTORB files
    integer*4  :: irec     ! record number

  ! PT170725: variables used in removing edge effects introduced when filtering prefit residuals
  INTEGER*4    :: edge1(100),edge2(100)    ! list of epochs of start and stop of gaps in kbr obs
  INTEGER*4    :: n_kbr_gaps               ! number of data gaps in the kbr observations
  LOGICAL      :: in_kbr_gap               ! flag for whether in a kbr data gap or not
  logical      :: invoke_NR_filt           ! flag to invoke (or not) the noise_robust_derivative filter

! PT170713: variables for reading the VCV file to enable the computation of postfit residuals
    logical      :: VCV_flag
    integer*4    :: nparam_vcv                                          ! number of parameters found in the VCV file
    integer*4    :: sTid,sScl, sBias, sEmp, sMasc, n_emp, nScl, nBias   ! we don't need these for the postfit residual computations
!    real(kind=8),allocatable :: soln(:)                                 ! solution vector from the VCV file
!    real(kind=8),allocatable :: apriori(:)
    character*30,allocatable :: prm_input(:)

! PT140416: debug partials for gracefit but needed in argument list here too !
    real(kind=8), allocatable :: kbrr_part(:,:,:),kbra_part(:,:,:)
  real (kind=8), allocatable :: Amat(:,:)
! PT140818: debug kbrr noise variables
    real(kind=8) :: kbrr_noise,noise_old

! RM210416: variables to add postfit to "simulate" noise
  character                  :: KBresid_file*100
  logical                    :: add_postfit = .false.
  real(kind=8)               :: KBresid(38),KBAresid(38)

! PT180125: variables used to determine whether GTORB files are binary or hd5
  LOGICAL :: is_h5
  CHARACTER(len=80) :: test
  INTEGER :: length

! PT200102: parameter variables
  integer*4 :: nmascons_t_truth, nmascons_t_apriori,nrvec_t_truth,nrvec_t_apriori

! RM191112: counter for counting every call to calc_beta_angle
  integer*4    :: nbeta

! PT210131: variable to store the first header line of the regularisation file
  character*70 :: header_line

! PT201109: 5-element integer to store yr/mo/day/hr/min
  integer*4    :: date(5)


!****************************************************************

! If there are no arguments, give instructions as to program usage
    if( iargc() == 0 ) then
        write(*,10)
        10 format(/,'###################################################################################################' &
        ,/,'                                  Program GRACESIM(new)                         ' &
        ,/,'Purpose: Fit integrated GRACE GTORB orbits to imported "truth" GTORB files' &
        ,/,'Usage: GRACESIM cmd-file rms-file plt-option GTORB_A GTORB_B GTORB_truth_A GTORB_truth_B' &
        ,' YYYY MM DD 0 0 0 0 0 0' &
        ,/ &
        ,/,'Example: GRACESIM_new GRACESIM_new.cmd thomas_new.rms 4 GTORB_2010-09-14_00-00-00_A_02.asc_iter3' &
        ,' GTORB_2010-09-14_00-00-00_B_02.asc_iter3 GTORB_2010-09-14_00-00-00_A_02.asc_truth ' &
        ,' GTORB_2010-09-14_00-00-00_B_02.asc_truth 2010 09 14 0 0 0 0 0 0' &
        ,/,'###################################################################################################')

        write(*,11)
        11 format(/,' a template GRACESIM.cmd file format is given in ~/gg/grace/tables/GRACESIM.cmd.template:',/, &
        'Simply add your preferred GPS and KBRR observation sigmas and you are ready to go.',/, &
        'As Tom Herring says "with perfect data, this should be easy !"',/)
        stop
    endif
!****************************************************************
    calling_prog = "GRACESIM"

! Get version number and write out status line
print*,'calling gversn'
    call gversn(version)
print*,'called it'
    write(*,'(a)')' '
    write(message,'(a,a120)') 'Started GRACESIM ',version
!    call status_update('STATUS',calling_prog,'gracesim/gracesim',' ',message,0)

! Assert that there are the correct number of command-line arguments
    if ( iargc() < NUM_EXPECTED_ARGS )then
        call status_update('FATAL',calling_prog,'gracesim',' ','Too few arguments in command line',0)
    endif

!************** PARSE COMMAND LINE AND COMMAND FILE *************

! Now read in the command-line and command file
    call getarg(1,cmdnam)      ! Command file
    call getarg(2,rmsnam)      ! rms file name
    call getarg(3,arg)         ! Ploting indicator. Placed into iop
    read(arg,'(i1)') iop
    call getarg(4,gt_fnam(1))  ! Integrated GRACE GTORB orbits for satellites A and B
    call getarg(5,gt_fnam(2))
    call getarg(6,truth_fnam(1))  ! Integrated GRACE GTORB orbits for satellites A and B
    call getarg(7,truth_fnam(2))
    call getarg(8,c_yr)        ! Year
    call getarg(9,c_month)     ! Month
    call getarg(10,c_day)       ! Day
    call getarg(11,msc_const_file) ! mascon constraint file name
    call getarg(12,arg) ! RM210416: option to add postfit residuals to prefit to simulate noise
    read(arg,'(a)') KBresid_file

! Get instructions
    call input_openCommand(cmdnam)
! PT221122: allocate the prmnam array. Was done in soln_mod but now needs to be done in calling program
    allocate(prmnam(prmnam_size))
    call command_readSetup()   ! PT191227: this contains "use soln_mod" which defines prmnam

! PT180125: determine whether files are hdf5 format or regular fortran binary format. For now, all files
!           must be of the same format!
  test = gt_fnam(1)
  length = LEN(TRIM(gt_fnam(1)))
  IF( test(length-2:length) == '.h5') THEN
     call status_update('STATUS','gracesim','gracesim',gt_fnam(1),'Input GTORB files are hdf5 format',0)
     is_h5 = .TRUE.
  ELSE
     call status_update('STATUS','gracesim','gracesim',gt_fnam(1),'Input GTORB files are NOT hdf5 format (assume bin format)',0)
     is_h5 = .FALSE.
  ENDIF

! PT131015: open the two "truth" GTORB files. They are given unit numbers LUGT93,4) which are 102, 103
  call status_update('STATUS',calling_prog,'gracesim',truth_fnam(1),'Opening truth GTORB files',0)
  IF (is_h5 ) THEN
     CALL input_openh5GTORB(truth_fnam,2)
     CALL header_readh5GTORB(irec,satics_t,apr_ScaleBias,GPSantoff,mcon_tides &
          ,combined_mascon_file,msc_hdr_code,ocean_mascon_file,msc_ocn_hdr_code,total_ocean_prim,2,prmnam,prmnam_size)
  ELSE
    call status_update('FATAL',calling_prog,'gracesim',' ',"Binary GTORB files no longer supported. Use .h5 instead",0)
  ENDIF

! PT200102: allocate variables required to read the truth GTORB files
  allocate(truthvec(6,nsat_t,nepochs_t)) !PT201006: Seb changed this !from nrvec_t to 6
  allocate(apr_prm(nparam_t))
  allocate(apr_wght(nobs_t))
  allocate(apr_const(nparam_t))
  allocate(sciframe_quat(4,nsat_t,nepochs_t))
  allocate(srf2trf_rotmat(3,3,nepochs_t,nsat_t))
  allocate(srf2trf_deriv_rotmat(3,3,nepochs_t,nsat_t))
  allocate(mcon_ocean(nmascons_t))
  allocate(msc_tide_amp(max_mcon_tides,2,nmascons_t))
! PT201013: replace "part" with the new Amat matrix
!  ALLOCATE(part(nobs_t,nparam_t,nepochs_t))
  allocate(Amat(nobs_t*nepochs_t,nparam_t))

! now read the truth GTORB files
  call status_update('STATUS','GRACESIM','gracesim',' ','Reading truth GTORB files',0)
  IF (is_h5 )THEN
   CALL input_readGTORBh5_Amat(satics_t,apr_ScaleBias,GPSantoff,truthvec,Amat,apr_prm &
          ,sciframe_quat,srf2trf_rotmat,srf2trf_deriv_rotmat,2, mcon_ocean &
          ,msc_tide_amp, mcon_tides)
  ELSE
    call status_update('FATAL','GRACESIM','gracesim',' ','code no longer works with binary fortran files. Must use h5 format',0)
  ENDIF

! and now deallocate these variables since we don't need them for the truth orbits (re-allocate below what is needed)
  deallocate(apr_prm)
  deallocate(apr_wght)
  deallocate(apr_const)
  deallocate(sciframe_quat)
  deallocate(srf2trf_rotmat)
  deallocate(srf2trf_deriv_rotmat)
  deallocate(mcon_ocean)
  deallocate(msc_tide_amp)


! PT200102: save the number of mascons and dimensioning for truthvec for the truth orbits
  nmascons_t_truth = nmascons_t
  nrvec_t_truth    = nrvec_t

! PT131017: reset nmascons_t and nparam_t so the code works properly for reading the next GTORB files
! PT140707: also need to subtract off the number of tidal amplitude parameters to reset this
    if(nmascons_t > 0 )then
      if(est_msc_tides > 0) then
        nparam_t = nparam_t - nmascons_t - nmsc_tid_constit_t*2
      else
        nparam_t = nparam_t - nmascons_t 
      endif
      nmascons_t = 1
    endif
    
! Open and read header of perturbed GTORB files
    call status_update('STATUS','GRACESIM','gracesim',gt_fnam(1),'Opening modelled GTORB files',0)
  IF (is_h5 ) THEN
     CALL input_openh5GTORB(gt_fnam,0)
     CALL header_readh5GTORB(irec,satics_t,apr_ScaleBias,GPSantoff,mcon_tides &
          ,combined_mascon_file,msc_hdr_code,ocean_mascon_file,msc_ocn_hdr_code,total_ocean_prim,0,prmnam,prmnam_size)
  ELSE
    call status_update('FATAL',calling_prog,'gracesim',' ',"Binary GTORB files no longer supported. Use .h5 instead",0)
  ENDIF
! PT200102: save the number of mascons and dimensioning for truthvec for the apriori orbits
  nmascons_t_apriori = nmascons_t
  nrvec_t_apriori    = nrvec_t

! PT140901: set end_neq to be nepochs_t if it wasn't set in the command file
  if(end_neq <= 0)end_neq = nepochs_t 
  call command_printSetup() ! Must come after reading GTORB file(s)

!************************ ALLOCATE ARRAYS ***********************

  call status_update('STATUS','GRACESIM','gracesim',' ','Allocating arrays',0)
! Variables storing data from command file and GTORB file
  allocate(apr_prm(nparam_t))
! PT201014: apr_wght now needs dimension of nobs_t*nepochs*t (not just nobs_t)
  ALLOCATE(apr_wght(nobs_t*nepochs_t))
  allocate(apr_const(nparam_t))
! PT150826: add tidal amplitude constraint variable
  allocate(apr_tide_ampl_const(max_mcon_tides,2))
! Variables used for calculations
  allocate(normeq(nparam_t,nparam_t))
  allocate(AtWb(nparam_t))
  allocate(adjust(nparam_t))
  ALLOCATE(pre_omc(nobs_t,nepochs_t))
  ALLOCATE(pre_omc_kb_unfilt(2,nepochs_t))   ! backup copy of unfiltered prefit KBRR/KBRA residuals. Used to create unfiltered postfit kbrr residuals for .kb file
  ALLOCATE(post_omc_kb_unfilt(2,nepochs_t))  ! backup copy of unfiltered prefit KBRR/KBRA residuals. Used to create unfiltered postfit kbrr residuals for .kb file
  ALLOCATE(post_omc(nobs_t,nepochs_t))
  allocate(bvec(nobs_t*nepochs_t))
! PT201014: array to indicate whether to use observations or not
  allocate(use_obs(nepochs_t,6))   ! 1: kbr; 2: kbrr; 3: kbra; 4: LR; 5: LRR; 6: LRA  (first 3: kbr, last 3: LRI)
  allocate(use_gvec(nepochs_t,2))  ! 1: kbr; 2: kbrr; 3: kbra; 4: LR; 5: LRR; 6: LRA  (first 3: kbr, last 3: LRI)
  use_obs(:,:) = .true.  ! PT201014: initialise the use of all inter-satellite observations (kbr and LRI) to be false
! Data read in or calculated from GTORB file(s)
  allocate(rvec(6,nsat_t,nepochs_t))                        !PT201007: redimensioned to (6,nsat_t,nepochs_t) to match Seb's changes
  allocate(sciframe_quat(4,nsat_t,nepochs_t))
  allocate(srf2trf_rotmat(3,3,nepochs_t,nsat_t))
  allocate(srf2trf_deriv_rotmat(3,3,nepochs_t,nsat_t))
! Variables used for rpy calculations and antenna offsets
  allocate(rpy(3,nsat_t,nepochs_t))
  allocate(AOC(9,3,nepochs_t))
  allocate(GPS_antoff_trf(6,nsat_t))
  allocate(GPSant_adj(3,nsat_t))
! Data read in from ACC1B file(s)
  allocate(acc(nepochs_t,4,nsat_t))
! Data read in from GNV1B files
  allocate(gvec(maxgprm,nsat_t,nepochs_t))
  allocate(num_gps_used(nsat_t))
!    allocate(uvec(maxgprm,nsat_t,nepochs_t))
! Data read in from KBR1B file
  allocate(kbrange(3,nepochs_t))
  allocate(kblt(3,nepochs_t))
  allocate(kbant(3,nepochs_t))
  allocate(kbrange_L1B(3,nepochs_t))
  allocate(kblt_L1B(3,nepochs_t))
  allocate(kbant_L1B(3,nepochs_t))
!    allocate(kbion_cor(nepochs_t))
!    allocate(kbsnr(2,nsat_t,nepochs_t))
! Data read in or calculated from THR1B file(s)
  allocate(thr(nepochs_t,16,nsat_t))
  allocate(thrusters_off(nepochs_t,nsat_t))
! Variables used for application of shadow condition
  allocate(shadow_counter(nsat_t))
  allocate(apply_shadow_cond(nepochs_t,nsat_t))
  allocate(shadow(nsat_t,nepochs_t))
! filtered omc for the kbrr and kbra observations
  allocate(filt_kbrr_omc(nepochs_t))
  allocate(filt_kbra_omc(nepochs_t))
  allocate(kb_1pr(nepochs_t))
! beta angle
  allocate(beta_angle(nepochs_t))
! mascon length-dependent constraints
  allocate(msc_const(nmascons_t,nmascons_t))
! mascon ocean/land identifier
  allocate(mcon_ocean(nmascons_t))
! mascon tidal amplitudes
  allocate(msc_tide_amp(max_mcon_tides,2,nmascons_t))
! partials of pos/vel wrt kbrr and kbra
  allocate(kbrr_part(6,2,nepochs_t))
  allocate(kbra_part(6,2,nepochs_t))
! variables for reading the VCV file
! PT210130: allocated further down when needed
!  allocate(soln(nparam_t))
!  allocate(apriori(nparam_t))
  allocate(prm_input(nparam_t))

!****************************************************************
! Read in the rest of the command file and command-line
  call command_readAprConst(apr_const,apr_tide_ampl_const)
  call command_readDataWeights(apr_wght)
  call command_readTol(kbrr_prefit_tol)
  call command_close()

! Write out date_str to be used for finding data files
  write(date_str,'(a4,a1,a2,a1,a2)')c_yr,"-",c_month,"-",c_day

! Read in the adjustments of CoM/GPSant for each sat from command line
! PT131015: increment this so that it reads from 11 onwards
! PT150812: turn this off !
!    do isat= 1, nsat_t
!      do i=1,3
!        call getarg(10+i+(isat-1)*3,arg)
!        read(arg,*)GPSant_adj(i,isat)
!      enddo
!    enddo
!****************************************************************
   GPSant_adj = 0.0 
! Open input and output files
  call input_openAllFiles(date_str)
  call output_openAllFiles(rmsnam)
! PT131015: read first the "truth" GTORB files and store the position/velocity information in gvec. The rest (scales, quaternions etc)
!           will be over-written by the subsequent call for reading the theoretical GTORB files (which is what we want to have happen)
! PT200102: now moved higher, so comment out here
!  call status_update('STATUS','GRACESIM','gracesim',' ','Reading truth GTORB files',0)
!  nmascons_t = nmascons_t_truth
!  nrvec_t = nrvec_t_truth
!  IF (is_h5 )THEN
!print*,'nparam_t,nmascons_t,nrvec_t',nparam_t,nmascons_t,nrvec_t
!   CALL input_readGTORBh5(satics_t,apr_ScaleBias,GPSantoff,truthvec,apr_prm &
!          ,sciframe_quat,srf2trf_rotmat,srf2trf_deriv_rotmat,2, mcon_ocean &
!          ,msc_tide_amp, mcon_tides)
!print*,'have read truth orbits'
!  ELSE
!    call input_readGTORB("direct   ",irec,satics_t,apr_ScaleBias,GPSantoff,truthvec,apr_prm &
!                          ,sciframe_quat,srf2trf_rotmat,srf2trf_deriv_rotmat,2,mcon_ocean,msc_tide_amp,mcon_tides)
!  ENDIF


!*********************************************************
!*********************************************************
!*                                                      **
!*  Compute the simulated inter-satellite observations  ** 
!*                                                      **
!*********************************************************
!*********************************************************
! PT131015: compute the range rate and store it and pos/vel in kbrange, gvec
  nrvec_t = nrvec_t_truth
  call kb_computeSim(truthvec,gvec,kbrange)
  nrvec_t = nrvec_t_apriori
! PT171129: need to generate a simulated range acceleration observation through differentiation of the range
  call kb_differentiate_kbr(.false.,kbrange,kblt,kbant)    ! .false. fixed spikes maybe? HM,SA 181215
! PT131015: Set kblt and kbant to zero, since our kbrange is derived from positions and velocities of the CoM
  kblt = 0.d0
  kbant = 0.d0
! PT210130: need to set up the gvec_epochs array to mimic having read when there are GNV1B epochs ...
  allocate( gvec_epochs(nepochs_t,2))
  do iepoch = 1,nepochs_t
    gvec_epochs(iepoch,:) = starting_epoch + (iepoch-1)*epoch_interval
  enddo
!****************************************************************


!*********************************************************
!*********************************************************
!*  PT181120                                            **
!*  read actual inter-satellite observations available  ** 
!*  Hardwired for GRACE RL02 data, just for one test    **
!*********************************************************
!*********************************************************
  use_kbr_epochs = .true.
  use_kbr_epochs = .false.
  if(use_kbr_epochs) then
! now we want to open an actual KBR1B file and read the actual epochs of KBR1B data for this day
    write(kb_fnam,'(a,a,a)')"KBR1B_",date_str,"_X_02.asc"
    call status_update('STATUS',calling_prog,'gracesim',kb_fnam,"Reading file to identify actual KBR epochs",0)
    call kbr_read_hdr(LUKB,calling_prog,kb_fnam,0,n_kbr)

! and read the actual data. Here, we pass different arrays to store the actual data so that they are not overwritten by the
! simulated observations later on. 
    call kbr_read_data(LUKB,calling_prog,kb_fnam,mission,nepochs_t,starting_epoch,epoch_interval &
                           ,kbrange_L1B,kblt_L1B,kbant_L1B,num_epochs)

! finally, setto zero any KBR simulated observations if there wasn't a measurement
    do iepoch = 1,nepochs_t
      if(kbrange_L1B(1,iepoch) == 0.d0)then
        write(message,'(a,i5)')"KBR data set to missing at epoch: ",iepoch
    !    call status_update('STATUS',calling_prog,'gracesim',' ',message,0)
        kbrange(1,iepoch) = 0.d0
      endif
    enddo

  endif
!****************************************************************


!****************************************************************
! Read in all needed data from input files.
  call status_update('STATUS','GRACESIM','gracesim',' ',"Reading GTORB files",0)
  nmascons_t = nmascons_t_apriori
  nrvec_t = nrvec_t_apriori
  IF (is_h5 )THEN
   ! PT201006: changed to v1 of subroutine and added "part" to passed arguments
   ! PT210130: updated to the "Amat" version here
   CALL input_readGTORBh5_Amat(satics_t,apr_ScaleBias,GPSantoff,rvec,Amat,apr_prm &
          ,sciframe_quat,srf2trf_rotmat,srf2trf_deriv_rotmat,0, mcon_ocean &
          ,msc_tide_amp, mcon_tides)


  ELSE
    call input_readGTORB("direct   ",irec,satics_t,apr_ScaleBias,GPSantoff,rvec,apr_prm &
                        ,sciframe_quat,srf2trf_rotmat,srf2trf_deriv_rotmat,0,mcon_ocean,msc_tide_amp,mcon_tides)
  ENDIF

! PT201006: the gracefit code now calls something to compute the kb partials, so add it here too.
  call status_update('STATUS','GRACESIM','gracesim',' ',"Calling kb_computePartials_Amat",0)
  call kb_computePartials_Amat(rvec, Amat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the L1B accelerometer and thrust data. Thrust data must be read first, since it is used in the fitting
! of the model to linearise the accelerometer data (inside input_readACC)
  call input_readACC1B(acc)
  call input_readTHR1B(thr,thrusters_off)
!****************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT161128: move the opening of the regularization file up higher in the program
  if(use_msc_len_const == 1)then
    open(243,file=msc_const_file(1:trimlen(msc_const_file)),status='old',iostat=ioerr)
    if(ioerr /= 0)then
      call status_update('FATAL','GRACESIM','gracesim',msc_const_file(1:trimlen(msc_const_file)) &
                           ,'Error opening msc regularization',0)
    endif
! PT210131: read the first header line of the regularisation file
    read(243,'(a)')header_line
! PT170206: verify that the mascon file code is compatible with that of the GTORB files
    read(header_line,'(a8)')msc_reg_code
! PT170713: check whether it is a VCV file to be used for postfit residual computations
    if(msc_reg_code(1:1) == "V")then
      VCV_flag = .true.
      call status_update('STATUS','GRACESIM','gracesim',msc_const_file &
                ,'Will compute postfit residuals using solution in VCV file',0)

    else if(msc_reg_code /= msc_hdr_code)then
      write(message,'(5a)')"Incompatible temporal mascon file (code ",msc_hdr_code,") and mascon regularization file (code " &
              ,msc_reg_code,"). Should not continue"
      call status_update('FATAL','GRACESIM','gracesim',msc_const_file,message,0)
    endif
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 

!****************************************************************

! precompute the yaw, pitch, roll angles and the antenna offset correction (AOC) of the KBR
  if(nsat_t == 2) call rpy_AOC(rvec(1:3,1,:),rvec(1:3,nsat_t,:),sciframe_quat(:,1,:),sciframe_quat(:,nsat_t,:),rpy,AOC) !@# condition...
! PT131015: for gracesim, we just make AOC = 0
  AOC = 0.d0

! Initialize necessary arrays
  call status_update('STATUS','GRACESIM','gracesim',' ','Initialising arrays',0)
  pre_omc = 0.d0
! PT201009: this should not be initialised here - it already contains the partials!
!  part    = 0.d0
  kbrr_misfit_num = 0
  shadow_counter = 0
  apply_shadow_cond = .false.

!********************************************************************************************************************************

!****************************************************************
!*  PT201109: read outlier removal file and flag unwanted epochs/obs      *
  read(c_yr,*)date(1)
  read(c_month,*)date(2)
  read(c_day,*)date(3)
  call mask_outliers(calling_prog,use_outlier_mask,nepochs_t,date(1:3),use_obs)
! DEBUG
!do iepoch=950,1050
!  print*,'iepoch, use_obs:',iepoch,use_obs(iepoch,:)
!enddo
!stop 'stop after mask_outlier'
!****************************************************************


!****************************************************************
!*                                                              *
!******                     EPOCH LOOP                  *********
!*                                                              *
!****************************************************************
! Report beginning of loop
  call status_update('STATUS','GRACESIM','gracesim',' ','Loop over epochs to prepare partials',0)
  last_epoch = nepochs_t

! Epoch loop
    do iepoch = 1, last_epoch ! nepochs_t

!******************** CALCULATIONS FOR KBAND ********************

! Only do calculations if there was a line at corresponding epoch (checked by making sure biased range is not 0.0)
      if (kbrange(iRANGE,iepoch) /= 0.d0) then
! Calculate pre_omc
        call kb_computeOMC(iepoch,pre_omc(:,iepoch),rvec(1:nrvecprm_t,:,iepoch),kbrange(:,iepoch),kblt(:,iepoch),kbant(:,iepoch),&
                           aoc(:,:,iepoch))
! PT140818: add 0.1 um/s random noise to the omc
       if(iepoch == 1) then
         kbrr_noise = (rand(0)-0.5d0)*1.d-7
         noise_old = kbrr_noise
       else
         kbrr_noise = noise_old * exp(-5/500.d0) + (rand(0)-0.5d0)*1.d-7
         noise_old = kbrr_noise
       endif
!       kbrr_noise = (rand(0)-0.5d0)*1.d-7
! PT171129: over-ride this and set the noise to zero
       kbrr_noise = 0.d0
       pre_omc(ikbrr,iepoch) = pre_omc(ikbrr,iepoch) + kbrr_noise

! DEBUG:
!print*,'pre_omc kbrr = ',iepoch,pre_omc(ikbrr,iepoch),nrvecprm_t


! Compute and place necessary partials into part matrix
! PT140511: pass back from kb_computePartials the dRR/dX, dRR/dY etc in new variable kbrr_part(6,2). Use this later in debug stuff.
! PT201006: I think the kb partials have now been computed higher up with kb_computePartials_vec
!        call kb_computePartials(iepoch,part(:,:,iepoch),rvec(:,:,iepoch),kbrr_part(:,:,iepoch),kbra_part(:,:,iepoch))

! Validate that the values calculated in omc are within our prefit tollerance !@# Make sure this does what we want it to do
! PT210129: turn this off - we don't need it for simulated data and it is now incompatible with Amat vs part code!
        !if(nkobs_t > 0)call kb_validateRR(part(ikbrr,:,iepoch),pre_omc(ikbrr,iepoch),kbrr_misfit_num,kbrr_prefit_tol,iepoch)
      endif
!****************************************************************
  
!********************* CALCULATIONS FOR GOBS ********************
     ! PT180730: need logic here to match the GTORB epoch (5sec ?) with GNV1B epochs (1s or 5s). Check each satellite separately
     do isatt = 1,nsat_t
       igvec =  1

       if(starting_epoch + (iepoch-1)*epoch_interval > gvec_epochs(igvec,isatt) ) then
         do while (gvec_epochs(igvec,isatt) < starting_epoch + (iepoch-1)*epoch_interval)
           igvec = igvec + 1
         enddo
         if(gvec_epochs(igvec,isatt) == starting_epoch + (iepoch-1)*epoch_interval)then
           use_gvec(iepoch, isatt) = .true.
         else
           use_gvec(iepoch, isatt) = .false.
         endif
       endif

       if( use_gvec(iepoch, isatt) ) then 
         ! Compute GPS antenna offset in terrestrial reference frame
         call quat_rot_vect(sciframe_quat(:,isatt,iepoch),GPSantoff(5:7,isatt),GPS_antoff_trf(1:3,isatt))
         call matmult(srf2trf_deriv_rotmat(:,:,iepoch,isatt),GPSant_adj(:,isatt),GPS_antoff_trf(4:6,isatt),3,3,1)

         ! Calculate pre_omc
         call gobs_computeOMC(pre_omc(:,iepoch),gvec(:,:,igvec),rvec(1:ngobs_t,:,iepoch),GPS_antoff_trf)
  
         ! Apply CoM conditions
         !CALL add_CoM_cond(part(:,:,iepoch),pre_omc(:,iepoch),rvec(1:nCoM_t,:,iepoch),gvec(:,:,igvec),GPS_antoff_trf,&
         !    srf2trf_rotmat(:,:,iepoch,:),srf2trf_deriv_rotmat(:,:,iepoch,:))

         ! Compute the beta angle
         !nbeta = nbeta+1 ! RM191112
         !CALL calc_beta_angle(iepoch,rvec(1:3,1,iepoch),rvec(4:6,1,iepoch),beta_angle(iepoch))
       endif

     enddo ! end of satellite do loop for GNV1B obs
  ENDDO  ! End of epoch loop
!********************************************************************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RM210416: determine whether we add postfit as noise to the prefit
  if(KBresid_file(:1) .ne. "0")then
    add_postfit = .true.
    call status_update('STATUS','GRACESIM','gracesim',KBresid_file,'Will add kb and kba postfit residuals to simulate noise',0)
  endif

! RM210416: read postfits from *.kb and *.kba file
  if(add_postfit)then
    !allocate(KBresid(nepochs_t,38))
    open(245,file=KBresid_file(1:trimlen(KBresid_file)),status='old',iostat=ioerr)
    if(ioerr /= 0)then
      call status_update('FATAL','GRACESIM','gracesim',KBresid_file(1:trimlen(msc_const_file)) &
                           ,'Error opening k-band residuals file',0)
    endif
    open(246,file=KBresid_file(1:trimlen(KBresid_file))//"a",status='old',iostat=ioerr)
    if(ioerr /= 0)then
      call status_update('FATAL','GRACESIM','gracesim',KBresid_file(1:trimlen(msc_const_file))//"a" &
                           ,'Error opening k-band residuals file',0)
    endif
    ! RM210416: skip header lines
    do i=1,6
      read(245,*)
      read(246,*)
    enddo
    do iepoch=1,nepochs_t
      read(245,*)KBresid(:)
      read(246,*)KBAresid(:)
      pre_omc(1:12,iepoch) = pre_omc(1:12,iepoch) + KBresid(25:36)
      pre_omc(ikbrr,iepoch) = pre_omc(ikbrr,iepoch) + (KBresid(2)/m_mm) ! these are filtered
      pre_omc(ikbra,iepoch) = pre_omc(ikbra,iepoch) + (KBAresid(2)/m_mm/m_mm) ! these are filtered and nonzero to being with?
      !pre_omc(1:12,iepoch) = pre_omc(1:12,iepoch) + KBresid(25:36)
      !pre_omc(ikbrr,iepoch) = pre_omc(ikbrr,iepoch) + (KBresid(37)/m_mm) ! these are unfiltered
      !pre_omc(ikbra,iepoch) = pre_omc(ikbra,iepoch) + (KBAresid(37)/m_mm/m_mm) ! these are unfiltered and nonzero to being with?

    enddo
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!******************* RANGE ACCELERATION COMPUTATIONS ******************
! PT170428: in this subroutine we will numerically differentiate the kbrr theoretical range
!           and partials to compute the information for the range accelerations
  if (use_kb_nr(2) > 0 .and. nkbra_t == 1) then
    nrvec_t = 6
    CALL kb_computeKBRA(use_obs,use_kb_nr(2),ikbrr,ikbra,nepochs_t,nobs_t,nparam_t,nrvec_t,nsat_t,nrvecprm_t,pre_omc,Amat,rvec &
       ,kbrr_part,kbra_part,n_kbr_gaps,kbrange, edge1,edge2)
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call status_update('STATUS','GRACESIM',calling_prog,'',"End of epoch loop that prepared the partial derivatives",0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! DEBUG
! PT210908: output the partial derivatives d_kbra/d_msc for a selection of mascons to assess what is going on
!do iepoch = 1,nepochs_t
!  print*,iepoch,Amat(ikbra+(iepoch-1)*nobs_t,24+4300:24+4581),' kbra partials'
!enddo
!stop 'stopped after partials debug'

!*********** ADD CONDITONS AND SET UP NORMAL EQUATIONS **********

! Condition on mean acceleration of the two satellites
! PT210131: turn this off. Not compatible with the "Amat" version
!    if(naobs_t > 0) call add_meanACC_cond(acc,apr_wght,apr_ScaleBias,part(:,:,1),pre_omc(:,1))

  !*********** FILTER KBRR/KBRA obs and/or remove a 1pr ?  **********
  IF(override_kbr_filt .AND. INT(use_kbrr_decomp(1)) == -2)THEN
     CALL status_update('STATUS','GRACESIM','gracesim',' ','Filtering of kbr obs will NOT occur because of missing obs',0)
  ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PT170717:   filtered KBRR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! here we filter the kbr obs. use_kbrr_decomp(1) = -2 in the case where we want to filter but only if there are no data gaps
  override_kbr_filt = .FALSE.
  IF( (INT(use_kbrr_decomp(1)) > 0 .OR. INT(use_kbrr_decomp(1)) == -2) .AND. .NOT. override_kbr_filt)THEN  ! we remove the high-frequency noise from all options apart from 0 (leave as is)
     ! DEBUG
     !print*,'pre_omc(12450:12550):'
     !do i=12450,12550
     !  print*,i,pre_omc(ikbrr,i)
     !enddo
     !stop
     IF(nkbrr_t > 0) then
      CALL kb_filtered(start_neq,end_neq,NINT(use_kbrr_decomp(2)), kbrange(iRANGE,:),nobs_t,pre_omc(iKBRR,:),filt_kbrr_omc)
    endif
     IF(nkbra_t > 0) then
      CALL kb_filtered(start_neq,end_neq,NINT(use_kbrr_decomp(2)), kbrange(iRANGE,:),nobs_t,pre_omc(iKBRA,:),filt_kbra_omc)
    endif

     write(1000,*) pre_omc(iKBRR,:)
     write(1001,*) filt_kbrr_omc
     !stop 'stopped after kb_filtered'

     ! generate a model with a 1/rev signal with linearly changing amplitude removed ?
     kb_1pr = 0.d0
     IF(use_kbrr_decomp(1) == 3)THEN
        CALL kb_sponge(nepochs_t,filt_kbrr_omc(:),kb_1pr)  ! remove kbrr sponge as well
     ELSE
        CALL status_update('STATUS','GRACESIM','gracefit',' ','Not removing kbrr empirical model',0)
        kb_1pr = 0.d0
     ENDIF

     ! replace the kbrr prefit omc with the filtered one, and with the 1/rev model removed
     IF(nkbrr_t >0 .AND. INT(use_kbrr_decomp(1)) > 0 .AND. .NOT. override_kbr_filt)THEN
     DO i=1, nepochs_t
           pre_omc(ikbrr,i) = filt_kbrr_omc(i) - kb_1pr(i)
     ENDDO
   ENDIF

     ! write out filtered PREFIT RESIDUALS
     CALL status_update('STATUS','GRACESIM','gracefit',' ','Writing out filtered prefit residuals',0)
     OPEN(397,file='prefit_sim.res',status='unknown')
     WRITE(397,'(a,a)')" iepoch :   kbrr  in mm/s  :    latitude     longitude  :"  &
          ,":     XYZ (GRACE A) in m      :      XYZvel (GRACE A) in mm/s      " &
          ,":     XYZ (GRACE B) in m  :      XYZvel (GRACE B)   in mm/s "
     DO iepoch = 1,nepochs_t
        satlon = datan2(rvec(2,1,iepoch),rvec(1,1,iepoch)) *180.d0/pi
        satlat = dasin(rvec(3,1,iepoch)/dsqrt(rvec(1,1,iepoch)**2 + rvec(2,1,iepoch)**2 + rvec(3,1,iepoch)**2))*180.d0/pi
        WRITE(397,'(i6,3f14.8,13f12.5)')iepoch,pre_omc(iKBRR,iepoch)*1.d3-kb_1pr(iepoch)*1.d3 &
             ,satlat, satlon &
             ,pre_omc(1:3,iepoch) &
             ,pre_omc(4:6,iepoch)*1.d3,pre_omc(7:9,iepoch),pre_omc(10:12,iepoch)*1.d3,kb_1pr(iepoch)*1.d6
     ENDDO
     CLOSE(397)


     ! PT170428: write out a kbra prefit file as well, if using acceleration observations
     IF(nkbra_t > 0)THEN
        OPEN(397,file='prefitRA_filt_sim.res',status='unknown')
        WRITE(397,'(a,a)')" iepoch :   kbrr  in mm/s  :    latitude     longitude  :"  &
             ,":     XYZ (GRACE A) in m      :      XYZvel (GRACE A) in mm/s      " &
             ,":     XYZ (GRACE B) in m  :      XYZvel (GRACE B)   in mm/s "
        DO iepoch = 1,nepochs_t
           satlon = datan2(rvec(2,1,iepoch),rvec(1,1,iepoch)) *180.d0/pi
           satlat = dasin(rvec(3,1,iepoch)/dsqrt(rvec(1,1,iepoch)**2 + rvec(2,1,iepoch)**2 + rvec(3,1,iepoch)**2))*180.d0/pi
           ! filtered kbra output
           WRITE(397,'(i6,3f14.8,13f12.5,f18.8)')iepoch,filt_kbra_omc(iepoch)*1.d6  &
                ,satlat, satlon &
                ,pre_omc(1:3,iepoch) &
                ,pre_omc(4:6,iepoch)*1.d3,pre_omc(7:9,iepoch),pre_omc(10:12,iepoch)*1.d3,kb_1pr(iepoch)*1.d6 &
                ,pre_omc(iKBRA,iepoch)*1.d6
           ! filtered kbra output

        ENDDO
        CLOSE(397)
     ENDIF


     !call kb_computeKBRA(ikbrr,ikbra,nepochs_t,nobs_t,nparam_t,nrvec_t,nsat_t,nrvecprm_t,pre_omc,part,rvec,kbrr_part,kbra_part &
     !   ,n_kbr_gaps,edge1,edge2)
     !if(nkbra_t >0) call kb_filtered(start_neq,end_neq, kbrange(iRANGE,:),nobs_t,pre_omc(iKBRA,:),filt_kbra_omc)


     ! replace the kbra prefit omc with the filtered one, and with the 1/rev model removed
     !FIXME: SA Shouldn't be beofre the writting part?
     IF(nkbra_t >0  ) THEN
     DO i=1, nepochs_t
         pre_omc(ikbra,i) = filt_kbra_omc(i)
     ENDDO
   ENDIF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! PT170717: unfiltered KBRR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! PT170717: also use unfiltered obs if the filter option of -2 was chosen (i.e. filter but not if there are data gaps)
  ELSE IF ( INT(use_kbrr_decomp(1)) == 0 .OR. (INT(use_kbrr_decomp(1)) == -2 .AND. override_kbr_filt) ) THEN

     kb_1pr = 0.d0

     ! write out the prefit residuals at each epoch, so that I don't have to run gracefit all the way through just to get this information
     CALL status_update('STATUS','GRACESIM','gracefit',' ','Writing out (unfiltered) prefit range rate residuals',0)
     OPEN(397,file='prefit_sim.res',status='unknown')
     WRITE(397,'(a,a,a)')" iepoch : kbrr  (mm/s):   latitude     longitude  "  &
          ,":         XYZ (GRACE A) in m      :      XYZvel (GRACE A) in mm/s       " &
          ,":         XYZ (GRACE B) in m       :      XYZvel (GRACE B)   in mm/s"
     DO iepoch = 1,nepochs_t
        IF(nkobs_t > 0) THEN
           WRITE(397,'(i6,3f14.8,13f12.5)')iepoch,pre_omc(13,iepoch)*1.d3  &
                ,satlat, satlon &
                ,pre_omc(1:3,iepoch) &
                ,pre_omc(4:6,iepoch)*1.d3,pre_omc(7:9,iepoch),pre_omc(10:12,iepoch)*1.d3,kb_1pr(iepoch)*1.d6
        ELSE
           WRITE(397,'(i6,3f14.8,13f12.5)')iepoch,0.d0  &
                ,satlat, satlon &
                ,pre_omc(1:3,iepoch) &
                ,pre_omc(4:6,iepoch)*1.d3,pre_omc(7:9,iepoch),pre_omc(10:12,iepoch)*1.d3,0.d0
        ENDIF
     ENDDO
     CLOSE(397)

     ! PT170428: write out a kbra prefit file as well, if using acceleration observations
     IF(nkbra_t > 0)THEN
     CALL status_update('STATUS','GRACESIM','gracefit',' ','Writing out (unfiltered) prefit range acceleration residuals',0)
        OPEN(397,file='prefitRA_unfilt_sim.res',status='unknown')
        WRITE(397,'(a,a)')" iepoch :   kbra  in um/s^2:    latitude     longitude  :"  &
             ,":     XYZ (GRACE A) in m      :      XYZvel (GRACE A) in mm/s      " &
             ,":     XYZ (GRACE B) in m  :      XYZvel (GRACE B)   in mm/s "
        DO iepoch = 1,nepochs_t
           satlon = datan2(rvec(2,1,iepoch),rvec(1,1,iepoch)) *180.d0/pi
           satlat = dasin(rvec(3,1,iepoch)/dsqrt(rvec(1,1,iepoch)**2 + rvec(2,1,iepoch)**2 + rvec(3,1,iepoch)**2))*180.d0/pi
           ! filtered kbra output
           WRITE(397,'(i6,3f14.8,13f12.5,f18.8)')iepoch,pre_omc(iKBRA,iepoch)*1.d6  &
                ,satlat, satlon &
                ,pre_omc(1:3,iepoch) &
                ,pre_omc(4:6,iepoch)*1.d3,pre_omc(7:9,iepoch),pre_omc(10:12,iepoch)*1.d3,kb_1pr(iepoch)*1.d6 &
                ,pre_omc(iKBRA,iepoch)*1.d6
           ! filtered kbra output

        ENDDO
        CLOSE(397)
     ENDIF



  ENDIF   ! end  of filtering KBR observations (or not)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! OUTPUT prefit residual statistics
  ! PT140901: only for the epochs that will be stacked into the normal equations
  sumsq = 0.d0
  DO iepoch = start_neq, end_neq
     DO j=1,ngobs_t*2 + nkobs_t
        sumsq(j) = sumsq(j) + pre_omc(j,iepoch)*pre_omc(j,iepoch)
     ENDDO
  ENDDO
  sumsq = dsqrt(sumsq / DBLE(end_neq-start_neq) )
  OPEN(4132,file='prefit.RMS',status='unknown')
  WRITE(4132,*)'Prefit RMS:',sumsq
  CLOSE(4132)
  ! PT140618: print this to the screen as well as to the file
  WRITE(message,'(a,3f8.2,3f10.4,a)')'Prefit RMS GRACE A obs (pos/vel):',(sumsq(j)*1.d3,j=1,6),' (mm, mm/s)'
  CALL status_update('STATUS','GRACESIM',' ',' ',message,0)
  WRITE(message,'(a,3f8.2,3f10.4,a)')'Prefit RMS GRACE B obs (pos/vel):',(sumsq(j+6)*1.d3,j=1,6),' (mm, mm/s)'
  CALL status_update('STATUS','GRACESIM',' ',' ',message,0)
  IF(nkbrr_t > 0)THEN
  WRITE(message,'(a,f10.4,a)')'Prefit RMS    KBRR obs          :',sumsq(13)*1.d6,' (um/s)'
  CALL status_update('STATUS','GRACESIM',' ',' ',message,0)
  ENDIF
  IF(nkbra_t > 0)THEN
     WRITE(message,'(a,f10.4,a)')'Prefit RMS    KBRA obs          :',sumsq(14)*1.d9,' (nm/s^2)'
     CALL status_update('STATUS','GRACESIM',' ',' ',message,0)
  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!! **** !!!! **** ************************************************************************************************************
!
! ADD RANDOM NOISE TO THE KBRR RPEFIT OBSERVATIONS
!
!!!!! **** !!!! **** ************************************************************************************************************

!!!!! **** !!!! **** ************************************************************************************************************
!
! DEBUG    KBRR PARTIALS
!
    call status_update('STATUS','GRACESIM','gracesim',' ','Writing out dKBRR/dP partial derivatives',0)
    open(398,file='kbrr.partials',status='unknown')
    write(398,'(a)')"dKBRR/dP where P = pos then vels for GRACE A (cols 2-7) and GRACE B (8-13). Column 1 is the epoch"
    do iepoch = 1,last_epoch
      write(398,*)iepoch,((kbrr_part(i,j,iepoch),i=1,6),j=1,2)
    enddo
    close(398)

!!!!! **** !!!! **** ************************************************************************************************************
! PT170713: we branch here to either read a VCV file to allow the computation of postfit residuals or we form the 
!           normal equations and perform the least squares inversion. The flag for one or the other comes from 
!           whether the first character in the so-called "mascon regularisation" file is a "V" (ie VCV file) or a "#" (ie regularisation file)
    if(VCV_flag)then   ! it is a VCV file
   
      call read_soln_v3("GRACESIM",243,nparam_t,nparam_vcv,apriori,soln,prm_input,normeq,.true. &
          ,sTid,sScl, sBias, sEmp, sMasc, n_emp, nScl, nBias)

! make the adjust vector here rather than computing it in a least squares solution
      adjust = soln - apriori
   



    else

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  N O R M A L    E Q U A T I O N    S T A C K I N G    L O O P
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Initialize the normal equations
     AtWb     = 0.d0
     normeq   = 0.d0

     ALLOCATE(soln(nparam_t))
     ALLOCATE(apriori(nparam_t))


     !********************************************************************************************************************************
     ! PT150519: derive the epochs when the solar radiation pressure force either starts or stops acting on the satellites. We use
     !           this information to correlate the y-scale with the x-scale and z-scale for the accelerometers
     ! PT150527: also add the relevant values into the normal equations at this step (from within the subroutine)
     ! PT201109: comment this out. We don't apply shadow constraints anymore.
     ! IF(INT(use_accel_scale_const) == 1) CALL shadow_SRP(thrusters_off,acc,rvec,srf2trf_rotmat,normeq,AtWb,apr_ScaleBias)
     !********************************************************************************************************************************


    CALL status_update('STATUS','GRACESIM','gracesim','',"Masking ",0) 
    ! Quick and Dirty check.    
    do concurrent (i = 1:nepochs_t)
        if (.not. use_obs(i,2 )) then 
               amat((i-1)*nobs_t+ikbrr,:) = 0.0
                pre_omc(ikbrr, i) = 0.0
                !print *, "KBRR", ikbrr, "epoch ", i , "set to 0"
        endif
        if (.not. use_obs(i, 3)) then
                 amat((i-1)*nobs_t+ikbra,:) = 0.0
                 !apr_wght((i-1)*nobs_t+ikbra) = 0.0
                 pre_omc(ikbra, i) = 0.0
                !print *, "KBRA", ikbra, "epoch ", i , "set to 0"
        endif

        if (.not. use_gvec(i, 1) ) then
                amat((i-1)*nobs_t + 1: (i-1)*nobs_t + 6, : )  = 0.0
                pre_omc(1:6, i ) = 0.0
        endif

         if (.not. use_gvec(i, 2) ) then
                 amat((i-1)*nobs_t + 7: (i-1)*nobs_t + 12, : )  = 0.0
                 pre_omc(7:12, i ) = 0.0
         endif
    enddo


    CALL status_update('STATUS','GRACESIM','gracesim',''," ... ADD weight",0)
    ! call output_2d("amat.nc", size(Amat, 1), size(Amat,2),  Amat)
    do i = 1,nobs_t*nepochs_t
       amat(i,:) = amat(i,:) * dsqrt(apr_wght(i))
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL status_update('STATUS','GRACESIM','gracesim','',"Compute KBRA done try Amat",0)
    call DSYRK('U', 'T', nparam_t , nobs_t*nepochs_t, 1.0d0, amat, nobs_t*nepochs_t, 0.0d0, normeq, nparam_t )

     CALL status_update('STATUS','GRACESIM','gracesim',' ','Finished stacking normal equations',0)

     ! build the b vector:
     do concurrent (i = 1:nepochs_t, j = 1: nobs_t)
       bvec(j+(i-1)*nobs_t) = pre_omc(j,i) * dsqrt(apr_wght(j+(i-1)*nobs_t))
     enddo

!     bvec = bvec * dsqrtapr_wght
     CALL DGEMV('T', nobs_t*nepochs_t, nparam_t, 1.0d0, amat, nobs_t*nepochs_t, bvec, 1, 0.0d0, Atwb, 1 )

     !****************************************************************
     ! output the normal equations to a file so that I can play with them without needing to restack all the time
     !    call status_update('STATUS','GRACEFIT','gracefit',' ',"NO LONGER writing out binary normal equations",0)
     CALL output_writeNORM(rmsnam, apr_prm, normeq, AtWb, apr_wght, mcon_tides,prmnam,prmnam_size)

     ! PT191126: stop here if requested by user via the "gracefit_step" command file option.
     if(gracefit_step == 2)then
       call status_update('STATUS',calling_prog,'gracesim',' ',"Program stopped by user after outputting normal equations",0)
       stop
     endif


!**********************************************************************
!
! A D D   C O N S T R A I N T S   T O   N O R M A L   E Q U A T I O N S
!
!**********************************************************************
    call status_update('STATUS','GRACESIM','gracesim',' ',"Adding constraints to normal equations",0)
! PT160215: here we add only up to (but not including) the mascon tidal parameters
    if(imsctide > 0) then
      do i = 1, imsctide - 1 
        normeq(i,i) = normeq(i,i) + apr_const(i)
      enddo
    else
      do i=1, nparam_t
        normeq(i,i) = normeq(i,i) + apr_const(i)
      enddo
    endif   

!open(unit=302,file="apr_const.txt")
!    write(302,*) apr_const
!    close(unit=302)

! PT150820:  constrain all ocean tide amplitudes 
    if(est_msc_tides == 1)then
      do i=imsctide,nparam_t
        if(prmnam(i)(16:18) == "M2s"    )then
          normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(1,1)**2
        elseif(prmnam(i)(16:18) == "M2c")then
          normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(1,2)**2
        elseif(prmnam(i)(16:18) == "O1s")then
          normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(2,1)**2
        elseif(prmnam(i)(16:18) == "O1c")then
          normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(2,2)**2
        elseif(prmnam(i)(16:18) == "S2s")then
          normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(3,1)**2
        elseif(prmnam(i)(16:18) == "S2c")then
          normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(3,2)**2
        elseif(prmnam(i)(16:18) == "K1s")then
          normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(4,1)**2
        elseif(prmnam(i)(16:18) == "K1c")then
          normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(4,2)**2
        elseif(prmnam(i)(16:18) == "K2s")then
          normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(5,1)**2
        elseif(prmnam(i)(16:18) == "K2c")then
          normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(5,2)**2
        endif

! kluge: leave only the amplitudes of the second mascon free to adjust
!      if(i < nparam_t-9)normeq(i,i) = normeq(i,i)+1.0/1e-12
      enddo
    endif
!****************************************************************
! Add mascon length-dependent constraint
    if(nmascons_t > 0 .and. use_msc_len_const == 1 )then
      call status_update('STATUS','GRACESIM','gracesim',msc_const_file &
                   ,"    Adding length-dependent constraints to mascon parameters",0)

! PT170307: read the second line in the regularisation file, which contains the length scales and constraint levels
      read(243,'(a)')message
! PT150721: check that it really was a header line
      if(message(1:7) == "# scale" .or. message(1:5) == "# lat")then
        call status_update('STATUS','GRACESIM','gracesim',' ',message,0)
      else
        call status_update('WARNING','GRACESIM','gracesim',' ',message,0)

        write(message,'(a)')'No header info found in mascon spatial constraint file '
        call status_update('WARNING','GRACESIM','gracesim',msc_const_file(1:trimlen(msc_const_file)),message,0)
        backspace(243)
      endif
! read the constraint matrix
      msc_const = 0.d0
      do i=1,nmascons_t
        ! PT210131: check whether the regularisation file is a vector of diagonal elements or a full matrix
        if(header_line(61:69) == "diag_only")then
          read(243,*,iostat=ioerr)msc_const(i,i)
        else
          read(243,*,iostat=ioerr)(msc_const(i,j),j=1,nmascons_t)
        endif
        if(ioerr /= 0)then
          write(message,'(a,i6,a,i6)')"Error reading line",i," of the mascon regularisation file which needs to have",nmascons_t
          call status_update('FATAL','GRACESIM','gracesim',msc_const_file,message,0)
        endif          
      enddo
      do i=imascons,imascons+nmascons_t-1
        do j=imascons,imascons+nmascons_t-1
          normeq(i,j) = normeq(i,j) + msc_const(i-imascons+1,j-imascons+1) ! *0.7 ! / 0.04 ! PT150618: scale it by 0.2 m sigma
        enddo
      enddo
    else
      call status_update('STATUS','GRACESIM','gracesim',' ',"NOT adding length-dependent constraints to mascon parameters",0)
    endif


!****************************************************************


     !**********************************************************************
     !**********************************************************************
     ! Add mascon mass conserving constraint if necessary           PT150721
     !**********************************************************************
     IF(use_msc_conserve_mass == 1)THEN
        CALL status_update('STATUS','GRACESIM','gracefit',' ',"    Adding mass conservation constraints to mascon parameters",0)

        ! the constraint equation is simply: m1+m2+m3+ ..... +m_n = 0
        ! therefore, the partials d/dm? = 1 for each mascon and the terms in the normal equations are just the variance of the observation
        DO i=imascons,imascons+nmascons_t-1
           DO j=imascons,imascons+nmascons_t-1
              normeq(i,j) = normeq(i,j) + 1.0/0.01**2  ! hardwire 1 cm as the sigma for this conditional equation. @#&& fix this later !
           ENDDO
        ENDDO
     ELSE
        CALL status_update('STATUS','GRACESIM','gracesim',' ',"NOT adding mass conservation constraints to mascon parameters",0)
     ENDIF
     !**********************************************************************

     Call status_update('STATUS', 'GRACESIM', 'gracesim', '', "Solving equations", 0 )
     !**********************************************************************
     !   Solve the normal equations for the estimated ICs and parameters
     IF(INT(use_kbrr_decomp(1)) == 0  .OR. use_kbrr_decomp(1) == 2.OR. use_kbrr_decomp(1) == 3 .or. use_kbrr_decomp(1) == 5)THEN  ! solve for all the observations and parameters
        CALL Chol_normSolve(nparam_t,normeq,AtWb,adjust)
     ELSE IF(INT(use_kbrr_decomp(1)) == 1)THEN  ! solve for only the mascons
        CALL Chol_normSolve(nmascons_t,normeq(imascons:imascons+nmascons_t-1,imascons:imascons+nmascons_t-1) &
             ,AtWb(imascons:imascons+nmascons_t-1),adjust(imascons:imascons+nmascons_t-1))
     ENDIF

     Call status_update('STATUS', 'GRACESIM', 'gracesim', '', "Solving equations done", 0 )
     !******************************************************************
  ENDIF   ! end of if statement as to whether to perform the LS inversion or read a solution from a VCV file
  ! Compute post fit. 
  ! Descale amat.
  do i = 1,nobs_t*nepochs_t
    amat(i,:) = amat(i,:)/dsqrt(apr_wght(i))
  enddo
  do i = 1, nepochs_t
    do j = 1, nobs_t
      bvec(j+(i-1)*nobs_t) = pre_omc(j,i)
    enddo
  enddo
  ! PT201113: this is the maths done in the call to DGEMV
  ! bvec(post_omc) = - amat * adjust + bvec(pre_omc) 
  CALL DGEMV('N', nobs_t*nepochs_t, nparam_t, -1.0d0, amat, nobs_t*nepochs_t, adjust, 1, 1.0d0, bvec, 1 )
  !write bvec into post_omc
  ! build the b vector:
  do i = 1, nepochs_t
    do j = 1, nobs_t
      post_omc(j,i) = bvec(j+(i-1)*nobs_t) 
    enddo
  enddo
  ! PT201113: compute  unfiltered kbrr postfit residuals
  IF(nkbrr_t > 0)post_omc_kb_unfilt(1,:) = post_omc(iKBRR,:) - pre_omc(iKBRR,:) + pre_omc_kb_unfilt(1,:)
  IF(nkbra_t > 0)post_omc_kb_unfilt(2,:) = post_omc(iKBRA,:) - pre_omc(iKBRA,:) + pre_omc_kb_unfilt(2,:)

  !********************** WRITE OUT SOLUTIONS *********************
  call output_writeVCV(gt_fnam,apr_prm,kbrr_misfit_num,adjust,normeq,prmnam,prmnam_size)
  CALL output_writeFIT(gt_fnam,apr_prm,apr_wght,adjust,kbrr_misfit_num,Amat, pre_omc,normeq,post_omc,prmnam,prmnam_size &
                       ,msc_const_file)  ! postfit residuals are computed in this subroutine
!****************************************************************

!   Plot the residuals
    if( iop > 0 )then
     IF(nkbrr_t /= 0) CALL plot_writeKB(iop,pre_omc,post_omc,post_omc_kb_unfilt(iKBRR,:),rvec(1:3,:,:),rpy,kbrr_prefit_tol,nepochs_t)
     IF(nkbra_t /= 0) CALL plot_writeKBA(iop,pre_omc,post_omc,post_omc_kb_unfilt(iKBRA,:),rvec(1:3,:,:),rpy,kbrr_prefit_tol,nepochs_t)
    endif
!****************************************************************

! ******************!  Output the mean beta angle *************************
   write(message,'(a,f8.2,a,f8.2)')"Mean beta angle = ",sum(beta_angle(:))/dble(num_gps_used(1))*180.d0/pi," degrees."
   call status_update('STATUS','GRACESIM','gracesim',' ',message,0)
!****************************************************************


!*********************** DEALLOCATE ARRAYS **********************

    call status_update('STATUS','GRACESIM','gracesim',' ','Deallocating arrays',0)
    deallocate(apr_prm)
    deallocate(apr_const)
    deallocate(apr_wght)
    deallocate(Amat)
    deallocate(normeq)
    deallocate(AtWb)
    deallocate(adjust)
    deallocate(pre_omc)
    deallocate(post_omc)
    deallocate(rvec)
    deallocate(sciframe_quat)
    deallocate(srf2trf_rotmat)
    deallocate(srf2trf_deriv_rotmat)
    deallocate(rpy)
    deallocate(AOC)
    deallocate(GPS_antoff_trf)
    deallocate(GPSant_adj)
    deallocate(acc)
    deallocate(gvec)
!    deallocate(uvec)
    deallocate(kbrange)
    deallocate(kblt)
    deallocate(kbant)
!    deallocate(kbion_cor)
!    deallocate(kbsnr)
    deallocate(thr)
    deallocate(thrusters_off)
    deallocate(shadow_counter)
    deallocate(apply_shadow_cond)
    deallocate(shadow)
    deallocate(filt_kbrr_omc)
!****************************************************************

!   That's all folks
    call status_update('STATUS','GRACESIM','gracesim',' ', 'Normal end of GRACESIM',0)

  end program

!********************************************************************************************************************************
