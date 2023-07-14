!********************************************************************************************************************************
!  File: gracefit.f90
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
!           subroutines in jplephem.f (which, incidentally, both names and opens the file !!)
!
! PT150721: add a conservation of mascon mass constraint to the inversion through a conditional equation
! PT180625: modified to read L1B ascii files for KBR, ACC, THR, LRI for GRACE/GRACE FO 
! PT109095: modified to re-invoke the capability to estimate only postfit residuals, given parameter estimates in a vcv file
! PT141020: added new logical array to indicate whether to use range/range rate/range accel observations
! PT220922: write out normal equations in netcdf format
!
! Output files:  .rms file
!                .svs file
!                .vcv file
!                .fit file
!                .resid file
!                .corr file
!                plot .kb file
!                plot .A/B file(s)
!                .nc file
!
!********************************************************************************************************************************

PROGRAM GRACEFIT

  ! PT170609: include the mod file to define the record length of the GTORB files
  USE gtorb_mod
  USE gracefit_mod
  use accred_mod      ! declares the variables for the accelerometer observations
  use soln_mod        ! PT190905: adds the declaration of the arrays for the solution (apriori, adjust, soln)
!  use norm_netcdf_mod ! PT221122: adds the variables for writing out a netcdf file of normal equations

  IMPLICIT NONE


  !********************  Variable declarations ********************

  ! Paramters
  ! PT180626: reduced by one so that we don't have to put a mascon reg file if not estimating mascons.
  INTEGER*4, PARAMETER :: NUM_EXPECTED_ARGS = 8  ! Number of expected arguments. 

  ! Status variables
  CHARACTER(120) :: version         ! Version of program being run
  CHARACTER(256) :: message         ! Message to be written out to standard output

  ! Command line related variables
  CHARACTER(80) :: cmdnam,rmsnam                         ! Files to be read from the command-line
!  CHARACTER(80) :: gt_fnam(2)                            ! GTORB file names to be used
  CHARACTER(80) :: arg                                   ! Dummy variable to be read in iop from the command-line
  INTEGER*4     :: iop                                   ! Plot resid indicator
  CHARACTER     :: date_str*10,c_yr*4,c_month*2,c_day*2  ! Variables for reading the date from the command-line
  real :: start, finish
  ! Variables storing data from command file and GTORB file(s)
  DOUBLE PRECISION, ALLOCATABLE :: apr_prm(:)   ! A priori values for parameters
  DOUBLE PRECISION, ALLOCATABLE :: apr_const(:) ! A priori constant values
  DOUBLE PRECISION, ALLOCATABLE :: apr_tide_ampl_const(:,:) ! A priori tidal amplitude constant values
  DOUBLE PRECISION, ALLOCATABLE :: apr_wght(:)  ! A priori data weights
  DOUBLE PRECISION :: kbrr_prefit_tol           ! Orbit misfit and kbrr prefit tolerances

  ! Data read in or calculated from GTORB file(s)
  integer*4,parameter           :: maxsats = 2
  INTEGER*4                     :: rec_len(2)                    ! record length of each GTORB file
  INTEGER*4                     :: irec                          ! record number for reading binary GTORB file
  DOUBLE PRECISION, ALLOCATABLE :: rvec(:,:,:)                   ! Vector of positions, velocities and partials
  DOUBLE PRECISION, ALLOCATABLE :: sciframe_quat(:,:,:)          ! Quaternions for each satellite for every epoch
  DOUBLE PRECISION, ALLOCATABLE :: srf2trf_rotmat(:,:,:,:)       ! SRF to TRF rotation matrix for GRACE A and B
  DOUBLE PRECISION, ALLOCATABLE :: srf2trf_deriv_rotmat(:,:,:,:) ! Differentiated SRF to TRF rotation matrix for GRACE A and B
  DOUBLE PRECISION :: satics_t(maxrvecprm+1,2)                   ! A priori values for IC's
  DOUBLE PRECISION :: apr_ScaleBias(max_SB,maxsats)               ! A priori values for scale and bias
  DOUBLE PRECISION :: GPSantoff(7,maxsats)                        ! Quaternion and offsets for antenna

  DOUBLE PRECISION:: tempsum
  integer:: k
  ! PT170206: add mascon filenames and codes
  CHARACTER*40  :: combined_mascon_file,ocean_mascon_file        ! name of temporal and ocean mascon files used in GRACEORB
  CHARACTER*8   :: msc_hdr_code,msc_ocn_hdr_code,msc_reg_code    ! mascon file codes of temporal and ocean mascon files used in GRACEORB
  CHARACTER*256 :: header_line                                   ! RM210415: contains header of reg file
  INTEGER*4     :: total_ocean_prim
  ! Variables used for rpy calculations and antenna offsets
  DOUBLE PRECISION, ALLOCATABLE :: rpy(:,:,:)          ! Roll, Pitch, Yaw of satellites for each epoch
  DOUBLE PRECISION, ALLOCATABLE :: AOC(:,:,:)          ! Antenna offset corrections
  DOUBLE PRECISION, ALLOCATABLE :: GPS_antoff_trf(:,:) ! Antenna offset in TRF
  DOUBLE PRECISION, ALLOCATABLE :: GPSant_adj(:,:)     ! Antenna offsets to be read in from command line

  ! Data read in from ACC1B file(s)
  DOUBLE PRECISION, ALLOCATABLE :: acc(:,:,:)    ! Acceleration data

  ! Data read in or calculated from THR1B file(s)
  DOUBLE PRECISION, ALLOCATABLE :: thr(:,:,:)          ! Thruster data  !@# Change index order (does it even need to be here?)
  LOGICAL, ALLOCATABLE          :: thrusters_off(:,:)  ! Data on whether thrusts are affecting accelerations measurements

  ! Data read in from GNV1B files
  DOUBLE PRECISION, ALLOCATABLE :: gvec(:,:,:)       ! Positions and velocities for each satellite at every epoch
  INTEGER*4       , ALLOCATABLE :: num_gps_used(:)   ! number of GPS observations that will be used (after decimating the 5sec data)
  logical         , ALLOCATABLE :: use_gvec(:,:)          ! logical to indicate whether there are GNV1B obs to match the GTORB obs
  integer*4                     :: igvec             ! counter through the GNV1B records
  ! PT190528: pass this variable back from input_readGNV1B
  integer*4                     :: n_gnv             ! number of GNV1B epochs read from the input GNV1B file

  !   Variables that might be used in the future but not needed in current version program
  !    double precision, allocatable :: uvec(:,:,:)  ! Uncertainties on positions and velocities for each satellite at every epoch
  CHARACTER(1)               :: lead
  ! Data read in from KBR1B file
  DOUBLE PRECISION, ALLOCATABLE :: kbrange(:,:)    ! Range value (biased), rate, and acceleration
  DOUBLE PRECISION, ALLOCATABLE :: kblt(:,:)       ! Light time corrections for range, rate and acceleration
  DOUBLE PRECISION, ALLOCATABLE :: kbant(:,:)      ! Antenna phase center corrections for range, rate and acceleration
  !   Variables that might be used in the future but not needed in current version program
  !    double precision, allocatable :: kbion_cor(:)  ! Ionospheric range correction for Ka frequencies
  !    double precision, allocatable :: kbsnr(:,:,:)  ! SNR K and Ka bands for both satellites
  INTEGER*4                     :: kbrr_misfit_num ! Number of kbrr misfits
  ! PT170703: 
  LOGICAL                       :: override_kbr_filt     ! flag to override the filtering of kbr data
  logical                       :: kb_differentiate_flag ! flag to indicate whether to differentiate the L1B data or not

  ! RM190724:
  LOGICAL                       :: override_kbr_cos_fft  ! flag to override the filtering of kbr data

  ! Counter variables
!  INTEGER*4 :: iepoch    ! Counter that runs through epochs (numbered)
  INTEGER*4 :: isatt     ! Counter that runs through satellites
  INTEGER*4 :: i,j       ! Counter variable

  ! Variables used for calculations
  DOUBLE PRECISION, ALLOCATABLE :: pre_omc(:,:)   ! Prefit Observed Minus Computed values to be used in LS algorithm
  DOUBLE PRECISION, ALLOCATABLE :: pre_omc_kb_unfilt(:,:)   ! KBRR unfiltered prefit  OMC, backup copy
  DOUBLE PRECISION, ALLOCATABLE :: post_omc_kb_unfilt(:,:)  ! KBRR unfiltered postfit OMC, backup copy
  DOUBLE PRECISION, ALLOCATABLE , target:: part(:,:,:)    ! Partials of observables with respect to the parameters
  double precision, dimension(:,:), pointer :: ptr_part
  DOUBLE PRECISION, ALLOCATABLE :: normeq(:,:)    ! The left side of the LS algorithm. P_t*W*P
! PT221122: now defined in norm_netcdf_mod
!  DOUBLE PRECISION, ALLOCATABLE :: AtWb(:)        ! The right side of the LS algorithm. P_t*W*OMC
  DOUBLE PRECISION, ALLOCATABLE :: adjust(:)      ! Solution to normal equations
  DOUBLE PRECISION, ALLOCATABLE :: post_omc(:,:)  ! Postfit Observed Minus Computed values
  DOUBLE PRECISION, ALLOCATABLE :: beta_angle(:)  ! beta angle at each GPS epoch
  DOUBLE PRECISION, ALLOCATABLE :: bvec(:)
  ! Variables used for application of shadow condition
  INTEGER*4, ALLOCATABLE :: shadow_counter(:)       ! Variable that counts the number of epochs the satellites are in shadow for
  LOGICAL, ALLOCATABLE   :: apply_shadow_cond(:,:)  ! Data on whether the shadow condition can be applied
  DOUBLE PRECISION, ALLOCATABLE :: shadow(:,:)      ! Variable used to indicate whether satellite is in shadow

  ! PT130813: variables used for inversion using just mascons and emf-decomposed kbrr omc
  DOUBLE PRECISION, ALLOCATABLE :: filt_kbrr_omc(:)      ! sum of the decomposed elements of the kbrr omc
  DOUBLE PRECISION, ALLOCATABLE :: filt_kbra_omc(:)      ! sum of the decomposed elements of the kbra omc
  DOUBLE PRECISION, ALLOCATABLE :: kb_1pr(:)          ! model to fit a 1/rev with linearly varying amplitude to the prefit kbrr residuals
  double precision              :: kb_sponge_model(7)    ! all the parameters for a linear plus sinusoid model

  ! RM190723: variables used for filtering kb with fft cosine
  INTEGER*4 :: nfft

  ! PT130826: variables to apply length-dependent mascon constraints
  DOUBLE PRECISION, ALLOCATABLE :: msc_const(:,:)
  CHARACTER                     :: msc_const_file*100
  INTEGER*4 :: trimlen        ! Function that returns the length of a string
  INTEGER*4 :: ioerr

  ! PT140527: mascon ocean/land identifier and tidal mascon information
  LOGICAL,          ALLOCATABLE :: mcon_ocean(:)
  INTEGER*4                     :: mcon_tides(4000,2)    ! we have to declare the dimensions of this array before calling header_readGTORB ....
  REAL(kind=8),     ALLOCATABLE :: msc_tide_amp(:,:,:)

  ! PT131202: variables for outputting prefit residuals before stacking normal equations
  DOUBLE PRECISION :: satlon, satlat
  REAL(kind=8),     ALLOCATABLE :: sumsq(:)

  ! PT130902: variables for openMP stacking of normal equations
  DOUBLE PRECISION, ALLOCATABLE :: normeq_tmp(:,:,:)    ! The left side of the LS algorithm. P_t*W*P     tmp for openMP implementation
  DOUBLE PRECISION, ALLOCATABLE :: AtWb_tmp(:,:)        ! The right side of the LS algorithm. P_t*W*OMC  tmp for openMP implementation

  DOUBLE PRECISION, ALLOCATABLE :: normeq_tmp2(:,:)    ! The left side of the LS algorithm. P_t*W*P     tmp for openMP implementation
  DOUBLE PRECISION, ALLOCATABLE :: AtWb_tmp2(:)    

  INTEGER*4 count_epoch_loop
  ! PT140617: last epoch that we want to process (hard-wired for now, but should set through the command file)
  INTEGER*4 :: last_epoch

  ! PT140523: range rate derived from GNV1B pos and vel
  REAL(kind=8), ALLOCATABLE :: kb_theor_gvec(:,:)

  ! PT151016: angles between sunsat vector and spacecraft XYZ axes
  REAL(kind=8) :: SRPdotprod

  ! PT/HMcQ160704 : variable to count the number of obs in shadow that are unaffected by thrusts
  INTEGER*4    :: n_shadow_obs=0

  ! PT170713: variables for reading the VCV file to enable the computation of postfit residuals
  LOGICAL      :: VCV_flag
  INTEGER*4    :: nparam_vcv, nfiles_vcv, nfile                       ! number of parameters found in the VCV file
  INTEGER*4    :: sTid,sScl, sBias, sEmp, sMasc, n_emp, nScl, nBias   ! we don't need these for the postfit residual computations
  CHARACTER*30,ALLOCATABLE :: prm_input(:)
  INTEGER*4,ALLOCATABLE    :: dates_vcv(:,:)

  ! PT140429: debug variables
  REAL(kind=8) :: kbtheoretical(3), IC_adj(12,2), prefit_adj(13),cmd_args(12),sumsq_for_rms(13),rms(13)
  REAL(kind=8) :: delta_kbrr
  REAL(kind=8) :: dot,amag3
  REAL(kind=8) :: tmp_LOS_A(3),tmp_LOS_B(3),tmp_LOS_mag,tmp_vel_A(3),tmp_vel_B(3),tmp_vel_mag
  REAL(kind=8), ALLOCATABLE :: kbrr_part(:,:,:),kbra_part(:,:,:)
  real (kind=8), allocatable :: Amat(:,:)
  ! PT170725: variables used in removing edge effects introduced when filtering prefit residuals
  INTEGER*4    :: edge1(10000),edge2(10000)    ! list of epochs of start and stop of gaps in kbr obs
  INTEGER*4    :: n_kbr_gaps               ! number of data gaps in the kbr observations
  LOGICAL      :: in_kbr_gap               ! flag for whether in a kbr data gap or not

  ! HM170912: variables for output of mascon partials
  INTEGER*4    :: im,imp1,imp2
  INTEGER*4, ALLOCATABLE  :: mascs(:)      ! mascons for which partials are to be output

  LOGICAL :: is_h5
  CHARACTER(len=80) :: test
  INTEGER :: length


  character (len=512) :: commandline
  integer :: lencmd, status_cmd

! PT190514: variables used to investigate the weird Arctic signals in Nov 2016 GRACE data
  real(kind=8) :: kbr_weights(2)

! RM191112: counter for counting every call to calc_beta_engle
  integer*4    :: nbeta
  
! PT221122: variables needed for writing out netcdf files
  integer*4    :: iwhere
  integer*4    :: sats(2)
  integer*4    :: epoch_netcdf(6)
  character*150:: netcdf_out
  logical      :: IC_solve,msc_solve,tmsc_solve
  integer*4    :: IC_offset,msc_offset
  real(kind=8),allocatable :: AtWb_flip(:),apriori_flip(:)
  integer*4   ,allocatable :: param_type_flip(:)

! PT200706: variables to store the inertial-efixed rotation and rotation rate information
  integer*4                 :: nepochs                  ! number of epochs of celestial-to-terrestrial rotations computed
  real(kind=8),allocatable  :: rot_mats(:,:,:)          ! array of inertial->efixed rotation matrices (one per epoch)                  
  real(kind=8),allocatable  :: rotdots(:,:,:)           ! array of time derivatives of rotation matrices
  real(kind=8),allocatable  :: rotaccs(:,:,:)           ! array of time derivatives of rotation matrices
  real(kind=8),allocatable  :: rpy_tmp(:,:)             ! temp array needed for calling generate_inert_efixed (but not used)
  real(kind=8),allocatable  :: rpydots(:,:)             ! roll/pitch/yaw angle rates for the rotation matrices
  real(kind=8),allocatable  :: rot_dates(:)             ! epochs for which rotation angles and rates have been computed
  real(kind=8)              :: mjd_start,mjdend         ! floating point start/end epoch of the orbit integration
  real(kind=8)  :: rot_i2e(3,3)           ! inert->efixedf rotation matrix for epoch "mjd"
  real(kind=8)  :: rotdot_i2e(3,3)        ! inert->efixed rotation rate matrix for epoch "mjd"
  real(kind=8)  :: rotacc_i2e(3,3)        ! inert->efixed rotation accel matrix for epoch "mjd"
  real(kind=8)  :: rot_e2i(3,3)           ! efixed->inert rotation matrix for epoch "mjd"
  real(kind=8)  :: rotdot_e2i(3,3)        ! efixed->inert rotation rate matrix for epoch "mjd"
  real(kind=8)  :: rotacc_e2i(3,3)        ! efixed->inert rotation accel matrix for epoch "mjd"
  integer*4     :: date_tmp(5)            ! used in the call ymdhms_to_jd           
  real(kind=8)  :: mjd,jd                 ! modified julian date and julian date for start of solution
  real(kind=8)  :: sunpos(6)              ! inertial sun coordinates (pos and vel)
  real(kind=8)  :: sunpos_efixed(6)       ! earth-fixed sun coordinates
  real(kind=8)  :: sec_tmp                ! temporary storage of seconds of day for the starting epoch
  integer*4     :: interp_flag            ! 1: only the rot_i2e and rot_e2i matrices; 2: all rotation rate and accel matrices

! RM
  integer*4 :: maxparam
  real(kind=8),allocatable :: VCV_obs_local(:,:)

! PT201109: 5-element integer to store yr/mo/day/hr/min
  integer*4    :: date(5)

  integer :: ompthread, ompmax

  INTEGER, EXTERNAL :: omp_get_thread_num, omp_get_max_threads
  call cpu_time(start)
  !****************************************************************
  calling_prog = "GRACEFIT  "

  ! If there are no arguments, give instructions as to program usage
  IF( iargc() == 0 ) THEN
     WRITE(*,10)
10   FORMAT(/,'###################################################################################################' &
          ,/,'                                  Program GRACEFIT(new)                         ' &
          ,/,'Purpose: Fit integrated GRACE GTORB orbits to imported GNV1B GPS tabular and KBR1B range rate files' &
          ,/,'Usage: gracefit cmd-file rms-file plt-option GTORB_A GTORB_B YR MO DD' &
          ,/ &
          ,/,'Example: gracefit_new gracefit_new.cmd thomas_new.rms 4 GTORB_2010-09-14_00-00-00_A_02.asc_iter3' &
          ,' GTORB_2010-09-14_00-00-00_B_02.asc_iter3 2010 09 14 0 0 0 0 0 0' &
          ,/,'###################################################################################################')

     WRITE(*,11)
11   FORMAT(/,' a template gracefit.cmd file format is given in ~/gg/grace/tables/gracefit.cmd.template:',/, &
          'Simply add your preferred GPS and KBRR observation sigmas and you are ready to go.',/, &
          'As Tom Herring says "with perfect data, this should be easy !"',/)
     STOP
  ENDIF
  !****************************************************************

  ! Get version number and write out status line
  CALL gversn(version)
  WRITE(*,'(a)')' '
  WRITE(message,'(a,a120)') 'Started GRACEFIT ',version
  !CALL status_update('STATUS',calling_prog,'gracefit/gracefit',' ',message,0)

  ! Assert that there are the correct number of command-line arguments
  IF ( iargc() < NUM_EXPECTED_ARGS )THEN
     CALL status_update('FATAL',calling_prog,'gracefit',' ','Too few arguments in command line',0)
  ENDIF

  !************** PARSE COMMAND LINE AND COMMAND FILE *************
  call get_command (commandline, lencmd, status_cmd)
  CALL status_update('STATUS', calling_prog, 'gracefit/cmdline', ' ', trim(commandline),0)


  ! Now read in the command-line and command file
  CALL getarg(1,cmdnam)         ! Command file
  CALL getarg(2,rmsnam)         ! rms file name
  CALL getarg(3,arg)            ! Ploting indicator. Placed into iop
  READ(arg,'(i1)') iop
  CALL getarg(4,gt_fnam(1))     ! Integrated GRACE GTORB orbits for satellites A and B
  CALL getarg(5,gt_fnam(2))
  CALL getarg(6,c_yr)           ! Year
  CALL getarg(7,c_month)        ! Month
  CALL getarg(8,c_day)          ! Day
  ! Get instructions
  CALL input_openCommand(cmdnam)
! PT221122: allocate the prmnam array. Was done in soln_mod but now needs to be done in calling program
  allocate(prmnam(prmnam_size))
  CALL command_readSetup()
  write(message,'(a,i2)')"   GRACE mission: ",mission
  call status_update('STATUS',calling_prog,'gracefit',' ',message,0)

! PT190905: we need the mascon regularisation file irrespective of whether we want to estimate mascons or not because
!           we use it to derive the postfit residuals only
  CALL getarg(9,msc_const_file) ! mascon constraint file name
  VCV_flag = .false.
  if(msc_const_file(1:1) /= "")then  ! a file was supplied on the command line
    OPEN(243,file=msc_const_file(1:trimlen(msc_const_file)),status='old',iostat=ioerr)
    IF(ioerr /= 0)THEN
      CALL status_update('FATAL',calling_prog,'gracefit',msc_const_file(1:trimlen(msc_const_file)) &
                       ,'Error opening msc regularization',0)
    ENDIF
    ! PT170206: verify that the mascon file code is compatible with that of the GTORB files
    !READ(243,'(a8)')msc_reg_code
    ! RM210415: store header line and extract code so that I can test reg file type later
    READ(243,'(a)')header_line
    msc_reg_code = header_line(1:8)

    ! PT170713: check whether it is a VCV file to be used for postfit residual computations
    IF(msc_reg_code(1:1) == "V")THEN
      VCV_flag = .TRUE.
      CALL status_update('STATUS',calling_prog,'gracefit',msc_const_file &
                       ,'Will compute postfit residuals using solution in VCV file',0)    
    endif   
    !CLOSE(243)  ! JP191122: close file, opened again with output_openAllFiles at L431! 
  endif

  test = gt_fnam(1)
  length = LEN(TRIM(gt_fnam(1)))
  IF( test(length-2:length) == '.h5') THEN
     call status_update('STATUS',calling_prog,'gracefit',' ','   Input GTORB files are hdf5 format',0)
     is_h5 = .TRUE.
  ELSE
     call status_update('FATAL',calling_prog,'gracefit',trim(test),'Input GTORB files do not end in ".h5" ',0)
     is_h5 = .FALSE.
  ENDIF
  IF (is_h5 ) THEN
     CALL input_openh5GTORB(gt_fnam,0)
     !NOTE: input_readGTORB MUST follow directly after command_readSetup (i.e. before anything else is set)
     CALL header_readh5GTORB(irec,satics_t,apr_ScaleBias,GPSantoff,mcon_tides &
          ,combined_mascon_file,msc_hdr_code,ocean_mascon_file,msc_ocn_hdr_code,total_ocean_prim,0,prmnam,prmnam_size)
  ELSE
     call status_update('FATAL','program_name','gracefit',' ',"Binary GTORB files no longer supported",0)
  ENDIF
  ! PT140901: set end_neq to be nepochs_t if it wasn't set in the command file
  IF(end_neq <= 0)end_neq = nepochs_t
  ! PT160223: reduce end_neq to the end of the orbit if it is too large in the command file
  IF(end_neq > nepochs_t) end_neq =  nepochs_t
  CALL command_printSetup() ! Must come after reading GTORB file(s)

! PT180625: move to here the reading of the mascon regularisation file. Only read the command line argument IF we want to estimate mascons
  if(use_msc_len_const == 1)then
    CALL getarg(9,msc_const_file) ! mascon constraint file name
  endif
  !************************ ALLOCATE ARRAYS ***********************

  CALL status_update('STATUS',calling_prog,'gracefit',' ','            Allocating arrays',0)
  ! Variables storing data from command file and GTORB file
  ALLOCATE(apr_prm(nparam_t))
  ALLOCATE(apr_const(nparam_t))
! PT201014: apr_wght now needs dimension of nobs_t*nepochs*t (not just nobs_t)
  ALLOCATE(apr_wght(nobs_t*nepochs_t))
  ALLOCATE(apr_tide_ampl_const(max_mcon_tides,2))
  ! Variables used for calculations
! PT201013: replace "part" with the new Amat matrix
!  ALLOCATE(part(nobs_t,nparam_t,nepochs_t))
  allocate(Amat(nobs_t*nepochs_t,nparam_t))
  ALLOCATE(normeq(nparam_t,nparam_t))
  ALLOCATE(AtWb(nparam_t))
  ALLOCATE(pre_omc(nobs_t,nepochs_t))
  ALLOCATE(pre_omc_kb_unfilt(2,nepochs_t))   ! backup copy of unfiltered prefit KBRR/KBRA residuals. Used to create unfiltered postfit kbrr residuals for .kb file
  ALLOCATE(post_omc_kb_unfilt(2,nepochs_t))  ! backup copy of unfiltered prefit KBRR/KBRA residuals. Used to create unfiltered postfit kbrr residuals for .kb file
  ALLOCATE(post_omc(nobs_t,nepochs_t))
  allocate(bvec(nobs_t*nepochs_t))
! PT201014: array to indicate whether to use observations or not
! PT220518: increased to 8 obs to include pos and vel
  allocate(use_obs(nepochs_t,8))   ! 1: kbr; 2: kbrr; 3: kbra; 4: LR; 5: LRR; 6: LRA  (first 3: kbr, last 3: LRI)
  allocate(use_gvec(nepochs_t,2))  ! 1: kbr; 2: kbrr; 3: kbra; 4: LR; 5: LRR; 6: LRA  (first 3: kbr, last 3: LRI)
  use_obs(:,:) = .true.  ! PT201014: initialise the use of all inter-satellite observations (kbr and LRI) to be false

  ! Data read in or calculated from GTORB file(s)
  ALLOCATE(rvec(6,nsat_t,nepochs_t))  !SAChange
  ALLOCATE(sciframe_quat(4,nsat_t,nepochs_t))
  ALLOCATE(srf2trf_rotmat(3,3,nepochs_t,nsat_t))
  ALLOCATE(srf2trf_deriv_rotmat(3,3,nepochs_t,nsat_t))
  ! Variables used for rpy calculations and antenna offsets
  ALLOCATE(rpy(3,nsat_t,nepochs_t))
  ALLOCATE(AOC(9,3,nepochs_t))
  ALLOCATE(GPS_antoff_trf(6,nsat_t))
  ALLOCATE(GPSant_adj(3,nsat_t))
  ! Data read in from ACC1B file(s)
  ALLOCATE(acc(nepochs_t,4,nsat_t))
  ! Data read in from GNV1B files
  ALLOCATE(gvec(maxgprm,nsat_t,nepochs_t))
  ALLOCATE(num_gps_used(nsat_t))
  !    allocate(uvec(maxgprm,nsat_t,nepochs_t))
  ! Data read in from KBR1B file
  ALLOCATE(kbrange(3,nepochs_t))
  ALLOCATE(kblt(3,nepochs_t))
  ALLOCATE(kbant(3,nepochs_t))
  !    allocate(kbion_cor(nepochs_t))
  !    allocate(kbsnr(2,nsat_t,nepochs_t))
  ! Data read in or calculated from THR1B file(s)
  ALLOCATE(thr(nepochs_t,16,nsat_t))
  ALLOCATE(thrusters_off(nepochs_t,nsat_t))
  ! Variables used for application of shadow condition
  ALLOCATE(shadow_counter(nsat_t))
  ALLOCATE(apply_shadow_cond(nepochs_t,nsat_t))
  ALLOCATE(shadow(nsat_t,nepochs_t))
  ! filtered omc for the kbrr and kbra observations
  ALLOCATE(filt_kbrr_omc(nepochs_t))
  ALLOCATE(filt_kbra_omc(nepochs_t))
  ALLOCATE(kb_1pr(nepochs_t))
  ! beta angle
  ALLOCATE(beta_angle(nepochs_t))
  ! mascon length-dependent constraints
  ALLOCATE(msc_const(nmascons_t,nmascons_t))
  ! mascon ocean/land identifier
  ALLOCATE(mcon_ocean(nmascons_t))
  ! range/range rate/range accel derived from GNV1B pos and vel
  ALLOCATE(kb_theor_gvec(3,nepochs_t))
  ! mascon tidal amplitudes
  ALLOCATE(msc_tide_amp(max_mcon_tides,2,nmascons_t))
  ALLOCATE(kbrr_part(6,2,nepochs_t))
  ALLOCATE(kbra_part(6,2,nepochs_t))
  ! varaibles for prefit residual statistics
  ALLOCATE(sumsq(ngobs_t*2 + nkobs_t))
  ! variables for reading the VCV file
  ALLOCATE(adjust(nparam_t))
  ALLOCATE(prm_input(nparam_t))
  ! variable for macon partial output
  ALLOCATE(mascs(nobs_t))


  !****************************************************************

  ! Read in the rest of the command file and command-line
  CALL command_readAprConst(apr_const,apr_tide_ampl_const)
  CALL command_readDataWeights(apr_wght)
  kbr_weights(1:2) = apr_wght(ikbrr:ikbra)
  CALL command_readTol(kbrr_prefit_tol)
  CALL command_close()
  ! PT140918: define GPSant_adj to be zero
  GPSant_adj = 0.d0

  ! Write out date_str to be used for finding data files
  WRITE(date_str,'(a4,a1,a2,a1,a2)')c_yr,"-",c_month,"-",c_day

  !****************************************************************

  ! Open input and output files
  ! PT180625: this includes the ACC1B, SCA1B, THR1B, GNV1B, LRI1B
  CALL input_openAllFiles(date_str)
  CALL output_openAllFiles(rmsnam)
  ! Read in all needed data from input files.
  CALL status_update('STATUS',calling_prog,'gracefit',' ',"Reading GTORB files",0)
  IF (is_h5 )THEN
     CALL input_readGTORBh5_Amat(satics_t,apr_ScaleBias,GPSantoff,rvec,Amat, apr_prm &
          ,sciframe_quat,srf2trf_rotmat,srf2trf_deriv_rotmat,0, mcon_ocean &
          ,msc_tide_amp, mcon_tides)
  ELSE
     CALL input_readGTORB("direct   ",irec,satics_t,apr_ScaleBias,GPSantoff,rvec,apr_prm &
          ,sciframe_quat,srf2trf_rotmat,srf2trf_deriv_rotmat,0,mcon_ocean,msc_tide_amp,mcon_tides)
  ENDIF

  
  ! PT201013: use the new Amat matrix in place of the old "part" matrix
  IF(nkobs_t > 0)call kb_computePartials_Amat(rvec, Amat) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the L1B accelerometer and thrust data. Thrust data must be read first, since it is used in the fitting
! of the model to linearise the accelerometer data (inside input_readACC)
  CALL input_readTHR1B(thrusters_off)
  CALL input_readACC1B(thrusters_off)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF(ngobs_t > 0) CALL input_readGNV1B(gvec,num_gps_used,n_gnv)  

  ! 190924 SA,HM. leadsat 4th args is nepochs_t not num_gps_used  
  CALL leadsat(mission,gvec, 6, nepochs_t, nsat_t, lead)

  ! PT140501: fix bug that kbrange was never initialised if no kbr kbrr or kbra obs were to be used
  IF(nkobs_t > 0) THEN
    CALL input_readKBR1B(kbrange,kblt,kbant) ! PT201014: pass in the array to fill out whether to use obs or not
    ! PT170831: perform numerical differentiation of kbrange to generate kbrr and kbra. Need to do that for kbant as well (kblt is correct and doesn't need fixing)
    ! SA Comment call
    ! SA to make achraf happy, set to false
    if (use_kb_nr(1) == 1)then
      kb_differentiate_flag = .true.
    else
      kb_differentiate_flag = .false.
    endif
    do iepoch = 1 , nepochs_t
      if(kbrange(2, iepoch) == 0)  then
         use_obs(iepoch, 1:3) = .false.
      endif
    enddo 
    call kb_differentiate_kbr(kb_differentiate_flag,kbrange,kblt,kbant)
  ELSE
    kbrange = 0.d0
  ENDIF

  ! PT150901: fix bug to set number of tidal amplitudes to zero if we don't want to estimate them, even if found in GTORB files
  IF(est_msc_tides == 0 )nmsc_tid_constit_t = 0

  ! PT161128: move the opening of the regularization file up higher in the program
  ! PT180723: only read the mascon reg file if a) we want to estimate mascons, b) we want to regularise the inversion
  IF(nmascons_t > 0 .and. use_msc_len_const == 1)THEN
    ! PT190905: need to rewind the first line that was read when the file was opened.
    rewind(243)  
    IF(ioerr /= 0)THEN
      CALL status_update('FATAL',calling_prog,'gracefit',msc_const_file(1:trimlen(msc_const_file)) &
             ,'Error opening msc regularization',ioerr)
    ENDIF
    ! PT170206: verify that the mascon file code is compatible with that of the GTORB files     
    !READ(243,'(a8)')msc_reg_code
    ! RM210415: store header line and extract code so that I can test reg file type later
    READ(243,'(a)')header_line
    msc_reg_code = header_line(1:8)

    ! PT170713: check whether it is a VCV file to be used for postfit residual computations
    IF(msc_reg_code(1:1) == "V")THEN
      VCV_flag = .TRUE.
      CALL status_update('STATUS',calling_prog,'gracefit',msc_const_file &
             ,'Will compute postfit residuals using solution in VCV file',0)
    ELSE IF(msc_reg_code /= msc_hdr_code)THEN
      WRITE(message,'(5a)')"Incompatible temporal mascon file (code ",msc_hdr_code,") and mascon regularization file (code " &
             ,msc_reg_code,"). Should not continue"
      CALL status_update('WARNING',calling_prog,'gracefit',msc_const_file,message,0)
    ENDIF
  ENDIF
  !****************************************************************

  ! precompute the yaw, pitch, roll angles and the antenna offset correction (AOC) of the KBR
  IF(nsat_t == 2) CALL rpy_AOC(rvec(1:3,1,:),rvec(1:3,nsat_t,:),sciframe_quat(:,1,:),sciframe_quat(:,nsat_t,:),rpy,AOC) !@# condition...

  ! Initialize necessary arrays
  CALL status_update('STATUS',calling_prog,'gracefit',' ','Initialising arrays',0)
  pre_omc = 0.d0
  kbrr_misfit_num = 0
  shadow_counter = 0
  apply_shadow_cond = .FALSE.


!****************************************************************
!*  PT201109: read outlier removal file and flag unwanted epochs/obs      *
  read(c_yr,*)date(1)
  read(c_month,*)date(2)
  read(c_day,*)date(3)
  call mask_outliers(calling_prog,use_outlier_mask,nepochs_t,date(1:3),use_obs)
! PT220524: also mask out obs that lie outside of start_new:end_neq
  if(start_neq > 1)then
    write(message,'(a,i5,a)')"Masking out obs before requested start epoch (1:",start_neq-1,")"
    CALL status_update('STATUS',calling_prog,'gracefit',' ',message,0)
    use_obs(1:start_neq-1,:) = .false.
  endif
  if(end_neq < nepochs_t)then
    write(message,'(a,i5,a,i5,a)')"Masking out obs after requested end epoch (",end_neq+1,":",nepochs_t,")"
    CALL status_update('STATUS',calling_prog,'gracefit',' ',message,0)
    use_obs(end_neq+1:nepochs_t,:) = .false.
  endif
!****************************************************************


!****************************************************************
!!  PT200706
!!  generate the inertial-to-earth_fixed matrices   
!!
!****************************************************************
  call gsec_to_ymdhms(starting_epoch,date_tmp,sec_tmp)
  call ymdhms_to_jd(date_tmp,sec_tmp,jd)
  mjd_start = jd - 2400000.5d0   
  mjdend = mjd_start + 1.d0
  nepochs =  1 + (mjdend - mjd_start ) * 24.d0*60.d0*12.d0 + 10.d0*12.d0  !(12x5sec epochs per minute for an extra 5 mins at each end of orbit)
  allocate(rot_mats(3,3,nepochs))
  allocate(rotdots(3,3,nepochs))
  allocate(rotaccs(3,3,nepochs))
  allocate(rpy_tmp(3,nepochs))
  allocate(rpydots(3,nepochs))
  allocate(rot_dates(nepochs))
  call generate_inert_efixed(calling_prog,.false.,mjd_start,mjdend,rot_mats,rotdots,rotaccs,rpy_tmp,rpydots,rot_dates,nepochs,0.d0,0.d0,0.d0)


!****************************************************************
!*                                                              *
!******                     EPOCH LOOP                  *********
!*                                                              *
!****************************************************************
  ! Report beginning of loop
  call status_update('STATUS',calling_prog,'gracefit',' ','Loop over epochs to prepare partials',0)

  sumsq = 0.d0

  ! PT170703: add an override of kbr filtering should there be any missing kbr observations. Set to false first.
  override_kbr_filt = .FALSE.
  override_kbr_cos_fft = .FALSE.

  ! Epoch loop
  nbeta = 0 ! RM191112: initialise nbeta counter
  n_kbr_gaps = 0




  do iepoch = 1, nepochs_t
     debug_iepoch = iepoch

     ! get the inertial-efixed rotation matrix for this epoch. Used in shadow and beta angle computations
     mjd = mjd_start + (dble(iepoch-1)*epoch_interval - 0)/86400.d0
     interp_flag = 1   ! 1: only the rot_i2e and rot_e2i; 2: all rotation matrices
     call inert_interp_v2(interp_flag,mjd,nepochs,rot_dates,rot_mats,rotdots,rotaccs,rot_i2e,rotdot_i2e &
                         ,rot_e2i,rotdot_e2i,rotacc_e2i,rotacc_i2e )

     !******************** CALCULATIONS FOR KBAND ********************

     ! Only do calculations if there was a line at corresponding epoch (checked by making sure biased range is not 0.0)
     ! PT201014: check this now using the KBRR value of the use_obs logical array
     if (use_obs(iepoch,2)) then

        ! PT170725: set flag to indicate that we are not in a kbr data gap
        in_kbr_gap = .false.

        ! Calculate pre_omc
        call kb_computeOMC(iepoch,pre_omc(:,iepoch),rvec(1:nrvecprm_t,:,iepoch),kbrange(:,iepoch),kblt(:,iepoch),kbant(:,iepoch),&
             aoc(:,:,iepoch))

     ELSE
        IF(.NOT. in_kbr_gap)THEN  ! the start of a gap in kbr observations
           n_kbr_gaps = n_kbr_gaps + 1
           edge1(n_kbr_gaps) = iepoch
           in_kbr_gap = .TRUE.
        ENDIF

        pre_omc(:,iepoch) = 0.d0
        ! PT170703: set the override to true so that filtering of kbr obs will NOT occur
        override_kbr_filt = .TRUE.
        override_kbr_cos_fft = .TRUE.
        ! push the value of the end of the gap along one more epoch
        edge2(n_kbr_gaps) = iepoch

     ENDIF  ! end of kbrange if statement
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

! PT220510: don't use pos/vel if use_obs mask says not to
       if( use_gvec(iepoch, isatt) .and. use_obs(iepoch,7) .and. use_obs(iepoch,8) ) then 
         ! Compute GPS antenna offset in terrestrial reference frame
         call quat_rot_vect(sciframe_quat(:,isatt,iepoch),GPSantoff(5:7,isatt),GPS_antoff_trf(1:3,isatt))
         call matmult(srf2trf_deriv_rotmat(:,:,iepoch,isatt),GPSant_adj(:,isatt),GPS_antoff_trf(4:6,isatt),3,3,1)

         ! Calculate pre_omc
         call gobs_computeOMC(pre_omc(:,iepoch),gvec(:,:,igvec),rvec(1:ngobs_t,:,iepoch),GPS_antoff_trf)
  
         ! Apply CoM conditions
         CALL add_CoM_cond(part(:,:,iepoch),pre_omc(:,iepoch),rvec(1:nCoM_t,:,iepoch),gvec(:,:,igvec),GPS_antoff_trf,&
             srf2trf_rotmat(:,:,iepoch,:),srf2trf_deriv_rotmat(:,:,iepoch,:))

         ! Compute the beta angle
         nbeta = nbeta+1 ! RM191112
         ! PT210401: add the inertial-to-efixed rotation matrix to the arguments.
         CALL calc_beta_angle(iepoch,rot_i2e,rvec(1:3,1,iepoch),rvec(4:6,1,iepoch),beta_angle(iepoch))
       endif

     enddo ! end of satellite do loop for GNV1B obs
     !****************************************************************

!    !******************* SHADOW COMPUTATIONS SETUP ******************
!
!     ! Set shadow based on whether the satellites are in shadow, partially in shadow or in full sunlight
!     CALL shadow_satelliteStatus(iepoch,beta_angle(iepoch),rvec(1:6,:,iepoch),sciframe_quat(:,:,iepoch),shadow(:,iepoch),sunpos)
!
!     ! Set whether the condition is applicable
!     IF(ncond_t > 0) CALL shadow_conditionApplicable(iepoch,shadow(:,iepoch),thrusters_off,shadow_counter,apply_shadow_cond &
!          , n_shadow_obs)
!
!    !****************************************************************
  ENDDO  ! End of epoch loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  CALL status_update('STATUS',calling_prog,'gracefit','',"Compute KBRA start",0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!******************* RANGE ACCELERATION COMPUTATIONS ******************
! PT170428: in this subroutine we will numerically differentiate the kbrr theoretical range
!           and partials to compute the information for the range accelerations
! PT170725: pass in information on the start/stop epochs of gaps in the kbr obs
!  PT180823: allow additional options for how to apply the NR, and pass the options into the subroutine
  if (use_kb_nr(2) > 0 .and. nkbra_t == 1) then
  nrvec_t = 6
  CALL kb_computeKBRA(use_obs,use_kb_nr(2),ikbrr,ikbra,nepochs_t,nobs_t,nparam_t,nrvec_t,nsat_t,nrvecprm_t,pre_omc,Amat,rvec &
       ,kbrr_part,kbra_part,n_kbr_gaps,kbrange, edge1,edge2)
  endif
   CALL status_update('STATUS',calling_prog,'gracefit','',"End of epoch loop that prepared the partial derivatives",0)
  ! PT/HmcQ160704: write out some shadow constraint information
  WRITE(message,'(a,2i8)')"Number of obs to constrain in shadow: ",n_shadow_obs,ncond_t
  CALL status_update('STATUS',calling_prog,'gracefit','',message,0)
  !********************************************************************************************************************************

! DEBUG
! PT210908: output the partial derivatives d_kbra/d_msc for a selection of mascons to assess what is going on
!do iepoch = 1,nepochs_t
!  print*,iepoch,Amat(ikbra+(iepoch-1)*nobs_t,24+4300:24+4500),' gracefit kbra partials 4300:4500'
!  !print*,iepoch,Amat(ikbra+(iepoch-1)*nobs_t,1:30),' gracefit kbra partials'
!enddo
!stop 'stopped after partials debug'




  !*********** FILTER KBRR/KBRA obs and/or remove a 1pr ?  **********
  IF(override_kbr_filt .AND. INT(use_kbrr_decomp(1)) == -2)THEN
     CALL status_update('STATUS',calling_prog,'gracefit',' ','Filtering of kbr obs will NOT occur because of missing obs',0)
  ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PT170717:   filtered KBRR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PT201113: save a copy of the prefit unfiltered kbrr residuals. Use it later to compute postfit unfiltered KBRR residuals
  IF(nkbrr_t > 0)pre_omc_kb_unfilt(1,:) = pre_omc(iKBRR,:)
  IF(nkbra_t > 0)pre_omc_kb_unfilt(2,:) = pre_omc(iKBRA,:)


  ! here we filter the kbr obs. use_kbrr_decomp(1) = -2 in the case where we want to filter but only if there are no data gaps
  override_kbr_filt = .FALSE.
  IF( (INT(use_kbrr_decomp(1)) > 0 .OR. INT(use_kbrr_decomp(1)) == -2) .AND. .NOT. override_kbr_filt)THEN  ! we remove the high-frequency noise from all options apart from 0 (leave as is)

     IF(nkbrr_t > 0) then

! PT190501: added an alternate way of using EMD for when the KBR data are segmented (2016/2017 GRACE data)
      if(INT(use_kbrr_decomp(1)) == 5)then
        CALL kb_filtered_v2(calling_prog,'kbrr',nepochs_t,start_neq,end_neq,NINT(use_kbrr_decomp(2)),pre_omc(iKBRR,:) &
                            ,filt_kbrr_omc)
      else
        CALL kb_filtered('kbrr',nepochs_t,start_neq,end_neq,NINT(use_kbrr_decomp(2)),pre_omc(iKBRR,:),filt_kbrr_omc)
      endif
     endif
     IF(nkbra_t > 0) then
      if(INT(use_kbrr_decomp(1)) == 5)then
        filt_kbra_omc = 0.d0
        CALL kb_filtered_v2(calling_prog,'kbra',nepochs_t,start_neq,end_neq,NINT(use_kbrr_decomp(2)),pre_omc(iKBRA,:) &
                            ,filt_kbra_omc)
      else
        CALL kb_filtered('kbra',nepochs_t,start_neq,end_neq,NINT(use_kbrr_decomp(2)),pre_omc(iKBRA,:),filt_kbra_omc)
      endif
     endif

     ! generate a model with a 1/rev signal with linearly changing amplitude removed ?
     kb_1pr = 0.d0
     IF(use_kbrr_decomp(1) == 3)THEN
        CALL kb_sponge(nepochs_t,filt_kbrr_omc(:),kb_1pr,kb_sponge_model)  ! remove kbrr sponge as well
     ELSE
        CALL status_update('STATUS',calling_prog,'gracefit',' ','Not removing kbrr empirical model',0)
        kb_1pr = 0.d0
     ENDIF

     ! replace the kbrr prefit omc with the filtered one, and with the 1/rev model removed
     IF(nkbrr_t >0 .AND. INT(use_kbrr_decomp(1)) > 0 .AND. .NOT. override_kbr_filt)THEN
       DO i=1, nepochs_t
           pre_omc(ikbrr,i) = filt_kbrr_omc(i) - kb_1pr(i)
       ENDDO
       write(1002,*) pre_omc(iKBRR,:)
     ENDIF

     ! write out filtered PREFIT RESIDUALS
     CALL status_update('STATUS',calling_prog,'gracefit',' ','Writing out filtered prefit residuals',0)
     OPEN(397,file='prefit.res',status='unknown')
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
        OPEN(397,file='prefitRA_filt.res',status='unknown')
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

     ! replace the kbra prefit omc with the filtered one, and with the 1/rev model removed
     !FIXME: SA Shouldn't be beofre the writting part?
     IF(nkbra_t >0  ) THEN
       DO i=1, nepochs_t
         pre_omc(ikbra,i) = filt_kbra_omc(i)
       ENDDO
     ENDIF

     !!! PT/RMcG190722: invoke the cosine fft on the kband obs if requested
  else if (nint(use_kb_cos_fft(1)) == 1) then   ! invoke the cosine fft filter on the kband observations

     kb_1pr = 0.d0
     nfft = 2**16
     IF(nkbrr_t > 0) then
       call kb_cos_fft(calling_prog,use_obs(:,2),'kbrr',nepochs_t,nfft,start_neq,end_neq,use_kb_cos_fft(2:3),pre_omc(iKBRR,:) &
                       ,filt_kbrr_omc)
       pre_omc(ikbrr,:) = filt_kbrr_omc(:)
     endif
     IF(nkbra_t > 0) then
       call kb_cos_fft(calling_prog,use_obs(:,3),'kbra',nepochs_t,nfft,start_neq,end_neq,use_kb_cos_fft(2:3),pre_omc(iKBRA,:) &
                       ,filt_kbra_omc)
       pre_omc(ikbra,:) = filt_kbra_omc(:)
     endif

     ! write out filtered PREFIT RESIDUALS
     CALL status_update('STATUS',calling_prog,'gracefit',' ','Writing out cosine fft filtered prefit residuals',0)
     OPEN(397,file='prefit.res',status='unknown')
     WRITE(397,'(a,a)')" iepoch :   kbrr  in mm/s  :    latitude     longitude  :"  &
          ,":     XYZ (GRACE A) in m      :      XYZvel (GRACE A) in mm/s      " &
          ,":     XYZ (GRACE B) in m  :      XYZvel (GRACE B)   in mm/s. Cosine FFT filtering."
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
        OPEN(397,file='prefitRA_filt.res',status='unknown')
        WRITE(397,'(a,a)')" iepoch :   kbrr  in mm/s  :    latitude     longitude  :"  &
             ,":     XYZ (GRACE A) in m      :      XYZvel (GRACE A) in mm/s      " &
             ,":     XYZ (GRACE B) in m  :      XYZvel (GRACE B)   in mm/s. Cosine FFT filtering. "
        DO iepoch = 1,nepochs_t
           satlon = datan2(rvec(2,1,iepoch),rvec(1,1,iepoch)) *180.d0/pi
           satlat = dasin(rvec(3,1,iepoch)/dsqrt(rvec(1,1,iepoch)**2 + rvec(2,1,iepoch)**2 + rvec(3,1,iepoch)**2))*180.d0/pi
           ! filtered kbra output
           WRITE(397,'(i6,3f14.8,13f12.5,f18.8)')iepoch,filt_kbra_omc(iepoch)*1.d6  &
                ,satlat, satlon &
                ,pre_omc(1:3,iepoch) &
                ,pre_omc(4:6,iepoch)*1.d3,pre_omc(7:9,iepoch),pre_omc(10:12,iepoch)*1.d3,kb_1pr(iepoch)*1.d6 &
                ,pre_omc(iKBRA,iepoch)*1.d6

        ENDDO
        CLOSE(397)
     ENDIF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! PT170717: unfiltered KBRR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! PT170717: also use unfiltered obs if the filter option of -2 was chosen (i.e. filter but not if there are data gaps)
  ELSE IF ( INT(use_kbrr_decomp(1)) == 0 .OR. (INT(use_kbrr_decomp(1)) == -2 .AND. override_kbr_filt) ) THEN

     kb_1pr = 0.d0

     ! write out the prefit residuals at each epoch, so that I don't have to run gracefit all the way through just to get this information
     CALL status_update('STATUS',calling_prog,'gracefit',' ','Writing out (unfiltered) prefit residuals',0)
     OPEN(397,file='prefit.res',status='unknown')
     WRITE(397,'(a,a,a)')" iepoch : kbrr  (um/s):   latitude     longitude  "  &
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
        OPEN(397,file='prefitRA_unfilt.res',status='unknown')
        WRITE(397,'(a,a)')" iepoch :   kbrr  in mm/s  :    latitude     longitude  :"  &
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
  CALL status_update('STATUS',calling_prog,' ',' ',message,0)
  WRITE(message,'(a,3f8.2,3f10.4,a)')'Prefit RMS GRACE B obs (pos/vel):',(sumsq(j+6)*1.d3,j=1,6),' (mm, mm/s)'
  CALL status_update('STATUS',calling_prog,' ',' ',message,0)
  IF(nkbrr_t > 0)THEN
  WRITE(message,'(a,f10.4,a)')'Prefit RMS    KBRR obs          :',sumsq(13)*1.d6,' (um/s)'
  CALL status_update('STATUS',calling_prog,' ',' ',message,0)
  ENDIF
  IF(nkbra_t > 0)THEN
     WRITE(message,'(a,f10.4,a)')'Prefit RMS    KBRA obs          :',sumsq(14)*1.d9,' (nm/s^2)'
     CALL status_update('STATUS',calling_prog,' ',' ',message,0)
  ENDIF

! PT191126: stop here if requested by user via the "gracefit_step" command file option.
  if(gracefit_step == 1)then
    call status_update('STATUS',calling_prog,'gracefit',' ',"Program stopped by user after computing prefit residuals",0)
    stop
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Second epoch loop
  WRITE(message,'(a,i7,a,i7)')'Start stacking normal equations: '
  CALL status_update('STATUS',calling_prog,'gracefit',' ',message,0)
!  call output_3d("gfit_mem.nc",  size(part, 1), size(part, 2), size(part,3),   part)
! DEBUG
!  call output_2d("preomc.nc",  size(pre_omc, 1), size(pre_omc, 2), pre_omc)
! write(*,*), ' KBR ', ikbr, ' KBRR ', ikbrr, ' KBRA ', ikbra

  ! PT170713: we branch here to either read a VCV file to allow the computation of postfit residuals or we form the
  !           normal equations and perform the least squares inversion. The flag for one or the other comes from
  !           whether the first character in the so-called "mascon regularisation" file is a "V" (ie VCV file) or a "#" (ie regularisation file)
!  IF(VCV_flag)THEN   ! it is a VCV file
! read the first mascon solution file
!    if(msc_reg_code(4:8) == "GRACE")then
      ! PT190313: reduced the number of arguments passed to read_soln_v3 to make compatible with the subroutine
!      CALL read_soln_v3("GRACEFIT",243,.false.,sTid,sScl, sBias, sEmp, sMasc, n_emp, nScl, nBias)
     !else if(msc_reg_code(4:6) == "ADD")then 
       !allocate(apriori(nparam_t))
       !allocate(soln(nparam_t))
       !allocate(VCV_obs_local(nparam_t,nparam_t))
       !call read_soln_addnorm(calling_prog," VCV ", 1,243, maxparam, nparam_t, apriori, soln, prm_input, VCV_obs_local, .false. , sTid, sScl, sBias, sEmp, sMasc, n_emp, nScl, nBias)
!     else
!       print*,"UH OH:",msc_reg_code(4:8),msc_reg_code(4:7)
!     endif
     ! make the adjust vector here rather than computing it in a least squares solution
!     adjust = soln - apriori
     !do i=nparam_t-50,nparam_t
     !  print*,"DEBUG! i,adjust",i,adjust(i)
     !enddo
  !ELSE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !
     !  N O R M A L    E Q U A T I O N    S T A C K I N G    L O O P
     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT220913: we only want to create and solve the normal equations if we are NOT using the VCV_flag for postfit residuals
  if(.not. VCV_flag)then
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


    CALL status_update('STATUS',calling_prog,'gracefit',''," ... Masking out observations to be ignored ",0) 
    ! Quick and Dirty check.    
    ! PT220124: fix this so it only uses obs between start_neq and end_neq
    do concurrent (i = 1:nepochs_t)
        if (i < start_neq .or. i > end_neq)then
               amat((i-1)*nobs_t+ikbrr,:) = 0.0
               pre_omc(ikbrr, i) = 0.0
               amat((i-1)*nobs_t+ikbra,:) = 0.0
               pre_omc(ikbra, i) = 0.0
               amat((i-1)*nobs_t + 1: (i-1)*nobs_t + 6, : )  = 0.0
               pre_omc(1:6, i ) = 0.0
               amat((i-1)*nobs_t + 7: (i-1)*nobs_t + 12, : )  = 0.0
               pre_omc(7:12, i ) = 0.0
        endif
        
        if ( .not. use_obs(i,2 )) then 
               amat((i-1)*nobs_t+ikbrr,:) = 0.0
                pre_omc(ikbrr, i) = 0.0
                !print *, "KBRR", ikbrr, "epoch ", i , "set to 0"
        endif
        if (  .not. use_obs(i, 3)) then
                 amat((i-1)*nobs_t+ikbra,:) = 0.0
                 pre_omc(ikbra, i) = 0.0
                !print *, "KBRA", ikbra, "epoch ", i , "set to 0"
        endif

        if ( .not. use_gvec(i, 1) ) then
                amat((i-1)*nobs_t + 1: (i-1)*nobs_t + 6, : )  = 0.0
                pre_omc(1:6, i ) = 0.0
        endif

         if ( .not. use_gvec(i, 2) ) then
                 amat((i-1)*nobs_t + 7: (i-1)*nobs_t + 12, : )  = 0.0
                 pre_omc(7:12, i ) = 0.0
         endif
    enddo


    CALL status_update('STATUS',calling_prog,'gracefit',''," ... ... ADD weights to diagonal elements",0)
    ! call output_2d("amat.nc", size(Amat, 1), size(Amat,2),  Amat)
    do i = 1,nobs_t*nepochs_t
       amat(i,:) = amat(i,:) * dsqrt(apr_wght(i))
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL status_update('STATUS',calling_prog,'gracefit',''," ... ... ... form the AtWA matrix",0)
    call DSYRK('U', 'T', nparam_t , nobs_t*nepochs_t, 1.0d0, amat, nobs_t*nepochs_t, 0.0d0, normeq, nparam_t )

     ! build the b vector:
    CALL status_update('STATUS',calling_prog,'gracefit',''," ... ... ... ... build the AtWb vector",0)
     do concurrent (i = 1:nepochs_t, j = 1: nobs_t)
       bvec(j+(i-1)*nobs_t) = pre_omc(j,i) * dsqrt(apr_wght(j+(i-1)*nobs_t))
     enddo



!     bvec = bvec * dsqrtapr_wght
     CALL DGEMV('T', nobs_t*nepochs_t, nparam_t, 1.0d0, amat, nobs_t*nepochs_t, bvec, 1, 0.0d0, Atwb, 1 )

     
     !****************************************************************
     ! PT221123: turn off the writing out of the binary .norm file. It is now replaced with the .nc version
     CALL output_writeNORM(rmsnam, apr_prm, normeq, AtWb, apr_wght, mcon_tides,prmnam,prmnam_size)

     ! Pt220922: write out normal equations in netcdf format
     allocate(duration(1))
     duration = dble((end_neq-start_neq)*epoch_interval)/86400.d0
     iwhere = index(rmsnam,'.')
     netcdf_out = rmsnam(1:iwhere)//"nc"
     epoch_netcdf(1:3) = date(1:3)
     epoch_netcdf(4) = 0
     epoch_netcdf(5) = 0
     epoch_netcdf(6) = 0   
     ! PT221122: convert prmnam into parameter types (the shorthand storage)
     allocate(param_type(nparam_t))
     if(norbprm_t > 0)IC_solve = .true.
     if(nmascons_t > 0)msc_solve = .true.
     if(est_msc_tides == 1)tmsc_solve = .true.
     msc_offset = 0                               ! PT221123: set to zero if the mascons come after the ICs
     IC_offset  = 0                               ! PT221123: set to zero if the ICs come before the mascons
     call output_codes2labels(1,1,IC_solve,msc_solve,tmsc_solve,nparam_t,nmascons_t,(nparam_t-imsctide),norbprm_t*2 &
          ,param_type,satnam(1),satnam(2),IC_offset,msc_offset,prmnam)

     ! PT221123: only write out the netcdf file if mascons have been estimated
     if(nmascons_t > 0)then
     
       ! need to reorder the AtWB, apriori and param_type to have mascons at the top, then ICs
       allocate(AtWb_flip(nparam_t))
       allocate(apriori_flip(nparam_t))
       allocate(param_type_flip(nparam_t))
       call reorder_normeq(calling_prog,nparam_t,nmascons_t,AtWb,apr_prm,param_type &
                            ,AtWb_flip,apriori_flip,param_type_flip)

              
       ! PT221123: norbprm_t is number of IC parameters per satellite, satnam is the char*1 name for the satellites (for GRACE A,B,C,D)
       call status_update('STATUS',calling_prog,'gracefit',' ','writing out netcdf version of normal equations',0)     
       call write_norm_netcdf_v2(calling_prog,netcdf_out,2.0,nparam_t,nmascons_t,norbprm_t*2,nmsc_tid_constit_t,1 &
                             , normeq(imascons:imascons+nmascons_t-1,imascons:imascons+nmascons_t-1) &
                             , normeq(iICparm:imascons-1,iICparm:imascons-1) &
                             , normeq(iICparm:imascons-1,imascons:imascons+nmascons_t-1)  &
                             , AtWb_flip, param_type_flip,apriori_flip &
                             , netcdf_out,epoch_netcdf,duration,satnam &
                             , .true.,norbprm_t*2,nmascons_t)
     else
       call status_update('STATUS',calling_prog,'gracefit',' ' &
                          ,"From 2022-11-23 onwards, normal equations NOT written out for just IC solutions",0)                                 
     endif                              
     !****************************************************************


     
     ! PT191126: stop here if requested by user via the "gracefit_step" command file option.
     if(gracefit_step == 2)then
       call status_update('STATUS',calling_prog,'gracefit',' ',"Program stopped by user after outputting normal equations",0)
       stop
     endif
  endif     ! end of the stacking of the normal equations
  
     ! RM210416: Move VCV_flag stuff to here, after the bvec is built and weights are incorporated to the amat, 
     ! this is important for computation of the postfit residuals
   IF(VCV_flag)THEN
     if(msc_reg_code(4:8) == "GRACE")then ! gracefit type vcv file
       CALL status_update('STATUS',calling_prog,'gracefit',msc_const_file,"Reading gracefit type VCV file",0)
       ! RM210416: deallocate apriori and soln, these are reallocated in read_soln_v3
       ! PT220913: they are not allocated, so no need to deallocate them
       !deallocate(apriori)
       !deallocate(soln)
       CALL read_soln_v3("GRACEFIT",243,.false.,sTid,sScl, sBias, sEmp, sMasc, n_emp, nScl, nBias)
       adjust = soln - apriori

     else if(msc_reg_code(4:6) == "ADD")then ! addnorm type vcv file
       CALL status_update('STATUS',calling_prog,'gracefit',msc_const_file,"Reading addnorm type VCV file",0)
       READ(243,'(79x,i2,27x)')nfiles_vcv
       READ(243,'(a)')message
       allocate(dates_vcv(nfiles_vcv,3))
       do i=1,nfiles_vcv
         read(243,'(19x,i4,3x,i2,3x,i2,a)')dates_vcv(i,1),dates_vcv(i,2),dates_vcv(i,3),message
         if(dates_vcv(i,3) == date(3))nfile = i
       enddo
       ! determine which no. file we need to extract IC's for...
       allocate(VCV_obs_local(nparam_t,nparam_t))
       maxparam = nparam_t
       call read_soln_addnorm(calling_prog," VCV ",nfile,243,maxparam, maxparam, apriori, soln, prm_input, VCV_obs_local, .false. , sTid, sScl, sBias, sEmp, sMasc, n_emp, nScl, nBias)
       adjust = soln - apriori
     else
       CALL status_update('FATAL',calling_prog,'gracefit',msc_const_file,"Unfamiliar VCV file type",0)
     endif
  ELSE 

     !**********************************************************************
     !
     ! A D D   C O N S T R A I N T S   T O   N O R M A L   E Q U A T I O N S
     !
     !**********************************************************************
     !**********************************************************************
     ! a priori parameter constraints
     !**********************************************************************
     ! PT150827: only do the ICs and mascons here. Mascon tidal amplitudes are
     !           done in the next block down
     IF(imsctide > 0) THEN
        DO i = 1, imsctide - 1
           normeq(i,i) = normeq(i,i) + apr_const(i)
        ENDDO
     ELSE
        DO i=1, nparam_t
           normeq(i,i) = normeq(i,i) + apr_const(i)
        ENDDO
     ENDIF
     !**********************************************************************
     ! a priori tidal amplitude parameter constraints
     !**********************************************************************
     ! PT150820:  constrain all ocean tide amplitudes, with different constraints
     !            for each constituent
     IF(est_msc_tides == 1)THEN
        DO i=imsctide,nparam_t
           IF(prmnam(i)(16:18) == "M2s"    )THEN
              normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(1,1)**2
           ELSEIF(prmnam(i)(16:18) == "M2c")THEN
              normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(1,2)**2
           ELSEIF(prmnam(i)(16:18) == "O1s")THEN
              normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(2,1)**2
           ELSEIF(prmnam(i)(16:18) == "O1c")THEN
              normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(2,2)**2
           ELSEIF(prmnam(i)(16:18) == "S2s")THEN
              normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(3,1)**2
           ELSEIF(prmnam(i)(16:18) == "S2c")THEN
              normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(3,2)**2
           ELSEIF(prmnam(i)(16:18) == "K1s")THEN
              normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(4,1)**2
           ELSEIF(prmnam(i)(16:18) == "K1c")THEN
              normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(4,2)**2
           ELSEIF(prmnam(i)(16:18) == "K2s")THEN
              normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(5,1)**2
           ELSEIF(prmnam(i)(16:18) == "K2c")THEN
              normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(5,2)**2
           ENDIF
        ENDDO
     ENDIF


     !**********************************************************************
     !**********************************************************************
     ! Add mascon length-dependent constraint if necessary
     !**********************************************************************
     IF(nmascons_t > 0 .and. use_msc_len_const == 1)THEN
        CALL status_update('STATUS',calling_prog,'gracefit',msc_const_file &
             ,"    Adding length-dependent constraints to mascon parameters",0)

        ! PT170307: read the second line in the regularisation file, which contains the length scales and constraint levels
        READ(243,'(a)')message
        ! PT150721: check that it really was a header line
        IF(message(1:7) == "# scale" .OR. message(1:5) == "# lat")THEN
           CALL status_update('STATUS',calling_prog,'gracefit',' ',message,0)
        ELSE
           WRITE(message,'(a)')'No (or wrong) header info found in mascon spatial constraint file '
           CALL status_update('WARNING',calling_prog,'gracefit',msc_const_file(1:trimlen(msc_const_file)),message,0)
           BACKSPACE(243)
        ENDIF
        ! read the constraint matrix
        ! RM210415: if reg file is diag only just fill in the diagonals of the constrait matrix
        msc_const = 0.d0
        DO i=1,nmascons_t
          if(header_line(61:69) == "diag_only")then
            read(243,*,iostat=ioerr)msc_const(i,i)
          else
            READ(243,*,iostat=ioerr)(msc_const(i,j),j=1,nmascons_t)
          endif
          IF(ioerr /= 0)THEN
            WRITE(message,'(a,i6,a)')"Error reading line",i," of the mascon regularisation file"
            CALL status_update('FATAL',calling_prog,'gracefit',msc_const_file,message,0)
          ENDIF
        ENDDO
        DO i=imascons,imascons+nmascons_t-1
           DO j=imascons,imascons+nmascons_t-1
              normeq(i,j) = normeq(i,j) + msc_const(i-imascons+1,j-imascons+1) 
           ENDDO
        ENDDO
     ELSE
        CALL status_update('STATUS',calling_prog,'gracefit',' ',"NOT adding length-dependent constraints to mascon parameters",0)
     ENDIF

     !**********************************************************************
     !**********************************************************************
     ! Add mascon mass conserving constraint if necessary           PT150721
     !**********************************************************************
     IF(use_msc_conserve_mass == 1)THEN
        CALL status_update('STATUS',calling_prog,'gracefit',' ',"    Adding mass conservation constraints to mascon parameters",0)

        ! the constraint equation is simply: m1+m2+m3+ ..... +m_n = 0
        ! therefore, the partials d/dm? = 1 for each mascon and the terms in the normal equations are just the variance of the observation
        DO i=imascons,imascons+nmascons_t-1
           DO j=imascons,imascons+nmascons_t-1
              normeq(i,j) = normeq(i,j) + 1.0/0.01**2  ! hardwire 1 cm as the sigma for this conditional equation. @#&& fix this later !
           ENDDO
        ENDDO
     ELSE
        CALL status_update('STATUS',calling_prog,'gracefit',' ',"NOT adding mass conservation constraints to mascon parameters",0)
     ENDIF
     !**********************************************************************

     Call status_update('STATUS', calling_prog, 'gracefit', '', "Solving equations", 0 )
     !**********************************************************************
     !   Solve the normal equations for the estimated ICs and parameters
     IF(INT(use_kbrr_decomp(1)) == 0  .OR. use_kbrr_decomp(1) == 2.OR. use_kbrr_decomp(1) == 3 .or. use_kbrr_decomp(1) == 5)THEN  ! solve for all the observations and parameters
        CALL Chol_normSolve(nparam_t,normeq,AtWb,adjust)
     ELSE IF(INT(use_kbrr_decomp(1)) == 1)THEN  ! solve for only the mascons
        CALL Chol_normSolve(nmascons_t,normeq(imascons:imascons+nmascons_t-1,imascons:imascons+nmascons_t-1) &
             ,AtWb(imascons:imascons+nmascons_t-1),adjust(imascons:imascons+nmascons_t-1))
     ENDIF

     !******************************************************************
  ENDIF   ! end of if statement as to whether to perform the LS inversion or read a solution from a VCV file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Compute post fit. 
  Call status_update('STATUS', calling_prog, 'gracefit', '', "Computing postfit residuals", 0 )
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
  CALL output_writeVCV(gt_fnam,apr_prm,kbrr_misfit_num,adjust,normeq,prmnam,prmnam_size)
  CALL output_writeFIT(gt_fnam,apr_prm,apr_wght,adjust,kbrr_misfit_num,Amat, pre_omc,normeq,post_omc,prmnam,prmnam_size &
                       ,msc_const_file)  ! postfit residuals are computed in this subroutine
  !      call output_writeCOR(normeq)
  !****************************************************************

  !   Plot the residuals
  IF( iop > 0 )THEN
     IF(nkbrr_t /= 0) CALL plot_writeKB(iop,pre_omc,post_omc,post_omc_kb_unfilt(iKBRR,:),rvec(1:3,:,:),rpy,kbrr_prefit_tol,nepochs_t)
     IF(nkbra_t /= 0) CALL plot_writeKBA(iop,pre_omc,post_omc,post_omc_kb_unfilt(iKBRA,:),rvec(1:3,:,:),rpy,kbrr_prefit_tol,nepochs_t)
     CALL plot_writePLT(iop,post_omc,rvec(1:6,:,:),gvec(4,:,:))
  ENDIF
  !****************************************************************

  ! ******************!  Output the mean beta angle *************************
! PT/HMcQ190301: divide by the number of epochs used, not the number of GPS epochs in the GNV1B file
  WRITE(message,'(a,f8.2,a,f8.2)')"Mean beta angle = ",SUM(beta_angle(:))/DBLE(nbeta/nsat_t)*180.d0/pi," degrees."
  CALL status_update('STATUS',calling_prog,'gracefit',' ',message,0)
  !****************************************************************


  !*********************** DEALLOCATE ARRAYS **********************

  CALL status_update('STATUS',calling_prog,'gracefit',' ','Deallocating arrays',0)
  DEALLOCATE(apr_prm)
  DEALLOCATE(apr_const)
  DEALLOCATE(apr_wght)
  DEALLOCATE(amat)
  DEALLOCATE(normeq)
  DEALLOCATE(AtWb)
  DEALLOCATE(adjust)
  DEALLOCATE(pre_omc)
  DEALLOCATE(post_omc)
  DEALLOCATE(rvec)
  DEALLOCATE(sciframe_quat)
  DEALLOCATE(srf2trf_rotmat)
  DEALLOCATE(srf2trf_deriv_rotmat)
  DEALLOCATE(rpy)
  DEALLOCATE(AOC)
  DEALLOCATE(GPS_antoff_trf)
  DEALLOCATE(GPSant_adj)
  DEALLOCATE(acc)
  DEALLOCATE(gvec)
  !    deallocate(uvec)
  DEALLOCATE(kbrange)
  DEALLOCATE(kblt)
  DEALLOCATE(kbant)
  !    deallocate(kbion_cor)
  !    deallocate(kbsnr)
  DEALLOCATE(thr)
  DEALLOCATE(thrusters_off)
  DEALLOCATE(shadow_counter)
  DEALLOCATE(apply_shadow_cond)
  DEALLOCATE(shadow)
  DEALLOCATE(filt_kbrr_omc)
  DEALLOCATE(filt_kbra_omc)
  !****************************************************************

  !   That's all folks
  call cpu_time(finish)
  WRITE(message,'(a,f15.2,a)')" Normal end of GRACEFIT, computing time  ", finish - start ," seconds."
  CALL status_update('STATUS',calling_prog,'gracefit',' ', message ,0)
END PROGRAM GRACEFIT

!********************************************************************************************************************************
