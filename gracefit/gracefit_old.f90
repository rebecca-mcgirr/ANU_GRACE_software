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
!
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
!********************************************************************************************************************************

PROGRAM GRACEFIT

  ! PT170609: include the mod file to define the record length of the GTORB files
  USE gtorb_mod
  USE gracefit_mod
  use accred_mod      ! declares the variables for the accelerometer observations
  use soln_mod        ! PT190905: adds the declaration of the arrays for the solution (apriori, adjust, soln)
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
  ! PT170206: add mascon filenames and codes
  CHARACTER*40  :: combined_mascon_file,ocean_mascon_file        ! name of temporal and ocean mascon files used in GRACEORB
  CHARACTER*8   :: msc_hdr_code,msc_ocn_hdr_code,msc_reg_code    ! mascon file codes of temporal and ocean mascon files used in GRACEORB
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
  logical                       :: use_gvec          ! logical to indicate whether there are GNV1B obs to match the GTORB obs
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
  DOUBLE PRECISION, ALLOCATABLE :: part(:,:,:)    ! Partials of observables with respect to the parameters
  DOUBLE PRECISION, ALLOCATABLE :: normeq(:,:)    ! The left side of the LS algorithm. P_t*W*P
  DOUBLE PRECISION, ALLOCATABLE :: AtWb(:)        ! The right side of the LS algorithm. P_t*W*OMC
  DOUBLE PRECISION, ALLOCATABLE :: adjust(:)      ! Solution to normal equations
  DOUBLE PRECISION, ALLOCATABLE :: post_omc(:,:)  ! Postfit Observed Minus Computed values
  DOUBLE PRECISION, ALLOCATABLE :: beta_angle(:)  ! beta angle at each GPS epoch

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
  DOUBLE PRECISION, ALLOCATABLE :: normeq_tmp(:,:)    ! The left side of the LS algorithm. P_t*W*P     tmp for openMP implementation
  DOUBLE PRECISION, ALLOCATABLE :: AtWb_tmp(:)        ! The right side of the LS algorithm. P_t*W*OMC  tmp for openMP implementation
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
  INTEGER*4    :: nparam_vcv                                          ! number of parameters found in the VCV file
  INTEGER*4    :: sTid,sScl, sBias, sEmp, sMasc, n_emp, nScl, nBias   ! we don't need these for the postfit residual computations
!  REAL(kind=8),ALLOCATABLE :: soln(:)                                 ! solution vector from the VCV file
!  REAL(kind=8),ALLOCATABLE :: apriori(:)
  CHARACTER*30,ALLOCATABLE :: prm_input(:)

  ! PT140429: debug variables
  REAL(kind=8) :: kbtheoretical(3), IC_adj(12,2), prefit_adj(13),cmd_args(12),sumsq_for_rms(13),rms(13)
  REAL(kind=8) :: delta_kbrr,sunpos
  REAL(kind=8) :: dot,amag3
  REAL(kind=8) :: tmp_LOS_A(3),tmp_LOS_B(3),tmp_LOS_mag,tmp_vel_A(3),tmp_vel_B(3),tmp_vel_mag
  REAL(kind=8), ALLOCATABLE :: kbrr_part(:,:,:),kbra_part(:,:,:)

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
  CALL status_update('STATUS','GRACEFIT','gracefit/gracefit',' ',message,0)

  ! Assert that there are the correct number of command-line arguments
  IF ( iargc() < NUM_EXPECTED_ARGS )THEN
     CALL status_update('FATAL','GRACEFIT','gracefit',' ','Too few arguments in command line',0)
  ENDIF

  !************** PARSE COMMAND LINE AND COMMAND FILE *************
  call get_command (commandline, lencmd, status_cmd)
  CALL status_update('STATUS', 'GRACEFIT', 'gracefit/cmdline', ' ', trim(commandline),0)


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
  CALL command_readSetup()
  write(message,'(a,i2)')"   GRACE mission: ",mission
  call status_update('STATUS','GRACEFIT','gracefit',' ',message,0)

! PT190905: we need the mascon regularisation file irrespective of whether we want to estimate mascons or not because
!           we use it to derive the postfit residuals only
  CALL getarg(9,msc_const_file) ! mascon constraint file name
  VCV_flag = .false.
  if(msc_const_file(1:1) /= "")then  ! a file was supplied on the command line
    OPEN(243,file=msc_const_file(1:trimlen(msc_const_file)),status='old',iostat=ioerr)
    IF(ioerr /= 0)THEN
      CALL status_update('FATAL','GRACEFIT','gracefit',msc_const_file(1:trimlen(msc_const_file)) &
                       ,'Error opening msc regularization',0)
    ENDIF
    ! PT170206: verify that the mascon file code is compatible with that of the GTORB files
    READ(243,'(a8)')msc_reg_code
    ! PT170713: check whether it is a VCV file to be used for postfit residual computations
    IF(msc_reg_code(1:1) == "V")THEN
      VCV_flag = .TRUE.
      CALL status_update('STATUS','GRACEFIT','gracefit',msc_const_file &
                       ,'Will compute postfit residuals using solution in VCV file',0)    
    endif   
    !CLOSE(243)  ! JP191122: close file, opened again with output_openAllFiles at L431! 
  endif


  test = gt_fnam(1)
  length = LEN(TRIM(gt_fnam(1)))
  IF( test(length-2:length) == '.h5') THEN
     call status_update('STATUS','GRACEFIT','gracefit',' ','   Input GTORB files are hdf5 format',0)
     is_h5 = .TRUE.
  ELSE
     call status_update('STATUS','GRACEFIT','gracefit',' ','Input GTORB files are binary format',0)
     is_h5 = .FALSE.
  ENDIF
  IF (is_h5 ) THEN
     CALL input_openh5GTORB(gt_fnam,0)
     !NOTE: input_readGTORB MUST follow directly after command_readSetup (i.e. before anything else is set)
     CALL header_readh5GTORB(irec,satics_t,apr_ScaleBias,GPSantoff,mcon_tides &
          ,combined_mascon_file,msc_hdr_code,ocean_mascon_file,msc_ocn_hdr_code,total_ocean_prim,0,prmnam,prmnam_size)
  ELSE
     CALL input_openGTORB(gt_fnam,0,"direct    ",rec_len)
     !NOTE: input_readGTORB MUST follow directly after command_readSetup (i.e. before anything else is set)
     CALL header_readbinGTORB(rec_len,irec,satics_t,apr_ScaleBias,GPSantoff,0,mcon_tides &
          ,combined_mascon_file,msc_hdr_code,ocean_mascon_file,msc_ocn_hdr_code,total_ocean_prim,prmnam,prmnam_size)
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

  CALL status_update('STATUS','GRACEFIT','gracefit',' ','            Allocating arrays',0)
  ! Variables storing data from command file and GTORB file
  ALLOCATE(apr_prm(nparam_t))
  ALLOCATE(apr_const(nparam_t))
  ALLOCATE(apr_wght(nobs_t))
  ALLOCATE(apr_tide_ampl_const(max_mcon_tides,2))
  ! Variables used for calculations
  ALLOCATE(part(nobs_t,nparam_t,nepochs_t))
  ALLOCATE(normeq(nparam_t,nparam_t))
  ALLOCATE(AtWb(nparam_t))
  ALLOCATE(pre_omc(nobs_t,nepochs_t))
  ALLOCATE(post_omc(nobs_t,nepochs_t))
  ALLOCATE(normeq_tmp(nparam_t,nparam_t))
  ALLOCATE(AtWb_tmp(nparam_t))
  ! Data read in or calculated from GTORB file(s)
  ALLOCATE(rvec(nrvec_t,nsat_t,nepochs_t))
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
!PT190828: these three (adjust, soln, apriori) should be defined in soln_mod
! PT190905: allocate these lower down, once we know if it is a full gracefit inversion or just computation of postfit residuals
  ALLOCATE(adjust(nparam_t))
!  ALLOCATE(soln(nparam_t))
!  ALLOCATE(apriori(nparam_t))
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

! DEBUG: PT190404 for a test ONLY, apply a 10 mm radial offset to the centre of mass of GRACE A - just to see what happens
!  GPSant_adj(1,1) = -1.100d0

  ! Write out date_str to be used for finding data files
  WRITE(date_str,'(a4,a1,a2,a1,a2)')c_yr,"-",c_month,"-",c_day


  !****************************************************************

  ! Open input and output files
  ! PT180625: this includes the ACC1B, SCA1B, THR1B, GNV1B, LRI1B
  CALL input_openAllFiles(date_str)
  CALL output_openAllFiles(rmsnam)
  ! Read in all needed data from input files.

  CALL status_update('STATUS','GRACEFIT','gracefit',' ',"Reading GTORB files",0)
  IF (is_h5 )THEN
     CALL input_readGTORBh5(satics_t,apr_ScaleBias,GPSantoff,rvec,apr_prm &
          ,sciframe_quat,srf2trf_rotmat,srf2trf_deriv_rotmat,0, mcon_ocean &
          ,msc_tide_amp, mcon_tides)
  ELSE
     CALL input_readGTORB("direct   ",irec,satics_t,apr_ScaleBias,GPSantoff,rvec,apr_prm &
          ,sciframe_quat,srf2trf_rotmat,srf2trf_deriv_rotmat,0,mcon_ocean,msc_tide_amp,mcon_tides)
  ENDIF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the L1B accelerometer and thrust data. Thrust data must be read first, since it is used in the fitting
! of the model to linearise the accelerometer data (inside input_readACC)
  CALL input_readTHR1B(thrusters_off)
  CALL input_readACC1B(thrusters_off)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF(ngobs_t > 0) CALL input_readGNV1B(gvec,num_gps_used,n_gnv)  
! DEBUG
!do i=1,100
!  print*,'gracefit gvec:',i,gvec_epochs(i,1),gvec(1:6,1,i)
!enddo

  ! 190924 SA,HM. leadsat 4th args is nepochs_t not num_gps_used  
  CALL leadsat(mission,gvec, 6, nepochs_t, nsat_t, lead)

  ! PT140501: fix bug that kbrange was never initialised if no kbr kbrr or kbra obs were to be used
  IF(nkobs_t > 0) THEN
     CALL input_readKBR1B(kbrange,kblt,kbant)
     ! PT170831: perform numerical differentiation of kbrange to generate kbrr and kbra. Need to do that for kbant as well (kblt is correct and doesn't need fixing)
     ! SA Comment call
     ! SA to make achraf happy, set to false
     if (use_kb_nr(1) == 1)then
       kb_differentiate_flag = .true.
     else
       kb_differentiate_flag = .false.
     endif 
     call kb_differentiate_kbr(kb_differentiate_flag,kbrange,kblt,kbant)
  ELSE
     kbrange = 0.d0
  ENDIF

  ! PT150901: fix bug to set number of tidal amplitudes to zero if we don't want to estimate them, even if found in GTORB files
  IF(est_msc_tides == 0 )nmsc_tid_constit_t = 0

  ! PT161128: move the opening of the regularization file up higher in the program
  ! PT180723: only read the mascon reg file if a) we want to estimate mascons, b) we want to regularise the inversion
  IF(nmascons_t > 0 .and. use_msc_len_const == 1)THEN
     ! PT190905: the file should already be open ... ? 
     !OPEN(243,file=msc_const_file(1:trimlen(msc_const_file)),status='old',iostat=ioerr)
     rewind(243)  ! PT190905: need to rewind the first line that was read when the file was opened.
     IF(ioerr /= 0)THEN
        CALL status_update('FATAL','GRACEFIT','gracefit',msc_const_file(1:trimlen(msc_const_file)) &
             ,'Error opening msc regularization',ioerr)
     ENDIF
     ! PT170206: verify that the mascon file code is compatible with that of the GTORB files     
     READ(243,'(a8)')msc_reg_code
     ! PT170713: check whether it is a VCV file to be used for postfit residual computations
     IF(msc_reg_code(1:1) == "V")THEN
        VCV_flag = .TRUE.
        CALL status_update('STATUS','GRACEFIT','gracefit',msc_const_file &
             ,'Will compute postfit residuals using solution in VCV file',0)

     ELSE IF(msc_reg_code /= msc_hdr_code)THEN
        WRITE(message,'(5a)')"Incompatible temporal mascon file (code ",msc_hdr_code,") and mascon regularization file (code " &
             ,msc_reg_code,"). Should not continue"
        CALL status_update('WARNING','GRACEFIT','gracefit',msc_const_file,message,0)
     ENDIF
  ENDIF
  !****************************************************************

  ! precompute the yaw, pitch, roll angles and the antenna offset correction (AOC) of the KBR
  IF(nsat_t == 2) CALL rpy_AOC(rvec(1:3,1,:),rvec(1:3,nsat_t,:),sciframe_quat(:,1,:),sciframe_quat(:,nsat_t,:),rpy,AOC) !@# condition...

  ! Initialize necessary arrays
  CALL status_update('STATUS','GRACEFIT','gracefit',' ','Initialising arrays',0)
  pre_omc = 0.d0
  part    = 0.d0
  kbrr_misfit_num = 0
  shadow_counter = 0
  apply_shadow_cond = .FALSE.


  !****************************************************************
  !*                                                              *
  !******                     EPOCH LOOP                  *********
  !*                                                              *
  !****************************************************************
  ! Report beginning of loop
  call status_update('STATUS','GRACEFIT','gracefit',' ','Loop over epochs to prepare partials',0)

  sumsq = 0.d0

  ! PT170703: add an override of kbr filtering should there be any missing kbr observations. Set to false first.
  override_kbr_filt = .FALSE.
  override_kbr_cos_fft = .FALSE.

  ! Epoch loop
  nbeta = 0 ! RM191112: initialise nbeta counter
  n_kbr_gaps = 0
  do iepoch = 1, nepochs_t
     debug_iepoch = iepoch
     !******************** CALCULATIONS FOR KBAND ********************

     ! Only do calculations if there was a line at corresponding epoch (checked by making sure biased range is not 0.0)
     if (kbrange(iRANGE,iepoch) /= 0.d0) then

        ! PT170725: set flag to indicate that we are not in a kbr data gap
        in_kbr_gap = .false.

        ! Calculate pre_omc
        call kb_computeOMC(iepoch,pre_omc(:,iepoch),rvec(1:nrvecprm_t,:,iepoch),kbrange(:,iepoch),kblt(:,iepoch),kbant(:,iepoch),&
             aoc(:,:,iepoch))

        ! Compute and place necessary partials into part matrix
        ! PT140511: pass back from kb_computePartials the dRR/dX, dRR/dY etc in new variable kbrr_part(6,2). Use this later in debug stuff.
        call kb_computePartials(iepoch,part(:,:,iepoch),rvec(:,:,iepoch),kbrr_part(:,:,iepoch),kbra_part(:,:,iepoch))
        do i=1,6
           IC_adj(i,1) = cmd_args(i)
           IC_adj(i,2) = cmd_args(i+6)
        enddo
        prefit_adj = 0.d0
        do i=1,3
           prefit_adj(ikbrr) = prefit_adj(ikbrr) + part(ikbrr,i,iepoch)*IC_adj(i,1) &
                +part(ikbrr,i+norbprm_t,iepoch)*IC_adj(i,2)  ! position contribution
           prefit_adj(ikbrr) = prefit_adj(ikbrr) + part(ikbrr,i+3,iepoch)*IC_adj(i+3,1) &
                +part(ikbrr,i+norbprm_t,iepoch)*IC_adj(i+3,2)
        enddo
        sumsq_for_rms(13) = sumsq_for_rms(13) + prefit_adj(13)**2

        ! Validate that the values calculated in omc are within our prefit tollerance !@# Make sure this does what we want it to do
        CALL kb_validateRR(part(ikbrr,:,iepoch),pre_omc(ikbrr,iepoch),kbrr_misfit_num,kbrr_prefit_tol,iepoch)
     ELSE
        IF(.NOT. in_kbr_gap)THEN  ! the start of a gap in kbr observations
           n_kbr_gaps = n_kbr_gaps + 1
           edge1(n_kbr_gaps) = iepoch
           in_kbr_gap = .TRUE.
        ENDIF

        WRITE(message,'(a,i6)')'No KBR obs for epoch',iepoch
! PT180921: turn this off - there are too many messages when GRACE is in eclipse at the end of the GRACE mission
!        CALL status_update('WARNING','GRACEFIT','gracefit',' ',message,0)
        pre_omc(:,iepoch) = 0.d0
        part(:,:,iepoch) = 0.d0
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
! debug
!if(isatt==1)print*,'igvec,iepoch,(iepoch-1)*epoch_interval',igvec,iepoch,(iepoch-1)*epoch_interval+1

       if(starting_epoch + (iepoch-1)*epoch_interval > gvec_epochs(igvec,isatt) ) then
         do while (gvec_epochs(igvec,isatt) < starting_epoch + (iepoch-1)*epoch_interval)
           igvec = igvec + 1
         enddo
         if(gvec_epochs(igvec,isatt) == starting_epoch + (iepoch-1)*epoch_interval)then
           use_gvec = .true.
         else
           use_gvec = .false.
         endif
       endif

       if(use_gvec)then 
         ! Compute GPS antenna offset in terrestrial reference frame
         call quat_rot_vect(sciframe_quat(:,isatt,iepoch),GPSantoff(5:7,isatt),GPS_antoff_trf(1:3,isatt))
         call matmult(srf2trf_deriv_rotmat(:,:,iepoch,isatt),GPSant_adj(:,isatt),GPS_antoff_trf(4:6,isatt),3,3,1)

         ! Calculate pre_omc
         call gobs_computeOMC(pre_omc(:,iepoch),gvec(:,:,igvec),rvec(1:ngobs_t,:,iepoch),GPS_antoff_trf)
  
         ! Calculate partials
         call gobs_computePartials(part(:,:,iepoch),rvec(:,:,iepoch),gvec(4,:,igvec))
  
         ! PT140509: compute prefit residuals corrected for parameter adjustments that are hardwired
         prefit_adj(1:12) = 0.d0
         do i=1,6
           do j=1,6
             prefit_adj(i)   = prefit_adj(i)   + part(i,j,iepoch)  *IC_adj(j,1)
             prefit_adj(i+6) = prefit_adj(i+6) + part(i+6,j,iepoch)*IC_adj(j,2)
           enddo
         enddo
         do i=1,12
           sumsq_for_rms(i) = sumsq_for_rms(i) + prefit_adj(i)**2
         enddo

         ! Apply CoM conditions
         CALL add_CoM_cond(part(:,:,iepoch),pre_omc(:,iepoch),rvec(1:nCoM_t,:,iepoch),gvec(:,:,igvec),GPS_antoff_trf,&
             srf2trf_rotmat(:,:,iepoch,:),srf2trf_deriv_rotmat(:,:,iepoch,:))

         ! Compute the beta angle
         nbeta = nbeta+1 ! RM191112
         CALL calc_beta_angle(iepoch,rvec(1:3,1,iepoch),rvec(4:6,1,iepoch),beta_angle(iepoch))
       endif

     enddo ! end of satellite do loop for GNV1B obs
     !****************************************************************

     !******************* SHADOW COMPUTATIONS SETUP ******************

     ! Set shadow based on whether the satellites are in shadow, partially in shadow or in full sunlight
     CALL shadow_satelliteStatus(iepoch,beta_angle(iepoch),rvec(1:6,:,iepoch),sciframe_quat(:,:,iepoch),shadow(:,iepoch),sunpos)
     ! DEBUG
     !if(iepoch > 0 .and. iepoch < 600 .and. shadow(1,iepoch) < 10.0)print*,'iepoch, shad',iepoch,shadow(1,iepoch),iepoch*5 &
     !             ,starting_epoch,starting_epoch+iepoch*5

     ! Set whether the condition is applicable
!print*,"iepoch,shadow(:,iepoch),thrusters_off(iepoch,:  )",iepoch,shadow(:,iepoch),thrusters_off(iepoch,:)
     IF(ncond_t > 0) CALL shadow_conditionApplicable(iepoch,shadow(:,iepoch),thrusters_off,shadow_counter,apply_shadow_cond &
          , n_shadow_obs)

     !****************************************************************
  ENDDO  ! End of epoch loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !******************* RANGE ACCELERATION COMPUTATIONS ******************
  ! PT170428: in this subroutine we will numerically differentiate the kbrr theoretical range
  !           and partials to compute the information for the range accelerations
  ! PT170725: pass in information on the start/stop epochs of gaps in the kbr obs
!  PT180823: allow additional options for how to apply the NR, and pass the options into the subroutine
  if (use_kb_nr(2) > 0 .and. nkbra_t == 1) then
  CALL kb_computeKBRA(use_kb_nr(2),ikbrr,ikbra,nepochs_t,nobs_t,nparam_t,nrvec_t,nsat_t,nrvecprm_t,pre_omc,part,rvec,kbrr_part &
       ,kbra_part,n_kbr_gaps,kbrange, edge1,edge2)
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  CALL status_update('STATUS','GRACEFIT','gracefit','',"End of epoch loop that prepared the partial derivatives",0)
  ! PT/HmcQ160704: write out some shadow constraint information
  WRITE(message,'(a,2i8)')"Number of obs to constrain in shadow: ",n_shadow_obs,ncond_t
  CALL status_update('STATUS','GRACEFIT','gracefit','',message,0)
  !********************************************************************************************************************************






  !*********** ADD ACCELEROMETER CONDITONS  **********

  ! Condition on mean acceleration of the two satellites
  IF(naobs_t > 0) THEN
     CALL add_meanACC_cond(apply_shadow_cond,acc,apr_wght,apr_ScaleBias,part(:,:,1),pre_omc(:,1))
  ENDIF
  !****************************************************



  !*********** FILTER KBRR/KBRA obs and/or remove a 1pr ?  **********
  IF(override_kbr_filt .AND. INT(use_kbrr_decomp(1)) == -2)THEN
     CALL status_update('STATUS','GRACEFIT','gracefit',' ','Filtering of kbr obs will NOT occur because of missing obs',0)
  ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PT170717:   filtered KBRR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! here we filter the kbr obs. use_kbrr_decomp(1) = -2 in the case where we want to filter but only if there are no data gaps
  override_kbr_filt = .FALSE.
  IF( (INT(use_kbrr_decomp(1)) > 0 .OR. INT(use_kbrr_decomp(1)) == -2) .AND. .NOT. override_kbr_filt)THEN  ! we remove the high-frequency noise from all options apart from 0 (leave as is)

     IF(nkbrr_t > 0) then
! PT190501: added an alternate way of using EMD for when the KBR data are segmented (2016/2017 GRACE data)
      if(INT(use_kbrr_decomp(1)) == 5)then
        CALL kb_filtered_v2(calling_prog,'kbrr',nepochs_t,start_neq,end_neq,NINT(use_kbrr_decomp(2)),pre_omc(iKBRR,:) &
                            ,filt_kbrr_omc)
!debug
!do i=1,nepochs_t
!  print*,i,filt_kbrr_omc(i),pre_omc(iKBRR,i)
!enddo
!stop 'stopped after debug'
      else
        CALL kb_filtered('kbrr',nepochs_t,start_neq,end_neq,NINT(use_kbrr_decomp(2)),pre_omc(iKBRR,:),filt_kbrr_omc)
      endif
     endif
     IF(nkbra_t > 0) then
      if(INT(use_kbrr_decomp(1)) == 5)then
        filt_kbra_omc = 0.d0
        CALL kb_filtered_v2(calling_prog,'kbra',nepochs_t,start_neq,end_neq,NINT(use_kbrr_decomp(2)),pre_omc(iKBRA,:) &
                            ,filt_kbra_omc)
!debug
!do i=1,nepochs_t
!  print*,i,filt_kbrr_omc(i),filt_kbra_omc(i), "kbra debug"
!enddo
!stop 'stopped after debug'
      else
        CALL kb_filtered('kbra',nepochs_t,start_neq,end_neq,NINT(use_kbrr_decomp(2)),pre_omc(iKBRA,:),filt_kbra_omc)
      endif
     endif

!do i=1,nepochs_t
!     write(1000,*) pre_omc(iKBRR,i)
!     write(1001,*) filt_kbrr_omc(i)
!enddo
     !stop 'stopped after kb_filtered'

     ! generate a model with a 1/rev signal with linearly changing amplitude removed ?
     kb_1pr = 0.d0
     IF(use_kbrr_decomp(1) == 3)THEN
        CALL kb_sponge(nepochs_t,filt_kbrr_omc(:),kb_1pr,kb_sponge_model)  ! remove kbrr sponge as well
     ELSE
        CALL status_update('STATUS','GRACEFIT','gracefit',' ','Not removing kbrr empirical model',0)
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
     CALL status_update('STATUS','GRACEFIT','gracefit',' ','Writing out filtered prefit residuals',0)
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

     !!! PT/RMcG190722: invoke the cosine fft on the kband obs if requested
  else if (nint(use_kb_cos_fft(1)) == 1) then   ! invoke the cosine fft filter on the kband observations

! PT190808: turn this off here. Deal with data gaps in the subroutine kb_cos_fft.
!     if ( override_kbr_cos_fft )then
!       call status_update('FATAL','GRACEFIT','gracefit',' ','kb_cos_fft not suitable when missing data exist',0)
!     endif

     kb_1pr = 0.d0
     nfft = 2**16
     IF(nkbrr_t > 0) then
       call kb_cos_fft(calling_prog,'kbrr',nepochs_t,nfft,start_neq,end_neq,use_kb_cos_fft(2:3),pre_omc(iKBRR,:),filt_kbrr_omc)
       pre_omc(ikbrr,:) = filt_kbrr_omc(:)
     endif
     IF(nkbra_t > 0) then
       call kb_cos_fft(calling_prog,'kbra',nepochs_t,nfft,start_neq,end_neq,use_kb_cos_fft(2:3),pre_omc(iKBRA,:),filt_kbra_omc)
       pre_omc(ikbra,:) = filt_kbra_omc(:)
     endif

     ! write out filtered PREFIT RESIDUALS
     CALL status_update('STATUS','GRACEFIT','gracefit',' ','Writing out cosine fft filtered prefit residuals',0)
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
           ! filtered kbra output

        ENDDO
        CLOSE(397)
     ENDIF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! PT170717: unfiltered KBRR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! PT170717: also use unfiltered obs if the filter option of -2 was chosen (i.e. filter but not if there are data gaps)
  ELSE IF ( INT(use_kbrr_decomp(1)) == 0 .OR. (INT(use_kbrr_decomp(1)) == -2 .AND. override_kbr_filt) ) THEN

     kb_1pr = 0.d0

!do i=1,nepochs_t
!  write(1000,*)pre_omc(iKBRR,i)
!enddo

     ! write out the prefit residuals at each epoch, so that I don't have to run gracefit all the way through just to get this information
     CALL status_update('STATUS','GRACEFIT','gracefit',' ','Writing out (unfiltered) prefit residuals',0)
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
  CALL status_update('STATUS','GRACEFIT',' ',' ',message,0)
  WRITE(message,'(a,3f8.2,3f10.4,a)')'Prefit RMS GRACE B obs (pos/vel):',(sumsq(j+6)*1.d3,j=1,6),' (mm, mm/s)'
  CALL status_update('STATUS','GRACEFIT',' ',' ',message,0)
  IF(nkbrr_t > 0)THEN
  WRITE(message,'(a,f10.4,a)')'Prefit RMS    KBRR obs          :',sumsq(13)*1.d6,' (um/s)'
  CALL status_update('STATUS','GRACEFIT',' ',' ',message,0)
  ENDIF
  IF(nkbra_t > 0)THEN
     WRITE(message,'(a,f10.4,a)')'Prefit RMS    KBRA obs          :',sumsq(14)*1.d9,' (nm/s^2)'
     CALL status_update('STATUS','GRACEFIT',' ',' ',message,0)
  ENDIF

! PT191126: stop here if requested by user via the "gracefit_step" command file option.
  if(gracefit_step == 1)then
    call status_update('STATUS',calling_prog,'gracefit',' ',"Program stopped by user after computing prefit residuals",0)
    stop
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  call output_3d("partials__2.nc",  size(part, 1), size(part, 2), size(part,3),   part)

!  call exit(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Second epoch loop
  WRITE(message,'(a,i7,a,i7)')'Start stacking normal equations: '
  CALL status_update('STATUS','GRACEFIT','gracefit',' ',message,0)

! DEBUG
!  call output_2d("preomc.nc",  size(pre_omc, 1), size(pre_omc, 2), pre_omc)
! write(*,*), ' KBR ', ikbr, ' KBRR ', ikbrr, ' KBRA ', ikbra

  ! PT170713: we branch here to either read a VCV file to allow the computation of postfit residuals or we form the
  !           normal equations and perform the least squares inversion. The flag for one or the other comes from
  !           whether the first character in the so-called "mascon regularisation" file is a "V" (ie VCV file) or a "#" (ie regularisation file)
  IF(VCV_flag)THEN   ! it is a VCV file
! PT190313: reduced the number of arguments passed to read_soln_v3 to make compatible with the subroutine
     CALL read_soln_v3("GRACEFIT",243,.false. &
          ,sTid,sScl, sBias, sEmp, sMasc, n_emp, nScl, nBias)

     ! make the adjust vector here rather than computing it in a least squares solution
     adjust = soln - apriori

! DEBUG: print out the first 10 mascon adjustment values
!     print*,'first 10 mascon adjusts:',adjust(1:35)
!     print*,'first 10 mascon soln:',soln(1:35)
!     print*,'first 10 mascon apriori:',apriori(1:35)

  ELSE
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
     IF(INT(use_accel_scale_const) == 1) CALL shadow_SRP(thrusters_off,acc,rvec,srf2trf_rotmat,normeq,AtWb,apr_ScaleBias)
     ! DEBUG
     !    do i=1,nparam_t
     !      write(*,*)i,(normeq(i,j),j=1,nparam_t)
     !    enddo
     !    do i=1,nparam_t
     !      print*,i,AtWb(i)
     !    enddo
     !    stop
     !********************************************************************************************************************************
     !    call shadow_SRPdotprod(thrusters_off,acc,rvec,srf2trf_rotmat,SRPdotprod)

     count_epoch_loop = start_neq  ! set this so that it prints out every 100 epochs from start_neq onwards
     normeq_tmp(:,:) = 0.d0
     AtWb_tmp (:)= 0.d0

! SA170907: put normeq, AtWb in "reduction" rather than shared (22% increase in speed). Remove OMP_CRITICAL which was a bottleneck.
!$OMP PARALLEL DO private (normeq_tmp,AtWb_tmp) shared(nobs_t,nparam_t,apr_wght,part,pre_omc,rvec,kbr_weights) reduction(+:normeq,AtWb)
! PT140901: include only the epochs between start_neq and end_neq (which can be defined in the input command file)
! PT190514: include latitude-dependent weighting on the kbrr/kbra observations so that I can downweight the Arctic obs during Nov 2016. For this, we require
!           the latitude of the satellites (can get away with using the Z coordinate, actually!)
     DO iepoch = start_neq,  end_neq, 1
        debug_iepoch = iepoch

        ! set the openMP variables to zero
        normeq_tmp(:,:) = 0.d0
        AtWb_tmp (:)= 0.d0

        ! Report status as to progress through epoch loop
! SA 200130 Preprocessor flags to remove this logging information in production mode         
#ifndef _NDEBUG
        IF(MOD(iepoch,100).EQ.0)THEN
           count_epoch_loop = count_epoch_loop + 100
           WRITE(message,'(a,i7,a,i7)')'Stack epoch loop: ',count_epoch_loop,' of',end_neq
           CALL status_update('STATUS','GRACEFIT','gracefit',' ',message,0)
        ENDIF
#endif
        !   Add Shadow condition
        CALL shadow_applyCondition_drag(iepoch,apply_shadow_cond,thrusters_off,apr_ScaleBias,part(:,:,iepoch) &
             ,pre_omc(:,iepoch), lead)

! DEBUG:
! PT190514: modify the uncertainty of the kbrr and kbra obs, subject to the location of the satellites
!        if(rvec(3,1,iepoch) > 6496700.d0)then   ! PT190514: just set it to downweight when GRACE A Z coord is > some value (ie the sat is at high latitude)
!          apr_wght(ikbrr:ikbra) = 1.d0/(10.d0**2)
!        else
!          apr_wght(ikbrr:ikbra) = kbr_weights(1:2)
!        endif
 
        !   Increment normal equation
        CALL LS_norminc(nobs_t,nparam_t,part(:,:,iepoch),pre_omc(:,iepoch),apr_wght,normeq_tmp,AtWb_tmp)
        normeq = normeq + normeq_tmp
        AtWb = AtWb + AtWb_tmp
     ENDDO  ! End of epoch loop
!$OMP END PARALLEL DO
     CALL status_update('STATUS','GRACEFIT','gracefit',' ','Finished stacking normal equations',0)

     !****************************************************************
     ! output the normal equations to a file so that I can play with them without needing to restack all the time
     !    call status_update('STATUS','GRACEFIT','gracefit',' ',"NO LONGER writing out binary normal equations",0)
     CALL output_writeNORM(rmsnam, apr_prm, normeq, AtWb, apr_wght, mcon_tides,prmnam,prmnam_size)

    call output_2d("normeq_fit.nc",  size(normeq, 1), size(normeq, 2), normeq)
    write(10002,*) AtWb

! write(*,*), ' KBR ', ikbr, ' KBRR ', ikbrr, ' KBRA ', ikbra
     ! PT191126: stop here if requested by user via the "gracefit_step" command file option.
     if(gracefit_step == 2)then
       call status_update('STATUS',calling_prog,'gracefit',' ',"Program stopped by user after outputting normal equations",0)
       stop
     endif

     !****************************************************************

     !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
     !&& !! DEBUG: write out the partials of KBRR and KBRA for a list of mascons, to compare them
     !&& !!    open(99,file="mascon_subset",status="old")
     !&& !!    im=0
     !&& !!100 read(99,*,end=110) mascs(im+1)
     !&& !!    im=im+1
     !&& !!    go to 100
     !&& !!1110 close(99)
     !&& ! PT170920: South American mascons
     !&& imp1=973
     !&& imp2=1417
     !&& !           Eurasian mascons
     !&& imp1=1967
     !&& imp2=3236
     !&& ! all mascons
     !&& imp1 = 1
     !&& imp2 = 7526

     !&& im=imp2-imp1+1
     !&& write(*,*) "logging partials for", im, " mascons"

     !&& open(97,file="mascon_partials_kbrr_allmsc",status="unknown")
     !&& open(98,file="mascon_partials_kbra_allmsc",status="unknown")
     !&& do iepoch = start_neq,end_neq
     !&&    !!      print*,iepoch,part(ikbrr:ikbra,imascons+1301,iepoch),part(ikbrr:ikbra,imascons+1364,iepoch) &
     !&&    !!        ,part(ikbrr:ikbra,imascons+1287,iepoch)
     !&&    write(97,*) iepoch,(part(ikbrr,imascons+i-1,iepoch),i=imp1,imp2)
     !&&    write(98,*) iepoch,(part(ikbra,imascons+i-1,iepoch),i=imp1,imp2)
     !&& enddo
     !&& close(97)
     !&& close(98)

     !&& !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


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
              !print*,'adding ',1.0/apr_tide_ampl_const(1,1)**2,' to M2s for mascon',i,imsctide,nparam_t
           ELSEIF(prmnam(i)(16:18) == "M2c")THEN
              normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(1,2)**2
              !print*,'adding ',1.0/apr_tide_ampl_const(1,2)**2,' to M2c for mascon',i,imsctide,nparam_t
           ELSEIF(prmnam(i)(16:18) == "O1s")THEN
              normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(2,1)**2
              !print*,'adding ',1.0/apr_tide_ampl_const(2,1)**2,' to O1s for mascon',i,imsctide,nparam_t
           ELSEIF(prmnam(i)(16:18) == "O1c")THEN
              normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(2,2)**2
              !print*,'adding ',1.0/apr_tide_ampl_const(2,2)**2,' to O1c for mascon',i,imsctide,nparam_t
           ELSEIF(prmnam(i)(16:18) == "S2s")THEN
              normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(3,1)**2
              !print*,'adding ',1.0/apr_tide_ampl_const(3,1)**2,' to S2s for mascon',i,imsctide,nparam_t
           ELSEIF(prmnam(i)(16:18) == "S2c")THEN
              normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(3,2)**2
              !print*,'adding ',1.0/apr_tide_ampl_const(3,2)**2,' to S2c for mascon',i,imsctide,nparam_t
           ELSEIF(prmnam(i)(16:18) == "K1s")THEN
              normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(4,1)**2
              !print*,'adding ',1.0/apr_tide_ampl_const(4,1)**2,' to K1s for mascon',i,imsctide,nparam_t
           ELSEIF(prmnam(i)(16:18) == "K1c")THEN
              normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(4,2)**2
              !print*,'adding ',1.0/apr_tide_ampl_const(4,2)**2,' to K1c for mascon',i,imsctide,nparam_t
           ELSEIF(prmnam(i)(16:18) == "K2s")THEN
              normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(5,1)**2
              !print*,'adding ',1.0/apr_tide_ampl_const(5,1)**2,' to K2s for mascon',i,imsctide,nparam_t
           ELSEIF(prmnam(i)(16:18) == "K2c")THEN
              normeq(i,i) = normeq(i,i) + 1.0/apr_tide_ampl_const(5,2)**2
              !print*,'adding ',1.0/apr_tide_ampl_const(5,2)**2,' to K2c for mascon',i,imsctide,nparam_t
           ENDIF
        ENDDO
     ENDIF


     !**********************************************************************
     !**********************************************************************
     ! Add mascon length-dependent constraint if necessary
     !**********************************************************************
     IF(nmascons_t > 0 .and. use_msc_len_const == 1)THEN
        CALL status_update('STATUS','GRACEFIT','gracefit',msc_const_file &
             ,"    Adding length-dependent constraints to mascon parameters",0)

        ! PT170307: read the second line in the regularisation file, which contains the length scales and constraint levels
        READ(243,'(a)')message
        ! PT150721: check that it really was a header line
        IF(message(1:7) == "# scale" .OR. message(1:5) == "# lat")THEN
           CALL status_update('STATUS','GRACEFIT','gracefit',' ',message,0)
        ELSE
           WRITE(message,'(a)')'No (or wrong) header info found in mascon spatial constraint file '
           CALL status_update('WARNING','GRACEFIT','gracefit',msc_const_file(1:trimlen(msc_const_file)),message,0)
           BACKSPACE(243)
        ENDIF
        ! read the constraint matrix
        DO i=1,nmascons_t
           READ(243,*,iostat=ioerr)(msc_const(i,j),j=1,nmascons_t)
           IF(ioerr /= 0)THEN
              WRITE(message,'(a,i6,a)')"Error reading line",i," of the mascon regularisation file"
              CALL status_update('FATAL','GRACEFIT','gracefit',msc_const_file,message,0)
           ENDIF
        ENDDO
        DO i=imascons,imascons+nmascons_t-1
           DO j=imascons,imascons+nmascons_t-1
              normeq(i,j) = normeq(i,j) + msc_const(i-imascons+1,j-imascons+1) ! *0.7 ! / 0.04 ! PT150618: scale it by 0.2 m sigma
           ENDDO
        ENDDO
     ELSE
        CALL status_update('STATUS','GRACEFIT','gracefit',' ',"NOT adding length-dependent constraints to mascon parameters",0)
     ENDIF

     !**********************************************************************
     !**********************************************************************
     ! Add mascon mass conserving constraint if necessary           PT150721
     !**********************************************************************
     IF(use_msc_conserve_mass == 1)THEN
        CALL status_update('STATUS','GRACEFIT','gracefit',' ',"    Adding mass conservation constraints to mascon parameters",0)

        ! the constraint equation is simply: m1+m2+m3+ ..... +m_n = 0
        ! therefore, the partials d/dm? = 1 for each mascon and the terms in the normal equations are just the variance of the observation
        DO i=imascons,imascons+nmascons_t-1
           DO j=imascons,imascons+nmascons_t-1
              normeq(i,j) = normeq(i,j) + 1.0/0.01**2  ! hardwire 1 cm as the sigma for this conditional equation. @#&& fix this later !
           ENDDO
        ENDDO
        ! RHS for each mascon parameter becomes the negative of the sum of all a priori mascons times the variance of the observation
        !      msc_sum_apr = sum( (imascons:imascons+nmascons_t-1))
        !      do i=imascons,imascons+nmascons_t-1
        !        AtWb(i) = AtWb(i) - msc_sum_apr * 0.005**2  ! &#$^% need to remove the hardwired value of 1 cm
        !      enddo
        !      enddo
     ELSE
        CALL status_update('STATUS','GRACEFIT','gracefit',' ',"NOT adding mass conservation constraints to mascon parameters",0)
     ENDIF
     !**********************************************************************

     !**********************************************************************
     !   Solve the normal equations for the estimated ICs and parameters
     IF(INT(use_kbrr_decomp(1)) == 0  .OR. use_kbrr_decomp(1) == 2.OR. use_kbrr_decomp(1) == 3 .or. use_kbrr_decomp(1) == 5)THEN  ! solve for all the observations and parameters
! PT/SA190820: replace with LaPack routines that do a Cholesky decomposition (faster by 30% or so)
        !CALL LS_normSolve(nparam_t,normeq,AtWb,adjust)
        CALL Chol_normSolve(nparam_t,normeq,AtWb,adjust)
     ELSE IF(INT(use_kbrr_decomp(1)) == 1)THEN  ! solve for only the mascons
        !CALL LS_normSolve(nmascons_t,normeq(imascons:imascons+nmascons_t-1,imascons:imascons+nmascons_t-1) &
        CALL Chol_normSolve(nmascons_t,normeq(imascons:imascons+nmascons_t-1,imascons:imascons+nmascons_t-1) &
             ,AtWb(imascons:imascons+nmascons_t-1),adjust(imascons:imascons+nmascons_t-1))
     ENDIF

     !******************************************************************
  ENDIF   ! end of if statement as to whether to perform the LS inversion or read a solution from a VCV file

  !********************** WRITE OUT SOLUTIONS *********************
  CALL output_writeVCV(gt_fnam,apr_prm,kbrr_misfit_num,adjust,normeq,prmnam,prmnam_size)
  CALL output_writeFIT(gt_fnam,apr_prm,apr_wght,adjust,kbrr_misfit_num,part,pre_omc,normeq,post_omc,prmnam,prmnam_size &
                       ,msc_const_file)  ! postfit residuals are computed in this subroutine
  !    call output_writeRMS(gt_fnam,post_omc,rvec(1:6,:,:))  ! Hard coded constant 6 for pos and vel (x y z). Assumed to always be in GTORB
  !    call output_writeSVS(apr_prm,adjust)
  !  print*,"**#$#%#%$ not outputting correlation file   *#&%#&% "
!      call output_writeCOR(normeq)
  !    call output_writeJPL(post_omc)
  !****************************************************************

  !   Plot the residuals
  IF( iop > 0 )THEN
     IF(nkbrr_t /= 0) CALL plot_writeKB(iop,pre_omc,post_omc,rvec(1:3,:,:),rpy,kbrr_prefit_tol,nepochs_t)
     IF(nkbra_t /= 0) CALL plot_writeKBA(iop,pre_omc,post_omc,rvec(1:3,:,:),rpy,kbrr_prefit_tol,nepochs_t)
     CALL plot_writePLT(iop,post_omc,rvec(1:6,:,:),gvec(4,:,:))
  ENDIF
  !****************************************************************

  ! ******************!  Output the mean beta angle *************************
! PT/HMcQ190301: divide by the number of epochs used, not the number of GPS epochs in the GNV1B file
!  WRITE(message,'(a,f8.2,a,f8.2)')"Mean beta angle = ",SUM(beta_angle(:))/DBLE(nepochs_t)*180.d0/pi," degrees."
  WRITE(message,'(a,f8.2,a,f8.2)')"Mean beta angle = ",SUM(beta_angle(:))/DBLE(nbeta/nsat_t)*180.d0/pi," degrees."
  CALL status_update('STATUS','GRACEFIT','gracefit/gracefit',' ',message,0)
  !****************************************************************


  !*********************** DEALLOCATE ARRAYS **********************

  CALL status_update('STATUS','GRACEFIT','gracefit',' ','Deallocating arrays',0)
  DEALLOCATE(apr_prm)
  DEALLOCATE(apr_const)
  DEALLOCATE(apr_wght)
  DEALLOCATE(part)
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
  CALL status_update('STATUS','GRACEFIT','gracefit',' ', 'Normal end of GRACEFIT',0)

END PROGRAM GRACEFIT

!********************************************************************************************************************************
