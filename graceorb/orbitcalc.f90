    subroutine orbitcalc(tin, efic, int_step, step)

! E-K Potter May 2010
!   This program will (for each arc)					
!       1. feed the earth-fixed satellite position vector into gravcalc 
!          subroutine to calculate acceleration due to gravity for 
!          each time step.
!       2. rotate the acceleration vector into inertial space
!       3. integrate the acceleration to determine new position vector 
!          for next time step.

!	jd	----  	julian day
!       t	----  	seconds of day (GPST)
!       tdtgpst	----  	TDT - GPST  (seconds)
!       iut1pol	----  	Binary coded control for diurnal pole (bit 1) and UT1 (bit 2)
!	idir 	----	flag to show if tranformation is ef to intertial (1) or vice versa (-1)
!	xpole 	----	N/A
!	ypole 	--	N/A
!	ingrav	----	N/A
!	efgrav 	----	N/A
!	efpos 	----	N/A
!	efcoor 	----	N/A
!	incoor 	----	N/A
!	efic 	----	6 position and velocity ICs in earth fixed coordinates
!	tin 	----	integration start time
!	efxin.. ---- 	N/A 
!     int_step  ----	integration step size (seconds)
!	step 	----	number of integration steps
!	rot_t 	----	transpose of matrix "rot" (declared in rotation_mod)
!	angvelten --	angular velocity tensor (used for general relativity calcs)
!	angvel   --	angular velocity
!       Other variables declared in modules listed below

! MODIFIED: APP 121004
! To allow user-defined names for the output file and the UT1, pole, and nutation tables
! PT131104: write the GTORB output file as a binary file

    use sat_mod
!    use rotation_mod
    use coeff_mod
    use spherhar_mod
    use gm_mod
    use evcoef_mod
    use timxtr_mod
    use gauleg_mod
    use lovenum_mod 
    use inmod_mod
    use accred_mod
    use gpsant_mod
    use accel_mod
    use dealias_mod
    use bsscl_mod
    use usno_mod     ! provides the max and actual number of eop values from file usno.final.data
    use gtorb_mod    ! PT170609: this declares the record length information for the GTORB file

    implicit none

    real(kind=8), dimension(10) :: efic           ! position, velocity and once- and twice-per-rev IC values
    real(kind=8)               :: tin            ! initial time (in seconds of day) of start of orbit
    integer*4 :: jd,jds,ioerr,i,j,iscrn
    integer*4 :: idir 
    integer*4 :: step
!    integer*4 :: num_mcon
    integer*4 :: date(5), day_of_year
    integer*4 mjd, PEPjd
    integer*4 :: dateend(5)
    integer :: time_array(8)
! PT170609: now declared in mod_subs/gtorb_mod
!    integer*4 :: nrec                            ! record length for binary file
    real(kind=8) ::  xpole, ypole
    real(kind=8), dimension(3) :: ingrav, efgrav, efpos
    real(kind=8) :: t  
    real(kind=8), dimension(6) :: efcoor, incoor
    real(kind=8) :: int_step
    real(kind=8), dimension(3,3) :: rot_t, angvelten
    real*8 jdate, sectag, grace_sec, grace_start_mjd
    real*8 mjdate, PEPt
    real*8 jdateend, mjdend
    real*8 :: tinend
    real*8 :: sectagend
    character*55 :: ef_filename
    character*10 :: agency, uname 
    data grace_start_mjd / 51544.5d0 /
    real(kind=8), dimension(3,3) :: SCArot_mati, rotSRF2e_init, SCArot_mat
    real(kind=8), dimension(3) :: GPSant_cos, GPSant, GPSant_iner, GPSant_efix
    real(kind=8) :: GPSant_mag
    real(kind=8), dimension(0:3) :: quatSRF2e_init, quatS2i, quati2e, quatS2e

! variables to store the inertial-efixed rotation and rotation rate information
  integer*4                 :: nepochs                  ! number of epochs of celestial-to-terrestrial rotations computed
  real(kind=8),allocatable  :: rot_mats(:,:,:)          ! array of inertial->efixed rotation matrices (one per epoch)                  
  real(kind=8),allocatable  :: rotdots(:,:,:)           ! array of time derivatives of rotation matrices
  real(kind=8),allocatable  :: rotaccs(:,:,:)           ! array of time derivatives of rotation matrices
  real(kind=8),allocatable  :: rpy(:,:)                 ! roll/pitch/yaw angles for the rotation matrices
  real(kind=8),allocatable  :: rpydots(:,:)             ! roll/pitch/yaw angle rates for the rotation matrices
  real(kind=8),allocatable  :: rot_dates(:)             ! epochs for which rotation angles and rates have been computed
  real(kind=8)              :: mjd_start                ! floating point start epoch of the orbit integration
  real(kind=8)              :: angvel(3),angvelrate(3)  ! angular velocity vector (roll, pitch, yaw) and its rate (rolldot, pitchdot, yawdot)
  real(kind=8) :: dot


! rotation and rotation rate matrices
  real(kind=8),dimension(3,3) :: rot_i2e,rotdot_i2e,rotacc_i2e,rot_e2i,rotdot_e2i,rotacc_e2i,rot_e2i_t
  real*8  tmpmat(3)

! variables for use in the IERS SOFA routines
    real(kind=8) :: Rnut_prec_framebias(3,3),casr
! interpolated xp, yp, ut1-utc values
    real(kind=8) :: xp_int,yp_int,ut1utc_int

! DEBUG variables
    logical :: debug


    casr = pi/180.d0/3600.d0

    call Date_and_Time(values=time_array)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call gmset ! to set constants for various earth parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   CALL subroutine to read in coefficents for evaluating spherical 
!   harmonics of static mean field
    call coeff_read
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   CALL subroutine to read in coefficients for Desai model for ocean pole tide
    call polcoeff_read
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   read in coefficients for everett interpolation, used by solred/lunred
    call evrtcf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   call gauleg to calculate the abscissa and weightings for use in the 
!   calpnm subroutine, which evaluates the legendre functions at each latitude
!   - these are then used to calculate the sph harm coefficients to represent a grid of data 
!   store values of associated Legendre functions at each quadrature point
!   Only half need to be stored since Pnm(-x)= (-1)^(n+m) Pnm(x)
    call gauleg(-1.d0, 1.d0, absc, wei, nlat)
    do i = nlatd2+1, nlat
      call calpnm(absc(i), pnmx(1, i-nlatd2)) !for each abscissa (latitude) calculate the legendre function array and store
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   unit numbers for table and screen
    iscrn = 40


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   set julian day for initial epoch (tin is in seconds past noon Jan 1 2000
!   and 2451545 is julian day of Jan 1 2000 12:00 hrs)
    jdate = tin/86400.d0+2451545.d0   ! actual julian *date* (non-integer)
    mjdate = jdate-2400000.5d0        ! actual modified julian *date* (non-integer)

    jd = aint(jdate)                  ! julian day (integer)
    t = (jdate-aint(jdate))*86400.d0  ! seconds of day for jdate
    mjd = aint(mjdate)                ! modified julian day (integer)

    tinend = tin+int_step*dble(step)
    mjdend = grace_start_mjd + tinend/86400.d0
    jdateend  = mjdend + 2400000.5d0
    call jd_to_ymdhmsnew(jdate, date, sectag)
    call jd_to_ymdhmsnew(jdateend, dateend, sectagend)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! allocate the variables to store the i2e rotation angles and rotation rate angles
    mjd_start = mjdate
    nepochs =  1+ (mjdend - mjd_start ) * 24.d0*60.d0*12.d0 + 10.d0*12.d0  !(12x5sec epochs per minute for an extra 5 mins at each end of orbit)

    allocate(rot_mats(3,3,nepochs))
    allocate(rotdots(3,3,nepochs))
    allocate(rotaccs(3,3,nepochs))
    allocate(rpy(3,nepochs))
    allocate(rpydots(3,nepochs))
    allocate(rot_dates(nepochs))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Generate the suite of inertial-efixed rotation matrices (and their time 
!   derivatives) using the SOFA routines that implement the IERS2010 standards
    debug = .false.
    call generate_inert_efixed('GRACEORB  ',debug,mjd_start,mjdend,rot_mats,rotdots,rotaccs,rpy,rpydots,rot_dates &
                               ,nepochs,0.d0,0.d0,0.d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   open output file with name format GTORB_date_time_sat_RL 
! PT170609: as of this date the record length is defined in mod_subs/gtorb_mod
    open (unit=15, file=output_file, status='unknown',access='direct',recl=GTORB_recl)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   PT140528: generate better a priori estimates for the biases (c0x c0y c0z) if a priori
!             values were not updated from a vcv file
    if(.not. vcv_bias)then
      if(gt_acc(1) == "Y")then
        call bias_from_ACC
      else
        call status_update('FATAL','GRACEORB','orbitcalc',' ',"cannot derive default biases because we haven't read the accel data",0)
      endif
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   position of the LC phase center of the satellite (hard-wired here)
!   set magnitude and direction cosines (taken from data)
!    if (sat.eq."A") then
!      GPSant_mag = 0.4514203544369704 
!      GPSant_cos(1) = -0.0008860920781892878
!      GPSant_cos(2) =-0.0008860920781892878 
!      GPSant_cos(3) = -0.9999992148405207
!    else
!      if (sat.eq."B") then
!        GPSant_mag = 0.4517310303930869d0
!        GPSant_cos(1) = 0.001332651421967077
!        GPSant_cos(2) = 0.001669134837480358
!        GPSant_cos(3) = -0.9999977190119395
!      else
!        stop "Sat must be A or B"
!      endif
!    endif
!    print*, "GPSant_mag, GPSant_cos", GPSant_mag, GPSant_cos 

!   calculate the components of the SRF vector from COM to GPS antenna phase centre for L1 and L2
!   *********************************************************************************************
!   I am really not convinced this is correct... but I'll leave it as is for now
!   *********************************************************************************************
    do i=1,3
      GPSant_L1(i)=GPSant_magL1*GPSant_cosL1(i) 
      GPSant_L2(i)=GPSant_magL2*GPSant_cosL2(i) 
! PT130405: try a different formula ...
!      GPSant(i)=GPSant_L2(i)-1227.6d0/1575.42d0*GPSant_L1(i)
      GPSant(i)=GPSant_L1(i)-1.984d0*(GPSant_L2(i)-0.779d0*GPSant_L1(i))
    enddo

!   calculate the magnitude and direction cosines for the antenna LC phase centre for output file header
    GPSant_mag = dsqrt(GPSant(1)**2+GPSant(2)**2+GPSant(3)**2)
    do i=1,3
      GPSant_cos(i)=GPSant(i)/GPSant_mag  
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   calculate the matrix to rotate from SRF into inertial space
! PT171003: only do this if we want accelerometer data
    if(gt_acc(1) == "Y")then
      call SCArot(jd,t,SCArot_mat,quatS2i)
    else
      SCArot_mat = 0.d0
    endif
!    call matmult(SCArot_mat,GPSant,GPSant_iner,3,3,1) now combined with next rotation

! PT140317:  calculate the rotation matrix from inertial to terrestrial using new routines
    call inert_interp(mjd_start,nepochs,rot_dates,rot_mats,rotdots,rotaccs,rot_i2e,rotdot_i2e &
                      ,rot_e2i,rotdot_e2i,rotacc_e2i,rotacc_i2e) !,rpy,rpydots,angvel,angvelrate)
  
! Combine matrices together then calculate quaternion!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   combine the rotation matrix from SRF-->inertial (apply first) with inertial-->efixed (apply second)   
    call matmult(rot_i2e, SCArot_mat, rotSRF2e_init, 3, 3, 3)  

!   calculate the SRF-->efixed quaternion to pass to and print out in adam moulton integrator
!    call rotmat2quat(rotSRF2e_init, quatSRF2e_init)
    call rotation_mat2quat_3d(rotSRF2e_init, quatSRF2e_init)

!   rotate the GPSant vector into efixed
    call matmult(rotSRF2e_init, GPSant, GPSant_efix, 3, 3, 1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Calculate angular velocity required for general relativity correction
! *************************************************************************
! Need to check - shouldn't the angvel be calculated at everytime step to feed through
! to the general rel calculation
! *************************************************************************
! PT140321: we now get this from the roll/pitch/yaw angles 
!    call trnsp(rot_e2i, rot_e2i_t, 3, 3)
!    call matmult(rotdot_e2i, rot_e2i_t, angvelten, 3, 3, 3) !angvel needed for general relativity correction 
!    angvel(1) = angvelten(3,2) 
!    angvel(2) = angvelten(1,3)
!    angvel(3) = angvelten(2,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Transform efix ICs into inertial space
    do i=1,3
      incoor(i) = rot_e2i(i,1)*efic(1)+rot_e2i(i,2)*efic(2)+rot_e2i(i,3)*efic(3)
    enddo
    do i=4,6
      incoor(i)=rot_e2i(i-3,1)*efic(4)+rot_e2i(i-3,2)*efic(5)+rot_e2i(i-3,3)*efic(6)+ &
               rotdot_e2i(i-3,1)*efic(1)+rotdot_e2i(i-3,2)*efic(2)+rotdot_e2i(i-3,3)*efic(3)
    enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! PT130321: calculate the accelerometer biases for this day so that they can be written to the file.
!!           They are stored in bs(3) which is defined in a common block. Need JD and tin as arguments.
    call biasscale(jd, tin, c0x,c0y,c0z,c1x,c1y,c1z,c2x,c2y,c2z,scl )

! PT130326: output the a priori scales and biases to the screen
    write(message,'(a,10f17.6,6f8.4)')"   Satellite apriori ICs ",efic,scl,bs(1)*1.d6,bs(2)*1.d6,bs(3)*1.d6
    call status_update('STATUS','GRACEORB','orbitcalc',' ',message,0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   CALL hdwrt to write header for output file

  if(ish5) then
    call orb_create(output_file,  orb_file)
    call hdwrt(orb_file,tin,date,sectag, tinend, dateend, sectagend, int_step, step, mjd, efic, &
           incoor, time_array, "J2000", GPSant_mag, GPSant_cos, gpsantoff)
  else
    call hdwrt_binary(tin,date,sectag, tinend, dateend, sectagend, int_step, step, mjd, efic, &
               incoor, time_array, "J2000", GPSant_mag, GPSant_cos, gpsantoff) ! to write header for GTORB output file
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   CALCULATE nd (the number of parameters plus 1 for the orbit set, pos and vel)
    nd = 13 + num_mcon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   CALL adam-moulton integrator
    call status_update('STATUS','GRACEORB','orbitcalc',' ','calling adam_moulton_part',0)
    call adam_moulton_part(incoor, efic(7:8), efic(9:10), jd, t, step, int_step, quatSRF2e_init &
! PT140318: pass the rotation matrices
                          ,nepochs,rot_dates,rot_mats,rotdots,rotaccs &
! PT140319: pass the Xp, Yp values into adam_moulton_part-> gravcalc -> sphharmfield -> poltid
! PT161005: now done with usno_mod.f90
                          ) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    return
    end

