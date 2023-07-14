  program pitch_correction

! program to calculate a few different things:
!
! 1. The roll/pitch/yaw angles between the spacecraft orientations and the nominal LOS orientation
! 2. The pitch correction required to be added to the accelerometer transplant of GRACE C to GRACE D
!
! The code here repeats aspects of code embedded inside GRACEORB and GRACEFIT. I am doing it again, deliberately,
! in an attempt to find where/why our computations of roll/pitch/yaw no longer seem to be correct ....
!
! P. Tregoning
! 1 November 2019


  use usno_mod     ! provides the max and actual number of eop values from file usno.final.data
  use rotation_mod
  use accred_mod   ! provides sca_obs/sca_obs_2 array definitions

  implicit none

  character*150 :: arg,message
  integer*4     :: year, month, day            ! epoch for which computations are required
  integer*4     :: mission                     ! 0: GRACE, 1: GRACE-FO
  character*150 :: GNV1B_1, GNV1B_2            ! filenames of GNV1B input data
  character*150 :: SCA1B_1, SCA1B_2            ! filenames of SCA1B input data
  integer*4,parameter :: lugnv1=10,lugnv2=11   ! unit numbers for the two GNV1B files

! PT200717: add unit numbers for the SCA1B files
  integer*4,parameter :: lusca1=12,lusca2=13   ! unit numbers for the two SCA1B files
  integer*4           :: n_sca1,n_sca2         ! number of star camera obs for each satellite

  ! Data read in from GNV1B files
  real(kind=8)    , allocatable :: gvec(:,:,:)       ! Positions and velocities for each satellite at every epoch
  real(kind=8)    , allocatable :: gvec_epochs(:,:)  ! Epochs of positions and velocities for each satellite at every epoch
  integer*4                     :: num_gps_used(2)   ! number of GPS observations that will be used (after decimating the 5sec data)
  logical                       :: use_gvec          ! logical to indicate whether there are GNV1B obs to match the GTORB obs
  integer*4                     :: igvec             ! counter through the GNV1B records
  integer*4                     :: epoch_interval    ! as it says
  ! PT190528: pass this variable back from input_readGNV1B
  integer*4                     :: nepochs_t         ! number of GNV1B epochs read from the input GNV1B file
  real(kind=8)                  :: gnv_step(2)       ! epoch step size for each GNV1B file
  integer*4                     :: gnv_span(2)       ! span of epochs of each GNV1B file ?

! the starting epoch in terms of GRACE seconds
  real(kind=8)                  :: starting_epoch

! variables to store the inertial-efixed rotation and rotation rate information
  real(kind=8),allocatable  :: rot_mats(:,:,:)          ! array of inertial->efixed rotation matrices (one per epoch)                  
  real(kind=8),allocatable  :: rotdots(:,:,:)           ! array of time derivatives of rotation matrices
  real(kind=8),allocatable  :: rotaccs(:,:,:)           ! array of time derivatives of rotation matrices
  real(kind=8),allocatable  :: rpy(:,:,:)                 ! roll/pitch/yaw angles for the rotation matrices
  real(kind=8),allocatable  :: rpydots(:,:)             ! roll/pitch/yaw angle rates for the rotation matrices
  real(kind=8),allocatable  :: rot_dates(:)             ! epochs for which rotation angles and rates have been computed
  real(kind=8)              :: mjd_start,mjdend         ! floating point start/end epoch of the orbit integration
  real(kind=8)              :: angvel(3),angvelrate(3)  ! angular velocity vector (roll, pitch, yaw) and its rate (rolldot, pitchdot, yawdot)
  real(kind=8) :: dot
  real(kind=8) :: mjd

! variables for computing unit vectors
  real(kind=8)  :: LOS(3)                      ! vector of the line-of-sight between the two satellites
  real(kind=8)  :: LOS_mag                     ! magnitude of the LOS vector
  real(kind=8)  :: vec1(3),vec2(3)

  integer*4     :: iepoch
  integer*4     :: ioerr
  character*20  :: calling_prog                ! name of program here. Used in subroutines to report what is happening
  integer*4     :: date(5)

! quaternion veriables
  real(kind=8)  :: quat1_rpy(0:3),quat2_rpy(0:3)               ! the roll/pitch/yaw computation for each satellite
  real(kind=8)  :: quat1_actual_ef(0:3),quat2_actual_ef(0:3)   ! quaternions for SRF to actual earth-fixed orientation of satellites
  real(kind=8)  :: quat_i2e(0:3)                               ! quaternion from inertial to earth-fixed
  real(kind=8)  :: quat1_srf_los_ef(0:3),quat2_srf_los_ef(0:3) ! quaternion from srf to nominal earth-fixed LOS

! AOC variable returned (but not used) from rpy_AOC.f90
  double precision, allocatable :: AOC(:,:,:)          ! Antenna offset corrections

!! PT200723: gamit rotsnp variables
!  real*8    :: fjd,fjd_dec,t, PEPt,sec
!  integer*4 :: PEPjd
!  integer*4,parameter        :: idir=1 ! Values to be input to rotsnp
!  integer*4  :: start_grace_seconds


  real*8        :: amag3
  integer*4     :: tmp_int
  logical       :: debug


  calling_prog = "pitch_correction"
  date = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! decode runstring
  call getarg(1,arg)
  if (arg(1:1) == " ")then
    write(message,'(a)')"Runstring: pitch_correction 2019 10 02 mission "
    call status_update('FATAL','UTIL','pitch_correction',' ',message,0)
  endif
  ! epoch
  read(arg,*)year
  call getarg(2,arg)
  read(arg,*)month
  call getarg(3,arg)
  read(arg,*)day
  ! mission
  call getarg(4,arg)
  read(arg,*)mission


  ! calculate (in GRACE seconds) the starting epoch at 00UT for this day
  date(1) = year
  date(2) = month
  date(3) = day
  call ymdhms_to_gracesec(date,0.d0,tmp_int)
  starting_epoch = dble(tmp_int)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! form up the names of the GNV1B files
  if (mission == 0)then
    write(GNV1B_1,'("GNV1B_",i4,"-",i2.2,"-",i2.2,"_A_02.asc")')year,month,day
    write(GNV1B_2,'("GNV1B_",i4,"-",i2.2,"-",i2.2,"_B_02.asc")')year,month,day
    write(SCA1B_1,'("SCA1B_",i4,"-",i2.2,"-",i2.2,"_A_02.asc")')year,month,day
    write(SCA1B_2,'("SCA1B_",i4,"-",i2.2,"-",i2.2,"_B_02.asc")')year,month,day
  else
    write(GNV1B_1,'("GNV1B_",i4,"-",i2.2,"-",i2.2,"_C_04.txt")')year,month,day
    write(GNV1B_2,'("GNV1B_",i4,"-",i2.2,"-",i2.2,"_D_04.txt")')year,month,day
    write(SCA1B_1,'("SCA1B_",i4,"-",i2.2,"-",i2.2,"_C_04.txt")')year,month,day
    write(SCA1B_2,'("SCA1B_",i4,"-",i2.2,"-",i2.2,"_D_04.txt")')year,month,day
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! open the GNV1B files
  open(lugnv1,file=GNV1B_1,status='old',iostat=ioerr)
  if(ioerr /= 0)then
    call status_update('FATAL','UTIL','pitch_correction',GNV1B_1,"Error opening GNV1B file",0)
  endif
  open(lugnv2,file=GNV1B_2,status='old',iostat=ioerr)
  if(ioerr /= 0)then
    call status_update('FATAL','UTIL','pitch_correction',GNV1B_2,"Error opening GNV1B file",0)
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the GNV1B files. Use the standard subroutines for reading L1B files
  ! first satellite
  call gnv_read_hdr(lugnv1,calling_prog ,GNV1B_1,mission,num_gps_used(1),gnv_step(1),gnv_span(1) )
  call gnv_read_hdr(lugnv2,calling_prog ,GNV1B_2,mission,num_gps_used(2),gnv_step(2),gnv_span(2) )

! allocate the array of epochs for gvec (unused in this subroutine but required when calling gnv_read_data)
  if(gnv_span(2) > gnv_span(1))nepochs_t = int(gnv_span(2) / gnv_step(2))
  if(gnv_span(1) >= gnv_span(2))nepochs_t = int(gnv_span(1) / gnv_step(1))

! Data read in from GNV1B files
  allocate(gvec(6,2,nepochs_t))
  allocate( gvec_epochs(nepochs_t,2))

  epoch_interval = 5.0
  call gnv_read_data(lugnv1,calling_prog ,GNV1B_1,mission,6,2,nepochs_t,1,starting_epoch &
                     ,gnv_step(1),num_gps_used(1),1.d0,gvec_epochs(:,1),nepochs_t,gvec)
  call gnv_read_data(lugnv2,calling_prog ,GNV1B_2,mission,6,2,nepochs_t,2,starting_epoch &
                     ,gnv_step(2),num_gps_used(2),1.d0,gvec_epochs(:,2),nepochs_t,gvec)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Generate the suite of inertial-efixed rotation matrices (and their time 
!   derivatives) using the SOFA routines that implement the IERS2010 standards
    allocate(rot_mats(3,3,nepochs_t))
    allocate(rotdots(3,3,nepochs_t))
    allocate(rotaccs(3,3,nepochs_t))
    allocate(rpy(3,2,nepochs_t))
    allocate(rpydots(3,nepochs_t))
    allocate(rot_dates(nepochs_t))

    debug = .false.
    call ymdhms_to_jd(date,0.d0,mjd_start)
    mjd_start = mjd_start - 2400000.5d0
    call ymdhms_to_jd(date,dble(nepochs_t*5),mjdend)
    mjdend = mjdend - 2400000.5d0
    call generate_inert_efixed(debug,mjd_start,mjdend,rot_mats,rotdots,rotaccs,rpy,rpydots,rot_dates,nepochs_t,0.d0,0.d0,0.d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the SCA1B files. Use the standard subroutines for reading L1B files.
! PT200717: get_sca_data opens, reads the file, infills the arrays and returns
!           in sca_obs and sca_obs_2 (defined in accred_mod)
  call get_sca_data(.true.,calling_prog,1,SCA1B_1)
  call get_sca_data(.true.,calling_prog,2,SCA1B_2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! do whatever is necessary to open the gamit nutabl., pole., ut1. files etc
!  frame = "J2000"
!  iut1 = 14
!  ipole = 15
!  inut = 16
!! PT200724: code from the old gracefit
!! Open the ut1, pole, and nutation files (error FATAL because they are needed for computation of sun position)
!!   open ut1 file
!    call input_openFile(iUT1,'ut1.','old','GRACEFIT','input_openAllFiles','FATAL','Error opening ut1 table: ')
!!   open pole file
!    call input_openFile(iPOLE,'pole.','old','GRACEFIT','input_openAllFiles','FATAL','Error opening pole table: ')
!!   open nutation file
!    call input_openFile(iNUT,'nutabl.','old','GRACEFIT','input_openAllFiles','FATAL','Error opening nutation table: ')
!
!  !call gsec_to_ymdhms (start_grace_seconds, date, sec)
!  call ymdhms_to_jd(date,0.d0,fjd)
!  fjd_dec = fjd - int(fjd)                  ! decimal part of fjd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



debug = .false.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! loop through all the epochs
  do iepoch = 1  ,nepochs_t
    if( amag3(gvec(1:3,1,iepoch)) > 0.d0 .and. amag3(gvec(1:3,2,iepoch))> 0.d0  )then   ! we have GNV1B data for both satellites


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! calculate the quaternion from SRF to LOS. Since these coords are Earth-fixed, it will be SRF to Earth-fixed LOS
      vec1 = gvec(1:3,1,iepoch)
      vec2 = gvec(1:3,2,iepoch)
!print*,'vec2-vec1',vec2-vec1,(vec2-vec1)/amag3(vec2-vec1)
      call srf_2_los(.false.,vec1,vec2,quat1_srf_los_ef,1)   ! 1 for SRF to LOS, -1 for LOS to SRF
      vec1 = gvec(1:3,1,iepoch)
      vec2 = gvec(1:3,2,iepoch)
      call srf_2_los(debug,vec2,vec1,quat2_srf_los_ef,1)   ! 1 for SRF to LOS, -1 for LOS to SRF
      quat1_srf_los_ef =   1.d0*quat1_srf_los_ef  ! change the sign so that it matches the computation in graceorb
      quat2_srf_los_ef =   1.d0*quat2_srf_los_ef  ! change the sign so that it matches the computation in graceorb
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! interpolate to derive the inertial-to-ef 3x3 rotation matrix for this epoch
      mjd = mjd_start + dble(iepoch-1)*5.d0/86400.d0 +0.d0/86400.d0
      call inert_interp(mjd,nepochs_t,rot_dates,rot_mats,rotdots,rotaccs,rot_i2e,rotdot_i2e &
                      ,rot_e2i,rotdot_e2i,rotacc_e2i,rotacc_i2e)

!call printmat(rot_i2e,3,3,'rot_i2e ')
! DEBUG:
! PT200723: code in the old gamit way of doing this, using rotsnp ...
!     call PEPtime(int(fjd),fjd_dec*86400.d0+(iepoch-1)*epoch_interval, PEPjd, PEPt)
!print*,'PEPjd,PEPt',int(sca_obs(iepoch,1)),int(fjd),fjd_dec*86400.d0+(iepoch-1)*epoch_interval,PEPjd,PEPt
!    print*, "efixinert"
!      call rotsnp_grace( -1,PEPjd,PEPt,tdtgpst,iut1pol,rot,rotdot,sidtm,iut1,ipole,inut,frame,precmod ) 
! replace the IERS2010 rotation matrix with the rotsnp one
!call transp(rot,rot_i2e,3,3)
!call printmat(rot,3,3,'rot_rotsnp ')
!rot_i2e = rot

      ! convert it to a quaternion
      call rotation_mat2quat_3d ( rot_i2e, quat_i2e )

      ! add the relevant SCA-to-inertial quaternion for each satellite, thus giving us the quaternion to 
      ! transform from SRF to earth-fixed frame for each satellite
      call quat_mul(quat_i2e,sca_obs(iepoch,2:5),quat1_actual_ef)
      call quat_mul(quat_i2e,sca_obs_2(iepoch,2:5),quat2_actual_ef)
      quat1_actual_ef(1:3) =  1.d0*quat1_actual_ef(1:3)  ! change the sign so that it matches the computation in graceorb
      quat2_actual_ef(1:3) =  1.d0*quat2_actual_ef(1:3)  ! change the sign so that it matches the computation in graceorb

      ! calculate the rpy differences of the actual orientation of each satellite wrt the nominal LOS orientation
      call quat_conj(quat1_actual_ef)
      call quat_mul(quat1_actual_ef,quat1_srf_los_ef,quat1_rpy)
      call quat_conj(quat2_actual_ef)
      call quat_mul(quat2_actual_ef,quat2_srf_los_ef,quat2_rpy)

      ! and convert the rpy angles to rotations around the along-track, cross-track and radial axes

! now computed using the "new" subroutine (new as of 2013!)
      call quat_to_rpy_new(quat1_rpy,rpy(1,1,iepoch),rpy(2,1,iepoch),rpy(3,1,iepoch))
      call quat_to_rpy_new(quat2_rpy,rpy(1,2,iepoch),rpy(2,2,iepoch),rpy(3,2,iepoch))



! DEBUG
!call printmat(rot_i2e,3,3,'rot_i2e ')
!print*,'PTC quat     ',int(sca_obs(iepoch,1)),sca_obs_2(iepoch,2:5)
!print*,'PTC quati_2_e',int(sca_obs(iepoch,1)),quat_i2e
!print*,'PTC qtsrf_2_e',int(sca_obs(iepoch,1)),quat2_srf_los_ef
!print*,'**pos sat1**:',int(sca_obs(iepoch,1)), vec1,vec2
!print*,'e-fixed LOS and unit LOS',int(sca_obs(iepoch,1)),vec2-vec1,(vec2-vec1)/amag3(vec2-vec1)
!print*,'quat2_srf_los_ef    ',int(sca_obs(iepoch,1)),quat2_srf_los_ef
!print*,'PTC rpy             ',int(sca_obs(iepoch,1)),rpy(:,1,iepoch)
print*,'rpy1:',iepoch,rpy(:,1,iepoch),rpy(:,2,iepoch),rpy(2,1,iepoch)-rpy(2,1,iepoch)
!print*," "

    endif


  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  end

!****************************************************************
  subroutine input_openFile(file_unit_num,file_name,file_status,program_name,subroutine_name,error_type,error_message)

    implicit none

!********************  Variable declarations ********************
    integer*4,    intent(in) :: file_unit_num   ! File unit number to be opened
    character(*), intent(in) :: file_name       ! Name of file to be opened
    character(*), intent(in) :: file_status     ! Status of file to be opened
    character(*), intent(in) :: program_name    ! Program name
    character(*), intent(in) :: subroutine_name ! Subroutine name
    character(*), intent(in) :: error_type      ! Type of error if any occurs
    character(*), intent(in) :: error_message   ! Message in case of error

    integer*4 :: ioerr      ! Standard error variable
!****************************************************************

! Open file and report stat
    open(unit=file_unit_num,file=file_name,status=file_status,iostat=ioerr)
    if (ioerr /= 0) call status_update(error_type,program_name,subroutine_name,file_name,error_message,ioerr)
    call status_update('STATUS',program_name,subroutine_name,file_name,'Opened file: ',ioerr)
!****************************************************************

    return
  end subroutine input_openFile
!****************************************************************



!Copyright (c) Massachusetts Institute of Technology and The University of
!California at San Diego, 1994. All rights reserved.

      subroutine rotsnp_grace( idir,jd,t,tdtgpst,iut1pol,rot,rotdot,sidtm,iut1,ipole,inut,frame,precmod )
!
!     Y. Bock and R. King 1986-1994.
!
!     Calculate the rotation matrix between earth-fixed and inertial
!
!     Input:
!       idir = 1   Earth-fixed to inertial
!             -1   Inertial to Earth-fixed
!       jd         Julian day (PEP JD, beginning at midnight, = MJD + 1)
!       t          Seconds-of-day  (GPST)
!       tdtgpst    TDT - GPST  (seconds)
!       iut1pol    Binary coded control for diurnal pole (bit 1) and UT1 (bit 2)
!       iut1, ipole, inut - unit numbers for tables 
!       frame,precmod - Inertial frame and precession model to be used in rotations

!      Output:
!        rot(3,3)     Rotation matix (rad)
!        rotdot(3,3)  Derivative of rotation matrix (rad/s)
!        sidtm        Sidereal time (rad)

      implicit none
!
      character*5 frame, precmod

      integer*4 idir,jd,jds,iut1pol,iut1,ipole,inut,i,j

      real*8 t,ts,tdtutc,tdtgpst,rot,rotdot,prec,rnut,srot,sdrot &
           , snpmat,pnsmat,snpdmt,pnsdmt,eqe,sidtm,xpole,ypole   &
           , gpstutc,taiutc,sidten,prot

      dimension prec(3,3),rnut(3,3),srot(3,3),sdrot(3,3),        &
                snpmat(3,3),pnsmat(3,3),snpdmt(3,3),pnsdmt(3,3)  &
              , rot(3,3),rotdot(3,3),sidten(3,3),prot(3,3)

! Precession and nutation transformations
      call pnrot_grace(inut,jd,t,tdtgpst,eqe,prec,rnut,frame,precmod)

! Sidereal rotation and polar motion
!     srotat (unique in GAMIT) expects UTC
      jds = jd
      ts = t
      gpstutc = taiutc(jds) - 19.d0
      call timinc(jds,ts,-gpstutc)
!     catch possible leap second
      if( jds.ne.jd ) then
        gpstutc = taiutc(jds) - 19.d0
        jds = jd
        ts = t
        call timinc(jds,ts,-gpstutc)
      endif
      tdtutc = taiutc(jds) + 32.184d0
!d      print *,'ROTSNP jd t jds ts gpstutc tdtutc '
!d     .               ,jd,t,jds,ts,gpstutc,tdtutc
      call srotat_grace( jds,ts,tdtutc,eqe,iut1pol,srot,sdrot,sidtm, xpole,ypole,sidten,prot,iut1,ipole,precmod)

! Combine the matrices
      if(idir.eq.-1 ) then
        call snp(prec,rnut,srot,snpmat)
        call snp(prec,rnut,sdrot,snpdmt)

        do 100 j=1,3
        do 100 i=1,3
        rot(I,J)=snpmat(I,J)
100     rotdot(I,J)= snpdmt(I,J)

      else
        call pns(prec,rnut,srot,pnsmat)
        call pns(prec,rnut,sdrot,pnsdmt)

         do 200 J=1,3
         do 200 I=1,3
         rot(I,J)=pnsmat(I,J)
200      rotdot(I,J)=pnsdmt(I,J)

      end if

      return
      end


!!!!!!!!!!******************!!!!!!!!!!!!!!!!!!!!*****************
!Copyright (c) Massachusetts Institute of Technology and the University of
!California at San Diego, 1994. All rights reserved.

      SUBROUTINE PNROT_grace(INUT,JD,T,TDTGPST,EQE,PREC,RNUT,frame,precmod)
!C
!C Written by Yehuda Bock (1987)
!C
!C Compute the precession and nutation matrices

!c     Input:
!c       inut       Unit number for nutation table
!c       jd         Julian day (PEP JD, beginning at midnight, = MJD + 1)
!c       t          Seconds-of-day  (GPST)
!c       tdtgpst    TDT - GPST  (seconds)
!c       frame      Inertial frame
!c       precmod    precession model
!c
!c      Output:
!c        EqE          Equation of the equinoxes (rad)
!c        prec         Precession matrix (rad)
!c        rnut         Nutation matrix (rad)

      implicit none

      character*5 frame, precmod

      integer*4 jd,inut

      real*8 t,tdtgpst,eqe,prec,rnut,deps,oblq,dpsi,fjdtdt
!C
      dimension prec(3,3),rnut(3,3)

!c Convert GPS time to Terrestrial Dynamical Time (still PEP JD = true JD + 0.5)

      FJDTDT= DBLE(JD) + T/86400.D0 + TDTGPST/86400.d0

!C Form the precession matrix


!c     temporary to get 1950 to 2000 matrix
!c      fjdsave = fjdtdt
!c      fjdtdt = 2433282.923d0
!c      call PRCES(FJDTDT,OBLQ,PREC,frame,precmod)
!c      print *,'B1950 PREC: ',prec
!c      fjdtdt = fjdsave

      call PRCES(FJDTDT,OBLQ,PREC,frame,precmod)

!C Form the nutation matrix

      Call NUTTAB(inut,FJDTDT,OBLQ,DPSI,DEPS,RNUT)

!C Compute the equation of the equinoxes for computation of GAST

      EqE=dpsi*dcos(oblq)

      return
      end



!!!!!!!!!!******************!!!!!!!!!!!!!!!!!!!!*****************
      SUBROUTINE NUTTAB( INUT,FJD,OBLQ,DPSI,DEPS,RNUT )
!C
!C     NUTTAB takes JD as an argument and return the nutations
!C     in longitude and obliquity (in radians) and the nutation matrix.
!C     Values of DEPS AND DPSI are interpolated from tables.
!C     R.I. Abbot - June 1984
!C     Merged orbits/model version moved to libary by R. King - January 1992
!c     Sense of matrix (from /orbits) opposite from original /model version
!C
!C     Input
!C
!C        INUT   Unit number of nutation table
!C        FJD    Julian date
!C        OBLQ   Mean obliquity of date
!C
!C      Output
!C
!C        DPSI   Delta psi
!C        DEPS   Delta epsilon
!C        RNUT   Nutation matrix

      implicit none

      integer*4 inut

      REAL*8 RNUT(3,3),RNA(3,3),RNB(3,3),RNC(3,3),WORK(3,3), pc,ps,pi,casr,deps,dpsi,oblq,fjd,cobliq,sobliq


      PI= 4.D0*DATAN(1.D0)
      CASR= PI/180.D0/3600.D0

!C
!C        Obtain nutation parameters from table
!C        PEP JD is required as time argument to nutred
!C
      CALL NUTRED( INUT,FJD,DPSI,DEPS )

!C
      DPSI=DPSI*CASR
      DEPS=DEPS*CASR
!C
!C
!C        Auxiliary quantities for the computation of the nutation matrix

      COBLIQ=DCOS(OBLQ)
      SOBLIQ=DSIN(OBLQ)
      PC=COBLIQ*DPSI
      PS=SOBLIQ*DPSI




!c       Calculate the nutation matrix (Reference:  Mueller, p. 75)

      Call ROTMAT(OBLQ,1,RNA)
      Call ROTMAT(-DPSI,3,RNB)
      Call ROTMAT(-OBLQ-DEPS,1,RNC)
      Call MATMPY(RNB,RNA,WORK,3,3,3)
      Call MATMPY(RNC,WORK,RNUT,3,3,3)

      RETURN
      END




!!!!!!!!!!******************!!!!!!!!!!!!!!!!!!!!*****************
      SUBROUTINE NUTRED( INUT,FJD,DPSI,DEPS )
!
!C     Read nutation values from an external data set
!C     R.I. Abbot - NOVEMBER 1984
!C     MODEL version into library by R. King  January 1992
!C
      implicit none

      integer*4 inut,itsav,nvtot,it,nr,nrec,npr,nvr , i,j,k,n,len,rcpar

      REAL*8 JD1,JD2,JDT1,JDT2,fjd , TAB(4),NUTVAL(8,2),Y1(2,2),Y2(2,2), fjdf,fjdn,deps,dpsi,rint,s,t,rjd,f2

      character*80 prog_name
      character*256 message

      logical first
!C
      SAVE NUTVAL,JDT1,JDT2,JD1,JD2,Y1,Y2,ITSAV,NVTOT,NVR,NREC,FIRST

      DATA NPR,RINT/4,0.5D0/,first/.true./

      FJDN = DINT(FJD)
      FJDF = FJD - FJDN
!C
      IF (.not.first) GO TO 130
!C
!C        Read the header records
!C
      REWIND INUT
      READ (INUT,50,end=520) JDT1,JDT2
   50 FORMAT (/,35X,F7.0,1X,F7.0)
      first = .false.
      NVTOT = 0
      NREC=0
      JDT2 = JDT2 - RINT
      ITSAV=-9999

!C
!C
!C        Read the data record into storage
!C!

  100 IF (NVTOT.GT.8-NPR) GO TO 120
!C        Keep reading until (8/NPR) records are in storage
      READ(INUT,110,end=530 )RJD,((NUTVAL(NVTOT+I,K),K=1,2),I=1,NPR),NVR
  110 FORMAT (1X,F5.0,1X,8(F7.0,1X),7X,I2)
      NREC = NREC + 1
!C        JD1 and JD2 are the limits of usable values in storage
!C        and depend on the interpolation scheme
      IF (NVTOT.EQ.0) JD1 = RJD + 2400000.D0 + RINT
      IF (NVR.EQ.0) NVR = NPR
      NVTOT = NVTOT + NVR
      ITSAV=-9999
      JD2 = RJD + 2400000.D0 + (NVR-2)*RINT
      IF (JD2.LT.JDT2) GO TO 100
!C        Save the first data just in case the header JDT1 is wrong
  120 IF ( first ) JDT1 = JD1
!C
!C
!C        Is JD within range of the table ?
!C
  130 IF (FJDN.LT.JDT1.OR.FJDN.GT.JDT2) GO TO 500
!C
!C        Is JD too close to JD1?
!C*      IF (FJDN.GE.JD1) GO TO 170
      IF (FJD.GE.JD1) GO TO 170
!C        If so, backspace and try again
!C**      N = (JD1-FJDN)/RINT/NPR
!c             round up any fractional part - rwk 92/2/27
      N = NINT ( (JD1-FJD)/RINT/NPR + .5d0 ) + 2
      IF (N.GT.NREC) GOTO 540
      DO 160 I=1,N
  160 BACKSPACE INUT
      NREC = NREC - N
      NVTOT = 0
      GO TO 100
!C
!C        Is JD too close to JD2 ?
!C*  170 IF (FJDN.LT.JD2) GO TO 200
  170 IF (FJD.LT.JD2) GO TO 200
!C        If so, shift storage and read another record
      NVTOT = NVTOT - NVR
      DO 180 I=1,NVTOT
      DO 180 K=1,2
  180 NUTVAL(I,K) = NUTVAL(NPR+I,K)
      JD1 = JD1 + NPR*RINT
      GO TO 100
!C
!C
!C        Calculate interpolation times and value of tabular points
!C!
  200 T = FJDN - JD1
      T = (T + FJDF)/RINT
      IT = T
      T = DMOD(T,1.D0)
      S = 1.D0 - T
      IF (IT.EQ.ITSAV) GO TO 230
      DO 225 K=1,2
      DO 210 I=1,4
      J = IT + I
  210 TAB(I) = NUTVAL(J,K)
!C
!C
!C        Calculate interpolation Y-vector
!C
      DO 220 I=1,2
      NR = I + 1
      F2 =     0.166666666666667D0 * (TAB(NR+1)+TAB(NR-1))
      Y1(I,K) =  1.333333333333333D0 * TAB(NR) - F2
  220 Y2(I,K) = -0.333333333333333D0 * TAB(NR) + F2
  225 CONTINUE
      ITSAV = IT
!C
!C
!C        Second difference interpolation
!C
  230 DPSI=(T*(Y1(2,1)+T*T*Y2(2,1))+S*(Y1(1,1)+S*S*Y2(1,1)))*1.D-4
      DEPS=(T*(Y1(2,2)+T*T*Y2(2,2))+S*(Y1(1,2)+S*S*Y2(1,2)))*1.D-4
      GOTO 999
!C
!C
!C        Out of range of table
!C
  500 WRITE(message,510) FJDN,JDT1,JDT2  
  510 FORMAT('JD= ',F8.0,' out of range of Nutation Table, JDT1= ',F8.0,' JDT2= ',F8.0)
!c     get calling module name for status_update
      len =  rcpar(0,prog_name)
      call status_update('FATAL',prog_name,'lib/nutred',' ',message,0)
  520 call status_update('FATAL',prog_name,'lib/nutred',' ','EOF on header',0)
  530 call status_update('FATAL',prog_name,'lib/nutred',' ','EOF in nutation table',0)
  540 call status_update('FATAL',prog_name,'lib/nutred',' ' ,'Backspacing ahead of start of nuttbl.',0)

  999 RETURN
      END


!*********************************************************************************
! Copyright (c) Massachusetts Institute of Technology and University of
! California at San Diego, 1994. All rights reserved.

      subroutine srotat_grace( jd,t,tdtutc,eqe,iut1pol,srot,sdrot,sidtm, xpole,ypole,sidten,prot,iut1,ipole,precmod )
! C
! C      Compute S rotation matrix from CRS true-of-date to CTRS

! c     Note that this routine (uniquely in GAMIT) uses UTC not GPST as
! c     a calling argument since the former is more natural for the Earth's
! c     rotation (rwk 940714)

! c     Input:
! c       jd         Julian day (PEP JD, beginning at midnight, = MJD + 1)
! c       t          Seconds-of-day
! c       tdtutc     TDT - UTC  (seconds)
! c       eqe        Equation of equinoxes
! c       iut1pol    Binary coded control for diurnal pole (bit 1) and UT1 (bit 2)
! c                   bit 3 (8) is Ray model; bit 4 (16) is libration terms for UT1
! c                  (iut1pol is isptide in 'model' and 'solve' 
! c       iut1, ipole, inut - unit numbers for tables
! c       precmod - precession model to use for rotations

! c      Output:
! c        srot(3,3)    Rotation matix - wobble and sideral (rad)
! c        sdrot(3,3)   Derivative of rotation matrix (rad/s)
! c        sidtm        Sidereal time (rad)
! c        xpole, ypole Pole position (rad)
! c        sidten(3,3)  Sideral time matrix |  These needed in MODEL for
! c        prot(3,3)    Wobble matrix       |  EOP partials

      implicit none

      logical bitmap

      character*5 precmod
      character*80 prog_name

      integer*4 jd,iuttyp,iut1,ipole,iut1pol,len,rcpar
! c      integer*4 i,j

      real*8 t,taiut1,taiutc,ut1utc,ut1dot,fract,tjd,xl,f,d,ascm &
           , xrot,yrot,prot,srot,sdrot,sidten,sdtmat             &
           , xpole,xpdot,ypole,ypdot,eqe,pi,casr,sidtm           &
           , tdtutc,fjdutc,rmjd,fjdtdt,dx,dy,dUT1,dUTlib,dLODlib 

! c MOD TAH 110207: Implemented IERS model
      real*8 cor_x, cor_y, cor_lod  ! returns from IERS routine
              ! pmut1_oceans and pm_gravi.  Units arc-sec and secs.

      dimension xrot(3,3),yrot(3,3),prot(3,3),srot(3,3),sdrot(3,3)
      dimension sidten(3,3),sdtmat(3,3)
                        
! c      get calling program name and m-file name for status_update
       len = rcpar(0,prog_name)

      fract = t/86400.d0

! c compute diurnal and seimdiurnal contributions to the pole position
! c and ut1 using VLBI derived values for tidal variations. dx, dy, are
! c output in milliarc seconds, and dut1 in millitime seconds.

! c **  this recalculation should be unncessary since tdtutc is passed in
! c **  add the following as a trap until we've verified this in all modules:
! c **  (tdt-utc was 52.184 in 1981 and has increased to 61.184 by 1994)
      if( dabs(tdtutc-56.d0).ge.20.d0 ) &
        call status_update('WARNING',prog_name,'lib/srotat',' ' &
                        , 'Input tdtutc undefined ',0)
      tdtutc = 32.184d0 + taiutc(jd)
! c     JD here is PEP JD (= true JD + 0.5  and MJD -1 + 2400000 )
      fjdtdt= dble(jd) + fract + tdtutc/86400.d0
      fjdutc= dble(jd) + fract 
      rmjd = fjdutc - 2400001.d0 
! cd     print *, 'jd fract fjdutc rmjd ',jd,fract,fjdutc,rmjd 
! c      if ((bitmap(iut1pol,1)).or.(kbit(iut1pol,2))) then
! c       call sd_comp(fjdutc-0.5d0,dx,dy,dut1) 
! c       write(*,'(a,f12.3,3f7.3)') 'SD_COMP: ',fjdutc,dx,dy,dut1
! c **    force use of new short-period ut1/pole terms -- tah/rwk 990320
! c         call ray( fjdutc-0.5d0,dx,dy,dut1 ) 
! c       write(*,'(a,f12.3,3f7.3)') 'RAY   P: ',fjdutc,dx,dy,dut1
! c      endif

!*      See if Ray or IERS model (bit 4 sets IERS)
! c      The Gipson model is implementbed (bit 6 sets Gipson model) -- lei 151007
       if( bitmap(iut1pol,4) ) then   ! IERS model
           call PMUT1_OCEANS (fjdutc-0.5d0,dX, dY, dUT1,cor_lod)             
!*          Solid Earth terms          
           call PM_GRAVI (fjdutc-0.5d0,cor_x,cor_y)
           dX = (dX + cor_x)*1000  ! Convert to mas
           dY = (dY + cor_y)*1000  ! Convert to mas
           dUT1 = dUT1*1000        ! Convert to ms
       elseif( bitmap(iut1pol,6) ) then   ! Gipson model 
           call gipson(fjdutc-0.5d0,dX, dY, dUT1)
           dX = dX*1000  ! Convert to mas
           dY = dY*1000  ! Convert to mas
           dUT1 = dUT1*1000        ! Convert to ms
       else 
           call ray( fjdutc-0.5d0,dx,dy,dut1 )  ! Ray Model
       end if
       

! c Read UT1 from the input file
      call ut1red_grace( iut1,jd,fract,taiut1,ut1dot,iuttyp )
      if( taiut1.ge.0.d0 ) ut1utc=taiutc(jd)-taiut1

! C Compute fundamental arguments for tidal effects on UT1
      tjd= jd + fract - 0.5D0
      call funarg( tjd,xl,f,d,ascm )
! C     Add tidal correction if UT1-UTC is regularized (UT1R)
      if( iuttyp.eq.2 ) then
        call ut1tid( ut1utc,xl,f,d,ascm )
      else if (iuttyp.eq.0 ) then
        call status_update('FATAL',prog_name,'lib/srotat',' ','UT1 type = 0, set = 2 or 4 in UT1. table',0)
      endif

! c Add in diurnal and semidiurnal tide contribution to UT1 if requested.
      if (bitmap(iut1pol,2)) then
        ut1utc = ut1utc + (dUT1/1000.d0)
      endif
! cd       print *,'Aft SE tides ut1utc ',ut1utc
! cd       print *,'iut1pol ',iut1pol

! c Add in semidiurnal libration contribution to UT1 if requested
! cd      print *,'bef lib ut1utc ',ut1utc
       if( bitmap(iut1pol,5) ) then 
         call UTLIBR( rmjd,dUTlib,dLODlib )
	 ut1utc = ut1utc + (dUTlib/1.d6) 
! cd       print *,'Aft lib iut1pol dUTlib ut1utc ',
! cd     .      iut1pol,dUTlib,ut1utc
       endif
! cd       stop

! c** Test of UTLIBR:
! cd      rmjd = 44239.1d0
! cd       call UTLIBR( rmjd,dUTlib,dLODlib )
! cd       print *,'TEST rmjd dUTlib,dLODlib ',rmjd,dutlib,dlodlib
! cd       rmjd = 55227.4d0
! cd       call UTLIBR( rmjd,dUTlib,dLODlib )
! cd       print *,'TEST rmjd dUTlib,dLODlib ',rmjd,dutlib,dlodlib
! cd       STOP
 
! c Read pole position from table
      call polred( ipole,jd,fract,xpole,ypole,xpdot,ypdot )
! cd    print *,'jd fract xpole ypole ',jd,fract,xpole,ypole
      

! c Add in diurnal and semi diurnal contributions to pole position if requested.
      if (bitmap(iut1pol,1)) then
        xpole = xpole + (dx/1000.d0)
        ypole = ypole + (dy/1000.d0)
      endif

! C     write(*,800) fjdutc-0.5d0, xpole, ypole, ut1utc, 
! C    .             dx, dy, dut1, iut1pol
! C800  format('EOP ',F14.6,1x,2F11.6,1x,F11.8,1x,2F8.3,1x,F8.5,1x,o4)

! c Convert pole position to radians
      pi= 4.D0*datan(1.D0)
      casr= pi/180.d0/3600.d0
      xpole = xpole*casr
      ypole = ypole*casr

! c Compute GAST and sidereal rotation matrix
      call sidmat(jd,fract,ut1utc,eqe,sidten,sdtmat,sidtm,precmod)
! c     print *,'for SIDMAT iut1pol ut1utc ',iut1pol,ut1utc
! c
! c      print*,'in SROTAT jd,t,tdtutc,eqe,iut1pol,srot,sdrot,sidtm '
! c     .                ,jd,t,tdtutc,eqe,iut1pol,srot,sdrot,sidtm
! c
! c      print*,'          xpole,ypole,sidten,prot,iut1,ipole'
! c     .                ,xpole,ypole,sidten,prot,iut1,ipole

! c ******** debug *************
! c      print*,' in srotat sidtm = ',sidtm*86400.d0/(2*pi)
! c      stop
! c ****************************

! c      write(*,700) ((SIDTEN(i,j),j=1,3),i=1,3)
! c 700  Format(' GAST rotation matrix :',/,3(1x,3D22.14,/))
! c      write(*,701) ((SDTMAT(i,j),j=1,3),i=1,3)
! c 701  format(' GAST dot rotation matrix :',/,3(1x,3D22.14,/))

! C Compute polar motion xp rotation
      Call rotmat(-xpole,2,xrot)
! C Compute polar motion yp rotation
      Call rotmat(-ypole,1,yrot)
! C Compute total polar motion rotation
      Call matmpy(xrot,yrot,prot,3,3,3)

! c      write(*,705) ((XROT(i,j),j=1,3),i=1,3)
! c 705  format(' X-Pole rotation matrix :',/,3(1x,3D22.14,/))
! c      write(*,706) ((YROT(i,j),j=1,3),i=1,3)
! c 706  format(' Y-Pole rotation matrix :',/,3(1x,3D22.14,/))
! c      write(*,702) ((PROT(i,j),j=1,3),i=1,3)
! c 702  format(' Pole rotation matrix :',/,3(1x,3D22.14,/))

! c Compute S-rotation matrix
      Call matmpy(prot,sidten,srot,3,3,3)
! c Compute SDOT-rotation matrix
      Call matmpy(prot,sdtmat,sdrot,3,3,3)
! c      write(*,703) ((srot(i,j),j=1,3),i=1,3)
! c 703  format(' S-rotation matrix :',/,3(1x,3D22.14,/))
! c      write(*,704) ((sdrot(i,j),j=1,3),i=1,3)
! c 704  format(' SDOT-rotation matrix :',/,3(1x,3D22.14,/))

      return
      end
!*********************************************************************************


!*********************************************************************************
      SUBROUTINE UT1RED_grace( IUT1,JD,FRACT,TAIUT1,UT1DOT,IUTTYP )
! C
! C     Read UT1 values from an external data set
! C     Output variable is TAI-UT1
! C     R. Abbot - Nov 1984 from PEP routines
! C     R. King -  Jun 1987
! C
      implicit none
! C
      character*32   afmt
      INTEGER*4 JD,MJD,JD1,JD2,JDT1,JDT2,J00001,iutval(12),itsav &
              , iut1,int,nvtot,nr,nrec,nrecs,ifirst,npr,nvr &
              , ioerr,itype,iuttyp,it,nback,nrbmin,lap,i,j,n

      REAL*8 TAB(4),UTVAL(12),Y1(2),Y2(2),taiut1,ut1dot,s,t,f2,fract, rmult

      integer*4 len,rcpar

      character*80 prog_name
      character*256 message


      SAVE UTVAL,JDT1,JDT2,JD1,JD2,Y1,Y2,ITSAV,NVTOT,NVR,NREC,IFIRST,NPR,INT,rmult,AFMT,itype

      DATA J00001/2400001/
      nrecs=0
      nback=0
      nrbmin=0
      lap=0

! c     get calling program name and m-file name for status_update
      len = rcpar(0,prog_name)

! c     handle problem when JD skips backwards, as in multiple ARC integrations
      if( jd.lt.jd1 ) ifirst = 0

      IF (IFIRST.NE.0) GO TO 130
! C
! C     Read the Header Records, assuming format goodies
! c     on the second line.
! C
      REWIND IUT1
! c     READ (IUT1,50) afmt,JDT1,JDT2,NPR,INT
      READ (IUT1,50) afmt,itype,JDT1,JDT2,NPR,INT,rmult
! c  50 FORMAT (/,a35,I7,1X,I7,2X,I1,2X,I1)
   50 FORMAT (/,a32,1x,i1,1x,i7,1X,I7,2X,I1,2X,I1,10x,g7.0)

! c     find the end of the format format
      i = index(afmt,' ')
      afmt = afmt(1:i-1)
! c     choke if you can't find it
      if (i.lt.1) goto 540

      IFIRST = 1
      NVTOT = 0
      NREC=0
      nrbmin= 12/npr + 1
      JDT2 = JDT2 - INT
      ITSAV=-9999

! C     Read data record into storage
! C     Keep reading until (12/NPR) records are in storage
! c     read them as integers, convert them to reals

  100 IF (NVTOT.GT.12-NPR) GO TO 120
      READ (IUT1,fmt=afmt,ERR=520,iostat=ioerr)MJD,(IUTVAL(I),I=1,NPR),NVR
      do 105 i = 1,npr
         utval(nvtot+i) = dble(iutval(i))
  105 continue
      NREC = NREC + 1

! C     JD1 and JD2 are the limits of usable values in storage
! C     and depend on the interpolation scheme
      IF (NVTOT.EQ.0) JD1 = MJD + J00001 + INT
      IF (NVR.EQ.0) NVR = NPR
      NVTOT = NVTOT + NVR
      ITSAV=-9999
      JD2 = MJD + J00001 + (NVR-2)*INT

      IF (JD2.LT.JDT2) GO TO 100
! C        Save the first date just in case the header JDT1 is wrong
  120 IF (IFIRST.EQ.0) JDT1 = JD1
! C
! C     Is JD within the range of the table?
! C
  130 IF (JD.LT.JDT1.OR.JD.GT.JDT2) GO TO 500
! C     Is JD too close to JD1 ?
      IF (JD.GE.JD1) GO TO 170
! C        If so, backspace and try again
      if( nback.gt.0 .and. nrec.ge.nrecs ) lap = nrbmin+1+nrec-nrecs
      nrecs = nrec
      nback = nback + 1
      if( nback.gt.100) then
! c       get calling program name and m-file name for status_update
        len = rcpar(0,prog_name)
        write(message,'(a,i8,a,2i8)') 'Infinite loop in UT1RED, jd=',jd,' jd1,jd2=',jd1,jd2
        call status_update('FATAL',prog_name,'lib/ut1red',' ',message,0)
      endif
      N = (JD1-JD)/INT/NPR + lap
      IF (N.GT.NREC) N = NREC
      DO 160 I=1,N
  160 BACKSPACE IUT1
      NREC = NREC - N
      NVTOT = 0
      GO TO 100
! C
! C        Is JD too close to JD2 ?
  170 IF (JD.LT.JD2) GO TO 200
! C        If so, shift storage and read another record
      NVTOT = NVTOT - NVR
      DO 180 I=1,NVTOT
  180 UTVAL(I) = UTVAL(NPR+I)
      JD1 = JD1 + NPR*INT
      GO TO 100
! C
! C        Calculate interpolation times and value of tabular points
  200 T = JD - JD1
      T = (T + FRACT)/INT
      IT = T
      T = DMOD(T,1.D0)
      S = 1.D0 - T
      IF (IT.EQ.ITSAV) GO TO 230
      DO 210 I=1,4
      J = IT + I
  210 TAB(I) = UTVAL(J)
! C
! C
! C        Calculate interpolation Y-vector
! C
      DO 220 I=1,2
      NR = I + 1
      F2 =     0.166666666666667D0 * (TAB(NR+1)+TAB(NR-1))
      Y1(I) =  1.333333333333333D0 * TAB(NR) - F2
  220 Y2(I) = -0.333333333333333D0 * TAB(NR) + F2
      ITSAV = IT
! C
! C
! C        Second difference interpolation
! C
! c 230 TAIUT1 = (T* (Y1(2)+T*T*Y2(2)) + S * (Y1(1)+S*S*Y2(1)))*1.D-5
  230 TAIUT1 = (T* (Y1(2)+T*T*Y2(2)) + S * (Y1(1)+S*S*Y2(1)))*rmult
      ut1dot = ( y1(2)+3.d0*t*t*y2(2) - y1(1)-3.d0*s*s*y2(1) )*rmult/int
      GO TO 999
! C
! C
! C        Out of range of table
! C
  500 WRITE(message,510) JD,JDT1,JDT2
  510 FORMAT('JD= ',I7,' out of range of ut1 table, JDT1= ',I7,' JDT2= ',I7)
      call status_update('FATAL',prog_name,'lib/ut1red',' ',message,0)
  520 WRITE(message,530) MJD,JD2,JDT2
  530 FORMAT('File error in UT1 table, MJD = ',I5,' JD2 = ',I7,' JDT2 = ',I7)
      call status_update('FATAL',prog_name,'lib/ut1red',' ',message,ioerr)
  540 call status_update('FATAL',prog_name,'lib/ut1red',' ','Cannot find format statement in ut1 file:',0)
  999 iuttyp = itype

      RETURN
      END


