  subroutine gravcalc(jd ,t, pos, vel, gtot, efpos, output_tide)

! Emma-Kate Potter, 10 May 2010   
! This program calls for the calculation of each contribution to 
! acceleration experienced by the satellite:

! 1) static field of earth
! 2) sun and moon point mass
! 3) solid earth tide
! 4) atmospheric mass contribution
! 5) ocean tides

! (6) Added mass cons: Initially just a single one at a set
!     location (stored in accel_mod)

    use dealias_mod
    use spherhar_mod
    use coeff_mod
    use gm_mod
    use accel_mod    ! Add accel_mod which contains mass-cons
    use inmod_mod
    use tide_mod     ! added for binary ocean tide height grid variables
    use sat_mod      ! added to have the satellite name to output status info
    use bsscl_mod    ! added to have the bias/scale variables to then pass into accelerom
    use mascon_mod   ! new mascon variable definitions
    use rotation_mod ! provides efixed-inertial matrices
    implicit none

    integer*4,    intent(in) :: jd                 ! julian day
    real(kind=8), intent(in) :: t                  ! seconds of julian day
    real(kind=8), intent(in) :: pos(3)             ! inertial position
    real(kind=8), intent(in) :: vel(3)             ! inertial velocity
    real(kind=8), intent(out):: gtot(3)            ! total inertial accelerations computed in this subroutine
    real(kind=8), intent(out):: efpos(3)           ! earth-fixed position
    logical     , intent(in) :: output_tide        ! flag as to whether to print ocean tide info to the screen

    integer*4 :: i, j, k, m, idir
    real(kind=8), dimension(3) :: gsphharm , gravsm, gravsolid
    real(kind=8), dimension(3) :: g_in, g_ef
    real(kind=8), dimension(3) :: efvel, efgstatic, efgsolid, accmtr, efgsphharm
    real(kind=8) :: rx, rlat, rlong, x, mcon_fact
    integer*4 :: jds
    integer :: lucoeff, ioerr, counter, ideg, iord, ierr
    real(kind=8) :: cval, sval, temp1, temp2 
    real(kind=8), dimension(3,3) :: rot_t, angvelten
    real(kind=8), dimension(3) :: angvel, grelacc
    real(kind=8), dimension(3) :: apulse         
    real(kind=8), dimension(6) :: sunpos
    integer :: tarray(3)
    real*8 :: PEPt
    integer*4 :: PEPjd 

! variables required in the call to rotsnp (although not used in gravcalc)
    real(kind=8) :: xpole,ypole

!  EWH Mascon variables
    real*8 mcon_dist     ! Distance from satellite to Mass con
    real*8 def_fact      ! scaling factor for the deformational components
    real*8 mcon_efacc(3) ! mcon earth fixed acceleration
    real*8 mcon_ifacc(3) ! mcon earth inertial acceleration
!  Ocean mascon variables
    real(kind=8) :: ocean_efacc(3)   ! sum of accelerations from the ocean tide

! PT130619: variables to output the eclipse information
    double precision shadow

! PT130917: variables for interpolating the ocean tide height value
    integer*4 point_num
    integer*4 irec
    integer*4 lat_band
    real(kind=8 ) dtime

! PT140428: variables for the partials of mascon acceleration wrt ewh height and position
    real(kind=8) mcon_vol, tmp_part(3), bias_perturb(3), rpy_perturb(3), scl_perturb(3)
    integer*4 l
    logical bitmap

! debug variables
    real*8 :: tmp_mcon_efacc(3),tmp_mcon_efacc_tide(3)
    real*8 :: debug_mcon_efacc(3),debug_otide_efacc(3)
    real*8 :: tide_volume
    logical :: debug

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   initialise contributions to gravity 
    gtot = 0.d0
    gsphharm = 0.d0
    gravsolid = 0.d0
    grelacc = 0.d0
    g_in = 0.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   CALCULATE THE SATELLITE POSITION IN EARTH FIXED COORDINATES
!   pos is ***inertial coordinates*** of satellite
    do i=1,3
      efpos(i) = rot_i2e(i,1)*pos(1)+rot_i2e(i,2)*pos(2)+rot_i2e(i,3)*pos(3)
      efvel(i) = rot_i2e(i,1)*vel(1)+rot_i2e(i,2)*vel(2)+rot_i2e(i,3)*vel(3) &
                +rotdot_i2e(i,1)*pos(1)+rotdot_i2e(i,2)*pos(2)+rotdot_i2e(i,3)*pos(3)
    enddo
    write(41,*)  efpos(1), efpos(2), efpos(3)

!   CALCULATE lat and long for earthfixed position vector
    call cart_sph_convert(rlat, rlong, rx, efpos)

! PT171003: put a trap here if the radius is less than the radius of the Earth (check on units of input ICs)
    if(rx < 6378137.d0)then
      write(message,'(a,f15.4,a)')"Error: satellite radius (",rx," m) is less than Earth radius. Is your satellite flying?"
      call status_update('FATAL','GRACEORB','gravcalc',' ',message,0)
    endif

!   COMPUTE the legendre polynomials for this spatial coordinate.
!   They are returned in Plm in spherhar_mod1 
!   set variable to feed to legendre calculation
    x = dsin(rlat)
    call legcalc_norm(x)
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   POINT EVALUATIONS OF GRAVITY CONTRITBUTION
!   These are calculations of the gravity contribution at a point (not as a
!   spherical harmonic representation of a field)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   CALCULATE  contribution of sun and moon and planets. This includes:
!   1) the direct point mass attraction of the planets on the satellite
!   2) the indirect effect of the point mass attraction of the planets on the
!      earth
!   3) the indirect J2 effect due to a point mass sun and moon on the
!      oblateness of the earth
    call planetfield(jd, t, pos, meancoefC(2,0), sunpos, gt_planet(1), gravsm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   CALCULATE solid earth tide contribution to gravity experienced by satellite
!   NOTE that the solidfield routine must be called after the planetfield
!   subroutine above because it uses the sun and moon coordinates determined
!   in that subroutine
!    call solidfield(rlat, rlong, rx, jd, t, rot_i2e, efgsolid)
    call solidfield_IERS2010(rlat, rlong, rx, jd, t, rot_i2e, efgsolid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    debug = .false.
    if(debug)print*,'efgsolid,lat,lon',efgsolid,rlat*180./pi,rlong*180./pi,jd,t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   CALCULATE contribution from non-gravitational accelerations (from 
!   accelerometer data).
    bias_perturb = 0.d0
    rpy_perturb = 0.d0
! DEBUG: PT220124: insert a +0.5 mrad bias in the pitch, just to see what happens!
!    rpy_perturb(2) = -0.5d-3

    scl_perturb = 0.d0
! uncomment to set a bias in the rpy orientations
!    rpy_perturb(1) = 2.d-3  
!    rpy_perturb(2) = 2.d-3  
!    rpy_perturb(3) = 2.d-3  
! PT171003: only call this if we want to use the accelerometer data
    if (gt_acc(1).eq."Y")then
      call accelerom(jd, t, bias_perturb, c0x,c0y,c0z,c1x,c1y,c1z,c2x,c2y,c2z,scl,scl_perturb,rpy_perturb, accmtr)
    else
      accmtr = 0.d0
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   CALCULATE the acceleration perturbation due to general relativistic effects
!   IERS 2003 expression is used here:
!   ******************* STILL ISSUES with general rel calculation **************
!   ***************related to the order of the cross products - need to check***
!   ******************but very very small effect on the final acc regardless****
!  PT140130: but it builds to 1.5 m oscillating effect after 12 hours of integration !!!
    call trnsp(rot_e2i, rot_t, 3, 3)
    call matmult(rotdot_e2i, rot_t, angvelten, 3, 3, 3)
! PT140131: angvelten is the angular velocity tensor, with the rotation angles about x (3,2), y(1,3) and z (2,1).
!           It also contains the -ives of them, and has zeroes down the diagonal.
    angvel(1) = angvelten(3,2)
    angvel(2) = angvelten(1,3)
    angvel(3) = angvelten(2,1)
    call generalrel(jd,t,-dsqrt(5.d0)*meancoefC(2,0) ,pos, vel, sunpos, angvel, grelacc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(debug)print*,'grelacc,pos,vel',grelacc,pos,vel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   CONTRIBUTIONS TO GRAVITY FROM SPHERICAL HARMONIC MODELS/RECONSTRUCTIONS
!   The sphharmfield subroutine includes:
!   1) static/mean gravity field
!   2) ocean tides
!   3) atmospheric tides
!   4) dealiasing (non-tidal atm and ocean)
    call sphharmfield(rlat, rlong,rx, jd, t, efgsphharm) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(debug)print*,'efgsphharm,lat,lon',efgsphharm,rlat*180./pi,rlong*180./pi
!DEBUG
! print the accelerations to date
if(dabs(efgsphharm(1)) > 100.d0)then
  write(message,'(a,3e16.9)')'static gravity acceleration too big',efgsphharm
  call status_update('FATAL','GRACEORB','gravcalc',' ',message,0)
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate Mascon accelerations and partial derivatives

! PT121221: merge of Tony's code that can exclude mascon computations
! PT210805: pass the flag into the subroutine to separate the need for the computation of accelerations and partials
  if (gt_mcon(1).eq."Y" .or. gt_mcon(1).eq."O") then
    call calc_mascon_acc_part(efpos,t,mcon_efacc)
  endif
  if(debug)print*,'mcon_efacc,efpos(1:6)',mcon_efacc,efpos
!print*,'returned in gravcalc from calc_mascon_acc_part'
!stop
! PT160909: check whether we are applying the ocean tide through the ocean mascon approach
  if(gt_oceantide(1).eq.'G') then
    call calc_ocean_acc_part(sat,rlat,rlong,efpos,jd,t,ocean_efacc,output_tide)
  endif
  if(debug)print*,'ocean_efacc,lat,lon',ocean_efacc ,rlat*180./pi,rlong*180./pi 

! END MASS CON Code          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ADD the e-fixed contributions to acceleration together then rotate
!   - all spherical harmonic models (efgsphharm)
!   - solid earth tide (efgsolid)
!   - mascons (mcon_efacc)
  g_ef = 0.d0
  g_ef(:) = g_ef(:) + efgsphharm(:)
  if (gt_solidtide(1).eq."Y") g_ef(:) = g_ef(:) + efgsolid(:)
  if (gt_mcon(1).eq."Y" .or. gt_mcon(1).eq."O")      g_ef(:) = g_ef(:) + mcon_efacc(:)
  if (gt_oceantide(1) == "G")      g_ef(:) = g_ef(:) + ocean_efacc(:)

! Now rotate this sum of efixed contributions into inertial space
  do i=1,3  
    g_in(i) = rot_e2i(i,1)*g_ef(1)+rot_e2i(i,2)*g_ef(2)+rot_e2i(i,3)*g_ef(3)  
  enddo
 
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ADD all gravity contributions to give total gravity experienced by satellite
    gtot = 0.d0

! DEBUG
!print*,'gravsm',jd,t,gravsm
!print*,'grelacc',jd,t,grelacc
!print*,'accmtr',jd,t,accmtr
!print*,'mcon_efacc',jd,t,mcon_efacc
!print*,'ocean_efacc',jd,t,ocean_efacc
!print*,'efgsphharm',jd,t,efgsphharm
!print*,'efgsolid',jd,t,efgsolid


    do i = 1, 3
      gtot(i) = gtot(i) + g_in(i)
      if (gt_sunmoon(1).eq."Y")    gtot(i) = gtot(i) + gravsm(i)
      if (gt_genrel(1).eq."Y")     gtot(i) = gtot(i) + grelacc(i)
      if (gt_acc(1).eq."Y")        gtot(i) = gtot(i) + accmtr(i)
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   return
   end



