  subroutine calc_perturbed_acc(norbs_perturbed,onepr,twopr,orbit_perturb,satpos,accel,jd,tin)

! Subroutine to adjust the acceleration computed at the satellite position for the difference in static gravity
! acceleration between the actual location of the satellite and the perturbed location of the satellite. To do
! this, we calculate the static gravity accelerations at both locations, difference them, rotate to inertial space
! and then adjust by the difference the acceleration that was computed by gravcalc.
!
! Initially we only do this for the Z coordinate, but it can eventually be expanded to XYZ and perhaps to other IC parameters.
!
! P. Tregoning
! 12 August 2014

  use bsscl_mod     ! required to get the a priori bias and scale information for accelerom.f90
  use usno_mod      ! provides the Xp, Yp etc information
  use rotation_mod  ! provides the rotation matrices efixed-inertial
  use inmod_mod     ! provides info on whether we want accelerometer obs to be used

  implicit none

  integer*4      , intent(in)    :: norbs_perturbed                  ! number of perturbed orbits that we need to keep track of
  real(kind=8)   , intent(in)    :: onepr(2)                         ! IC amplitudes of twice-per-rev along-track accelerations
  real(kind=8)   , intent(in)    :: twopr(2)                         ! IC amplitudes of twice-per-rev along-track accelerations
  real(kind=8)   , intent(in)    :: orbit_perturb(norbs_perturbed)   ! perturbations of each IC parameter 
  real(kind=8)   , intent(in   ) :: satpos(3+norbs_perturbed*3)      ! first three are inertial XYZ of sat, then X-perturbed XYZ, then Y-perturbed XYZ then Z-perturbed XYZ
  real(kind=8)   , intent(inout) :: accel(3+norbs_perturbed*3)       ! first three are inertial XYZ accelerations from GRAVCALC, then the acceleration partials for each of 
                                                    ! the position perturbations. We want to replace the partials with the accelerations at the perturbed locations

  integer*4   ,intent(in)    :: jd                 !  julian day
  real(kind=8),intent(in)    :: tin                !  seconds of day

! local variables
  real(kind=8) :: sat_ef(6)                        ! earth-fixed position of the satellite
  real(kind=8) :: sat_perturb_ef(6)                ! earth-fixed perturbed position of the satellite
  real(kind=8) :: satlat_ef,satlon_ef              ! earth-fixed lat/lon of satellite position
  real(kind=8) :: pertlat_ef,pertlon_ef            ! earth-fixed lat/lon of the perturbed satellite position
  real(kind=8) :: radius_ef                        ! satellite radius (we don't use it here)
  real(kind=8) :: statgrav_sat(3),statgrav_pert(3) ! static gravity at the satellite and perturbed satellite locations
  real(kind=8) :: accel_diff_ef(3)                 ! difference in earth-fixed static gravity acceleration at the two locations
  real(kind=8) :: accel_diff_if(3)                 ! difference in inertial    static gravity acceleration at the two locations
  real(kind=8) :: accel_diff_scale_if(3)           ! computed (inertial) differences in non-gravitational accelerations because of scale perturbations
  real(kind=8) :: accel_diff_bias_if(3)            ! computed (inertial) differences in non-gravitational accelerations because of bias perturbations
  real(kind=8) :: accel_diff_rpy_if(3)             ! computed (inertial) differences in non-gravitational accelerations because of roll,pitch,yaw perturbations
  real(kind=8) :: x
  real(kind=8) :: amag3
  real(kind=8) :: tmp_pos(6)

! accelerometer values
  real(kind=8) :: acc_calib_if(3)                  ! calibrated non-gravitational accelerations in inertial XYZ directions
  real(kind=8) :: acc_calib_no2pr_if(3)            ! calibrated non-gravitational accelerations in inertial XYZ directions wihtout twice-per-rev
  real(kind=8) :: acc_pert_scale_if(3)             ! non-gravitational accelerations in inertial XYZ directions perturbed by a scale error
  real(kind=8) :: acc_pert_bias_if(3)              ! non-gravitational accelerations in inertial XYZ directions perturbed by a bias error
  real(kind=8) :: acc_pert_rpy_if(3)               ! non-gravitational accelerations in inertial XYZ directions perturbed by a roll,pitch,yaw error
  real(kind=8) :: tmp_bias(2)                      ! parameter to allow us to perturb the relevant bias temporarily
  real(kind=8) :: tmp_scale(3)                     ! array of temporary scales to allow us to perturb the relevant bias temporarily
  integer*4    :: bias_epoch,bias_ptr
  real(kind=8) :: bias_perturb(3)                  ! array to pass a bias perturbation into accelerom.f90
  real(kind=8) :: rpy_perturb(3)                   ! array to pass a rpy  perturbation into accelerom.f90
  real(kind=8) :: scl_perturb(3)                   ! array to pass a scl  perturbation into accelerom.f90

! variables for a once-per-revolution orbit perturbation in along-track
  real(kind=8) :: acc_pert_1pr_if(3)               ! array of accelerations perturbed in the along-track direction by a once-per-rev   
  real(kind=8) :: accel_diff_1pr_if(3)             ! non-gravitational accelerations in inertial XYZ directions perturbed by a 1pr along-track error
  real(kind=8) :: oneperrev_ampl                   ! amplitude of the once-per-rev empirical function for this epoch

! variables for a twice-per-revolution orbit perturbation in along-track
  real(kind=8) :: acc_pert_2pr_if(3)               ! array of accelerations perturbed in the along-track direction by a twice-per-rev   
  real(kind=8) :: accel_diff_2pr_if(3)             ! non-gravitational accelerations in inertial XYZ directions perturbed by a 2pr along-track error
  real(kind=8) :: twoperrev_ampl                   ! amplitude of the twice-per-rev empirical function for this epoch
  real(kind=8) :: secs_per_rev                     ! temporary means of getting the period of the satellite (ie hardwire it)
  real(kind=8) :: pi

! PT170821: variable required in the call of efixed_inert (but not used in computations from this subroutine)
    real(kind=8) :: ifacc(3)

! counters
  integer*4    :: i,j

  pi = 4.d0*datan(1.d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!  *************                     actual satellite orbit                       ************* 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ok, convert the satellite location to earth-fixed. We ignore the values 4:6 because we don't need to use the velocities here.
  call efix_inert(rot_i2e, rotdot_i2e, rotacc_i2e ,satpos(1:6), sat_ef(1:6), ifacc, .false. )
! convert to lat/lon
  call cart_sph_convert(satlat_ef,satlon_ef,radius_ef,sat_ef(1:3) )
!  print*,'sat_ef',sat_ef

! compute the static gravity field accelerations at the satellite location
  x = dsin(satlat_ef)
  call legcalc_norm(x)
  call sphharmfield(satlat_ef, satlon_ef, radius_ef, jd, tin,  statgrav_sat ) 
!  print*,'statgrav_pert',statgrav_pert

! compute the default calibrated inertial accelerometer accelerations in XYZ
  bias_perturb = 0.d0
  rpy_perturb = 0.d0
  scl_perturb = 0.d0
! PT140820: just use tin (seconds since 2000/01/01 1200UT) as our time for the twice-per-rev. Hardwire the period to 93 mins for now !!!!
  secs_per_rev = 93.d0*60.d0  ! 93 minute period  x 60 seconds/minute
! PT150814: the 1/rev and 2/rev values in the gracefit output file are already in um/s^2, so don't multiply here by 1e-6
  oneperrev_ampl = onepr(1)*1.d-6 * dsin(2.d0*pi*tin/secs_per_rev) &
                                  + onepr(2)*1.d-6 * dcos(2.d0*pi*tin/secs_per_rev)
  twoperrev_ampl = twopr(1)*1.d-6 * dsin(4.d0*pi*tin/secs_per_rev) &
                                  + twopr(2)*1.d-6 * dcos(4.d0*pi*tin/secs_per_rev)
!  print*,'twoperrev_ampl = ',twoperrev_ampl
  bias_perturb(1) = oneperrev_ampl + twoperrev_ampl
  acc_calib_if = 0.d0
  if(gt_acc(1) == "Y")call accelerom(jd, tin, bias_perturb,c0x,c0y,c0z,c1x,c1y,c1z,c2x,c2y,c2z,scl,scl_perturb,rpy_perturb &
                                    ,acc_calib_if)
! PT140820: we now update the three accelerations for the twice-per-rev along-track accelerations via the code in accelerom.f90
  bias_perturb = 0.d0
  acc_calib_no2pr_if = 0.d0
  scl_perturb = 0.d0
  if(gt_acc(1) == "Y")call accelerom(jd, tin, bias_perturb,c0x,c0y,c0z,c1x,c1y,c1z,c2x,c2y,c2z,scl,scl_perturb,rpy_perturb &
                                    ,acc_calib_no2pr_if)
  accel(1:3) = accel(1:3) + acc_calib_if(1:3) - acc_calib_no2pr_if(1:3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!  *************                   perturbed orbits                   ************* 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do i=4,3+norbs_perturbed*3,3         ! 4,7,10 are XYZ of position perturbations, 13,16,19 are XYZ of velocity perturbations, 22,25,28 are scale, 31,34,37 are bias
    tmp_pos = 0.d0
    tmp_pos(1:3) = satpos(i:i+2)
!print*,'** tin,i,tmp_pos',tin,i,tmp_pos(1:3),tmp_pos(1:3) - satpos(1:3),amag3(tmp_pos(1:3) - satpos(1:3))
    call efix_inert(rot_i2e, rotdot_i2e, rotacc_i2e ,tmp_pos, sat_perturb_ef(1:6), ifacc, .false. )
    call cart_sph_convert(pertlat_ef,pertlon_ef,radius_ef,sat_perturb_ef(1:3) )
! get the static gravity accelerations at the perturbed location
    x = dsin(pertlat_ef)
    call legcalc_norm(x)
    call sphharmfield(pertlat_ef, pertlon_ef, radius_ef, jd, tin,  statgrav_pert ) 
!  print*,'statgrav_sat ',statgrav_sat
!  print*,'statgrav_pert',statgrav_pert

! compute the difference in static gravity accelerations at the satellite and perturbed satellite locations
    accel_diff_ef = statgrav_pert - statgrav_sat
! rotate the acceleration difference to inertial space
    do j=1,3  
        accel_diff_if(j) = rot_e2i(j,1)*accel_diff_ef(1)+rot_e2i(j,2)*accel_diff_ef(2)+rot_e2i(j,3)*accel_diff_ef(3)  
    enddo
   

!!!!!!!!!!!!!!!!!!!!
! perturbed scales !
!!!!!!!!!!!!!!!!!!!!
    accel_diff_scale_if = 0.d0
    rpy_perturb = 0.d0
    bias_perturb = 0.d0
    bias_perturb(1) = oneperrev_ampl + twoperrev_ampl
    scl_perturb(1:3) = 0.d0

    if (i>= 22 .and. i <= 28)then
      if(i == 22 )then
        tmp_scale = scl
        scl_perturb(1) = orbit_perturb(7)
      elseif(i == 25)then
        tmp_scale = scl
        scl_perturb(2) = orbit_perturb(8)
      elseif(i == 28)then
        tmp_scale = scl
        scl_perturb(3) = orbit_perturb(9)
      endif
      acc_pert_scale_if = 0.d0
      if(gt_acc(1) == "Y")call accelerom(jd, tin, bias_perturb,c0x,c0y,c0z,c1x,c1y,c1z,c2x,c2y,c2z,tmp_scale &
         ,scl_perturb,rpy_perturb,acc_pert_scale_if)
! now calculate the difference in non-gravitational acceleration between the actual and the perturbed non-grav accelerations
      accel_diff_scale_if = acc_pert_scale_if - acc_calib_if
!print*,'jd,tin,i,acc_pert_scale_if',jd,tin,i,acc_pert_scale_if
    endif

!!!!!!!!!!!!!!!!!!!!
! perturbed biases !
!!!!!!!!!!!!!!!!!!!!
    accel_diff_bias_if = 0.d0
    rpy_perturb = 0.d0
    bias_perturb = 0.d0
    bias_perturb(1) = oneperrev_ampl + twoperrev_ampl
    scl_perturb = 0.d0
    if (i >= 31 .and. i <= 37)then
      if(i == 31 )then
        bias_perturb(1) = orbit_perturb(10)*1.d-6
      elseif(i == 34)then
        bias_perturb(2) = orbit_perturb(11)*1.d-6
      elseif(i == 37)then
        bias_perturb(3) = orbit_perturb(12)*1.d-6
      endif
      acc_pert_bias_if = 0.d0 
      if(gt_acc(1) == "Y")call accelerom(jd, tin, bias_perturb,c0x,c0y,c0z,c1x,c1y,c1z,c2x,c2y,c2z,scl,scl_perturb,rpy_perturb &
                                         ,acc_pert_bias_if)

! now calculate the difference in non-gravitational acceleration between the actual and the perturbed non-grav accelerations
     accel_diff_bias_if = acc_pert_bias_if - acc_calib_if
    endif

!!!!!!!!!!!!!!!!!!!!
!  once-per-rev   !
!!!!!!!!!!!!!!!!!!!!
    accel_diff_1pr_if = 0.d0
    rpy_perturb = 0.d0
    scl_perturb = 0.d0
    if (i >= 40 .and. i <= 43)then
      bias_perturb = 0.d0
      if(i == 40 )then
        bias_perturb(1) = oneperrev_ampl + twoperrev_ampl + orbit_perturb(13)*1.d-6 * dsin(2.d0*pi*tin/secs_per_rev)
      else if(i == 43)then
        bias_perturb(1) = oneperrev_ampl + twoperrev_ampl + orbit_perturb(14)*1.d-6 * dcos(2.d0*pi*tin/secs_per_rev)
      endif
      acc_pert_1pr_if = 0.d0
      if(gt_acc(1) == "Y")call accelerom(jd, tin, bias_perturb,c0x,c0y,c0z,c1x,c1y,c1z,c2x,c2y,c2z,scl,scl_perturb,rpy_perturb &
                                        ,acc_pert_1pr_if)

! now calculate the difference in non-gravitational acceleration between the actual and the perturbed non-grav accelerations
     accel_diff_1pr_if = acc_pert_1pr_if - acc_calib_if
    endif

!!!!!!!!!!!!!!!!!!!!
!  twice-per-rev   !
!!!!!!!!!!!!!!!!!!!!
    accel_diff_2pr_if = 0.d0
    rpy_perturb = 0.d0
    scl_perturb = 0.d0
    if (i >= 46 .and. i <= 49)then
      bias_perturb = 0.d0
      if(i == 46 )then
        bias_perturb(1) = oneperrev_ampl + twoperrev_ampl + orbit_perturb(15)*1.d-6 * dsin(4.d0*pi*tin/secs_per_rev)
      else if(i == 49)then
        bias_perturb(1) = oneperrev_ampl + twoperrev_ampl + orbit_perturb(16)*1.d-6 * dcos(4.d0*pi*tin/secs_per_rev)
      endif
      acc_pert_2pr_if = 0.d0
      if(gt_acc(1) == "Y")call accelerom(jd, tin, bias_perturb,c0x,c0y,c0z,c1x,c1y,c1z,c2x,c2y,c2z,scl,scl_perturb,rpy_perturb &
                                        ,acc_pert_2pr_if)

! now calculate the difference in non-gravitational acceleration between the actual and the perturbed non-grav accelerations
     accel_diff_2pr_if = acc_pert_2pr_if - acc_calib_if
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!
!  orientation errors   !
!!!!!!!!!!!!!!!!!!!!!!!!!
    accel_diff_rpy_if = 0.d0
    rpy_perturb = 0.d0
    scl_perturb = 0.d0
    if (i >= 52 .and. i <= 58)then
      bias_perturb = 0.d0
      if(i == 52 )then
! roll
        rpy_perturb(1)  = orbit_perturb(17)*1.d0
      else if(i == 55 )then
! pitch
        rpy_perturb(2)  = orbit_perturb(18)*1.d0
      else if(i == 58 )then
! yaw
        rpy_perturb(3)  = orbit_perturb(19)*1.d0
      endif
      acc_pert_rpy_if = 0.d0
      if(gt_acc(1) == "Y")call accelerom(jd, tin, bias_perturb,c0x,c0y,c0z,c1x,c1y,c1z,c2x,c2y,c2z,scl,scl_perturb,rpy_perturb &
                                        ,acc_pert_rpy_if)

! now calculate the difference in non-gravitational acceleration between the actual and the perturbed non-grav accelerations
     accel_diff_rpy_if = acc_pert_rpy_if - acc_calib_if
!print*,'jd,tin,rpy_perturb,acc_diff',jd,tin,rpy_perturb,acc_pert_rpy_if-acc_calib_if
    endif

! Now, the accelerations at the location of the perturbed orbit are found by adding the difference to the GRAVCALC accelerations
    accel(i:i+2) = accel(1:3) + accel_diff_if(1:3) + accel_diff_scale_if(1:3) + accel_diff_bias_if(1:3) &
                   + accel_diff_1pr_if(1:3) + accel_diff_2pr_if(1:3) + accel_diff_rpy_if(1:3)


  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  stop 'stopped at bottom of calc_perturbed_acc'
  return

  end subroutine calc_perturbed_acc

    
    
