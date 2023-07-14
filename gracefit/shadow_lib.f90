!********************************************************************************************************************************
!  File: shadow_lib.f90
!
!  Purpose: Set of subroutines used concerning eclipse conditions 
!
!  Author: Thomas Greenspan
!
!  API:
!       shadow_satelliteStatus      :  Sets indicator array, shadow, to show whether each satellite is in 
!                                      complete shadow (calls shadow_isInShadow)
!       shadow_isInShadow           :  Sets indicator array, shadow, to show whether each satellite is in
!                                      complete shadow (calls get_lambda)
!       shadow_conditionApplicable  :  Sets whether the shadow condition can be applied for a given epoch
!       shadow_applyCondition       :  Applies the shadow condition (as indicated by input)
!       shadow_SRP                  :  Calculate the start and stop times of when the solar radiation pressure is 
!                                      either applied or removed.
!       shadow_findAccObs           :  find a non-thrust-affected acc obs that is roughly half a revolution away from a shadow observation
!       leadsat                     :  determine which is the lead and which is the trail satellite
!
!   July 23, 2013
!
!********************************************************************************************************************************
!********************************************************************************************************************************
! shadow_satelliteStatus: Sets the status array shadow to show whether the satellites are in
!                         shadow given the julian date and the satellite positions.
! Author: Paul Tregoning
!         (Modified by Thomas Greenspan)
!********************************************************************************************************************************

subroutine shadow_satelliteStatus(iepoch,rot_i2e,beta,sat_pos,sciframe_quat,shadow,sunpos_efixed)
  use gracefit_mod
  implicit none

  include 'input.h90'

  !********************  Variable declarations ********************

  integer*4, intent(in)         :: iepoch             ! Epoch in loop at which the subroutine is called
  real(kind=8), intent(in)      :: rot_i2e(3,3)       ! inertial to earth-fixed rotation matrix
  real(kind=8), intent(in)      :: beta               ! current beta angle (for debug)
  double precision, intent(in)  :: sat_pos(6,nsat_t)  ! Vector with satellite positions and velocities
  real(kind=8),     intent(in)  :: sciframe_quat(4,2) ! quaternion to get from srf to inertial space
  double precision, intent(out) :: shadow(nsat_t)     ! Array containing information on whether the satellites are in shadow
  double precision, intent(out) :: sunpos_efixed(6)   ! Position of the sun in efixed frame

  integer*4        :: isat          ! Counter variable for satellites
  integer*4        :: date(5)       ! Varaible needed to calculate the julian date
  double precision :: sec           ! Variable needed to calculate the julian date
  double precision :: fjd           ! Julian date
  double precision :: fjd_dec       ! Integer and decimal parts of fjd
  double precision :: ephtime(2)    ! Integer part of epoch and fraction of a day
  double precision :: sunpos(6)     ! Position of the sun 
  double precision :: satpos(3,nsat_t)   ! inertial position of a satellite
  double precision :: satvel(3,nsat_t)   ! inertial velocity of a satellite
!  double precision :: PEPt                 
!  integer*4        :: PEPjd
!  integer*4        :: idir, iut1pol ! Values to be input to rotsnp
!  double precision :: rot(3,3)      ! Rotation matrix to be used in satellite position calculation
!  double precision :: rotderiv(3,3) ! Derivative of above matrix (Not used here but needed as input for rotsnp)
!  double precision :: sidtm(3,3)    ! Sidereal time (rad) (Not used here but needed as input for rotsnp)

  real(kind=8)     :: dot_angle(3)    ! dot products between the inertial satellite XYZ axes vectors and the sun-sat vector
  real(kind=8)     :: vect_srf(3)     ! unit vector of the satellite bus axes (SRF frame)
  real(kind=8)     :: vect_efixed(3)  ! unit vector of the satellite bus axes (inertial frame)
  real(kind=8)     :: sun_sat_vect(3) ! sun to satellite vector
  real(kind=8)     :: mag_sunsat      ! magnitude of sun to satellite vector
  real*8           :: dot             ! dot product function
  integer*4        :: i
  real(kind=8)     :: orb_norm(3),amag3,beta_angle
!****************************************************************

! Compute Julian date
! PT150612: this time computation is now done in shadow_getSun
!    call gsec_to_ymdhms (starting_epoch, date, sec)
!    call ymdhms_to_jd(date,sec,fjd)

! Set Values
!    fjd_dec = fjd - int(fjd)                  ! decimal part of fjd
!    ephtime(1) = dble(int(fjd))               ! integer part of fjd
!    ephtime(2) = fjd_dec + ((iepoch-1)*epoch_interval+ tdtoff)/SEC_IN_DAY  ! fraction of a day   (tdtoff included as per graceorb/planetfield.f90)
!****************************************************************

! PT150525: get the sun coordinates wrt in earth-fixed frame
  call shadow_getsun(starting_epoch,iepoch,1,sunpos)  ! sunpos (in metres)  and is returned as sun relative to earth in inertial coords

! PT200706: now rotate inertial Sun coords to earth-fixed 
  call matmult(rot_i2e,sunpos(1:3),sunpos_efixed(1:3),3,3,1)
  do isat = 1, nsat_t
     ! PT150525: change this to using the sun in e-fixed rather than satellite in inertial
     call shadow_isInShadow(sunpos_efixed(1:3),sat_pos(:,isat),shadow(isat))
  enddo

  !****************************************************************


  return
end subroutine shadow_satelliteStatus

!********************************************************************************************************************************
!********************************************************************************************************************************
! shadow_isInShadow: Calculates whether a GRACE satellite is in the Earth's shadow
!                    given a vector of the sun's position. a vector of the satellite's
!                    position and the value that will be set
!                    NOTE: shadow is set between 0.0 (shadow) and 1.0 (full sunlight)
! Author: Paul Tregoning
!********************************************************************************************************************************

subroutine shadow_isInShadow(sunpos,satpos,shadow)

  implicit none

  !********************  Variable declarations ********************
  ! PT150522: fixed up the definitions of which coordinates are with respect to what
  double precision, intent(in)  :: sunpos(3)   ! Position of the sun (m) in earth-fixed frame              ! PT150522: this is sun wrt earth
  double precision, intent(in)  :: satpos(3)   ! Position of a satellite (m)                               ! PT150522: this is sat wrt earth
  double precision, intent(out) :: shadow      ! Variable set to indicate whether satellite is in shadow

  double precision, parameter :: EARTH_RAD = 6378136.3d0  ! value taken from arc/egm08.f 
  double precision, parameter :: SUN_RAD   = 696000.D+03  ! value taken from arc/egm08.f
  double precision :: dot                      ! Function that returns the dot product of two vectors
  double precision :: earthsat_vec(3)          ! Vector from the satellite to the earth
  double precision :: mag_earthsat             ! Magnitude of earth to satellite vector
  double precision :: sunsat_vec(3)            ! Vector from the sun to the satellite
  double precision :: mag_sunsat               ! Magnitude of sun to satellite vector
  double precision :: earth_ang                ! Angle between satellite and earth
  double precision :: sun_ang                  ! Angle between satellite and sun
  double precision :: sunsatearth_ang          ! Angle between sun-sat and earth-sat vectors
  integer*4 :: idbprt,idebug                   ! Constants needed as input for get_lambda subroutine
  !****************************************************************

  mag_sunsat = 0.d0
  mag_earthsat = 0.d0 

  ! calculate the satellite-sun vector, the magnitude of the earth-sat and sun-sat vectors
  ! sunpos (in km) : sun wrt earth.   satpos: satellite wrt earth. Therefore, sunsat_vec (sun to sat) = satpos - sunpos
  sunsat_vec = satpos - sunpos
  mag_sunsat = dsqrt(sum(sunsat_vec(1:3)**2))
  mag_earthsat = dsqrt(sum(satpos(1:3)**2))
  ! change the sign of satpos so that it is the earth to satellite vector
  ! PT150522: no, it is already the Earth to satellite vector
  !  earthsat_vec = -satpos
  earthsat_vec = satpos

  ! get the angles subtended at the satellite by half the disc of the Earth. Use a different formula from that in arc/shadow, because
  ! the GRACE satellites are so much closer to the Earth.
  earth_ang = dasin(EARTH_RAD/mag_earthsat)

  ! we can assume that the sun angle will be small, so this arc approximation is ok
  sun_ang = SUN_RAD/mag_sunsat

  ! calculate the angle at the satellite between the sat-sun and sat-earth vectors
  sunsatearth_ang = dacos( dot(earthsat_vec,sunsat_vec)/(mag_sunsat*mag_earthsat) )

  ! now, the satellite will be in eclipse if the angle subtended by half of the earth disc is greater than the sun-sat-earth angle plus the 
  ! angle subtended by half of the sun disc. But, we'll just use the gamit/arc/get_lambda.f routine to compute it.
  call get_lambda( sun_ang,earth_ang,sunsatearth_ang,shadow,idbprt,idebug )   

  !****************************************************************

  return
end subroutine shadow_isInShadow

!********************************************************************************************************************************
!********************************************************************************************************************************
! shadow_conditionApplicable: Determines whether the shadow condition should be applied for each
!                             satellite at a given epoch given the epoch number, the indicator
!                             as to whether the satellite is in shadow, the array on whether the
!                             thrusters are having any effect on acceleration measurements, the
!                             array to be set and a counter for how long each satellite has been
!                             in shadow (for current eclipse)
! Author: Thomas Greenspan
!********************************************************************************************************************************

subroutine shadow_conditionApplicable(iepoch,shadow,thrusters_off,shadow_counter,apply_shadow_cond, n_shadow_obs)
  use gracefit_mod
  implicit none


  !********************  Variable declarations ********************

  integer*4, intent(in)        :: iepoch                              ! Epoch number in which subroutine is being called
  double precision, intent(in) :: shadow(nsat_t)                      ! Variable set to indicate whether satellite is in shadow
  logical, intent(in)          :: thrusters_off(nepochs_t,nsat_t)     ! Data on whether thrusters are having any effect on accelerations
  integer*4, intent(inout)     :: shadow_counter(nsat_t)              ! Counter as to how many epochs the satellites have been in shadow
  logical, intent(out)         :: apply_shadow_cond(nepochs_t,nsat_t) ! Data on whether the shadow condition can be applied
  integer*4, intent(inout)     :: n_shadow_obs                        ! number of shadow obs to be used to apply shadow constraints

  double precision, parameter :: ECLIPSE = 0.0  ! Constant indicating the satellite is in shadow
  integer*4 :: isat  ! Counter variable for satellites
  integer*4 :: i     ! Counter variable
  !****************************************************************

  do isat = 1, nsat_t
     if(shadow(isat) == ECLIPSE)then  ! If the satellite is in complete shadow
        shadow_counter(isat) = shadow_counter(isat) + 1
        if( (shadow_counter(isat) > BUTTERWORTH_OFFSET).and.thrusters_off(iepoch,isat) ) then
           apply_shadow_cond(iepoch,isat) = .true.
           n_shadow_obs = n_shadow_obs + 1
        endif


     else
        if(shadow_counter(isat) > 0)then  ! This only happens when the satellite has just left shadow.
           do i = max(iepoch - BUTTERWORTH_OFFSET,1), iepoch - 1  ! max function in case satellite comes out of shadow in first few epochs
              apply_shadow_cond(i,isat) = .false.
           enddo
           shadow_counter(isat) = 0
        endif

     endif  ! End of "in shadow" if statement
  enddo  ! End of satellite loop

  !****************************************************************

  return
end subroutine shadow_conditionApplicable

!********************************************************************************************************************************
!********************************************************************************************************************************
! shadow_applyCondition: Places partials appropriate partials and observed minues calculated
!                        values into part and pre_omc respectively given the array determining
!                        whether to apply the conditions or not, the acceleration of each
!                        satellite, the a priori values and the two arrays mentioned above
! Author: Thomas Greenspan
!********************************************************************************************************************************

subroutine shadow_applyCondition(iepoch,apply_shadow_cond,thrusters_off,acc,apr_ScaleBias,part,pre_omc)

  ! MODS
  ! PT150609: add a new data point, where we use the a priori scale and the difference in y-acc between shadow and point separated by 1/2 revolution
  !           to define a new "known" calibrated point. We need this in order to stop the scale factor going to zero when we tightly constraint these
  !           condition equations.
  ! PT150609: pass the entire acc and apply_shadow_cond arrays into the subroutine, not just one epoch's worth of values.
  use gracefit_mod
  implicit none


  !********************  Variable declarations ********************

  integer*4 , intent(in)        :: iepoch                                  ! epoch number
  logical, intent(in)           :: apply_shadow_cond(nepochs_t,nsat_t)     ! Data on whether the shadow condition can be applied
  logical,          intent(out) :: thrusters_off(nepochs_t,nsat_t)         ! Data on whether thrusters are having any effect on accelerations
  double precision, intent(in)  :: acc(nepochs_t,4,nsat_t)                 ! Accelerations for each sattelite at every epoch
  double precision, intent(in)  :: apr_ScaleBias(max_SB,maxsat)            ! A priori values for parameters
  double precision, intent(out) :: part(nobs_t,nparam_t)                   ! Partials of observables with respect to the parameters
  double precision, intent(out) :: pre_omc(nobs_t)                         ! Observed Minus Computed values to be used in LS algorithm

  double precision, parameter :: ECLIPSE = 0.0  ! Constant indicating the satellite is in shadow
  integer*4 :: isat  ! Counter variable for satellites
  integer*4 :: i     ! Counter variable
  character*100 :: message  ! for status_update calls
  integer*4 :: second_epoch ! epoch number of an epoch half a revolution from a shadow epoch
  real(kind=8) :: y_acc_tmp ! difference between two y-acc observations, calibrated by the a priori scale factor
  real(kind=8) :: z_acc_tmp ! difference between two z-acc observations, calibrated by the a priori scale factor


  !****************************************************************
  ! NOTE:   The conditions are as follows  
  !  Y - When a satellite is in shadow, the acceleration in the
  !      y-direction is close to 0 since there is no force from
  !      the sun. The conditional equation is:
  !      Scl_y * ACC_y + Bias_y = 0.0
  !      The partials for scale and bias are thus ACC_y and 1
  !      respectively and observed is 0.0
  !  Z - When a satellite is in shadow, the acceleration in the
  !      z-direction is slightly negative since the albedo effect
  !      is (assumed) stronger than the radiation off the top of
  !      the satellite. The conditional equation is:
  !      Scl_y * ACC_y + Bias_y = -10.0
  !      The partials for scale and bias are thus ACC_y and 1
  !      respectively and observed is -10.0 (value may change)
  !    As the top of the satellite cools off, the radiation it
  !      gives off will decrease and the observed acceleration
  !      will decrease (become more negative). This is not modeled
  !      here. We use a constant with a great enough uncertainty.
  !    ->If this condition proves very useful we can refine our
  !      observed value (make time dependent, more accurate, etc...)
  !
  ! PT160405: Z - because the satellites are pitched to align with the
  !               LOS vector, the lead satellite has a -ive bias and
  !               the trailing satellite has a +ive bias because of 
  !               a small component of atmospheric drag in the SRF Z
  !               direction. So we can't expect each of them to be
  !               -10 nm/s^2 during eclipse. Rather, the mean of both
  !               will have this value.
  !****************************************************************

  ! ! ! ! ! ! ************************
  ! PT150612: the conditional shadow equations are as follows:
  !           row icond  : yacc on GRACE A during shadow
  !               icond+1: zacc on GRACE A during shadow
  !               icond+2: yacc on GRACE A 1/2 revolution from shadow obs
  !               icond+3: zacc on GRACE A 1/2 revolution from shadow obs
  !
  !               icond+ncond_t  : yacc on GRACE B during shadow
  !               icond+ncond_t+1: zacc on GRACE B during shadow
  !               icond+ncond_t+2: yacc on GRACE B 1/2 revolution from shadow obs
  !               icond+ncond_t+3: zacc on GRACE B 1/2 revolution from shadow obs

  do isat = 1, nsat_t
     if(apply_shadow_cond(iepoch,isat))then ! Make sure it is ok to apply condition
        do i = 1, 3   ! run through coordinates
           if(i == 2 )then    ! apply conditional equations on "y"
              if(nScale_t > 0 .and. nBias_t > 0)then
                 ! first, the accelerometer observation when the satellite is in shadow
                 !         Scale
                 ! PT140922: "calibrate" the accelerometer value to use in the partial (ie add the bias) so that the scale partial varies with time
                 ! PT150611: again, remove this pre-calibration
                 part(icond+ncond_t*(isat-1),iScale-1+i+nsatprm_t*(isat-1)) =  acc(iepoch,i+1,isat) ! + apr_ScaleBias(i+3,isat))   ! d[equ]/dScl_y  = ACC_y
                 !         bias
                 part(icond+ncond_t*(isat-1),iBias-1+i+nsatprm_t*(isat-1)) = 1.d0            ! d[equ]/dBias_y = 1.0
                 !         Set pre_omc
                 pre_omc(icond+(isat-1)*ncond_t) = 0.d0 - ( apr_ScaleBias(i,isat)*acc(iepoch,i+1,isat) + apr_ScaleBias(i+3,isat)) ! Observed is 0.0

                 ! PT150609: add another conditional observation here, being half a revolution away from the shadow observation. Use the a priori scale value to
                 !           scale the difference between the two obs, which then gives the "calibrated" value of the second (since the shadow ob should be zero)
                 ! find a non-thrust observation that is half a revolution away from this epoch
                 ! PT1501015: only do this if it was requested in the gracefit command file
                 if(int(use_accel_scale_const) == 1)then
                    call shadow_findAccObs(iepoch,i,isat,acc,thrusters_off,second_epoch)
                    if(second_epoch /= -999)then
                       y_acc_tmp = apr_ScaleBias(i,isat)*(acc(second_epoch,i+1,isat) - acc(iepoch,i+1,isat))
                       print*, y_acc_tmp, isat
!!!           Scale
                       part(icond+ncond_t*(isat-1)+2,iScale-1+i+nsatprm_t*(isat-1)) = acc(second_epoch,i+1,isat)  !+apr_ScaleBias(i+3,isat))      ! d[equ]/dScl_y  = ACC_y
                       !!           Bias
                       part(icond+ncond_t*(isat-1)+2,iBias-1+i+nsatprm_t*(isat-1)) = 1.d0            ! d[equ]/dBias_y = 1.0
                       !!           set pre-omc
                       pre_omc(icond+(isat-1)*ncond_t+2) = y_acc_tmp - ( apr_ScaleBias(i,isat)*acc(second_epoch,i+1,isat) &
                            + apr_ScaleBias(i+3,isat)) ! Observed is 0.0
                    endif
                 endif

              else
                 write(message,'(a)')'Error: cannot apply shadow scale conditions unless estimating bias and scale'
                 call status_update('FATAL','GRACEFIT','shadow_ApplyCondition',' ',message,0)

              endif

           endif

           if(i == 3 )then    ! apply conditional equations on "z"
              !         Set part
              if(nScale_t > 0 .and. nBias_t > 0) then
                 ! PT140922: "calibrate" the accelerometer value to use in the partial (ie add the bias) so that the scale partial varies with time
                 ! PT150611: no, don't
                 !           ! Scale 
                 part(icond+ncond_t*(isat-1)+1,iScale-1+i+nsatprm_t*(isat-1)) = acc(iepoch,i+1,isat) !  + apr_ScaleBias(i+4,isat)  ! d[equ]/dScl_z  = ACC_z
                 !           ! Bias
                 part(icond+ncond_t*(isat-1)+1,iBias-1+i+nsatprm_t*(isat-1)) = 1.d0            ! d[equ]/dBias_z = 1.0
                 !         Set pre_omc
                 !! PT140801: try reducing the -10 to -5 nm/s^2 (just to see what happens)
                 pre_omc(icond+(isat-1)*ncond_t+1)= -5.d-3 - (apr_ScaleBias(i,isat)*acc(iepoch,i+1,isat)+apr_ScaleBias(i+3,isat)) ! Observed is -10.0

                 ! PT150610: add another conditional observation here, being half a revolution away from the shadow observation. Use the a priori scale value to
                 !           scale the difference between the two obs, which then gives the "calibrated" value of the second (since the shadow ob should be zero)
                 if(int(use_accel_scale_const) == 1)then
                    if(second_epoch /= -999)then
                       z_acc_tmp = apr_ScaleBias(i,isat)*(acc(second_epoch,i+1,isat) - acc(iepoch,i+1,isat)) - 5.d-3  ! use -5nm/s^2 as the value during eclipse

                       !!           Scale
                       part(icond+ncond_t*(isat-1)+3,iScale-1+i+nsatprm_t*(isat-1)) = acc(second_epoch,i+1,isat)   ! +apr_ScaleBias(i+3,isat))      ! d[equ]/dScl_y  = ACC_y
                       !!           Bias
                       part(icond+ncond_t*(isat-1)+3,iBias-1+i+nsatprm_t*(isat-1)) = 1.d0            ! d[equ]/dBias_y = 1.0
                       !!           set pre-omc
                       pre_omc(icond+(isat-1)*ncond_t+3) = z_acc_tmp - ( apr_ScaleBias(i,isat)*acc(second_epoch,i+1,isat) &
                            + apr_ScaleBias(i+3,isat)) ! Observed is 0.0
                    endif
                 endif

              else
                 write(message,'(a)')'Error: cannot apply shadow bias conditions unless estimating bias and scale'
                 call status_update('FATAL','GRACEFIT','shadow_ApplyCondition',' ',message,0)

              endif
           endif
        enddo  ! End of coordinate loop
     endif  ! End of "apply condition" if statment
  enddo  ! End of satellite loop
  !****************************************************************

  return
end subroutine shadow_applyCondition

!********************************************************************************************************************************


!********************************************************************************************************************************
! shadow_applyCondition: Places partials appropriate partials and observed minues calculated
!                        values into part and pre_omc respectively given the array determining
!                        whether to apply the conditions or not, the acceleration of each
!                        satellite, the a priori values and the two arrays mentioned above
! Author: Thomas Greenspan
!********************************************************************************************************************************

subroutine shadow_applyCondition_drag(iepoch,apply_shadow_cond,thrusters_off,apr_ScaleBias,part,pre_omc, lead)

  ! MODS
  ! PT150609: add a new data point, where we use the a priori scale and the difference in y-acc between shadow and point separated by 1/2 revolution
  !           to define a new "known" calibrated point. We need this in order to stop the scale factor going to zero when we tightly constraint these
  !           condition equations.
  ! PT150609: pass the entire acc and apply_shadow_cond arrays into the subroutine, not just one epoch's worth of values.
  ! PT160405: use a different condition equation that accounts for the fact that the atmospheric drag does affect the z-acc obs, and
  !           with a different sign for leading and trailing satellite. Thus, the mean of both satellites has a (slightly) negative
  !           value rather than each satellite having that value.
  ! PT180626: modifications to use the "acc_obs" array rather than the old "acc" array.

  use gracefit_mod
  use accred_mod      ! brings in the accelerometer observations

  implicit none


  !********************  Variable declarations ********************

  integer*4 , intent(in)        :: iepoch                                  ! epoch number
  logical, intent(in)           :: apply_shadow_cond(nepochs_t,nsat_t)     ! Data on whether the shadow condition can be applied
  logical,          intent(out) :: thrusters_off(nepochs_t,nsat_t)         ! Data on whether thrusters are having any effect on accelerations
!  double precision, intent(in)  :: acc(nepochs_t,4,nsat_t)                 ! Accelerations for each sattelite at every epoch
  double precision, intent(in)  :: apr_ScaleBias(max_SB,maxsat)            ! A priori values for parameters
  double precision, intent(out) :: part(nobs_t,nparam_t)                   ! Partials of observables with respect to the parameters
  double precision, intent(out) :: pre_omc(nobs_t)                         ! Observed Minus Computed values to be used in LS algorithm

  double precision, parameter :: ECLIPSE = 0.0  ! Constant indicating the satellite is in shadow
  integer*4 :: isat  ! Counter variable for satellites
  integer*4 :: i     ! Counter variable
  integer*4 :: ahead, behind ! Which Satellite is in front and behind
  character*100 :: message  ! for status_update calls
  integer*4 :: second_epoch ! epoch number of an epoch half a revolution from a shadow epoch
  real(kind=8) :: y_acc_tmp ! difference between two y-acc observations, calibrated by the a priori scale factor
  real(kind=8) :: z_acc_tmp ! difference between two z-acc observations, calibrated by the a priori scale factor
  real(kind=8) :: f1, f2    ! weighting factors for calculating the mean of the calibrated Z-acc observations
  real(kind=8) :: acc_trailing, acc_leading   ! temporary computations of the observed acc_Z values, calibrated using apriori scl/bsz
  character(1) :: lead      ! Which Satellite is ahead

! local acc array, filled out based on the available acc_obs data
  real(kind=8),allocatable :: acc(:,:,:)
  integer*4                :: iacc
  integer*4                :: acc_ptr(2)
  allocate(acc(nepochs_t,4,nsat_t))

  !****************************************************************
  ! NOTE:   The conditions are as follows  
  !  Y - When a satellite is in shadow, the acceleration in the
  !      y-direction is close to 0 since there is no force from
  !      the sun. The conditional equation is:
  !      Scl_y * ACC_y + Bias_y = 0.0
  !      The partials for scale and bias are thus ACC_y and 1
  !      respectively and observed is 0.0
  !  Z - When a satellite is in shadow, the acceleration in the
  !      z-direction is slightly negative since the albedo effect
  !      is (assumed) stronger than the radiation off the top of
  !      the satellite. The conditional equation is:
  !      Scl_y * ACC_y + Bias_y = -10.0
  !      The partials for scale and bias are thus ACC_y and 1
  !      respectively and observed is -10.0 (value may change)
  !    As the top of the satellite cools off, the radiation it
  !      gives off will decrease and the observed acceleration
  !      will decrease (become more negative). This is not modeled
  !      here. We use a constant with a great enough uncertainty.
  !    ->If this condition proves very useful we can refine our
  !      observed value (make time dependent, more accurate, etc...)
  !
  ! PT160405: Z - because the satellites are pitched to align with the
  !               LOS vector, the lead satellite has a -ive bias and
  !               the trailing satellite has a +ive bias because of 
  !               a small component of atmospheric drag in the SRF Z
  !               direction. So we can't expect each of them to be
  !               -10 nm/s^2 during eclipse. Rather, the mean of both
  !               will have this value.
  !****************************************************************


  ! ! ! ! ! ! ************************
!  PT180626: this needs to be fixed so that it is GRACE FO compatible !!!
  ! PT170702: Adding in Differing lead satellites
  if (lead == 'A' .or. lead == 'C') THEN
     ahead = 1
     behind = 2
  else
     ahead = 2
     behind = 1
  endif
  ! ! ! ! ! ! ************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT180626: we need to convert iepoch into a row number within the acc_obs array of accelerometer
! observations. We can do this through knowing 
! a) the start epoch of the accelerometer obs (in gracesec)
! b) the starting_epoch of the orbit integration
! c) the "iepoch" of the epoch entered into this subroutine
  do isat = 1,nsat_t
    acc_ptr(isat) = nint(starting_epoch - acc_obs(isat,1,1))+iepoch ! +1 ?
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PT150612: the conditional shadow equations are as follows:
  !           row icond  : yacc on GRACE A during shadow
  !               icond+1: zacc on GRACE A during shadow
  !               icond+2: yacc on GRACE A 1/2 revolution from shadow obs
  !               icond+3: zacc on GRACE A 1/2 revolution from shadow obs
  !
  !               icond+ncond_t  : yacc on GRACE B during shadow
  !               icond+ncond_t+1: zacc on GRACE B during shadow
  !               icond+ncond_t+2: yacc on GRACE B 1/2 revolution from shadow obs
  !               icond+ncond_t+3: zacc on GRACE B 1/2 revolution from shadow obs

  ! PT160405: reorganise this. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Apply first for the Y axis   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do isat = 1, nsat_t
     if(apply_shadow_cond(iepoch,isat))then ! Make sure it is ok to apply condition
!print*,'Apply shadow condition isat,iepoch',isat,iepoch
        !       apply conditional equations on "y"
        i=2
	!print*, nScale_t, nBias_t
        if(nScale_t > 0 .and. nBias_t > 0)then
           ! first, the accelerometer observation when the satellite is in shadow
           !         Scale
           part(icond+ncond_t*(isat-1),iScale-1+i+nsatprm_t*(isat-1)) =  acc_obs(isat,acc_ptr(isat),i+1) ! + apr_ScaleBias(i+3,isat))   ! d[equ]/dScl_y  = ACC_y
           !         bias
           part(icond+ncond_t*(isat-1),iBias-1+i+nsatprm_t*(isat-1)) = -1.d0            ! d[equ]/dBias_y = 1.0
           !         Set pre_omc
           pre_omc(icond+(isat-1)*ncond_t) = 0.d0 - ( apr_ScaleBias(i,isat)*acc_obs(isat,acc_ptr(isat),i+1) + apr_ScaleBias(i+3,isat)) ! Observed is 0.0

! DEBUG
!print*,"applyCondition: iepoch,acc_ptr,acc_obs",isat,iepoch,acc_ptr(isat),acc_obs(isat,acc_ptr(isat),i+1)

        else
           write(message,'(a)')'Error: cannot apply shadow scale conditions unless estimating bias and scale'
           call status_update('FATAL','GRACEFIT','shadow_ApplyCondition',' ',message,0)

        endif

     endif
  enddo  ! end of satellite loop applying the condition on Y-acc observations

  ! PT160405: Apply the condition to the mean value of Z-acc, but only when applicable obs exist for both satellites.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   Apply next for the Z axis   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! our conditional observation equation is that the (perhaps weighted?) mean value is slightly negative
  !
  !     -5 =  f1*(Scl_Za*ACC_Za + Bias_Za) + f2*(Scl_Zb*ACC_Zb + Bias_Zb)] 
  !
  ! to begin with, make f1 = f2 = 0.5  We can perhaps improve this later to introduce asymmetry should it be realistic (it will 
  ! depend on the pitch of each satellite wrt the tangent of the satellite orbit.

  if(apply_shadow_cond(iepoch,1).and.apply_shadow_cond(iepoch,2))then ! Make sure it is ok to apply condition
     i = 3    ! apply conditional equations on "z"
     ! PT160405: we may want to come back and define the proportion of drag affecting each satellite (i.e. does one tilt away
     !           from nominal more than the other?). For now, just assign them as being weighted 50% each.
     f1 = 0.5
     f2 = 1.0 - f1

     if(nScale_t > 0 .and. nBias_t > 0)then
        !     ! Scale Lead 
        part(icond+1,iScale-1+i) = f2*acc_obs(ahead,acc_ptr(ahead),i+1)               !! d[equ]/dScl_z  = f1*ACC_z
        !     ! Bias  Lead
        part(icond+1,iBias-1+i) = f2*1.d0                                         ! d[equ]/dBias_z = f1*1.0
        !     ! Scale Tail
        part(icond+1,iScale-1+i+nsatprm_t) = f1*acc_obs(behind,acc_ptr(behind),i+1)   !! d[equ]/dScl_z  = f2*ACC_z
        !     ! Bias  Tail
        part(icond+1,iBias-1+i+nsatprm_t) = f1*1.d0                               ! d[equ]/dBias_z = f2*1.0
        !     ! Set pre_omc
        ! PT140801: try reducing the -10 to -5 nm/s^2 (just to see what happens)
        pre_omc(icond+1)= -5.d-3 - (f1*(apr_ScaleBias(i,1)*acc_obs(1,acc_ptr(1),i+1)+apr_ScaleBias(i+3,1)) &
             + f2*(apr_ScaleBias(i,2)*acc_obs(2,acc_ptr(2),i+1)+apr_ScaleBias(i+3,2))) ! Observed is -5.0 nm/s^2

        ! PT160420: the condition making the mean of the accZ obs of both satellites during shadow equal to -5 works, but it 
        !           doesn't impose that the trailing satellite should have calibrated values that are > leading satellite. So we
        !           need to add another conditional equation:
        !
        !    accZ_trailing = accZ_leading + abs(d_acc_Z)     ! ie that the trailing satellite value is greater
        !
        !           to implement this, we need a new condition equation. Since we have reduced the acc_Z down to one equation for two satellites
        !           (it had been one equation per satellite) we can make use of the GRACE B acc_Z row to add the new equation.
        !
        !               0  =  (acc_Z_obs_trailing * sclz_trailing + bsz_trailing) - (acc_Z_obs_leading * sclz_leading + bsz_leading) - (acc_Z_obs_leading - acc_Z_obs_trailing)
        !
        !  which has partials of 
        !
        !               d/dsclZ_trailing =  acc_Z_obs_trailing ,   d/dbsZ_trailing =  1
        !               d/dsclZ_leading  = -acc_Z_obs_leading  ,   d/dbsZ_leading  = -1

        ! PT160421: GRACE A trails from Dec 2005 to into 2011 ..... hardwire it accordingly at present but need to make the code better here !!!
        !     ! Scale Lead
        part(icond+ncond_t+1,iScale-1+i) = 1.d0*acc_obs(ahead,acc_ptr(ahead),i+1)                  !! d[equ]/dScl_z  = f1*ACC_z
        !     ! Bias  Lead
        part(icond+ncond_t+1,iBias-1+i) = 1.d0                                         ! d[equ]/dBias_z = f1*1.0
        !     ! Scale Tail
        part(icond+ncond_t+1,iScale-1+i+nsatprm_t) = -1.*acc_obs(behind,acc_ptr(behind),i+1)        !! d[equ]/dScl_z  = f2*ACC_z
        !     ! Bias  Tail
        part(icond+ncond_t+1,iBias-1+i+nsatprm_t) = -1.d0                              ! d[equ]/dBias_z = f2*1.0
        !     ! Set pre_omc
        acc_trailing = apr_ScaleBias(i,behind)*acc_obs(behind,acc_ptr(behind),i+1)+apr_ScaleBias(i+3,behind)
        acc_leading  = apr_ScaleBias(i,ahead)*acc_obs(ahead,acc_ptr(ahead),i+1)+apr_ScaleBias(i+3,ahead)
        pre_omc(icond+ncond_t+1)= 0 - ( acc_trailing - acc_leading - abs(acc_trailing - acc_leading))

     else
        write(message,'(a)')'Error: cannot apply Z shadow bias conditions unless estimating bias and scale'
        call status_update('FATAL','GRACEFIT','shadow_ApplyCondition',' ',message,0)
     endif
  endif   ! end of check whether it is appropriate to apply shadow condition
  !****************************************************************

  return
end subroutine shadow_applyCondition_drag

!********************************************************************************************************************************


!********************************************************************************************************************************
! shadow_SRP: Determines when each satellite goes through the penumbra of the Earth, which is when the solar radiation pressure
!             force is either applied or removed. We can use this information to relate the scale of the accelerometers in each
!             direction (with knowledge of the satellite orientation, surface area and a little bit of geometry ...)
! 
! Author: Paul Tregoning
!         19 May 2015
!        
!********************************************************************************************************************************

subroutine shadow_SRP(thrusters_off,acc,rvec,srf2trf_rotmat,normeq,AtWb,apr_ScaleBias)

  use gracefit_mod
  implicit none


  logical, intent(in)            :: thrusters_off(nepochs_t,nsat_t)      ! Data on whether thrusters are having any effect on accelerations
  double precision, intent(in)   :: acc(nepochs_t,4,nsat_t)              ! Accelerations for each sattelite at every epoch
  real(kind=8), intent(in)       :: rvec(nrvec_t,nsat_t,nepochs_t)       ! Earth-fixed satellite positions
  real(kind=8), intent(in)       :: srf2trf_rotmat(3,3,nepochs_t,nsat_t) ! rotation matrices from SRF to TRF
  real(kind=8), intent(inout)    :: normeq(nparam_t,nparam_t)            ! the normal equations
  real(kind=8), intent(inout)    :: AtWb(nparam_t)                       ! the RHS
  real(kind=8), intent(in)       :: apr_ScaleBias(max_SB,maxsat)         ! A priori values for scale and bias parameters

  ! local variables
  integer*4 :: iepoch, counter, i,j,isat
  integer*4 :: change_shad                             ! counter of how many epochs since the shadow derivative has been non-zero
  integer*4 :: n_transition                            ! number of usable epochs (unaffected by thrusts) during a shadow transition 
  real(kind=8) :: clean_acc(nepochs_t,7)
  real(kind=8) :: shadow(nepochs_t,2)
  real(kind=8) :: mean_value(3)
  real(kind=8) :: sunpos_efixed(6)
  real(kind=8) :: amag3,dot
  real(kind=8) :: sun_to_sat(3,nsat_t)
  real(kind=8) :: ang_X(nepochs_t), ang_Y(nepochs_t), ang_Z(nepochs_t), sat_X_trf(3),sat_Y_trf(3),sat_Z_trf(3)
  real(kind=8) :: dx_acc,dy_acc,dz_acc
  real(kind=8) :: sigma_sclYZ,sigma_sclYX,sigma_sclXZ,YZ_ratio,YX_ratio,XZ_ratio
  integer*4    :: param_x,param_y,param_z              ! parameter numbers for the accelerometer scale parameters
  real(kind=8) :: tmp_omc(2)                           ! obs - comp for the accel scale conditional observations

  call status_update('STATUS','GRACEFIT','shadow_SRP',' ',"Applying accelerometer scale constraints",0)

  ! loop over satellites
  do isat = 1,nsat_t

     ! we need to make an array that contains only the accelerometer values 
     ! Y-axis accelerometer values
     counter = 0
     do iepoch=1,nepochs_t
        if(thrusters_off(iepoch,isat) )then
           counter = counter + 1
           clean_acc(counter,1) = iepoch * 5.0
           call shadow_getsun(starting_epoch,iepoch,-1,sunpos_efixed) ! sunpos in units of m (I think!!) and is sun relative to earth in earth-fixed coords
           call shadow_isInShadow(sunpos_efixed(1:3),rvec(1:3, isat, iepoch),shadow(counter,1))

           do j=1,3
              clean_acc(counter,j+1) = acc(iepoch,j+1,isat)
           enddo
        endif
     enddo

     ! remove the mean value
     do j=1,3
        mean_value(j) = sum(clean_acc(1:counter,j+1)) / dble(counter)
        clean_acc(1:counter,j+1) = clean_acc(1:counter,j+1) - mean_value(j)
     enddo

     ! now calculate the derivative of the curve. Here we hope that the edges of the entry/exit to shadow stick out !!!
     do j=1,3
        call noise_robust_deriv(shadow(:,1), shadow(:,2), epoch_interval, counter, 5)
     enddo

     change_shad = 0
     do i=1,counter
        iepoch = abs(clean_acc(i,1))/5.0
        ! PT150525: if the derivative of the shadow factor is not zero then the satellite is going into or out of eclipse
        if(abs(shadow(i,2)) > 0.d0)then
           change_shad = change_shad + 1
           ! get the earth-fixed sun coordinates
           call shadow_getsun(starting_epoch,iepoch,-1,sunpos_efixed)  ! sunpos (in metres)  and is returned as sun relative to earth in earth-fixed coords

           ! compute the sun-to-satellite unit vector. Both sunpos and rvec are e-fixed (I think) so just difference them.
           sun_to_sat(:,isat) = rvec(1:3, isat, iepoch) - sunpos_efixed(1:3)
           sun_to_sat(:,isat) = sun_to_sat(:,isat)/amag3(sun_to_sat(:,isat))
           ! extract the unit vectors of the SRF X,Y,Z directions from the full matrix of values
           sat_X_trf(1:3) = srf2trf_rotmat(:,1,iepoch,isat)
           sat_Y_trf(1:3) = srf2trf_rotmat(:,2,iepoch,isat)
           sat_Z_trf(1:3) = srf2trf_rotmat(:,3,iepoch,isat)

           ! compute the angle between the Z-axis and the sun-sat vector
           ang_Z(i) = dot(sun_to_sat(:,isat),sat_Z_trf)

           ! compute the angle between the Y-axis and the sun-sat vector
           ang_Y(i) = dot(sun_to_sat(:,isat),sat_Y_trf)

           ! compute the angle between the X-axis and the sun-sat vector
           ang_X(i) = dot(sun_to_sat(:,isat),sat_X_trf)
        elseif(change_shad > 0)then
           ! satellite has finished changing through the Earth's penumbra. We want to use the past few obs to get our condition equation values
           ! number of points in transition
           n_transition = change_shad

           ! check that there were continuous obs through the transition
           if (i-n_transition-1 > 0 .and. n_transition > 5)then
              if( (clean_acc(i-1,1) - clean_acc(i-n_transition-1,1) ) > n_transition * epoch_interval )then
                 !              print*,'cannot use this shadow transition: epochs missing'
              else
                 !            print*,'have a good transition here',clean_acc(i-1,1),clean_acc(i-n_transition,1),n_transition
                 ! calculate the difference in the uncalibrated accelerometer observations in the three orthogonal directions
                 dx_acc =    (clean_acc(i-1,2) - clean_acc(i-n_transition-1,2) )           
                 dy_acc =    (clean_acc(i-1,3) - clean_acc(i-n_transition-1,3) )          
                 dz_acc =    (clean_acc(i-1,4) - clean_acc(i-n_transition-1,4) ) 

                 !print*,"Xacc values",clean_acc(i-1,2)*1.d3 ,clean_acc(i-n_transition-1,2)*1.d3,dx_acc*1.d3
                 !print*,"Yacc values",clean_acc(i-1,3)*1.d3 ,clean_acc(i-n_transition-1,3)*1.d3,dy_acc*1.d3
                 !print*,"Zacc values",clean_acc(i-1,4)*1.d3 ,clean_acc(i-n_transition-1,4)*1.d3,dz_acc*1.d3


                 ! see if I can make some estimate of the uncertainty of this obervation, based on the accuracy of picking the dx_acc, dy_acc, dz_acc
                 ! Let the y-scale = 1.0 and the uncertainty of each accelerometer obs = 1 nm/s^2 (therefore sigma of 1.414 nm/s^2 for each d?_acc)
                 ! scale Z

                 sigma_sclYZ =  sqrt(  sigma_acc**2*(ang_Z(i-1)/(dz_acc*ang_Y(i-1)))**2  &
                      + sigma_acc**2*(  dy_acc*ang_Z(i-1)*(-dz_acc)/(ang_Y(i-1)*(dz_acc**2)) )**2) 

                 ! scale X
                 sigma_sclYX =  sqrt(  sigma_acc**2*(ang_X(i-1)/(dx_acc*ang_Y(i-1)))**2  &
                      + sigma_acc**2*(  dy_acc*ang_X(i-1)*(-dx_acc)/(ang_Y(i-1)*(dx_acc**2)) )**2) 

                 ! scale Z wrt X
                 sigma_sclXZ =  sqrt(  sigma_acc**2*(ang_Z(i-1)/(dz_acc*ang_X(i-1)))**2  &
                      + sigma_acc**2*(  dx_acc*ang_Z(i-1)*(-dz_acc)/(ang_X(i-1)*(dz_acc**2)) )**2) 

                 ! now calculate the relation between sclY and the other two scales
                 YZ_ratio = dy_acc/dz_acc * ang_Z(i-1)/ang_Y(i-1)
                 YX_ratio = dy_acc/dx_acc * ang_X(i-1)/ang_Y(i-1)
                 XZ_ratio = dx_acc/dz_acc * ang_Z(i-1)/ang_X(i-1)

                 write(*,'(2i6,3(a,f9.3,a,e10.3))') isat,iepoch,' YZ ratio',YZ_ratio,' +/-',sigma_sclYZ,'    YX ratio',YX_ratio,' +/-' &
                      ,sigma_sclYX ,'    XZ ratio',XZ_ratio,' +/-',sigma_sclXZ


                 ! PT150527: now, increment these condition equations into the normal equations
                 param_x = iScale+nsatprm_t*(isat-1)
                 param_y = iScale+nsatprm_t*(isat-1) + 1
                 param_z = iScale+nsatprm_t*(isat-1) + 2
                 !             relation of sclZ to sclY
                 normeq(param_z,param_z) = normeq(param_z,param_z) + 1.0/(sigma_sclYZ)**2
                 normeq(param_y,param_z) = normeq(param_y,param_z) - YZ_ratio/(sigma_sclYZ)**2
                 normeq(param_z,param_y) = normeq(param_y,param_z)
                 normeq(param_y,param_y) = normeq(param_y,param_y) + YZ_ratio**2/(sigma_sclYZ)**2
                 !             relation of sclZ to sclX
                 normeq(param_z,param_z) = normeq(param_z,param_z) + 1.0/(sigma_sclXZ)**2
                 normeq(param_x,param_z) = normeq(param_x,param_z) - XZ_ratio/(sigma_sclXZ)**2
                 normeq(param_z,param_x) = normeq(param_x,param_z)
                 normeq(param_x,param_x) = normeq(param_x,param_x) + XZ_ratio**2/(sigma_sclXZ)**2

                 ! PT150527: and now into the RHS
                 tmp_omc(1) = 0.0 - (apr_ScaleBias(3,isat) - apr_ScaleBias(2,isat)*YZ_ratio)
                 tmp_omc(2) = 0.0 - (apr_ScaleBias(3,isat) - apr_ScaleBias(1,isat)*XZ_ratio)

                 ! PT150602: these are wrong: I need to divide by sigma^2. Apply the fix after the tests of 2 June 2015 will have finished.
                 !             scale X
                 AtWb(param_x) = AtWb(param_x) - (XZ_ratio * tmp_omc(2))/(sigma_sclXZ)**2
                 !             scale Y
                 AtWb(param_y) = AtWb(param_y) - (YZ_ratio * tmp_omc(1))/(sigma_sclYZ)**2
                 !             scale Z
                 AtWb(param_z) = AtWb(param_z) + tmp_omc(2)/(sigma_sclXZ)**2 + tmp_omc(1)/(sigma_sclYZ)**2
              endif
              change_shad = 0
           endif
        else
           !  satellite is either in eclipse or not in eclipse, but it is not going through the transition from one to the other
        endif
     enddo

  enddo

  ! stop 'stopped in shadow_SRP'




  !****************************************************************

  return
end subroutine shadow_SRP

!********************************************************************************************************************************

!********************************************************************************************************************************
! shadow_getsun: Calculates the sun vector (relative to the Earth) in earth-fixed coordinates
!
! Author: Paul Tregoning
!         21 May 2015
!********************************************************************************************************************************

subroutine shadow_getsun(start_grace_seconds,iepoch,idir,sunpos)
  use gracefit_mod
  implicit none

  include 'input.h90'

  real(kind=8),intent(in)     :: start_grace_seconds    ! GRACE seconds of the first epoch of the orbit
  integer*4   ,intent(in)     :: iepoch
  integer*4   ,intent(in)     :: idir
  real(kind=8),intent(out)    :: sunpos(6)

  integer*4        :: date(5)          ! Varaible needed to calculate the julian date
  double precision :: sec              ! Variable needed to calculate the julian date
  double precision :: fjd              ! Julian date
  double precision :: fjd_dec          ! Integer and decimal parts of fjd
  double precision :: ephtime(2)       ! Integer part of epoch and fraction of a day
  double precision :: sunpos_efixed(3) ! Position of the sun in efixed frame
  double precision :: PEPt                 
  integer*4        :: PEPjd
  integer*4        :: iut1pol       ! Values to be input to rotsnp
  double precision :: rot(3,3)      ! Rotation matrix to be used in satellite position calculation
  double precision :: rotderiv(3,3) ! Derivative of above matrix (Not used here but needed as input for rotsnp)
  double precision :: sidtm(3,3)    ! Sidereal time (rad) (Not used here but needed as input for rotsnp)


  ! Compute Julian date
  call gsec_to_ymdhms (start_grace_seconds, date, sec)
  call ymdhms_to_jd(date,sec,fjd)

  ! Set Values
  fjd_dec = fjd - int(fjd)                  ! decimal part of fjd
  ephtime(1) = dble(int(fjd))               ! integer part of fjd
  ephtime(2) = fjd_dec + ((iepoch-1)*epoch_interval+ tdtoff)/SEC_IN_DAY  ! fraction of a day   (tdtoff included as per graceorb/planetfield.f90)
  !****************************************************************

  ! Set sunpos_i
  call dpleph (ephtime, 11, 3, sunpos)
  sunpos = 1.d0*sunpos*1.d3    ! puts sunpos in units of m (I think!!) and is sun relative to earth in inertial coordinates
  ! sunpos from the JPL ephemeris are coordinates in the intertial frame. Rotate intertial to e-fixed for sunpos if so requested.
  if(idir == 1) then
     return  ! we wanted - and we already have - sunpos in inertial coordinates !
  else
     call status_update('FATAL','GRACEFIT','shadow_getSun',' ',"Requests for Sun coords in earth-fixed no longer supported",0)
!     ! need to convert sunpos into earth-fixed coordinates
!     iut1pol = 7  ! The value used in graceorb
!
!     ! Get rotation matrix between inert and efixed
!     ! PT140731: include (iepoch-1)*5 seconds of day into the time epoch here
!     call PEPtime(int(fjd),fjd_dec*SEC_IN_DAY+(iepoch-1)*epoch_interval, PEPjd, PEPt)
!     call rotsnp(idir,PEPjd,PEPt,tdtoff,iut1pol,rot,rotderiv,sidtm,iUT1,iPOLE,iNUT,'J2000','IAU76')  ! PEPtime test
!
!     ! convert sunpos to earth-fixed coordinates
!     call matmult(rot,sunpos(1:3),sunpos_efixed(1:3),3,3,1)
!     sunpos(1:3) = sunpos_efixed(1:3)

  endif



  !****************************************************************
  return

end subroutine shadow_getsun
!********************************************************************************************************************************


!********************************************************************************************************************************
! shadow_findAccObs: Find an accelerometer observation that is ~half an orbit away from another one and that is not affected by thrusts.
!                    Returns either an epoch number or -999 if no appropriate epoch was found.
!
! Author: Paul Tregoning
!         09 June 2015
!********************************************************************************************************************************

subroutine shadow_findAccObs(iepoch,iaxis,isat,acc,thrusters_off,second_epoch)
  use gracefit_mod
  implicit none


  integer*4, intent(in)     :: iepoch                            ! current epoch
  integer*4, intent(in)     :: iaxis                             ! which acc observation are we looking for (2: y-acc, 3: z-acc)
  integer*4, intent(in)     :: isat                              ! satellite number
  real(kind=8), intent(in)  :: acc(nepochs_t,4,nsat_t)           ! raw accelerometer observations
  logical, intent(in)       :: thrusters_off(nepochs_t,nsat_t)   ! Data on whether thrusters are having any effect on accelerations
  integer*4, intent(out)    :: second_epoch                      ! epoch of the other accelerometer observation half a revolution away

  ! local variables
  integer*4  :: half_rev     ! the number of epochs in half a revolution 
  integer*4  :: idir         ! forward search or backward search, depending on whether we are close to the start/end of the orbit integration
  logical found
  integer*4  :: tmp_epoch
  character*100 :: message

  half_rev = 558    ! (93/2 mins * 12 epochs/min = 558 epochs)


  ! do we go forwards or backwards
  if (iepoch < 558) then 
     idir = 1
  else
     idir = -1
  endif

  ! start searching
  found = .false.
  tmp_epoch = 0
  second_epoch = -999
  do while (.not. found .and. tmp_epoch < half_rev/2 .and.iepoch+idir*half_rev+idir*tmp_epoch > 0 )
     if(thrusters_off(iepoch+idir*half_rev+idir*tmp_epoch,isat))then
        second_epoch = iepoch+idir*half_rev+tmp_epoch
        found = .true.
        ! DEBUG
        !print*,'iepoch,second_epoch,idir,tmp_epoch,acc',iepoch,second_epoch,idir,tmp_epoch,acc(iepoch,iaxis+1,isat),isat,iaxis
     else
        tmp_epoch = tmp_epoch + 1
     endif
  enddo

  if(second_epoch == -999) then
     write(message,'(a,i6)')"No available non-thrust accel observation found ~1/2 revolution from epoch: ",iepoch
     call status_update('WARNING','GRACEFIT','shadow_findAccObs',' ',message,0)
  endif
  ! that should be all we need to do
  !****************************************************************
  return

end subroutine shadow_findAccObs
!********************************************************************************************************************************


!********************************************************************************************************************************
! shadow_SRPdotprod: Calculate the dot products of the Sun-Sat vector into the along-track, cross-track and radial directions of
!                      the accelerometers. Returns the three angles
!
! Author: Paul Tregoning
!         16 October 2015
!********************************************************************************************************************************
subroutine shadow_SRPdotprod(thrusters_off,acc,rvec,srf2trf_rotmat,SRPdotprod)

  use gracefit_mod
  implicit none


  logical, intent(in)            :: thrusters_off(nepochs_t,nsat_t)      ! Data on whether thrusters are having any effect on accelerations
  double precision, intent(in)   :: acc(nepochs_t,4,nsat_t)              ! Accelerations for each sattelite at every epoch
  real(kind=8), intent(in)       :: rvec(nrvec_t,nsat_t,nepochs_t)       ! Earth-fixed satellite positions
  real(kind=8), intent(in)       :: srf2trf_rotmat(3,3,nepochs_t,nsat_t) ! rotation matrices from SRF to TRF
  real(kind=8), intent(out)      :: SRPdotprod(3)                        ! the three dot products, one for each axis

  ! local variables
  integer*4 :: iepoch, counter, i,j,isat
  real(kind=8) :: clean_acc(nepochs_t,7)
  real(kind=8) :: shadow(nepochs_t,2)
  real(kind=8) :: mean_value(3)
  real(kind=8) :: sunpos_efixed(6)
  real(kind=8) :: amag3,dot
  real(kind=8) :: sun_to_sat(3,nsat_t)
  real(kind=8) :: tmpshadow, sat_X_trf(3),sat_Y_trf(3),sat_Z_trf(3)

  call status_update('STATUS','GRACEFIT','shadow_SRP',' ',"Breaking SRP into three orthogonal components for GRACE A",0)
  isat = 1  ! only GRACE A at this stage while I'm testing whether this even works!!!

  ! loop through all the epochs
  do iepoch=1,nepochs_t
     ! get the earth-fixed sun coordinates
     call shadow_getsun(starting_epoch,iepoch,-1,sunpos_efixed)  ! sunpos (in metres)  and is returned as sun relative to earth in earth-fixed coords
     call shadow_isInShadow(sunpos_efixed(1:3),rvec(1:3, isat, iepoch),tmpshadow)

     ! compute the sun-to-satellite unit vector. Both sunpos and rvec are e-fixed (I think) so just difference them.
     sun_to_sat(:,isat) = rvec(1:3, isat, iepoch) - sunpos_efixed(1:3)
     sun_to_sat(:,isat) = sun_to_sat(:,isat)/amag3(sun_to_sat(:,isat))
     ! extract the unit vectors of the SRF X,Y,Z directions from the full matrix of values
     sat_X_trf(1:3) = srf2trf_rotmat(:,1,iepoch,isat)
     sat_Y_trf(1:3) = srf2trf_rotmat(:,2,iepoch,isat)
     sat_Z_trf(1:3) = srf2trf_rotmat(:,3,iepoch,isat)

     ! compute the angle between the X-axis and the sun-sat vector
     SRPdotprod(1) = acos(dot(sun_to_sat(:,isat),sat_X_trf))

     ! compute the angle between the Y-axis and the sun-sat vector
     SRPdotprod(2) = acos(dot(sun_to_sat(:,isat),sat_Y_trf))

     ! compute the angle between the Z-axis and the sun-sat vector
     SRPdotprod(3) = dot(sun_to_sat(:,isat),sat_Z_trf)

  enddo

  return
end subroutine shadow_SRPdotprod

!********************************************************************************************************************************

!********************************************************************************************************************************
! leadsat: Finds the satellite in the lead and returns it
!
! Author: Erin Gray
!         2 August 2017
!********************************************************************************************************************************

subroutine leadsat(mission,gvec, nparams,nepochs_t, nsat_t, lead)

  implicit none

  !********************  Variable declarations ********************
  integer*4, intent(in)   :: mission                       ! GRACE: 0 , GRACE FO: 1, GRACE II: 2
  integer*4, intent(in)   :: nepochs_t                     ! number of epochs
  integer*4, intent(in)   :: nparams                       ! number of parameters? I think it means pos/vel ... 6
  integer*4, intent(in)   :: nsat_t                        ! number of satellites
  real(kind=8),intent(in) :: gvec(nparams,nsat_t,nepochs_t)! Array from GNV1B files
  character*1,intent(out) :: lead                          ! satellite in the lead

  integer*4     :: iepoch, a, b  !Counters for A, B and epoch number
  
  !****************************************************************
  !********************     INITIAL VALUES     ********************
  a = 0
  b = 0
  !****************************************************************
  do iepoch = 1, nepochs_t                                                   ! Run Through All Epochs
     if (gvec(6,1, iepoch) > 0 .and. gvec(6,2, iepoch) > 0) THEN       ! If Heading North
        if(gvec(3, 1, iepoch) > gvec(3, 2, iepoch)) THEN               
           a = a + 1
        else if (gvec(3, 1, iepoch) < gvec(3, 2, iepoch)) THEN         
           b = b + 1
        endif
     else if (gvec(6,1, iepoch) < 0 .and. gvec(6,2, iepoch) < 0) THEN  ! If Heading South
        if(gvec(3, 1, iepoch) > gvec(3, 2, iepoch)) THEN
           b = b + 1
        else if(gvec(3, 1, nepochs_t) < gvec(3, 2, iepoch)) THEN
           a = a + 1
        endif
     endif
  enddo

  !print*, 'A, B', a, b

  if (a > b) THEN                                                           ! Set Lead to the correct satellite
     if(mission == 0)lead = 'A'
     if(mission == 1)lead = 'C'
  else if (b > a) THEN
     if(mission == 0)lead = 'B'
     if(mission == 1)lead = 'D'
  endif

  return
end subroutine leadsat









































