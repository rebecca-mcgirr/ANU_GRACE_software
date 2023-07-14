subroutine calc_beta_angle(iepoch,rot_i2e,satpos,satvel,beta_angle)

!  subroutine to calculate the beta angle, being the elevation of the sun from the orbital plane. Both the sunpos and satpos
!  variables are assumed to be in the same reference frame (which, in the case of gracefit, is the inertial frame)
!
!  P. Tregoning
!  21 June 2013
!
! MODS
! PT200707: passed in the rot_i2e matrix and removed all information related to rotsnp

  use gracefit_mod
  implicit none

   
  include 'input.h90'

  !********************  Variable declarations ********************

  integer*4       , intent(in)  :: iepoch          ! epoch number
  real(kind=8)    , intent(in)  :: rot_i2e(3,3)    ! inertial to efixed rotation matrix for this epoch
  double precision, intent(in)  :: satpos(3),satvel(3) ! Vectors lying in the orbital plane
  double precision, intent(out) :: beta_angle      ! Beta angle for this epoch

  double precision :: sunpos(6)     ! Position of the sun in inertial frame
  double precision :: sunpos_efixed(6)  ! Position of the sun in earth-fixed frame
!  double precision :: satpos(3)     ! position of the satellite
!  double precision :: satvel(3)     ! Velocity of the satellite
  double precision :: orb_norm(3)   ! Vector normal to position-velocity plane
  double precision :: dot           ! Function that returns the dot product of vectors of length 3 or less
  double precision :: amag3         ! Function that returns the magnigtude of a vector of length 3
  integer*4        :: date(5)       ! Varaible needed to calculate the julian date
  double precision :: sec           ! Variable needed to calculate the julian date
  double precision :: fjd           ! Julian date
  double precision :: fjd_dec       ! Integer and decimal parts of fjd
  double precision :: ephtime(2)    ! Integer part of epoch and fraction of a day
  character(128)   :: message       ! Message to be written out
!****************************************************************

! Compute Julian date
  call gsec_to_ymdhms (starting_epoch, date, sec)
  call ymdhms_to_jd(date,sec,fjd)
! Set Values
  fjd_dec = fjd - int(fjd)                                   ! decimal part of fjd
  ephtime(1) = dble(int(fjd))                                ! integer part of fjd
  ephtime(2) = fjd_dec + ((iepoch-1)*5 + tdtoff)/SEC_IN_DAY  ! fraction of a day   (tdtoff included as per graceorb/planetfield.f90)

! Set sunpos_i
  call dpleph (ephtime, 11, 3, sunpos)    ! this returns the position of the sun relative to the Earth, in km (I think)
  sunpos = 1.d0*sunpos*1.d3               ! puts sunpos in units of m (I think!!) 

! sunpos from the JPL ephemeris are coordinates in the inertial frame, but rvec is e-fixed frame. Rotate e-fixed to intertial for satpos and satvel
! PT200707: in fact, I will rotate the inertial Sun position to earth-fixed
!  call matmult(rot_i2e,vec1,satpos,3,3,1)
!  call matmult(rot_i2e,vec2,satvel,3,3,1)
  call matmult(rot_i2e,sunpos(1:3),sunpos_efixed(1:3),3,3,1)
!****************************************************************

! calculate the normal to the orbital plane (cross product of position and velocity vectors)
  call cross(satpos,satvel,orb_norm)

! calculate the dot product of this and the earth-to-sun vector
  beta_angle = pi/2.d0 - dacos( dot(orb_norm,sunpos_efixed(1:3))/(amag3(sunpos_efixed(1:3))*amag3(orb_norm)) )
!****************************************************************

!    write(message,'(a,f8.2,a)')"Beta angle = ",beta_angle*rad_deg," degrees."
!    call status_update('STATUS','GRACEFIT','gracefit/gracefit',' ',message,0)

  return
end subroutine calc_beta_angle

