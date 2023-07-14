  subroutine calc_beta_angle(sunpos,satpos,satvel,beta_angle)

!  subroutine to calculate the beta angle, being the elevation of the sun from the orbital plane. Both the sunpos and satpos
!  variables are assumed to be in the same reference frame (which, in the case of gracefit, is the inertial frame)
!
!  P. Tregoning
!  21 June 2013

  implicit none

  double precision sunpos(3),satpos(3),satvel(3),orb_norm(3),beta_angle
  double precision dot,amag3,pi

  pi = 4.d0*datan(1.d0)

! sunpos as passed in is the earth wrt sun. We want the other sign.
  sunpos = -sunpos

! calculate the normal to the orbital plane (cross product of position and velocity vectors)
  call cross(satpos,satvel,orb_norm)
!  print*,'amag3(sunpos),amag3(satpos),amag3(orb_norm)',amag3(sunpos),amag3(satpos),amag3(orb_norm)

! calculate the dot product of this and the earth-to-sun vector
  beta_angle = pi/2.d0 - dacos( dot(orb_norm,sunpos)/(amag3(sunpos)*amag3(orb_norm)) )

  return
  end

