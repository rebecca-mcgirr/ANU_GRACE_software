    subroutine generalrel(jd,t,J2, pos, vel, sunpos, angvel, grelacc)

! Emma-Kate Potter, 16 Sept 2010
! This program calculates the general relativistic contributions to acceleration

! The relativistic correction to the acceleration of an artificial Earth satellite taken from
! IERS conventions. 
! NOTE that the equation given in the IERS conventions 2003 documents has an error and the 
! updated draft documents (2010) are used instead.
!
! MODS
! 
! PT160201: fixed bug in a1 computation of J2 effect, plus error in equation 4.2.25 for a3 effect

    use gm_mod

     implicit none

     integer :: i
     real(kind=8) :: r, v, c, beta, gam, rs , Ie 
     real(kind=8) :: vdotv, rdotv, rdotJ 
     real(kind=8), dimension(3) ::  pos, vel, angvel, J
     real(kind=8), dimension(3) ::  posS, velS
     real(kind=8), dimension(3) ::  rcrossv, vcrossJ, bigcross
     real(kind=8), dimension(3) ::  VScrossRS, RScrossv
     real(kind=8), dimension(3) ::  grelacc
     real(kind=8), dimension(6) ::  sunpos
! PT140205: general relativity accelerations due to J2
     real(kind=8), dimension(3) ::  a1,a3
     real(kind=8) :: J2,a2
     real*8 factor

! PT140407: pass in the epoch to use in debug
     integer*4 jd
     real*8  t

!    CALCULATE distance from earth to satellite
     r = dsqrt(pos(1)*pos(1)+pos(2)*pos(2)+pos(3)*pos(3))
     v = dsqrt(vel(1)*vel(1)+vel(2)*vel(2)+vel(3)*vel(3))

!    SET speed of light, c
     c = 299792458.d0

!    CALCULATE Moment of inertia of earth around polar axis, Ie
     Ie = 8.0365d37

!    ASSIGN VALUES to relativity constants, beta and gamma
     beta = 1.d0
     gam  = 1.d0

!    CALCULATE dot product of satellite velocity with itself (vdotv)
     vdotv = vel(1)*vel(1)+vel(2)*vel(2)+vel(3)*vel(3)
!    CALCULATE dot product of satellite position with velocity (rdotv)
     rdotv = pos(1)*vel(1)+pos(2)*vel(2)+pos(3)*vel(3)

!    CALCULATE Angular momentum of the earth per unit mass, J
     J = angvel*Ie/Me 
!  print*,'generalrel: J = ',J,dsqrt(j(1)*j(1)+j(2)*j(2)+j(3)*j(3))  ! |J| has a value of 981034718.66715....

!    CONVERT and assign Sun's position and velocity (m/day to m/s) vectors 
     do i=1,3
       posS(i)=sunpos(i)
! PT140130: EK's code incremented sunpos by i+1 to get the velocities. Should it be i+3 ... ?
       velS(i)=sunpos(i+3)/(24.d0*60.d0*60.d0)
     enddo

!    CALCULATE the distance between the sun and the earth
     RS = dsqrt(posS(1)*posS(1)+posS(2)*posS(2)+posS(3)*posS(3))
!    CALCULATE the dot product of position and angular momentum J
     rdotJ = pos(1)*J(1)+pos(2)*J(2)+pos(3)*J(3)
!    CALCULATE the cross product of satellite position and velocity
     rcrossv = crossproduct(pos,vel)
!    CALCULATE the cross product of satellite velocity and the Earth's angular momentum
     vcrossJ = crossproduct(vel, J)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    CALCULATE the cross product of the Sun's velocity and the scaled Sun's position
     VScrossRS = crossproduct(velS, (-gm(3)*posS/(c*c*RS*RS*RS)))
!    CALCULATE the cross product of the above VScrossRS value and the satellite velocity
     bigcross = crossproduct(VScrossRS, vel)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   CALCULATE  the general relativity contribution to acceleration 
! DEBUG: according to Charlie Lineweaver, beta and gamma are related and, in special relativity, don't have the values of 1
!        beta = v/c   and gamma = 1/sqrt(1-beta^2)
!     beta = v/c
!     gam = 1.d0/dsqrt(1.d0-beta*beta)
!
! PT140404: turn the gamma*vdotv term into just gamma*1.d0
!     beta = 1.25d0
     do i = 1, 3
       grelacc(i) = gm(1)/(c*c*r*r*r)*( (2.d0*(beta+gam)*gm(1)/r - gam*vdotv)*pos(i) + 2.d0*(1+gam)*(rdotv)*vel(i)) &    ! Schwartzchild term    (1e-8  m/s^2)
                     + (1.d0+gam) * gm(1)/(c*c*r*r*r) * ((3.d0/(r*r)*rcrossv(i)*rdotJ)+vcrossJ(i)) &                 ! Lens-Thirring effect  (1e-10 m/s^2)
                     + ((1.d0+2.d0*gam)*(bigcross(i)))                                                               ! de Sitter effect      (1e-16 m/s^2)
! DEBUG:
  factor = gm(1)/(c*c*r*r*r)
!  beta = 1.25d0
!  if (mod(t,300.d0) == 0.d0 )print*,'t,i,first term, second term,ratio', t,i, &
!        2.d0*beta*GM(1)*factor,factor*gam*vdotv*pos(i), gam*vdotv*pos(i)/(2.d0*beta*GM(1) )
!  stop
!  print*,'i,v,vdotv,2.d0*(beta+gam)*gm(1)/r*pos(i),gam*v*pos(i),2.d0*(1+gam)*(rdotv)*vel(i)',i,v,vdotv &
!           ,factor*2.d0*(beta+gam)*gm(1)/r*pos(i) &
!           ,factor*gam*vdotv*pos(i) &
!           ,factor*2.d0*(1+gam)*(rdotv)*vel(i)
!  print*,'v*pos(i),vdotv*pos(i),1.d0*pos(i)',v*pos(i)/1.d6,vdotv*pos(i)/1.d6,1.d0*pos(i)/1.d6
!  print*,'v, vdotv',v,vdotv
! 
!  print*,' '
!  print*,


! PT140205: add the acceleration perturbations caused by the Earth oblation (the J2 effect). Use
!           equations 4.2.23, 4.2.24 and 4.2.25 from Soffel (1989) "Relativity in Astrometry, Celestial
!           mechanics and Geodesy"
! PT160201: should be 5 in a1(3) equation, not 2
       if(i < 3)then
         a1(i) = 2.d0*(beta+gam)*gm(1)/(c*c*r)*J2*(Ae/r)**2 * ( (pos(i)/r)*(2.d0-(9.d0*pos(3)**2/(r*r)) )*gm(1)/(r*r) )
       else
         a1(i) = 2.d0*(beta+gam)*gm(1)/(c*c*r)*J2*(Ae/r)**2 * ( (pos(i)/r)*(5.d0-(9.d0*pos(3)**2/(r*r)) )*gm(1)/(r*r) )
       endif
       a2    = 3.d0*(gam+1)   *gm(1)/(c*c*r)*J2*(Ae/r)**2 * ( (1.d0-5.d0*(pos(3)/r)**2) * (pos(1)*vel(1)+pos(2)*vel(2))/r  &
                                                             +(3.d0 - 5.d0*(pos(3)/r)**2) * pos(3)*vel(3)/r ) * v/r  
! PT160201: according to Bob King, there is one too many 1/r in the a3 equations. Change c^2 r   to be  c^2
       if(i < 3)then
         a3(i) = -3.d0/2.d0*gam *gm(1)/(c*c)*J2*(Ae/r)**2 * ( (pos(i)/r)*(1.d0-(5.d0*pos(3)**2/(r*r)) )  *v*v/(r*r) )
       else
         a3(i) = -3.d0/2.d0*gam *gm(1)/(c*c)*J2*(Ae/r)**2 * ( (pos(i)/r)*(3.d0-(5.d0*pos(3)**2/(r*r)) )  *v*v/(r*r) )
       endif
       grelacc(i) = grelacc(i) + a1(i) + a2 + a3(i)
     enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
     return

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       contains
!      function to calculate the cross product of two vectors
       function crossproduct(a,b)
       implicit none
       real*8, dimension(3) :: a, b
       real*8, dimension(3) :: crossproduct

!      a × b = (a2b3 − a3b2) i + (a3b1 − a1b3) j + (a1b2 − a2b1) k
       crossproduct(1) = a(2)*b(3) - a(3)*b(2)
       crossproduct(2) = a(3)*b(1) - a(1)*b(3)
       crossproduct(3) = a(1)*b(2) - a(2)*b(1)

       end function crossproduct
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     end

