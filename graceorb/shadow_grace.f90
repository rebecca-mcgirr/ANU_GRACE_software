  subroutine shadow_grace(sunpos,sbcor,shadow)

! subroutine to calculate whether the GRACE satellite(s?) are in the Earth's shadow or not.
! Code adapted from gg/gamit/arc/shadow.f
!
! P. Tregoning
! 19 June 2013
!
! In:
!      sunpos:  sun position
!      sbcor :  Earth-satellite vector
! OUT:
!      shadow value (0 to 1)
!
!
!
!  bcor  :  sun-satellite vector
! ubcor  :  normalised sun-satellite vector

  implicit none

  double precision sunpos(3), sbcor(3), shadow
  double precision mag_sunsat,mag_earthsat,rsbx,rs,rp,sunrad,ertrad,sep,dot
  double precision sbcor_norm(3),bcor(3),sepp(3),ubcor(3)
  integer i,j,idbprt,idebug
  double precision earth_ang,sun_ang,sunsatearth_ang,satearth(3)

  mag_sunsat = 0.d0
  mag_earthsat = 0.d0 
  ertrad = 6378136.3d0   ! value taken from arc/egm08.f
  sunrad = 696000.D+03   ! value taken from arc/egm08.f 

! calculate the satellite-sun vector, the magnitude of the earth-sat and sun-sat vectors
! sunpos (in km) : earth wrt sun, sbcor: satellite wrt earth. Therefore, bcor (sat wrt sun) = -(sunpos+sbcor)
  bcor = -(sunpos*1.d3 - sbcor)
  mag_sunsat = dsqrt(bcor(1)**2+bcor(2)**2+bcor(3)**2)
  mag_earthsat = dsqrt(sbcor(1)**2+sbcor(2)**2+sbcor(3)**2)

! change the sign of sbcor so that it is the sat-earth vector
  satearth = -sbcor

! get the angle subtended at the satellite by half the disc of the Earth
  earth_ang = ertrad/mag_earthsat

! get the angle subtended at the satellite by half the disc of the sun 
  sun_ang = sunrad/mag_sunsat

! calculate the angle at the satellite between the sat-sun and sat-earth vectors
  sunsatearth_ang = dacos( dot(satearth,bcor)/(mag_sunsat*mag_earthsat) )

! now, the satellite will be in eclipse if the angle subtended by half of the earth disc is greater than the sun-sat-earth angle plus the angle subtended by half of the sun disc.
  if (earth_ang > sunsatearth_ang + sun_ang)then
!    print*,"Eclipse: earth_ang,sun_ang,sunsatearth_ang", earth_ang,sun_ang+sunsatearth_ang
!     shadow = 1.d0
  else
!     shadow = 0.d0
!    print*,"No eclipse: earth_ang,sun_ang,sunsatearth_ang", earth_ang,sun_ang+sunsatearth_ang
  endif
! normalise the satellite-sun vector and earth-satellite vector. Why??
!  ubcor = bcor/mag_sunsat
!  sbcor_norm = sbcor/mag_earthsat
!
! rsbx is the projection of sbcor along bcor (earth-sat vector along sun-sat vector)
!  rsbx = dot(sbcor_norm,ubcor)
!  print*,'rsbx = ',rsbx,sbcor_norm,ubcor
!
!! cross of earth-sat and sun-sat vectors to generate "sepp" (why? not sure why!)
! call cross(sbcor_norm,ubcor,sepp) 
!  print*,'sepp=',sepp
!
! rs, rp are apparent (from satellite) radii of sun and earth
! sep is apparent separation of their centers
!  rs=sunrad/mag_sunsat
!  rp=ertrad/rsbx
!  sep=dsqrt(sepp(1)**2+sepp(2)**2+sepp(3)**2)/rsbx    

!  print*,'calling get_lambda with',  rs,rp,sep  

  call get_lambda( sun_ang,earth_ang,sunsatearth_ang,shadow,idbprt,idebug )   
!  print*,"from get_lambda: earth_ang,sun_ang+sunsatearth_ang", earth_ang,sun_ang+sunsatearth_ang,shadow,sunsatearth_ang,sun_ang

!  if (shadow .eq. 0 ) stop "stopped in shadow_grace"
!
!  print*,"shadow_grace: shadow = ",shadow

!  stop "stopped in shadow_grace"
  return
  end
