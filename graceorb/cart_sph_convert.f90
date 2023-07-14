  subroutine cart_sph_convert(lat, long, R, cart) 

! This subroutine converts  
! cartesian coordinates cart(1), (2), (3) into 
! spherical coordinates R, lat, long

  real (kind=8) :: lat, long, R
  real (kind=8), dimension(3) :: cart 
  real (kind=8) :: pi

  pi = 4.0d0*datan(1.d0) 

  R  = dsqrt(cart(1)**2+cart(2)**2+cart(3)**2)
  lat = dasin(cart(3)/R)
  long = datan2(cart(2),cart(1))

  return
  END

