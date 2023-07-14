  subroutine cart_to_sph(xyz,sph)

! convert cartesian coordinates into spherical lat, lon, radius
!
! P. Tregoning
! 26 June 2017

  implicit none

  real(kind=8),  intent(in)  :: xyz(3)      ! cartesian coordinates to be transformed
  real(kind=8),  intent(out) :: sph(3)      ! lat/lon/radius (in radians and metres)
  real(kind=8) :: pi
  real(kind=8) amag3

  pi = 4.d0*datan(1.d0)

! first, the radius
  sph(3) = dsqrt(xyz(1)**2+xyz(2)**2+xyz(3)**2)

! now the latitude
 sph(1) = dasin(xyz(3)/sph(3))

! now the longitude
  if(xyz(1) == 0.d0 .and. xyz(2) > 0.d0)then
    sph(2) = 0.d0
  elseif (xyz(1) == 0.d0 .and. xyz(2) < 0.d0)then
    sph(2) = pi
  else if (xyz(1) > 0.d0 .and. xyz(2) == 0.d0)then
    sph(2) = pi/2.d0
  else if (xyz(1) < 0.d0 .and. xyz(2) == 0.d0)then
    sph(2) = 3.d0*pi/2.d0
  else
    sph(2) = datan(xyz(2)/xyz(1))
    if(sph(2) < 0.d0 .and. xyz(1) < 0.d0)then
      sph(2) = pi+sph(2)
    elseif(sph(2) > 0.d0 .and. xyz(1) < 0.d0)then
      sph(2) = pi+sph(2)
    elseif(sph(2) < 0.d0 .and. xyz(2) < 0.d0)then
      sph(2) = 2.d0*pi+sph(2)
    endif
  endif

  return
  end

