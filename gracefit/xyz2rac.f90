subroutine xyz2rac( posvel,dxyz,drac )

  !   Convert delta x,y,z to radial, along-track, cross-track

  !   posvel(6)    position,velocity of SV
  !   dxyz(3)      cartesian difference
  !   drac(3)      rad,along,cross difference

  implicit none

  real*8 :: posvel(6),dxyz(3),drac(3)
  real*8 :: radlen,sprodr,vellen,c(3),sproda
  real*8 :: dot

  !   Length of position vector
  radlen = dsqrt( dot(posvel,posvel) )

  !   Projection of dxyz onto position vector
  sprodr = dot(dxyz,posvel)

  !   Radial difference
  drac(1) = sprodr/radlen

  !   Length of velocity vector
  vellen = dsqrt( dot(posvel(4),posvel(4)) )

  !   Projection of dxyz onto velocity vector
  sproda = dot(dxyz,posvel(4))

  !   Along-track difference
  drac(2) = sproda/vellen

  !   Cross-track difference
  call cross(posvel(1),posvel(4),c)
  drac(3) = dot(dxyz,c)


  return
end subroutine xyz2rac
