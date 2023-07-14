subroutine los_to_srf(vec1,vec2,quat_los_srf)

  ! subroutine to generate the quaternion that will rotate the LOS vector into the nominal satellite reference frame
  !
  ! IN:  
  !     vec1, vec2 : the XYZ positions of the two satellites
  !
  ! OUT:
  !     quat       : the quaternion for the satellite of coordinates vec1.
  !
  ! P. Tregoning
  ! 11 April 2013

  implicit none

  double precision vec1(3),vec2(3),quat_los_srf(4)
  double precision LOS(3),cross_LOS(3),rad_LOS(3),cross_in_SRF(3),rad_in_SRF(3)
  double precision srf_X(3),srf_Y(3),srf_Z(3),quat_x(4),quat_y(4),quat_z(4),quat_xy(4)
  double precision quat_tmp(4),quat_c(4),tmpvec(3)
  integer i
  logical debug

  debug = .true.

  ! first compute the line of sight vector
  do i = 1,3
     LOS(i) = vec2(i) - vec1(i)
     print*,"los_to_srf LOS(i)",LOS(i),vec2(i),vec1(i)
  enddo

  ! find the vector orthogonal to the LOS vector (the cross-track direction). To do this, we use the cross product of the two position vectors.
  call cross (vec2,vec1,cross_LOS)
  if(debug)print*,'cross_LOS',cross_LOS

  ! and the third orthogonal axis
  call cross (LOS,cross_LOS,rad_LOS)
  if(debug)print*,'rad_LOS',rad_LOS

  ! find the quaternion that aligns LOS with SRF X
  srf_X = 0.d0
  srf_X(1) = 1.d0
  call findquat(LOS,srf_X,quat_x)
  if(debug)print*,'quat_x ',quat_x

  ! apply quat_x to cross_LOS
  call quat_rot_vect(quat_x,cross_LOS,cross_in_SRF)

  ! now find the quaternion that aligns it with SRF Y
  srf_Y = 0.d0
  srf_Y(2) = 1.d0
  call findquat(cross_in_SRF,srf_Y,quat_y)

  ! combine quat_x and quat_y
  call quat_mul(quat_y,quat_x,quat_xy)

  ! apply quat_xy to rad_LOS
  call quat_rot_vect(quat_xy,rad_LOS,rad_in_SRF)

  ! and the quaternion that aligns it with SRF Z
  srf_Z = 0.d0
  srf_Z(3) = 1.d0
  call findquat(rad_in_SRF,srf_Z,quat_z)

  ! combine quat_z and quat_xy
  call quat_mul(quat_z,quat_xy,quat_los_srf)

  ! we're done.
  return

end subroutine los_to_srf

