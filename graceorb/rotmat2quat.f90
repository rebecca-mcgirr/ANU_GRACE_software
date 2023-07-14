    subroutine rotmat2quat(rotmat, quat)

! E-K Potter May 2011
! This subroutine returns the quaternion values for a given rotation matrix

    implicit none

    real (kind=8), dimension(3,3) :: rotmat
    real (kind=8), dimension(0:3) :: quat 
    
    if (1.d0+rotmat(1,1)+rotmat(2,2)+rotmat(3,3).lt.0.d0) then
      print*, "STOPPING because quaternion cannot be calculated"
    endif
    quat(0) = dsqrt(1.d0+rotmat(1,1)+rotmat(2,2)+rotmat(3,3))/2
    quat(1) = (rotmat(2,3)-rotmat(3,2))/(4.d0*quat(0))
    quat(2) = (rotmat(3,1)-rotmat(1,3))/(4.d0*quat(0))
    quat(3) = (rotmat(1,2)-rotmat(2,1))/(4.d0*quat(0))

    return
    end

