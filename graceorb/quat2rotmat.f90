    subroutine quat2rotmat(quat, rotmat)

! E-K Potter May 2011
! This subroutine returns the quaternion values for a given rotation matrix

    implicit none

    real (kind=8), dimension(3,3) :: rotmat
    real (kind=8), dimension(0:3) :: quat 
    
    rotmat(1, 1) = quat(0)*quat(0) + quat(1)*quat(1) - quat(2)*quat(2) - quat(3)*quat(3)
    rotmat(1, 2) = 2.d0*quat(1)*quat(2) + 2.d0*quat(0)*quat(3)
    rotmat(1, 3) = 2.d0*quat(1)*quat(3) - 2.d0*quat(0)*quat(2)
    rotmat(2, 1) = 2.d0*quat(1)*quat(2) - 2.d0*quat(0)*quat(3)
    rotmat(2, 2) = quat(0)*quat(0) - quat(1)*quat(1) + quat(2)*quat(2) - quat(3)*quat(3)
    rotmat(2, 3) = 2.d0*quat(2)*quat(3) + 2.d0*quat(0)*quat(1)
    rotmat(3, 1) = 2.d0*quat(1)*quat(3) + 2.d0*quat(0)*quat(2)
    rotmat(3, 2) = 2.d0*quat(2)*quat(3) - 2.d0*quat(0)*quat(1)
    rotmat(3, 3) = quat(0)*quat(0) - quat(1)*quat(1) - quat(2)*quat(2) + quat(3)*quat(3)

    return
    end

