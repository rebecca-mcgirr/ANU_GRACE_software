    subroutine rpy2rotmat(roll, pitch, yaw, rotmat, rotmat_deriv_roll, rotmat_deriv_pitch, rotmat_deriv_yaw)

! Tony Purcell April 2013
! This subroutine calculates the rotation matrix defined by a series of three angles.
! These angles represent rotations around the x-axis, y-axis and z-axis respectively,
! corresponding to roll, pitch and yaw.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The formulae used for the rotation matrices are:
!
!             Rx               |            Ry              |              Rz
!
!      1      0         0      |   cos(th)    0    sin(th)  |   cos(th)  -sin(th)    0        
!      0    cos(th)  -sin(th)  |     0        1      0      |   sin(th)   cos(th)    0
!      0    sin(th)   cos(th)  |  -sin(th)    0    cos(th)  |     0         0        1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The formulae used for the derivatives of the rotation matrices are:
!
!             R'x              |            R'y             |              Rz
!
!      0      0         0      |  -sin(th)    0    cos(th)  |  -sin(th)  -cos(th)    0        
!      0   -sin(th)  -cos(th)  |     0        0      0      |   cos(th)  -sin(th)    0
!      0    cos(th)  -sin(th)  |  -cos(th)    0   -sin(th)  |     0         0        0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    implicit none

    real (kind=8), dimension(3, 3) :: a, b, rotmat, rotmat_roll, rotmat_pitch, rotmat_yaw, rottmp
    real (kind=8), dimension(3, 3) :: rotmat_deriv_roll, rotmat_deriv_pitch, rotmat_deriv_yaw
    real (kind=8), dimension(3, 3) :: deriv_roll, deriv_pitch, deriv_yaw
    real (kind=8) :: roll, pitch, yaw

    integer*4 :: i, j, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set rotation matrix entries

    rotmat_roll = 0.0d0
    rotmat_pitch = 0.0d0
    rotmat_yaw = 0.0d0

    rotmat_roll(1, 1) = 1.0d0
    rotmat_roll(2, 2) = dcos(roll)
    rotmat_roll(2, 3) = -dsin(roll)
    rotmat_roll(3, 2) = dsin(roll)
    rotmat_roll(3, 3) = dcos(roll)

! DEBUG:
!    write(6,*) ""
!    write(6,*) "rotmat_roll:"
!    do i = 1, 3
!      write(6,*) (rotmat_roll(i,j), j = 1, 3)
!    enddo
!    write(6,*) ""

    rotmat_pitch(1, 1) = dcos(pitch)
    rotmat_pitch(1, 3) = dsin(pitch)
    rotmat_pitch(2, 2) = 1.0d0
    rotmat_pitch(3, 1) = -dsin(pitch)
    rotmat_pitch(3, 3) = dcos(pitch)

! DEBUG:
!    write(6,*) ""
!    write(6,*) "rotmat_pitch:"
!    do i = 1, 3
!      write(6,*) (rotmat_pitch(i,j), j = 1, 3)
!    enddo
!    write(6,*) ""

    rotmat_yaw(1, 1) = dcos(yaw)
    rotmat_yaw(1, 2) = -dsin(yaw)
    rotmat_yaw(2, 1) = dsin(yaw)
    rotmat_yaw(2, 2) = dcos(yaw)
    rotmat_yaw(3, 3) = 1.0d0

! DEBUG:
!    write(6,*) ""
!    write(6,*) "rotmat_yaw:"
!    do i = 1, 3
!      write(6,*) (rotmat_yaw(i,j), j = 1, 3)
!    enddo
!    write(6,*) ""

! Now combine to form the complete rotational matrix

    call matmult(rotmat_roll, rotmat_pitch, rottmp, 3, 3, 3)
    call matmult(     rottmp,   rotmat_yaw, rotmat, 3, 3, 3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now calculate the derivatives of the rotational matrices

    deriv_roll = 0.0d0
    deriv_pitch = 0.0d0
    deriv_yaw = 0.0d0

    deriv_roll(2, 2) = -dsin(roll)
    deriv_roll(2, 3) = -dcos(roll)
    deriv_roll(3, 2) = dcos(roll)
    deriv_roll(3, 3) = -dsin(roll)

    deriv_pitch(1, 1) = -dsin(pitch)
    deriv_pitch(1, 3) = dcos(pitch)
    deriv_pitch(3, 1) = -dcos(pitch)
    deriv_pitch(3, 3) = -dsin(pitch)

    deriv_yaw(1, 1) = -dsin(yaw)
    deriv_yaw(1, 2) = -dcos(yaw)
    deriv_yaw(2, 1) = dcos(yaw)
    deriv_yaw(2, 2) = -dsin(yaw)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now calculate the derivatives of the combined rotational matrix with respect to each of the
! rotation corrections 

    call matmult( deriv_roll, rotmat_pitch,             rottmp, 3, 3, 3)
    call matmult(     rottmp,   rotmat_yaw,  rotmat_deriv_roll, 3, 3, 3)

    call matmult(rotmat_roll,  deriv_pitch,             rottmp, 3, 3, 3)
    call matmult(     rottmp,   rotmat_yaw, rotmat_deriv_pitch, 3, 3, 3)

    call matmult(rotmat_roll, rotmat_pitch,            rottmp, 3, 3, 3)
    call matmult(     rottmp,    deriv_yaw,  rotmat_deriv_yaw, 3, 3, 3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    return
    end
