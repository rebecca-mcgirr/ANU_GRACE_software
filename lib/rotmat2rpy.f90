  subroutine rotmat2rpy(matrix,rpy)

! given a 3 x 3 rotation matrix, calculate the three Euler angles which are the rotations about the XYZ axes (roll, pitch, yaw)
!
! P. Tregoning
! 21 March 2014

  implicit none

  real(kind=8), intent(in ) :: matrix(3,3)             ! input rotation matrix
  real(kind=8), intent(out) :: rpy(3)                  ! roll, pitch, yaw euler angles


  rpy(2) = -dasin(matrix(3,1))                     ! pitch, or rotation about Y
  rpy(3) =  dasin( matrix(2,1) / dcos(rpy(2)) )    ! yaw  , or rotation about Z
  rpy(1) =  dasin( matrix(3,2) / dcos(rpy(2)) )    ! roll , or rotation about X

  return
  end subroutine rotmat2rpy

