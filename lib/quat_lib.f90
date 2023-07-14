subroutine quat_conj ( q )
!
!*******************************************************************************
!
!! QUAT_CONJ conjugates a quaternion.
!
!
!  Discussion:
!
!    A quaternion is a quadruplet (A,B,C,D) of real numbers, which
!    may be written as
!
!      Q = A + Bi + Cj + Dk.
!
!    The conjugate of Q is
!
!      conj ( Q ) = A - Bi - Cj - Dk.
!
!  Modified:
!
!    29 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real Q(4).  On input, the quaternion to be conjugated.
!    On output, the conjugated quaternion.
!
  implicit none
!
  real(kind=8) :: q(4)
!
  q(2) = - q(2)
  q(3) = - q(3)
  q(4) = - q(4)

  return
end

subroutine quat_inv ( q )
!
!*******************************************************************************
!
!! QUAT_INV inverts a quaternion.
!
!
!  Discussion:
!
!    A quaternion is a quadruplet (A,B,C,D) of real numbers, which
!    may be written as
!
!      Q = A + Bi + Cj + Dk.
!
!    The inverse of Q is
!
!      inverse ( Q ) = conjugate ( Q ) / ( norm ( Q ) )**2.
!
!  Modified:
!
!    16 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real Q(4).  On input, the quaternion to be inverted.
!    On output, the inverse of the input quaternion.
!
  implicit none
!
  real(kind=8) :: q(4)
!
  q(1:4) = q(1:4) / sum ( q(1:4)**2 ) 
  q(2:4) = - q(2:4)

  return
end

subroutine quat_mul ( q1, q2, q3 )
!
!*******************************************************************************
!
!! QUAT_MUL multiplies two quaternions.
!
!
!  Discussion:
!
!    A quaternion is a quadruplet (A,B,C,D) of real numbers, which
!    may be written as
!
!      Q = A + Bi + Cj + Dk.
!
!    To multiply two quaternions, use the relationships:
!
!      ij = -ji = k
!      jk = -kj = i
!      ki = -ik = j
!      ii = jj = kk = -1
!
!  Modified:
!
!    29 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real Q1(4), Q2(4), the two quaternions to be multiplied.
!
!    Output, real Q3(4), the product of the two quaternions.
!
  implicit none
!
  real(kind=8) :: q1(4)
  real(kind=8) :: q2(4)
  real(kind=8) :: q3(4)
!
  q3(1) = q1(1) * q2(1) - q1(2) * q2(2) - q1(3) * q2(3) - q1(4) * q2(4)
  q3(2) = q1(1) * q2(2) + q1(2) * q2(1) + q1(3) * q2(4) - q1(4) * q2(3)
  q3(3) = q1(1) * q2(3) - q1(2) * q2(4) + q1(3) * q2(1) + q1(4) * q2(2)
  q3(4) = q1(1) * q2(4) + q1(2) * q2(3) - q1(3) * q2(2) + q1(4) * q2(1)

  return
end
function quat_norm ( q )
!
!*******************************************************************************
!
!! QUAT_NORM computes the norm of a quaternion.
!
!
!  Discussion:
!
!    A quaternion is a quadruplet (A,B,C,D) of real numbers, which
!    may be written as
!
!      Q = A + Bi + Cj + Dk.
!
!    The norm of Q is
!
!      norm(Q) = sqrt ( A**2 + B**2 + C**2 + D**2 ).
!
!  Modified:
!
!    29 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real Q(4), the quaternion.
!
!    Output, real QUAT_NORM, the norm of the quaternion.
!
  implicit none
!
  real(kind=8) :: q(4)
  real(kind=8) :: quat_norm
!
  quat_norm = sqrt ( sum ( q(1:4)**2 ) )

  return
end

subroutine rotation_axis2quat_3d ( axis, angle, q )
!
!*******************************************************************************
!
!! ROTATION_AXIS2QUAT_3D converts a rotation from axis to quaternion format in 3D.
!
!
!  Definition:
!
!    A rotation quaternion Q has the form:
!
!      Q = A + Bi + Cj + Dk
!
!    where A, B, C and D are real numbers, and i, j, and k are to be regarded
!    as symbolic constant basis vectors, similar to the role of the "i"
!    in the representation of imaginary numbers.
!
!    A is the cosine of half of the angle of rotation.  (B,C,D) is a
!    unit vector pointing in the direction of the axis of rotation.
!    Rotation multiplication and inversion can be carried out using
!    this format and the usual rules for quaternion multiplication
!    and inversion.
!
!  Modified:
!
!    24 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real AXIS(3), the axis vector which remains unchanged by
!    the rotation.
!
!    Input, real ANGLE, the angular measurement of the rotation about
!    the axis, in radians.
!
!    Output, real Q(4), the quaternion representing the rotation.
!
  implicit none
!
  real(kind=8) :: axis(3)
  real(kind=8) :: angle
  real(kind=8) :: norm
  real(kind=8) :: q(4)
!
  norm = sqrt ( axis(1) * axis(1) + axis(2) * axis(2) + axis(3) * axis(3) )

  if ( norm == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ROTATION_AXIS2QUAT_3D - Fatal error!'
    write ( *, '(a)' ) '  The axis vector is null.'
  end if

  q(1) = cos ( 0.5E+00 * angle )

  q(2) = axis(1) * sin ( 0.5E+00 * angle ) / norm
  q(3) = axis(2) * sin ( 0.5E+00 * angle ) / norm
  q(4) = axis(3) * sin ( 0.5E+00 * angle ) / norm

  return
end

subroutine rotation_quat2axis_3d ( q, axis, angle )
!
!*******************************************************************************
!
!! ROTATION_QUAT2AXIS_3D converts a rotation from quaternion to axis format in 3D.
!
!
!  Definition:
!
!    A rotation quaternion Q has the form:
!
!      Q = A + Bi + Cj + Dk
!
!    where A, B, C and D are real numbers, and i, j, and k are to be regarded
!    as symbolic constant basis vectors, similar to the role of the "i"
!    in the representation of imaginary numbers.
!
!    A is the cosine of half of the angle of rotation.  (B,C,D) is a
!    vector pointing in the direction of the axis of rotation.
!    Rotation multiplication and inversion can be carried out using
!    this format and the usual rules for quaternion multiplication
!    and inversion.
!
!  Modified:
!
!    02 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real Q(4), the quaternion representing the rotation.
!
!    Output, real AXIS(3), the axis vector which remains unchanged by
!    the rotation.
!
!    Output, real ANGLE, the angular measurement of the rotation about
!    the axis, in radians.
!
  implicit none
!
  real(kind=8) :: axis(3)
  real(kind=8) :: angle
  real(kind=8) :: cos_phi
  real(kind=8) :: q(4)
  real(kind=8) :: sin_phi
!
  sin_phi = sqrt ( sum ( q(2:4)**2 ) )

  cos_phi = q(1)

  angle = 2.0E+00 * atan2 ( sin_phi, cos_phi )

  if ( sin_phi == 0.0E+00 ) then
    axis(1:3) = (/ 1.0E+00, 0.0E+00, 0.0E+00 /)
  else
    axis(1:3) = (/ q(2), q(3), q(4) /) / sin_phi
  end if

  return
end

subroutine rotation_mat2quat_3d ( a, q )
!
!*******************************************************************************
!
!! ROTATION_MAT2QUAT_3D converts a rotation from matrix to quaternion format in 3D.
!
!
!  Discussion:
!
!    The computation is based on the fact that a rotation matrix must
!    have an eigenvector corresponding to the eigenvalue of 1, hence:
!
!      ( A - I ) * v = 0.
!
!    The eigenvector V is the axis of rotation.
!
!  Reference:
!
!    Jack Kuipers
!    Quaternions and Rotation Sequences,
!    Princeton, 1998.
!
!  Modified:
!
!    27 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A(3,3), the rotation matrix.
!
!    Output, real Q(4), the quaternion representing the rotation.
!
  implicit none
!
  real(kind=8) :: a(3,3)
  real(kind=8) :: angle
  real(kind=8) :: cos_phi
  real(kind=8) :: norm
  real(kind=8) :: q(4)
  real(kind=8) :: sin_phi
!
  norm = sqrt ( ( a(3,2) - a(2,3) )**2 + ( a(1,3) - a(3,1) )**2 + ( a(2,1) - a(1,2) )**2 )

  if ( norm == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ROTATION_MAT2AXIS_3D - Fatal error!'
    write ( *, '(a)' ) '  A is not a rotation matrix,'
    write ( *, '(a)' ) '  or there are multiple axes of rotation.'
  end if

  angle = acos( 0.5E+00 * ( a(1,1) + a(2,2) + a(3,3) - 1.0E+00 ) )

  cos_phi = cos ( 0.5E+00 * angle )

  sin_phi = sqrt ( 1.0E+00 - cos_phi**2 )

  q(1) = cos_phi
  q(2) = sin_phi * ( a(3,2) - a(2,3) ) / norm
  q(3) = sin_phi * ( a(1,3) - a(3,1) ) / norm
  q(4) = sin_phi * ( a(2,1) - a(1,2) ) / norm

  return
end

subroutine rotation_quat2mat_3d ( q, a )
!
!*******************************************************************************
!
!! ROTATION_QUAT2MAT_3D converts a rotation from quaternion to matrix format in 3D.
!
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Modified:
!
!    27 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real AXIS(3), the axis vector which remains unchanged by
!    the rotation.
!
!    Input, real ANGLE, the angular measurement of the rotation about
!    the axis, in radians.
!
!    Output, real A(3,3), the rotation matrix.
!
  implicit none
!
  real(kind=8) :: a(3,3)
  real(kind=8) :: angle
  real(kind=8) :: ca
  real(kind=8) :: cos_phi
  real(kind=8) :: q(4)
  real(kind=8) :: sa
  real(kind=8) :: sin_phi
  real(kind=8) :: v1
  real(kind=8) :: v2
  real(kind=8) :: v3
!
  sin_phi = sqrt ( sum ( q(2:4)**2 ) )

  cos_phi = q(1)

  angle = 2.0E+00 * atan2 ( sin_phi, cos_phi )

  if ( sin_phi == 0.0E+00 ) then
    v1 = 1.0E+00
    v2 = 0.0E+00
    v3 = 0.0E+00
  else
    v1 = q(2) / sin_phi
    v2 = q(3) / sin_phi
    v3 = q(4) / sin_phi
  end if

  ca = cos ( angle )
  sa = sin ( angle )

  a(1,1) =                    v1 * v1 + ca * ( 1.0E+00 - v1 * v1 )
  a(1,2) = ( 1.0E+00 - ca ) * v1 * v2 - sa * v3
  a(1,3) = ( 1.0E+00 - ca ) * v1 * v3 + sa * v2

  a(2,1) = ( 1.0E+00 - ca ) * v2 * v1 + sa * v3
  a(2,2) =                    v2 * v2 + ca * ( 1.0E+00 - v2 * v2 )
  a(2,3) = ( 1.0E+00 - ca ) * v2 * v3 - sa * v1

  a(3,1) = ( 1.0E+00 - ca ) * v3 * v1 - sa * v2
  a(3,2) = ( 1.0E+00 - ca ) * v3 * v2 + sa * v1
  a(3,3) =                    v3 * v3 + ca * ( 1.0E+00 - v3 * v3 )

  return
end

subroutine quat_rot_vect ( q, v, w )
!
!*******************************************************************************
!
!! QUAT_ROT_VECT applies a quaternion rotation to a vector in 3d.
!
!
!  Discussion:
!
!    If Q is a unit quaternion that encodes a rotation of ANGLE
!    radians about the vector AXIS, then for an arbitrary real
!    vector V, the result W of the rotation on V can be written as:
!
!      W = Q * V * Conj(Q)
!
!  Modified:
!
!    29 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real Q(4), the quaternion defining the rotation.
!
!    Input, real V(3), the vector to be rotated.
!
!    Output, real W(3), the rotated vector.
!
  implicit none
!
  real(kind=8) :: q(4)
  real(kind=8) :: v(3)
  real(kind=8) :: w(3)
!
  w(1) = &
         ( 2.0E+00 * ( q(1) * q(1) + q(2) * q(2) ) - 1.0E+00 ) * v(1) &
       +   2.0E+00 * ( q(2) * q(3) - q(1) * q(4) )             * v(2) &
       +   2.0E+00 * ( q(2) * q(4) + q(1) * q(3) )             * v(3)

  w(2) = &
           2.0E+00 * ( q(2) * q(3) + q(1) * q(4) )             * v(1) &
       + ( 2.0E+00 * ( q(1) * q(1) + q(3) * q(3) ) - 1.0E+00 ) * v(2) &
       +   2.0E+00 * ( q(3) * q(4) - q(1) * q(2) )             * v(3)

  w(3) = &
           2.0E+00 * ( q(2) * q(4) - q(1) * q(3) )             * v(1) &
       +   2.0E+00 * ( q(3) * q(4) + q(1) * q(2) )             * v(2) &
       + ( 2.0E+00 * ( q(1) * q(1) + q(4) * q(4) ) - 1.0E+00 ) * v(3)

  return
end

!
!*******************************************************************************

subroutine quat_to_rpy( quat, roll, pitch, yaw )
!
! subroutine to compute from a quaternion the corresponding yaw, pitch and roll angles
!
! P. Tregoning
! 10 April 2013

  implicit none

  double precision quat(4),yaw,pitch,roll

  roll=atan2(2.d0*(quat(3)*quat(4)+quat(1)*quat(2)),(quat(1)*quat(1)-quat(2)*quat(2)-quat(3)*quat(3)+quat(4)*quat(4)))

  yaw=atan2(2.d0*(quat(2)*quat(3)+quat(1)*quat(4)),(quat(1)*quat(1)+quat(2)*quat(2)-quat(3)*quat(3)-quat(4)*quat(4)))

  pitch = asin(-2.d0*( quat(2)*quat(4) - quat(1)*quat(3) ) )


  return
end


!
!*******************************************************************************

subroutine quat_to_rpy_new( quat, roll, pitch, yaw )
!
! subroutine to compute from a quaternion the corresponding yaw, pitch and roll angles
!
! P. Tregoning
! 10 April 2013
!
! Modified: APP170313
! original form was based on the rotation matrix for Euler angles not yaw, pitch, roll

  implicit none

  double precision quat(4),yaw,pitch,roll

  roll=-atan2(2.d0*(quat(3)*quat(4)-quat(1)*quat(2)),(quat(1)*quat(1)-quat(2)*quat(2)-quat(3)*quat(3)+quat(4)*quat(4)))

  yaw=-atan2(2.d0*(quat(2)*quat(3)-quat(1)*quat(4)),(quat(1)*quat(1)+quat(2)*quat(2)-quat(3)*quat(3)-quat(4)*quat(4)))

  pitch = asin(2.d0*( quat(2)*quat(4) + quat(1)*quat(3) ) )


  return
end


!
!*******************************************************************************

subroutine quat_to_ypr_new( quat, roll, pitch, yaw )
!
! subroutine to compute from a quaternion the corresponding yaw, pitch and roll angles
!
! P. Tregoning
! 10 April 2013
!
! Modified: APP170313
! original form was based on the rotation matrix for Euler angles not yaw, pitch, roll

  implicit none

  double precision quat(4),yaw,pitch,roll

  roll=-atan2(2.d0*(quat(3)*quat(4)-quat(1)*quat(2)),(quat(1)*quat(1)-quat(2)*quat(2)-quat(3)*quat(3)+quat(4)*quat(4)))

  yaw=-atan2(2.d0*(quat(2)*quat(3)-quat(1)*quat(4)),(quat(1)*quat(1)+quat(2)*quat(2)-quat(3)*quat(3)-quat(4)*quat(4)))

  pitch = asin(2.d0*( quat(2)*quat(4) + quat(1)*quat(3) ) )


  return
end


!
!*******************************************************************************

subroutine quat_to_ypr( quat, roll, pitch, yaw )
!
! subroutine to compute from a quaternion the corresponding yaw, pitch and roll angles
!
! P. Tregoning
! 10 April 2013

  implicit none

  double precision quat(4),yaw,pitch,roll

  roll=atan2(2.d0*(quat(3)*quat(4)+quat(1)*quat(2)),(quat(1)*quat(1)-quat(2)*quat(2)-quat(3)*quat(3)+quat(4)*quat(4)))

  yaw=atan2(2.d0*(quat(2)*quat(3)+quat(1)*quat(4)),(quat(1)*quat(1)+quat(2)*quat(2)-quat(3)*quat(3)-quat(4)*quat(4)))

  pitch = asin(-2.d0*( quat(2)*quat(4) - quat(1)*quat(3) ) )


  return
end


