    subroutine findquat_ISRF2EF(u, v, quat)

! Tony Purcell April 11 2013
! This subroutine returns the quaternion defining a rotation from an idealised satellite reference frame (SRF)
! to an inertial reference frame. In the idealised satellite reference frame the satellites are assumed to be 
! perfectly aligned with the line of sight vector between them with zero yaw and roll

    implicit none

    real (kind=8), dimension(3) :: u, v            ! the position vectors of the two spacecraft in EF coords
    real (kind=8), dimension(0:3) :: quat, quat_2  ! the quaternion describing the rotation from SRF to EF

    real (kind=8), dimension(3) :: x, y, z, LOS, x_axis, y_axis, z_axis
    real (kind=8), dimension(0:3) :: quat_x, quat_y, quat_z, quat_xy, quat_tmp
    real (kind=8), dimension(0:3) :: quat_x_2, quat_y, quat_z, quat_xy_2, quat_xyz_2

    real (kind=8) :: mag_u, mag_v, mag_factor, mag_quat, mag_cross, mag_LOS, mag_y, mag_z, mag_LOS
    integer :: k

    quat = 0.0d0

    mag_u = dsqrt(u(1)**2 + u(2)**2 + u(3)**2)
    mag_v = dsqrt(v(1)**2 + v(2)**2 + v(3)**2)
 
! Check that neither of the input vectors is small enough to cause numerical problems
! if either of them is, return an error message and stop run

    if ( min(mag_u, mag_v) .lt. 1.0d-16 ) then
      write(6, *) "ERROR: Subroutine findquat_ISRF2EF - zero vector input"
      stop
    endif

! convert to unit vectors

    u = u/mag_u
    v = v/mag_v

! calculate unit line of sight vector (corresponds to x-axis in idealised SRF

    LOS(1) = u(1) - v(1) 
    LOS(2) = u(2) - v(2) 
    LOS(3) = u(3) - v(3) 

    mag_LOS = dsqrt(LOS(1)**2 + LOS(2)**2 + LOS(3)**2)

! check that the difference between the vectors is non-zero
! if not, return an error message and stop run

    if ( mag_LOS .lt. 1.0d-16 ) then
      write(6, *) "ERROR: Subroutine findquat_ISRF2EF - identical vectors input"
      stop
    endif

    LOS = LOS/mag_LOS

! Calculate cross product of the unit vectors and store in quaternion vector entries to define
! the axis of rotation. Order here is chosen so that the z-axis will point toward the satellite's nadir
! (that is, toward the Earth)

    call cross(LOS, u, y)

!    y(1) = (u(2) * LOS(3) - u(3) * LOS(2))
!    y(2) = (u(3) * LOS(1) - u(1) * LOS(3))
!    y(3) = (u(1) * LOS(2) - u(2) * LOS(1))

! test if vectors are parallel, if so return an error message and stop run

    mag_y = dsqrt(y(1)**2 + y(2)**2 + y(3)**2)

    if ( mag_y .lt. 1.0d-16 ) then
      write(6, *) "ERROR: Subroutine findquat_ISRF2EF - parallel vectors input"
      stop
    endif

    y = y/mag_y

! calculate unit z-axis vector (obtained from LOS x y)

    call cross(LOS, y, z)
    mag_z = dsqrt(z(1)**2 + z(2)**2 + z(3)**2)
    z = z/mag_z

!    z(1) = (LOS(2) * y(3) - LOS(3) * y(2))
!    z(2) = (LOS(3) * y(1) - LOS(1) * y(3))
!    z(3) = (LOS(1) * y(2) - LOS(2) * y(1))

! Find quaternion that rotates the SRF x-axis onto the line of sight vector

    x_axis = 0.0d0
    x_axis(1) = 1.0d0

    call findquat(x_axis, LOS, quat_x)

! Now find quaternion describing rotation from rotated y-axis to EF y-axis

    y_axis = 0.0d0
    y_axis(2) = 1.0d0

    call quat_rot_vect(quat_x, y_axis, rot_y_axis)
    call findquat(rot_y_axis, y, quat_y)

! combine the two quaternions

    call quat_mul(quat_y, quat_x, quat_xy)

! Now find quaternion describing rotation from rotated z-axis to EF z-axis

    z_axis = 0.0d0
    z_axis(3) = 1.0d0

    call quat_rot_vect(quat_xy, z_axis, rot_z_axis)
    call findquat(rot_z_axis, z, quat_z)

! combine the two quaternions

    call quat_mul(quat_z, quat_xy, quat)

! We have calculated quaternion from the idealised SRF to inertial coordinate system
! To calculate the quaternion describing the opposite rotation take the conjugate

!    call quat_conj(quat)

    return
    end
