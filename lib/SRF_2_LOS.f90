    subroutine SRF_2_LOS(debug,u, v, quat, flag)

! Tony Purcell April 11 2013
! This subroutine returns the quaternion defining a rotation from an idealised satellite reference frame (SRF)
! to an inertial reference frame. In the idealised satellite reference frame the satellites are assumed to be 
! perfectly aligned with the line of sight vector between them with zero yaw and roll

    implicit none

    logical      , intent(in)                  :: debug     ! flag to output (or not) debug in this subroutine
    real (kind=8), intent(inout), dimension(3) :: u, v      ! the position vectors of the two spacecraft in EF coords
    real (kind=8), intent(out), dimension(0:3) :: quat      ! the quaternion describing the rotation from SRF to EF
    integer*4,     intent(in)                  :: flag      ! indicates either SRF to LOS (1) or LOS to SRF (-1)

    real (kind=8), dimension(3) :: x, y, z, LOS, x_axis, y_axis, z_axis
    real (kind=8), dimension(3) :: rot_y_axis, rot_z_axis
    real (kind=8), dimension(0:3) :: quat_x, quat_y, quat_z, quat_xy, quat_tmp

    real (kind=8) :: mag_u, mag_v, mag_factor, mag_quat, mag_cross, mag_LOS, mag_y, mag_z
    integer :: k
    real*8  :: amag3,dot

! DEBUG
    real*8 :: u_tmp(3),v_tmp(3)

    quat = 0.0d0

    mag_u = dsqrt(u(1)**2 + u(2)**2 + u(3)**2)
    mag_v = dsqrt(v(1)**2 + v(2)**2 + v(3)**2)
 
! Check that neither of the input vectors is small enough to cause numerical problems
! if either of them is, return an error message and stop run

    if ( min(mag_u, mag_v) .lt. 1.0d-16 ) then
      write(6, *) "ERROR: Subroutine SRF_2_LOS - zero vector input"
      stop
    endif

! convert to unit vectors
! PT200724: do NOT normalise these - it introduces rounding errors that make the pitch incorrect!
!    u = u/mag_u
!    v = v/mag_v

! calculate unit line of sight vector (corresponds to x-axis in idealised SRF
    LOS = v - u 
    mag_LOS = dsqrt(LOS(1)**2 + LOS(2)**2 + LOS(3)**2)

! check that the difference between the vectors is non-zero
! if not, return an error message and stop run
    if ( mag_LOS .lt. 1.0d-16 ) then
      write(6, *) "ERROR: Subroutine SRF_2_LOS - identical vectors input",v,u
      stop
    endif

    LOS = LOS/mag_LOS
    if(debug)print*,'LOS = ',LOS,amag3(LOS)
    if(debug)print*,'u   = ',u  ,amag3(u  )
    if(debug)print*,'v   = ',v  ,amag3(v  )

! Calculate cross product of the unit vectors and store in quaternion vector entries to define
! the axis of rotation. Order here is chosen so that the z-axis will point toward the satellite's nadir
! (that is, toward the Earth)
    call cross(LOS, u, y)
    if(debug)print*,'LOS x u',y/amag3(y)  , amag3(y)



! test if vectors are parallel, if so return an error message and stop run
    mag_y = amag3(y)

    if ( mag_y .lt. 1.0d-16 ) then
      write(6, *) "ERROR: Subroutine SRF_2_LOS - parallel vectors input"
      stop
    endif

!    write(6,*) "SRF_2_LOS: y vector:", y
    y = y/mag_y

! calculate unit z-axis vector (obtained from LOS x y)

    call cross(LOS, y, z)
    if(debug)print*,'LOS x y',z/amag3(z)  , amag3(z)
    mag_z = amag3(z)
    z = z/mag_z

! Find quaternion that rotates the SRF x-axis onto the line of sight vector
    x_axis = 0.0d0
    x_axis(1) = 1.0d0

!    call findquat(x_axis, LOS, quat_x)
! PT191108: because x_axis is (1 0 0), quat_x can be written explicitly without computation:
    quat_x(1) =  0.d0
    quat_x(2) = -1.d0*LOS(3)
    quat_x(3) =       LOS(2)
    quat_x(0) = 1.d0 + dot(x_axis,LOS)
    mag_quat = dsqrt(quat_x(0)**2 + quat_x(1)**2 + quat_x(2)**2 + quat_x(3)**2)
    quat_x = quat_x/mag_quat
    if(debug)print*,'quat-x',quat_x

! Now find quaternion describing rotation from rotated y-axis to EF y-axis
    y_axis = 0.0d0
    y_axis(2) = 1.0d0
    rot_y_axis = 0.d0

    call quat_rot_vect(quat_x, y_axis, rot_y_axis)
    if(debug)call printmat_R8(rot_y_axis,3,1,"rot_y_axis ")
    call findquat(rot_y_axis, y, quat_y)
    if(debug)print*,'quat-y',quat_y

! combine the two quaternions
    call quat_mul(quat_y, quat_x, quat_xy)
    if(debug)print*,'quat-xy',quat_xy

! Now find quaternion describing rotation from rotated z-axis to EF z-axis
    z_axis = 0.0d0
    z_axis(3) = 1.0d0
    rot_z_axis = 0.d0

    call quat_rot_vect(quat_xy, z_axis, rot_z_axis)
    if(debug)call printmat_R8(rot_z_axis,3,1,"rot_z_axis ")
    call findquat(rot_z_axis, z, quat_z)
    if(debug)print*,'quat-z',quat_z

! combine the two quaternions
    call quat_mul(quat_z, quat_xy, quat)

    if ( flag .eq. -1 ) then

! We have calculated quaternion from the idealised SRF to inertial coordinate system
! To calculate the quaternion describing the opposite rotation take the conjugate
      call quat_conj(quat)
    else
      if ( flag .ne. 1 ) write(6,*) "WARNING: Subroutine SRF_2_LOS - flag should be 1 or -1"
    endif

    return
    end
