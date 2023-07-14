    subroutine srf_2_los_quadprec(debug,u, v, quat, flag)

! Tony Purcell April 11 2013
! This subroutine returns the quaternion defining a rotation from an idealised satellite reference frame (SRF)
! to an inertial reference frame. In the idealised satellite reference frame the satellites are assumed to be 
! perfectly aligned with the line of sight vector between them with zero yaw and roll
!
! PT140826: made some of this double precision to see why we are getting some errors somewhere ....

    implicit none

    logical :: debug                               ! flag to output (or not) debug in this subroutine
    real (kind=8), dimension(3) :: u, v            ! the position vectors of the two spacecraft in EF coords
    real (kind=16), dimension(0:3) :: quat         ! the quaternion describing the rotation from SRF to EF
    integer :: flag                                ! indicates either SRF to LOS (1) or LOS to SRF (-1)

    real (kind=16), dimension(3) :: x, y, z, LOS, x_axis, y_axis, z_axis
    real (kind=16), dimension(3) :: rot_y_axis, rot_z_axis
    real (kind=16), dimension(0:3) :: quat_x, quat_y, quat_z, quat_xy, quat_tmp

    real (kind=16) :: mag_u, mag_v, mag_factor, mag_quat, mag_cross, mag_LOS, mag_y, mag_z
    integer :: k
    real(kind=8)  :: amag3
    real(kind=16)  :: amag3_quadprec
    real(kind=16)  :: u_quad(3),v_quad(3),tmpvec_quad(3)

    quat = 0.0

    u_quad = u
    v_quad = v
    mag_u = amag3_quadprec(u_quad)
    mag_v = amag3_quadprec(v_quad)
 
! Check that neither of the input vectors is small enough to cause numerical problems
! if either of them is, return an error message and stop run

    if ( min(mag_u, mag_v) .lt. 1.0e-16 ) then
      call status_update('FATAL','LIB','SRF_2_LOS_quadprec',' ', "zero vector input",0)
    endif

! convert to unit vectors
!    u_quad = u_quad/mag_u
!    v_quad = v_quad/mag_v
!    print*,'u_quad =',u_quad,amag3_quadprec(u_quad)
!    print*,'v_quad =',v_quad,amag3_quadprec(v_quad)

! calculate unit line of sight vector (corresponds to x-axis in idealised SRF
    LOS = v_quad - u_quad 
    mag_LOS = amag3_quadprec(LOS)

! check that the difference between the vectors is non-zero
! if not, return an error message and stop run
    if ( mag_LOS .lt. 1.0e-16 ) then
      write(6, *) "ERROR: Subroutine SRF_2_LOS - identical vectors input",v,u
      stop
    endif

    LOS = LOS/mag_LOS
    if(debug)print*,'LOS = ',LOS,amag3_quadprec(LOS)
    if(debug)print*,'u   = ',u  ,amag3_quadprec(u_quad  )
    if(debug)print*,'v   = ',v  ,amag3_quadprec(v_quad  )

! Calculate cross product of the unit vectors and store in quaternion vector entries to define
! the axis of rotation. Order here is chosen so that the z-axis will point toward the satellite's nadir
! (that is, toward the Earth)
    call cross_quadprec(LOS, u_quad, y)
    if(debug)print*,'LOS x u',y  , amag3_quadprec(y)



! test if vectors are parallel, if so return an error message and stop run
    mag_y = amag3_quadprec(y)

    if ( mag_y .lt. 1.0e-16 ) then
      write(6, *) "ERROR: Subroutine SRF_2_LOS - parallel vectors input"
      stop
    endif

!    write(6,*) "SRF_2_LOS: y vector:", y
    y = y/mag_y
    if(debug)print*,'y-axis',y  

! calculate unit z-axis vector (obtained from LOS x y)

    call cross_quadprec(LOS, y, z)
    if(debug)print*,'LOS x y',z  , amag3_quadprec(z)
    mag_z = amag3_quadprec(z)
    z = z/mag_z
    if(debug)print*,'z-axis',z  

! Find quaternion that rotates the SRF x-axis onto the line of sight vector
    x_axis = 0.0e0
    x_axis(1) = 1.0e0

    call findquat_quad(x_axis, LOS, quat_x)
    if(debug)print*,'quat-x',quat_x

! Now find quaternion describing rotation from rotated y-axis to EF y-axis
    y_axis = 0.0d0
    y_axis(2) = 1.0d0

    call quat_rot_vect_quad(quat_x, y_axis, rot_y_axis)
    call findquat_quad(rot_y_axis, y, quat_y)
    if(debug)print*,'quat-y',quat_y

! combine the two quaternions
    call quat_mul_quad(quat_y, quat_x, quat_xy)
    if(debug)print*,'quat-xy',quat_xy

! Now find quaternion describing rotation from rotated z-axis to EF z-axis
    z_axis = 0.0d0
    z_axis(3) = 1.0d0

    call quat_rot_vect_quad(quat_xy, z_axis, rot_z_axis)
    call findquat_quad(rot_z_axis, z, quat_z)

! combine the two quaternions
    call quat_mul_quad(quat_z, quat_xy, quat)
    if(debug)then
      print*,'quat  ',quat
    endif


    if ( flag .eq. -1 ) then

! We have calculated quaternion from the idealised SRF to inertial coordinate system
! To calculate the quaternion describing the opposite rotation take the conjugate
      call quat_conj_quad(quat)
    else
      if ( flag .ne. 1 ) write(6,*) "WARNING: Subroutine SRF_2_LOS - flag should be 1 or -1"
    endif

    return
    end
