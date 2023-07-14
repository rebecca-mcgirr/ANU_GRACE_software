    subroutine findquat(u, v, quat)

! Tony Purcell
! April 10, 2013
! This subroutine returns the quaternion defining a rotation from vector u to vector v

    implicit none

    real (kind=8), dimension(3) :: u, v, cross_prod
    real (kind=8), dimension(0:3) :: quat 
    real (kind=8) :: mag_u, mag_v, mag_factor, mag_quat, mag_cross, dot_prod
    real (kind=8) :: dot
    integer :: i, k
    real*8 :: amag3

    external cross,dot

    quat = 0.0d0

    u = u/amag3(u)
    v = v/amag3(v)
    mag_u = amag3(u)
    mag_v = amag3(v)
    mag_factor = 1.0d0/(mag_u * mag_v)
!print*,"mag_u,mag_v,mag_factor",mag_u,mag_v,mag_factor

! Check that neither of the input vectors is small enough to cause numerical problems
! if either of them is, return a warning message and a trivial rotation of 0 radians

    if ( min(mag_u, mag_v) .lt. 1.0d-16 ) then
      write(6, *) "WARNING: Subroutine findquat - zero vector input"
      quat(0) = 1.0d0
      return
    endif

! Calculate cross product of the unit vectors and store in quaternion vector entries to define
! the axis of rotation

!    call cross(u, v, cross_prod)
    call cross(u, v, cross_prod)
!print*,'after cross, U',u,amag3(u)
!print*,'after cross, V', v,amag3(v)
!print*,'after cross, U x V', cross_prod,amag3(cross_prod)

    do i = 1, 3
      quat(i) = mag_factor * cross_prod(i)
    enddo

!    quat(1) = mag_factor * (u(2) * v(3) - u(3) * v(2))
!    quat(2) = mag_factor * (u(3) * v(1) - u(1) * v(3))
!    quat(3) = mag_factor * (u(1) * v(2) - u(2) * v(1))

! test if vectors are parallel, if so return warning message and a trivial rotation of 0 radians

    mag_cross = dsqrt(quat(1)**2 + quat(2)**2 + quat(3)**2)

    dot_prod = dot(u, v)
!print*,'mag_cross and dot_prod',mag_cross,dot_prod,mag_u,mag_v,mag_factor
!print*,'cross_prod:',cross_prod
    if ( mag_cross .lt. 1.0d-15 ) then
      quat = 0.0d0
      !write(6, *) "WARNING: Subroutine findquat - parallel vectors input"
      if ( dot_prod .gt. 0.0d0) then
        quat(0) = 1.0d0
      else
        quat(1) = 1.0d0
      endif
      return
    endif

! set the magnitude of the scalar entry to match that of the cross product calculated above
! Note that by incorporating the mag_factor term we have effectively used unit vectors for
! the cross and dot products so that:
! scale factor = sin(theta)/sin(theta/2) = 2 * cos(theta/2)
! Thus we have:
! quat(0) = 2 * (cos(theta/2))^2 = 1 + cos(theta) = 1 + mag_factor * u.v

!    quat(0) = 1.0d0 + mag_factor * (u(1) * v(1) + u(2) * v(2) + u(3) * v(3))

    quat(0) = 1.0d0 + mag_factor * dot(u, v)

! normalise quaternion entries
    mag_quat = dsqrt(quat(0)**2 + quat(1)**2 + quat(2)**2 + quat(3)**2)
    do k = 0,3
      quat(k) = quat(k)/mag_quat
    enddo
!print*,'findquat: quat = ',quat

    return
    end
