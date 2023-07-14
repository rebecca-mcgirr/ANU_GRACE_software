    subroutine findquat_quad(u, v, quat)

! Tony Purcell
! April 10, 2013
! This subroutine returns the quaternion defining a rotation from vector u to vector v

    implicit none

    real (kind=16), dimension(3) :: u, v, cross_prod
    real (kind=16), dimension(0:3) :: quat 
    real (kind=16) :: mag_u, mag_v, mag_factor, mag_quat, mag_cross, dot_prod
    real (kind=16) :: dot_quad,amag3_quadprec
    integer :: i, k

    external cross_quadprec,dot_quadprec

    quat = 0.0e0

    mag_u = sqrt(u(1)**2 + u(2)**2 + u(3)**2)
    mag_v = sqrt(v(1)**2 + v(2)**2 + v(3)**2)
    mag_factor = 1.0e0/(mag_u * mag_v)

! Check that neither of the input vectors is small enough to cause numerical problems
! if either of them is, return a warning message and a trivial rotation of 0 radians

    if ( min(mag_u, mag_v) .lt. 1.0e-16 ) then
      write(6, *) "WARNING: Subroutine findquat - zero vector input"
      quat(0) = 1.0e0
      return
    endif

! Calculate cross product of the unit vectors and store in quaternion vector entries to define
! the axis of rotation

!    call cross(u, v, cross_prod)
    call cross_quadprec(u, v, cross_prod)

    do i = 1, 3
      quat(i) = mag_factor * cross_prod(i)
    enddo


! test if vectors are parallel, if so return warning message and a trivial rotation of 0 radians
    mag_cross = amag3_quadprec(quat(1:3))

    dot_prod = dot_quad(u, v)

    if ( mag_cross .lt. 1.0e-20 ) then
      quat = 0.0e0
      write(6, *) "WARNING: Subroutine findquat - parallel vectors input"
      if ( dot_prod .gt. 0.0e0) then
        quat(0) = 1.0e0
      else
        quat(1) = 1.0e0
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

    quat(0) = 1.0e0 + mag_factor * dot_quad(u, v)

! normalise quaternion entries

    mag_quat = sqrt(quat(0)**2 + quat(1)**2 + quat(2)**2 + quat(3)**2)
    do k = 0,3
      quat(k) = quat(k)/mag_quat
    enddo

    return
    end
