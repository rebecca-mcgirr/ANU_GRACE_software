    subroutine cart_convert(lat, long, sp, cart)

! This subroutine carries out a coordinate transformation on a vector
! in spherical coordinates (sp) into cartesian coordinates (cart)

    real*8, dimension(3,3) :: T
    real*8, dimension(3) :: cart, sp
    real*8 :: pi, colat, lat, long

    pi = 4.0d0*datan(1.d0)

    colat =pi/2.d0 - lat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The elements of the transformation matrix from spherical to cartesian
   T(1,1) = dcos(lat)*dcos(long)
   T(1,2) = -1.d0*dsin(lat)*dcos(long)
   T(1,3) = -1.d0*dsin(long)
   T(2,1) = dcos(lat)*dsin(long)
   T(2,2) = -1.d0*dsin(lat)*dsin(long)
   T(2,3) = dcos(long)
   T(3,1) = dsin(lat)
   T(3,2) = dcos(lat)
   T(3,3) = 0.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO i = 1,3
      cart(i) = 0.d0
      DO j = 1,3
        cart(i) = cart(i) + T(i,j)*sp(j)
      ENDDO
    ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    return
    END



