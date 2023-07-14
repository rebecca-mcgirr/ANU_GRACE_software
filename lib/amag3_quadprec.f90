      function amag3_quadprec(a)
!
! computes the magnitude of a vector using a dot product
! P Tregoning 3/95
!
! PT140826: turned into a quadruple precision fortran90 subroutine

      implicit none
      real(kind=16) :: amag3_quadprec,a(3)

      amag3_quadprec = sqrt(a(1)*a(1) +a(2)*a(2) + a(3)*a(3))

      return
      end
