      SUBROUTINE CROSS_quadprec(A,B,C)
!
! COMPUTE CROSS PRODUCT OF TWO VECTORS A AND B WITH RESULTS IN C
! S. A. GOUREVITCH	6/81
!
! PT140826: turned into a quadruple precision fortran90 routine

      implicit none

      real(kind=16) :: a,b,c
      DIMENSION A(3),B(3),C(3)
!
      C(1)=A(2)*B(3)-A(3)*B(2)
      C(2)=A(3)*B(1)-A(1)*B(3)
      C(3)=A(1)*B(2)-A(2)*B(1)
!

     RETURN
      END
