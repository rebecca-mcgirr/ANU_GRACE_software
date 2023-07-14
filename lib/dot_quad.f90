      FUNCTION DOT_quad(A,B)
!
! COMPUTE DOT PRODUCT OF VECTORS A AND B
! S. A. GOUREVITCH	6/81
!                      
      implicit none

      real(kind=16) a,b,dot_quad
      DIMENSION A(3),B(3)

      DOT_quad=A(1)*B(1)+A(2)*B(2)+A(3)*B(3)

      RETURN
      END
