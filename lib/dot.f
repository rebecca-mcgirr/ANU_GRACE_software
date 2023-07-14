      FUNCTION DOT(A,B)
C
C COMPUTE DOT PRODUCT OF VECTORS A AND B
C S. A. GOUREVITCH	6/81
C                      
      implicit none

      real*8 a,b,dot
      DIMENSION A(3),B(3)

      DOT=A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
      RETURN
      END
