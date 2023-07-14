      SUBROUTINE GAULEG(X1,X2,X,W,N)

!c This subroutine uses Legendre polynomials to calculate the abscissa and
!c weights so that an integration can be performed over the interval X1, X2
!c using N integration points.

! x1, x2  : define range over which abscissa should be computed
! X       : output values
! W       : output weights
! N       : max degree of legendre polynomial

!c This subroutine is taken verbatim from Numerical Recipes page 145

      use gm_mod
      implicit none
!      IMPLICIT REAL*8 (A-H, O-Z)
      integer*4 :: N, M, I, J
      real*8 x(n), w(n), x1, x2
      REAL*8, PARAMETER :: EPS = 3.D-14
      real*8 :: XM, XL, Z, P1, P2, P3, PP, Z1

!      pi = 4.d0*datan(1.d0)
      M = (N + 1)/2
      XM = dble((X2 + X1)/2)
      XL = dble((X2 - X1)/2)
      DO 12 I = 1, M
        Z = dCOS(pi * (dble(I) - .25D0)/(dble(N) + .5D0))

!c first guess as to the location of the I-th pole of the Legendre polynomial
!c of degree N

1       CONTINUE
        P1 = 1.D0
        P2 = 0.D0
        DO 11 J = 1, N
          P3 = P2
          P2 = P1
          P1 = ((2.D0 * J - 1.D0) * Z * P2 - (J - 1.D0) * P3)/dble(J)
11        CONTINUE

!c recursion formula for Legendre polynomials of degree 1 up to N

        PP = dble(N) * (Z * P1 - P2)/(Z*Z - 1.D0)

!c recursion formula for the derivative of the Legendre polynomial of degree N

        Z1 = Z
        Z = Z1 - P1/PP

!c iterate the location of the current pole using Newton-Raphson

        IF ( ABS(Z - Z1) .GT. EPS ) GO TO 1

!c If necessary do further iterations

        X(I) = dble(XM - XL * Z)
        X(N + 1 - I) = dble(XM + XL * Z)
        W(I) = dble(2 * XL/((1 - Z * Z) * PP * PP))
        W(N + 1 - I) = W(I)

!c calculate corresponding weights
12      CONTINUE

      RETURN
      END
