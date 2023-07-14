    subroutine calpnm(xx,pnms)
!c
!c     calculate the associated Legendre function values at x up to degree nmax
!c
!c     The normalisation is such that the integral of Pnm^2 from -1 to 1 is 2.
!c     This is slightly different from the normalisation used in the older
!c     sea-level programs.  In those programs, if m is not equal to 0, the
!c     normalisation had an extra factor sqrt(2), so that integral over the
!c     sphere of Pnm cos or sin m phi is 4 pi.  With this normalisation
!c     the integral over the sphere of Pnm exp(i m phi) is equal to 4 pi.
!c
    use gm_mod
    use gauleg_mod
    implicit none
    real*8 aa, as, s, sum, c, d, x
    real*8 :: pnm(maxnm), pnms(maxnm)
    real*8 :: xx
    integer*4 :: i, i0, i1, i2, n, n0, m, j

      x = dble(xx)
      pnm(1) = 1.d0
      pnm(2) = x
      i1 = 1
      i2 = 2
      do 10 i = 2, nmax
        i0 = i1
        i1 = i2
        i2 = i * (i + 1)/2 + 1
        pnm(i2) = ((2.d0 * dble(i) - 1.d0) * x * pnm(i1) - (dble(i) - 1.d0) * pnm(i0))/dble(i)
10      continue

!c classic recurrence relation for order 0 Legendre polynomials

      do 20 i = 1, nmax
        i2 = i * (i + 1)/2 + 1
        pnm(i2) = dsqrt(2.d0 * dble(i) + 1.d0) * pnm(i2)
20      continue

!c apply normalisation as discussed above

      AA = 1 - x * x
      IF ( AA .NE. 0 ) GO TO 100
      DO 200 n = 0, nmax
        n0 = n * (n + 1)/2 + 1
        do 200 m = 1, n
200       pnm(m + n0) = 0.d0
      RETURN

!c at the poles set the associated Legendre polynomials identically to zero

100   CONTINUE
      AS = dSQRT(AA)
      pnm(3) = AS * dSQRT(1.5D0)
      pnm(5) = 3 * AS * x * dSQRT(2.5D0/3.d0)
      pnm(6) = 3 * AA * dSQRT(5.D0/24.d0)

!c calculate associated Legendre polynomials for degrees 1 and 2

      DO 70 N = 3, nmax
        SUM = 1.d0
        DO 72 J = 1, N
          S = (2.d0 * dble(J) - 1.d0)/(2.d0 * dble(J))
72        SUM = SUM * dSQRT(AA* S)
        n0 = n * (n + 1)/2 + 1
        pnm(n + n0) = SUM * dSQRT(2.D0 * dble(N) + 1.d0)

!c calculate order n associated Legendre polynomial for current degree

        SUM = 1.d0
        DO 74 J = 2, N
          S = (2.d0 * dble(J) - 1.d0)/(2.d0 * dble(J) - 2.d0)
74        SUM = SUM * dSQRT(AA * S)
        pnm(n - 1 + n0) = x * dSQRT(2.D0 * dble(N) + 1.d0) * SUM

!c calculate order n - 1 associated Legendre polynomial for current degree

        DO 76 M = n - 2, 1, -1
          C = 2.d0 * (dble(M) + 1.d0)/dSQRT((dble(N) + dble(M) + 1.D0) * (dble(N) - dble(M)))
          C = C * x/AS
          D = (dble(N) + dble(M) + 2.d0) * (dble(N) - dble(M) - 1.d0)/((dble(N) + dble(M) + 1.D0) * (dble(N) - dble(M)))
          D = dSQRT(D)
76        pnm(M + n0) = C * pnm(M + 1 + n0) - D * pnm(M + 2 + n0)
70        CONTINUE

!c calculate associated Legendre polynomials of order 1, ..., n - 2 for the
!c current degree

      do i = 1, maxnm
        pnms(i) = dble(pnm(i))
      enddo

      return
      end
