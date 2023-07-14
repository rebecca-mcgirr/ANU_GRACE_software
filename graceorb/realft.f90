      SUBROUTINE REALFT(DATA,N,ISIGN)

!c Subroutine to calculate FFT of real function whose entries are given in the
!c array DATA.
!c
!c The program is taken verbatim from Numerical Recipes page 507

      use gm_mod
!      dimension DATA(2 * N)
      implicit none
      integer*4 :: N, N2P3, I, I1, I2, I3, I4, isign
      real*8 :: H1R, H1I, H2R, H2I 
      real*8 :: C1, C2
      real*8 :: data(2*n)
      real*8 wr, wi, wpr, wpi, wtemp, theta

!      print*, "realft", pi
!      pi = 4.d0 * datan(1.d0)
      THETA = pi/dble(n)
      C1 = 0.5d0
      IF ( ISIGN .EQ. 1 ) THEN
        C2 = -0.5d0
        CALL FOUR1(DATA, N, +1)
      ELSE
        C2 = 0.5d0
        THETA = -THETA
      ENDIF
      WPR = -2.d0 * DSIN(0.5D0 * THETA)**2
      WPI = DSIN(THETA)
      WR = 1.d0 + WPR
      WI = WPI
      N2P3 = 2 * dble(N) + 3
      DO 11 I = 2, N/2 + 1
        I1 = 2 * I - 1
        I2 = I1 + 1
        I3 = N2P3 - I2
        I4 = I3 + 1
        H1R = C1 * (DATA(I1) + DATA(I3))
        H1I = C1 * (DATA(I2) - DATA(I4))
        H2R = -C2 * (DATA(I2) + DATA(I4))
        H2I = C2 * (DATA(I1) - DATA(I3))
        DATA(I1) = H1R + dble(wr) * H2R - dble(wi) * H2I
        DATA(I2) = H1I + dble(wr) * H2I + dble(wi) * H2R
        DATA(I3) = H1R - dble(wr) * H2R + dble(wi) * H2I
        DATA(I4) = -H1I + dble(wr) * H2I + dble(wi) * H2R
        WTEMP = WR
        WR = WR * WPR - WI * WPI + WR
        WI = WI * WPR + WTEMP * WPI + WI
11    CONTINUE
      IF ( ISIGN .EQ. 1) THEN
        H1R = DATA(1)
        DATA(1) = H1R + DATA(2)
        DATA(2) = H1R - DATA(2)
      ELSE
        H1R = DATA(1)
        DATA(1) = C1 * (H1R + DATA(2))
        DATA(2) = C1 * (H1R - DATA(2))
        CALL FOUR1(DATA, N, -1)
      ENDIF

      RETURN
      END


      SUBROUTINE FOUR1(DATA, NN, ISIGN)

!c This subroutine calculates the FFT or reverse FFT of a given function. It
!c is called from REALFT above.

!c This program is taken verbatim from Numerical Recipes page 501

!      DIMENSION DATA(2 * NN)

      implicit none
      integer*4 :: NN, N, J, I, M, MMAX, istep, isign
      real*8 :: data(2*nn)
      real*8 wr, wi, wpr, wpi, wtemp, theta
      real*8 :: tempr, tempi

      N = 2 * NN
      J = 1
      DO 11 I = 1, N, 2
        IF ( J .GT. I ) THEN
          TEMPR = DATA(J)
          TEMPI = DATA(J + 1)
          DATA(J) = DATA(I)
          DATA(J + 1) = DATA(I + 1)
          DATA(I) = TEMPR
          DATA(I + 1) = TEMPI
        ENDIF
        M = N/2
1       IF ( ( M .GE. 2 ) .AND. ( J .GT. M ) ) THEN
          J = J - M
          M = M/2
          GO TO 1
        ENDIF
        J = J + M
11    CONTINUE
      MMAX = 2
2     IF ( N .GT. MMAX ) THEN
        ISTEP = 2 * MMAX
        THETA = 6.28318530717959D0/(ISIGN * MMAX)
        WPR = -2.D0 * DSIN(0.5D0 * THETA)**2
        WPI = DSIN(THETA)
        WR = 1.D0
        WI = 0.D0
        DO 13 M = 1, MMAX, 2
          DO 12 I = M, N, ISTEP
            J = I + MMAX
            TEMPR = dble(wr) * DATA(J) - dble(wi) * DATA(J + 1)
            TEMPI = dble(wr) * DATA(J + 1) + dble(wi) * DATA(J)
            DATA(J) = DATA(I) - TEMPR
            DATA(J + 1) = DATA(I + 1) - TEMPI
            DATA(I) = DATA(I) + TEMPR
            DATA(I + 1) = DATA(I + 1) + TEMPI
12        CONTINUE
          WTEMP = WR
          WR = WR * WPR - WI * WPI + WR
          WI = WI * WPR + WTEMP * WPI + WI
13      CONTINUE
        MMAX = ISTEP
      GO TO 2
      ENDIF


      RETURN
      END
