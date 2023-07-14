      SUBROUTINE nutred_grace( INUT,FJD,DPSI,DEPS )
C
C     Read nutation values from an external data set
C     R.I. Abbot - NOVEMBER 1984
C     MODEL version into library by R. King  January 1992
C
c PT140220: modified to make use of the format statement that is written into the second header line!!!

      implicit none

      integer*4 inut,itsav,nvtot,it,nr,nrec,npr,nvr
     .        , i,j,k,n,len,rcpar

      REAL*8 JD1,JD2,JDT1,JDT2,fjd
     .     , TAB(4),NUTVAL(8,2),Y1(2,2),Y2(2,2)
     .     , fjdf,fjdn,deps,dpsi,rint,rmult,s,t,rjd,f2

      character*80 prog_name
      character*256 message
      character*50 afmt

      logical first
C
      SAVE NUTVAL,JDT1,JDT2,JD1,JD2,Y1,Y2,ITSAV,NVTOT,NVR,NREC,FIRST
     .      ,rint

c      DATA NPR,RINT/4,0.5D0/,first/.true./
      data first,rmult/.true.,1.d-4/

      FJDN = DINT(FJD)
      FJDF = FJD - FJDN
C
      IF (.not.first) GO TO 130
C
C        Read the header records
C
      REWIND INUT
c      READ (INUT,50,end=520) JDT1,JDT2
c   50 FORMAT (/,35X,F7.0,1X,F7.0)
      READ (inut,50) afmt,JDT1,JDT2,NPR,rint,rmult
   50 format (/,a29,6x,F7.0,1X,F7.0,1x,i2,f3.1,9x,e7.4)
      first = .false.
      NVTOT = 0
      NREC=0
      JDT2 = JDT2 - RINT
      ITSAV=-9999
C
C
C        Read the data record into storage
C
  100 IF (NVTOT.GT.8-NPR) GO TO 120
C        Keep reading until (8/NPR) records are in storage
      READ(INUT,fmt=afmt,end=530 )RJD
     .           ,((NUTVAL(NVTOT+I,K),K=1,2),I=1,NPR),NVR
c      print*,rjd,((NUTVAL(NVTOT+I,k),K=1,2),I=1,NPR)
c      READ(INUT,110,end=530 )RJD
c     .           ,((NUTVAL(NVTOT+I,K),K=1,2),I=1,NPR),NVR
c  110 FORMAT (1X,F5.0,1X,8(F7.0,1X),7X,I2)
      NREC = NREC + 1
C        JD1 and JD2 are the limits of usable values in storage
C        and depend on the interpolation scheme
      IF (NVTOT.EQ.0) JD1 = RJD + 2400000.D0 + RINT
      IF (NVR.EQ.0) NVR = NPR
      NVTOT = NVTOT + NVR
      ITSAV=-9999
      JD2 = RJD + 2400000.D0 + (NVR-2)*RINT
      IF (JD2.LT.JDT2) GO TO 100
C        Save the first data just in case the header JDT1 is wrong
  120 IF ( first ) JDT1 = JD1
C
C
C        Is JD within range of the table ?
C
  130 IF (FJDN.LT.JDT1.OR.FJDN.GT.JDT2) GO TO 500
C
C        Is JD too close to JD1?
C*      IF (FJDN.GE.JD1) GO TO 170
      IF (FJD.GE.JD1) GO TO 170
C        If so, backspace and try again
C**      N = (JD1-FJDN)/RINT/NPR
c             round up any fractional part - rwk 92/2/27
!
      N = NINT ( (JD1-FJD)/RINT/NPR + .5d0 ) + 2
      IF (N.GT.NREC) GOTO 540
      DO 160 I=1,N
  160 BACKSPACE INUT
      NREC = NREC - N
      NVTOT = 0
      GO TO 100
C
C        Is JD too close to JD2 ?
  170 continue
      IF (FJD.LT.JD2) GO TO 200
C        If so, shift storage and read another record
      NVTOT = NVTOT - NVR
      DO   I=1,NVTOT
        DO   K=1,2
          NUTVAL(I,K) = NUTVAL(NPR+I,K)
        enddo
      enddo
      JD1 = JD1 + NPR*RINT
      GO TO 100
C
C
C        Calculate interpolation times and value of tabular points
C
  200 continue
      T = FJDN - JD1
      T = (T + FJDF)/RINT
      IT = T
      T = DMOD(T,1.D0)
      S = 1.D0 - T
      IF (IT.EQ.ITSAV) GO TO 230
      DO  K=1,2
        DO  I=1,4
          J = IT + I
          TAB(I) = NUTVAL(J,K)
        enddo
C
C
C        Calculate interpolation Y-vector
C
          DO  I=1,2
           NR = I + 1
            F2 =     0.166666666666667D0 * (TAB(NR+1)+TAB(NR-1))
            Y1(I,K) =  1.333333333333333D0 * TAB(NR) - F2
            Y2(I,K) = -0.333333333333333D0 * TAB(NR) + F2
          enddo

      enddo
      ITSAV = IT
C
C
C        Second difference interpolation
C
  230 DPSI=(T*(Y1(2,1)+T*T*Y2(2,1))+S*(Y1(1,1)+S*S*Y2(1,1)))*rmult
      DEPS=(T*(Y1(2,2)+T*T*Y2(2,2))+S*(Y1(1,2)+S*S*Y2(1,2)))*rmult
c      print*,'nutred_grace: dpsi,deps',dpsi,deps
      GOTO 999
C
C
C        Out of range of table
C
  500 WRITE(message,510) FJDN,JDT1,JDT2  
  510 FORMAT('JD= ',F8.0,' out of range of Nutation Table, JDT1= '
     .      ,F8.0,' JDT2= ',F8.0)
c     get calling module name for report_stat
      len =  rcpar(0,prog_name)
      call report_stat('FATAL',prog_name,'lib/nutred',' ',message,0)
  520 call report_stat('FATAL',prog_name,'lib/nutred',' '
     .                ,'EOF on header',0)
  530 call report_stat('FATAL',prog_name,'lib/nutred',' '
     .                ,'EOF in nutation table',0)
  540 call report_stat('FATAL',prog_name,'lib/nutred',' '
     .                ,'Backspacing ahead of start of nuttbl.',0)

  999 RETURN
      END
