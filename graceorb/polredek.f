Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      SUBROUTINE POLREDEK( IPOLE,JD,FRACT,XPOLE,YPOLE,XPDOT,YPDOT )
C
C     Read pole values from an external data set

      implicit none

      INTEGER*4 JD,MJD,JD1,JD2,JDT1,JDT2,J00001,ipolvl(12,2),itsav
     .        , ipole,nvtot,it,nr,nrec,nrecs,ifirst,int,npr,nvr
     .        , nback,nrbmin,lap,i,j,k,n,len,rcpar


      character*80 prog_name
      character*40   afmt
      character*256  message

      REAL*8 TAB(4),POLVAL(12,2),Y1(2,2),Y2(2,2),t,s,fract,rmult,f2
     .     , xpole,ypole,xpdot,ypdot

      SAVE POLVAL,JDT1,JDT2,JD1,JD2,Y1,Y2,
     $     ITSAV,NVTOT,NVR,NREC,IFIRST,npr,int,rmult,afmt

C     DATA NPR,INT/6,5/
      DATA J00001/2400001/
      nrecs=0
      nback=0
      nrbmin=0
      lap=0

c     get calling module name for report_stat
      len = rcpar(0,prog_name)

      IF (IFIRST.NE.0) GO TO 130
C
C
C       Read the header records

      REWIND IPOLE
c     READ (IPOLE,50) afmt,JDT1,JDT2
      READ (ipole,50) afmt,JDT1,JDT2,npr,int,rmult

c  50 FORMAT (/,a35,I7,1X,I7)
   50 FORMAT (/,a35,i7,1X,i7,1x,i2,1x,i2,9x,g7.0)


c     find the end of the format statement
      i = index(afmt,' ')
      afmt = afmt(1:i-1)

c     choke if you can't find a format statement
      if (i .lt. 1) goto 540

      IFIRST = 1
      NVTOT = 0
      NREC=0
      nrbmin= 12/npr + 1
      JDT2 = JDT2 - INT
      ITSAV=-9999
C
C       Read the data records into storage
C       Keep reading until (12/NPR) records are in storage
c       read values as integers, convert them to reals

  100 IF (NVTOT.GT.12-NPR) GO TO 120
       READ(IPOLE,fmt=afmt,err=520,END=520)
     .   MJD,((IPOLVL(I,K),K=1,2),I=1,NPR),NVR
         do 117 i=1,npr
            do 115 k=1,2
               polval(nvtot+i,k) = dble(ipolvl(i,k))
  115       continue
  117    continue

       NREC = NREC + 1
C       JD1 and JD2 are the limits of usable values in storage
C       and depend on the interpolation scheme
      IF (NVTOT.EQ.0) JD1 = MJD + J00001 + INT
      IF (NVR.EQ.0) NVR = NPR
      NVTOT = NVTOT + NVR
      ITSAV=-9999
      JD2 = MJD + J00001 + (NVR-2)*INT
      IF (JD2.LT.JDT2) GO TO 100
C       Save the first date just in case the header JDT1 is wrong
  120 IF (IFIRST.EQ.0) JDT1 = JD1
C
C
C       Is JD within the range of the table ?
C
  130 IF (JD.LT.JDT1 .OR. JD.GT.JDT2) GO TO 500
C
C
C       Is JD too close to JD1 ?
C
      IF ( JD.GE.JD1 ) GO TO 170
C       If so, backspace and try again
      if( nback.gt.0 .and. nrec.ge.nrecs ) lap = nrbmin+1+nrec-nrecs
      nrecs = nrec
      nback = nback + 1
      if( nback.gt.100) then
         write(message,'(a,i8,a,2i8)')
     .     'Infinite loop reading pole. for jd=',jd,' jd1,jd2=',jd1,jd2
         call report_stat('FATAL',prog_name,'lib/polred',' ',message,0)
      endif
      N = (JD1-JD)/INT/NPR + lap
      IF (N.GT.NREC) N = NREC
      DO 160 I=1,N
  160 BACKSPACE IPOLE
      NREC = NREC - N
      NVTOT = 0
      GO TO 100
C
C       Is JD too close to JD2 ?
  170 IF ( JD.LT.JD2 ) GO TO 200
C       If so, shift storage and read another record
      NVTOT = NVTOT - NVR
      DO 180 I=1,NVTOT
      DO 180 K=1,2
  180 POLVAL(I,K) = POLVAL(NPR+I,K)
      JD1 = JD1 + NPR*INT
      GO TO 100
C
C
C       Calculate interpolation times and value of tabular points
C
  200 T = JD - JD1
      T = (T + FRACT)/INT
      IT = T
      T = DMOD(T,1.D0)
      S = 1.D0 - T
      IF (IT.EQ.ITSAV) GO TO 230
      DO 225 K=1,2
      DO 210 I=1,4
      J = IT + I
  210 TAB(I) = POLVAL(J,K)
C
C
C       Calculate interpolation Y-vector
C
      DO 220 I=1,2
      NR = I + 1
      F2 =     0.166666666666667D0 * (TAB(NR+1)+TAB(NR-1))
      Y1(I,K) =  1.333333333333333D0 * TAB(NR) - F2
  220 Y2(I,K) = -0.333333333333333D0 * TAB(NR) + F2
  225 CONTINUE
      ITSAV = IT
C
C
C       Second difference interpolation
C
c 230 XPOLE=(T*(Y1(2,1)+T*T*Y2(2,1))+S*(Y1(1,1)+S*S*Y2(1,1)))*1.D-3
c     YPOLE=(T*(Y1(2,2)+T*T*Y2(2,2))+S*(Y1(1,2)+S*S*Y2(1,2)))*1.D-3
  230 XPOLE=(T*(Y1(2,1)+T*T*Y2(2,1))+S*(Y1(1,1)+S*S*Y2(1,1)))*rmult
      YPOLE=(T*(Y1(2,2)+T*T*Y2(2,2))+S*(Y1(1,2)+S*S*Y2(1,2)))*rmult
      xpdot=( y1(2,1)+3.d0*t*t*y2(2,1) - y1(1,1)-3.d0*s*s*y2(1,1) )
     .        *rmult/int
      ypdot=( y1(2,2)+3.d0*t*t*y2(2,2) - y1(1,2)-3.d0*s*s*y2(1,2) )
     .        *rmult/int
      GOTO 999
C
C
C        Out of range of table
C
  500 WRITE(message,510) JD,JDT1,JDT2
  510 FORMAT('JD= ',I7,' out of range of pole table, JDT1= '
     .,I7,' JDT2= ',I7)
      call report_stat('FATAL',prog_name,'lib/polred',' ',message,0)
  520 WRITE(message,530) MJD,JD2,JDT2
  530 FORMAT('File error in pole table, MJD = '
     .,I5,'  JD2 = ',I7,'  JDT2 = ',I7)
      call report_stat('FATAL',prog_name,'lib/polred',' ',message,0)
  540 call report_stat('FATAL',prog_name,'lib/polred',' '
     .,'Cannot find format statement in pole file:',0)

  999 RETURN
      END
