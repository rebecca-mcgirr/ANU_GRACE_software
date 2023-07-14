!Copyright (c) Massachusetts Institute of Technology and the University of
!California at San Diego, 1994. All rights reserved.

      SUBROUTINE PNROT_grace(INUT,JD,T,TDTGPST,EQE,PREC,RNUT,frame,precmod)
!
! Written by Yehuda Bock (1987)
!
! Compute the precession and nutation matrices

!     Input:
!       inut       Unit number for nutation table
!       jd         PEP Julian day 
!       t          Seconds-of-day  of PEP julian day (GPST)
!       tdtgpst    TDT - GPST  (seconds)
!       frame      Inertial frame
!       precmod    precession model
!
!      Output:
!        EqE          Equation of the equinoxes (rad)
!        prec         Precession matrix (rad)
!        rnut         Nutation matrix (rad)

      implicit none

      character*5 frame, precmod

      integer*4 jd,inut

      real*8 t,tdtgpst,eqe,prec,rnut,deps,oblq,dpsi,fjdtdt
      real*8 oblm,oblt,eqeq,pi
      real*8 matrix1(3,3),matrix2(3,3),matrix3(3,3),tmp(3,3)
! precession variables
      real*8 R_pb(3,3),Xr,Yr,Zr,X,Y,Z,phi,S,C,F &
             ,tt,t1,t2,t3,t4,t5,oblq_pt,casr &
             ,pzeta(3,3),ptheta(3,3),pz(3,3),prec_pt(3,3),work(3,3)
! nutation variables
      real*8 dpsi_pt,deps_pt,zeta,theta,RNA(3,3),RNB(3,3),RNC(3,3)

      dimension prec(3,3),rnut(3,3)

      pi = 4.d0*datan(1.d0)
      casr = pi/180.d0/3600.d0

! Convert GPS time to Terrestrial Dynamical Time (still PEP julian day)

      FJDTDT= DBLE(JD) + T/86400.D0 + TDTGPST/86400.d0 

! Form the precession matrix
! PT140219: replace old gamit call with my own computation of the precession matrix
!      call prces(FJDTDT,OBLQ,PREC,frame,precmod)

! PT140218: code for zeta theta and z. Ref: Capitaine et al, Astronomy and Astrophysics, 412, 567–586 (2003)
      tt = (dble(JD) +T/86400.d0 - 0.5D0 - 2451545.d0)/36525.d0
      t1 = tt
      t2 = tt*tt
      t3 = t2*tt
      t4 = t3*tt
      t5 = t4*tt

! zeta, theta and z are in arcseconds. Multiplying by casr (which is pi/180/3600) converts to radians.
      zeta =    2.650545d0    + 2306.0832270d0*t1 +0.29884990d0*t2 +0.018018280d0*t3 -0.0000059710d0*t4 -0.0000003173d0*t5
      theta= 2004.191903d0*t1     -0.4294934d0*t2 -0.04182264d0*t3 -0.000007089d0*t4 -0.0000001274d0*t5
      z    =    2.650545d0    + 2306.077181d00*t1 +1.0927348d00*t2 +0.018268370d0*t3 -0.0000285960d0*t4 -0.0000002904d0*t5

!  Calculate the precession matrix from the reference to the current epoch
      Call ROTMAT(-zeta*casr,3,pzeta)
      Call ROTMAT(theta*casr,2,ptheta)
      Call ROTMAT(-z*casr,3,pz)
      Call MATMPY(ptheta,pzeta,work,3,3,3)
      Call MATMPY(pz,work,prec_pt,3,3,3)

! PT140217: compute my own obliquity of the ecliptic (Capitaine et al, 2003)
      oblq_pt = (84381.406d0-46.836769d0*t1-0.0001831d0*t2+0.00200340d0*t3 - 0.000000576d0*t4 - 0.0000000434d0*t5) *casr
       prec=prec_pt
       OBLQ = oblq_pt

! Form the nutation matrix
! PEP JD is required as time argument to nutred
      CALL nutred_grace( INUT,FJDTDT,DPSI,DEPS )
      DPSI=DPSI*CASR
      DEPS=DEPS*CASR

! Calculate the nutation matrix (Reference:  Mueller, p. 75)
      Call ROTMAT(OBLQ,1,RNA)
      Call ROTMAT(-DPSI,3,RNB)
      Call ROTMAT(-OBLQ-DEPS,1,RNC)
      Call MATMPY(RNB,RNA,WORK,3,3,3)
      Call MATMPY(RNC,WORK,RNUT,3,3,3)

!      call printmat(rnut,3,3,'gamit_rnut ')

! Compute the equation of the equinoxes for computation of GAST
      EqE=dpsi*dcos(oblq)

      return
      end
