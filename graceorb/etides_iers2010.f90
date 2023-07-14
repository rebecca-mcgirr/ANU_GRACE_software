  subroutine etides_iers2010

! define all the in-phase and out-of-phase tidal amplitudes pertaining to the effect on low-degree
! coefficients for the solid body tide computations
!
! all coefficients are defined in arrays declared in etide_mod. Nothing passed in/out as arguments
!
! P. Tregoning
! 29 October 2018

  use etides_mod

  implicit none


! allocate the arrays for deg 20
  num_deg20 = 21
  allocate(ztid_doodson(num_deg20))   ! Doodson number
  allocate(ztidip(num_deg20))         ! deg 2 in-phase amplitudes
  allocate(ztidop(num_deg20))         ! deg 2 out-of-phase amplitudes

! allocate the arrays for deg 21
  num_deg21 = 48
  allocate(dtid_doodson(num_deg21))   ! Doodson number
  allocate(dtidip(num_deg21))         ! deg 2 in-phase amplitudes
  allocate(dtidop(num_deg21))         ! deg 2 out-of-phase amplitudes

! allocate the arrays for deg 22
  num_deg22 = 2
  allocate(sdtid_doodson(num_deg21))   ! Doodson number
  allocate(sdtidip(num_deg22))         ! deg 2 in-phase amplitudes
  allocate(sdtidop(num_deg22))         ! deg 2 out-of-phase amplitudes


! corrections to C20
! unnamed 
  ztid_doodson(1) = ' 55.565'
  ztidip(1) = 16.6d-12
  ztidop(1) = -6.7d-12
! unnamed
  ztid_doodson(2) = ' 55.575'
  ztidip(2) = -0.1d-12
  ztidop(2) =  0.1d-12
! Sa
  ztid_doodson(3) = ' 56.554'  
  ztidip(3) = -1.2d-12
  ztidop(3) =  0.8d-12
! Ssa
  ztid_doodson(4) = ' 57.555'  
  ztidip(4) = -5.5d-12
  ztidop(4) =  4.3d-12
! unnamed
  ztid_doodson(5) = ' 57.565'  
  ztidip(5) =  0.1d-12
  ztidop(5) = -0.1d-12
! unnamed
  ztid_doodson(6) = ' 58.554'  
  ztidip(6) = -0.3d-12
  ztidop(6) =  0.2d-12
! Msm
  ztid_doodson(7) = ' 63.655'  
  ztidip(7) = -0.3d-12
  ztidop(7) =  0.7d-12
! unnamed
  ztid_doodson(8) = ' 65.445'  
  ztidip(8) =  0.1d-12
  ztidop(8) = -0.2d-12
! Mm
  ztid_doodson(9) = ' 65.455'  
  ztidip(9) = -1.2d-12
  ztidop(9) =  3.7d-12
! unnamed
  ztid_doodson(10) = ' 65.465'  
  ztidip(10) =  0.1d-12
  ztidop(10) = -0.2d-12
! unnamed
  ztid_doodson(11) = ' 65.655'  
  ztidip(11) =  0.1d-12
  ztidop(11) = -0.2d-12
! Msf
  ztid_doodson(12) = ' 73.555'  
  ztidip(12) =  0.0d-12
  ztidop(12) =  0.6d-12
! unnamed
  ztid_doodson(13) = ' 75.355'  
  ztidip(13) =  0.0d-12
  ztidop(13) =  0.3d-12
! Mf
  ztid_doodson(14) = ' 75.555'  
  ztidip(14) =  0.6d-12
  ztidop(14) =  6.3d-12
! unnamed
  ztid_doodson(15) = ' 75.565'  
  ztidip(15) =  0.2d-12
  ztidop(15) =  2.6d-12
! unnamed
  ztid_doodson(16) = ' 75.575'  
  ztidip(16) =  0.0d-12
  ztidop(16) =  0.2d-12
! Mstm
  ztid_doodson(17) = ' 83.655'  
  ztidip(17) =  0.1d-12
  ztidop(17) =  0.2d-12
! Mtm
  ztid_doodson(18) = ' 85.455'  
  ztidip(18) =  0.4d-12
  ztidop(18) =  1.1d-12
! unnamed
  ztid_doodson(19) = ' 85.465'  
  ztidip(19) =  0.2d-12
  ztidop(19) =  0.5d-12
! Msqm
  ztid_doodson(20) = ' 93.555'  
  ztidip(20) =  0.1d-12
  ztidop(20) =  0.2d-12
! Mqm
  ztid_doodson(21) = ' 95.355'  
  ztidip(21) =  0.1d-12
  ztidop(21) =  0.1d-12

! Corrections to C21and S21 (diurnal)
! 2Q1
  dtid_doodson(1) = '125.755'
  dtidip(1)       =   -0.1e-12
  dtidop(1)       =    0.0e-12
! sigma_1
  dtid_doodson(2) = '127.555'
  dtidip(2)       =   -0.1e-12
  dtidop(2)       =    0.0e-12
! unnamed
  dtid_doodson(3) = '135.645'
  dtidip(3)       =   -0.1e-12
  dtidop(3)       =    0.0e-12
! Q1
  dtid_doodson(4) = '135.655'
  dtidip(4)       =   -0.7e-12
  dtidop(4)       =    0.1e-12
! rho1
  dtid_doodson(5) = '137.455'
  dtidip(5)       =   -0.1e-12
  dtidop(5)       =    0.0e-12
! unnamed
  dtid_doodson(6) = '145.545'
  dtidip(6)       =   -1.3e-12
  dtidop(6)       =    0.1e-12
! O1
  dtid_doodson(7) = '145.555'
  dtidip(7)       =   -6.8e-12
  dtidop(7)       =    0.6e-12
! tau1
  dtid_doodson(8) = '147.555'
  dtidip(8)       =    0.1e-12
  dtidop(8)       =    0.0e-12
! Ntau1
  dtid_doodson(9) = '153.655'
  dtidip(9)       =    0.1e-12
  dtidop(9)       =    0.0e-12
! unnamed
  dtid_doodson(10) = '155.445'
  dtidip(10)       =    0.1e-12
  dtidop(10)       =    0.0e-12
! Lk1
  dtid_doodson(11) = '155.455'
  dtidip(11)       =    0.4e-12
  dtidop(11)       =    0.0e-12
! No1
  dtid_doodson(12) = '155.655'
  dtidip(12)       =    1.3e-12
  dtidop(12)       =   -0.1e-12
! unnamed
  dtid_doodson(13) = '155.655'
  dtidip(13)       =    0.3e-12
  dtidop(13)       =    0.0e-12
! chi1
  dtid_doodson(14) = '157.455'
  dtidip(14)       =    0.3e-12
  dtidop(14)       =    0.0e-12
! unnamed
  dtid_doodson(15) = '157.465'
  dtidip(15)       =    0.1e-12
  dtidop(15)       =    0.0e-12
! pi1
  dtid_doodson(16) = '162.556'
  dtidip(16)       =   -1.9e-12
  dtidop(16)       =    0.1e-12
! unnamed
  dtid_doodson(17) = '163.545'
  dtidip(17)       =    0.5e-12
  dtidop(17)       =    0.0e-12
! P1
  dtid_doodson(18) = '163.555'
  dtidip(18)       =  -43.4e-12
  dtidop(18)       =    2.9e-12
! unnamed
  dtid_doodson(19) = '164.554'
  dtidip(19)       =    0.6e-12
  dtidop(19)       =    0.0e-12
! S1
  dtid_doodson(20) = '164.556'
  dtidip(20)       =    1.6e-12
  dtidop(20)       =   -0.1e-12
! unnamed
  dtid_doodson(21) = '165.345'
  dtidip(21)       =    0.1e-12
  dtidop(21)       =    0.0e-12
! unnamed
  dtid_doodson(22) = '165.535'
  dtidip(22)       =    0.1e-12
  dtidop(22)       =    0.0e-12
! unnamed
  dtid_doodson(23) = '165.545'
  dtidip(23)       =   -8.8e-12
  dtidop(23)       =    0.5e-12
! K1
  dtid_doodson(24) = '165.555'
  dtidip(24)       =  470.9e-12
  dtidop(24)       =  -30.2e-12
! unnamed
  dtid_doodson(25) = '165.565'
  dtidip(25)       =   68.1e-12
  dtidop(25)       =   -4.6e-12
! unnamed
  dtid_doodson(26) = '165.575'
  dtidip(26)       =   -1.6e-12
  dtidop(26)       =    0.1e-12
! unnamed
  dtid_doodson(27) = '166.455'
  dtidip(27)       =    0.1e-12
  dtidop(27)       =    0.0e-12
! unnamed
  dtid_doodson(28) = '166.544'
  dtidip(28)       =   -0.1e-12
  dtidop(28)       =    0.0e-12
! psi1
  dtid_doodson(29) = '166.554'
  dtidip(29)       =  -20.6e-12
  dtidop(29)       =   -0.3e-12
! unnamed
  dtid_doodson(30) = '166.556'
  dtidip(30)       =    0.3e-12
  dtidop(30)       =    0.0e-12
! unnamed
  dtid_doodson(31) = '166.564'
  dtidip(31)       =   -0.3e-12
  dtidop(31)       =    0.0e-12
! unnamed
  dtid_doodson(32) = '167.355'
  dtidip(32)       =   -0.2e-12
  dtidop(32)       =    0.0e-12
! unnamed
  dtid_doodson(33) = '167.365'
  dtidip(33)       =   -0.1e-12
  dtidop(33)       =    0.0e-12
! phi1
  dtid_doodson(34) = '167.555'
  dtidip(34)       =   -5.0e-12
  dtidop(34)       =    0.3e-12
! unnamed
  dtid_doodson(35) = '167.565'
  dtidip(35)       =    0.2e-12
  dtidop(35)       =    0.0e-12
! unnamed
  dtid_doodson(36) = '168.554'
  dtidip(36)       =   -0.2e-12
  dtidop(36)       =    0.0e-12
! theta1
  dtid_doodson(37) = '173.655'
  dtidip(37)       =   -0.5e-12
  dtidop(37)       =    0.0e-12
! unnamed
  dtid_doodson(38) = '173.665'
  dtidip(38)       =   -0.1e-12
  dtidop(38)       =    0.0e-12
! unnamed
  dtid_doodson(39) = '175.445'
  dtidip(39)       =    0.1e-12
  dtidop(39)       =    0.0e-12
! J1
  dtid_doodson(40) = '175.455'
  dtidip(40)       =   -2.1e-12
  dtidop(40)       =    0.1e-12
! unnamed
  dtid_doodson(41) = '175.465'
  dtidip(41)       =   -0.4e-12
  dtidop(41)       =    0.0e-12
! So1
  dtid_doodson(42) = '183.555'
  dtidip(42)       =   -0.2e-12
  dtidop(42)       =    0.0e-12
! unnamed
  dtid_doodson(43) = '185.355'
  dtidip(43)       =   -0.1e-12
  dtidop(43)       =    0.0e-12
! Oo1
  dtid_doodson(44) = '185.555'
  dtidip(44)       =   -0.6e-12
  dtidop(44)       =    0.0e-12
! unnamed
  dtid_doodson(45) = '185.565'
  dtidip(45)       =   -0.4e-12
  dtidop(45)       =    0.0e-12
! unnamed
  dtid_doodson(46) = '185.575'
  dtidip(46)       =   -0.1e-12
  dtidop(46)       =    0.0e-12
! v1
  dtid_doodson(47) = '195.455'
  dtidip(47)       =   -0.1e-12
  dtidop(47)       =    0.0e-12
! unnamed
  dtid_doodson(48) = '195.465'
  dtidip(48)       =   -0.1e-12
  dtidop(48)       =    0.0e-12

! degree 2,2 semi-diurnal
! N2
  sdtid_doodson(1) = '245.655'
  sdtidip(1)       =   -0.3e-12
  sdtidop(1)       =    0.0e-12
! M2
  sdtid_doodson(2) = '255.555'    ! PT211013: fixed bug (it said 245.655 but should be 255.555)
  sdtidip(2)       =   -1.2e-12
  sdtidop(2)       =    0.0e-12

  return
  end
















