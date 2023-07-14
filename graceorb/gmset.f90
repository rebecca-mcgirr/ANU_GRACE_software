    subroutine gmset

!   subroutine to set various earth parameters
!   mass data used to calculate gmpln for planets (other than the earth) 
!   were taken from: http://ssd.jpl.nasa.gov/?planet_phys_par

    use gm_mod

    implicit none

    Gconst = 6.674e-11 ! Gravitational constant
! PT150603: change from 4415 to the IERS 2010 standard value of 4418
    gm(1)= 398600.4418D9 !GM constant for the earth
    gm(2) = gm(1)/81.3005883D0 !GM constant for the moon
    gm(3) = (gm(1)+gm(2))*328900.55523496d0 ! GM constant for the Sun
    gmpln(1) = 0.330104d0*1.d24*Gconst ! Mercury 
    gmpln(2) = 4.86732d0*1.d24*Gconst ! Venus
    gmpln(4) = 0.641693d0*1.d24*Gconst ! Mars
    gmpln(5) = 1898.13d0*1.d24*Gconst ! Jupiter
    gmpln(6) = 568.319d0*1.d24*Gconst ! Saturn
    gmpln(7) = 86.8103d0*1.d24*Gconst ! Uranus
    gmpln(8) = 102.410d0*1.d24*Gconst ! Neptune
! PT130829: what happens if I increase this by 0.1m?
!           or decrease by 10 m
!   Ae = 6.3781363d6 ! radius of the earth
    Ae = 6.37813646d6 ! radius of the earth  PT140801: increased this from 6.37813646d6 to 6.37813666d6 for a test
    pi = 4.0d0*datan(1.d0) ! pi
    Me = 5.9736e24 ! mass of the earth
    rhow = 1025.d0 ! density of water
    Om = 7.292115d-5 ! Earth's angular velocity

    return
    end
