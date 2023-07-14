  subroutine solidfield(rlat, rlong, rx, jd, t, rot_i2e, gravsolid)

! Emma-Kate Potter, 11 May, 2010
! This program calculates the perturbation to the gravity field due to 
! solid earth tides.

! Calculations are lifted from GAMIT arc subroutines sbfn and sbfn1 and modified

! Determine quantities for solid-Earth tides
! Effect for C20, C21, S21, C22, and S22, according to IERS Standards
! (1992), Chapt. 7, pp 52-55.  Use degree & order Love numbers and
! frequency-dependent corrections for n=2,m=1 constituents from Richard
! Eanes (06 April 1994) instead of IERS Standards Table 7.1.
! Initialize parameters and variables
! (cztid, cctid, and cstid arrays zeroed out in init.f)
! nominal Love numbers for degree 2, order 0, 1, 2
! (values from Richard Eanes, 06 April 1994)

! values required for  
! and fjd (which is the actual time calculated from jd and t)
    use coor_mod
!    use rotation_mod
    use spherhar_mod
    use gm_mod

    implicit none

    real(kind=8), intent(inout) :: rlat, rlong, rx   ! earth-fixed lat/lon/radius (is that what rx is?)
    integer*4,    intent(in)    :: jd
    real(kind=8), intent(in)    :: t
    real(kind=8), intent(in)    :: rot_i2e(3,3)      ! inertial to rfixed rotation matrix

    integer*4 :: nm, iord, ideg
    real(kind=8) :: klove20, klove21, klove22
    real(kind=8) :: twopi, tdtoff, gmst  
    real(kind=8) :: slatm, slats, xeqm, xeqs, rpc2, rpc, rc, ertrad, erad3 
    real(kind=8) :: slonm, slons, clonm, clons 
    real(kind=8) :: gmor3m, gmor3s, rc2, rpc3, rc3 
    real(kind=8) :: coef2, coef20, coef21, coef22
    real(kind=8) :: sumrt, sumlatt,sumlont, multl, sumr, sumlat, sumlon
    real(kind=8), dimension(3) :: gravsolid
    real(kind=8), dimension(3) :: grav
    real(kind=8), dimension(2) :: cztid, cctid, cstid 
    real(kind=8), dimension(6) :: fund_arg
    real(kind=8), dimension(6) :: delc21, dels21

    tdtoff = 51.184d0
    twopi = 2.d0*pi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! setting love number values (see GAMIT code)
    klove20 = 0.299d0
    klove21 = 0.300d0
    klove22 = 0.302d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   calculating the earth-sun & earth-moon distances 
!   these coordinates are currently returned from the GAMIT PEP tables
!   need to change to use the same (JPL) coordinates as for the planetfield.f90
    rpc2 = pccor(1)**2+pccor(2)**2+pccor(3)**2
    rpc   =   dsqrt(rpc2)
    rpc3  =   rpc2*rpc
    rc2=ccor(1)**2+ccor(2)**2+ccor(3)**2
    rc =dsqrt(rc2)
    rc3=rc2*rc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  sine of geocentric latitudes of Moon, Sun in body-fixed coordinates
!  (EKP - from inertial to bf)
    slatm = (rot_i2e(3,1)*pccor(1)+rot_i2e(3,2)*pccor(2)+rot_i2e(3,3)*pccor(3) )/ rpc
    slats = -(rot_i2e(3,1)*ccor(1)+rot_i2e(3,2)*ccor(2)+rot_i2e(3,3)*ccor(3) )/ rc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  sine, cosine of geocentric east longitudes of Moon, Sun in body-fixed coordinates
    xeqm = rpc * dsqrt( 1.0d0 - slatm**2 )
    xeqs = rc * dsqrt( 1.0d0 - slats**2 )
    slonm =   (rot_i2e(2,1)*pccor(1) + rot_i2e(2,2)*pccor(2) +rot_i2e(2,3)*pccor(3) )/ xeqm
    slons = - (rot_i2e(2,1)*ccor(1)  + rot_i2e(2,2)*ccor(2)  +rot_i2e(2,3)*ccor(3)  )/ xeqs
    clonm =   (rot_i2e(1,1)*pccor(1) + rot_i2e(1,2)*pccor(2) +rot_i2e(1,3)*pccor(3) )/ xeqm
    clons = - (rot_i2e(1,1)*ccor(1)  + rot_i2e(1,2)*ccor(2)  +rot_i2e(1,3)*ccor(3)  )/ xeqs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  intermediate products
    erad3  = Ae * Ae * Ae
    gmor3m = gm(2) / rpc3
    gmor3s = gm(3) / rc3

    coef2  = (erad3) / (gm(1) * dsqrt(5.d0))
    coef20 = coef2 * klove20
    coef21 = coef2 * klove21 / dsqrt( 3.d0)
    coef22 = coef2 * klove22 / dsqrt(12.d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  changes to normalized degree 2 coefficients -- "Step 1"
!  (see eqns. 1a-1c, p. 52, IERS Standards '92)
    cztid(1) = coef20 * ( &
                 (gmor3m * (3.d0*slatm*slatm - 1.d0) / 2.d0) &
               + (gmor3s * (3.d0*slats*slats - 1.d0) / 2.d0) )
    cctid(1) = coef21 * ( &
                (gmor3m * 3.d0*slatm*dsqrt(1.d0 - slatm*slatm) * clonm) &
              + (gmor3s * 3.d0*slats*dsqrt(1.d0 - slats*slats) * clons) )
    cstid(1) = coef21 * ( &
                (gmor3m * 3.d0*slatm*dsqrt(1.d0 - slatm*slatm) * slonm) &
              + (gmor3s * 3.d0*slats*dsqrt(1.d0 - slats*slats) * slons) )
    cctid(2) = coef22 * ( &
                (gmor3m * 3.d0*(1.d0 - slatm*slatm) * (clonm**2 - slonm**2)) &
              + (gmor3s * 3.d0*(1.d0 - slats*slats) * (clons**2 - slons**2)) )
    cstid(2) = coef22 * ( &
               (gmor3m * 3.d0*(1.d0 - slatm*slatm) * (2.d0*clonm*slonm)) &
             + (gmor3s * 3.d0*(1.d0 - slats*slats) * (2.d0*clons*slons)) )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  changes to normalized (2,1) coefficients -- "Step 2",
!  due to frequency-dependent Love numbers (see eqn. 2, p. 53,
!  IERS Standards '92); coefficients from Richard Eanes, including
!  effect of revised FCN frequency (06 April 95).  Corrections
!  for (2,2) coefficients not needed due to the use of order-
!  dependent Love numbers above.

!  first, get Brown's astronomical angular arguments
!  plus GMST+pi
!  **Note: CALC-8 and MODEST use GAST instead of GMST here!!!!!

!  tide_angles subroutine in GAMIT takes the conventional JD 
 ! call tide_angles( fjd-0.5d0+tdtoff/86400.d0, fund_arg ) !- is the call in ARC's
!                                                            sbfn when the fjd is
!                                                            PEP time
   call tide_angles_gamit( jd+(t+tdtoff)/86400.d0, fund_arg )
!   call fund_arguments_iers1992( jd+(t+tdtoff)/86400.d0, fund_arg )

!  this routine from CALC is GAMIT/lib as part of sd_comp.f
    gmst = fund_arg(6) - pi 

!                                    l l' F D Om GMST+pi
!   (1) psi1 (Doodson = 166.554; Brown = 0 1 0 0 0 1)
!
    delc21(1) = 20.7496d-12 * &
                dsin( gmst + fund_arg(2) )
    dels21(1) = 20.7496d-12 * &
                dcos( gmst + fund_arg(2) )

!               (2) unnamed (Doodson = 165.565; Brown = 0 0 0 0 -1 1)

    delc21(2) = -68.4248d-12 * &
                 dsin( gmst - fund_arg(5) )
    dels21(2) = -68.4248d-12 * &
                 dcos( gmst - fund_arg(5) )

!               (3) K1 (Doodson = 165.555; Brown = 0 0 0 0 0 1)

    delc21(3) = -473.2131d-12 * &
                  dsin( gmst )
    dels21(3) = -473.2131d-12 * &
                 dcos( gmst )

!               (4) unnamed (Doodson = 165.545; Brown = 0 0 0 0 1 1)

    delc21(4) = 8.8203d-12 * &
                 dsin( gmst + fund_arg(5) )
    dels21(4) = 8.8203d-12 * &
                 dcos( gmst + fund_arg(5) )

!               (5) P1 (Doodson = 163.555; Brown = 0 0 -2 2 -2 1)

    delc21(5) = 48.7141d-12 * &
                 dsin( gmst - 2.d0*fund_arg(3) &
                + 2.d0*fund_arg(4) - 2.d0*fund_arg(5) )
    dels21(5) = 48.7141d-12 * &
                 dcos( gmst - 2.d0*fund_arg(3) &
                + 2.d0*fund_arg(4) - 2.d0*fund_arg(5) )

!               (6) O1 (Doodson = 145.555; Brown = 0 0 -2 0 -2 1)

    delc21(6) = 16.4005d-12 * &
                 dsin( gmst - 2.d0*fund_arg(3) - 2.d0*fund_arg(5) )
    dels21(6) = 16.4005d-12 * &
                 dcos( gmst - 2.d0*fund_arg(3) - 2.d0*fund_arg(5) )
!
    cctid(1) = cctid(1) + delc21(1) + delc21(2) + delc21(3) &
              + delc21(4) + delc21(5) + delc21(6)
    cstid(1) = cstid(1) + dels21(1) + dels21(2) + dels21(3) &
              + dels21(4) + dels21(5) + dels21(6)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The following is the same as gravstatic subroutine and modified for the 
! solid earth tide contribution (only deg/ord 20 21 22)
  nm = 4
  sumrt = 0.d0
  sumlatt = 0.d0
  sumlont = 0.d0

  do ideg = 2, 2 
!   compute the sectorial term (order = 0). Here, m=0, therefore the
!   C * cos(m*lambda) = C (and the S coefficient is zero)

    multl = (Ae/rx)**ideg

    sumrt = sumrt+Plm(nm)*cztid(1)*dble(ideg+1)*multl   ! * cos(m*lambda)
    sumlatt = sumlatt+dcos(rlat)*PDlm(nm)*cztid(1)*multl
    sumlont = sumlont+0.d0

    nm = nm +1

    do iord = 1,2 
      sumr =  Plm(nm)*multl*dble(ideg+1)*(cctid(iord)* &
              dcos(iord*rlong)+cstid(iord)*dsin(iord*rlong) )
      sumlat = dcos(rlat)*( PDlm(nm)*multl*(cctid(iord)* &
               dcos(iord*rlong)+cstid(iord)*dsin(iord*rlong) ))
      sumlon = (dble(iord)*Plm(nm)*multl*(-cctid(iord)* &
               dsin(iord*rlong)+cstid(iord)*dcos(iord*rlong) ))/dcos(rlat)

      sumrt = sumrt + sumr
      sumlatt = sumlatt + sumlat
      sumlont = sumlont + sumlon

      nm = nm + 1

    enddo
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set long range to -180/180
  if (rlong.gt.pi)rlong = rlong-2.d0*pi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! scale to give gravitational potential at earth's surface by multiplying
! by GM/R values of gravitational constants and earth radius
! grav() array gives r, theta, phi elements of gravity vector in spherical coordinates
  grav(1) = -gm(1)/((rx)**2)*(sumrt)
  grav(2) = gm(1)/((rx)**2)*(sumlatt)
  grav(3) = gm(1)/((rx)**2)*(sumlont)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! convert calculated gravity into cartesian coordinates:
  call cart_convert(rlat,rlong,grav,gravsolid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  return 
  end
