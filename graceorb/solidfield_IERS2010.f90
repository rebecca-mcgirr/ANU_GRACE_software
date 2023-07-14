!> Subroutine to estimate the effect of the solid earth tide on the gravitational potential fields.
!! The method follows the IERS standards 2010 (Petit et al 2010, IERS Technical Note No. 36 Chap 6)
!! @param [in] rlat, rlong, rx coordinates of the sattelite
!! @param [in] jd (int), t(float) julian date in 2 parts (for effect of the tide)
!! @param [in] rot_i2e rotation matrix
!! @param [out] gravsolid acceleration due the solid earth tide tide component on thye satelite 
!!
!! @todo Later output spherical harmonic coefficients instead of acceleration on the satellite..
!!
!! @Author S. Allgeyer
!! @Date 14 July 2021
!!
!! @note rewritten from scratch.  

subroutine solidfield_iers2010(rlat, rlong, rx, jd, t, rot_i2e, gravsolid)



use coor_mod
use spherhar_mod
use gm_mod
use etides_mod           ! IERS 2010 solid body tides variables 

implicit none

	real(kind=8), intent(inout) :: rlat, rlong, rx   ! earth-fixed lat/lon/radius (is that what rx is?)
	integer*4,    intent(in)    :: jd
	real(kind=8), intent(in)    :: t
	real(kind=8), intent(in)    :: rot_i2e(3,3)      ! inertial to rfixed rotation matrix

! coefficients for the body tide potential perturbations
	real(kind=8) :: cztid(4),cctid(9),cstid(9)       

! intermediate computations for Sun, Moon coords
	real(kind=8) :: C_mlon_m(4),S_mlon_m(4),C_mlon_s(4),S_mlon_s(4)
	real(kind=8) :: clatm,slatm,clats,slats,slonm, slons, clonm, clons
	real(kind=8) :: cmoon2,csun2,cmoon3,csun3

! counters
	integer*4    :: i,i1,idx

! Doodson tide angle and time variables
	real(kind=8) :: theta, tjd

! Legendre polynomials
! PT210918: commented out. It conflicts with the entries in etides_mod
	!real(kind=8) :: Pn_mon(4), Pn_Sun(4)
	!real(kind=8) :: Pnm_mon(9), Pnm_Sun(9)

	integer*4 :: nm, iord, ideg
	real(kind=8) :: twopi, tdtoff, gmst  
	real(kind=8) :: xeqm, xeqs, rpc2, rpc, rc, ertrad, erad3 
	real(kind=8) :: gmor3m, gmor3s, rc2, rpc3, rc3 
	real(kind=8) :: sumrt, sumlatt,sumlont, multl, sumr, sumlat, sumlon
	real(kind=8), dimension(3) :: gravsolid
	real(kind=8), dimension(3) :: grav
	real(kind=8), dimension(6) :: fund_arg
	real(kind=8), dimension(6) :: beta
	double precision :: teph 
	real(kind=8) :: delc20, delc21, dels21, delc22, dels22
	
	double precision, dimension(2) :: plm_moon(0:4,0:4), plm_sun(0:4,0:4) 
	double precision :: radius_earth
	double precision, dimension(2) :: knm_elastic(0:4, 0:4)

	double precision, dimension (2) :: C_coef (0:4, 0:4), S_coef (0:4, 0:4)
	double precision :: cmoon, csun
	double precision	:: lon_moon, lon_sun
	double precision, dimension(2), allocatable :: pnm_dummy(:,:)
	!!sofa calls:
	double precision :: iau_FAL03, iau_FALP03, iau_FAF03, iau_FAD03, iau_FAOM03, iau_GMST82, iau_GST94, iau_GMST06

	integer :: IY, IM, ID,  J_flag
	double precision :: TAI_UTC, FD


	C_coef(:,:) = 0.0
	S_coef(:,:) = 0.0

	radius_earth = 6.371229D06
	tdtoff = 51.184d0
	twopi = 2.d0*pi

	cztid = 0.d0
	cctid = 0.d0
	cstid = 0.d0

 

	knm_elastic(:,:) = 0.0
	! fill the array
	knm_elastic(2,0) = 0.29525d0
	knm_elastic(2,1) = 0.29470d0
	knm_elastic(2,2) = 0.29801d0
	knm_elastic(3,0) = 0.093d0
	knm_elastic(3,1) = 0.093d0
	knm_elastic(3,2) = 0.093d0
	knm_elastic(3,3) = 0.094d0
	!special case for deg4
	knm_elastic(4,0) = -0.00087d0
	knm_elastic(4,1) = -0.00079d0
	knm_elastic(4,2) = -0.00057d0
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
	! print*, slatm, slats
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  sine, cosine of geocentric east longitudes of Moon, Sun in body-fixed coordinates
	xeqm = rpc * dsqrt( 1.0d0 - slatm**2 )
	xeqs = rc * dsqrt( 1.0d0 - slats**2 )
    slonm =   (rot_i2e(2,1)*pccor(1) + rot_i2e(2,2)*pccor(2) +rot_i2e(2,3)*pccor(3) )/ xeqm
    slons = - (rot_i2e(2,1)*ccor(1)  + rot_i2e(2,2)*ccor(2)  +rot_i2e(2,3)*ccor(3)  )/ xeqs
    clonm =   (rot_i2e(1,1)*pccor(1) + rot_i2e(1,2)*pccor(2) +rot_i2e(1,3)*pccor(3) )/ xeqm
    clons = - (rot_i2e(1,1)*ccor(1)  + rot_i2e(1,2)*ccor(2)  +rot_i2e(1,3)*ccor(3)  )/ xeqs

	! lon_moon = asin(slonm)
	! lon_sun  = asin(slons)
	lon_sun = atan2(slons, clons)
	lon_moon = atan2(slonm, clonm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	plm_moon = 0.0
	call legendre_matrix(4,slatm,plm_moon)
	call legendre_matrix(4,slats,plm_sun)

	! degree 2, 3 case (6.6)
	do ideg = 2,3
		cmoon = gm(2)/ gm(1) * (radius_earth/rpc)**(ideg+1)
		csun = gm(3)/ gm(1) * (radius_earth/rc)**(ideg+1)
		do iord = 0,ideg
			C_coef(ideg, iord) = C_coef(ideg, iord) +  knm_elastic(ideg, iord) / (2*ideg +1) * cmoon * plm_moon(ideg, iord) * dcos(iord*lon_moon)
			S_coef(ideg, iord) = S_coef(ideg, iord) +  knm_elastic(ideg, iord) / (2*ideg +1) * cmoon * plm_moon(ideg, iord) * dsin(iord*lon_moon)
			C_coef(ideg, iord) = C_coef(ideg, iord) +  knm_elastic(ideg, iord) / (2*ideg +1) * csun  * plm_sun (ideg, iord) * dcos(iord*lon_sun )
			S_coef(ideg, iord) = S_coef(ideg, iord) +  knm_elastic(ideg, iord) / (2*ideg +1) * csun  * plm_sun (ideg, iord) * dsin(iord*lon_sun )
		enddo
	enddo
	! degree 4 special case (6.7)
	ideg = 4
	cmoon = gm(2)/ gm(1) * (radius_earth/rpc)**(3)
	csun = gm(3)/ gm(1) * (radius_earth/rc)**(3)
	do iord = 0, 2
		C_coef(ideg, iord) = C_coef(ideg, iord) + knm_elastic(ideg, iord) / (5) * cmoon * plm_moon(2, iord) * dcos(iord*lon_moon)
		S_coef(ideg, iord) = S_coef(ideg, iord) + knm_elastic(ideg, iord) / (5) * cmoon * plm_moon(2, iord) * dsin(iord*lon_moon)

		C_coef(ideg, iord) = C_coef(ideg, iord) + knm_elastic(ideg, iord) / (5) * csun * plm_sun(2, iord) * dcos(iord*lon_sun)
		S_coef(ideg, iord) = S_coef(ideg, iord) + knm_elastic(ideg, iord) / (5) * csun * plm_sun(2, iord) * dsin(iord*lon_sun)
	enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! Now start step 2. 
	tjd = jd + (t+tdtoff)/86400.d0
	teph = (tjd - 2451545.0 ) / 36525.0

	fund_arg(1) = iau_FAL03(teph)
	fund_arg(2) = iau_FALP03(teph)
	fund_arg(3) = iau_FAF03(teph)
	fund_arg(4) = iau_FAD03(teph)
	fund_arg(5) = iau_FAOM03(teph)

	CALL iau_JD2CAL ( 0.0d0, jd+t/86400.0, IY, IM, ID, FD, J_flag )
	CALL iau_DAT ( IY, IM, ID, FD, TAI_UTC, J_flag )

	fund_arg(6) = iau_GMST06(0.0d0, tjd, 0.0d0 , tjd-(TAI_UTC+32.184)/8600.)+ pi
	gmst = iau_GMST06(0.0d0, tjd, 0.0d0 , tjd-(TAI_UTC+32.184)/8600.) ! - pi

	beta(1) = gmst + pi - fund_arg(3) - fund_arg(5)
	beta(2) = fund_arg(3) + fund_arg(5)
	beta(3) = beta(2) - fund_arg(4)
	beta(4) = beta(2) - fund_arg(1)
	beta(5) = -1.0* fund_arg(5)
	beta(6) = beta(2) - fund_arg(4) - fund_arg(2)

	do i = 1, 6
		beta(i) = norm_angle(beta(i))
	enddo


! Zonal tidal constituents. There are "num_deg20" of them defined in etides_IERS2010.f90
  do i=1,num_deg20
	theta = 0.0
!print*,"deg20, i,ztid_doodson(i)",i,num_deg20,ztid_doodson(i)
	theta= theta_tide ( ztid_doodson(i) , beta)
	C_coef(2, 0) = C_coef(2, 0) + ztidip(i)*dcos(theta) - ztidop(i)*dsin(theta)


  enddo

! Tesseral tides. There are "num_deg21" of them defined in etides_IERS2010.f90
  do i=1,num_deg21 
	theta = 0.0
	theta = theta_tide ( dtid_doodson(i) , beta)
	C_coef(2, 1) = C_coef(2, 1) + dtidip(i)*dsin(theta) + dtidop(i)*dcos(theta)
        ! SA/PT 211013: changed dtidip to dtidop to fix a bug
	S_coef(2, 1) = S_coef(2, 1) + dtidip(i)*dcos(theta) - dtidop(i)*dsin(theta)

  enddo

! Semi-diurnal tides. There are "num_deg22" of them defined in etides_IERS2010.f90
  do i=1,num_deg22 
	theta = 0.0
	theta = theta_tide ( sdtid_doodson(i) , beta)
	C_coef(2, 2) = C_coef(2, 2) + sdtidip(i)*dcos(theta) 
	S_coef(2, 2) = S_coef(2, 2) - sdtidip(i)*dsin(theta) 
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The following is the same as gravstatic subroutine and modified for the 
! solid earth tide contribution (only deg/ord 20 21 22 and 40 41 42 have non-zero values)
  nm = 4
  sumrt = 0.d0
  sumlatt = 0.d0
  sumlont = 0.d0

  do ideg = 2, 4 
!   compute the sectorial term (order = 0). Here, m=0, therefore C * cos(m*lambda) = C (and the S coefficient is zero)

	multl = (Ae/rx)**ideg

	sumrt = sumrt+Plm(nm) * C_coef(ideg, 0) *dble(ideg+1)*multl   
	sumlatt = sumlatt+dcos(rlat)*PDlm(nm) * C_coef(ideg, 0)*multl
	sumlont = sumlont+0.d0

	nm = nm +1

	do iord = 1,ideg 
	  sumr =  Plm(nm)*multl*dble(ideg+1)*(C_coef(ideg, iord)* &
			  dcos(iord*rlong) + S_coef(ideg, iord) *dsin(iord*rlong) )
	  sumlat = dcos(rlat)*( PDlm(nm)*multl*( C_coef(ideg, iord) * &
			   dcos(iord*rlong)+S_coef(ideg, iord)*dsin(iord*rlong) ))
	  sumlon = (dble(iord)*Plm(nm)*multl*(- C_coef(ideg, iord)* &
			   dsin(iord*rlong)+ S_coef(ideg, iord) *dcos(iord*rlong) ))/dcos(rlat)

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

  contains
  	function theta_tide ( Doodson_number , beta)
		double precision :: theta_tide
		double precision, dimension(6) :: beta
		character, dimension(1) :: Doodson_number(7)
		integer, dimension(6) :: Darg

		integer :: i
		theta_tide = 0.0d0
!print*,'theta_tide doodson_number = xx',Doodson_number,'yy'
		read(Doodson_number(1),'(i1)') Darg(1)
		read(Doodson_number(2),'(i1)') Darg(2)
		read(Doodson_number(3),'(i1)') Darg(3)

		read(Doodson_number(5),'(i1)') Darg(4)
		read(Doodson_number(6),'(i1)') Darg(5)
		read(Doodson_number(7),'(i1)') Darg(6)

		do i=2,6
		  Darg(i) = Darg(i) - 5
		enddo  

		! do i =1, 6
			! theta_tide = theta_tide + Darg(i) * beta(i)
		! enddo
		theta_tide = sum(Darg * beta)
		return 
 	end function

	function norm_angle(angle)

		double precision norm_angle, angle, twopi
  														  
		twopi = 2.d0*datan(1.0d0) * 4.0d0
		norm_angle = dmod(angle,twopi)  
		if( norm_angle.lt.0.d0 ) norm_angle = norm_angle + twopi
  
		return
	end function 
  end
