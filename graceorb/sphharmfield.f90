   subroutine sphharmfield(rlat, rlong, rx, jd, tin,  g_cart )  

! Emma-Kate Potter, 23 April, 2010
! This program outputs the acceleration due to gravity in a cartesian coordinate 
! system in earth fixed frame at the input earth-fixed coordinates.
! The acceleration uses a spherical harmonic model of degree, ndeg, to determine 
! the gravitational potential field. The gradient of the potential field is then 
! calculated in a spherical coordinate system before being transformed to cartesian.

! The contributions to the spherical harmonic coefficients are:
! 1) the mean/static field model of the earth (staticfield)
! 2) the ocean tide model (oceanfield)
! 3) the atmospheric tide model (atmfield)
! 4) the non-tidal atmosphere and ocean (dealias)

! This code is modified from a spherical harmonic evaluation from
! P Tregonning, 18 June, 2007

! APP130226: This subroutine is called from gravcalc which is often executed
!            for the same time step (but different satellite positions).
!            To reduce superfluous operations a check has been introduced
!            to see if the current time step is the same as the previous
!            time step for which we generated a gravity field. If so, the
!            previously generated field is re-used, otherwise generate a new
!            field using standard formulation

! use modules for common variables
  use spherhar_mod  ! common module to store legendre functions
                    ! common to legcalc subroutine
  use coeff_mod     ! common module to store spher harm coefficients, common also 
                    ! to adam_moulton subroutine
  use gm_mod        ! stores earth and other constants
  use inmod_mod     ! stores the input model flags and details
  use dealias_mod   ! stores the input model flags and details
  use gauleg_mod
  use omtides_mod   ! brings in information for the long-period ocean tides

  implicit none

  real(kind=8),intent(inout) :: rlat, rlong, rx    !  efixed lat/lon/radius (spherical coords?)
  integer*4   ,intent(in)    :: jd                 !  julian day
  real(kind=8),intent(in)    :: tin                !  seconds of day
  real(kind=8),intent(out)   :: g_cart(3)          !  cartesian accelerations (in e-fixed frame)

  integer  lucoeff, ioerr, ierr, counter
  integer :: ideg,iord,nm,nord
  integer :: i,j,k,l,m,n, icount
  real(kind=8) :: x,cval,sval
  real(kind=8) :: temp1, temp2, error_check, t
  real(kind=8) :: multl, sumrt, sumlatt, sumlont, sumr, sumlat, sumlon 
  real(kind=8), dimension(3) :: grav, pos 
  integer:: tarray(3)
  integer*4 :: nvals    ! dimensioning for coefficients in the dealiasing file (see mod_subs/dealias_mod.f90)

! tidal acceleration computations
    real(kind=8) :: g_tide(3)

! PT161007: individual computations of each spherical harmonic field to get the individual effects
    real(kind=8) :: static_acc(3),delias_acc(3),opole_acc(3),atm_tide_acc(3),dealias_acc(3),atide_acc(3),total_acc(3)
    real(kind=8) :: factor

! Debug flag
    logical :: debug

! PT161007: save the g_tide computation for when we don't move to a new epoch
    save g_tide

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check if current time step is the same as the previous

  if ( ( jd .eq. jd_old ) .and. ( abs(tin - t_old) .le. 0.01 ) ) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! APP130226:
! if current time step is the same as the previous timestep for which we calculated
! a spherical harmonic field, then re-use the previously generated field
    do ideg = 0 , ndeg
      do iord = 0 , ideg
        coefC(ideg,iord) = coefC_old(ideg,iord)
        coefS(ideg,iord) = coefS_old(ideg,iord)
      enddo
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  else

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! APP130226: Otherwise we generate the coefficients of the spherical harmonic field
! initialise the set of spherical harmonic coefficients
    coefC = 0.d0
    coefS = 0.d0 
    dCdealias = 0.d0
    dSdealias = 0.d0
    dCotide = 0.d0
    dSotide = 0.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! STATICFIELD
! Already have spherical harmonic coefficients of mean field in the correct format
! which were read in by coeff_read (called in graceorb)
    if (gt_statfield(1).eq.'Y') then
! PT130923: the GOCE static gravity field, go_cons_gcf_2_dir_r4, starts at degree 2
      do ideg = 2, max_statgrav_degree
        do iord = 0, ideg
          coefC(ideg, iord) = meancoefC(ideg, iord)
          coefS(ideg, iord) = meancoefS(ideg, iord)
        enddo
      enddo
    endif  
! DEBUG: evaluate the spherical harmonics at this location
! calculate the associated Legendre functions
  call legcalc_norm(sin(rlat))
factor = 1.0
debug = .false.
if(debug) then
  print*,'static_acc computation for time',tin
  call compute_sph_acc(rlat,rlong,rx,coefC,coefS,max_statgrav_degree,Plm,PDlm,"stokes",static_acc)
  print*,'tin,static_acc: ',tin,static_acc
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! OCEANFIELD
! Call for oceanfield subroutine to call for the calculation of the tidal height grid
! for current epoch, convert this into C&S coefficients for tidal height and scale
! to be added to mean field C&S coefficents for calculating potential
    if (gt_oceantide(1).eq.'Y') then
! PT160909: no longer allow the code to use spotl routines to calculate the grid on the fly
!      call oceanfield(jd, tin, rlat,rlong)
      call status_update('FATAL','GRACEORB','sphharmfield',' ',"Code will no longer convert tide heights to spherical harmonics",0)
    endif

! but we will allow spherical harmonic representations to be read in and used 
    if (gt_oceantide(1).eq.'S') then
      call oceanfield_sph(jd,tin,rlat,rlong, rx,g_tide)
!
!! PT130923: we should use at least from degree 1 onwards (code was from deg 2 to deg 100
      do ideg = 1 , 200  ! ocean tide model goes up to 256, but use 200 here (????)
        do iord = 0 , ideg
! PT130926: turn off the addition of the spherical harmonic representation of the ocean tide. Truncation errors
!           can be as large as 0.7 m when a 3.5 m tide height occurs at the ocean/coast interface
! PT161006: don't use the coefficients anymore - use the tidal accelerations computed in oceanfield_sph
!          coefC(ideg,iord) = coefC(ideg,iord) + dCotide(ideg,iord)
!          coefS(ideg,iord) = coefS(ideg,iord) + dSotide(ideg,iord)
        enddo
      enddo
    else
      g_tide = 0.d0
    endif  
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Seb's code: add for ocean tides
    call OCT_CS(octides_data, jd, tin)
	do ideg = 0 , 200
		do iord = 0, ideg
			coefC(ideg,iord) = coefC(ideg,iord) + dCotide(ideg,iord)
			coefS(ideg,iord) = coefS(ideg,iord) + dSotide(ideg,iord)
		enddo
	enddo

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ATMFIELD
! Call for atmfield subroutine to call for the calculation of a pressure grid
! for current epoch (by interpolation?), convert this into C&S coefficients for 
! pressure and then scale to be added to mean field C&S coefficents for calculating potential
    if (gt_atmtide(1).eq.'Y') then
!      call status_update('FATAL','GRACEORB','sphharmfield',' ',"Atmospheric model code not yet implemented in integrator",0)
      call atm_tide_sph(jd, tin, rlat,rlong, rx, atm_tide_acc)
      do ideg = 2 , 100 ! atm tide model goes up to 180, but use 100 here (????) 
        do iord = 0 , ideg
          coefC(ideg,iord) = coefC(ideg,iord) + dCatide(ideg,iord)
          coefS(ideg,iord) = coefS(ideg,iord) + dSatide(ideg,iord)
        enddo
      enddo
    endif
if(debug) then
  call compute_sph_acc(rlat,rlong,rx,dCatide,dSatide,maxsize,Plm,PDlm,"stokes",atide_acc)
  print*,'tin,atide_acc: ',tin,atide_acc
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DEALIAS
! Call for dealias subroutine which take the C&S coefficients for non-tidal ocean and atm
! variations at 6hourly intervals and interpolates to determine the contribution at 
! current epoch
    if (gt_dealias(1).eq.'Y') then
! PT150528: the call to dealias was commented out (back in 2013?) because all the coefficients are now read in inside graceorb.f90
!      call dealias(jd, tin)
!     calculate seconds past epoch (spe) Jan 1 2000 from jd and tin
      t = (dble(jd)+tin/86400.0d0-2451545.0d0)*86400.0d0
      nvals = (aod1b_degree+1)*(aod1b_degree+2)/2
      call lininterp_r8(max_num_epoch, nvals, epoch_time, CT, t, num_dealias_epochs, nvals, dCdealias, error_check)
      call lininterp_r8(max_num_epoch, nvals, epoch_time, ST, t, num_dealias_epochs, nvals, dSdealias, error_check)

!      write(6, *) t, " Degree 2 Dealias: ", dCdealias(1), dSdealias(1)
      icount = 0
! PT130923: we should use at least from degree 1 onwards (code was from deg 2 to deg 100)
! PT131108: no, in graceorb the dealiasing is read in from deg2 onwards, so have to be consistent in both places!
!           either change both or leave both alone !!!
      do ideg = 2 , 100  ! dealiasing files only go to 100
        do iord = 0 , ideg
          icount = icount + 1
          coefC(ideg,iord) = coefC(ideg,iord) + dCdealias(icount)
          coefS(ideg,iord) = coefS(ideg,iord) + dSdealias(icount)
        enddo
      enddo
    endif
!factor = 1.0
if(debug) then
  call compute_sph_acc(rlat,rlong,rx,dCdealias,dSdealias,maxsize,Plm,PDlm,"stokes",dealias_acc)
  print*,'tin,dealias_acc: ',tin,dealias_acc,rlat,rlong
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! POLE TIDES (SOLID EARTH AND OCEAN)
! Call for poltid to calculate the contribution of the solid earth pole tide to the C21 and S21 coefficients
! PT150528: changed this from gt_dealias(1) to gt_poletide(1)
    if (gt_poletide(1).eq.'Y') then
      poltype = "both" ! replace with "sold" for solid earth or "ocen" for ocean pole tide
      call poltid(jd,tin,poltype)
if(debug)print*,'dCsolpol21,dSsolpol21',dCsolpol21,dSsolpol21
      coefC(2,1) = coefC(2,1) + dCsolpol21 
      coefS(2,1) = coefS(2,1) + dSsolpol21 
! PT160121: why aren't we adding the ocean pole tide ? Add it on here
      do ideg = 2,ndeg
        do iord = 0,ideg
          coefC(ideg,iord) = coefC(ideg,iord) + dCocepol(ideg,iord) 
          coefS(ideg,iord) = coefS(ideg,iord) + dSocepol(ideg,iord)
        enddo
      enddo
    endif
!factor = 1.0
if(debug) then
  call compute_sph_acc(rlat,rlong,rx,dCocepol,dCocepol,maxsize,Plm,PDlm,"stokes",opole_acc)
  print*,'tin,opole_acc: ',tin,opole_acc
  print*,'solpol21 = ',dCsolpol21,dSsolpol21
  factor = 1.0
  call compute_sph_acc(rlat,rlong,rx,coefC,coefS,maxsize,Plm,PDlm,"stokes",opole_acc)
  print*,'tin,ALL_acc: ',tin,opole_acc
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! APP130226: Store spherical harmonic coefficients in case they need to be re-used
    do ideg = 0 , ndeg
      do iord = 0 , ideg
        coefC_old(ideg,iord) = coefC(ideg,iord)
        coefS_old(ideg,iord) = coefS(ideg,iord)
      enddo
    enddo
! Store current time to check in subsequent call
    t_old = tin
    jd_old = jd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DEBUG: print out the total gravitational potential
if(debug) then
  call compute_sph_acc(rlat,rlong,rx,coefC,coefS,maxsize,Plm,PDlm,"stokes",total_acc)
  print*,'tin,total_acc: ',tin,total_acc
endif

  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE the components of gravitational acceleration in spherical coordinates due
! to all spherical harmonic models (above)

! APP130226: 
! The form of the gravitational potential is:
!    U(lambda, phi, r) = sum_n sum_m U_nm
! where:
!    Unm = P_nm*(cos(phi)) *[ C_nm*cos(m*lambda) + S_nm*sin(m*lambda) ] * (Ae/r)^(n+1)
! Gravitational acceleration is found by taking the gradient of the potential:
!    grav = grad U
! gravity therefore has the following components:
!   radial:          dU_nm/dr
!                       = -(n + 1)/Ae * P_nm*(cos(phi)) *[ C_nm*cos(m*lambda) + S_nm*sin(m*lambda) ] * (Ae/r)^(n+2)
!   co-latitudinal:  1/r * dU_nm/dphi
!                       = -sin(phi)/Ae * P'_nm*(cos(phi)) *[ C_nm*cos(m*lambda) + S_nm*sin(m*lambda) ] * (Ae/r)^(n+2)
!   longitudinal:    1/(r * sin(phi)) * dU_nm/dlambda
!                       = m/(Ae * sin(phi)) * P_nm*(cos(phi)) *[ S_nm*cos(m*lambda) - C_nm*sin(m*lambda) ] * (Ae/r)^(n+2)
!   it should be noted that rlat = pi/2 - phi, so changing from co-latitudinal to latitudinal changes the sign
!   of the corresponding component

! set error flag to zero
  ioerr = 0

! starting at degree/order of 2/0 means that the first legendre 
! function required is element 4
  nm = 4
! PT130923: but for deg 1 it is 2 (and the deg1 terms for the tide are important
!  nm = 2
  sumrt = 0.d0
  sumlatt = 0.d0
  sumlont = 0.d0 

! PT130920: evaluation of the spherical harmonic model takes ~7 nanoseconds
  do ideg = 2,  ndeg 
!  do ideg = 2,  11    ! To match arc integration
!  do ideg = 2,  2     ! To match Matlab/arc integration
!   do ideg = 2,  0     ! To match central force        

!   compute the sectorial term (order = 0). Here, m=0, therefore the
!   C * cos(m*lambda) = C (and the S coefficient is zero)

    multl = (Ae/rx)**ideg

    sumrt = sumrt+Plm(nm)*coefC(ideg,0)*dble(ideg+1)*multl   ! * cos(m*lambda) 
    sumlatt = sumlatt+dcos(rlat)*PDlm(nm)*coefC(ideg,0)*multl
    sumlont = sumlont+0.d0

    nm = nm +1

   do iord = 1, ideg    
      sumr =  Plm(nm)*multl*dble(ideg+1)*(coefC(ideg,iord)* &
              dcos(iord*rlong)+coefS(ideg,iord)*dsin(iord*rlong) )
      sumlat = dcos(rlat)*( PDlm(nm)*multl*(coefC(ideg,iord)* &
               dcos(iord*rlong)+coefS(ideg,iord)*dsin(iord*rlong) ))
      sumlon = (dble(iord)*Plm(nm)*multl*(-coefC(ideg,iord)* &
               dsin(iord*rlong)+coefS(ideg,iord)*dcos(iord*rlong) ))/dcos(rlat)

      sumrt = sumrt + sumr
      sumlatt = sumlatt + sumlat
      sumlont = sumlont + sumlon

      nm = nm + 1
    enddo
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set longitude range to -180/180
  if (rlong.gt.pi)rlong = rlong-2.d0*pi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SCALE to give gravitational potential at earth's surface by multiplying 
! by GM/R values of gravitational constants and earth radius 
! grav() array gives r, theta, phi elements of gravity vector in spherical coordinates         
  grav(1) = -gm(1)/((rx)**2)*(sumrt+1.d0) 
  grav(2) = gm(1)/((rx)**2)*(sumlatt)
  grav(3) = gm(1)/((rx)**2)*(sumlont)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!if(debug)print*,'tin, grav',tin,grav
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! convert calculated gravity into cartesian coordinates:
  call cart_convert(rlat,rlong,grav,g_cart)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(debug)print*,'**** **** tin, grav,g_cart,g_tide',tin, grav,g_cart,g_tide

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! add on the ocean tide acceleration from spherical harmonics if necessary
  if(gt_oceantide(1).eq.'S')g_cart = g_cart + g_tide
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  return
  end

