   subroutine poltid(jd, tin,  poltype)  

! Emma-Kate Potter, 17, Nov, 2010
! This subroutine calculates the change in C21 and S21 due to the effect of the solid earth pole tide
! It uses the pole position file pole. also used by GAMIT
! Based on IERS 2003, section 6.2 and 7.1.4
!
! MOD: PT140819 updated to IERS2010 (and fixed bug in the sign of the m2 term in the C21 correction)
!      PT190722 updated to use the Wahr et al (2016) formulation (eq 22)

    use gm_mod
    use lovenum_mod
    use coeff_mod
    use usno_mod      ! PT161005: Xp, Yp now passed in via this mod file

    implicit none

    integer*4   ,intent(in)  :: jd                 !  julian day
    real(kind=8),intent(in)  :: tin                !  seconds of day
    character*4 ,intent(in ) :: poltype            ! "sold" for solid earth, "ocen" for ocean pole tide, "both" for both

    real*8 :: yr
    real*8 :: xpolav, ypolav, xpole, ypole, xpdot, ypdot
    real*8 :: fract, m1, m2
    real*8 :: jdr, sec 
    integer*4 :: ipole
    integer*4 :: yro, doy 
    integer*4, dimension(5) :: date
    character*4 :: bothstrg, ocenstrg, soldstrg
    real*8 :: Rn
    integer :: n, m, j
    integer, dimension(3) :: tarray
! PT131030: add a local jd_for_polredek variable
    integer*4 jd_for_polredek
! PT140318: added time variable for IERS2010 computations
    real(kind=8) :: dt
    real(kind=8) :: dummy_int

    integer*4  :: i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   calculate julian date from jd and seconds of day
    jdr = dble(jd)+tin/86400.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   define the estimates of the mean pole from IERS 2003, section 7.1.4
    xpolav = 0.054+(jdr-2451544.5d0)/365.25d0*0.00083 
    ypolav = 0.357+(jdr-2451544.5d0)/365.25d0*0.00395
! PT140318: we need to update this to IERS2010 standards, which has a cubic pre-2010 and a linear
!           model post-2010 (p115 of IERS Conventions 2010)
! PT190724: the 2018 update of the IERS pole tide model is just a linear model all the time
    dt = (jdr - 2451544.5d0)/365.25d0
    xpolav =  55.0d-3 + 1.677d-3*dt
    ypolav = 320.5d-3 + 3.460d-3*dt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   call polredek (based on GAMIT subroutine polred) to read in the polar motion values and return 
!   the x and y values at time jd, tin
    fract=tin/86400.d0
! PT131030: I think we need to increment JD if fract becomes 1.0 (ie tin is 86400)
    if(tin == 86400.d0)then
      fract = 0.d0
      jd_for_polredek = jd+1
    else
      jd_for_polredek = jd
    endif
!   interpolate the Xp, Yp to the required epoch. Using iau_interp will also add in the sub-daily signals
!   NOTES: 
!     1. Xp(:,2) are the mjdates to which the Xp, Yp values refer
!     2. We don't care about UT1-UTC here (and haven't passed the information through) so just pass Xp again where ut1utc would be
    call iau_INTERP (Xp(:,2),Xp(:,1),Yp(:,1),Xp(:,1),n_usnovals,jdr-2400000.5d0,xpole,ypole,dummy_int)

!   PRINT *, 'DEBBUGpoltid', jdr-2400000.5d0,xpole,ypole 
!    ipole = 32
!    close(ipole)
!    open(unit=ipole,file='pole.',status='old')
!    call polredek(ipole,jd_for_polredek,fract,xpole,ypole,xpdot,ypdot)
    m1 = xpole-xpolav
    m2 = -1.d0*(ypole-ypolav)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   calculate contribution to C21 and S21 if including both ocean and solid or solid earth only
    if (poltype.eq.'both'.or.poltype.eq.'sold') then
! PT140319: fixed bug in the sign of the m2 term for the Csolpol21 computation (should be +ive. P94 of IERS2010 standards)
      dCsolpol21 = -1.333d-9*(m1+0.0115*m2)
      dSsolpol21 = -1.333d-9*(m2-0.0115*m1)
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (poltype.eq."both" .or. poltype.eq."ocen") then
      m1=m1/3600.d0*pi/180.d0 !convert seconds to radians
      m2=m2/3600.d0*pi/180.d0 !convert seconds to radians

!     doing the calculation on the fly takes about 0.025 seconds (no need for 10 min interpolation)
! PT160121: this equation doesn't match the one in IERS2010 p94 !
!      Rn = Om*Om*(Ae**6)*4.d0*pi*rhow/Gconst/Me/Me
      Rn = Om*Om*(Ae**4)/GM(1) *4.0*pi*Gconst*rhow/9.7803278    ! from Eq 6.23b of IERS2010 p94. 
   
!     calculate the contribution to the C & S coefficients 
      do n=0,100
        do m= 0,n
! PT160121: why is it "+m1*0.0036" for the AI and BI terms? The IERS2010 says "-m1*0.0036" ..... (IERS 2010 p94)
        dCocepol(n,m)= Rn*(1.d0+klove(n))/(2.d0*dble(n)+1.d0)*(AR(n,m)*(m1*0.6870d0+m2*0.0036d0)+AI(n,m)*(m2*0.6870d0-m1*0.0036d0))
        dSocepol(n,m)= Rn*(1.d0+klove(n))/(2.d0*dble(n)+1.d0)*(BR(n,m)*(m1*0.6870d0+m2*0.0036d0)+BI(n,m)*(m2*0.6870d0-m1*0.0036d0))
        enddo
      enddo
      
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    return
    end

