   subroutine planetfield(jdGPS, tGPS, sbcor, C20, sunpos, & 
                          plnflag, gravsm)  
!   subroutine planetfield(jd, t, jdGPS, tGPS, sbcor, C20, sunpos, & 
!                          plnflag, gravsm)  

! Emma-Kate Potter, 18 May, 2011
! This program calculates the point mass contribution to gravity for a satellite due to
! sun and the moon and planets.
!   1) the direct point mass attraction of the planets on the satellite
!   2) the indirect effect of the point mass attraction of the planets on the earth
!   3) the indirect J2 effect due to a point mass sun and moon on the oblateness of the earth

! NOTE THAT THE TIME INPUT IS PEP TIME

! GAMIT subroutines solred and lunred are called to calculate the position of the sun
! and moon respectively.

    use coor_mod     ! store pccor and ccor
    use tabsm_mod    ! store lunar and solar ephemeris data
    use gm_mod       ! store earth and other constants

    implicit none

    integer*4 :: i, ispeed, jd,j, jdGPS
    real(kind=8) ::  time, t, s, tdtoff, tGPS
    real(kind=8), dimension(2) ::  ephtime
    real(kind=8), dimension(3) :: sbcor,  ccor3, pccor3, bcor, pbcor
    real(kind=8), dimension(3) :: gravsm, gravm, gravs
    real(kind=8), dimension(6) :: sunpos
    real(kind=8) :: rc, rc2, rc3 ! earth-sun dist, dist^2 and dist^3 
    real(kind=8) :: rpc, rpc2, rpc3 ! moon-earth dist, dist^2 and dist^3 
    real(kind=8) :: rb, rb2, rb3 ! satellite-sun dist, dist^2 and dist^3 
    real(kind=8) :: rpb, rpb2, rpb3 ! moon-satellite dist, dist^2 and dist^3 
    real(kind=8) :: sums, sump 
    character*5 :: arc_frame
    character*1 :: plnflag 
    real(kind=8), dimension(1:6, 1:11) :: plncor
    real(kind=8), dimension(1:3, 1:11) :: plnbcor, plncor3
    real(kind=8), dimension(1:3, 1:11) :: gravpln
    real(kind=8), dimension(1:3) :: sum_inJ2, acc_inJ2
    real(kind=8), dimension(1:11) :: sumpln  
    real(kind=8), dimension(1:11) :: rpln, rpln2, rpln3
    real(kind=8), dimension(1:11) :: rplnb, rplnb2, rplnb3
    real(kind=8) :: C20     

!  sbcor(i) i  coordinates of earth satellite relative to earth
!  bcor(i)  i  coordinates of earth satellite relative to the sun
!  ccor(i)  i  coordinates of earth relative to the sun (should be i=1,6)
!  pccor(i) i  coordinates of moon relative to earth
!  pbcor(i) i  coordinates of moon relative to earth satellite
!  plncor(i,j) i  coordinates of j planets relative to earth
!  plnbcor(i,j) i  coordinates of j planets relative to earth satellite

!The following is copied from the JPL planetary ephemeris subroutine PLEPH 
!C            THE NUMBERING CONVENTION FOR 'NTARG' AND 'NCENT' IS:
!C
!C                1 = MERCURY           8 = NEPTUNE
!C                2 = VENUS             9 = PLUTO
!C                3 = EARTH            10 = MOON
!C                4 = MARS             11 = SUN
!C                5 = JUPITER          12 = SOLAR-SYSTEM BARYCENTER
!C                6 = SATURN           13 = EARTH-MOON BARYCENTER
!C                7 = URANUS           14 = NUTATIONS (LONGITUDE AND OBLIQ)
!C                            15 = LIBRATIONS, IF ON EPH FILE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    arc_frame = "J2000"

    tdtoff = 51.184d0

    s = 0.d0

! the following removed - not necessary because the PEPtime is not required for the
! JPL ephemeris call
!    call timinc(jd, t, s) ! to make sure day is increased at t=86400
!    time = jd +t/86400.d0 ! PEP time passed to this subroutine

    ephtime(1) = dble(jdGPS)
    ephtime(2) = (tGPS+tdtoff)/86400.d0 ! adding tdtoff converts from GPS time to Ephemeris time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   DETERMINE coor of earth relative to sun
    ispeed=0
!    call solred(ispeed,time+tdtoff/86400.d0,ccor) ! PEP time
!    the above call for GAMIT subroutine to calculate solar coordinates
!    has been replaced with JPL ephemeris call below
    call dpleph (ephtime, 11, 3, plncor(1:6,11)) ! pln = sun = 11
!   The dpelph call returns the position of the sun relative to the earth, 
!   but we need the earth relative to the sun (hence the negative in the 
!   following expression):
    ccor(1:3) = -1.d0 * plncor(1:3,11)*1000.d0 ! in metres
    sunpos(1:6) = -1.d0 * plncor(1:6,11)*1000.d0  ! sunpos contains the pos and vel of the 
                                                  ! earth relative to the sun for the 
                                                  ! subroutine "generalrel" 
    rc2=ccor(1)**2+ccor(2)**2+ccor(3)**2
    rc =dsqrt(rc2)
    rpln(11) = rc
    rc3=rc2*rc
    do i=1,3
      ccor3(i)=ccor(i)/rc3
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   DETERMINE coor of moon relative to earth
    ispeed=0
!    call lunred(ispeed,time+tdtoff/86400.d0,pccor) ! PEP time
!    the above call for GAMIT subroutine to calculate lunar coordinates
!    has been replaced with JPL ephemeris call below
    call dpleph (ephtime, 10, 3, plncor(1:6,10)) ! pln = moon = 10 
    pccor(1:3) = plncor(1:3,10)*1000.d0  ! convert km into metres
    rpc2  =   pccor(1)**2+pccor(2)**2+pccor(3)**2
    rpc   =   dsqrt(rpc2)
    rpln(10) = rpc
    rpc3  =   rpc2*rpc
    do i=1,3
      pccor3(i)=pccor(i)/rpc3
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   DETERMINE coor of the planets relative to the earth
!    print*, "ephtime", ephtime
    do j = 1,8
      if (j.ne.3) then
        call dpleph (ephtime, j , 3, plncor(1:6,j)) 
        plncor(1:6,j)=plncor(1:6,j)*1000.d0 ! convert km into metres
        rpln2(j) = plncor(1,j)**2+plncor(2,j)**2+plncor(3,j)**2
        rpln(j) = dsqrt(rpln2(j))
        rpln3(j) = rpln(j)*rpln2(j)
        do i = 1,3
          plncor3(i,j)=plncor(i,j)/rpln3(j)
        enddo
      endif
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   DETERMINE coor of earth satellite relative to sun
    do i=1,3
!     bcor = (earth wrt sun) + (satellite wrt earth)
      bcor(i)=ccor(i)+sbcor(i)
    enddo
    rb2=bcor(1)**2+bcor(2)**2+bcor(3)**2
    rb=dsqrt(rb2)
    rb3=rb2*rb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   DETERMINE coor of moon relative to earth satellite
    do i=1,3
!     pbcor = (moon wrt earth) - (satellite wrt earth)
      pbcor(i)=pccor(i)-sbcor(i)
    enddo
    rpb2  =   pbcor(1)**2+pbcor(2)**2+pbcor(3)**2
    rpb   =   dsqrt(rpb2)
    rpb3  =   rpb2*rpb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   DETERMINE coor of planets relative to earth satellite
    do j = 1,8 
      if (j.ne.3) then
        do i=1,3
!         pbcor = (planet rel earth) - (satellite rel earth)
          plnbcor(i,j)=plncor(i,j)-sbcor(i)
        enddo

        rplnb2(j)  =   plnbcor(1,j)**2+plnbcor(2,j)**2+plnbcor(3,j)**2
        rplnb(j)   =   dsqrt(rplnb2(j))
        rplnb3(j)  =   rplnb2(j)*rplnb(j)
      endif
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   CALCULATE the indirect J2 effect - see the GOCE standards document:
!   Doc. No: GO-TN-HPF-GS-0111 
!   Assume that this formulation represents the effective acceleration on the
!   satellite wrt the earth 

    sum_inJ2 = 0.d0 ! initialise sum

!   SUM UP contributions to indirect J2 effect from moon (10) and sun (11)
    do j=10, 11 ! for the moon and sun respectively
! PT130414: bounds checking when compiling showed up that gm(j+1) was beyond the array.
!           gm(1) is the Earth, gm(2) moon and gm(3) sun. Changed to gm(j-8)
! PT130830: according to the GOCE document (p33), the plncor(3,j)/rpln(j) should be squared !!
      sum_inJ2(1) = sum_inJ2(1) + (gm(j-8)/rpln(j)**3) * (Ae/rpln(j))**2 &
                      * (5.d0*(plncor(3,j)/rpln(j))**2 -1) * plncor(1,j)

      sum_inJ2(2) = sum_inJ2(2) + (gm(j-8)/rpln(j)**3) * (Ae/rpln(j))**2 &
                      * (5.d0*(plncor(3,j)/rpln(j))**2 -1) * plncor(2,j)

      sum_inJ2(3) = sum_inJ2(3) + (gm(j-8)/rpln(j)**3) * (Ae/rpln(j))**2 &
                      * (5.d0*(plncor(3,j)/rpln(j))**2 -3) * plncor(3,j)
    enddo

!   SCALE indirect J2 contribution (from moon + sun) by the fully normalised
!   second zonal harmonic coefficient (C20)
    do i=1,3
      acc_inJ2(i) = -3.d0 * dsqrt(5.d0) / 2.d0 * C20 * sum_inJ2(i)
    enddo 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 1, 3
!     Effect of the Moon on the motion of the Earth satellite
!     sump = (sat wrt moon)-(moon wrt earth)
      sump = pbcor(i)/rpb3-pccor3(i)
      gravm(i) = gm(2)*sump

!     Effect of the Sun on the motion of the satellite
!     sums = (earth wrt sun)-(satellite wrt sun)
      sums = ccor3(i)-bcor(i)/rb3
      gravs(i) =gm(3)*sums

!     Effect of the planets on the motion of the Earth satellite
!     sump = (sat wrt moon)-(moon wrt earth)
      do j = 1,8
        if (j.ne.3) then
          sumpln(j) = plnbcor(i,j)/rplnb3(j)-plncor3(i,j)
          gravpln(i,j) = gmpln(j)*sumpln(j)
        endif
      enddo

!     add all the contributions due to the moon, sun 
      gravsm(i) = gravm(i) + gravs(i) 

!     add the contributions due to the rest of the planets
      if (plnflag.eq.'Y') then
      do j = 1,8
        if (j.ne.3) then
          gravsm(i) = gravsm(i) + gravpln(i,j)
        endif
      enddo
      endif

!     adjust for the indirect J2 effect
      gravsm(i) = gravsm(i) - acc_inJ2(i)

    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    return 
    end











