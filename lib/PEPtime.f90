   subroutine PEPtime(jd, t, PEPjd, PEPt)

! Emma-Kate Potter, 23 March , 2011
! This subroutine uses the julian day (integer) and seconds of day to calculate 
! equivalent PEP jd and time 

    integer*4 :: jd, mjd, PEPjd
    real*8 :: t, jdate, mjdate, PEPt,PEPjdate

    jdate = dble(jd) +t/86400.d0      ! actual julian *date* (non-integer)
    mjdate = jdate-2400000.5d0        ! actual modified julian *date* (non-integer)

    mjd = aint(mjdate)                ! modified julian day (integer)
    PEPjdate = jdate + 0.5d0          ! PEP julian day as a real
    PEPjd = int(PEPjdate)
    PEPt = (PEPjdate-aint(PEPjdate))*86400.d0   ! seconds of day for PEP time (same as s-o-d for MJD)

! DEBUG:
!    print*,'PEPtime: jd,t,jdate,mjdate,mjd,PEPjd,PEPt',jd,t,jdate,mjdate,mjd,PEPjd,PEPt
   return
   end
