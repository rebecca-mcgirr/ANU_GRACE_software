      subroutine gsec_to_ymdhms (epoch,date,seconds)

!----------------------------------------------------------------------
!
!     MOD SCM 110328: Created GRACE SEC to YMDHMS (cloned from JD_TO_YMDHMS).
!
!     Author: T.Herring           11:28 AM  THU., 10  JULY, 1986
!
!----------------------------------------------------------------------
!     Routine to convert a GRACE epoch in seconds to a calender data with
!     hours and minutes, and floating-point seconds (Real*8).
!
!     CALLING SEQUENCE:
!     =================
!     CALL GSEC_to_YMDHMS ( epoch, date, seconds )
!
!     WHERE:
!     epoch   GRACE eopch in seconds (REAL*8 INPUT)
!     date    is an array containing the calender date with (full) year,
!             month of year, day of month, hours, minutes. Note the year
!             will be returned with the centuries added.
!             (I*4 5 element array OUTPUT)
!     seconds is the floating point seconds part of the MJD
!             (REAL*8 OUTPUT)
!
!----------------------------------------------------------------------
!
!     century     - Number of century from 1600 Jan 0.
!     date(5)     - the calender date corresponding to
!                 - 'epoch' resolution to the minute
!     day         - day of month
!     day_of_year - Number of days from start of year
!     days_to_month(13)  - number of days to start of each month
!                 - in a non-leap-year
!     
!     month       - month of year
!     
!     year        - years since start of century
!     years_from_1600 - Number of years since 1600 Jan 0.
!     
!     days_from_1600  - Number of days elapsed since 1600 Jan 0.
!               - (MJD -94554.0 Julian date 2305447.0)

    integer*4 :: century,date(5),day,day_of_year,days_to_month(13),month,year,years_from_1600,days_from_1600

!     
!      epoch    - the julian date or modified julian date
!               - to be converted to calender date
!   fraction    - the fraction of a day part of MJD
!     
!   mjd         - epoch converted to a MJD
!   mjd_day     - the whole number days in the mjd
!     
!   seconds     - the seconds part of the MJD (<60)
     
    real*8 :: epoch, fraction, mjd, mjd_day, seconds, grace_epoch_mjd

! New variables 090923: Dates and second computed + 1usec when
!     seconds > 59.

    real*8 :: fracp,dsec
     
!   leap_year   - Indicates that this days is a leap year
     
    logical :: leap_year
     
    data days_to_month /    0,  31,  59,  90, 120, 151, 181, 212, 243, 273, 304, 334, 365 /

    data grace_epoch_mjd / 51544.5d0 /
     
!**** START, convert full jd to mjd if we have too.

    mjd = epoch/86400.d0 + grace_epoch_mjd

!**** Remove the fractional part of the mjd
     
    mjd_day  = aint ( mjd )
    fraction =  mod ( mjd,1.d0 )
     
    if( mjd < 0 .and. fraction /= 0.d0 ) then
        mjd_day  = mjd_day - 1
        fraction = fraction + 1.d0
    end if
     
!**** Now convert MJD (even day) to date (year, month, day )
!     Get number of days since 1600.
     
    days_from_1600 = mjd_day - (-94554.0d0)
    years_from_1600 = days_from_1600/365.d0
     
!**** Now compute day_of_year and years_from_1600 accounting for leap years
!     Just to get us into loop
    day_of_year = 0
    do while ( day_of_year.le.0 )
         
        century = years_from_1600/100
         
        day_of_year =  days_from_1600 - years_from_1600*365 &
        - (years_from_1600 -   1)/  4 &
        + (years_from_1600 +  99)/100 &
        - (years_from_1600 + 399)/400 - 1
         
!**** If we are 1600 then add one day
        if( years_from_1600 == 0 ) then
            day_of_year = day_of_year + 1
        end if
         
!**** See if the leap days have taken us to a earlier year
        if( day_of_year <= 0 ) then
            years_from_1600 = years_from_1600 - 1
        end if
    end do
     
!**** We now have number of days from start of year and the year
!     Convert years back to start of century    
    year = mod( years_from_1600, 100)
     
!**** See if this is a leap year     
    leap_year = .false.

!**** We are at beginning of century
    if( year == 0 ) then
        if( mod(century,4) == 0 ) leap_year = .TRUE. 
    else
        if( mod(year,4)   == 0 ) leap_year = .TRUE. 
    end if
     
!**** If the day of year is less than 60 then the leap years do no not
!     matter,  if we are greater than or equal to 60, need to account for the leap years
     
!**** Dont worry about leap years
    if( day_of_year < 60 ) then
        if( day_of_year <= 31 ) then
            month = 1
            day   = day_of_year

!**** We are in February
        else
            month = 2
            day   = day_of_year - 31
        end if

!**** Need to account for leap years
    else
        if( leap_year .and. day_of_year == 60 ) then
            month  = 2
            day    = 29
        else
            if( leap_year ) day_of_year = day_of_year - 1
             
!**** Now find month
            month = 2
            do while ( day_of_year.gt. days_to_month(month) )
                month = month + 1
            end do
            month = month - 1
            day   = day_of_year - days_to_month(month)
        end if
    end if
     
!**** Now save the date
     
    date(1) = years_from_1600 + 1600
    date(2) = month
    date(3) = day
     
!**** Now convert the fraction of a day to hours, minutes and seconds     
    date(4) = fraction*24.d0
    date(5) = fraction*1440.d0 - date(4)*60.d0     
    seconds = 86400.d0*fraction - date(4)*3600.d0 - date(5)*60.d0

    if( seconds >= 59.0d0 ) then
!**** MOD TAH 090923: Try calc with +1 usec (0.01d-9 days)
        dsec = 1.d-6   ! Seconds
        fracp = fraction + dsec/86400.d0

        date(4) = fracp*24.d0
        date(5) = fracp*1440.d0 - date(4)*60.d0         
        seconds = 86400.d0*fracp - date(4)*3600.d0 - date(5)*60.d0 - dsec
    end if
  
!**** Thats all
    return
  end subroutine gsec_to_ymdhms
