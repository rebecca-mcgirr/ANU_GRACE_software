  subroutine ymdhms_to_gracesec(date,sec,grace_sec)

! subroutine to convert from ymdhms to grace seconds (starting at 1200 on 01 01 2000)
!
! P. Tregoning
! 27 March 2017

  implicit none

  integer*4, intent(in)  :: date(5)      ! year, month, day, hr, min
  real*8   , intent(in)  :: sec          ! seconds of day
  integer*4, intent(out) :: grace_sec    ! grace seconds

! MJD of gracesec=0
  real*8, parameter :: grace_start_jd = 2451545.0d0  !  1200 on 1 Jan 2000
  real*8            :: jd


! first, compute the MJD of the input epoch
  call ymdhms_to_jd(date,sec,jd)  

! now get the grace seconds
  grace_sec = nint((jd - grace_start_jd)*86400.d0)

  return
  end


