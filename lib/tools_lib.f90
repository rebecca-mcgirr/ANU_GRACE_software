!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! useful tool-type subroutines
!
! status_update   ::  prints output status information 
! bitmap          ::  decodes a bit-mapped number to see whether a particular bit is present
! trimlen       ::  trims down a character variable to the leading part before the first blank character
! get_keyword        ::  get the value associated with a command in an input command file
! amag3         ::  function to compute the magnitude of a three-dimensional vector
! rotmat        ::  create a rotation matrix for a rotation of "theta" radians around one of the three axes
! lwr_to_upper      ::  turn lower-case string into upper-case string
! last_nonblank         ::  return the position of the last non-blank character in a string
! xyz_to_geod   ::  converts cartesian XYZ into geodetic lat/lon/h
! ymdhms_to_jd  ::  year/month/day/hr/min/sec to julian day 
! ymd_to_jd     ::  convert year/month/day to day-of-year
! datetime        ::  form a runtime string of date and time
! first_word       ::  returns the first word of a string
! fund_arguments_iers1992 :: Brown's astronomical arguments for solid Earth tide computations
! jd_to_ymdhms  ::  convert julian day to ymdhms

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine status_update(state,program_name,subr_name,file_name,char_message,err_status)

! subroutine to print out status information. This is called by many programs and subroutines
!
! P. Tregoning
! 7 June 2021

! passed variables
  character(*),  intent(in)    :: state               ! status of reported message: "STATUS", "WARNING", "FATAL"
  character(*), intent(in)    :: program_name         ! name of program that is being run
  character(*), intent(in)    :: subr_name            ! name of subroutine within program that is calling for output info to be printed
  character(*), intent(in)    :: file_name            ! name of file related to output info to be printed (can be blank)
  character(*),  intent(in)    :: char_message        ! character message to be printed
  integer*4,     intent(inout) :: err_status          ! error status after doing whatever was required

! local variables
  integer*4       :: luout                            ! unit number of output file
  character*100   :: filename                         ! constructed output file name
  logical         :: exist
  character*500   :: output_message
  character*19    :: date_string
  character*150    :: file_string
  character*8     :: status_string                    ! 'STATUS'/'FATAL'/'WARNING' plus one blank character

! set blank the status string
  status_string = " "

! first, establish the output file name
  if(state(1:5) == "STATU")then
    filename = trim(program_name)//".status"
    luout = 101
  else if(state(1:5) == "WARNI")then
    filename = trim(program_name)//".warning"
    luout = 102
  else if(state(1:5) == "FATAL")then
    filename = trim(program_name)//".fatal"
    luout = 103
  else
    stop "unknown status_update print output state"
  endif

! now open the file
  open(luout, file=trim(filename), status="unknown", position="append")

! get the date and time
  call datetime(date_string)

! make a filename string, if needed
  if(len(trim(file_name)) > 0)then
    file_string = " (file "//trim(file_name)//")"
  else 
    file_string = " "
  endif

! create the status string
  status_string = trim(state)

! create the output message
  write(output_message,'(a,a,a2,a,a1,a,a2,a,a)')status_string,trim(date_string),": ",trim(program_name),"/",trim(subr_name) &
                                                 ,": ",trim(char_message),trim(file_string)
! write the message to the file and the screen
  write(*,'(a)')trim(output_message)
  write(luout,'(a)')trim(output_message)


! close the file
  close(luout)

  if(trim(status_string) == "FATAL")stop

  return
  end subroutine status_update
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical function bitmap(int_value, bit_value)

! logical function to decode a bit-mapped value into the relevant components and then decide whether 
! the requested one was present.
!
! The lowest bit is 1
!
! P. Tregoning
! 7 June 2021
  implicit none

! passed arguments
  integer*4    :: int_value(*)        ! bit-mapped input number
  integer*4    :: bit_value           ! required "bit" to test for


! local variables
  integer*4  :: pwr2_array(32)                   ! the array of powers of 2     
  integer*4 ia, ib                               ! array index and bit index

  data pwr2_array /    1,      2,       4,       8,      16,    32 , 64,    128,     256,     512,    1024,  2048,     &
                  4096,   8192,   16384,   32768,   65536, 131072, 262144,     524288,   1048576,  2097152, 4194304,   &
               8388608,   16777216,  33554432, 67108864, 134217728,  268435456, 536870912 ,  1073741824,-2147483648   /

!     1. Decompose IBIT into an array index IA and a bit index IB.
  ia = (bit_value+31)/32
  ib = bit_value - (ia-1)*32

! Test whether the required bit was present
  bitmap = IAND(int_value(ia),pwr2_array(ib)) .NE. 0

  return
  end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer function trimlen(char_string)

! returns the value of the last leading character of a string before the first blank character
!
! P. Tregoning
! 7 June 2021

  implicit none

  character(*) :: char_string
  trimlen = index(char_string," ")-1

  return
  end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_keyword(ifile,word_key,word_cmd,found_key,job)

! search for a character match in a line read from an input command file
!
! P. Tregoning
! 7 June 2021
  implicit none

! passed arguments
  integer*4     :: luin                            ! unit number of input file
  character*(*) :: word_key,word_cmd               ! keyword command, command value found
  integer*4     :: found_key                       ! -1: not found, 0: found but no value, >0 : value returned    
  integer*4     :: job                             ! 1: rewind and search for key word, 2: rewind and stop on key word, 3: search without rewind
 
! local variables
  character*256 :: line
  integer*4     :: cmd_len                          ! number of characters in the command string to search for
  integer*4    :: i,length,ifile,mchkey,strjst
  integer      :: last_nonblank,ioerr
  logical      :: blank_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ioerr = 0
!  set word_cmd length to -1 to indicate default wkey not found
  found_key = -1
  word_cmd = " "
  cmd_len = index(word_key," ")-1
  length = 0

  ! search key word from the begining
   rewind(ifile)
   line = " "
   do while (found_key == -1 .and. ioerr == 0)
     read(ifile,'(a)',end=200,iostat=ioerr) line  
     if(ioerr == 0 .and. line(1:1) == " ")then
       cmd_len = index(line,":")
       ! find the first non-blank character
       blank_char = .true.
       i = 2
       do while (blank_char .and. i < cmd_len)
         if(line(i:i) /= " ")then
           blank_char = .false.
         else
           i=i+1
         endif
           
       enddo

       ! is it the right command?
       if(line(i:cmd_len-1) == word_key)then
         ! find the last non-blank character in the line
         do i=256,cmd_len,-1
           if(line(i:i) /= " " .and. length == 0)then
             length = i
           endif
         enddo
         if(length == cmd_len) then  ! the line was blank after the command and colon
           found_key = 0
           return
         else
           word_cmd = line(cmd_len+1:length)
           found_key = length
           return
         endif
       endif

    endif
  enddo
 
200 return
    end subroutine get_keyword
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  double precision function amag3(vector)

! computes the magnitude of a R*8 three-dimensional vector
!
! P. Tregoning
! 8 June 2021

  implicit none

! passed variables
  real(kind=8), intent(in)  :: vector(3)   ! input vector

  amag3 = dsqrt(vector(1)**2+vector(2)**2+vector(3)**2)
  return
  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rotmat(theta,axis,rot_matrix)

! subroutine to establish a rotation matrix for a rotation of theta (in radians) around axis "axis"
!
! P. Tregoning
! 8 June 2021

  implicit none

  real(kind=8), intent(in)  :: theta                   ! rotation angle (in radians)
  integer*4,    intent(in)  :: axis                    ! axis about which rotation is required
  real(kind=8), intent(out) :: rot_matrix(3,3)         ! the derived rotation matrix

! local variables
  integer*4 i,j
  real*8 angle,r,sinang,cosang

! set all elements to zero
  rot_matrix = 0.d0

! now fill out the required values
  rot_matrix(axis,axis) = 1.d0      ! 1 on the diagonal for this axis

  i=MOD(axis,3) + 1
  j=MOD(j,3) + 1
  rot_matrix(i,i)=dcos(theta)
  rot_matrix(j,j)=dcos(theta)
  rot_matrix(i,j)=dsin(theta)
  rot_matrix(j,i)=-dsin(theta)

  return
  end subroutine rotmat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine lwr_to_upper(char_string)

! routine to convert a character string to upper case. Based on routine by Clive Page
!
! P. Tregoning
! 8 June 2021

  implicit none

! passed variables
  character*(*), intent(inout) :: char_string

! local variables
  integer*4  :: string_length
  integer*4  :: i,j

! get length of input string
  string_length = index(char_string," ")-1

! now convert lower to upper, if necessary
  do i = 1, string_length
    j = iachar(char_string(i:i))
    if (j>= iachar("a") .and. j<=iachar("z") ) then
      char_string(i:i) = achar(iachar(char_string(i:i))-32)
    end if
  end do
  return

  end subroutine lwr_to_upper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer function last_nonblank(string_in)

! function to return the position of the last non-blank character in a string
!
! P. Tregoning
! 9 June 2021

  implicit none

! passed variables
  character*(*) :: string_in

! local variables
  integer*1 string_length,i

  string_length = len(string_in)

! work backwards until the first non-blank character is found
  do i=string_length,1,-1
    if(string_in(i:i) /= " ")cycle
  enddo

  last_nonblank = i

  return
  end 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ymdhms_to_jd(date_array,sec,jd)

! subroutine to convert ymdhms to julian day
!
! P. Tregoning
! 9 June 2021

  implicit none

! passed variables
  integer*4, intent(in)    :: date_array(5)       ! yr/mo/dd/hr/min
  real(kind=8), intent(in) :: sec
  real(kind=8), intent(out):: jd                  ! output julian day (date + seconds)

! local variables
  integer*4  :: days_before_month(12)                ! number of days in the year before the 1st of each month
  integer*4  :: year,month,day,hr,minute
  integer*4  :: total_leap_days                      ! number of leap days since 1600
  integer*4  :: days_since_1600                      ! self explanatory

  data  days_before_month /0,31,59,90,120,151,181,212,243,273,304,334 /
  year = date_array(1)
  month = date_array(2)
  day = date_array(3)
  hr = date_array(4)
  minute = date_array(5)

! number of leap days since 1600 to the year before the requested year 
  total_leap_days = ((year-1600)-1)/4 - ((year-1600) + 99)/100 + ((year-1600) + 399)/400  + 1
  if( (year-1600).eq.0 ) total_leap_days = total_leap_days - 1
  
! number of days since 1600
  days_since_1600 = (year-1600)*365  + total_leap_days + days_before_month(month) + day
 
! add another day for the current year after February if it is a leap year
  if(mod((year-1600),4) == 0 .and. (mod((year-1600),100) /= 0 .or. mod((year-1600),400) == 0) .and. month > 2) then
    days_since_1600 = days_since_1600 + 1
  endif
 
! The midified Julian Day of 1600 Jan 0 is -94554
  jd = 2400000.5d0 - 94554.d0 + days_since_1600 + hr/24.d0 + minute/1440.d0 + sec/86400.d0 


  return
  end subroutine ymdhms_to_jd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ymd_to_doy(input_ymd,doy)

! convert year/month/day to day of year, taking into account leap years
!
! P. Tregoning 
! 17 June 2021

  implicit none

! passed variables
  integer*4, intent(in)   :: input_ymd(3)      ! year, month, day of month
  integer*4, intent(out)  :: doy               ! day of year

! local variables
  integer*4   :: days_to_month(12)
  integer*4   :: leap_day

! number of days of the year prior to the first day of each month (without leap year consideration)
  data  days_to_month /0,31,59,90,120,151,181,212,243,273,304,334 /

! do we need a leap year correction
  if (mod(input_ymd(1),4) == 0)then
    if(input_ymd(2) > 2 .or. ( input_ymd(2) == 2 .and. input_ymd(3) == 29))then
      leap_day = 1
    else
      leap_day = 0
    endif
  else
    leap_day = 0
  endif

! calculate day of year
  doy = days_to_month(input_ymd(2)) + input_ymd(3) + leap_day


  return
  end subroutine ymd_to_doy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine xyz_to_geod(rotmat,XYZ,llh)

! subroutine to convert cartesian XYZ into geodetic lat/lon/ell_height
!
! P. Tregoning
! 10 June 2021

  use  gm_mod   ! provides the WGS84 Earth radius and flattening values, and declaration of pi

  implicit none

! passed variables
  real(kind=8), intent(out) :: rotmat(3,3)     ! rotation matrix from XYZ to NEU
  real(kind=8), intent(in ) :: XYZ(3)          ! input cartesian XYZ coordinates
  real(kind=8), intent(out) :: llh(3)          ! output geodetic lat/lon/ell_height

! local variables
  real(kind=8),parameter :: e2 = 2.d0*earthflat_wgs84 - earthflat_wgs84**2   
  logical      :: converged  
  real(kind=8) :: rad_lat,rad_curv,rad_equatorial,lat_tmp,ell_h_tmp
  real(kind=8) :: twopi
  twopi = 8.d0*datan(1.d0)              

!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!
! compute the longitude
  llh(2) = datan2(XYZ(2),XYZ(1))
  if(llh(2) < 0)llh(2) = twopi + llh(2)

! now, iterate the lat/ell_h until convergence
  ! starting values
  rad_equatorial = dsqrt(XYZ(1)**2 + XYZ(2)**2)
  llh(1) = datan2(XYZ(3),rad_equatorial)
  llh(3) = 0.d0

  ! loop until convergence
  converged = .false.
  do while (.not. converged)
    ! radius of curvature
    rad_curv = earthrad_wgs84 / dsqrt(1.d0 - e2*dsin( llh(1) )**2)
    rad_lat  = rad_equatorial * (1.d0 - e2*rad_curv/(rad_curv+llh(3)) )

    ! new lat and ell_h
    lat_tmp = datan2(XYZ(3),rad_lat)
    if(dabs(lat_tmp) < twopi/8.d0)then
      ell_h_tmp = rad_lat/dcos(lat_tmp)
    else
      ell_h_tmp = XYZ(3)/dsin(lat_tmp) - (1.d0 - e2) * rad_curv
    endif

    ! check for convervence
    if(dabs(lat_tmp-llh(1)) < 1.d-4 .and. dabs(ell_h_tmp-llh(3)) < 1.d-4)then  ! convergence is < 0.1 mm
      converged = .true.
    endif
    llh(1) = lat_tmp
    llh(3) = ell_h_tmp
  enddo

! lat seems to be co-latitude. Convert to latitude
  llh(1) = twopi/4.d0 - llh(1)

  return
  end subroutine xyz_to_geod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine datetime(string)

! form a 19-char string of date and time
!
! P. Tregoning
! 10 June 2021

  implicit none

! passed variables
  character*19, intent(out) :: string

! local variables
  integer*4 :: tmpdate(3),tmptime(3)

! first, get the date
  call idate(tmpdate)
  call itime(tmptime)

! now form the YYMMDD:HHMM:SS.S character string
  string = " "
  write(string,'(3i2.2,a1,2i2.2,a1,f4.1)')tmpdate(3)-2000,tmpdate(2),tmpdate(1),":",tmptime(1:2),":",dble(tmptime(3))

! ensure that the seconds has no blank
  if (string(13:13) == " ")write(string(13:13),'(i1.1)')0


  return
  end subroutine datetime
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine first_word(string, word, next_blank)

! subroutine to extract the next word in a line, starting from "next_blank". 
!
! P. Tregoning
! 15 June 2021

  implicit none

! passed variables
  character(*), intent(in)    :: string        ! input character string
  integer*4,    intent(inout) :: next_blank    ! input: place to start in the string; Output: first blank after next word
  character(*), intent(out)   :: word          ! word found

! local variables
  integer*4     :: string_len                     ! length of input string
  integer*4     :: i

! get the length of the input string
  string_len = len(string)

! get the first non-blank character on or after the input "next_blank"
  i = next_blank
  do while (string(i:i) == " ")
    i = i+1
  enddo

! now find the next occurrence of a blank character, which will be after the end of the word
  next_blank = i
  do while (string(next_blank:next_blank) /= " ")
    next_blank = next_blank + 1
  enddo

! finally, the word is found as string(i:next_blank - 1)
  word = string(i:next_blank - 1)

  return 
  end subroutine first_word

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fund_arguments_iers1992( epoch, fund_arg )
!
! subroutine to compute the fundamental arguments (l, l', F, D, omega, L). Equations from page 32 of IERS1992 standards
!
! P. Tregoning
! 15 June 2021

  implicit none

! passed variables
  real(kind=8), intent(in)   :: epoch       ! required julian day (including fractional part of day)
  real(kind=8), intent(out)  :: fund_arg(6) ! computed values. fund_arg(6) is GMT + pi. All in radians

! local variables
  real(kind=8)    ::  t                     ! Julian centuries of 36525 days of 86400 seconds since J2000
  real(kind=8)    :: moon_anomaly(5),Sun_anomaly(5),F(5),D(5),omega(5)
  real(kind=8)    :: secs_in_360            ! number of seconds in 360 degrees
  real(kind=8)    :: pi
  real(kind=8)    :: univ_to_sidereal       ! ratio of universal to sidereal time
  real(kind=8)    :: fractional_day         ! fractional part of the Julian day

! provide precomputed  coefficients for the equations
  data moon_anomaly    /0.064d0, 31.310d0, 715922.633d0, 485866.733d0, 1325.0d0 /
  data Sun_anomaly   /-0.012d0, -0.577d0, 1292581.224d0, 1287099.804d0, 99.0d0 /
  data F     /0.011d0, -13.257d0, 295263.137d0, 335778.877d0, 1342.0d0/
  data D     /0.019d0, -6.891d0, 1105601.328d0, 1072261.307d0, 1236.0d0/
  data omega    / 0.008d0, 7.455d0, -482890.539d0, 450160.280d0, -5.0d0/
  
  pi = 4.d0*datan(1.d0)
  secs_in_360 = 360.d0*60.d0*60.d0

! convert "epoch" to "t"
  t = (epoch - 2451545.d0)/36525.d0

! l:     mean anomaly of the moon
  fund_arg(1) =  moon_anomaly(1)*t**3 + moon_anomaly(2)*t**2 + moon_anomaly(3)*t &
               + moon_anomaly(4) + mod( moon_anomaly(5) * t, 1.d0 ) * secs_in_360
  fund_arg(1) = mod(fund_arg(1),secs_in_360)

! l':    mean anomaly of the Sun
  fund_arg(2) = Sun_anomaly(1)*t**3 + Sun_anomaly(2)*t**2 + Sun_anomaly(3)*t &
               + Sun_anomaly(4) + mod( Sun_anomaly(5) * t, 1.d0 ) * secs_in_360
  fund_arg(2) = mod(fund_arg(2),secs_in_360)

! F:     mean longitude of the Moon - omega
  fund_arg(3) = F(1)*t**3 + F(2)*t**2 + F(3)*t + F(4) + mod( F(5) * t, 1.d0 ) * secs_in_360
  fund_arg(3) = mod(fund_arg(3),secs_in_360)

! D:     mean elongation of the Moon from the Sun
  fund_arg(4) = D(1)*t**3 + D(2)*t**2 + D(3)*t + D(4) + mod( D(5) * t, 1.d0 ) * secs_in_360
  fund_arg(4) = mod(fund_arg(4),secs_in_360)

! omega: mean longitude of the ascending node of the Moon
  fund_arg(5) = omega(1)*t**3 + omega(2)*t**2 + omega(3)*t + omega(4) + mod( omega(5) * t, 1.d0 ) * secs_in_360
  fund_arg(5) = mod(fund_arg(5),secs_in_360)

! convert from arcseconds to radians
  fund_arg(1:5) = fund_arg(1:5) / (3600.d0 * pi/180.d0)


! now compute Apparent Greenwich Sidereal Time (GST) (page 30 of IERS1992 standards)
  fractional_day = epoch - aint(epoch-0.5d0) + 0.5d0
  t = ( aint(epoch-0.5d0) + 0.5d0 - 2451545.d0)/36525.d0
  

  fund_arg(6) = ( 24110.54841d0  + 8640184.81266d0*t + 0.093104d0*t**2 - 6.2d-6*t**3 ) /86400.d0

  ! ratio of universal to sidereal time (p30 of IERS1992
  univ_to_sidereal = 1.002737909350795d0 + 5.9006d-11*t - 5.9d-15*t**2

! add the fractional part of the day and then convert to radians 
  fund_arg(6) = (fund_arg(6) + univ_to_sidereal*fractional_day) * 2.d0*pi


  return
  end subroutine fund_arguments_iers1992
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!TITLE TIDE_ANGLES

      subroutine tide_angles_gamit( epoch, fund_arg )

!     Routine to compute the value of the fundamental argument
!     for Brown's arguments.  The sixth entry is returned as GST
!     plus pi.  The additional pi is used for compatability with
!     Doodson's Tide argument.

! PHYSICAL CONSTANTS NEEDED FOR SD_COMP

!   pi          - Define here to full precision
!   rad_to_deg  - Conversion from radians to degs.
!   DJ2000      - Julian date of J2000
!   sec360      - number of seconds in 360 degreees.

      real*8 pi, rad_to_deg, DJ2000, sec360

      parameter ( pi            = 3.1415926535897932D0 )
      parameter ( DJ2000        = 2451545.d0           )
      parameter ( sec360        = 1296000.d0           )

!     Computed quanities
      parameter ( rad_to_deg    = 180.d0   /pi         )

!-------------------------------------------------------------------

! PASSED VARIABLES

! INPUT
! epoch  - Julian date for arguments (jd + fraction of day)

! OUTPUT
! fund_arg(6) -  Brown's arguments plus GST+pi (rads)
!
      real*8 epoch, fund_arg(6)

! LOCAL VARIABLES
!      cent             - Julian centuries to DJ2000.
!      el,eld           - Mean longitude of moon minus mean
!                       - longitude of moon's perigee (arcsec)
!      elc(5)           - Coefficients for computing el
!      elp,elpd         - Mean longitude of the sun minus mean
!                       - longitude of sun perigee (arcsec)
!      elpc(5)          - Coeffiecents for computing elp
!      f,fd             - Moon's mean longitude minus omega (sec)
!      fc(5)            - Coefficients for computing f
!      d,dd             - Mean elongation of the moon from the
!                       - sun (arcsec)
!      dc(5)            - coefficients for computing d
!      om,omd           - longitude of the ascending node of the
!                       - moon's mean orbit on the elliptic
!                       - measured from the mean equinox of date
!      omc(5)           - Coefficients for computing om.
!      gst              - Greenwich mean sidereal time (rad)

      real*8 cent, el,eld, elc(5), elp, elpd, elpc(5), &
         f,fd, fc(5), d,dd, dc(5), om,omd, omc(5), gst

!      fract        - fraction of a day from 0:00 hrs UT.
!      Jd_0hr       - Julian date at zero hours UT
!      t_0hr        - Days since DJ2000 at 0:00 hrs UT
!      gstd         - GMST at 0:00 hrs UT1 of day being evaluated
!      diurnv       - Ratio of solar days to sidreal days on
!                     day of evalution.

      real*8 fract, t_0hr, gstd, diurnv, jd_0hr

!***  DATA statements for the fundamental arguments.

      data elc    /     0.064d0,    31.310d0,    715922.633d0, &
                  485866.733d0,    1325.0d0 /
      data elpc   /    -0.012d0,    -0.577d0,   1292581.224d0, &
                 1287099.804d0,      99.0d0 /
      data fc     /     0.011d0,   -13.257d0,    295263.137d0, &
                  335778.877d0,    1342.0d0/
      data dc     /     0.019d0,    -6.891d0,    1105601.328d0, &
                 1072261.307d0,    1236.0d0/
      data omc    /     0.008d0,     7.455d0,    -482890.539d0, &
                  450160.280d0,      -5.0d0/

!***  Get the number of centuries to current time

      cent = (epoch-dj2000) / 36525.d0

!***  Compute angular arguments
      el = elc(1) * cent**3 + elc(2) * cent**2 + elc(3) * cent &
               + elc(4) + mod( elc(5) * cent, 1.d0 ) * sec360
      el = mod( el, sec360 )
      eld = 3.d0 * elc(1) * cent**2 + 2.d0 * elc(2) * cent + elc(3) &
           + elc(5) * sec360
!
      elp = elpc(1) * cent**3 + elpc(2) * cent**2 + elpc(3) * cent &
          + elpc(4) + mod( elpc(5) * cent, 1.d0 ) * sec360
      elp = mod( elp, sec360 )
      elpd = 3.d0 * elpc(1) * cent**2 + 2.d0 * elpc(2) * cent + elpc(3) &
            + elpc(5) * sec360
!
      f = fc(1) * cent**3 + fc(2) * cent**2 + fc(3) * cent &
          + fc(4) + mod( fc(5) * cent, 1.d0 ) * sec360 
      f = mod( f, sec360 )
      fd = 3.d0 * fc(1) * cent**2 + 2.d0 * fc(2) * cent + fc(3) &
          + fc(5) * sec360
!
      d = dc(1) * cent**3 + dc(2) * cent**2 + dc(3) * cent &
          + dc(4) + mod( dc(5) * cent, 1.d0 ) * sec360
      d = mod( d, sec360 )
      dd = 3.d0 * dc(1) * cent**2 + 2.d0 * dc(2) * cent + dc(3) &
          + dc(5) * sec360
!
      om = omc(1) * cent**3 + omc(2) * cent**2 + omc(3) * cent &
          + omc(4) + mod( omc(5) * cent, 1.d0 ) * sec360
      om = mod( om, sec360 )
      omd = 3.d0 * omc(1) * cent**2 + 2.d0 * omc(2) * cent + omc(3) &
           + omc(5) * sec360
!

!**** Now compute GMST.  (CALC 7.1 Algorithm)
!     Remove the fractional part of the julian date
!     Get jd at 0:00 UT
      jd_0hr = aint(epoch-0.5d0) + 0.5d0
!                         ! Days since J2000.0
      t_0hr = jd_0hr - dj2000
!                         ! 0:00 hrs at start of day
      cent = t_0hr / 36525.d0

!                         ! Fraction of a day
      fract = epoch - jd_0hr

      diurnv = ( 1.002737909350795d0 + 5.9006d-11*cent &
                                    - 5.9d-15*cent**2 )
!
!**** COMPUTE GST in cycles
      gstd = ( 24110.54841d0  + 8640184.81266d0*cent &
                             + 0.093104d0*cent**2   &
                             - 6.2d-6*cent**3 ) /86400.d0

      gstd = mod(gstd,1.d0)
!                                             ! Rads
      gst = (gstd + diurnv*fract) * 2.d0*pi

!***  Now save the values.  Convert values from arcseconds to radians

      fund_arg(1) = el / (3600.d0*rad_to_deg)
      fund_arg(2) = elp/ (3600.d0*rad_to_deg)
      fund_arg(3) = f  / (3600.d0*rad_to_deg)
      fund_arg(4) = d  / (3600.d0*rad_to_deg)
      fund_arg(5) = om / (3600.d0*rad_to_deg)
      fund_arg(6) = gst + pi

!**** Thats all
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!TITLE 'JD_TO_YMDHMS'
 
      SUBROUTINE JD_to_YMDHMS ( epoch, date, seconds )
!     ------------------------
!    .,      10JUL86                <871210.1317>

! MOD TAH 090923: Try to fix 60.0 s problem.
 
!     Author: T.Herring           11:28 AM  THU., 10  JULY, 1986
!
!$S "JD_to_YMDHMS"
!----------------------------------------------------------------------
!     Routine to convert a Julian date (with fractions of a day)
!     or a full Julian date, to a calender data with hours and minutes,
!     and floating-point seconds (Real*8).
!
!     NOTE: If a full Julian date is used the resolution of the seconds
!     will only be about 10 microseconds,  a MJD should yield a resolution
!     of about 0.1 microseconds in the seconds.
!
!     This routine is only valid for dates after 1600 Jan 0 (MJD -94554)
!
!     CALLING SEQUENCE:
!     =================
!     CALL JD_to_YMDHMS ( epoch, date, seconds )
!
!     WHERE:
!     epoch   is a modified or full julian date with possibly a fractional
!             part of a day.  The two type of input are destinguished by
!             the full Julian date being greater than 2,000,000.
!             (REAL*8 INPUT)
!     date    is an array containing the calender date with (full) year,
!             month of year, day of month, hours, minutes. Note the year
!             will be returned with the centuries added.
!             (I*4 5 element array OUTPUT)
!     seconds is the floating point seconds part of the MJD
!             (REAL*8 OUTPUT)
!
!----------------------------------------------------------------------
!
!$E
 
 
 
!         century     - Number of century from 1600 Jan 0.
!         date(5)     - the calender date corresponding to
!                     - 'epoch' resolution to the minute
!   day         - day of month
!   day_of_year - Number of days from start of year
!   days_to_month(13)   - number of days to start of each month
!               - in a non-leap-year
 
!   month       - month of year
 
!   year        - years since start of century
!   years_from_1600 - Number of years since 1600 Jan 0.
 
!   days_from_1600  - Number of days elapsed since 1600 Jan 0.
!               - (MJD -94554.0 Julian date 2305447.0)
      integer*4 century, date(5), day, day_of_year, days_to_month(13),month, year, years_from_1600, days_from_1600
 
!      epoch    - the julian date or modified julian date
!               - to be converted to calender date
!   fraction    - the fraction of a day part of MJD
 
!   mjd         - epoch converted to a MJD
!   mjd_day     - the whole number days in the mjd
 
!   seconds     - the seconds part of the MJD (<60)
 
      real*8 epoch, fraction, mjd, mjd_day, seconds

! New variables 090923: Dates and second computed + 1usec when
!     seconds > 59.
      real*8 fracp, dsec
 
!       leap_year   - Indicates that this days is a leap year
 
      logical leap_year
 
      data days_to_month /   0,  31,  59,  90, 120, 151, 181, 212, 243, 273, 304, 334, 365       /
 
 
!**** START, convert full jd to mjd if we have too.
!                                      ! epoch must be Full julian date
      if( epoch.gt.2000000.d0 ) then
!                                      ! convert to MJD.
 
          mjd = epoch - 2400000.5d0
      else
          mjd = epoch
      end if
 
!**** Remove the fractional part of the mjd
 
      mjd_day  = aint ( mjd )
      fraction =  mod ( mjd,1.d0 )
 
      if( mjd.lt.0 .and. fraction.ne.0.d0 ) then
          mjd_day  = mjd_day - 1
          fraction = fraction + 1.d0
      end if
 
!**** Now convert MJD (even day) to date (year, month, day )
!     Get number of days since 1600.
 
      days_from_1600 = mjd_day - (-94554.0d0)
      years_from_1600 = days_from_1600/365.d0
 
!**** Now compute day_of_year and years_from_1600 accounting for
!     leap years
 
!                        ! Just to get us into loop
      day_of_year = 0
      do while ( day_of_year.le.0 )
 
          century = years_from_1600/100
 
          day_of_year =  days_from_1600 - years_from_1600*365 - (years_from_1600-1)/4 + (years_from_1600 + 99)/100 &
                         - (years_from_1600 + 399)/400 - 1
 
!         If we are 1600 then add one day
          if( years_from_1600.eq.0 ) then
              day_of_year = day_of_year + 1
          end if
 
!         See if the leap days have taken us to a earlier year
          if( day_of_year.le.0 ) then
              years_from_1600 = years_from_1600 - 1
          end if
      end do
 
!**** We now have number of days from start of year and the year
!     Convert years back to start of century
 
      year    = mod( years_from_1600, 100)
 
!**** See if this is a leap year
 
      leap_year = .false.
!                             ! we are at beginning of century
      if( year.eq.0 ) then
          if( mod(century,4).eq.0 ) leap_year = .true.
      else
          if( mod(year,4)   .eq.0 ) leap_year = .true.
      end if
 
!**** If the day of year is less than 60 then the leap years do no not
!     matter,  if we are greater than or equal to 60, need to account
!     for the leap years
 
!                                     ! Dont worry about leap years
      IF( day_of_year.lt.60 ) THEN
          if( day_of_year.le.31 ) then
              month = 1
              day   = day_of_year
!                 ! we are in February
          else
              month = 2
              day   = day_of_year - 31
          end if
!                                     ! Need to account for leap years
      ELSE
          if( leap_year .and. day_of_year.eq.60 ) then
              month  = 2
              day    = 29
          else
              if( leap_year ) day_of_year = day_of_year - 1
 
!****         Now find month
              month = 2
              do while ( day_of_year.gt. days_to_month(month) )
                  month = month + 1
              end do
              month = month - 1
              day   = day_of_year - days_to_month(month)
          end if
      END IF
 
!**** Now save the date
 
      date(1) = years_from_1600 + 1600
      date(2) = month
      date(3) = day
 
!**** Now convert the fraction of a day to hours, minutes and seconds
 
      date(4) = fraction*24.d0
      date(5) = fraction*1440.d0 - date(4)*60.d0
 
      seconds = 86400.d0*fraction - date(4)*3600.d0 - date(5)*60.d0

      if( seconds.ge. 59.0d0 ) then
! MOD TAH 090923: Try calc with +1 usec (0.01d-9 days)
          dsec = 1.d-6   ! Seconds 
          fracp = fraction + dsec/86400.d0 

          date(4) = fracp*24.d0
          date(5) = fracp*1440.d0 - date(4)*60.d0
 
          seconds = 86400.d0*fracp - date(4)*3600.d0 - date(5)*60.d0 - dsec 
      end if     
        
 
!**** Thats all
      RETURN
      end subroutine JD_to_YMDHMS
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




