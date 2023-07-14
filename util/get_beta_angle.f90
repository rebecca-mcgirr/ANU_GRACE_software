  program get_beta_angle

! program to calculate the beta angle given a satellite position and velocity for a given date. Input command line
! arguments are year, day-of-year, x, y, z, vx, vy, vz, frame
!
! P. Tregoning
! 26 June 2013
!
! program will use similar code to gracefit, but does so as a stand-alone entity

  implicit none

  character arg*80,gnv1b_file*26,message*250,char1*1,char2*1,inframe*1,string*8
  integer*4 i,j,yr,doy,date(5),ioerr,iut1,ipole,inut,ncarg,imonth,iday,ihr,imin,isec
  double precision fjd,ephtime(2),sunpos_i(6),epoch,junk(3),PEPt,beta_ang_calc,pi
  double precision satpos_i(3),satvel_i(3),satpos_ef(3),satvel_ef(3),tmp(6)

! variables for the ef-inertial rotation
  integer*4 jd,PEPjd,sec
  double precision idir,iut1pol,tdtoff,rot(3,3),rotdot(3,3),sidtm,sod,seconds,decyrs

  pi = 4.d0*datan(1.d0)
  tdtoff = 51.184d0

! get the year and day-of-year from command line
  ncarg = command_argument_count()
  if (ncarg .eq. 0 .or. (ncarg .ne. 11 .and. ncarg .ne. 14) ) then
    print*,'Runstring: get_beta_angle 2013 177 43200 I x y z vx vy vz string'
    print*,'           get_beta_angle 2013 01 05 12 00 00 I x y z vx vy vz string'
    stop
  endif 
  call getarg(1,arg)
  read(arg,*)yr
  date(1) = yr
  if ( ncarg == 11 ) then
    call getarg(2,arg)
    read(arg,*)doy
    call getarg(3,arg)
    read(arg,*)sec
! convert the year, doy and seconds of day entered into fjd.
    call yds_to_jd(yr,doy,sec,fjd)
! convert fjd to ymdhms.
    call jd_to_ymdhms(fjd,date,seconds) 
    call getarg(4,arg)
    read(arg,'(a)')inframe
    do i=1,6
      j = i+4
      call getarg(j,arg)
      read(arg,*)tmp(i)
    enddo
    call getarg(11,arg)
    read(arg,'(a)')string     
  else
    call getarg(2,arg)
    read(arg,*)date(2)
    call getarg(3,arg)
    read(arg,*)date(3)
    call getarg(4,arg)
    read(arg,*)date(4)
    call getarg(5,arg)
    read(arg,*)date(5)
    call getarg(6,arg)
    read(arg,*)seconds
! convert the year, month day hr min and seconds of day entered into fjd.
    call ymdhms_to_jd(date,seconds,fjd)
! convert fjd to year, doy ,sec.
    call jd_to_yds ( fjd, yr, doy, sec )
    call fix_y2k(yr)
    call getarg(7,arg)
    read(arg,'(a)')inframe
    do i=1,6
      j = i+7
      call getarg(j,arg)
      read(arg,*)tmp(i)
    enddo
    call getarg(14,arg)
    read(arg,'(a)')string
  endif
  if (inframe .eq. 'E') then
    satpos_ef = tmp(1:3)
    satvel_ef = tmp(4:6)
  else
    satpos_i = tmp(1:3)
    satvel_i = tmp(4:6)
 endif
!  print*,'satpos_i satvel_i:',satpos_i,satvel_i
!  print*,'satpos_ef satvel_ef:',satpos_ef,satvel_ef
!  print*,'yr doy sec inframe:',yr,doy, sec, inframe
!  print*,'yr month day hr min sec inframe:',date, seconds, inframe
!  print*,'string',string

! open the pole., ut1. and nutabl. files.
! Unit=12 is used for the JPLEPH file
  iut1 = 11
  ipole = 13
  inut = 14
  open (unit=iut1,file="ut1.",status='old',iostat=ioerr)
  if (ioerr /= 0) then
      call status_update('FATAL','UTIL','beta_angle ','ut1.','Error opening ut1 table: ',ioerr)
  endif
  open (unit=ipole,file="pole.",status='old',iostat=ioerr)
  if (ioerr /= 0) then
      call status_update('FATAL','GRACEFIT','beta_angle ','pole.','Error opening pole table: ',ioerr)
  endif
  open (unit=inut,file="nutabl.",status='old',iostat=ioerr)
  if (ioerr /= 0) then
      call status_update('FATAL','UTIL','beta_angle ','nutabl.','Error opening nutation table: ',ioerr)
  endif

! get integer jd
 jd = int(fjd)

! get seconds of day 
  sod = (fjd - dble(jd)) * 86400.d0

! get date in decimal years 
 call jd_to_decyrs( fjd, decyrs )

! get the sun position for this epoch
  ephtime(1) = dble(jd)                         ! integer part of the epoch
  ephtime(2) = fjd - dble(jd) + tdtoff/86400.d0 ! fraction of a day
  call dpleph (ephtime, 11, 3, sunpos_i)
  sunpos_i = -1.d0*sunpos_i*1.d3    ! puts sunpos in units of m (I think!!) and with the same sign as in graceorb
!  print*,'sunpos_i:',sunpos_i,ephtime

! to define the orbital plane of the satellites, we need the position and velocity of the satellite in inertial coords.

  if ( inframe .eq. 'E' ) then
! now we need to rotate the earth-fixed pos/vel into inertial space. Get the rotation matrix
     idir = 1  ! for ef to inertial
     iut1pol = 7  ! the value used in graceorb
     call PEPtime(jd,sod, PEPjd, PEPt)
     call rotsnp(idir,PEPjd,PEPt,tdtoff,iut1pol,rot,rotdot,sidtm, &
                 iut1,ipole,inut,'J2000','IAU76') !PEPtime test
!    print*,"PEPjd,PEPt,fjd,sec",PEPjd,PEPt,jd,sod
!    call printmat(rot,3,3,"rot ")

! rotate the position and velocity of the satellite into inertial space
     call matmult(rot,satpos_ef,satpos_i,3,3,1)
     call matmult(rot,satvel_ef,satvel_i,3,3,1)
  endif

! now calculate the beta angle
  call calc_beta_angle(sunpos_i(1:3),satpos_i,satvel_i,beta_ang_calc)
!  print*,'sunpos_i,satpos_i,satvel_i,ephtime',sunpos_i,satpos_i,satvel_i,ephtime

  write(message,'(a8,1x,f14.9,1x,i4,1x,i3.3," :  ",i4,2("-",i2.2),1x,2(i2.2,":"),f5.2,a,f6.2)')string,decyrs,yr,doy,&
  (date(i),i=1,5),seconds," beta angle: ",beta_ang_calc*180.d0/pi
  call status_update('STATUS','UTIL','get_beta_angle',' ',message,0)

  stop

1000 write(message,'(a)')"No END OF HEADER line found in the GNV1B file."
     call status_update('FATAL','UTIL','util/get_beta_angle',' ',message,0)

  end
