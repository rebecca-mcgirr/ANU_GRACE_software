  program dist_deriv

! ************************************************************************************************
! program to calculate the time derivative of the distance function between two satellites
!
! A. Purcell
! 27 January 2022
! ************************************************************************************************

! ************************************************************************************************
  use gracefit_mod
 
  implicit none

  include '../gracefit/input.h90'

  double precision   :: dot_angle_calc           ! internal function

  character(len=200) :: line, tmp_char
  character(len=150) :: arg, message
  character(len=20)  :: out_angle, out_deriv, out_dist
  character(len=1)   :: sat, frame

  character(len=150), dimension(2) :: gnv        ! filenames of GNV1B input data


  double precision :: fact, tmp0, tmp1, tmp2, tmp3

  double precision, parameter :: tol = 1.0d-3

  double precision, dimension(60) :: angle, deriv, dist           ! output arrays
  double precision, dimension(3)  :: dvec,  dvec0, dvel, dvel0    ! delta-position and delta-velocity arrays

!------------------------------------------------------------------------------------
  ! Data read in from GNV1B files
  double precision, dimension(2) :: gnv_step     ! epoch step size for each GNV1B file

  double precision, allocatable :: epochs(:)     ! matrix of angles between LOS vectors

  double precision, allocatable :: gvec(:,:,:)     ! Positions and velocities for each satellite at every epoch
  double precision, allocatable :: gvec_epoch(:,:) ! Epochs of positions and velocities for each satellite at every epoch

  integer*4, dimension(2) :: num_gps_used        ! number of GPS observations that will be used (after decimating the 5sec data)
  integer*4, dimension(2) :: gnv_span            ! span of epochs of each GNV1B file ?
!------------------------------------------------------------------------------------

  integer*4, parameter :: lu_angle=12, lu_deriv=13, lu_dist=14

  integer*4, dimension(5) :: date

  integer*4 :: iepoch                            ! epoch counter
  integer*4 :: i, iline, irec, isat, j, tmp_int
  integer*4 :: ioerr


  logical   :: debug, tst1, tst2
! ************************************************************************************************


! ************************************************************************************************
  debug = .false.

  calling_prog = "dist_deriv"

  out_angle = "output.angle.txt"
  out_deriv = "output.deriv.txt"
  out_dist  = "output.dist.txt"

  fact = 1.8000d2/pi
! ************************************************************************************************


! ******************************************************************
! decode runstring
  call getarg(1,arg)

  if ( arg(1:1) == " ") then
    write(message,'(a)')"Runstring: dist_deriv 2019 10 02 mission "
    call report_stat('FATAL','UTIL','dist_deriv',' ',message,0)
  else
!   Do Nothing
  end if

! epoch: date(1) = year, date(2) = month, date(3) = day
  date = 0

  read(arg,*) date(1)

  call getarg(2,arg)
  read(arg,*) date(2)

  call getarg(3,arg)
  read(arg,*) date(3)

! mission
  call getarg(4,arg)
  read(arg,*) mission

  write(6,'(a,i5,2i3,/)') "Input date: ", date(1), date(2), date(3)
! ******************************************************************


! ******************************************************************
! form up the names of the GNV1B files
  if ( mission == 0 ) then
    write(gnv(1),'("data/GNV1B_",i4,"-",i2.2,"-",i2.2,"_A_02.asc")') date(1), date(2), date(3)
    write(gnv(2),'("data/GNV1B_",i4,"-",i2.2,"-",i2.2,"_B_02.asc")') date(1), date(2), date(3)
  else
    write(gnv(1),'("data/GNI1B_",i4,"-",i2.2,"-",i2.2,"_C_04.txt")') date(1), date(2), date(3)
    write(gnv(2),'("data/GNI1B_",i4,"-",i2.2,"-",i2.2,"_D_04.txt")') date(1), date(2), date(3)
  end if

  write(6,'(3(a,/))') "Input files: ", gnv(1), gnv(2)
! ******************************************************************


! ******************************************************************
! open the GNV1B files
  open(LUGN(1), file=gnv(1), status='old', iostat=ioerr)
  if ( ioerr /= 0 ) &
    call report_stat('FATAL','UTIL','dist_deriv',gnv(1),"Error opening GNV1B file",0)

  open(LUGN(2), file=gnv(2), status='old', iostat=ioerr)
  if ( ioerr /= 0 ) &
    call report_stat('FATAL','UTIL','dist_deriv',gnv(2),"Error opening GNV1B file",0)
! ******************************************************************


! ******************************************************************
  call ymdhms_to_gracesec(date, 0.d0, tmp_int)
  starting_epoch = dble(tmp_int)

  write(6,'(a,i5,2i3,f14.2)') "date: ", date(1), date(2), date(3), starting_epoch

! Data read in from GNV1B files
  allocate(gvec(6,2,nepochs_t))
  allocate(gvec_epoch(nepochs_t,2))

  do isat = 1, 2
    ioerr = 0
    iline = 0
    line  = "test"

    do while ( line(1:1) /= "#" )
      iline = iline + 1
      read(LUGN(isat),'(a200)',end=999,iostat=ioerr) line

      if ( ioerr /= 0 ) then
        write(6,'(2(a,/),a,i1,2a,/,2(a,i4,/),a,/)') &
          "++++++++++++++++++++++++++++++++++++++++++++", &
          "ERROR: dist_deriv",                            &
          "       error reading file: ", trim(gnv(isat)), &
          "       invalid header entry at line: ", iline, &
          "       error code: ", ioerr,                   &
          "++++++++++++++++++++++++++++++++++++++++++++"
        stop
      else
!       Do Nothing
      end if

      if ( index(line,"num_records: ") == 0 ) cycle

      irec = index(line,": ") + 1
      tmp_char = trim(line(irec:200))
      read(tmp_char,*) num_gps_used(isat)

    end do

    if ( num_gps_used(isat) == 0 ) then
      write(6,'(2(a,/),a,i1,2a,/,a,i1,a,/,a,/)') &
        "++++++++++++++++++++++++++++++++++++++++++++", &
        "ERROR: dist_deriv",                            &
        "       error reading file ", isat, ": ", trim(gnv(isat)), &
        "       num_gps_used(",isat,") = 0",         &
        "++++++++++++++++++++++++++++++++++++++++++++"
      stop
    else
!     Do Nothing
    end if

    write(6,'(2(a,/),a,i1,2a,/,a,/,a,i4,/,a,i1,a,i6,/)') &
      "============================================", &
      "REPORT: dist_deriv",                           &
      "        reading file ", isat, ": ", trim(gnv(isat)), &
      "        header successfully read",             &
      "        lines in header: ", iline,             &
      "        num_gps_used(", isat, ") = ", num_gps_used(isat)

    irec = 0
    do while ( ioerr == 0 )
      irec = irec + 1
      read(LUGN(isat),*,end=998,iostat=ioerr) gvec_epoch(irec,isat), sat, frame, &
        gvec(1,isat,irec), gvec(2,isat,irec), gvec(3,isat,irec), tmp1, tmp2, tmp3, &
        gvec(4,isat,irec), gvec(5,isat,irec), gvec(6,isat,irec)

      if ( ioerr /= 0 ) then
        write(6,'(2(a,/),a,i1,2a,/,2(a,i4,/),a,/)') &
          "++++++++++++++++++++++++++++++++++++++++++++", &
          "ERROR: dist_deriv",                            &
          "       error reading file: ", trim(gnv(isat)), &
          "       invalid data entry at record: ", irec,  &
          "       error code: ", ioerr,                   &
          "++++++++++++++++++++++++++++++++++++++++++++"
        stop
      else
!       Do Nothing
      end if

    end do

    if ( irec /= num_gps_used(isat) ) then
      write(6,'(2(a,/),a,i1,2a,/,2(a,i4,/),a,/)') &
        "++++++++++++++++++++++++++++++++++++++++++++", &
        "ERROR: dist_deriv",                            &
        "       error reading file: ", trim(gnv(isat)), &
        "       no. of records read: ", irec,           &
        "       no. of records expected: ", num_gps_used(isat), &
        "++++++++++++++++++++++++++++++++++++++++++++"
      stop
    else
      gnv_step(isat) = gvec_epoch(2,isat)    - gvec_epoch(1,isat)
      gnv_span(isat) = gvec_epoch(irec,isat) - gvec_epoch(1,isat)

      write(6,'(a,i5,a,/,a,/)') &
        "        ", irec, " GNV records read ", &
        "============================================"
    end if

  end do
! ******************************************************************

! ******************************************************************
! If dealing with GRACE data we need to convert to inertial coordinates
! Leave this for later
  if ( mission == 0 ) then
    write(6,'(4(a,/))') &
      "**************************************************************", &
      "PROGRAM dist_deriv WARNING",  &
      "Conversion of GRACE data to inertial coords not yet functional", &
      "Execution terminated", &
      "**************************************************************"
    stop
!    call earth_fixed2inert(,gvec)
  else
!   Do Nothing
  end if
! ******************************************************************


! ******************************************************************
! loop through all epochs
! lead satellite is assumed to be satellite 1
  angle = 0.00d0

  open(unit=lu_angle, file=out_angle)
  open(unit=lu_dist,  file=out_dist)
  open(unit=lu_deriv, file=out_deriv)

  do iepoch = 1, nepochs_t
    if ( ( sqrt(gvec(1,1,iepoch)**2 + gvec(2,1,iepoch)**2 + gvec(3,1,iepoch)**2) < tol ) .or. &
         ( sqrt(gvec(1,2,iepoch)**2 + gvec(2,2,iepoch)**2 + gvec(3,2,iepoch)**2) < tol ) ) cycle

    do i = 1, 3
      dvec0(j) = gvec(i,1,iepoch)   - gvec(i,2,iepoch)
      dvel0(j) = gvec(3+i,1,iepoch) - gvec(3+i,2,iepoch)
    end do

    do j = 1, 60
      angle(j) =  1.80d2
      deriv(j) =  9.99d6
      dist(j)  = -1.00d0

      if ( sqrt(gvec(1,2,iepoch+j)**2 + gvec(2,2,iepoch+j)**2 + gvec(3,2,iepoch+j)**2) < tol ) cycle

      do i = 1, 3
        dvec(j) = gvec(i,1,iepoch)   - gvec(i,2,iepoch+j)
        dvel(j) = gvec(3+i,1,iepoch) - gvec(3+i,2,iepoch+j)

        dist(j)  = dist(j)  + dvec(j) * dvec(j)
        deriv(j) = deriv(j) + dvel(j) * dvec(j)
        angle(j) = dot_angle_calc(dvec,dvec0) * fact
      end do

    end do

    write(lu_angle,'(i6,60f12.3)') iepoch, (angle(j), j=1,60)
    write(lu_deriv,'(i6,60f12.3)') iepoch, (deriv(j), j=1,60)
    write(lu_dist, '(i6,60f12.3)') iepoch, (dist(j),  j=1,60)

  end do
! ******************************************************************

  close(lu_angle)
  close(lu_deriv)
  close(lu_dist)


  stop

999 write(6,'(2(a,/),a,i1,2a,/,a,/,a,i6,/,a,/)') &
      "++++++++++++++++++++++++++++++++++++++++++++", &
      "ERROR: dist_deriv",                            &
      "       error reading file ", isat, ": ", trim(gnv(isat)), &
      "       unexpected end of file no. ",           &
      "       found at line: ", iline, " of file header", &
      "++++++++++++++++++++++++++++++++++++++++++++"
    stop

998 write(6,'(2(a,/),a,i1,2a,/,a,/,a,i6,/,a,/)') &
      "++++++++++++++++++++++++++++++++++++++++++++", &
      "ERROR: dist_deriv",                            &
      "       error reading file ", isat, ": ", trim(gnv(isat)), &
      "       unexpected end of file no. ",           &
      "       found at record: ", irec, " of file",   &
      "++++++++++++++++++++++++++++++++++++++++++++"
    stop

  end


  double precision function dot_angle_calc (v1, v2)

! ************************************************************************************************
!
!     Author: A.Purcell           2:30 PM  FRI., 21  JANUARY, 2022
!
! ************************************************************************************************
! Subroutine to calculate the angle between two vectors (in radians) using the normalised dot product
!
!    v1, v2    : double precision 3 dimensional Cartesian vectors 
! ************************************************************************************************


! **********************************************************
  double precision, dimension(3), intent(in)  :: v1, v2


  double precision :: dot, mag1, mag2


  integer :: i


  logical, parameter :: debug_flag = .false.


  character(len=10) :: prog_name = "dot_angle "
! **********************************************************


! **********************************************************
  if ( debug_flag ) write(6,'(2(a,/),2(a,3f12.5,/))') &
    "***********************************", &
    "FUNCTION dot_angle_calc ",  &
    "v1: ", v1(1), v1(2), v1(3), &
    "v2: ", v2(1), v2(2), v2(3)
! **********************************************************


! **********************************************************
  dot  = 0.0d0
  mag1 = 0.0d0
  mag2 = 0.0d0

  do i = 1, 3
    dot  = dot  + v1(i) * v2(i)
    mag1 = mag1 + v1(i)**2
    mag2 = mag2 + v2(1)**2
  end do

  dot_angle_calc = acos(dot/(sqrt(mag1)*sqrt(mag2)))
! **********************************************************

! **********************************************************
  if ( debug_flag ) write(6,'(a,f12.5,/,2(a,/))') &
    "angle: ", dot_angle_calc,          &
    "EXITING FUNCTION dot_angle_calc ", &
    "***********************************"
! **********************************************************

  return
  end function dot_angle_calc
