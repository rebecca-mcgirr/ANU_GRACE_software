!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  Routines related to inertial-efixed computations
!!
!!  generate_inert_efixed  ::   generates the complete celestial-to-terrestrial 3x3 rotation matrix and its time derivative. 
!!  read_eop_ut1           ::   reads the usno pole/ut1-utc file and returns all the values (either bull_a or bull_b)
!!  iau_interp             ::   interpolates pole/ut1-utc information to epoch required and adds the sub-daily tidal signals
!!  inert_interp           ::   interpolates the celestial-terrestrial rotation and rotation rate matrices to the required epoch
!!
!!
!!  P. Tregoning
!!  17 March 2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine generate_inert_efixed(calling_prog,debug,mjd_start,mjd_end,rot_mats,rotdots,rotaccs,rpy,rpydots,rot_dates,nepochs &
                                   ,dXp,dYp,dUT1)

! subroutine to interact with the SOFA IERS2010 routines to calculate all the precession/nutation/eop
! information (and time derivatives) for the entire interval of the required orbit.
!
! P. Tregoning
! 17 March 2014
!
! PT170327: moved to grace/lib and added dXp, dYp, dUT1 (three daily bias values to the interpolated EOP values)

  use usno_mod    ! provides variables of max and actual number of eop values in usno.finals.data

  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  character(*),intent(in)   :: calling_prog                   ! name of calling program
  logical     ,intent(in)   :: debug                          ! flag to print output statements or not
  real(kind=8),intent(in )  :: mjd_start,mjd_end              ! start and end julian dates for the required orbit integration

  integer*4   ,intent(in)   :: nepochs                        ! number of epochs of celestial-to-terrestrial rotations computed
  real(kind=8),intent(in)   :: dXp,dYp,dUT1                   ! daily bias values to be added to interpolated Xp,Yp,UT1 (in mas)
  real(kind=8),intent(out)  :: rot_mats(3,3,nepochs)          ! vector of rotations about X-axis                  
  real(kind=8),intent(out)  :: rotdots(3,3,nepochs)           ! array of time derivatives of rotation matrices
  real(kind=8),intent(out)  :: rotaccs(3,3,nepochs)            ! array of double time derivatives of rotation matrices (ie the acceleration of the rotation)
  real(kind=8),intent(out)  :: rpy(3,nepochs)                 ! roll/pitch/yaw angles for the rotation matrices
  real(kind=8),intent(out)  :: rpydots(3,nepochs)             ! roll/pitch/yaw angle rates for the rotation matrices
  real(kind=8),intent(out)  :: rot_dates(nepochs)             ! vector of epochs at which the rotation matrices have been computed


! time variables
  real(kind=8)   ::  TT1, TT2                 !  time variables for the SOFA IAU routines

! local variables
  real(kind=8)   ::  tmp_mjd                        ! the temporary counter of the mjd 
  real(kind=8)   ::  xp_int, yp_int,ut1utc_int      ! X-pole, Y-pole, UT1-UTC interpolated to the required epoch (in the loop over the orbit span)
  integer*4      ::  status_flag                    ! something required for an iau_SOFA routine (but not used here)
  integer*4      ::  iepoch                         ! epoch counter
  character      ::  usnofile*16                    ! (hardwired) name of input usno eop/ut1-utc file
  real(kind=8)   ::  pi
  character      :: message*256

! SA add
  real (kind=16):: T =0, ERA=0, F=0
  double precision :: RC2I(3,3)=0, SP=0, RPOM(3,3)=0, RC2T(3,3)=0
  double precision :: iau_SP00
! SA END ADD
  integer*4  :: i,j,iusno
  integer :: file1, file2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pi = 4.d0*datan(1.d0)
  iusno = 40
  usnofile = "usno.finals.data"
  rot_mats = 0.d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first, get the usno pole and ut1-utc information
  call read_eop_ut1(calling_prog,"bull_b",iusno,usnofile)
! store the mjdates within the Xp and Yp arrays (for later use in poltid to correct C21, S21)
  Xp(1:n_usnovals,2) = mjdates(1:n_usnovals)
  Yp(1:n_usnovals,2) = mjdates(1:n_usnovals)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now derive the cellestial-to-terrestrial rotation matrix for each required epoch. Start and end 5 minutes before the orbit times to 
! allow for interpolation.
  tmp_mjd = mjd_start - 5.d0 / (60.d0*24.d0)
  iepoch = 0
! PT140514: stop the loop if we reach the maximum epochs (this fixes a datetimee bug)
  do while ( tmp_mjd < mjd_end + 5.d0 / (60.d0*24.d0) .and. iepoch < nepochs)
    iepoch = iepoch + 1
    rot_dates(iepoch) = tmp_mjd

    if(mod (iepoch,1000) == 0 .and. debug)then 
      write(message,'(a,i6,a,f12.6,a)')  ' Generate rotation matrices for epoch',iepoch,' (mjd ',dble(tmp_mjd),')'
      call status_update('STATUS',calling_prog,'generate_inert_efixed',' ',message,0)
    endif

! interpolate pole/ut1-utc to this epoch
    call iau_INTERP (mjdates,Xp(:,1),Yp(:,1),UT1utc,n_usnovals,tmp_mjd,xp_int,yp_int,ut1utc_int)
!print*,'xp_int,yp_int,ut1utc_int',xp_int,yp_int,ut1utc_int

! PT170327: add the daily bias values to the interpolated values. These may well be zero, but not necessarily.
    xp_int = xp_int + dXp/1.d3
    yp_int = yp_int + dYp/1.d3
    ut1utc_int = ut1utc_int + dUT1/1.d3

! convert the epoch to the required time units for SOFA precession/nutation/eop routine
    call iau_taitt(2400000.5d0,tmp_mjd,TT1,TT2,status_flag)
! generate the complete rotation matrix for this epoch
! SA Change 1/(180/3600.0) with * 20.0
!    call iau_c2t06a(TT1, TT2, 2400000.5d0, tmp_mjd+(ut1utc_int+0.01d0)/dble(86400.0), XP_int*pi/dble(180)/dble(3600.0) & 
!            ,Yp_int*pi/dble(180)/dble(3600.0) &
!            , rot_mats(:,:,iepoch) )
!  SA We "explode" the c2t06a subrouting, below to convert the ERA in quad
!  precision
    CALL iau_C2I06A ( TT1, TT2, RC2I )
    ! T,F, and ERA are in quad precision here
    ! Contrary to the ERA we don't check which of tmp_mjd + ut1utc_in_day or
    ! 2400000.5 is the biggest in  T = smallest - (bigest - DJ00) we will
    ! need to chage the order after the 8th November 8429 .... Assumimg a Grace
    ! mission every 10years, we will need to change it for GRACE#842
    ! BTW, 2400000.5 - 2451545D0 = -51544.5 (If my computation are correct)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! I modified this part to quadruple precision (REAL*16) due to some
    ! unstabilities to compute the normal equation with respect to error in
    ! UT1-UTC
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! quad precision version
! PT170323: add 5 mas error to UT1-UTC just to see the effect on the orbit
!    ut1utc_int = ut1utc_int + 0.005d0
    T=tmp_mjd+ut1utc_int/86400.q0 -51544.5q0
    F=MOD(2400000.5q0,1q0)+ MOD (tmp_mjd+ut1utc_int/86400.q0, 1q0)
    ERA=MOD ( 2.0q0*pi  * ( F + 0.7790572732640q0 + 0.00273781191135448q0 * T ) ,   2.0q0*pi)
! double precision version
!    T=tmp_mjd+ut1utc_int/86400.d0 -51544.5d0
!    F=MOD(2400000.5d0,1d0)+ MOD (tmp_mjd+ut1utc_int/86400.d0, 1d0)
!    ERA=MOD ( 2.0d0*pi  * ( F + 0.7790572732640d0 + 0.00273781191135448d0 * T ) ,   2.0d0*pi)
    SP = iau_SP00 ( TT1, TT2 )
!PT170323: add 2 mas to XP for a test
!    YP_int = YP_int + 0.002d0
    CALL iau_POM00( XP_int*pi/dble(180)/dble(3600.0), Yp_int*pi/dble(180)/dble(3600.0), SP, RPOM) 
    CALL iau_C2TCIO(RC2I, dble(ERA), RPOM, RC2T)
    rot_mats(:,:,iepoch) = RC2T
! DEBUG
!    call printmat(rot_mats(:,:,iepoch),3,3,'rot_mats ')
!    stop 'stopped in generate_inert_efixd'

! PT140321: compute and store the angular velocity vector (from roll/pitch/yaw angles)
    call rotmat2rpy( rot_mats(:,:,iepoch),rpy(:,iepoch) )

    tmp_mjd = tmp_mjd + 5.d0 / 86400.d0
  enddo

! now, generate the time derivatives of the rotation matrices
  call status_update('STATUS',calling_prog,'generate_inert_efixed',' ','generate the time derivatives of rotation matrices',0)
  do i=1,3
    do j=1,3
      call noise_robust_deriv(rot_mats(i,j,:),rotdots(i,j,:),5.d0,nepochs,5)
! PT170821: and again to get the rotation acceleration matrices
      call noise_robust_deriv(rotdots(i,j,:),rotaccs(i,j,:),5.d0,nepochs,5)
    enddo
! and again to get the time derivatives of roll/pitch/yaw
    call noise_robust_deriv( rpy(i,:),rpydots(i,:),5.d0,nepochs,5 )
  enddo  
! SA190702 Comment those line which wrote the Rotation matrices (it was debugging)
!  open (unit=50, file = "RMATdot.txt")
!  do i = 1 , iepoch
!          write (50,*) rotdots(:,:,i)
!  enddo
!  close(50)
!
!  open(unit=51, file = "RMAT.txt")
!  do i = 1 , iepoch
!          write (51,*) rot_mats(:,:,i)
!  enddo
!  close(51)

! close the USNO file
  close(iusno)
  
  return
  end subroutine generate_inert_efixed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine inert_interp(mjd,nepochs,rot_dates,rot_mats,rotdots,rotaccs,rot_i2e,rotdot_i2e &
                         ,rot_e2i,rotdot_e2i,rotacc_e2i,rotacc_i2e ) !,rpy,rpydots,angvel,angvelrates)

! subroutine to interpolate the celestial-efixed rotation and rotation rate information to the required epoch.
! For the rotation, we first interpolate the three angles, then reconstruct the rotation matrix.
! For the rotation rate, this doesn't work, so we interpolate each of the 3x3 matrix elements to the epoch
!
! Routine returns the inertial->efixed rotation and rotation rate matrices and also the efixed->inertial
!
! P. Tregoning
! 18 March 2014
!
! PT170821: added rotaccs, rotacc_e2i, rotacc_i2e

  implicit none

  real(kind=8) , intent(in)  :: mjd                    ! mjd of required epoch
  integer*4    , intent(in)  :: nepochs                ! number of epochs in rotation matrix arrays
  real(kind=8) , intent(in)  :: rot_dates(nepochs)     ! vector of epochs in rotation matrix arrays
  real(kind=8) , intent(in)  :: rot_mats(3,3,nepochs)  ! array of rotation matrices
  real(kind=8) , intent(in)  :: rotdots(3,3,nepochs)   ! array of rotation rate matrices
  real(kind=8) , intent(in)  :: rotaccs(3,3,nepochs)   ! array of rotation acceleration matrices
!  real(kind=8) , intent(in)  :: rpy(3,nepochs)         ! roll/pitch/yaw angles for the rotation matrices
!  real(kind=8) , intent(in)  :: rpydots(3,nepochs)     ! roll/pitch/yaw angle rates for the rotation matrices
  real(kind=8) , intent(out) :: rot_i2e(3,3)           ! inert->efixed rotation matrix for epoch "mjd"
  real(kind=8) , intent(out) :: rotdot_i2e(3,3)        ! inert->efixed rotation rate matrix for epoch "mjd"
  real(kind=8) , intent(out) :: rotacc_i2e(3,3)        ! inert->efixed rotation accel matrix for epoch "mjd"
  real(kind=8) , intent(out) :: rot_e2i(3,3)           ! efixed->inert rotation matrix for epoch "mjd"
  real(kind=8) , intent(out) :: rotdot_e2i(3,3)        ! efixed->inert rotation rate matrix for epoch "mjd"
  real(kind=8) , intent(out) :: rotacc_e2i(3,3)        ! efixed->inert rotation accel matrix for epoch "mjd"
!  real(kind=8) , intent(out) :: angvel(3)              ! roll/pitch/yaw angles interpolated to the epoch
!  real(kind=8) , intent(out) :: angvelrates(3)         ! roll/pitch/yaw angle rates interpolated to the epoch

! local variables
  integer*4  :: i,j
  real(kind=8) :: tmpmat(3,3)

! interpolate each element of the rotation and rotation rate matrices to the required epoch.
! Do this with LAGINT, the lagrangian interpolation routine that is part of IAU_interp
  do i=1,3
    do j=1,3
      call lagint(rot_dates,rot_mats(i,j,:),nepochs,mjd,   rot_i2e(i,j) )
      call lagint(rot_dates, rotdots(i,j,:),nepochs,mjd,rotdot_i2e(i,j) )
      call lagint(rot_dates, rotaccs(i,j,:),nepochs,mjd,rotacc_i2e(i,j) )
    enddo
    enddo

! now, the rotation matrix for the other direction is the inverse. Because it's orthogonal, this is just the transpose
  call transp(rot_i2e,rot_e2i,3,3)
! it seems that we can do the same for the rotation rate matrix (**** check against rotsnp !!!)
  call transp(rotdot_i2e,rotdot_e2i,3,3)
  call transp(rotacc_i2e,rotacc_e2i,3,3)

! c'est tout 
  return
  end subroutine inert_interp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine inert_interp_v2(interp_flag,mjd,nepochs,rot_dates,rot_mats,rotdots,rotaccs,rot_i2e,rotdot_i2e &
                         ,rot_e2i,rotdot_e2i,rotacc_e2i,rotacc_i2e ) !,rpy,rpydots,angvel,angvelrates)

! subroutine to interpolate the celestial-efixed rotation and rotation rate information to the required epoch.
! For the rotation, we first interpolate the three angles, then reconstruct the rotation matrix.
! For the rotation rate, this doesn't work, so we interpolate each of the 3x3 matrix elements to the epoch
!
! Routine returns the inertial->efixed rotation and rotation rate matrices and also the efixed->inertial
!
! P. Tregoning
! 18 March 2014
!
! PT170821: added rotaccs, rotacc_e2i, rotacc_i2e

  implicit none

  integer*4    , intent(in)  :: interp_flag            ! 1: only the rot_i2e and rot_e2i; 2: all rotation matrices 
  real(kind=8) , intent(in)  :: mjd                    ! mjd of required epoch
  integer*4    , intent(in)  :: nepochs                ! number of epochs in rotation matrix arrays
  real(kind=8) , intent(in)  :: rot_dates(nepochs)     ! vector of epochs in rotation matrix arrays
  real(kind=8) , intent(in)  :: rot_mats(3,3,nepochs)  ! array of rotation matrices
  real(kind=8) , intent(in)  :: rotdots(3,3,nepochs)   ! array of rotation rate matrices
  real(kind=8) , intent(in)  :: rotaccs(3,3,nepochs)   ! array of rotation acceleration matrices
!  real(kind=8) , intent(in)  :: rpy(3,nepochs)         ! roll/pitch/yaw angles for the rotation matrices
!  real(kind=8) , intent(in)  :: rpydots(3,nepochs)     ! roll/pitch/yaw angle rates for the rotation matrices
  real(kind=8) , intent(out) :: rot_i2e(3,3)           ! inert->efixed rotation matrix for epoch "mjd"
  real(kind=8) , intent(out) :: rotdot_i2e(3,3)        ! inert->efixed rotation rate matrix for epoch "mjd"
  real(kind=8) , intent(out) :: rotacc_i2e(3,3)        ! inert->efixed rotation accel matrix for epoch "mjd"
  real(kind=8) , intent(out) :: rot_e2i(3,3)           ! efixed->inert rotation matrix for epoch "mjd"
  real(kind=8) , intent(out) :: rotdot_e2i(3,3)        ! efixed->inert rotation rate matrix for epoch "mjd"
  real(kind=8) , intent(out) :: rotacc_e2i(3,3)        ! efixed->inert rotation accel matrix for epoch "mjd"
!  real(kind=8) , intent(out) :: angvel(3)              ! roll/pitch/yaw angles interpolated to the epoch
!  real(kind=8) , intent(out) :: angvelrates(3)         ! roll/pitch/yaw angle rates interpolated to the epoch

! local variables
  integer*4  :: i,j
  real(kind=8) :: tmpmat(3,3)

! interpolate each element of the rotation and rotation rate matrices to the required epoch.
! Do this with LAGINT, the lagrangian interpolation routine that is part of IAU_interp
  do i=1,3
    do j=1,3
      call lagint(rot_dates,rot_mats(i,j,:),nepochs,mjd,   rot_i2e(i,j) )
      if(interp_flag /= 1)call lagint(rot_dates, rotdots(i,j,:),nepochs,mjd,rotdot_i2e(i,j) )
      if(interp_flag /= 1)call lagint(rot_dates, rotaccs(i,j,:),nepochs,mjd,rotacc_i2e(i,j) )
    enddo
    enddo

! now, the rotation matrix for the other direction is the inverse. Because it's orthogonal, this is just the transpose
  call transp(rot_i2e,rot_e2i,3,3)
! it seems that we can do the same for the rotation rate matrix (**** check against rotsnp !!!)
  if(interp_flag /= 1)call transp(rotdot_i2e,rotdot_e2i,3,3)
  if(interp_flag /= 1)call transp(rotacc_i2e,rotacc_e2i,3,3)

! c'est tout 
  return
  end subroutine inert_interp_v2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_eop_ut1(calling_prog,bulletin,iusno,infile)

! subroutine to read the usno.finals.data file and extract the xp, yp, ut1-utc values
!
! P. Tregoning
! 25 February 2014
!
! S. Allgeyer
! 19 April 2014
! => Read the pole positon as a string, check the length of the string to
! determine the EOF
!
! PT190701: add the reading of Bulletin A values (both predicted and estimated). 

  use usno_mod    ! provides max and actual number of eop values in file

  implicit none

  character(*),intent(in)    ::  calling_prog           ! name of calling program
  character, intent(in)      ::  bulletin*6             ! either "bull_a" or "bull_b"
  integer*4, intent(in)      ::  iusno                  ! unit number of input usno.finals.data file
  character*16, intent(in)   ::  infile                 ! name of input file (usno.finals.data)
  real(kind=8) :: tmpdPSI, tmpdEPSILON
! local variables
  integer*4    :: ioerr  
  character    :: message*200
  integer*4    :: iyr, imon,iday
  real(kind=8) :: sigx,sigy,sigut1
  character*31 :: readval

! PT190701: variables used when reading both Bulletin A and Bulletin B
  character*185 :: line_usnofile               ! read the whole line of the usno.finals.data file (to be decoded subsequently)
  real(kind=8)  :: Xp_a,Yp_a,ut1utc_a          ! Bulletin A values read from a single line
  real(kind=8)  :: Xp_b,Yp_b,ut1utc_b          ! Bulletin B values read from a single line
  
  
! open the file
  open(iusno,file=infile,status='old',iostat=ioerr)
  if(ioerr /= 0)then
    call status_update('FATAL',calling_prog,'read_eop_ut1',infile,'Error opening eop/ut1 file',0)
  endif

  ioerr = 0
  n_usnovals = 1
! loop through and read to the end of the file, storing the values as we go
! PT190701: modify this to read bull_a and bull_b. Store and return whichever one has been requested.
  do while (ioerr == 0 )
    line_usnofile = " "
    read(iusno,'(a)')line_usnofile
!print*,line_usnofile
!    read(iusno,101,iostat=ioerr,end=1000)mjdates(n_usnovals),xp_a,yp_a,ut1utc_a,xp_b,yp,ut1utc_b
!101 format(6x,f9.2,2x,2(1x,f9.6,9x,f9.6),3x,f10.7,66x,2f10.6,f11.7)
    read(line_usnofile,101)mjdates(n_usnovals),xp_a,yp_a ,ut1utc_a,xp_b,yp_b,ut1utc_b
101 format(6x,f9.2,3x,f9.6,10x,f9.6,12x,f10.7,66x,2f10.6,f11.7)
!print*,'after reading the line ioerr=',mjdates(n_usnovals),xp_a,yp_a ,ut1utc_a,xp_b,yp_b,ut1utc_b

! check whether there are bulletin B values available
    if(line_usnofile(138:144) /= "       ")then  
    ! there are interpolated bulletin B values. Transfer them to the output array
      xp(n_usnovals,1) = xp_b
      yp(n_usnovals,1) = yp_b
      ut1utc(n_usnovals) = ut1utc_b
      eop_type(n_usnovals) = "B"
! if not, what about bulletin A interpolated?
    else if (line_usnofile(17:17) == "I")then
      xp(n_usnovals,1) = xp_a
      yp(n_usnovals,1) = yp_a
      ut1utc(n_usnovals) = ut1utc_a
      eop_type(n_usnovals) = "A"
! if not, what about Bulletin A predicted?
    else  if (line_usnofile(17:17) == "P")then
!print*,'predicted Bulletin A values',mjdates(n_usnovals)
      xp(n_usnovals,1) = xp_a
      yp(n_usnovals,1) = yp_a
      ut1utc(n_usnovals) = ut1utc_a
      eop_type(n_usnovals) = "P"
    else if (line_usnofile(17:17) == " ")then ! no values read. Stop reading the file
!print*,'run out of Bulletin A values',mjdates(n_usnovals)
      ioerr = -1
    endif
!print*,xp(n_usnovals,1),yp(n_usnovals,1),ut1utc(n_usnovals)

    if(ioerr == 0 )then
      n_usnovals = n_usnovals + 1
    else  
      n_usnovals = n_usnovals - 1
      ioerr = -1          ! sets the flag to end the do while loop. We have reached the end of Bulletin B values.
    endif
    if(n_usnovals > max_usnovals) then
      call status_update('FATAL',calling_prog,'inert_lib/read_eop_ut1',infile,'EOP/UT1 file has too many records',0)
    endif
  enddo

1000 continue
  write(message,'(a,i6,a)')'Have read ',n_usnovals,' EOP/UT1-UTC values from file '
  call status_update('STATUS',calling_prog,'inert_lib/read_eop_ut1',infile,message,0)

  return
  end subroutine read_eop_ut1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine adjust_R(rot_matrix,rpy_angles)

! subroutine to add small angles to the roll/pitch/yaw angles of a rotation matrix, then recompute and replace
! the rotation matrix
!
! P. Tregoning
! 28 March 2014

  use rpy_mod

  implicit none

  integer*4 i,j
  real(kind=8), intent(inout) :: rot_matrix(3,3)    ! rotation matrix to be adjusted
  real(kind=8), intent(inout) :: rpy_angles(3)      ! roll/pitch/yaw angles by which it should be adjusted

  real(kind=8),dimension(3,3) :: Rx, Ry, Rz, Rtmp1, Rtmp2
  real(kind=8) :: pi

  pi = 4.d0*datan(1.d0)

! add the adjustment angles (passed in through rpy_mod) to the rpy angles for this rotation matrix
  rpy_angles = rpy_angles + rot_ang_errors * pi/180.d0  / 3600.d3   ! conversion from milliarcseconds to radians

! convert into a rotation matrix
  call rotmat(-rpy_angles(1),1,Rx)
  call rotmat(-rpy_angles(2),2,Ry)
  call rotmat(-rpy_angles(3),3,Rz)

  call matmult(Ry,Rx,Rtmp1,3,3,3)
  call matmult(Rz,Rtmp1,Rtmp2,3,3,3)

! DEBUG
!  print*,'adjustment angles:',rot_ang_errors, rot_ang_errors * pi/180.d0  / 3600.d3
!  call printmat(Rx,3,3,'Rx ')
!  call printmat(Ry,3,3,'Ry ')
!  call printmat(Rz,3,3,'Rz ')
!
!  call printmat (rot_matrix,3,3,'rotmat_orig ')
!  call printmat(Rtmp2,3,3,'rotmat_adjusted ')
!
!  stop


  return
  end subroutine adjust_R


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------------------------------------
      SUBROUTINE iau_INTERP (RJD,X,Y,UT1,N,rjd_int,x_int,y_int,ut1_int)
!
!      This subroutine takes a series of x, y, and UT1-UTC values
!      and interpolates them to an epoch of choice. This routine
!      assumes that the values of x and y are in seconds of
!      arc and that UT1-UTC is in seconds of time. At least
!      one point before and one point after the epoch of the
!      interpolation point are necessary in order for the
!      interpolation scheme to work. 
!
!      parameters are :
!      RJD     - array of the epochs of data (given in mjd)
!      X       - array of x polar motion (arcsec)
!      Y       - array of y polar motion (arcsec)
!      UT1     - array of UT1-UTC (sec)
!      n       - number of points in arrays
!      rjd_int - epoch for the interpolated value
!      x_int   - interpolated value of x
!      y_int   - interpolated value of y
!      ut1_int - interpolated value of ut1-utc
!
!      CALLED SUBROUTINE : LAGINT (Lagrange interpolation)
!                          PMUT1_OCEANS (Diurnal and semidiurnal oceanic effects)
!                          PM_GRAVI (Diurnal and semidiurnal lunisolar effects)
!
!       coded by Ch. BIZOUARD (Observatoire de Paris) : November 2002
!                                           Corrected : September 2007   
  
! SA modify declaration for f90
      implicit none
      integer i,n
      double precision  RJD(N), X(N), Y(N), UT1(N), rjd_int,  x_int, y_int, ut1_int,  cor_x, cor_y, cor_ut1, cor_lod
      

      CALL LAGINT (RJD,X,n,rjd_int,x_int)
  
      CALL LAGINT (RJD,Y,n,rjd_int,y_int)
      
      CALL LAGINT (RJD,UT1,n,rjd_int,ut1_int)
!  --------------
!  Oceanic effect      
!  --------------
! PT200707: renamed to avoid a clash of names with ~/gg/libraries/comlib   
      CALL PMUT1_OCEANS_anu (rjd_int,cor_x,cor_y,cor_ut1,cor_lod)

      x_int = x_int + cor_x
      y_int = y_int + cor_y
      ut1_int = ut1_int + cor_ut1

!  Lunisolar effect.
!  PT200707: rename to avoid a clash with ~/gg/libraries/comlib
      CALL PM_GRAVI_anu (rjd_int,cor_x,cor_y)
      
      x_int   = x_int + cor_x
      y_int   = y_int + cor_y
     
      RETURN

      END
!
!----------------------------------------------------------------
! ***************************************************************
      SUBROUTINE LAGINT (X,Y,n,xint,yout)
!  
!      This subroutine performs lagrangian interpolation
!      within a set of (X,Y) pairs to give the y
!      value corresponding to xint. This program uses a
!      window of 4 data points to perform the interpolation.
!      if the window size needs to be changed, this can be
!      done by changing the indices in the do loops for
!      variables m and j.
!
!      PARAMETERS ARE :
!      X     - array of values of the independent variable
!      Y     - array of function values corresponding to x
!      n     - number of points
!      xint  - the x-value for which estimate of y is desired
!      yout  - the y value returned to caller

      implicit none
      REAL*8 X(n),Y(n),xint,yout,term
      INTEGER i,j,k,m,n

      yout = 0.d0
      do  i = 1,n-1
        if ( xint .ge. X(i) .and. xint .lt. X(i+1) ) k = i
      enddo
    
      if ( k .lt. 2 ) k = 2
      if ( k .gt. n-2 ) k = n-2
!      if(k < 6800)print*,'lagint: k,X(k),X(k+1),Y(k),Y(k+1),xint',k,X(k),X(k+1),Y(k),Y(k+1),xint   
      do m = k-1,k+2
        term = y(m)
        do  j = k-1,k+2
          if ( m .ne. j ) then
            term = term * (xint - X(j))/(X(m) - X(j))
          end if
        enddo 
        yout = yout + term
      enddo
      
      return
      end
! ***************************************************************

!----------------------------------------------------------------
      SUBROUTINE PMUT1_OCEANS_anu (rjd,cor_x,cor_y,cor_ut1,cor_lod)
!
!     This subroutine provides, in time domain, the diurnal/subdiurnal
!     tidal effets on polar motion ("), UT1 (s) and LOD (s). The tidal terms,
!     listed in the program above, have been extracted from the procedure   
!     ortho_eop.f coed by Eanes in 1997.
!     
!     N.B.:  The fundamental lunisolar arguments are those of Simon et al.  
!
!     These corrections should be added to "average"
!     EOP values to get estimates of the instantaneous values.
!
!      PARAMETERS ARE :
!      rjd      - epoch of interest given in mjd
!      cor_x    - tidal correction in x (sec. of arc)
!      cor_y    - tidal correction in y (sec. of arc)
!      cor_ut1  - tidal correction in UT1-UTC (sec. of time)
!      cor_lod  - tidal correction in length of day (sec. of time)
!
!      coded by Ch. Bizouard (2002), initially coded by McCarthy and 
!      D.Gambis(1997) for the 8 prominent tidal waves.  
      
      IMPLICIT NONE
      
      INTEGER nlines
      PARAMETER(nlines=71)
      DOUBLE PRECISION ARG(6)     ! Array of the tidal arguments   
      DOUBLE PRECISION DARG(6)    ! Array of their time derivative 
      
      REAL*4 XCOS(nlines),XSIN(nlines),YCOS(nlines),YSIN(nlines),UTCOS(nlines),UTSIN(nlines)
     
      REAL*8 t,ag,dag,rjd,halfpi,secrad,cor_x,cor_y,cor_ut1,cor_lod
      INTEGER NARG(nlines,6),i,j
      
      halfpi = 1.5707963267948966d0
      secrad=2.d0*halfpi/(180.d0*3600.d0)	

!   Oceani!  tidal terms present in x (microas),y(microas),ut1(microseconds)       
!   NARG(j,6) : Multipliers of GMST+pi and Delaunay arguments. 
	
       data( &
       NARG(j,1),NARG(j,2),NARG(j,3),NARG(j,4),NARG(j,5),NARG(j,6),      &
       XSIN(j),XCOS(j),YSIN(j),YCOS(j),UTSIN(j),UTCOS(j),j=1,nlines)/    &
      1,-1, 0,-2,-2,-2,  -0.05,   0.94,  -0.94,  -0.05,  0.396, -0.078,  &
      1,-2, 0,-2, 0,-1,   0.06,   0.64,  -0.64,   0.06,  0.195, -0.059,  &
      1,-2, 0,-2, 0,-2,   0.30,   3.42,  -3.42,   0.30,  1.034, -0.314,  &
      1, 0, 0,-2,-2,-1,   0.08,   0.78,  -0.78,   0.08,  0.224, -0.073,  &
      1, 0, 0,-2,-2,-2,   0.46,   4.15,  -4.15,   0.45,  1.187, -0.387,  &
      1,-1, 0,-2, 0,-1,   1.19,   4.96,  -4.96,   1.19,  0.966, -0.474,  &
      1,-1, 0,-2, 0,-2,   6.24,  26.31, -26.31,   6.23,  5.118, -2.499,  &
      1, 1, 0,-2,-2,-1,   0.24,   0.94,  -0.94,   0.24,  0.172, -0.090,  &
      1, 1, 0,-2,-2,-2,   1.28,   4.99,  -4.99,   1.28,  0.911, -0.475,  &
      1, 0, 0,-2, 0, 0,  -0.28,  -0.77,   0.77,  -0.28, -0.093,  0.070,  &
      1, 0, 0,-2, 0,-1,   9.22,  25.06, -25.06,   9.22,  3.025, -2.280,  &
      1, 0, 0,-2, 0,-2,  48.82, 132.91,-132.90,  48.82, 16.020,-12.069,  &
      1,-2, 0, 0, 0, 0,  -0.32,  -0.86,   0.86,  -0.32, -0.103,  0.078,  &
      1, 0, 0, 0,-2, 0,  -0.66,  -1.72,   1.72,  -0.66, -0.194,  0.154,  &
      1,-1, 0,-2, 2,-2,  -0.42,  -0.92,   0.92,  -0.42, -0.083,  0.074,  &
      1, 1, 0,-2, 0,-1,  -0.30,  -0.64,   0.64,  -0.30, -0.057,  0.050,  &
      1, 1, 0,-2, 0,-2,  -1.61,  -3.46,   3.46,  -1.61, -0.308,  0.271,  &
      1,-1, 0, 0, 0, 0,  -4.48,  -9.61,   9.61,  -4.48, -0.856,  0.751,  &
      1,-1, 0, 0, 0,-1,  -0.90,  -1.93,   1.93,  -0.90, -0.172,  0.151,  &
      1, 1, 0, 0,-2, 0,  -0.86,  -1.81,   1.81,  -0.86, -0.161,  0.137,  &
      1, 0,-1,-2, 2,-2,   1.54,   3.03,  -3.03,   1.54,  0.315, -0.189,  &
      1, 0, 0,-2, 2,-1,  -0.29,  -0.58,   0.58,  -0.29, -0.062,  0.035,  &
      1, 0, 0,-2, 2,-2,  26.13,  51.25, -51.25,  26.13,  5.512, -3.095,  &
      1, 0, 1,-2, 2,-2,  -0.22,  -0.42,   0.42,  -0.22, -0.047,  0.025,  &
      1, 0,-1, 0, 0, 0,  -0.61,  -1.20,   1.20,  -0.61, -0.134,  0.070,  &
      1, 0, 0, 0, 0, 1,   1.54,   3.00,  -3.00,   1.54,  0.348, -0.171,  &
      1, 0, 0, 0, 0, 0, -77.48,-151.74, 151.74, -77.48,-17.620,  8.548,  &
      1, 0, 0, 0, 0,-1, -10.52, -20.56,  20.56, -10.52, -2.392,  1.159,  &
      1, 0, 0, 0, 0,-2,   0.23,   0.44,  -0.44,   0.23,  0.052, -0.025,  &
      1, 0, 1, 0, 0, 0,  -0.61,  -1.19,   1.19,  -0.61, -0.144,  0.065,  &
      1, 0, 0, 2,-2, 2,  -1.09,  -2.11,   2.11,  -1.09, -0.267,  0.111,  &
      1,-1, 0, 0, 2, 0,  -0.69,  -1.43,   1.43,  -0.69, -0.288,  0.043,  &
      1, 1, 0, 0, 0, 0,  -3.46,  -7.28,   7.28,  -3.46, -1.610,  0.187,  &
      1, 1, 0, 0, 0,-1,  -0.69,  -1.44,   1.44,  -0.69, -0.320,  0.037,  &
      1, 0, 0, 0, 2, 0,  -0.37,  -1.06,   1.06,  -0.37, -0.407, -0.005,  &
      1, 2, 0, 0, 0, 0,  -0.17,  -0.51,   0.51,  -0.17, -0.213, -0.005,  &
      1, 0, 0, 2, 0, 2,  -1.10,  -3.42,   3.42,  -1.09, -1.436, -0.037,  &
      1, 0, 0, 2, 0, 1,  -0.70,  -2.19,   2.19,  -0.70, -0.921, -0.023,  &
      1, 0, 0, 2, 0, 0,  -0.15,  -0.46,   0.46,  -0.15, -0.193, -0.005,  &
      1, 1, 0, 2, 0, 2,  -0.03,  -0.59,   0.59,  -0.03, -0.396, -0.024,  &
      1, 1, 0, 2, 0, 1,  -0.02,  -0.38,   0.38,  -0.02, -0.253, -0.015,  &
      2,-3, 0,-2, 0,-2,  -0.49,  -0.04,   0.63,   0.24, -0.089, -0.011,  &
      2,-1, 0,-2,-2,-2,  -1.33,  -0.17,   1.53,   0.68, -0.224, -0.032,  &
      2,-2, 0,-2, 0,-2,  -6.08,  -1.61,   3.13,   3.35, -0.637, -0.177,  &
      2, 0, 0,-2,-2,-2,  -7.59,  -2.05,   3.44,   4.23, -0.745, -0.222,  &
      2, 0, 1,-2,-2,-2,  -0.52,  -0.14,   0.22,   0.29, -0.049, -0.015,  &
      2,-1,-1,-2, 0,-2,   0.47,   0.11,  -0.10,  -0.27,  0.033,  0.013,  &
      2,-1, 0,-2, 0,-1,   2.12,   0.49,  -0.41,  -1.23,  0.141,  0.058,  &
      2,-1, 0,-2, 0,-2, -56.87, -12.93,  11.15,  32.88, -3.795, -1.556,  &
      2,-1, 1,-2, 0,-2,  -0.54,  -0.12,   0.10,   0.31, -0.035, -0.015,  &
      2, 1, 0,-2,-2,-2, -11.01,  -2.40,   1.89,   6.41, -0.698, -0.298,  &
      2, 1, 1,-2,-2,-2,  -0.51,  -0.11,   0.08,   0.30, -0.032, -0.014,  &
      2,-2, 0,-2, 2,-2,   0.98,   0.11,  -0.11,  -0.58,  0.050,  0.022,  &
      2, 0,-1,-2, 0,-2,   1.13,   0.11,  -0.13,  -0.67,  0.056,  0.025,  &
      2, 0, 0,-2, 0,-1,  12.32,   1.00,  -1.41,  -7.31,  0.605,  0.266,  &
      2, 0, 0,-2, 0,-2,-330.15, -26.96,  37.58, 195.92,-16.195, -7.140,  &
      2, 0, 1,-2, 0,-2,  -1.01,  -0.07,   0.11,   0.60, -0.049, -0.021,  &
      2,-1, 0,-2, 2,-2,   2.47,  -0.28,  -0.44,  -1.48,  0.111,  0.034,  &
      2, 1, 0,-2, 0,-2,   9.40,  -1.44,  -1.88,  -5.65,  0.425,  0.117,  &
      2,-1, 0, 0, 0, 0,  -2.35,   0.37,   0.47,   1.41, -0.106, -0.029,  &
      2,-1, 0, 0, 0,-1,  -1.04,   0.17,   0.21,   0.62, -0.047, -0.013,  &
      2, 0,-1,-2, 2,-2,  -8.51,   3.50,   3.29,   5.11, -0.437, -0.019,  &
      2, 0, 0,-2, 2,-2,-144.13,  63.56,  59.23,  86.56, -7.547, -0.159,  &
      2, 0, 1,-2, 2,-2,   1.19,  -0.56,  -0.52,  -0.72,  0.064,  0.000,  &
      2, 0, 0, 0, 0, 1,   0.49,  -0.25,  -0.23,  -0.29,  0.027, -0.001,  &
      2, 0, 0, 0, 0, 0, -38.48,  19.14,  17.72,  23.11, -2.104,  0.041,  &
      2, 0, 0, 0, 0,-1, -11.44,   5.75,   5.32,   6.87, -0.627,  0.015,  &
      2, 0, 0, 0, 0,-2,  -1.24,   0.63,   0.58,   0.75, -0.068,  0.002,  &
      2, 1, 0, 0, 0, 0,  -1.77,   1.79,   1.71,   1.04, -0.146,  0.037,  &
      2, 1, 0, 0, 0,-1,  -0.77,   0.78,   0.75,   0.45, -0.064,  0.017,  &
      2, 0, 0, 2, 0, 2,  -0.33,   0.62,   0.65,   0.19, -0.049,  0.018/   
  
      T = (rjd - 51544.5D0)/36525.0D0  ! julian century

!  Arguments in the following order : chi=GMST+pi,l,lp,F,D,Omega
!  et leur derivee temporelle 

      ARG(1) = (67310.54841d0 + (876600d0*3600d0 + 8640184.812866d0)*T + 0.093104d0*T**2 - 6.2d-6*T**3)*15.0d0 + 648000.0d0
      ARG(1)=dmod(ARG(1),1296000d0)*secrad 
   
      DARG(1) = (876600d0*3600d0 + 8640184.812866d0 + 2.d0 * 0.093104d0 * T - 3.d0 * 6.2d-6*T**2)*15.d0
      DARG(1) = DARG(1)* secrad / 36525.0D0   ! rad/day


      ARG(2) = -0.00024470d0*T**4 + 0.051635d0*T**3 + 31.8792d0*T**2 + 1717915923.2178d0*T + 485868.249036d0
      ARG(2) = DMOD(ARG(2),1296000d0)*secrad
      
      DARG(2) = -4.d0*0.00024470d0*T**3 + 3.d0*0.051635d0*T**2 + 2.d0*31.8792d0*T + 1717915923.2178d0 
      DARG(2) = DARG(2)* secrad / 36525.0D0   ! rad/day

      ARG(3) = -0.00001149d0*T**4 - 0.000136d0*T**3 -  0.5532d0*T**2 + 129596581.0481d0*T + 1287104.79305d0
      ARG(3) = DMOD(ARG(3),1296000d0)*secrad

      DARG(3) = -4.D0*0.00001149d0*T**3 - 3.d0*0.000136d0*T**2 -  2.D0*0.5532d0*T + 129596581.0481d0
      DARG(3) = DARG(3)* secrad / 36525.0D0   ! rad/day
          
      ARG(4) = 0.00000417d0*T**4 - 0.001037d0*T**3 - 12.7512d0*T**2 + 1739527262.8478d0*T + 335779.526232d0
      ARG(4) = DMOD(ARG(4),1296000d0)*secrad

      DARG(4) = 4.d0*0.00000417d0*T**3 - 3.d0*0.001037d0*T**2 - 2.d0 * 12.7512d0*T + 1739527262.8478d0 
      DARG(4) = DARG(4)* secrad / 36525.0D0   ! rad/day
    
      ARG(5) = -0.00003169d0*T**4 + 0.006593d0*T**3 - 6.3706d0*T**2 + 1602961601.2090d0*T + 1072260.70369d0
      ARG(5) = DMOD(ARG(5),1296000d0)*secrad

      DARG(5) = -4.d0*0.00003169d0*T**3 + 3.d0*0.006593d0*T**2 - 2.d0 * 6.3706d0*T + 1602961601.2090d0
      DARG(5) = DARG(5)* secrad / 36525.0D0   ! rad/day

      ARG(6) = -0.00005939d0*T**4 + 0.007702d0*T**3 + 7.4722d0*T**2 - 6962890.2665d0*T + 450160.398036d0
      ARG(6) = DMOD(ARG(6),1296000d0)*secrad

      DARG(6) = -4.d0*0.00005939d0*T**3 + 3.d0 * 0.007702d0*T**2 + 2.d0 * 7.4722d0*T - 6962890.2665d0
      DARG(6) = DARG(6)* secrad / 36525.0D0   ! rad/day

!  CORRECTIONS

	cor_x  = 0.d0
	cor_y  = 0.d0
	cor_ut1= 0.d0
	cor_lod= 0.d0

 	do j=1,nlines
 	
	ag  = 0.d0
 	dag = 0.d0
		do i=1,6
 		ag  = ag  + dble(narg(j,i))*ARG(i)
 		dag = dag + dble(narg(j,i))*DARG(i)
		enddo
	ag=dmod(ag,4.d0*halfpi)

        cor_x =cor_x + dble(XCOS(j))*dcos(ag) + dble(XSIN(j)) *dsin(ag)
        cor_y =cor_y + dble(YCOS(j))*dcos(ag) + dble(YSIN(j)) *dsin(ag)
        cor_ut1=cor_ut1+dble(UTCOS(j))*dcos(ag)+dble(UTSIN(j))*dsin(ag)
        cor_lod=cor_lod -(-dble(UTCOS(j)) * dsin(ag) &
                          + dble(UTSIN(j)) * dcos(ag) ) * dag   	 

        enddo
  
       cor_x   = cor_x * 1.0d-6   ! arcsecond (")
       cor_y   = cor_y * 1.0d-6   ! arcsecond (")
       cor_ut1 = cor_ut1 * 1.0d-6 ! second (s)
       cor_lod = cor_lod * 1.0d-6 ! second (s)
 
      RETURN
      END
      	
!----------------------------------------------------------------
      SUBROUTINE PM_GRAVI_anu (rjd,cor_x,cor_y)
!
!     This subroutine provides, in time domain, the diurnal
!     lunisolar effet on polar motion (")
!     
!     N.B.:  The fundamental lunisolar arguments are those of Simon et al.  
!
!     These corrections should be added to "average"
!     EOP values to get estimates of the instantaneous values.
!
!      PARAMETERS ARE :
!      rjd      - epoch of interest given in mjd
!      cor_x    - tidal correction in x (sec. of arc)
!      cor_y    - tidal correction in y (sec. of arc)
!
!      coded by Ch. Bizouard (2002)
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      INTEGER nlines
      PARAMETER(nlines=10)
      DOUBLE PRECISION ARG(6)    ! Array of the tidal arguments   
      REAL*4 XCOS(nlines),XSIN(nlines),YCOS(nlines),YSIN(nlines)
      INTEGER NARG(nlines,6)
      
      halfpi = 1.5707963267948966d0
      secrad=2.d0*halfpi/(180.d0*3600.d0)	

!   Diurnal lunisolar tidal terms present in x (microas),y(microas)      
!   NARG(j,6) : Multipliers of GMST+pi and Delaunay arguments. 
	
       data( NARG(j,1),NARG(j,2),NARG(j,3),NARG(j,4),NARG(j,5),NARG(j,6), &
      XSIN(j),XCOS(j),YSIN(j),YCOS(j),j=1,nlines)/     &
      1,-1, 0,-2, 0,-1,    -.44,   .25,   -.25,  -.44, &
      1,-1, 0,-2, 0,-2,   -2.31,  1.32,  -1.32, -2.31, &
      1, 1, 0,-2,-2,-2,    -.44,   .25,   -.25,  -.44, &
      1, 0, 0,-2, 0,-1,   -2.14,  1.23,  -1.23, -2.14, &
      1, 0, 0,-2, 0,-2,  -11.36,  6.52,  -6.52,-11.36, &
      1,-1, 0, 0, 0, 0,     .84,  -.48,    .48,   .84, &
      1, 0, 0,-2, 2,-2,   -4.76,  2.73,  -2.73, -4.76, &
      1, 0, 0, 0, 0, 0,   14.27, -8.19,   8.19, 14.27, &
      1, 0, 0, 0, 0,-1,    1.93, -1.11,   1.11,  1.93, &
      1, 1, 0, 0, 0, 0,     .76,  -.43,    .43,   .76/
 
      T = (rjd - 51544.5D0)/36525.0D0  ! julian century

!  Arguments in the following order : chi=GMST+pi,l,lp,F,D,Omega
!  et leur derivee temporelle 

      ARG(1) = (67310.54841d0 + (876600d0*3600d0 + 8640184.812866d0)*T + 0.093104d0*T**2 - 6.2d-6*T**3)*15.0d0 + 648000.0d0
      ARG(1)=dmod(ARG(1),1296000d0)*secrad 
   

      ARG(2) = -0.00024470d0*T**4 + 0.051635d0*T**3 + 31.8792d0*T**2 + 1717915923.2178d0*T + 485868.249036d0
      ARG(2) = DMOD(ARG(2),1296000d0)*secrad
      

      ARG(3) = -0.00001149d0*T**4 - 0.000136d0*T**3 -  0.5532d0*T**2 + 129596581.0481d0*T + 1287104.79305d0
      ARG(3) = DMOD(ARG(3),1296000d0)*secrad

          
      ARG(4) = 0.00000417d0*T**4 - 0.001037d0*T**3 - 12.7512d0*T**2 + 1739527262.8478d0*T + 335779.526232d0
      ARG(4) = DMOD(ARG(4),1296000d0)*secrad

    
      ARG(5) = -0.00003169d0*T**4 + 0.006593d0*T**3 - 6.3706d0*T**2 + 1602961601.2090d0*T + 1072260.70369d0
      ARG(5) = DMOD(ARG(5),1296000d0)*secrad

  
      ARG(6) = -0.00005939d0*T**4 + 0.007702d0*T**3 + 7.4722d0*T**2 - 6962890.2665d0*T + 450160.398036d0
      ARG(6) = DMOD(ARG(6),1296000d0)*secrad


!  CORRECTIONS

	cor_x  = 0.d0
	cor_y  = 0.d0

 	do j=1,nlines
 	
	ag  = 0.d0
		do i=1,6
 		ag  = ag  + dble(narg(j,i))*ARG(i)
		enddo
	ag=dmod(ag,4.d0*halfpi)

        cor_x =cor_x+dble(XCOS(j))*dcos(ag)+dble(XSIN(j))*dsin(ag)
        cor_y =cor_y+dble(YCOS(j))*dcos(ag)+dble(YSIN(j))*dsin(ag) 

        enddo
  
      cor_x = cor_x * 1.0d-6   ! arcsecond (")
      cor_y = cor_y * 1.0d-6   ! arcsecond (")
 
      RETURN

      END
