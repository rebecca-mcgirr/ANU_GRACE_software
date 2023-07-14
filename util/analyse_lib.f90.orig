!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!  calc_quadratic_annual     :  computes a 5-parameter LS fit of quadratic+annual
!!!  read_fitfile_to_msc       :  determines number of mascons in a fit file, then leaves it at the line of the first mascon
!!!  read_fitfile_msc_vals     :  reads the fit file EWH values
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calc_quadratic_annual(n_obs,obs,quad_model)

! subroutine to perform a least squares inversion to fit a quadratic+annual to a time series
!
! P. Tregoning
! 17 January 2020

  implicit none

! passed variables
  integer*4    :: n_obs          ! number of observations in the time series  
  real(kind=8) :: obs(n_obs,3)   ! the time series of observations (value, sigma, epoch)
  real(kind=8) :: quad_model(5)  ! offset, rate, ampl_sin, ampl_cos, accel

! local variables
  integer*4                 :: nparam = 5
  real(kind=8), allocatable :: A(:,:),At(:,:),W(:,:),AtW(:,:),AtWA(:,:),B(:,:), AtWb(:,:), VCV(:,:),soln(:,:)
  integer*4                 :: iobs, iparam,iter
  real(kind=8)              :: dt,pi
  logical                   :: debug

  pi = 4.d0*datan(1.d0)
  debug = .false.

! allocate the LS arrays
  allocate(A(n_obs,nparam))
  allocate(W(n_obs,n_obs))
  allocate(At(nparam,n_obs))
  allocate(AtW(nparam,n_obs))
  allocate(AtWA(nparam,nparam))
  allocate(VCV(nparam,nparam))
  allocate(B(n_obs,1))
  allocate(AtWB(nparam,1))
  allocate(soln(nparam,1))

! define a priori values
  quad_model = 0.d0

! iterate the LS (to validate the code. Shouldn't be needed afterwards)
  do iter=1,1
    W = 0.d0
    A = 0.d0
    B = 0.d0

! loop through the observations
    do iobs = 1,n_obs

      dt = obs(iobs,3) - obs(1,3)
      A(iobs,1) = 1.d0                    ! offset parameter
      A(iobs,2) = dt                      ! rate parameter
      A(iobs,3) = dsin(2.d0*pi*dt)        ! sin amplitude
      A(iobs,4) = dcos(2.d0*pi*dt)        ! cosine amplitude
      A(iobs,5) = dt**2                   ! quadratic parameter

      W(iobs,iobs) = 1.d0/obs(iobs,2)**2  ! set the observation weights based on the uncertainty of the observation

      B(iobs,1) = obs(iobs,1)-(quad_model(1) + quad_model(2)*dt + quad_model(3)*dsin(2.d0*pi*dt) &
                                   + quad_model(4)*dcos(2.d0*pi*dt)+quad_model(5)*dt**2)
    enddo

    ! LS solution
    call transp(A,At,n_obs,nparam)
    call matmult(At,W,AtW,nparam,n_obs,n_obs)
    call matmult(AtW,A,AtWA,nparam,n_obs,nparam)
    call invert(AtWA,VCV, nparam)
    call matmult(AtW,B,AtWB,nparam,n_obs,1)
    call matmult(VCV,AtWB,soln,nparam,nparam,1)

    if(debug)print*,"iter ",iter," soln:", soln
    ! update the a priori
    quad_model = quad_model + soln(:,1)

  enddo

! deallocate arrays
  deallocate(A)
  deallocate(W)
  deallocate(At)
  deallocate(AtW)
  deallocate(AtWA)
  deallocate(VCV)
  deallocate(B)
  deallocate(AtWB)
  deallocate(soln)

  end subroutine calc_quadratic_annual
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_fitfile_to_msc(lu_addnorm,n_msc,epoch,flag)

! subroutine to read the header of the solution file, leaving the file at the first line of 
! the mascon values. Flag permits to also determine how many mascons there are in the file (flag=1)
! or not (flag = 2)
!
! P. Tregoning
! 15 January 2020

  implicit none

  integer*4    :: lu_addnorm,n_msc,flag
  real(kind=8) :: epoch                     ! decimal year of the first file found in the addnorm fit file

! local variables
  character*34 :: line
  integer*4    :: ioerr,n_lines_hdr,ihdr
  integer*4    :: ymd(3),doy


! read through to the first mascon line
  n_lines_hdr = 0
  line = " "
  ymd = 0
  do while (line(7:9) /= " MC")
    read(lu_addnorm,'(a)',iostat=ioerr)line

    ! get the yr/mo/day for the first file in the addnorm solution. This will be close enough to use as the time variable
    ! for fitting a quadratic+annual to the EWH time series
    if (line(1:9) == "File:   1")then
      read(line(19:34),*)ymd(1:3)
      call ymd_to_doy(ymd,doy)
      epoch = dble(ymd(1)) + dble(doy)/365.d0
    endif

    if(line(7:9) /= " MC")n_lines_hdr = n_lines_hdr + 1
  enddo

! determine how many mascons there are, if required
  if(flag == 1)then
    n_msc = 1
    do while (line(7:9) == " MC")
      read(lu_addnorm,'(a)',iostat=ioerr,end=1000)line
      if(ioerr == 0 .and. line(7:9) == " MC")n_msc = n_msc + 1      
    enddo
1000 continue

! rewind the file, then step back through to the line of the first mascon
    rewind(lu_addnorm)
    do ihdr = 1,n_lines_hdr
      read(lu_addnorm,'(a)')line
    enddo

  else
    backspace(lu_addnorm)

  endif

  return
  end subroutine read_fitfile_to_msc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_fitfile_msc_vals(lu_addnorm,n_msc,msc_vals)

! subroutine to read the mascon estimates
! 
! P. Tregoning
! 15 January 2020

  implicit none

! passed variables
  integer*4    :: lu_addnorm,n_msc,n_files
  real(kind=8) :: msc_vals(n_msc,3)

! local variables
  character*50 :: tmpchar(3)
  real*8       :: tmpvals(2)
  integer*4    :: imsc
  character    :: line*90

! file should already be at the first line of the mascons, so we simply need to read them
  do imsc=1,n_msc
    read(lu_addnorm,'(a90)')line
    read(line(67:90),*)msc_vals(imsc,1:2)
  enddo

  return
  end subroutine read_fitfile_msc_vals
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





