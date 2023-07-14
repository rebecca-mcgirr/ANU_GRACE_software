  program calc_range_accel

! program to calculate analytically the range acceleration, given position and velocity of each GRACE satellite
!
! P. Tregoning
! 17 August 2017

  implicit none

! command line arguments
  character    :: infile_A*100,infile_B*100     ! input files containing GRACE A and GRACE B position/velocity
  character    :: outfile*100                   ! output range acceleration computations

! position/velocity variables
  real(kind=8), allocatable :: pos_A(:,:), pos_B(:,:) ! cartesian positions
  real(kind=8), allocatable :: vel_A(:,:), vel_B(:,:) ! cartesian velocities
  real(kind=8), allocatable :: acc_A(:,:), acc_B(:,:) ! cartesian accelerations
  real(kind=8), allocatable :: gracesec(:)
  real(kind=8), allocatable :: gracesec_acc(:)

! KBR arrays
  real(kind=8), allocatable :: KBR(:),KBRR(:),KBRA(:)
  real(kind=8), allocatable :: dRdt(:), dRRdt(:)

! spline variables
  real(kind=8), allocatable :: b(:,:,:),c(:,:,:),d(:,:,:)
  real(kind=8), allocatable :: acc_A_interp(:,:),acc_B_interp(:,:)
  real(kind=8) :: ispline

! counters
  integer*4  :: n_obs                                 ! number of observations in the input pos/vel files
  integer*4  :: i,j,iobs

! unit numbers
  integer*4  :: luin_A,luin_B,luout

! other stuff
  character*200 :: message,line
  integer*4     :: ioerr,junk
  real*8        :: tmp_kbra,amag3

! debug vectors
  real(kind=8), allocatable :: acc_A_deriv(:,:),acc_B_deriv(:,:)


! decode command line
  call getarg(1,infile_A)
  if(infile_A(1:1) == "")then
    print*,"Runstring: calc_range_accel truth_A.posvel truth_B.posvel truth_AB.KBRA"
    stop
  endif
  call getarg(2,infile_B)
  call getarg(3,outfile)

! open the files
  luin_A = 10
  luin_B = 11
  luout  = 12
  open(luin_A,file=infile_A,status='old')
  open(luin_B,file=infile_B,status='old')
  open(luout ,file=outfile,status='unknown')
! read header lines
  read(luin_A,'(a)')line
  read(luin_B,'(a)')line

! find out how many values there are in the file
  ioerr = 0
  n_obs = 0
  do while (ioerr == 0)
    read(luin_A,*,iostat=ioerr,end=1000)junk
    if(ioerr == 0)n_obs = n_obs+1
  enddo
1000  print*,'There are ',n_obs,' position/velocity epochs in file ',infile_A

! allocate arrays
  ! pos/vel
  allocate(pos_A(n_obs,3))
  allocate(vel_A(n_obs,3))
  allocate(acc_A(n_obs,3))
  allocate(pos_B(n_obs,3))
  allocate(vel_B(n_obs,3))
  allocate(acc_B(n_obs,3))

! spline interpolated values
  allocate(acc_A_interp(n_obs,3))
  allocate(acc_B_interp(n_obs,3))
  allocate(b(n_obs,2,3))
  allocate(c(n_obs,2,3))
  allocate(d(n_obs,2,3))

  ! range acceleration
  allocate(KBR(n_obs))
  allocate(KBRR(n_obs))
  allocate(KBRA(n_obs))
! time derivatives
  allocate(dRdt(n_obs))
  allocate(dRRdt(n_obs))
  allocate(acc_A_deriv(n_obs,3))
  allocate(acc_B_deriv(n_obs,3))

  ! gracesec
  allocate(gracesec(n_obs))
  allocate(gracesec_acc(n_obs))

! rewind the GRACE A file, then re-read the header line
  rewind(luin_A)
  read(luin_A,'(a)')line

! now read in the position/velocity information
  do iobs = 1,n_obs
    read(luin_A,*)gracesec(iobs),pos_A(iobs,:),vel_A(iobs,:),acc_A(iobs,:)
    read(luin_B,*)gracesec(iobs),pos_B(iobs,:),vel_B(iobs,:),acc_B(iobs,:)
  enddo
  print*,'Have read in the position/velocity data'

! PT170820: calculate the accelerations with time derivatives of velocity .. just to see ...
  do i=1,3
    call noise_robust_deriv(vel_A(:,i) ,  acc_A_deriv(:,i), 5.d0, n_obs, 5)
    call noise_robust_deriv(vel_B(:,i) ,  acc_B_deriv(:,i), 5.d0, n_obs, 5)
  enddo


! PT170825: 
  do iobs = 1,n_obs
    do j=1,3
! take the average of the point and the previous one .... in other words, assume the acceleration points are midway through the epoch
!      acc_A(iobs,j) = (acc_A(iobs,j)+acc_A(iobs+1,j))/2.d0
!      acc_B(iobs,j) = (acc_B(iobs,j)+acc_B(iobs+1,j))/2.d0


    enddo
! create an array of epochs that are 2.5 seconds earlier, so that the accelerations from GRACEORB are midway through the epoch step (and not assigned to the i+1th epoch)
    gracesec_acc(iobs) = gracesec(iobs) - 2.5
  enddo

! move it back one epoch
  do iobs = 2, n_obs-1
    do j=1,3
       acc_A(iobs-1,j) = acc_A(iobs,j)
       acc_B(iobs-1,j) = acc_B(iobs,j)
    enddo
  enddo

! PT170828: fit a cubic spline to each of the acceleration time series, so that we can interpolate accelerations to the i, i+1, i+2 ... epochs of pos/vel
  do j=1,3
    call spline(gracesec_acc, acc_A(:,j), b(:,1,j), c(:,1,j), d(:,1,j), n_obs)
    call spline(gracesec_acc, acc_B(:,j), b(:,2,j), c(:,2,j), d(:,2,j), n_obs)
  enddo

!! now interpolate the cubic splines onto the same epochs as the pos/vel values
!  do iobs = 1,n_obs
!    do j=1,3
!      acc_A_interp(iobs,j) = ispline(gracesec(iobs), gracesec_acc, acc_A(:,j), b(:,1,j), c(:,1,j), d(:,1,j), n_obs)
!      acc_B_interp(iobs,j) = ispline(gracesec(iobs), gracesec_acc, acc_B(:,j), b(:,2,j), c(:,2,j), d(:,2,j), n_obs)
!    enddo
!! print*,iobs,acc_A(iobs,1),acc_A_interp(iobs,1)
!  enddo  


! DEBUG
! PT170825: uncomment here to use dV/dt rather than accelerations from graceorb
!  acc_A = acc_A_deriv
!  acc_B = acc_B_deriv
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                            !!
!! compute range acceleration !!
!!                            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  kbr = 0.d0
  kbrr = 0.d0
  kbra = 0.d0
  do iobs = 1,n_obs
    ! calculate the range
!    kbr(iobs) = dsqrt( (pos_A(iobs,1)-pos_B(iobs,1))**2 + (pos_A(iobs,2)-pos_B(iobs,2))**2 &
!                         + (pos_A(iobs,3)-pos_B(iobs,3))**2 )
    kbr(iobs) = amag3(pos_A(iobs,:)-pos_b(iobs,:))


    ! calculate the range rate by the analytical expression
    ! kbrr = sum( dP_{xyz} * dV_{xyz} ) / range
    do i=1,3
      kbrr(iobs) = kbrr(iobs) + ( (pos_A(iobs,i) - pos_B(iobs,i)) * (vel_A(iobs,i)-vel_B(iobs,i)) ) /kbr(iobs)
!print*,iobs,i,kbrr(iobs),pos_A(iobs,:),pos_B(iobs,:)
    enddo

    ! analytical expression for range acceleration derived in sympy, but similar to Thomas (1999, page 3-16)
      ! kbra = 1/R ( -Rdot^2 + sum(delta_p*delta_A) + sum(delta_V^2) )
    tmp_kbra = 0.d0
    do i=1,3
      tmp_kbra = tmp_kbra   &
                 + (pos_A(iobs,i)-pos_B(iobs,i))*(acc_A(iobs,i)-acc_B(iobs,i))  &
                 + (vel_A(iobs,i)-vel_B(iobs,i) )**2
! PT170825: swap the signs of v^2 to match the formula of Chen et al (2008)
!      tmp_kbra = tmp_kbra   &
!                 + (pos_A(iobs,i)-pos_B(iobs,i))*(acc_A(iobs,i)-acc_B(iobs,i))  &
!                 - (vel_A(iobs,i)-vel_B(iobs,i) )**2
    enddo
    kbra(iobs) = 1/kbr(iobs) * (-1.d0*kbrr(iobs)**2 + tmp_kbra)  
! PT170825: swap the sign of RR^2 to match Chen et al (2008)
!    kbra(iobs) = 1/kbr(iobs) * (1.d0*kbrr(iobs)**2 + tmp_kbra)  

  enddo


!!!!!!!!!!!!!!!!!!!!!!
!!                  !!
!! time derivatives !!
!!                  !!
!!!!!!!!!!!!!!!!!!!!!!
  call noise_robust_deriv(kbr ,  dRdt, 5.d0, n_obs, 7)
  call noise_robust_deriv(dRdt, dRRdt, 5.d0, n_obs, 7)

! output the computed range, range rate and range acceleration + the time derivatives for comparison
  do iobs = 1,n_obs
    write(luout,'(f11.1,5e24.15,i8,6e24.15)')gracesec(iobs),kbr(iobs),kbrr(iobs),kbra(iobs),dRdt(iobs),dRRdt(iobs),iobs &
           ,acc_A(iobs,:),acc_A_deriv(iobs,:)
  enddo


  end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine spline (x, y, b, c, d, n)
!======================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function fspline can be used for interpolation
!======================================================================
implicit none
integer n
double precision x(n), y(n), b(n), c(n), d(n)
integer i, j, gap
double precision h

gap = n-1
! check input
if ( n < 2 ) return
if ( n < 3 ) then
  b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
  c(1) = 0.
  d(1) = 0.
  b(2) = b(1)
  c(2) = 0.
  d(2) = 0.
  return
end if
!
! step 1: preparation
!
d(1) = x(2) - x(1)
c(2) = (y(2) - y(1))/d(1)
do i = 2, gap
  d(i) = x(i+1) - x(i)
  b(i) = 2.0*(d(i-1) + d(i))
  c(i+1) = (y(i+1) - y(i))/d(i)
  c(i) = c(i+1) - c(i)
end do
!
! step 2: end conditions 
!
b(1) = -d(1)
b(n) = -d(n-1)
c(1) = 0.0
c(n) = 0.0
if(n /= 3) then
  c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
  c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
  c(1) = c(1)*d(1)**2/(x(4)-x(1))
  c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
end if
!
! step 3: forward elimination 
!
do i = 2, n
  h = d(i-1)/b(i-1)
  b(i) = b(i) - h*d(i-1)
  c(i) = c(i) - h*c(i-1)
end do
!
! step 4: back substitution
!
c(n) = c(n)/b(n)
do j = 1, gap
  i = n-j
  c(i) = (c(i) - d(i)*c(i+1))/b(i)
end do
!
! step 5: compute spline coefficients
!
b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
do i = 1, gap
  b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
  d(i) = (c(i+1) - c(i))/d(i)
  c(i) = 3.*c(i)
end do
c(n) = 3.0*c(n)
d(n) = d(n-1)
end subroutine spline

  function ispline(u, x, y, b, c, d, n)
!======================================================================
! function ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! ispline = interpolated value at point u
!=======================================================================
implicit none
double precision ispline
integer n
double precision  u, x(n), y(n), b(n), c(n), d(n)
integer i, j, k
double precision dx

! if u is ouside the x() interval take a boundary value (left or right)
if(u <= x(1)) then
  ispline = y(1)
  return
end if
if(u >= x(n)) then
  ispline = y(n)
  return
end if

!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
i = 1
j = n+1
do while (j > i+1)
  k = (i+j)/2
  if(u < x(k)) then
    j=k
    else
    i=k
   end if
end do
!*
!  evaluate spline interpolation
!*
dx = u - x(i)
ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
end function ispline     
