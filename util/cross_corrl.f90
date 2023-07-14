  program cross_corrl

! front-end program to calculate the cross correlation of two time series and, in particular, the time lag of the second wrt the first
!
! P. Tregoning
! 12 September 2018

  implicit none

  real(kind=8),allocatable :: series1(:),series2(:)
  integer*4                :: nvals,i,j,ioerr
  integer*4                :: delay
  character*100            :: line,file1,file2

! blah
  character*1 :: sat
  integer*4 :: gracesec
  real*8 :: acc(3)

! decode command line
  call getarg(1,file1)
  open(10,file=file1,status='old')
  call getarg(2,file2)
  open(11,file=file2,status='old')

! determine the number of elements in the first series
  line= " "
  nvals = 0
  ioerr = 0
  do while (ioerr == 0)
    read(10,*,iostat=ioerr,end=1000)acc(3)
    if(ioerr == 0 .and. dabs(acc(3)) > 520. .and. dabs(acc(3)) < 540.)nvals = nvals + 1
  enddo
1000 print*,"There are ",nvals," values in file ",file1
  rewind(10)

! allocate arrays
  allocate(series1(nvals))
  allocate(series2(nvals))

! read both series
  do i=1,nvals
    read(10,*)series1(i)
    read(11,*)series2(i)
  enddo


! ok, so now calculate the cross-correlation between the two series
  print*,"Compute the cross correlation ...."
  call xcorr(series1,series2,nvals,delay)
  print*,"delay, ccmax and ccpol",delay

  end
!----------------------------------------------------------------------------------------------


  subroutine xcorr(x, y, n,  delay)
  ! Cross-correlation of array x and y.
  ! Return time shift
  implicit none
  real(8) :: crosscorr
  integer :: n, shift, delay, ccpol
  real(8), dimension(0:n-1) :: x, y
  real(8) :: cc, ccmax, ccmin
  integer :: k, kmin

  shift = 1
  ccmax = 0
  ccmin = 0
  kmin = 0
!  do k = -n+1,n-1,shift
  do k = 1,6000
      cc = dot_product(x(1+k:6000+k), y(1:6000))
      cc = cc/sqrt(cc)
print*,'k,cc',k,cc,ccmax
     if (cc.gt.ccmax) then
         ccmax = cc
         delay = k
print*,'updated delay           : k,cc',k,cc
     endif
     if (cc.lt.ccmin) then
         ccmin = cc
         kmin = k
     endif
  enddo
  if (ccmax.gt.-ccmin) then
     ccpol = 1
  else
     ccmax = -ccmin
     delay = kmin
     ccpol = -1
  endif
  ccmax = ccmax/sqrt( dot_product(x,x) * dot_product(y,y) ) 

  end subroutine xcorr




