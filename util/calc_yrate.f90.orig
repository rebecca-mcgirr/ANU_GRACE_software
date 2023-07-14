  subroutine calc_yrate(nepochs,accobs,axis,rate)

! subroutine to calculate the rate of the y-axis accelerations using least squares
!
! P. Tregoning
! 27 November 2015

  implicit none

! passed variables
  integer*4,   intent(in)  :: nepochs
  integer*4,   intent(in)  :: axis                 ! axis for which to estimate the rate (1: X, 2: Y, 3: Z)
  real(kind=8),intent(in)  :: accobs(nepochs,5)   
  real(kind=8),intent(out) :: rate

! local variables
  integer*4     :: i,j,counter,iepoch
! local variables for LS
  integer*4                :: nobs,nparm,iter
  real(kind=8),allocatable :: A(:,:),At(:,:),AtA(:,:),VCV(:,:),B(:,:),AtB(:,:),soln(:,:)
  real(kind=8)             :: offset

! we have two parameters: the rate and the offset
  nparm = 2

! read through the accobs and extract out the non-thrust observations. We will use them all !!!
  nobs = 0
  do iepoch = 1,nepochs
    if(accobs(iepoch,5) /= 1.) then    ! it is a non-thrust observation
      nobs = nobs + 1
    endif
  enddo

! allocate the arrays
  allocate(A(nobs,nparm))
  allocate(At(nparm,nobs))
  allocate(AtA(nparm,nparm))
  allocate(VCV(nparm,nparm))
  allocate(B(nobs,1))
  allocate(AtB(nparm,1))
  allocate(soln(nparm,1))

! form up the arrays for LS
  rate = 0.
  offset = 0.
  do iter = 1,2
    counter = 0
    do iepoch = 1,nepochs
      if(accobs(iepoch,5) /= 1.  ) then    ! it is a non-thrust observation
        counter = counter + 1
        A(counter,1) = 1.
        A(counter,2) = iepoch*5.0
        B(counter,1) = accobs(iepoch,1+axis) - (offset + rate*(iepoch*5))
      endif
    enddo

! LS soln
    call transp(A,At,nobs,nparm)
    call matmult(At,A,AtA,nparm,nobs,nparm)
    call invert(AtA,VCV,nparm)
  
    call matmult(At,B,AtB,nparm,nobs,1)
    call matmult(VCV,AtB,soln,nparm,nparm,1)

    offset = offset + soln(1,1)
    rate = rate + soln(2,1)
  enddo

  return
  end


