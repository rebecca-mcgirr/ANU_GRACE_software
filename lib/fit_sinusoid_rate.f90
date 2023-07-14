  subroutine fit_sinusoid_rate(calling_prog,nobs,omega,values,model_params)

! subroutine to generate model values for a fit of a sinusoid, offset and rate to a vector of data
!
! P. Tregoning
! 8 August 2019

  implicit none

  character(*), intent(in)  :: calling_prog
  integer*4,    intent(in)  :: nobs                  ! number of values in the vector of data
  real(kind=8), intent(in)  :: omega                 ! period of sinusoidal function
  real(kind=8), intent(in)  :: values(nobs)          ! the actual data
  real(kind=8), intent(out) :: model_params(4)       ! parameter estimates (offset, rate, AmpC,AmpS

! local arrays
  real(kind=8),allocatable  :: Amat(:,:), Bmat(:,:),At(:,:),AtA(:,:),VCV(:,:),AtB(:,:),xhat(:,:)  ! LS matrices
  real(kind=8) :: offset, rate, AmpC,AmpS     ! parameters            
  integer*4    :: nparam,i,iter

  nparam = 4
! allocate LS arrays
  allocate(Amat(nobs,nparam))
  allocate(Bmat(nobs,1))
  allocate(At(nparam,nobs))
  allocate(AtA(nparam,nparam))
  allocate(VCV(nparam,nparam))
  allocate(AtB(nparam,1))
  allocate(xhat(nparam,1))

! define a priori parameter values
  offset = 0.d0
  rate = 0.d0
  AmpC = 0.d0
  AmpS = 0.d0

! form the Amat and Bmat
  do iter = 1,2
    do i=1,nobs
      Amat(i,1) = 1.d0
      Amat(i,2) = dble(i)
      Amat(i,3) = dcos(omega*dble(i))
      Amat(i,4) = dsin(omega*dble(i))

      Bmat(i,1) = values(i) - (offset + rate*dble(i) + AmpC*dcos(omega*dble(i)) + AmpS*dsin(omega*dble(i)))
    enddo

    ! LS solution
    call transp(Amat,At,nobs,nparam)
    call matmult(At,Amat,AtA,nparam,nobs,nparam)
    call invert(AtA,VCV,nparam)
    call matmult(At,Bmat,AtB,nparam,nobs,1)
    call matmult(VCV,AtB,xhat,nparam,nparam,1)

    offset = offset + xhat(1,1)
    rate = rate + xhat(2,1)
    AmpC = AmpC + xhat(3,1)
    AmpS = AmpS + xhat(4,1)
  enddo

  model_params(1) = offset
  model_params(2) = rate
  model_params(3) = AmpC
  model_params(4) = AmpS

  return
  end subroutine fit_sinusoid_rate


