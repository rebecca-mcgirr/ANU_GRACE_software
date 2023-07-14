program main

! program to filter the low frequency components of the uncalibrated accelerometer data, especially in the 
! cross-track. This is achieved with the forward discrete FFT to transform the accelerometer data to the 
! frequency domain (analysis), the cosine high-pass fileter is then used to filter the low frequencies,
! finally, the backward FFT is used to transform the filtered frequencies back to the time domain (synthesis).
!
! This test program is based on model_acc_v2
!
! R. McGirr
! 26 June 2019

  implicit none

! command line arguments
  character*150               :: acc_file               ! input ACC1B file
  character*20                :: model_type(3)          ! character string to define the model to fit
  integer*4                   :: start_ep(3),end_ep(3)  ! start and end epochs over which to fit the model

! unit numbers
  integer*4                   :: luacc_in

! accelerometer arrays and start/stop times
  integer*4                   :: acc_span(2)            ! start/stop time (in grace seconds) of accelerometer data
  integer*4                   :: nvals_acc
  real(kind=8),allocatable    :: acc_obs(:,:,:)         ! array to store the ACC1B accelerometer observations
  real(kind=8),allocatable    :: acc_obs_fixed(:,:,:)   ! array to store the fixed/filtered/detrended accelerometer

! variables related to the rfft and filtering
  integer*4                   :: nfft
  real(kind=4)                :: flo
  real(kind=4)                :: fhi
  real(kind=4),allocatable    :: p(:)
  real(kind=8),allocatable    :: f(:)
  character*20                :: ftype  

! other variables
  integer*4                   :: narg,i,icomp,iepoch,j
  integer*4                   :: mission                ! mission identifier
  character*20                :: calling_prog                
  character*100               :: arg
  character*250               :: message

! debug
  logical                     :: debug

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  luacc_in = 11
  calling_prog = "RFFT_ACC_TEST"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! decode command line
  call getarg(1,acc_file)
  if(acc_file(1:1) == "")then
    print*,"rfft_acc_test acc_file  model_X start end model_Y start end model_Z start end"
    print*,"e.g. rfft_acc_test ACC1B_2016-07-28_A_02.asc none start end rfft start end rfft start end"
    print*," "
    stop
  endif

! character string for which model to use
  narg = 2
  do i=1,3
    call getarg(narg,model_type(i))
    call getarg(narg+1,arg)
    read(arg,*)start_ep(i)
    call getarg(narg+2,arg)
    read(arg,*)end_ep(i)
    narg = narg + 3
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! open files
  call level1B_open(luacc_in,calling_prog,acc_file)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read files
! read the ACC1B header
  mission = -99
  call acc_read_hdr(luacc_in,calling_prog,acc_file,mission,acc_span)

! now we can dimension the arrays to hold the accelerometer observations.
  nvals_acc = acc_span(2) - acc_span(1) + 1
  allocate(acc_obs(1,nvals_acc,5))              ! gracesec, XYZ accelerations, 5th component is a thrust flag
  allocate(acc_obs_fixed(1,nvals_acc,5))        ! gracesec, XYZ accelerations, 5th component is a thrust flag

! set all the acc obs flags to good a priori
  acc_obs(1,:,5) = 1

! read the ACC1B linear accelerations
  call acc_read_data(luacc_in,calling_prog,acc_file,mission,nvals_acc,3,acc_span,acc_obs)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! remove some sort of model
  do icomp = 1,3
! do nothing
    if(model_type(icomp)(1:4) == "none")then
      write(message,'(a,i2,a)')" No model applied to component",icomp,"."
      call status_update('STATUS','UTIL',calling_prog,' ',message,0)
      acc_obs_fixed(1,:,icomp+1) = acc_obs(1,:,icomp+1)

! rfft/cosine highpass filter
    else if (model_type(icomp)(1:4) == "rfft")then
      debug = .true.
      nfft = nvals_acc*1 
      ftype = 'high'
      flo = 0.0
      fhi = 1.0/5670.0
      allocate(p(nvals_acc/2 + 1))
      allocate(f(nvals_acc/2 + 1))      
      call rfft(debug,calling_prog,icomp,nvals_acc,nfft,acc_obs(1,:,icomp+1),acc_obs_fixed(1,:,icomp+1),ftype,flo,fhi,p,f)
      write(message,'(a,i1,a)')"  Axis ",icomp,": filtered high frequency components"
      call status_update('STATUS','UTIL',calling_prog,' ',message,0)

! model not coded
    else
      write(message,'(a,a,a)')'Model "',model_type,'" not coded. You will have to write it yourself ...'
      call status_update('FATAL','UTIL','model_acc_v2',' ',message,0)
    endif
  enddo

! print out filtered ACC1B
  print*,"GRACE SECONDS        Acc_X              Acc_Y                Acc_Z     (um/s^2)      Models applied: ",model_type
  do iepoch = 1,10
    write(*,'(i12,6e25.15)')nint(acc_obs(1,iepoch,1)),(acc_obs_fixed(1,iepoch,j),j=2,4),(acc_obs(1,iepoch,j),j=2,4)
    !write(*,'(2e25.15)')acc_obs(1,iepoch,3),acc_obs_fixed(1,iepoch,3)
  enddo

  stop
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rfft(debug,calling_prog,icomp,n,nfft,r,r_out,ftype,flo,fhi,p,f)

! Perform Fourier analysis on accelerometer data, then filter low frequencies with a
! cosine highpass filter, finally, perform Fourier synthesis to transform frquencies
! to accelerations
!
! R. McGirr
! 26 June 2019

  implicit none

! passed variables
  logical     ,intent(in)     :: debug                   ! debug flag
  character*20,intent(in)     :: calling_prog            ! name of calling program
  integer*4   ,intent(in)     :: icomp                   ! component (X,Y or Z accelerations)
  integer*4   ,intent(in)     :: n
  integer*4   ,intent(in)     :: nfft
  real(kind=4),intent(in)     :: flo
  real(kind=4),intent(in)     :: fhi
  real(kind=8),intent(in)     :: r(n)
  real(kind=8),intent(out)    :: r_out(n)
  real(kind=4),intent(out)    :: p(n/2 + 1)
  real(kind=8),intent(out)    :: f(n/2 + 1)

! local variables
  integer*4                   :: lensav
  integer*4                   :: lenwrk
  integer*4                   :: ier
  integer*4                   :: inc
  integer*4                   :: lenr
  real(kind=8)                :: mean
  real(kind=4)                :: r_tmp(n)
  real(kind=4)                :: r_filt(n)
  real(kind=4),allocatable    :: filt(:)
  real(kind=4),allocatable    :: work(:)
  real(kind=4),allocatable    :: wsave(:)
  character*20                :: ftype  

  integer*4                   :: i,j

  write(*,'(a)')'rfft_test'
  write(*,'(a,i8)')'  The number of data items is N = ', n

! calculate and remove a mean value so that de-meaned data are sent to the forward rfft
  mean = sum(r,n)/dble(n)
  r_tmp = r - mean
  if(debug)print*,"mean value is",mean,sum(r,n),dble(n)
  !do i = 1,10
  !  write(*,'(2e25.15)')r(i),r_tmp(i)
  !enddo

!  Set work vectors.
  lensav = n + int(log(real(n))/log( 2.0E+00 )) + 4
  lenwrk = n

  write(*,'(a,i8)')'  LENSAV = ', lensav
  write(*,'(a,i8)')'  LENWRK = ', lenwrk

  allocate(work(lenwrk))
  allocate(wsave(lensav))

  call rfft1i(nfft, wsave, lensav, ier)

! Compute the FFT coefficients.
  inc = 1
  lenr = nfft
  !print*,"data = "
  !do i = 1,10
  !  write(*,'(e25.15)')r_tmp(i)
  !enddo
  call rfft1f(nfft, inc, r_tmp, lenr, wsave, lensav, work, lenwrk, ier)   
  !print*,"FFT coeffs = "
  !do i = 1,10
  !  write(*,'(e25.15)')r_tmp(i)
  !enddo

! Fill in frequency array (cycles/second)
  print*,"Length of sample frequency array",n/2 + 1
  do i = 1,n/2 + 1
    f(i) = (i-1) * 1./(n*inc)
  enddo

! Cosine high pass filter
  allocate(filt(n/2 + 1))
  flo = 0.0
  fhi = 1.0/6000.0
  ftype = 'high'
  call cos_filt(debug,calling_prog,ftype,f,n,flo,fhi,filt)

! Multiply spectrum by filter
  r_filt(1) = filt(1) * r_tmp(1)
  do i = 2,n/2
    r_filt(2*i-1) = filt(i) * r_tmp(2*i-1)
    r_filt(2*i) = filt(i) * r_tmp(2*i)
  enddo
  r_filt(n) = filt(n/2 + 1) * r_tmp(n)

! Convert to power spectra for plotting
  j = 0
  p(1) = r_filt(1)**2
  do i = 2,n/2
    j = j + 2
    p(i) = r_filt(j)**2 + r_filt(j+1)**2
  enddo
  p(n/2 + 1) = r_filt(n)**2
  
  !do i = 1,n/2 + 1
    !write(*,'(2e25.15)')f(i),p(i)
  !enddo

!  Compute inverse FFT of coefficients.  Should get back the
!  original data.

  call rfft1b(n, inc, r_filt, lenr, wsave, lensav, work, lenwrk, ier)

  deallocate(work)
  deallocate(wsave)

  r_out = r_filt + mean

  !do i = 1,10
  !  write(*,'(3e25.15)')r(i),r_filt(i),r_out(i)
  !enddo

  return
end

