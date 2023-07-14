 program calc_daily_mean_acc

! program to compute a model to calculate the daily mean of the biased acc
! observations once the thrusts have been removed via a low pass filter
!
! R. McGirr
! 12 August 2020

  implicit none

! command line arguments
  character*150 :: acc_file         ! input ACC1B file (format can be GRACE, GRACE FO or ???? )
  character*20  :: model_type(3)       ! character string to define the model to fit ("quadratic", "exponential", [add more] )

! unit numbers
  integer*4     :: luacc_in

! space gravity mission identifier (0: GRACE, 1: GRACE FO, 2: GRACE II)
  integer*4     :: mission

! accelerometer arrays and start/stop times
  integer*4                :: acc_span(2)                     ! start/stop time (in grace seconds) of accelerometer data
  real(kind=8),allocatable :: acc_obs(:,:,:)       ! array to store the ACC1B accelerometer observations
  real(kind=8),allocatable :: acc_obs_fixed(:,:,:) ! array to store the fixed/filtered/detrended accelerometer observations
  integer*4                :: nvals_acc
  integer*4                :: n_extend                   ! n epochs of extended data

! variables related to the rfft and cos filter
  integer*4                :: nfft
  real(kind=4)             :: flo             ! frequencies below flo are removed
  real(kind=4)             :: fhi             ! frequencies above fhi are kept
  real(kind=8)             :: mean
! other variables
  integer*4     :: icomp,iepoch,i
  character     :: calling_prog*17
  character*250 :: message

! debug
  logical debug

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  luacc_in = 10
  calling_prog = "UTIL/CALC_DAILY_MEAN_ACC"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! decode command line
  call getarg(1,acc_file)
  if(acc_file(1:1) == "")then
    print*,"calc_dail_mean_acc acccelerometer_file model_X model_Y model_Z"
    print*,"e.g. calc_dail_mean_acc ACC1B_2016-07-28_A_02.asc none lowpass none"
    print*," "
    stop
  endif

! character string for which model to use
  call getarg(2,model_type(1))
  call getarg(3,model_type(2))
  call getarg(4,model_type(3))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  A C C 1 B   D A T A
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! open files
  call level1B_open(luacc_in,calling_prog,acc_file)
! read the ACC1B header
  mission = -99
  call acc_read_hdr(luacc_in,calling_prog,acc_file,mission,acc_span)

! now we can dimension the arrays to hold the accelerometer observations.
  nvals_acc = acc_span(2) - acc_span(1) + 1
  allocate(acc_obs(1,nvals_acc,4))              ! gracesec, XYZ accelerations
  allocate(acc_obs_fixed(1,nvals_acc,4))        ! gracesec, XYZ accelerations

! read the ACC1B linear accelerations
  call acc_read_data(luacc_in,calling_prog,acc_file,mission,nvals_acc,3,acc_span,acc_obs)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   R E M O V E    H I G H    F R E Q U E N C Y    C O M P O N E N T S
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do icomp = 1,3
! *** do nothing ***
    if(model_type(icomp)(1:4) == "none")then
      write(message,'(a,i2,a)')" No model applied to component",icomp,"."
      call status_update('STATUS','UTIL',calling_prog,' ',message,0)
      acc_obs_fixed(1,:,icomp+1) = acc_obs(1,:,icomp+1)
! *** fft cosine filter ***
    else if (model_type(icomp)(1:7) == "lowpass")then
      debug = .false.
      nfft = 2**20

      ! set low and high frequency cutoffs
      flo = 0.2d-3 ! close to 1/rev frequency
      fhi = 0.3d-3
      acc_obs_fixed(1,:,icomp+1) = acc_obs(1,:,icomp+1)
      call acc_low_pass(debug,calling_prog,icomp,nvals_acc,nfft,acc_obs(1,:,icomp+1),acc_obs_fixed(1,:,icomp+1),flo,fhi)
      write(message,'(a,i1,a,e12.5,a)')"  Axis ",icomp,": filtered high frequency components above ",fhi," cycles/sec"
      call status_update('STATUS','UTIL',calling_prog,' ',message,0)

! *** model not coded. ***
    else
      write(message,'(a,a,a)')'Model "',model_type(icomp),'" not coded. You will have to write it yourself ...'
      call status_update('FATAL','UTIL','calc_daily_mean_acc',' ',message,0)
    endif

  enddo

! RM190815: don't write extended part of accelerometer obs
  if(acc_file(1:6) .eq. "extend")then
    n_extend = (nvals_acc - 86400) / 2
  else 
    n_extend = 0
  endif

! Find mean of low passed observations
  do icomp = 1,3
    if (model_type(icomp)(1:7) == "lowpass")then
      mean = sum(acc_obs_fixed(1,n_extend+1:nvals_acc-n_extend,icomp+1))/dble(nvals_acc-(2*n_extend))
      write(message,'(a,i1,a,e12.5)')"  Mean of component ",icomp,": ",mean*1.d9
      call status_update('STATUS','UTIL',calling_prog,' ',message,0)
    endif
  enddo

! PT180702: we need this program to output something ..... how about in ACC1B format without the header?
  print*,"GRACE SECONDS        Acc_X              Acc_Y                Acc_Z     (um/s^2)      Models applied: ",model_type
  do iepoch = 1+n_extend,nvals_acc-n_extend
    write(*,'(i12,6e25.15)')nint(acc_obs(1,iepoch,1)),(acc_obs_fixed(1,iepoch,i),i=2,4),(acc_obs(1,iepoch,i),i=2,4)
  enddo

  end



