!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  acc_read_hdr         : read the header information for an ACC1B file
!!  acc_read_data        : read the actual linear acceleration information 
!!  acc_form_quadratic   : form A and B matrices for fitting a quadratic expression to the data
!!  acc_form_exponential : fit an exponential model (including a quadratic term) instead of a quadratic
!!  acc_EMD1             : remove long-wavelength components from ACC1B using the EMD filter (gracefit version)
!!  acc_EMD2             : as per acc_EMD1 but uses a stopping criteria s number in the sifting process
!!  acc_inflexions       : determine whether EMD components have inflexion points within 1/rev (or are long wavelength instead)
!!  acc_zero_crossings   : determine the frequency content of each EMD component by how many times the time series crosses zero
!!  acc_emd_psd          : compute psd of each EMD component and, from that, decide whether to keep the component
!!  acc_cos_fft          : remove long-wavelength frequencies via FFT analysis/synthesis and cos highpass filtering
!!  acc_low_pass         : use the same method as acc_cos_fft but use a low pass filter (i.e for thrust removal)
! 
!   P. Tregoning
!   31 May 2018
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine acc_low_pass(debug,calling_prog,icomp,nvals,nfft,acc_in,acc_out,flo,fhi)

! Perform Fourier analysis on accelerometer data, then filter low frequencies
! with a cosine high pass filter, finally, perform Fourier synthesis to transform 
! frequencies back to a time series of accelerations
!
! note: nfft should be set to a value equal to or greater than nvals, the larger
! the number the 'cleaner' the spectrum will appear, but the calculations will be 
! slower. The FFT is most efficient if nfft is a power of two. 
!
! R. McGirr
! 12 August 2020
!

  implicit none

! passed variables
  logical     ,intent(in)     :: debug                   ! debug flag
  character(*),intent(in)     :: calling_prog            ! name of calling program
  integer*4   ,intent(in)     :: icomp                   ! acc component
  integer*4   ,intent(in)     :: nvals                   ! lengh of input data
  integer*4   ,intent(in)     :: nfft                    ! see note in subroutine header
  real(kind=4),intent(in)     :: flo                     ! low frequency cutoff
  real(kind=4),intent(in)     :: fhi                     ! high frequency cutoff
  real(kind=8),intent(in)     :: acc_in(nvals)           ! input vector time series
  real(kind=8),intent(out)    :: acc_out(nvals)          ! input vector time series

! local variables
  integer*4                   :: i
  real(kind=8)                :: mean
  character*9                 :: spec_outfile            ! soectrum outfile
  character*10                :: ftype                   ! set to 'high'
  character*10                :: wtype                   ! set to 'hann'
  real(kind=4),allocatable    :: pspec_raw(:)            ! contains power spectrum
  real(kind=4),allocatable    :: pspec_filt(:)           ! contains filtered power spectrum
  real(kind=4),allocatable    :: fspec(:)                ! contains frequencies
  real(kind=4),allocatable    :: filt(:)                 ! contains filter
  real(kind=8),allocatable    :: acc_tmp(:)              ! output filtered vector time series
  real(kind=8),allocatable    :: acc_tmp1(:)             ! input vector time series with slope between 1st/last point removed
  real(kind=4),allocatable    :: window(:)               ! contains window coefficients
  character*13                :: filter_type             ! currently either "hann_extended" or "slope"
  real(kind=8)                :: slope                   ! slope between first and last point

  call status_update("STATUS",calling_prog,"ACC_LOW_PASS","","Using ACC_LOW_PASS",0)

! allocate arrays
  allocate(pspec_raw(nfft/2 + 1))
  allocate(pspec_filt(nfft/2 + 1))
  allocate(fspec(nfft/2 + 1))
  allocate(filt(nfft/2 + 1))
  allocate(acc_tmp(nvals))
  allocate(acc_tmp1(nvals))
  allocate(window(nvals)) 

  !filter_type = 'hann_extended'
  filter_type = 'uniform'

  if(filter_type(1:13) == "hann_extended")then
    ! compute mean
    mean = sum(acc_in,nvals)/dble(nvals)

    ! compute window using hann window function. Requires extended ACC data before/after the day to be processed
    window = 0.0
    wtype = 'hann'
    call window_function(calling_prog,nvals,wtype,window)

    ! call compute_fft subroutine
    ftype = 'low'
    call compute_fft(debug,calling_prog,nvals,nfft,(acc_in-mean)*window,acc_tmp,ftype,flo,fhi,fspec,pspec_raw,pspec_filt,filt)

  else if (filter_type(1:7) == "uniform")then
    mean = sum(acc_in,nvals)/dble(nvals)

    ! use a uniform window
    window = 0.0
    wtype = 'uniform'
    call window_function(calling_prog,nvals,wtype,window)

    ! call compute_fft subroutine
    ftype = 'low'
    call compute_fft(debug,calling_prog,nvals,nfft,(acc_in-mean)*window,acc_tmp,ftype,flo,fhi,fspec,pspec_raw,pspec_filt,filt)
  endif

! remove window add in mean
  acc_out = (acc_tmp / window) + mean

! write out frequencies and power spectrum
  write(spec_outfile,'(a,i1)')"fft.spec",icomp
  open(1000,file=spec_outfile,status='unknown')
  write(1000,*)flo,fhi
  do i = 1,nfft/2 + 1
    write(1000,*)fspec(i),pspec_raw(i),pspec_filt(i),filt(i)
  enddo
  close(1000)

! deallocate arrays
  deallocate(pspec_raw)
  deallocate(pspec_filt)
  deallocate(fspec)
  deallocate(filt)
  deallocate(acc_tmp)
  deallocate(window) 

  end subroutine acc_low_pass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine acc_cos_fft(debug,calling_prog,icomp,nvals,nfft,acc_in,acc_out,flo,fhi)
! Perform Fourier analysis on accelerometer data, then filter low frequencies
! with a cosine high pass filter, finally, perform Fourier synthesis to transform 
! frequencies back to a time series of accelerations
!
! note: nfft should be set to a value equal to or greater than nvals, the larger
! the number the 'cleaner' the spectrum will appear, but the calculations will be 
! slower. The FFT is most efficient if nfft is a power of two. 
!
! R. McGirr
! 26 June 2019
!
! MODS:
! PT190819: added a second filtering option, being to remove the slope between the first/last point. The data then shouldn't need to 
!           be extended beyond the day boundaries.

  implicit none

! passed variables
  logical     ,intent(in)     :: debug                   ! debug flag
  character(*),intent(in)     :: calling_prog            ! name of calling program
  integer*4   ,intent(in)     :: icomp                   ! acc component
  integer*4   ,intent(in)     :: nvals                   ! lengh of input data
  integer*4   ,intent(in)     :: nfft                    ! see note in subroutine header
  real(kind=4),intent(in)     :: flo                     ! low frequency cutoff
  real(kind=4),intent(in)     :: fhi                     ! high frequency cutoff
  real(kind=8),intent(in)     :: acc_in(nvals)           ! input vector time series
  real(kind=8),intent(out)    :: acc_out(nvals)          ! input vector time series

! local variables
  integer*4                   :: i
  real(kind=8)                :: mean
  character*9                 :: spec_outfile            ! soectrum outfile
  character*10                :: ftype                   ! set to 'high'
  character*10                :: wtype                   ! set to 'hann'
  real(kind=4),allocatable    :: pspec_raw(:)            ! contains power spectrum
  real(kind=4),allocatable    :: pspec_filt(:)           ! contains filtered power spectrum
  real(kind=4),allocatable    :: fspec(:)                ! contains frequencies
  real(kind=4),allocatable    :: filt(:)                 ! contains filter
  real(kind=8),allocatable    :: acc_tmp(:)              ! output filtered vector time series
  real(kind=8),allocatable    :: acc_tmp1(:)             ! input vector time series with slope between 1st/last point removed
  real(kind=4),allocatable    :: window(:)               ! contains window coefficients
  character*13                :: filter_type             ! currently either "hann_extended" or "slope"
  real(kind=8)                :: slope                   ! slope between first and last point

  call status_update("STATUS",calling_prog,"ACC_FFT","","Using ACC_FFT",0)

! allocate arrays
  allocate(pspec_raw(nfft/2 + 1))
  allocate(pspec_filt(nfft/2 + 1))
  allocate(fspec(nfft/2 + 1))
  allocate(filt(nfft/2 + 1))
  allocate(acc_tmp(nvals))
  allocate(acc_tmp1(nvals))
  allocate(window(nvals)) 

  filter_type = 'hann_extended'
!  filter_type = 'slope'

  if(filter_type(1:13) == "hann_extended")then
    ! compute mean
    mean = sum(acc_in,nvals)/dble(nvals)

    ! compute window using hann window function. Requires extended ACC data before/after the day to be processed
    window = 0.0
    wtype = 'hann'
    call window_function(calling_prog,nvals,wtype,window)

    ! call compute_fft subroutine
    ftype = 'high'
    call compute_fft(debug,calling_prog,nvals,nfft,(acc_in-mean)*window,acc_tmp,ftype,flo,fhi,fspec,pspec_raw,pspec_filt,filt)

  else if (filter_type(1:5) == "slope")then
    ! PT190819: try just removing the slope between the first/last points, which would set the first/last to zero by default. Then, using
    !           a "uniform" window - with no padding - we can perform the high-pass filtering
    slope = (acc_in(nvals) - acc_in(1))/dble(nvals)
    do i=1,nvals
      acc_tmp1(i) = acc_in(i) - (acc_in(1)+dble(i)*slope)
    enddo
    mean = sum(acc_tmp1,nvals)/dble(nvals)

    ! use a uniform window
    window = 0.0
    wtype = 'uniform'
    call window_function(calling_prog,nvals,wtype,window)

    ! call compute_fft subroutine
    ftype = 'high'
    call compute_fft(debug,calling_prog,nvals,nfft,(acc_tmp1-mean)*window,acc_tmp,ftype,flo,fhi,fspec,pspec_raw,pspec_filt,filt)
  endif

! remove window add in mean
  acc_out = (acc_tmp / window) + mean

! write out frequencies and power spectrum
  write(spec_outfile,'(a,i1)')"fft.spec",icomp
  open(1000,file=spec_outfile,status='unknown')
  write(1000,*)flo,fhi
  do i = 1,nfft/2 + 1
    write(1000,*)fspec(i),pspec_raw(i),pspec_filt(i),filt(i)
  enddo
  close(1000)

! deallocate arrays
  deallocate(pspec_raw)
  deallocate(pspec_filt)
  deallocate(fspec)
  deallocate(filt)
  deallocate(acc_tmp)
  deallocate(window) 

  end subroutine acc_cos_fft

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine acc_emd_psd(debug,calling_prog,nvals,ndecomp,start_epoch,end_epoch,acc_decomp,stop_comp)

! subroutine to calculate the psd of each EMD component and then decide
! whether the component is signal or long-wavelength, non-linear stuff
! that we don't want to include
!
! P. Tregoning
! 26 September 2018

  implicit none

! passed variables
  logical       , intent(in)    :: debug                     ! debug flag  
  character*(*) , intent(in)    :: calling_prog              ! name of calling program
  integer*4     , intent(in)    :: nvals                     ! dimensioning of EMD component vector
  integer*4     , intent(in)    :: ndecomp                   ! number of EMD components
  integer*4     , intent(in)    :: start_epoch               ! epoch in the acc_decomp array for which we start the psd
  integer*4     , intent(in)    ::   end_epoch               ! epoch in the acc_decomp array for which we end   the psd
  real(kind=8)  , intent(in)    :: acc_decomp(ndecomp,nvals) ! EMD component observations
  integer*4     , intent(inout) :: stop_comp                 ! set to 1 if sub-revolution inflexions found

! local variables
  integer*4    :: iepoch,icomp,i
  integer*4    :: counter
  logical      :: pos_val
  character*100:: message

! variables for the psd computation
  integer*4,parameter :: m = 1024                          ! numbers plucked out of the air!
  integer*4,parameter :: k = 40              
  real*4              :: w1(1:4*m),w2(1:m),psd(m)
  real(kind=8)        :: sum_pwr,pwr_tol

  pwr_tol = 0.d0

! for each component, calculate a psd. Use the Numerical Recipes "spctrm". It reads data from a file
! with unit number "9"
  do icomp=1,ndecomp
    ! write the values to a file for the Numerical Recipies subroutine to read
    open(9,file='decomp.tmp',status='unknown')
    do iepoch = start_epoch,end_epoch
      write(9,*)acc_decomp(icomp,iepoch)*1.e9
    enddo
    rewind(9)
    ! now calculate the PSD
print*,'calling spctrm icomp=',icomp
!    call spctrm(psd,m,k,.true.,w1,w2)
! debug
do i=2,22,2
  print*,i,psd(i),icomp,' period psd icomp'
enddo
    ! calculate the integrated power for periods < 2 revolutions (~10000 seconds)
    sum_pwr = sum(psd(1:100))
    if(sum_pwr > pwr_tol .and. icomp >= stop_comp)then
      write(message,'(a,e12.5,a)')"Component has enough high frequency power ",sum_pwr," to be included"
      call status_update('STATUS',calling_prog,'acc_emd_psd',' ',message,0)
      stop_comp = icomp+1
    endif
    ! close the file
    close(9)
  enddo
        
stop 'stopped in acc_emd_psd'
  return
  end subroutine acc_emd_psd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine acc_zero_crossings(debug,calling_prog,icomp,nvals,n_crossings,stop_comp,acc_data)

! subroutine to determine whether the acc component is long-wavelength or
! whether it has high frequency crossings of zero (and is therefore signal)
!
! P. Tregoning
! 26 September 2018

  implicit none

! passed variables
  logical       , intent(in)    :: debug           ! debug flag  
  character*(*) , intent(in)    :: calling_prog    ! name of calling program
  integer*4     , intent(in)    :: icomp           ! actual decomposed component being tested 
  integer*4     , intent(in)    :: nvals           ! dimensioning of EMD component vector
  integer*4     , intent(inout) :: stop_comp       ! set to 1 if sub-revolution inflexions found
  real(kind=8)  , intent(in)    :: acc_data(nvals) ! EMD component observations
  integer*4     , intent(in)    :: n_crossings     ! the minimum number of zero crossings for which we include the component

! local variables
  integer*4    :: iepoch,start_epoch,end_epoch
  integer*4    :: counter
  logical      :: pos_val
  character*100:: message

! determine the start/end points in (potentially) extended data so that we count zero crossings
  if(nvals > 86400)then
    start_epoch = 1+ nint(dble(nvals - 86400))/2.d0
    end_epoch =  nvals - nint(dble((nvals - 86400))/2.d0) 
  else
    start_epoch = 1
    end_epoch = nvals
  endif
!print*,nvals,start_epoch,dble(nvals - 86400),nint(dble(nvals - 86400))

!! PT181003: trap occasions when start or end epochs have been calculated to be negative (ie < 24 hrs of data has been passed in)
!  if(start_epoch < 0)start_epoch = 1
!  if(end_epoch < 0)end_epoch = nvals

! set the pos_val flag based on the first value
  if(acc_data(start_epoch) > 0.d0)then
    pos_val = .true.
  else
    pos_val = .false.
  endif

  counter = 0
  do iepoch = start_epoch+1,end_epoch
    if(pos_val)then
      if(acc_data(iepoch) < 0.d0)then
        pos_val = .false.
        counter = counter + 1
      endif
    else
      if(acc_data(iepoch) > 0.d0)then
        pos_val = .true.
        counter = counter + 1
      endif
    endif
  enddo


  if (counter >= n_crossings)then
    write(message,'(a,i4,a,i6)')"Found ",counter," zero crossings. Keep EMD component",icomp
    call status_update('STATUS',calling_prog,'acc_zero_crossings',' ',message,0)
    stop_comp = 0
  else
    write(message,'(a,i4,a,i6)')"Found only ",counter," zero crossings. Discard EMD component",icomp
    call status_update('STATUS',calling_prog,'acc_zero_crossings',' ',message,0)
    stop_comp = icomp
  endif

  return
  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine acc_inflexions(debug,calling_prog,icomp,nvals,max_tol,n_inflexions,stop_comp,acc_data)

! subroutine to determine whether the acc component is long-wavelength or
! whether it has inflexion points within 1 revolution (and is therefore signal)
!
! P. Tregoning
! 6 August 2018

  implicit none

! passed variables
  logical       , intent(in)    :: debug           ! debug flag  
  character*(*) , intent(in)    :: calling_prog    ! name of calling program
  integer*4     , intent(in)    :: icomp           ! actual decomposed component being tested 
  integer*4     , intent(in)    :: nvals           ! dimensioning of EMD component vector
  integer*4     , intent(inout) :: stop_comp       ! set to 1 if sub-revolution inflexions found
  real(kind=8)  , intent(in)    :: acc_data(nvals) ! EMD component observations
  real(kind=8)                  :: max_tol         ! tolerance for the peak-to-peak test of the component
  integer*4     , intent(in)    :: n_inflexions    ! the minimum number of inflexions for which we include the component

! local variables
  integer*4    :: iepoch
  integer*4    :: nobs_per_rev                     ! number of 1 second observations in a revolution
  integer*4    :: dt                               ! number of epochs since last inflexion point
  integer*4    :: acc_extend                       ! number of epochs that acc_data have been extended by
  integer*4    :: end_epoch                        ! end epoch for extended acc data
  real(kind=8) :: dval                             ! change in value since last epoch
  integer*4    :: last_inflexion                   ! epoch counter of epoch of last  inflexion 
  integer*4    :: first_inflexion                  ! epoch counter of epoch of first inflexion 
  logical      :: ascending                        ! true if time series is increasing
  character*90 :: message
  logical      :: not_sure                         ! logical as to whether to keep EMD component or not
  integer*4    :: n_cases                          ! number of times a 1/rev signal is found in a component
  
! set the number of 1s observations per revolution (say 90 mins/rev)
  nobs_per_rev = 5400

! loop over all observations in this EMD component
  ascending = .false.
  not_sure = .true.
  n_cases = 0
  stop_comp = 1
  
! RM180914 if acc has been extended we want to ignore added data possibly corrupted by mode-mixing/edge-effects
  if (nvals > 86400) then
    acc_extend = nint((dble(nvals - 86400))/2.d0)
    iepoch = (acc_extend + 2)
    end_epoch = (nvals - acc_extend)
    !acc_extend = 0
    !iepoch = 2
    !end_epoch = nvals
  else
    acc_extend = 0
    iepoch = 2
    end_epoch = nvals
  end if 
  do while (iepoch <= end_epoch .and. not_sure)
    if(iepoch == acc_extend+2)then
      dt = 0
      dval = acc_data(iepoch) - acc_data(iepoch-1)
      last_inflexion = acc_extend+2
      !last_inflexion = 2
      first_inflexion = acc_extend
    else if (iepoch == acc_extend+3)then
      dval = acc_data(iepoch) - acc_data(iepoch-1)
      if(dval > 0.d0)ascending = .true.
    else
      dval = acc_data(iepoch) - acc_data(iepoch-1)
      if(dabs(acc_data(iepoch)-acc_data(last_inflexion)) > max_tol &
               .and. (dval > 0.d0 .and. .not. ascending .or. dval < 0.d0 .and. ascending))then ! change of sign from last epoch
if(debug)print*,'inflexion found. iepoch, last_inflexion,acc_data(iepoch),acc_data(last_inflexion):' &
               ,iepoch, last_inflexion,acc_data(iepoch),acc_data(last_inflexion)
        if(first_inflexion == acc_extend)then ! record the first inflexion point
if(debug)print*,"iepoch,first_inflexion",iepoch,first_inflexion
          first_inflexion = iepoch
          last_inflexion = iepoch
          not_sure = .true.
        else
!print*,'inflexion found',iepoch,last_inflexion,nobs_per_rev,acc_data(iepoch) ,acc_data(iepoch-1)
          if(iepoch - last_inflexion < nobs_per_rev/2 .and. n_cases >= n_inflexions)then  ! change within one revolution. Bail out of subroutine
            not_sure = .false.  ! we know that we want to keep this 
            stop_comp = 0
            write(message,'(a,i4,a,i6,a,2i7)')"Found ",n_inflexions," inflexion points within ",iepoch-last_inflexion &
                                         ," epochs. Keep EMD component",icomp,iepoch
            call status_update('STATUS',calling_prog,'acc_inflexions',' ',message,0)
          else
if(debug)print*,"iepoch-last_inflexion,n_cases",iepoch-last_inflexion,n_cases+1
            last_inflexion = iepoch
            n_cases = n_cases + 1
          endif
        endif
      endif
    endif
!print*,'iepoch,last_inflexion,dval',iepoch,last_inflexion,dval
    iepoch = iepoch + 1

  enddo

  return


  end subroutine acc_inflexions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine acc_EMD1(debug,calling_prog,icomp,nvals,ncomp,acc_secs,acc_data_in,acc_decomp,end_comp)

! apply the EMD filter (used in gracefit to remove high frequency noise from KBR)
! to remove the low-frequency components from the ACC1B data
!
! P. Tregoning
! 6 August 2018
!
! PT180817: modify the EMD process to try to reduce mode mixing that sees some of the 1/rev signal 
!           appearing in a lower frequency mode. We don't want that!

  implicit none

! passed variables
  logical     ,intent(in)     :: debug                   ! debug flag
  character(*),intent(in)     :: calling_prog            ! name of calling program
  integer*4   ,intent(in)     :: icomp                   ! component (X,Y or Z accelerations) passed in (used to name output file)
  integer*4   ,intent(in)     :: nvals                   ! dimensioning of rows in acc_data for one component
  integer*4   ,intent(in)     :: ncomp                   ! number of components requested from the EMD filter
  integer*4   ,intent(in)     :: acc_secs                ! number of seconds from start to end of acc data
  real(kind=8),intent(in)     :: acc_data_in(nvals)      ! storage of actual linear acceleration data for one component
  real(kind=8),intent(out)    :: acc_decomp(ncomp,nvals) ! decomposition of the input time series
  integer*4   ,intent(out)    :: end_comp                ! the last component to include in the reconstruction

! local variables
  real(kind=8)                :: mean
  real(kind=8),allocatable    :: acc_data_tmp(:)
  integer*4                   :: iepoch,i,j
  character*11                :: emd_outfile
  integer*4                   :: stop_comp
  real(kind=8)                :: max_tol                 ! peak-to-peak range that the component must exceed to be included in recomposition
  integer*4                   :: n_inflexions            ! minimum number of inflexions for which we keep a component (try ~10 per 24-hours)
  integer*4                   :: n_crossings             ! minimum number of zero crossings for which we keep a component (try ~10 per 24-hours)
  real(kind=8)                :: pi

! variables for the second run of the EMD
  real(kind=8),allocatable :: acc_decomp1(:,:)
  real(kind=8),allocatable :: acc_data_tmp2(:)
  call status_update("STATUS",calling_prog,"ACC_EMD1","","Using ACC_EMD1",0)
  pi = 4.d0*datan(1.d0)

  if(debug)print*,"inside acc_EMD1"
! allocate temporary storate
  allocate(acc_data_tmp(nvals))
  allocate(acc_data_tmp2(nvals))
  allocate(acc_decomp1(ncomp,nvals))

! calculate and remove a mean value so that de-meaned data are sent to the EMD filter
  mean = sum(acc_data_in,nvals)/dble(nvals)
  acc_data_tmp = acc_data_in - mean
  if(debug)print*,"mean value is",mean

! PT180820: try adding in a once-per-tworev signal to see whether the mode mixing of the once-per-rev signal
!           can be mitigated. The assumption is that the EMD will identify the 1//tworev and then the code below
!           will deem it - and, therefore, all others of lower wavelength - to be excluded from the reconstruction
!  do iepoch = 1,nvals
!    acc_data_tmp(iepoch) = acc_data_tmp(iepoch) + 10.d-9*dsin(dble(iepoch)*1.d0/5400.d0 * (0.88d0*pi) )
!  enddo

! call the EMD filter
  acc_decomp = 0.d0
  if(debug)print*,"Calling EMD filter",ncomp,nvals
  call emd(acc_data_tmp,nvals,ncomp,acc_decomp)
  if(debug)print*,"Back from EMD filter"


! starting from the highest frequency component, determine whether there is a point of inflexion within
! one revolution. If so, the component contains high-frequency signal and we want to keep it. If not, it
! is a long-wavelength component and it - and all others below it - should be excluded from the reconstruction
  stop_comp = 0
  i = 1

! set the tolerance based on whether it is X, Y or Z accelerometer obs
  if(icomp == 1)then
    max_tol = 100.e-9
  else if (icomp == 2)then
    max_tol = 15.e-9
  else
    max_tol = 5.e-9
  endif

! how many inflexions must occur before we keep a component?
  n_inflexions = nint(15*dble(acc_secs)/86400.d0)  ! allow 10 per day
  do while (stop_comp == 0 .and. i <= ncomp)
! PT180926: replace this with a routine that counts zero crossings of the data instead of inflexion points
!    call acc_inflexions(debug,calling_prog,i,nvals,max_tol,n_inflexions,stop_comp,acc_decomp(i,:))
    n_crossings = 10
    call acc_zero_crossings(debug,calling_prog,i,nvals,n_crossings,stop_comp,acc_decomp(i,:))
    if(stop_comp == 0)i = i + 1
  enddo

! save and pass back the value of the last component to include in the reconstruction
  end_comp = i

! debug: write out to a file all the components of the EMD decomposition at each epoch
  write(emd_outfile,'(a,i1)')"emd.decomp",icomp
  open(1000,file=emd_outfile,status='unknown')
  do iepoch = 1,nvals
    write(1000,*)iepoch,(acc_decomp(j,iepoch),j=1,ncomp)
  enddo
  close(1000)


! add back the mean value
  acc_decomp(1,:) = acc_decomp(1,:) + mean

  deallocate(acc_data_tmp)
  deallocate(acc_data_tmp2)
  deallocate(acc_decomp1)

  return

  end subroutine acc_EMD1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine acc_EMD2(debug,calling_prog,icomp,nvals,ncomp,acc_secs,acc_data_in,acc_decomp,end_comp)

! apply the EMD filter (used in gracefit to remove high frequency noise from KBR)
! to remove the low-frequency components from the ACC1B data
!
! P. Tregoning
! 6 August 2018
!
! PT180817: modify the EMD process to try to reduce mode mixing that sees some of the 1/rev signal 
!           appearing in a lower frequency mode. We don't want that!

  implicit none

! passed variables
  logical     ,intent(in)     :: debug                   ! debug flag
  character(*),intent(in)     :: calling_prog            ! name of calling program
  integer*4   ,intent(in)     :: icomp                   ! component (X,Y or Z accelerations) passed in (used to name output file)
  integer*4   ,intent(in)     :: nvals                   ! dimensioning of rows in acc_data for one component
  integer*4   ,intent(in)     :: ncomp                   ! number of components requested from the EMD filter
  integer*4   ,intent(in)     :: acc_secs                ! number of seconds from start to end of acc data
  real(kind=8),intent(in)     :: acc_data_in(nvals)      ! storage of actual linear acceleration data for one component
  real(kind=8),intent(out)    :: acc_decomp(ncomp,nvals) ! decomposition of the input time series
  integer*4   ,intent(out)    :: end_comp                ! the last component to include in the reconstruction

! local variables
  real(kind=8)                :: mean
  real(kind=8),allocatable    :: acc_data_tmp(:)
  integer*4                   :: iepoch,i,j
  character*11                :: emd_outfile
  integer*4                   :: stop_comp
  real(kind=8)                :: max_tol                 ! peak-to-peak range that the component must exceed to be included in recomposition
  integer*4                   :: n_inflexions            ! minimum number of inflexions for which we keep a component (try ~10 per 24-hours)
  integer*4                   :: n_crossings             ! criterion for the number of zero crossings of the data
  real(kind=8)                :: pi
  integer*4                   :: nvals_tmp
  logical                     :: use_data

! variables for the second run of the EMD
  real(kind=8),allocatable :: acc_decomp1(:,:)
  real(kind=8),allocatable :: acc_data_tmp2(:)
  call status_update("STATUS",calling_prog,"ACC_EMD2","","Using ACC_EMD2",0)
  pi = 4.d0*datan(1.d0)

  if(debug)print*,"inside acc_EMD2"
! PT181003: limit the end of the data to be filtered in the case that there are records with zero values
  use_data = .true.
  nvals_tmp = 0
  do while (use_data .and. nvals_tmp < nvals)
!print*,nvals_tmp+1,acc_data_in(nvals_tmp+1)
    if(acc_data_in(nvals_tmp+1) == 0.d0)then
print*,'zero acc obs at epoch',nvals_tmp+1
      use_data = .false.
    else
      nvals_tmp = nvals_tmp + 1
    endif
  enddo

! allocate temporary storate
  allocate(acc_data_tmp(nvals_tmp))
  allocate(acc_data_tmp2(nvals_tmp))
  allocate(acc_decomp1(ncomp,nvals_tmp))
  acc_decomp1 = 0.d0


! calculate and remove a mean value so that de-meaned data are sent to the EMD filter
  mean = sum(acc_data_in(1:nvals_tmp))/dble(nvals_tmp)
  acc_data_tmp(1:nvals_tmp) = acc_data_in(1:nvals_tmp) - mean
  if(debug)print*,"mean value is",mean


! call the EMD filter
  acc_decomp = 0.d0
  if(debug)print*,"Calling EMD filter",ncomp,nvals
  call emd_v2(acc_data_tmp(1:nvals_tmp),nvals_tmp,ncomp,acc_decomp1)
  if(debug)print*,"Back from EMD filter"


! starting from the highest frequency component, determine whether there is a point of inflexion within
! one revolution. If so, the component contains high-frequency signal and we want to keep it. If not, it
! is a long-wavelength component and it - and all others below it - should be excluded from the reconstruction
  stop_comp = 0
  i = 1

! set the tolerance based on whether it is X, Y or Z accelerometer obs
  if(icomp == 1)then
    max_tol = 100.e-9
  else if (icomp == 2)then
    max_tol = 7.e-9
  else
    max_tol = 5.e-9
  endif

! how many inflexions must occur before we keep a component?
  n_inflexions = nint(15*dble(acc_secs)/86400.d0)  ! allow 10 per day
  do while (stop_comp == 0 .and. i <= ncomp)
! PT180926: replace this with a routine that counts zero crossings of the data instead of inflexion points
!    call acc_inflexions(debug,calling_prog,i,nvals,max_tol,n_inflexions,stop_comp,acc_decomp(i,:))
    n_crossings = 17
    call acc_zero_crossings(debug,calling_prog,i,nvals_tmp,n_crossings,stop_comp,acc_decomp1(i,:))
    if(stop_comp == 0)i = i + 1
  enddo

! save and pass back the value of the last component to include in the reconstruction
  end_comp = i

! debug: write out to a file all the components of the EMD decomposition at each epoch
  write(emd_outfile,'(a,i1)')"emd.decomp",icomp
  open(1000,file=emd_outfile,status='unknown')
  do iepoch = 1,nvals_tmp
    write(1000,*)iepoch,(acc_decomp1(j,iepoch),j=1,ncomp)
  enddo
  close(1000)


! add back the mean value and transfer from acc_decomp1 to acc_decomp
  do iepoch = 1,nvals_tmp
    acc_decomp(1,iepoch)    = acc_decomp1(1,iepoch) + mean
    acc_decomp(2:13,iepoch) = acc_decomp1(2:13,iepoch)
  enddo

  deallocate(acc_data_tmp)
  deallocate(acc_data_tmp2)
  deallocate(acc_decomp1)

  return

  end subroutine acc_EMD2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine acc_EMD3(debug,calling_prog,icomp,nvals,ncomp,acc_secs,acc_data_in,acc_decomp,end_comp)

! apply the EMD filter (used in gracefit to remove high frequency noise from KBR)
! to remove the low-frequency components from the ACC1B data
!
! P. Tregoning
! 6 August 2018
!
! PT180817: modify the EMD process to try to reduce mode mixing that sees some of the 1/rev signal 
!           appearing in a lower frequency mode. We don't want that!
! RM181203: subroutine calls emd_v3

  implicit none

! passed variables
  logical     ,intent(in)     :: debug                   ! debug flag
  character(*),intent(in)     :: calling_prog            ! name of calling program
  integer*4   ,intent(in)     :: icomp                   ! component (X,Y or Z accelerations) passed in (used to name output file)
  integer*4   ,intent(in)     :: nvals                   ! dimensioning of rows in acc_data for one component
  integer*4   ,intent(in)     :: ncomp                   ! number of components requested from the EMD filter
  integer*4   ,intent(in)     :: acc_secs                ! number of seconds from start to end of acc data
  real(kind=8),intent(in)     :: acc_data_in(nvals)      ! storage of actual linear acceleration data for one component
  real(kind=8),intent(out)    :: acc_decomp(ncomp,nvals) ! decomposition of the input time series
  integer*4   ,intent(out)    :: end_comp                ! the last component to include in the reconstruction

! local variables
  real(kind=8)                :: mean
  real(kind=8),allocatable    :: acc_data_tmp(:)
  integer*4                   :: iepoch,i,j
  character*11                :: emd_outfile
  integer*4                   :: stop_comp
  real(kind=8)                :: max_tol                 ! peak-to-peak range that the component must exceed to be included in recomposition
  integer*4                   :: n_inflexions            ! minimum number of inflexions for which we keep a component (try ~10 per 24-hours)
  integer*4                   :: n_crossings             ! criterion for the number of zero crossings of the data
  real(kind=8)                :: pi
  integer*4                   :: nvals_tmp
  logical                     :: use_data

! variables for the second run of the EMD
  real(kind=8),allocatable :: acc_decomp1(:,:)
  real(kind=8),allocatable :: acc_data_tmp2(:)
  call status_update("STATUS",calling_prog,"ACC_EMD3","","Using ACC_EMD3",0)
  pi = 4.d0*datan(1.d0)

  if(debug)print*,"inside acc_EMD3"
! PT181003: limit the end of the data to be filtered in the case that there are records with zero values
  use_data = .true.
  nvals_tmp = 0
  do while (use_data .and. nvals_tmp < nvals)
!print*,nvals_tmp+1,acc_data_in(nvals_tmp+1)
    if(acc_data_in(nvals_tmp+1) == 0.d0)then
print*,'zero acc obs at epoch',nvals_tmp+1
      use_data = .false.
    else
      nvals_tmp = nvals_tmp + 1
    endif
  enddo

! allocate temporary storate
  allocate(acc_data_tmp(nvals_tmp))
  allocate(acc_data_tmp2(nvals_tmp))
  allocate(acc_decomp1(ncomp,nvals_tmp))
  acc_decomp1 = 0.d0


! calculate and remove a mean value so that de-meaned data are sent to the EMD filter
  mean = sum(acc_data_in(1:nvals_tmp))/dble(nvals_tmp)
  acc_data_tmp(1:nvals_tmp) = acc_data_in(1:nvals_tmp) - mean
  if(debug)print*,"mean value is",mean


! call the EMD filter
  acc_decomp = 0.d0
  if(debug)print*,"Calling EMD filter",ncomp,nvals
  call emd_v3(acc_data_tmp(1:nvals_tmp),nvals_tmp,ncomp,acc_decomp1)
  if(debug)print*,"Back from EMD filter"


! starting from the highest frequency component, determine whether there is a point of inflexion within
! one revolution. If so, the component contains high-frequency signal and we want to keep it. If not, it
! is a long-wavelength component and it - and all others below it - should be excluded from the reconstruction
  stop_comp = 0
  i = 1

! set the tolerance based on whether it is X, Y or Z accelerometer obs
  if(icomp == 1)then
    max_tol = 100.e-9
  else if (icomp == 2)then
    max_tol = 7.e-9
  else
    max_tol = 5.e-9
  endif

! how many inflexions must occur before we keep a component?
  n_inflexions = nint(15*dble(acc_secs)/86400.d0)  ! allow 10 per day
  do while (stop_comp == 0 .and. i <= ncomp)
! PT180926: replace this with a routine that counts zero crossings of the data instead of inflexion points
!    call acc_inflexions(debug,calling_prog,i,nvals,max_tol,n_inflexions,stop_comp,acc_decomp(i,:))
    n_crossings = 3
    call acc_zero_crossings(debug,calling_prog,i,nvals_tmp,n_crossings,stop_comp,acc_decomp1(i,:))
    if(stop_comp == 0)i = i + 1
  enddo

! save and pass back the value of the last component to include in the reconstruction
  end_comp = i

! debug: write out to a file all the components of the EMD decomposition at each epoch
  write(emd_outfile,'(a,i1)')"emd.decomp",icomp
  open(1000,file=emd_outfile,status='unknown')
  do iepoch = 1,nvals_tmp
    write(1000,*)iepoch,(acc_decomp1(j,iepoch),j=1,ncomp)
  enddo
  close(1000)


! add back the mean value and transfer from acc_decomp1 to acc_decomp
  do iepoch = 1,nvals_tmp
    acc_decomp(1,iepoch)    = acc_decomp1(1,iepoch) + mean
    acc_decomp(2:13,iepoch) = acc_decomp1(2:13,iepoch)
  enddo

  deallocate(acc_data_tmp)
  deallocate(acc_data_tmp2)
  deallocate(acc_decomp1)

  return

  end subroutine acc_EMD3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine acc_fit_exponential(debug,calling_prog,nvals,acc_data,thr_obs,start_ep,end_ep,params)

! subroutine to fit a model y = Ae^(Bx+C) + De^(Ex^2 + F) to the data. 
! Not sure how much of the model we then want to remove .... we need to leave the "bias" intact!
!
! P. Tregoning
! 6 July 2018
!

  implicit none

! passed variables
  logical     ,intent(in)      :: debug             ! debug flag
  character(*),intent(in)      :: calling_prog      ! name of calling program
  integer*4   ,intent(in)      :: nvals             ! dimensioning of rows in acc_data for one component
  real(kind=8),intent(inout)   :: acc_data(nvals)    ! storage of actual linear acceleration data for one component
  real(kind=8),intent(out)     :: thr_obs(nvals)    ! thrust flags indicating whether acc_obs are affected by thrusts
  integer*4   ,intent(in)      :: start_ep,end_ep   ! start and end of epochs to use in the fit of the model
  real(kind=8),intent(inout)   :: params(4)         ! actual estimated parameter values

! local variables
  integer*4,parameter      :: nparams=4
  real(kind=8),allocatable :: Amat(:,:),Bmat(:,:),At(:,:),AtA(:,:),AtB(:,:),VCV(:,:),soln(:,:)
  real(kind=8)             :: apr(nparams)            ! array to store the iterated a priori parameter values
  real(kind=8)             :: dt
  integer*4                :: i,iter,next_row,next_col,iepoch ! various counters
  integer*4                :: n_affected              ! number of observations affected by thrusts
  real(kind=8)             :: meanval                 ! mean value of the input obs unaffected by thrusts
  integer*4                :: n_good                  ! number of observations unaffected by thrusts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! first, allocate the least squares arrays
  allocate(Amat(nvals,nparams))
  allocate(Bmat(nvals,1))           
  allocate(At(nparams,nvals))
  allocate(AtA(nparams,nparams))
  allocate(VCV(nparams,nparams))
  allocate(AtB(nparams,1))
  allocate(soln(nparams,1))

! set the a priori values
  apr(1:4) = 0.5d0
  apr(3:4) = 1.0d0
  meanval = 0.d0
  n_good = 0
! calculate a mean value of the obs
  do iepoch = start_ep, end_ep
    if(thr_obs(iepoch) == 1.d0)then               ! epoch is not affected by thrusts
      n_good = n_good + 1
      meanval = meanval + acc_data(iepoch)
    endif
  enddo
  if(n_good > 0)meanval = meanval/dble(n_good)

! iterate the least squares 10 times
  do iter = 1,20
    next_row = 0
    next_col = 0
    Amat = 0.d0
    Bmat = 0.d0
    n_affected = 0
      do iepoch = start_ep,end_ep
        if(thr_obs(iepoch) == 1.d0)then               ! epoch is not affected by thrusts
          dt = dble(iepoch)/86400.d0
          Amat(iepoch,1) = exp(apr(2)*dt)       
          Amat(iepoch,2) = dt*apr(1)*exp( apr(2)*dt)     
!          Amat(iepoch,3) = apr(1)*exp( apr(2)*dt+apr(3))

          Amat(iepoch,3) = exp(apr(4)*dt**2)       
          Amat(iepoch,4) = dt**2*apr(3)*exp( apr(4)*dt**2)     
!          Amat(iepoch,6) = apr(4)*exp( apr(5)*dt**2+apr(6))


          Bmat(iepoch,1) = acc_data(iepoch) - (apr(1)*exp( apr(2)*dt) + apr(3)*exp( apr(4)*dt**2) )  !- meanval 
!print*,'iter,iepoch,Amat(iepoch,:),Bmat(iepoch,1)',iter,iepoch,Amat(iepoch,:),Bmat(iepoch,1)
!print*,'start_ep,end_ep,iter,iepoch,dt,Bmat(iepoch,1)',start_ep,end_ep,iter,iepoch,dt,Bmat(iepoch,1)
        else
!print*,'epoch ',iepoch,' is thrust-affected'
          n_affected = n_affected + 1
        endif
    enddo

! Least squares solution
    ! build normal equations
    call transp(Amat,At,nvals,nparams)
    call matmult(At,Amat,AtA,nparams,nvals,nparams)
!call printmat(Ata,nparams,nparams,"AtA")
    call invert(AtA,VCV,nparams)
!call printmat(VCV,nparams,nparams,"AtA")

    ! RHS
    call matmult(At,Bmat,AtB,nparams,nvals,1)

    ! solution
    call matmult(VCV,AtB,soln,nparams,nparams,1)

! update the parameter values
    next_row = 0
    do i=1,nparams
      apr(i) = apr(i) + soln(i,1)
    enddo
    if(debug)print*,iter,apr,soln(:,1)
  enddo

! remove the modelled rate and acceleration from the observations, but leave the offset (the bias) uncorrected
  do iepoch = 1,nvals
    dt = dble(iepoch)/86400.d0
!print*,iepoch,acc_obs(iepoch,:)," uncorrected"
      acc_data(iepoch) = acc_data(iepoch) - (apr(1)*exp( apr(2)*dt) + apr(3)*exp( apr(4)*dt**2) )
  enddo  

! return the final parameter values
  params = apr

  return
  end subroutine acc_fit_exponential
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine acc_fit_quadratic(debug,calling_prog,nvals,acc_data,thr_obs,start_ep,end_ep,params)

! subroutine to form up the A and B matrices for the least squares inversion of 
! a quadratic fit to the data
!
! P. Tregoning
! 31 May 2018
!
! PT180702: modified to work on only one component of accelerometer data, rather than all three at once.

  implicit none

! passed variables
  logical     ,intent(in)      :: debug
  character(*),intent(in)      :: calling_prog      ! name of calling program
  integer*4   ,intent(in)      :: nvals             ! dimensioning of rows in acc_data for one component
  real(kind=8),intent(inout)   :: acc_data(nvals)    ! storage of actual linear acceleration data for one component
  real(kind=8),intent(in)      :: thr_obs(nvals)    ! thrust flags indicating whether acc_obs are affected by thrusts
  integer*4   ,intent(in)      :: start_ep,end_ep   ! start and end of epochs to use in the fit of the model
  real(kind=8),intent(inout)   :: params(3)         ! actual estimated parameter values

! local variables
  integer*4,parameter      :: nparams=3
  real(kind=8),allocatable :: Amat(:,:),Bmat(:,:),At(:,:),AtA(:,:),AtB(:,:),VCV(:,:),soln(:,:)
  real(kind=8)             :: apr(3)               ! array to store the iterated a priori parameter values
  real(kind=8)             :: dt
  integer*4                :: i,iter,next_row,next_col,iepoch ! various counters
  integer*4                :: n_affected              ! number of observations affected by thrusts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! first, allocate the least squares arrays
  allocate(Amat(nvals,nparams))
  allocate(Bmat(nvals,1))           
  allocate(At(nparams,nvals))
  allocate(AtA(nparams,nparams))
  allocate(VCV(nparams,nparams))
  allocate(AtB(nparams,1))
  allocate(soln(nparams,1))

! set the a priori values
  apr = 0.d0

! iterate the least squares 10 times
  do iter = 1,4
    next_row = 0
    next_col = 0
    Amat = 0.d0
    Bmat = 0.d0
    n_affected = 0
      do iepoch = start_ep,end_ep
        if(thr_obs(iepoch) == 1.d0)then               ! epoch is not affected by thrusts
! PT190424: subtract the first epoch to reduce correlations.... has no effect if start_epoch = 1
          dt = dble(iepoch-start_ep+1)/86400.d0
          Amat(iepoch,1) = dt**2    ! quadratic term 
          Amat(iepoch,2) = dt       ! linear term
          Amat(iepoch,3) = 1.d0     ! offset term

          Bmat(iepoch,1) = acc_data(iepoch) - (apr(1)*dt**2 + apr(2)*dt + apr(3))
        else
          n_affected = n_affected + 1
        endif
    enddo

! Least squares solution
    ! build normal equations
    call transp(Amat,At,nvals,nparams)
    call matmult(At,Amat,AtA,nparams,nvals,nparams)
    if(debug)call printmat_R8(VCV,nparams,nparams,"acc_fit_quad:AtA ")
    call invert(AtA,VCV,nparams)
    if(debug)call printmat_R8(VCV,nparams,nparams,"acc_fit_quad:VCV ")

    ! RHS
    call matmult(At,Bmat,AtB,nparams,nvals,1)

    ! solution
    call matmult(VCV,AtB,soln,nparams,nparams,1)

! update the parameter values
    next_row = 0
    do i=1,3
      apr(i) = apr(i) + soln(i,1)
    enddo
    if(debug)print*,"acc_fit_quadratic,iter,apr,soln:",iter,apr,soln(:,1)
  enddo

! remove the modelled rate and acceleration from the observations, but leave the offset (the bias) uncorrected
  do iepoch = 1,nvals
    dt = dble(iepoch-start_ep + 1)/86400.d0
    acc_data(iepoch) = acc_data(iepoch) - (apr(1)*dt**2 + apr(2)*dt )  ! subtract the accel and rate only
  enddo  

! return the final parameter values
  params = apr

! deallocate arrays
  deallocate(Amat)
  deallocate(Bmat)           
  deallocate(At)
  deallocate(AtA)
  deallocate(VCV)
  deallocate(AtB)
  deallocate(soln)

  return
  end subroutine acc_fit_quadratic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine acc_read_data(luin,calling_prog,acc_file,mission,nvals,ncols,acc_span,acc_obs )

! subroutine to read all the linear accelerations of an ACC1B file. 
!
! P. Tregoning
! 31 May 2018
!
! PT180810: made variable the number of observations per epoch read from the file

  implicit none

! passed variables
  integer*4   ,intent(in)    :: luin              ! unit number of file
  character(*),intent(in)    :: calling_prog      ! name of calling program
  character(*),intent(in)    :: acc_file          ! name of accelerometer file
  integer*4   ,intent(inout) :: mission           ! 0: GRACE, 1: GRACE FO, 2: GRACE II  
  integer*4   ,intent(in)    :: acc_span(2)       ! start/stop epochs (in gracesec) of accelerometer data
  integer*4   ,intent(in)    :: nvals             ! dimensioning of rows in acc_data array
  integer*4   ,intent(in)    :: ncols             ! number of columns of observations to read from the ACC1B file
  real(kind=8),intent(out)   :: acc_obs(nvals,ncols+2)  ! gracesec, XYZ accel obs, accel_flag

! local variables
  integer*4     :: ioerr,irow,gracesec
  character*1   :: sat
  character*250 :: line,message
  real(kind=8),allocatable  :: tmp_acc_obs(:)
  allocate(tmp_acc_obs(ncols))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the accelerometer data. The header should already have been read, 
! so we just need to read in the data
  ioerr = 0
  gracesec = 0
  do while (ioerr == 0 .and. gracesec <= acc_span(2))
    read(luin,*,iostat=ioerr,end=1000)gracesec,sat,tmp_acc_obs(1:ncols)
    if(ioerr == 0 .and. gracesec >= acc_span(1) .and. gracesec <= acc_span(2))then
      ! work out which row it should go into, being the offset from the start epoch
        irow = 1 + (gracesec - acc_span(1) )  
      ! store the data
      ! PT180917: only store if less than the dimensioning of the array
        if(irow <= nvals)then
          acc_obs(irow,1)   = gracesec
          acc_obs(irow,2:ncols+1) = tmp_acc_obs(1:ncols)
        endif
    endif
  enddo

1000 continue
  write(message,'(a,i7,a)')"  read",irow," accelerometer obs from ACC1B file"
  call status_update('STATUS',calling_prog,'acc_read_data',acc_file,message,0)

  return
  end subroutine acc_read_data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine acc_read_hdr(luin,calling_prog,acc_file,mission,acc_span )

! subroutine to read through the header information of an ACC1B file
!
! P. Tregoning
! 31 May 2018

  implicit none

! passed variables
  integer*4   ,intent(in)    :: luin         ! unit number of file
  character(*),intent(in)    :: calling_prog ! name of calling program
  character(*),intent(in)    :: acc_file     ! name of accelerometer file
  integer*4   ,intent(inout) :: mission      ! 0: GRACE, 1: GRACE FO, 2: GRACE II  
  integer*4   ,intent(out)   :: acc_span(2)  ! start/stop epochs (in gracesec) of accelerometer data

! local variables
  real(kind=8)  :: seconds
  integer*4     :: ioerr
  character*250 :: line,message

! check whether the mission variable has been set
  if(mission < 0) then

  ! not set. Determine the mission (GRACE, GRACE FO, GRACE II) from the first line of the header
    read(luin,'(a)')line
    if(line(1:7) == "PRODUCE")then
      mission = 0    ! it is a GRACE accelerometer file
      call status_update('STATUS',calling_prog,'acc_read_hdr',' ','GRACE format found ',0)
    elseif(line(1:7) == "header:")then
      mission = 1    ! it is a GRACE FO accelerometer file
      call status_update('STATUS',calling_prog,'acc_read_hdr',' ','GRACE FO format found',0)
    endif
    backspace(luin)
  endif

! now, read the header information.
  line = ""
  ioerr = 0
  if(mission == 0)then
    do while (line(1:13) /= "END OF HEADER" .and. ioerr == 0)
      read(luin,'(a)',iostat=ioerr)line
      if(line(1:31) == "TIME FIRST OBS(SEC PAST EPOCH):")then
         read(line(32:42),*) seconds              ! HM190608 read seconds as real in case it is <10^8 
         acc_span(1)=int(seconds)
      elseif(line(1:31) == "TIME LAST OBS(SEC PAST EPOCH) :")then
         read(line(32:42),*) seconds
         acc_span(2)=int(seconds)
      endif
    enddo

  else if (mission == 1)then
    do while (line(1:20) /= "# End of YAML header" .and. ioerr == 0)
      read(luin,'(a)',iostat=ioerr)line
      if(line(1:25)     == "    start_time_epoch_secs")then
        read(line(27:36),*)acc_span(1)
      elseif(line(1:25) == "    stop_time_epoch_secs:")then
        read(line(26:35),*)acc_span(2)
      endif
    enddo    

  endif

  write(message,'(a,i12,a,i11)')"Span of ACC1B data: ",acc_span(1),"   to",acc_span(2)
  call status_update('STATUS',calling_prog,'acc_read_hdr',acc_file,message,0)

  return
  end subroutine acc_read_hdr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









