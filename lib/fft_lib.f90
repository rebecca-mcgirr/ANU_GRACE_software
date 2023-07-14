!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  This library contains the functions necessary to perform a Fast Fourier 
!!  Transform on real single precision 1D data and filter the frequencies using 
!!  a raised cosine low/high pass filter. This library contains a "lite" version  
!!  of FFTPACK5.1 containing routines needed for real FFT synthesis and analysis.
!!  Originally written in FORTRAN77 and translated to FORTRAN90. 
!! 
!!  compute_fft      : Fourier analysis/synthesis and filtering
!!  window_function  : Returns window/taper coefficients
!!  cos_filt         : returns the sample frequencies 
!!  rfft1b           : initialisation routine for rfftb1
!!  rfft1f           : initialisation routine for rfftf1
!!  rfft1i           : initialisation routine for rffti1
!!  rfftb1           : backward transform of a real coefficient array
!!  rfftf1           : forward transform of a real periodic sequence
!!  rffti1           : initialise rfftf1 and rfftb1
!!  r1f2kb           : auxilliary function
!!  r1f2kf           : auxilliary function
!!  r1f3kb           : auxilliary function
!!  r1f3kf           : auxilliary function
!!  r1f4kb           : auxilliary function
!!  r1f4kf           : auxilliary function
!!  r1f5kb           : auxilliary function
!!  r1f5kf           : auxilliary function
!!  r1fgkb           : auxilliary function
!!  r1fgkf           : auxilliary function
!!  xerfft           : error handler
!!
!!  R. McGirr
!!  19 June 2019
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_fft(debug,calling_prog,nvals,nfft,obs_in,obs_out,ftype,flo,fhi,fspec,pspec_raw,pspec_filt,filt)

! Perform Fourier analysis on a time series, then filter low frequencies with a
! cosine high/low pass filter, finally, perform Fourier synthesis to transform 
! frequencies back to a time series
!
! note: nfft should be set to a value equal to or greater than nvals, the larger
! the number the 'cleaner' the spectrum will appear, but the calculations will be 
! slower. The FFT is most efficient if nfft is a power of two. 
!
! R. McGirr
! 26 June 2019

  implicit none

! passed variables
  logical     ,intent(in)     :: debug                   ! debug flag
  character(*),intent(in)     :: calling_prog            ! name of calling program
  integer*4   ,intent(in)     :: nvals                   ! lengh of input data
  integer*4   ,intent(in)     :: nfft                    ! see note in subroutine header
  real(kind=4),intent(in)     :: flo                     ! low frequency cutoff
  real(kind=4),intent(in)     :: fhi                     ! high frequency cutoff
  real(kind=8),intent(in)     :: obs_in(nvals)           ! input vector time series
  real(kind=8),intent(out)    :: obs_out(nvals)          ! input vector time series
  real(kind=4),intent(out)    :: pspec_raw(nfft/2 + 1)   ! contains power spectrum
  real(kind=4),intent(out)    :: pspec_filt(nfft/2 + 1)  ! contains filtered power spectrum
  real(kind=4),intent(out)    :: fspec(nfft/2 + 1)       ! contains frequencies
  real(kind=4),intent(out)    :: filt(nfft/2 + 1)
  character*10,intent(in)     :: ftype                   ! set to 'high' or 'low'

! local variables
  integer*4                   :: lensav
  integer*4                   :: lenwrk
  integer*4                   :: ier
  integer*4                   :: inc
  integer*4                   :: lenr
  integer*4                   :: i,j
  real(kind=4)                :: fft_tmp(nfft)
  real(kind=4)                :: fft_filt(nfft)
  real(kind=4),allocatable    :: work(:)
  real(kind=4),allocatable    :: wsave(:) 
  character*100               :: message

  call status_update("STATUS",calling_prog,"COMPUTE_FFT","","Using COMPUTE_FFT",0)
  if(debug)print*,"obs_in",obs_in(1:10)

! calculate and remove a mean value so that de-meaned data are sent to the forward rfft
  do i = 1,nvals
    fft_tmp(i) = obs_in(i)
  enddo

! set length of work arrays for use in rfft1f and rfft1b
  lensav = nfft + int(log(real(nfft))/log( 2.0E+00 )) + 4
  lenwrk = nfft
  if(debug)print*,"lensav is",lensav,"lenwrk is",lenwrk

! allocate work arrays
  allocate(work(lenwrk))
  allocate(wsave(lensav))

! call rfft initialisation routine
  call rfft1i(nfft, wsave, lensav, ier)
  if(debug)print*,"wsave",wsave(1:10)

! compute the FFT coefficients via the rfft forward routine
  inc = 1
  lenr = nfft
  call rfft1f(nfft, inc, fft_tmp, lenr, wsave, lensav, work, lenwrk, ier) 
  if(debug)print*,"fft_tmp",fft_tmp(1:10)

! fill in frequency array (cycles/second)
  if(debug)print*,"length of sample frequency array",nfft/2 + 1
  do i = 1,nfft/2 + 1
    fspec(i) = (i-1) * 1./(nfft*inc)
  enddo
  if(debug)print*,"fspec",fspec(1:10)

! convert FFT coefficients to power spectra for plotting
  j = 0
  pspec_raw(1) = fft_tmp(1)**2
  do i = 2,nfft/2
    j = j + 2
    pspec_raw(i) = fft_tmp(j)**2 + fft_tmp(j+1)**2
  enddo
  pspec_raw(nfft/2 + 1) = fft_tmp(nfft)**2

! filter frequencies via the high pass cosine filter
  if(debug)print*,"ftype ",ftype,"flo",flo,"fhi",fhi
  call cos_filt(debug,calling_prog,ftype,fspec,nfft,flo,fhi,filt)
  if(debug)print*,"filt",filt(1:10)

! multiply FFT coefficients by filter
  fft_filt(1) = filt(1) * fft_tmp(1)
  do i = 2,nfft/2
    fft_filt(2*i-2) = filt(i) * fft_tmp(2*i-2)
    fft_filt(2*i-1) = filt(i) * fft_tmp(2*i-1)
  enddo
  fft_filt(nfft) = filt(nfft/2 + 1) * fft_tmp(nfft)

! convert FFT coefficients to power spectra for plotting
  j = 0
  pspec_filt(1) = fft_filt(1)**2
  do i = 2,nfft/2
    j = j + 2
    pspec_filt(i) = fft_filt(j)**2 + fft_filt(j+1)**2
  enddo
  pspec_filt(nfft/2 + 1) = fft_filt(nfft)**2

! compute inverse FFT coefficients via the rfft backward routine
  call rfft1b(nfft, inc, fft_filt, lenr, wsave, lensav, work, lenwrk, ier)

! deallocate arrays
  deallocate(work)
  deallocate(wsave)

  obs_out = fft_filt(:nvals)
  if(debug)print*,"obs_out",obs_out(1:10)

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine window_function(calling_prog,n,wtype,w)

! Returns an array of weighting coefficients given length of array and window
! type ('hann' or 'uniform'). Windows or tapering functions are used to convolve
! signals to smooth discontinuities and reduce spectral leaage when performing 
! a forward FFT.  
!
! The Hann window is a taper formed by using a weighted cosine. The returned 
! window has a maximum value normalized to one. The Hann window, in general,
! is appropriate for most cases due to its good frequency resolution and 
! reduced spectral leakage. 
!
! The uniform window (rectangular/no window) simply returns an array filled with
! ones. No window is appropriate for signals with an integer number of revolutions
! or if the signal spectrum is relatively flat or broadband in frquency content.
!
! Other window/tapering functions that may be of interest include Hamming,
! Blackman, Tukey, Blackman and Gaussian
!
! R. McGirr
! 15 July 2019

  implicit none

! passed variables
  character(*),intent(in)     :: calling_prog          ! name of calling program
  integer*4   ,intent(in)     :: n                     ! length of vector
  character*8 ,intent(in)     :: wtype                 ! set to 'hann' or 'uniform'
  real(kind=4),intent(out)    :: w(n)        	       ! ouput window vector


! local variables
  integer*4                   :: i
  real(kind=4)                :: alpha = 0.5       	       
  real(kind=4)                :: a(2)        	        
  real(kind=4)                :: fac(n)        	        
  character*250               :: message
  real(kind=8),parameter      :: pi = 4.d0*atan(1.d0)  ! approximation of pi


  if( wtype .eq. "hann" )then
    a = (/ alpha, 1.-alpha /)
    do i = 0,n-1
      fac(i+1) = -pi + i*(pi+pi)/(n-1) 
    enddo
    do i = 0,1
      w = w + a(i+1) * cos(i*fac)
    enddo
    write(message,'(a,i6,a)')" Applying Hann window to input data of length",n," epochs"
  else if( wtype .eq. "uniform" )then
    w = 1.0
    write(message,'(a,i6,a)')" Applying uniform window to input data of length",n," epochs"
  else
    write(message,'(a,i6,a)')" Window type",wtype," not recognised"
  endif
  call status_update('STATUS','LIB',calling_prog,' ',message,0)

end subroutine window_function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cos_filt(debug,calling_prog,ftype,fs,n,flo,fhi,filt)

! Built for use in conjunction with rfft1f
! Raised cosine low pass filter
!          1 if fs < flo
!   filt = 0.5 * (cos(pi*(f - flo)/(flo - fhi)) + 1) if flo < fs < fhi
!          0 if fs > fhi
!
! Raised cosine high pass filter
!          0 if fs < flo
!   filt = 0.5 * (-cos(pi*(f - flo)/(flo - fhi)) + 1) if flo < fs < fhi
!          1 if fs > fhi
!
! R. McGirr
! 27 June 2019

  implicit none

! passed variables
  integer*4   ,intent(in)     :: n                     ! length of vector
  real(kind=4),intent(in)     :: flo                   ! low frequency cutoff
  real(kind=4),intent(in)     :: fhi                   ! high frequency cutoff
  real(kind=4),intent(in)     :: fs(n/2 + 1)           ! input frequencies
  real(kind=4),intent(out)    :: filt(n/2 + 1)         ! ouput filter
  character*10,intent(in)     :: ftype                 ! filter type (high or low)
  character(*),intent(in)     :: calling_prog          ! name of calling program
  logical     ,intent(in)     :: debug                 ! debug flag

! local variables
  integer*4                   :: i
  character*250               :: message
  real(kind=8),parameter      :: pi = 4.d0*atan(1.d0)  ! approximation of pi

  write(message,'(3a)')" Applying cosine ",ftype," pass filter to input frequencies "
  call status_update('STATUS','UTIL',calling_prog,' ',message,0)

! raised cosine high pass filter
  if( ftype .eq. "high" )then
    do i = 1,n/2 + 1
      if( fs(i) < flo )then
        filt(i) = 0.0
      else if( fs(i) > fhi )then
        filt(i) = 1.0
      else
        filt(i) = 0.5*(-cos(pi*((fs(i) - flo)/(flo - fhi))) + 1)
      endif
    enddo
! raised cosine low pass filter
  else if( ftype .eq. "low" )then
    do i = 1,n/2 + 1
      if( fs(i) < flo )then
        filt(i) = 1.0
      else if( fs(i) > fhi )then
        filt(i) = 0.0
      else
        filt(i) = 0.5*(cos(pi*((fs(i) - flo)/(flo - fhi))) + 1)
      endif
    enddo
  else
    write(message,'(a)')" Input not understood, set ftype of cosine filter to 'high' or 'low' "
    call status_update('FATAL','UTIL',calling_prog,' ',message,0)
  endif
end

subroutine rfft1b ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

!*****************************************************************************80
!
!! RFFT1B: real single precision backward fast Fourier transform, 1D.
!
!  Discussion:
!
!    RFFT1B computes the one-dimensional Fourier transform of a periodic
!    sequence within a real array.  This is referred to as the backward
!    transform or Fourier synthesis, transforming the sequence from 
!    spectral to physical space.  This transform is normalized since a 
!    call to RFFT1B followed by a call to RFFT1F (or vice-versa) reproduces
!    the original array within roundoff error. 
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.  The transform is most efficient when N is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, 
!    in array R, of two consecutive elements within the sequence. 
!
!    Input/output, real ( kind = 4 ) R(LENR), on input, the data to be 
!    transformed, and on output, the transformed data.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the R array.  
!    LENR must be at least INC*(N-1) + 1. 
!
!    Input, real ( kind = 4 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to RFFT1I before the first call to routine
!    RFFT1F or RFFT1B for a given transform length N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least N + INT(LOG(REAL(N))) + 4. 
!
!    Workspace, real ( kind = 4 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.  
!    LENWRK must be at least N. 
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENR not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough.
!
  implicit none

  integer ( kind = 4 ) lenr
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) n
  real ( kind = 4 ) r(lenr)
  real ( kind = 4 ) work(lenwrk)
  real ( kind = 4 ) wsave(lensav)

  ier = 0

  if (lenr < inc*(n-1) + 1) then
    ier = 1
    call xerfft ('rfft1b ', 6)
  else if (lensav < n + int(log( real ( n, kind = 4 ) )/log( 2.0E+00 )) +4) then
    ier = 2
    call xerfft ('rfft1b ', 8)
  else if (lenwrk < n) then
    ier = 3
    call xerfft ('rfft1b ', 10)
  end if

  if (n == 1) then
    return
  end if

  call rfftb1 (n,inc,r,work,wsave,wsave(n+1))

  return
end
subroutine rfft1f ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

!*****************************************************************************80
!
!! RFFT1F: real single precision forward fast Fourier transform, 1D.
!
!  Discussion:
!
!    RFFT1F computes the one-dimensional Fourier transform of a periodic
!    sequence within a real array.  This is referred to as the forward 
!    transform or Fourier analysis, transforming the sequence from physical
!    to spectral space.  This transform is normalized since a call to 
!    RFFT1F followed by a call to RFFT1B (or vice-versa) reproduces the 
!    original array within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.  The transform is most efficient when N is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, in 
!    array R, of two consecutive elements within the sequence. 
!
!    Input/output, real ( kind = 4 ) R(LENR), on input, contains the sequence 
!    to be transformed, and on output, the transformed data.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the R array.  
!    LENR must be at least INC*(N-1) + 1.
!
!    Input, real ( kind = 4 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to RFFT1I before the first call to routine RFFT1F
!    or RFFT1B for a given transform length N. 
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 4 ) WORK(LENWRK). 
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array. 
!    LENWRK must be at least N. 
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENR not big enough:
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough.
!
  implicit none

  integer ( kind = 4 ) lenr
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) n
  real ( kind = 4 ) work(lenwrk)
  real ( kind = 4 ) wsave(lensav)
  real ( kind = 4 ) r(lenr)

  ier = 0

  if (lenr < inc*(n-1) + 1) then
    ier = 1
    call xerfft ('rfft1f ', 6)
  else if (lensav < n + int(log( real ( n, kind = 4 ) )/log( 2.0E+00 )) +4) then
    ier = 2
    call xerfft ('rfft1f ', 8)
  else if (lenwrk < n) then
    ier = 3
    call xerfft ('rfft1f ', 10)
  end if

  if (n == 1) then
    return
  end if

  call rfftf1 (n,inc,r,work,wsave,wsave(n+1))

  return
end
subroutine rfft1i ( n, wsave, lensav, ier )

!*****************************************************************************80
!
!! RFFT1I: initialization for RFFT1B and RFFT1F.
!
!  Discussion:
!
!    RFFT1I initializes array WSAVE for use in its companion routines 
!    RFFT1B and RFFT1F.  The prime factorization of N together with a
!    tabulation of the trigonometric functions are computed and stored
!    in array WSAVE.  Separate WSAVE arrays are required for different
!    values of N. 
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.  The transform is most efficient when N is a product of 
!    small primes.
!
!    Output, real ( kind = 4 ) WSAVE(LENSAV), containing the prime factors of
!    N and also containing certain trigonometric values which will be used in
!    routines RFFT1B or RFFT1F.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least N + INT(LOG(REAL(N))) + 4.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    2, input parameter LENSAV not big enough.
!
  implicit none

  integer ( kind = 4 ) lensav

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) n
  real ( kind = 4 ) wsave(lensav)

  ier = 0

  if (lensav < n + int(log( real ( n, kind = 4 ) )/log( 2.0E+00 )) +4) then
    ier = 2
    call xerfft ('rfft1i ', 3)
  end if

  if (n == 1) then
    return
  end if

  call rffti1 (n,wsave(1),wsave(n+1))

  return
end
subroutine rfftf1 ( n, in, c, ch, wa, fac )

!*****************************************************************************80
!
!! RFFTF1 is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) in
  integer ( kind = 4 ) n

  real ( kind = 4 ) c(in,*)
  real ( kind = 4 ) ch(*)
  real ( kind = 4 ) fac(15)
  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) ix2
  integer ( kind = 4 ) ix3
  integer ( kind = 4 ) ix4
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) kh
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nl
  real ( kind = 4 ) sn
  real ( kind = 4 ) tsn
  real ( kind = 4 ) tsnm
  real ( kind = 4 ) wa(n)

      nf = int ( fac(2) )
      na = 1
      l2 = n
      iw = n

      do 111 k1=1,nf
         kh = nf-k1
         ip = fac(kh+3)
         l1 = l2/ip
         ido = n/l2
         idl1 = ido*l1
         iw = iw-(ip-1)*ido
         na = 1-na
         if (ip /= 4) go to 102
         ix2 = iw+ido
         ix3 = ix2+ido
         if (na /= 0) go to 101
	     call r1f4kf (ido,l1,c,in,ch,1,wa(iw),wa(ix2),wa(ix3))
         go to 110
  101    call r1f4kf (ido,l1,ch,1,c,in,wa(iw),wa(ix2),wa(ix3))
         go to 110
  102    if (ip /= 2) go to 104
         if (na /= 0) go to 103
	 call r1f2kf (ido,l1,c,in,ch,1,wa(iw))
         go to 110
  103    call r1f2kf (ido,l1,ch,1,c,in,wa(iw))
         go to 110
  104    if (ip /= 3) go to 106
         ix2 = iw+ido
         if (na /= 0) go to 105
	 call r1f3kf (ido,l1,c,in,ch,1,wa(iw),wa(ix2))
         go to 110
  105    call r1f3kf (ido,l1,ch,1,c,in,wa(iw),wa(ix2))
         go to 110
  106    if (ip /= 5) go to 108
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (na /= 0) go to 107
         call r1f5kf (ido,l1,c,in,ch,1,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
  107    call r1f5kf (ido,l1,ch,1,c,in,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
  108    if (ido == 1) na = 1-na
         if (na /= 0) go to 109
	 call r1fgkf (ido,ip,l1,idl1,c,c,c,in,ch,ch,1,wa(iw))
         na = 1
         go to 110
  109    call r1fgkf (ido,ip,l1,idl1,ch,ch,ch,1,c,c,in,wa(iw))
         na = 0
  110    l2 = l1
  111 continue

      sn = 1.0E+00 / real ( n, kind = 4 )
      tsn =  2.0E+00  / real ( n, kind = 4 )
      tsnm = -tsn
      modn = mod(n,2)
      nl = n-2
      if(modn /= 0) nl = n-1
      if (na /= 0) go to 120
      c(1,1) = sn*ch(1)
      do 118 j=2,nl,2
	 c(1,j) = tsn*ch(j)
	 c(1,j+1) = tsnm*ch(j+1)
  118 continue
      if(modn /= 0) return
      c(1,n) = sn*ch(n)
      return
  120 c(1,1) = sn*c(1,1)
      do 122 j=2,nl,2
	 c(1,j) = tsn*c(1,j)
	 c(1,j+1) = tsnm*c(1,j+1)
  122 continue
      if(modn /= 0) return
      c(1,n) = sn*c(1,n)

  return
end
subroutine rffti1 ( n, wa, fac )

!*****************************************************************************80
!
!! RFFTI1 is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number for which factorization 
!    and other information is needed.
!
!    Output, real ( kind = 4 ) WA(N), trigonometric information.
!
!    Output, real ( kind = 4 ) FAC(15), factorization information.  
!    FAC(1) is N, FAC(2) is NF, the number of factors, and FAC(3:NF+2) are the 
!    factors.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) arg
  real ( kind = 8 ) argh
  real ( kind = 8 ) argld
  real ( kind = 4 ) fac(15)
  real ( kind = 4 ) fi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ib
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ipm
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) ld
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nfm1
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nq
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) ntry
  integer ( kind = 4 ) ntryh(4)
  real ( kind = 8 ) tpi
  real ( kind = 4 ) wa(n)

  save ntryh

  data ntryh / 4, 2, 3, 5 /

      nl = n
      nf = 0
      j = 0
  101 j = j+1
      if (j-4) 102,102,103
  102 ntry = ntryh(j)
      go to 104
  103 ntry = ntry+2
  104 nq = nl/ntry
      nr = nl-ntry*nq
      if (nr) 101,105,101
  105 nf = nf+1
      fac(nf+2) = ntry
      nl = nq
      if (ntry /= 2) go to 107
      if (nf == 1) go to 107
      do 106 i=2,nf
         ib = nf-i+2
         fac(ib+2) = fac(ib+1)
  106 continue
      fac(3) = 2
  107 if (nl /= 1) go to 104
      fac(1) = n
      fac(2) = nf
      tpi = 8.0D+00 * atan ( 1.0D+00 )
      argh = tpi / real ( n, kind = 8 )
      is = 0
      nfm1 = nf-1
      l1 = 1
      if (nfm1 == 0) return
      do 110 k1=1,nfm1
         ip = fac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = n/l2
         ipm = ip-1
         do 109 j=1,ipm
            ld = ld+l1
            i = is
            argld = real ( ld, kind = 8 ) * argh
            fi = 0.0E+00
            do 108 ii=3,ido,2
               i = i+2
               fi = fi + 1.0E+00
               arg = fi*argld
	       wa(i-1) = cos ( arg )
	       wa(i) = sin ( arg )
  108       continue
            is = is+ido
  109    continue
         l1 = l2
  110 continue

  return
end
subroutine rfftb1 ( n, in, c, ch, wa, fac )

!*****************************************************************************80
!
!! RFFTB1 is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) in
  integer ( kind = 4 ) n

  real ( kind = 4 ) c(in,*)
  real ( kind = 4 ) ch(*)
  real ( kind = 4 ) fac(15)
  real ( kind = 4 ) half
  real ( kind = 4 ) halfm
  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) ix2
  integer ( kind = 4 ) ix3
  integer ( kind = 4 ) ix4
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nl
  real ( kind = 4 ) wa(n)

      nf = fac(2)
      na = 0
      do 10 k1=1,nf
      ip = fac(k1+2)
      na = 1-na      
      if(ip <= 5) go to 10
      if(k1 == nf) go to 10
      na = 1-na
   10 continue 
      half = 0.5E+00
      halfm = -0.5E+00
      modn = mod(n,2)
      nl = n-2
      if(modn /= 0) nl = n-1
      if (na == 0) go to 120

      ch(1) = c(1,1)
      ch(n) = c(1,n)
      do j=2,nl,2
	    ch(j) = half*c(1,j)
	    ch(j+1) = halfm*c(1,j+1)
      end do

      go to 124
  120 do 122 j=2,nl,2
	 c(1,j) = half*c(1,j)
	 c(1,j+1) = halfm*c(1,j+1)
  122 continue
  124 l1 = 1
      iw = 1
      do 116 k1=1,nf
         ip = fac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idl1 = ido*l1
         if (ip /= 4) go to 103
         ix2 = iw+ido
         ix3 = ix2+ido
         if (na /= 0) go to 101
         call r1f4kb (ido,l1,c,in,ch,1,wa(iw),wa(ix2),wa(ix3))
         go to 102
  101    call r1f4kb (ido,l1,ch,1,c,in,wa(iw),wa(ix2),wa(ix3))
  102    na = 1-na
         go to 115
  103    if (ip /= 2) go to 106
         if (na /= 0) go to 104
	 call r1f2kb (ido,l1,c,in,ch,1,wa(iw))
         go to 105
  104    call r1f2kb (ido,l1,ch,1,c,in,wa(iw))
  105    na = 1-na
         go to 115
  106    if (ip /= 3) go to 109
         ix2 = iw+ido
         if (na /= 0) go to 107
	 call r1f3kb (ido,l1,c,in,ch,1,wa(iw),wa(ix2))
         go to 108
  107    call r1f3kb (ido,l1,ch,1,c,in,wa(iw),wa(ix2))
  108    na = 1-na
         go to 115
  109    if (ip /= 5) go to 112
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (na /= 0) go to 110
         call r1f5kb (ido,l1,c,in,ch,1,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 111
  110    call r1f5kb (ido,l1,ch,1,c,in,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  111    na = 1-na
         go to 115
  112    if (na /= 0) go to 113
	 call r1fgkb (ido,ip,l1,idl1,c,c,c,in,ch,ch,1,wa(iw))
         go to 114
  113    call r1fgkb (ido,ip,l1,idl1,ch,ch,ch,1,c,c,in,wa(iw))
  114    if (ido == 1) na = 1-na
  115    l1 = l2
         iw = iw+(ip-1)*ido
  116 continue

  return
end
subroutine r1f2kb (ido,l1,cc,in1,ch,in2,wa1)

!*****************************************************************************80
!
!! R1F2KB is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 4 ) cc(in1,ido,2,l1)
  real ( kind = 4 ) ch(in2,ido,l1,2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = 4 ) wa1(ido)

  do k=1,l1
    ch(1,1,k,1) = cc(1,1,1,k)+cc(1,ido,2,k)
    ch(1,1,k,2) = cc(1,1,1,k)-cc(1,ido,2,k)
  end do

      if (ido-2) 107,105,102
 102  idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            
            ch(1,i-1,k,1) = cc(1,i-1,1,k)+cc(1,ic-1,2,k)
            ch(1,i,k,1) = cc(1,i,1,k)-cc(1,ic,2,k)
            
            ch(1,i-1,k,2) = wa1(i-2)*(cc(1,i-1,1,k)-cc(1,ic-1,2,k)) &
                 -wa1(i-1)*(cc(1,i,1,k)+cc(1,ic,2,k))
            ch(1,i,k,2) = wa1(i-2)*(cc(1,i,1,k)+cc(1,ic,2,k))+wa1(i-1) &
                 *(cc(1,i-1,1,k)-cc(1,ic-1,2,k))

 103     continue
 104  continue
      if (mod(ido,2) == 1) return
 105  do 106 k=1,l1
         ch(1,ido,k,1) = cc(1,ido,1,k)+cc(1,ido,1,k)
         ch(1,ido,k,2) = -(cc(1,1,2,k)+cc(1,1,2,k))
 106  continue
 107  continue

  return
end
subroutine r1f2kf (ido,l1,cc,in1,ch,in2,wa1)

!*****************************************************************************80
!
!! R1F1KF is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 4 ) ch(in2,ido,2,l1)
  real ( kind = 4 ) cc(in1,ido,l1,2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = 4 ) wa1(ido)

  do k=1,l1
    ch(1,1,1,k) = cc(1,1,k,1)+cc(1,1,k,2)
    ch(1,ido,2,k) = cc(1,1,k,1)-cc(1,1,k,2)
  end do

      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do i=3,ido,2
            ic = idp2-i
            ch(1,i,1,k) = cc(1,i,k,1)+(wa1(i-2)*cc(1,i,k,2) &
              -wa1(i-1)*cc(1,i-1,k,2))
            ch(1,ic,2,k) = (wa1(i-2)*cc(1,i,k,2) &
              -wa1(i-1)*cc(1,i-1,k,2))-cc(1,i,k,1)
            ch(1,i-1,1,k) = cc(1,i-1,k,1)+(wa1(i-2)*cc(1,i-1,k,2) &
              +wa1(i-1)*cc(1,i,k,2))
            ch(1,ic-1,2,k) = cc(1,i-1,k,1)-(wa1(i-2)*cc(1,i-1,k,2) &
              +wa1(i-1)*cc(1,i,k,2))
          end do
  104 continue
      if (mod(ido,2) == 1) return
  105 do 106 k=1,l1
         ch(1,1,2,k) = -cc(1,ido,k,2)
         ch(1,ido,1,k) = cc(1,ido,k,1)
  106 continue
  107 continue

  return
end
subroutine r1f3kb (ido,l1,cc,in1,ch,in2,wa1,wa2)

!*****************************************************************************80
!
!! R1F3KB is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 4 ) arg
  real ( kind = 4 ) cc(in1,ido,3,l1)
  real ( kind = 4 ) ch(in2,ido,l1,3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = 4 ) taui
  real ( kind = 4 ) taur
  real ( kind = 4 ) wa1(ido)
  real ( kind = 4 ) wa2(ido)

  arg = 2.0E+00 * 4.0E+00 * atan ( 1.0E+00 ) / 3.0E+00 
  taur = cos ( arg )
  taui = sin ( arg )

  do k = 1, l1
    ch(1,1,k,1) = cc(1,1,1,k) + 2.0E+00 * cc(1,ido,2,k)
    ch(1,1,k,2) = cc(1,1,1,k) + ( 2.0E+00 * taur ) * cc(1,ido,2,k) &
      - ( 2.0E+00 *taui)*cc(1,1,3,k)
    ch(1,1,k,3) = cc(1,1,1,k) + ( 2.0E+00 *taur)*cc(1,ido,2,k) &
      + 2.0E+00 *taui*cc(1,1,3,k)
  end do

  if (ido == 1) then
    return
  end if

  idp2 = ido+2

      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
        ch(1,i-1,k,1) = cc(1,i-1,1,k)+(cc(1,i-1,3,k)+cc(1,ic-1,2,k))
        ch(1,i,k,1) = cc(1,i,1,k)+(cc(1,i,3,k)-cc(1,ic,2,k))

        ch(1,i-1,k,2) = wa1(i-2)* &
       ((cc(1,i-1,1,k)+taur*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)))- &
       (taui*(cc(1,i,3,k)+cc(1,ic,2,k)))) &
                         -wa1(i-1)* &
       ((cc(1,i,1,k)+taur*(cc(1,i,3,k)-cc(1,ic,2,k)))+ &
       (taui*(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))) 

            ch(1,i,k,2) = wa1(i-2)* &
       ((cc(1,i,1,k)+taur*(cc(1,i,3,k)-cc(1,ic,2,k)))+ &
       (taui*(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))) &
                        +wa1(i-1)* &
       ((cc(1,i-1,1,k)+taur*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)))- &
       (taui*(cc(1,i,3,k)+cc(1,ic,2,k))))

              ch(1,i-1,k,3) = wa2(i-2)* &
       ((cc(1,i-1,1,k)+taur*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)))+ &
       (taui*(cc(1,i,3,k)+cc(1,ic,2,k)))) &
         -wa2(i-1)* &
       ((cc(1,i,1,k)+taur*(cc(1,i,3,k)-cc(1,ic,2,k)))- &
       (taui*(cc(1,i-1,3,k)-cc(1,ic-1,2,k))))

            ch(1,i,k,3) = wa2(i-2)* &
       ((cc(1,i,1,k)+taur*(cc(1,i,3,k)-cc(1,ic,2,k)))- &
       (taui*(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))) &
                       +wa2(i-1)* &
       ((cc(1,i-1,1,k)+taur*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)))+ &
       (taui*(cc(1,i,3,k)+cc(1,ic,2,k))))

  102    continue
  103 continue

  return
end
subroutine r1f3kf (ido,l1,cc,in1,ch,in2,wa1,wa2)

!*****************************************************************************80
!
!! R1F3KF is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 4 ) arg
  real ( kind = 4 ) cc(in1,ido,l1,3)
  real ( kind = 4 ) ch(in2,ido,3,l1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = 4 ) taui
  real ( kind = 4 ) taur
  real ( kind = 4 ) wa1(ido)
  real ( kind = 4 ) wa2(ido)

  arg= 2.0E+00 * 4.0E+00 * atan( 1.0E+00 )/ 3.0E+00 
  taur=cos(arg)
  taui=sin(arg)

  do k=1,l1
    ch(1,1,1,k) = cc(1,1,k,1)+(cc(1,1,k,2)+cc(1,1,k,3))
    ch(1,1,3,k) = taui*(cc(1,1,k,3)-cc(1,1,k,2))
    ch(1,ido,2,k) = cc(1,1,k,1)+taur*(cc(1,1,k,2)+cc(1,1,k,3))
  end do

  if (ido == 1) then
    return
  end if

      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i

            ch(1,i-1,1,k) = cc(1,i-1,k,1)+((wa1(i-2)*cc(1,i-1,k,2)+ &
             wa1(i-1)*cc(1,i,k,2))+(wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
             cc(1,i,k,3)))

            ch(1,i,1,k) = cc(1,i,k,1)+((wa1(i-2)*cc(1,i,k,2)- &
             wa1(i-1)*cc(1,i-1,k,2))+(wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
             cc(1,i-1,k,3)))

            ch(1,i-1,3,k) = (cc(1,i-1,k,1)+taur*((wa1(i-2)* &
             cc(1,i-1,k,2)+wa1(i-1)*cc(1,i,k,2))+(wa2(i-2)* &
             cc(1,i-1,k,3)+wa2(i-1)*cc(1,i,k,3))))+(taui*((wa1(i-2)* &
             cc(1,i,k,2)-wa1(i-1)*cc(1,i-1,k,2))-(wa2(i-2)* &
             cc(1,i,k,3)-wa2(i-1)*cc(1,i-1,k,3))))

            ch(1,ic-1,2,k) = (cc(1,i-1,k,1)+taur*((wa1(i-2)* &
             cc(1,i-1,k,2)+wa1(i-1)*cc(1,i,k,2))+(wa2(i-2)* &
             cc(1,i-1,k,3)+wa2(i-1)*cc(1,i,k,3))))-(taui*((wa1(i-2)* &
             cc(1,i,k,2)-wa1(i-1)*cc(1,i-1,k,2))-(wa2(i-2)* &
             cc(1,i,k,3)-wa2(i-1)*cc(1,i-1,k,3))))

            ch(1,i,3,k) = (cc(1,i,k,1)+taur*((wa1(i-2)*cc(1,i,k,2)- &
             wa1(i-1)*cc(1,i-1,k,2))+(wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
             cc(1,i-1,k,3))))+(taui*((wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
             cc(1,i,k,3))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
             cc(1,i,k,2))))

            ch(1,ic,2,k) = (taui*((wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
             cc(1,i,k,3))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
             cc(1,i,k,2))))-(cc(1,i,k,1)+taur*((wa1(i-2)*cc(1,i,k,2)- &
             wa1(i-1)*cc(1,i-1,k,2))+(wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
             cc(1,i-1,k,3))))
  102    continue
  103 continue

  return
end
subroutine r1f4kb (ido,l1,cc,in1,ch,in2,wa1,wa2,wa3)

!*****************************************************************************80
!
!! R1F4KB is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 4 ) cc(in1,ido,4,l1)
  real ( kind = 4 ) ch(in2,ido,l1,4)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = 4 ) sqrt2
  real ( kind = 4 ) wa1(ido)
  real ( kind = 4 ) wa2(ido)
  real ( kind = 4 ) wa3(ido)

  sqrt2=sqrt( 2.0E+00 )

  do k=1,l1
    ch(1,1,k,3) = (cc(1,1,1,k)+cc(1,ido,4,k)) &
      -(cc(1,ido,2,k)+cc(1,ido,2,k))
    ch(1,1,k,1) = (cc(1,1,1,k)+cc(1,ido,4,k)) &
      +(cc(1,ido,2,k)+cc(1,ido,2,k))
    ch(1,1,k,4) = (cc(1,1,1,k)-cc(1,ido,4,k)) &
      +(cc(1,1,3,k)+cc(1,1,3,k))
    ch(1,1,k,2) = (cc(1,1,1,k)-cc(1,ido,4,k)) &
      -(cc(1,1,3,k)+cc(1,1,3,k))
  end do

      if (ido-2) 107,105,102

  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
        ch(1,i-1,k,1) = (cc(1,i-1,1,k)+cc(1,ic-1,4,k)) &
        +(cc(1,i-1,3,k)+cc(1,ic-1,2,k))
        ch(1,i,k,1) = (cc(1,i,1,k)-cc(1,ic,4,k)) &
        +(cc(1,i,3,k)-cc(1,ic,2,k))
        ch(1,i-1,k,2)=wa1(i-2)*((cc(1,i-1,1,k)-cc(1,ic-1,4,k)) &
        -(cc(1,i,3,k)+cc(1,ic,2,k)))-wa1(i-1) &
        *((cc(1,i,1,k)+cc(1,ic,4,k))+(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))
        ch(1,i,k,2)=wa1(i-2)*((cc(1,i,1,k)+cc(1,ic,4,k)) &
        +(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))+wa1(i-1) &
        *((cc(1,i-1,1,k)-cc(1,ic-1,4,k))-(cc(1,i,3,k)+cc(1,ic,2,k)))
        ch(1,i-1,k,3)=wa2(i-2)*((cc(1,i-1,1,k)+cc(1,ic-1,4,k)) &
        -(cc(1,i-1,3,k)+cc(1,ic-1,2,k)))-wa2(i-1) &
        *((cc(1,i,1,k)-cc(1,ic,4,k))-(cc(1,i,3,k)-cc(1,ic,2,k)))
        ch(1,i,k,3)=wa2(i-2)*((cc(1,i,1,k)-cc(1,ic,4,k)) &
        -(cc(1,i,3,k)-cc(1,ic,2,k)))+wa2(i-1) &
        *((cc(1,i-1,1,k)+cc(1,ic-1,4,k))-(cc(1,i-1,3,k) &
        +cc(1,ic-1,2,k)))
        ch(1,i-1,k,4)=wa3(i-2)*((cc(1,i-1,1,k)-cc(1,ic-1,4,k)) &
        +(cc(1,i,3,k)+cc(1,ic,2,k)))-wa3(i-1) &
       *((cc(1,i,1,k)+cc(1,ic,4,k))-(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))
        ch(1,i,k,4)=wa3(i-2)*((cc(1,i,1,k)+cc(1,ic,4,k)) &
        -(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))+wa3(i-1) &
        *((cc(1,i-1,1,k)-cc(1,ic-1,4,k))+(cc(1,i,3,k)+cc(1,ic,2,k)))
  103    continue
  104 continue
      if (mod(ido,2) == 1) return
  105 continue
      do 106 k=1,l1
         ch(1,ido,k,1) = (cc(1,ido,1,k)+cc(1,ido,3,k)) &
         +(cc(1,ido,1,k)+cc(1,ido,3,k))
         ch(1,ido,k,2) = sqrt2*((cc(1,ido,1,k)-cc(1,ido,3,k)) &
         -(cc(1,1,2,k)+cc(1,1,4,k)))
         ch(1,ido,k,3) = (cc(1,1,4,k)-cc(1,1,2,k)) &
         +(cc(1,1,4,k)-cc(1,1,2,k))
         ch(1,ido,k,4) = -sqrt2*((cc(1,ido,1,k)-cc(1,ido,3,k)) &
         +(cc(1,1,2,k)+cc(1,1,4,k)))
  106 continue
  107 continue

  return
end
subroutine r1f4kf (ido,l1,cc,in1,ch,in2,wa1,wa2,wa3)

!*****************************************************************************80
!
!! R1F4KF is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 4 ) cc(in1,ido,l1,4)
  real ( kind = 4 ) ch(in2,ido,4,l1)
  real ( kind = 4 ) hsqt2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = 4 ) wa1(ido)
  real ( kind = 4 ) wa2(ido)
  real ( kind = 4 ) wa3(ido)

  hsqt2=sqrt( 2.0E+00 )/ 2.0E+00 

  do k=1,l1
    ch(1,1,1,k) = (cc(1,1,k,2)+cc(1,1,k,4))+(cc(1,1,k,1)+cc(1,1,k,3))
    ch(1,ido,4,k) = (cc(1,1,k,1)+cc(1,1,k,3))-(cc(1,1,k,2)+cc(1,1,k,4))
    ch(1,ido,2,k) = cc(1,1,k,1)-cc(1,1,k,3)
    ch(1,1,3,k) = cc(1,1,k,4)-cc(1,1,k,2)
  end do

      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            ch(1,i-1,1,k) = ((wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
             cc(1,i,k,2))+(wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
             cc(1,i,k,4)))+(cc(1,i-1,k,1)+(wa2(i-2)*cc(1,i-1,k,3)+ &
             wa2(i-1)*cc(1,i,k,3)))
            ch(1,ic-1,4,k) = (cc(1,i-1,k,1)+(wa2(i-2)*cc(1,i-1,k,3)+ &
             wa2(i-1)*cc(1,i,k,3)))-((wa1(i-2)*cc(1,i-1,k,2)+ &
             wa1(i-1)*cc(1,i,k,2))+(wa3(i-2)*cc(1,i-1,k,4)+ &
             wa3(i-1)*cc(1,i,k,4)))
            ch(1,i,1,k) = ((wa1(i-2)*cc(1,i,k,2)-wa1(i-1)* &
             cc(1,i-1,k,2))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
             cc(1,i-1,k,4)))+(cc(1,i,k,1)+(wa2(i-2)*cc(1,i,k,3)- &
             wa2(i-1)*cc(1,i-1,k,3)))
            ch(1,ic,4,k) = ((wa1(i-2)*cc(1,i,k,2)-wa1(i-1)* &
             cc(1,i-1,k,2))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
             cc(1,i-1,k,4)))-(cc(1,i,k,1)+(wa2(i-2)*cc(1,i,k,3)- &
             wa2(i-1)*cc(1,i-1,k,3)))
            ch(1,i-1,3,k) = ((wa1(i-2)*cc(1,i,k,2)-wa1(i-1)* &
             cc(1,i-1,k,2))-(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
             cc(1,i-1,k,4)))+(cc(1,i-1,k,1)-(wa2(i-2)*cc(1,i-1,k,3)+ &
             wa2(i-1)*cc(1,i,k,3)))
            ch(1,ic-1,2,k) = (cc(1,i-1,k,1)-(wa2(i-2)*cc(1,i-1,k,3)+ &
             wa2(i-1)*cc(1,i,k,3)))-((wa1(i-2)*cc(1,i,k,2)-wa1(i-1)* &
             cc(1,i-1,k,2))-(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
             cc(1,i-1,k,4)))
            ch(1,i,3,k) = ((wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
             cc(1,i,k,4))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
             cc(1,i,k,2)))+(cc(1,i,k,1)-(wa2(i-2)*cc(1,i,k,3)- &
             wa2(i-1)*cc(1,i-1,k,3)))
            ch(1,ic,2,k) = ((wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
             cc(1,i,k,4))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
             cc(1,i,k,2)))-(cc(1,i,k,1)-(wa2(i-2)*cc(1,i,k,3)- &
             wa2(i-1)*cc(1,i-1,k,3)))
  103    continue
  104 continue
      if (mod(ido,2) == 1) return
  105 continue
      do 106 k=1,l1
            ch(1,ido,1,k) = (hsqt2*(cc(1,ido,k,2)-cc(1,ido,k,4)))+cc(1,ido,k,1)
            ch(1,ido,3,k) = cc(1,ido,k,1)-(hsqt2*(cc(1,ido,k,2)-cc(1,ido,k,4)))
            ch(1,1,2,k) = (-hsqt2*(cc(1,ido,k,2)+cc(1,ido,k,4)))-cc(1,ido,k,3)
            ch(1,1,4,k) = (-hsqt2*(cc(1,ido,k,2)+cc(1,ido,k,4)))+cc(1,ido,k,3)
  106 continue
  107 continue

  return
end
subroutine r1f5kb (ido,l1,cc,in1,ch,in2,wa1,wa2,wa3,wa4)

!*****************************************************************************80
!
!! R1F5KB is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 4 ) arg
  real ( kind = 4 ) cc(in1,ido,5,l1)
  real ( kind = 4 ) ch(in2,ido,l1,5)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = 4 ) ti11
  real ( kind = 4 ) ti12
  real ( kind = 4 ) tr11
  real ( kind = 4 ) tr12
  real ( kind = 4 ) wa1(ido)
  real ( kind = 4 ) wa2(ido)
  real ( kind = 4 ) wa3(ido)
  real ( kind = 4 ) wa4(ido)

  arg= 2.0E+00 * 4.0E+00 * atan( 1.0E+00 ) / 5.0E+00
  tr11=cos(arg)
  ti11=sin(arg)
  tr12=cos( 2.0E+00 *arg )
  ti12=sin( 2.0E+00 *arg )

  do k=1,l1
    ch(1,1,k,1) = cc(1,1,1,k)+ 2.0E+00 *cc(1,ido,2,k)+ 2.0E+00 *cc(1,ido,4,k)
    ch(1,1,k,2) = (cc(1,1,1,k)+tr11* 2.0E+00 *cc(1,ido,2,k) &
      +tr12* 2.0E+00 *cc(1,ido,4,k))-(ti11* 2.0E+00 *cc(1,1,3,k) &
      +ti12* 2.0E+00 *cc(1,1,5,k))
    ch(1,1,k,3) = (cc(1,1,1,k)+tr12* 2.0E+00 *cc(1,ido,2,k) &
      +tr11* 2.0E+00 *cc(1,ido,4,k))-(ti12* 2.0E+00 *cc(1,1,3,k) &
      -ti11* 2.0E+00 *cc(1,1,5,k))
    ch(1,1,k,4) = (cc(1,1,1,k)+tr12* 2.0E+00 *cc(1,ido,2,k) &
      +tr11* 2.0E+00 *cc(1,ido,4,k))+(ti12* 2.0E+00 *cc(1,1,3,k) &
      -ti11* 2.0E+00 *cc(1,1,5,k))
    ch(1,1,k,5) = (cc(1,1,1,k)+tr11* 2.0E+00 *cc(1,ido,2,k) &
      +tr12* 2.0E+00 *cc(1,ido,4,k))+(ti11* 2.0E+00 *cc(1,1,3,k) &
      +ti12* 2.0E+00 *cc(1,1,5,k))
  end do

  if (ido == 1) return

      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
        ch(1,i-1,k,1) = cc(1,i-1,1,k)+(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +(cc(1,i-1,5,k)+cc(1,ic-1,4,k))
        ch(1,i,k,1) = cc(1,i,1,k)+(cc(1,i,3,k)-cc(1,ic,2,k)) &
        +(cc(1,i,5,k)-cc(1,ic,4,k))
        ch(1,i-1,k,2) = wa1(i-2)*((cc(1,i-1,1,k)+tr11* &
        (cc(1,i-1,3,k)+cc(1,ic-1,2,k))+tr12 &
        *(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))-(ti11*(cc(1,i,3,k) &
        +cc(1,ic,2,k))+ti12*(cc(1,i,5,k)+cc(1,ic,4,k)))) &
        -wa1(i-1)*((cc(1,i,1,k)+tr11*(cc(1,i,3,k)-cc(1,ic,2,k)) &
        +tr12*(cc(1,i,5,k)-cc(1,ic,4,k)))+(ti11*(cc(1,i-1,3,k) &
        -cc(1,ic-1,2,k))+ti12*(cc(1,i-1,5,k)-cc(1,ic-1,4,k))))

        ch(1,i,k,2) = wa1(i-2)*((cc(1,i,1,k)+tr11*(cc(1,i,3,k) &
        -cc(1,ic,2,k))+tr12*(cc(1,i,5,k)-cc(1,ic,4,k))) &
        +(ti11*(cc(1,i-1,3,k)-cc(1,ic-1,2,k))+ti12 &
        *(cc(1,i-1,5,k)-cc(1,ic-1,4,k))))+wa1(i-1) &
        *((cc(1,i-1,1,k)+tr11*(cc(1,i-1,3,k) &
        +cc(1,ic-1,2,k))+tr12*(cc(1,i-1,5,k)+cc(1,ic-1,4,k))) &
        -(ti11*(cc(1,i,3,k)+cc(1,ic,2,k))+ti12 &
        *(cc(1,i,5,k)+cc(1,ic,4,k))))

        ch(1,i-1,k,3) = wa2(i-2) &
        *((cc(1,i-1,1,k)+tr12*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +tr11*(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))-(ti12*(cc(1,i,3,k) &
        +cc(1,ic,2,k))-ti11*(cc(1,i,5,k)+cc(1,ic,4,k)))) &
       -wa2(i-1) &
       *((cc(1,i,1,k)+tr12*(cc(1,i,3,k)- &
        cc(1,ic,2,k))+tr11*(cc(1,i,5,k)-cc(1,ic,4,k))) &
        +(ti12*(cc(1,i-1,3,k)-cc(1,ic-1,2,k))-ti11 &
        *(cc(1,i-1,5,k)-cc(1,ic-1,4,k))))

        ch(1,i,k,3) = wa2(i-2) &
       *((cc(1,i,1,k)+tr12*(cc(1,i,3,k)- &
        cc(1,ic,2,k))+tr11*(cc(1,i,5,k)-cc(1,ic,4,k))) &
        +(ti12*(cc(1,i-1,3,k)-cc(1,ic-1,2,k))-ti11 &
        *(cc(1,i-1,5,k)-cc(1,ic-1,4,k)))) &
        +wa2(i-1) &
        *((cc(1,i-1,1,k)+tr12*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +tr11*(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))-(ti12*(cc(1,i,3,k) &
        +cc(1,ic,2,k))-ti11*(cc(1,i,5,k)+cc(1,ic,4,k))))

        ch(1,i-1,k,4) = wa3(i-2) &
        *((cc(1,i-1,1,k)+tr12*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +tr11*(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))+(ti12*(cc(1,i,3,k) &
        +cc(1,ic,2,k))-ti11*(cc(1,i,5,k)+cc(1,ic,4,k)))) &
        -wa3(i-1) &
       *((cc(1,i,1,k)+tr12*(cc(1,i,3,k)- &
        cc(1,ic,2,k))+tr11*(cc(1,i,5,k)-cc(1,ic,4,k))) &
        -(ti12*(cc(1,i-1,3,k)-cc(1,ic-1,2,k))-ti11 &
        *(cc(1,i-1,5,k)-cc(1,ic-1,4,k))))

        ch(1,i,k,4) = wa3(i-2) &
       *((cc(1,i,1,k)+tr12*(cc(1,i,3,k)- &
        cc(1,ic,2,k))+tr11*(cc(1,i,5,k)-cc(1,ic,4,k))) &
        -(ti12*(cc(1,i-1,3,k)-cc(1,ic-1,2,k))-ti11 &
        *(cc(1,i-1,5,k)-cc(1,ic-1,4,k)))) &
        +wa3(i-1) &
        *((cc(1,i-1,1,k)+tr12*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +tr11*(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))+(ti12*(cc(1,i,3,k) &
        +cc(1,ic,2,k))-ti11*(cc(1,i,5,k)+cc(1,ic,4,k))))

        ch(1,i-1,k,5) = wa4(i-2) &
        *((cc(1,i-1,1,k)+tr11*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +tr12*(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))+(ti11*(cc(1,i,3,k) &
        +cc(1,ic,2,k))+ti12*(cc(1,i,5,k)+cc(1,ic,4,k)))) &
        -wa4(i-1) &
        *((cc(1,i,1,k)+tr11*(cc(1,i,3,k)-cc(1,ic,2,k)) &
        +tr12*(cc(1,i,5,k)-cc(1,ic,4,k)))-(ti11*(cc(1,i-1,3,k) &
        -cc(1,ic-1,2,k))+ti12*(cc(1,i-1,5,k)-cc(1,ic-1,4,k))))

        ch(1,i,k,5) = wa4(i-2) &
        *((cc(1,i,1,k)+tr11*(cc(1,i,3,k)-cc(1,ic,2,k)) &
        +tr12*(cc(1,i,5,k)-cc(1,ic,4,k)))-(ti11*(cc(1,i-1,3,k) &
        -cc(1,ic-1,2,k))+ti12*(cc(1,i-1,5,k)-cc(1,ic-1,4,k)))) &
        +wa4(i-1) &
        *((cc(1,i-1,1,k)+tr11*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +tr12*(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))+(ti11*(cc(1,i,3,k) &
        +cc(1,ic,2,k))+ti12*(cc(1,i,5,k)+cc(1,ic,4,k))))

  102    continue
  103 continue

  return
end
subroutine r1f5kf (ido,l1,cc,in1,ch,in2,wa1,wa2,wa3,wa4)

!*****************************************************************************80
!
!! R1F5KF is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 4 ) arg
  real ( kind = 4 ) cc(in1,ido,l1,5)
  real ( kind = 4 ) ch(in2,ido,5,l1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = 4 ) ti11
  real ( kind = 4 ) ti12
  real ( kind = 4 ) tr11
  real ( kind = 4 ) tr12
  real ( kind = 4 ) wa1(ido)
  real ( kind = 4 ) wa2(ido)
  real ( kind = 4 ) wa3(ido)
  real ( kind = 4 ) wa4(ido)

  arg= 2.0E+00 * 4.0E+00 * atan( 1.0E+00 ) / 5.0E+00
  tr11=cos(arg)
  ti11=sin(arg)
  tr12=cos( 2.0E+00 *arg)
  ti12=sin( 2.0E+00 *arg)

  do k=1,l1
    ch(1,1,1,k) = cc(1,1,k,1)+(cc(1,1,k,5)+cc(1,1,k,2))+ &
      (cc(1,1,k,4)+cc(1,1,k,3))
    ch(1,ido,2,k) = cc(1,1,k,1)+tr11*(cc(1,1,k,5)+cc(1,1,k,2))+ &
      tr12*(cc(1,1,k,4)+cc(1,1,k,3))
    ch(1,1,3,k) = ti11*(cc(1,1,k,5)-cc(1,1,k,2))+ti12* &
      (cc(1,1,k,4)-cc(1,1,k,3))
    ch(1,ido,4,k) = cc(1,1,k,1)+tr12*(cc(1,1,k,5)+cc(1,1,k,2))+ &
      tr11*(cc(1,1,k,4)+cc(1,1,k,3))
    ch(1,1,5,k) = ti12*(cc(1,1,k,5)-cc(1,1,k,2))-ti11* &
      (cc(1,1,k,4)-cc(1,1,k,3))
  end do

  if (ido == 1) return

      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i

            ch(1,i-1,1,k) = cc(1,i-1,k,1)+((wa1(i-2)*cc(1,i-1,k,2)+ &
            wa1(i-1)*cc(1,i,k,2))+(wa4(i-2)*cc(1,i-1,k,5)+wa4(i-1)* &
            cc(1,i,k,5)))+((wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
            cc(1,i,k,3))+(wa3(i-2)*cc(1,i-1,k,4)+ &
            wa3(i-1)*cc(1,i,k,4))) 

            ch(1,i,1,k) = cc(1,i,k,1)+((wa1(i-2)*cc(1,i,k,2)- &
             wa1(i-1)*cc(1,i-1,k,2))+(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)* &
             cc(1,i-1,k,5)))+((wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
             cc(1,i-1,k,3))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
             cc(1,i-1,k,4)))

            ch(1,i-1,3,k) = cc(1,i-1,k,1)+tr11* &
            ( wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)*cc(1,i,k,2) &
             +wa4(i-2)*cc(1,i-1,k,5)+wa4(i-1)*cc(1,i,k,5))+tr12* &
            ( wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)*cc(1,i,k,3) &
             +wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)*cc(1,i,k,4))+ti11* &
            ( wa1(i-2)*cc(1,i,k,2)-wa1(i-1)*cc(1,i-1,k,2) &
             -(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)*cc(1,i-1,k,5)))+ti12* &
            ( wa2(i-2)*cc(1,i,k,3)-wa2(i-1)*cc(1,i-1,k,3) &
             -(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)*cc(1,i-1,k,4)))

            ch(1,ic-1,2,k) = cc(1,i-1,k,1)+tr11* &
            ( wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)*cc(1,i,k,2) &
             +wa4(i-2)*cc(1,i-1,k,5)+wa4(i-1)*cc(1,i,k,5))+tr12* &
           ( wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)*cc(1,i,k,3) &
            +wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)*cc(1,i,k,4))-(ti11* &
            ( wa1(i-2)*cc(1,i,k,2)-wa1(i-1)*cc(1,i-1,k,2) &
             -(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)*cc(1,i-1,k,5)))+ti12* &
            ( wa2(i-2)*cc(1,i,k,3)-wa2(i-1)*cc(1,i-1,k,3) &
             -(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)*cc(1,i-1,k,4))))

            ch(1,i,3,k) = (cc(1,i,k,1)+tr11*((wa1(i-2)*cc(1,i,k,2)- &
             wa1(i-1)*cc(1,i-1,k,2))+(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)* &
             cc(1,i-1,k,5)))+tr12*((wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
             cc(1,i-1,k,3))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
             cc(1,i-1,k,4))))+(ti11*((wa4(i-2)*cc(1,i-1,k,5)+ &
             wa4(i-1)*cc(1,i,k,5))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
             cc(1,i,k,2)))+ti12*((wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
             cc(1,i,k,4))-(wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
             cc(1,i,k,3))))

            ch(1,ic,2,k) = (ti11*((wa4(i-2)*cc(1,i-1,k,5)+wa4(i-1)* &
             cc(1,i,k,5))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
             cc(1,i,k,2)))+ti12*((wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
             cc(1,i,k,4))-(wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
             cc(1,i,k,3))))-(cc(1,i,k,1)+tr11*((wa1(i-2)*cc(1,i,k,2)- &
             wa1(i-1)*cc(1,i-1,k,2))+(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)* &
             cc(1,i-1,k,5)))+tr12*((wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
             cc(1,i-1,k,3))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
             cc(1,i-1,k,4))))

            ch(1,i-1,5,k) = (cc(1,i-1,k,1)+tr12*((wa1(i-2)* &
             cc(1,i-1,k,2)+wa1(i-1)*cc(1,i,k,2))+(wa4(i-2)* &
             cc(1,i-1,k,5)+wa4(i-1)*cc(1,i,k,5)))+tr11*((wa2(i-2)* &
             cc(1,i-1,k,3)+wa2(i-1)*cc(1,i,k,3))+(wa3(i-2)* &
             cc(1,i-1,k,4)+wa3(i-1)*cc(1,i,k,4))))+(ti12*((wa1(i-2)* &
             cc(1,i,k,2)-wa1(i-1)*cc(1,i-1,k,2))-(wa4(i-2)* &
             cc(1,i,k,5)-wa4(i-1)*cc(1,i-1,k,5)))-ti11*((wa2(i-2)* &
             cc(1,i,k,3)-wa2(i-1)*cc(1,i-1,k,3))-(wa3(i-2)* &
             cc(1,i,k,4)-wa3(i-1)*cc(1,i-1,k,4))))

            ch(1,ic-1,4,k) = (cc(1,i-1,k,1)+tr12*((wa1(i-2)* &
             cc(1,i-1,k,2)+wa1(i-1)*cc(1,i,k,2))+(wa4(i-2)* &
             cc(1,i-1,k,5)+wa4(i-1)*cc(1,i,k,5)))+tr11*((wa2(i-2)* &
             cc(1,i-1,k,3)+wa2(i-1)*cc(1,i,k,3))+(wa3(i-2)* &
             cc(1,i-1,k,4)+wa3(i-1)*cc(1,i,k,4))))-(ti12*((wa1(i-2)* &
             cc(1,i,k,2)-wa1(i-1)*cc(1,i-1,k,2))-(wa4(i-2)* &
             cc(1,i,k,5)-wa4(i-1)*cc(1,i-1,k,5)))-ti11*((wa2(i-2)* &
             cc(1,i,k,3)-wa2(i-1)*cc(1,i-1,k,3))-(wa3(i-2)* &
             cc(1,i,k,4)-wa3(i-1)*cc(1,i-1,k,4))))

            ch(1,i,5,k) = (cc(1,i,k,1)+tr12*((wa1(i-2)*cc(1,i,k,2)- &
             wa1(i-1)*cc(1,i-1,k,2))+(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)* &
             cc(1,i-1,k,5)))+tr11*((wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
             cc(1,i-1,k,3))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
             cc(1,i-1,k,4))))+(ti12*((wa4(i-2)*cc(1,i-1,k,5)+ &
             wa4(i-1)*cc(1,i,k,5))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
             cc(1,i,k,2)))-ti11*((wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
             cc(1,i,k,4))-(wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
             cc(1,i,k,3))))

            ch(1,ic,4,k) = (ti12*((wa4(i-2)*cc(1,i-1,k,5)+wa4(i-1)* &
             cc(1,i,k,5))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
             cc(1,i,k,2)))-ti11*((wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
             cc(1,i,k,4))-(wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
             cc(1,i,k,3))))-(cc(1,i,k,1)+tr12*((wa1(i-2)*cc(1,i,k,2)- &
             wa1(i-1)*cc(1,i-1,k,2))+(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)* &
             cc(1,i-1,k,5)))+tr11*((wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
             cc(1,i-1,k,3))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
             cc(1,i-1,k,4))))

  102    continue
  103 continue

  return
end
subroutine r1fgkb (ido,ip,l1,idl1,cc,c1,c2,in1,ch,ch2,in2,wa)

!*****************************************************************************80
!
!! R1FGKB is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) l1

  real ( kind = 4 ) ai1
  real ( kind = 4 ) ai2
  real ( kind = 4 ) ar1
  real ( kind = 4 ) ar1h
  real ( kind = 4 ) ar2
  real ( kind = 4 ) ar2h
  real ( kind = 4 ) arg
  real ( kind = 4 ) c1(in1,ido,l1,ip)
  real ( kind = 4 ) c2(in1,idl1,ip)
  real ( kind = 4 ) cc(in1,ido,ip,l1)
  real ( kind = 4 ) ch(in2,ido,l1,ip)
  real ( kind = 4 ) ch2(in2,idl1,ip)
  real ( kind = 4 ) dc2
  real ( kind = 4 ) dcp
  real ( kind = 4 ) ds2
  real ( kind = 4 ) dsp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idij
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) ipp2
  integer ( kind = 4 ) ipph
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) nbd
  real ( kind = 4 ) tpi
  real ( kind = 4 ) wa(ido)

  tpi= 2.0E+00 * 4.0E+00 * atan( 1.0E+00 )
  arg = tpi / real ( ip, kind = 4 )
  dcp = cos(arg)
  dsp = sin(arg)
  idp2 = ido+2
  nbd = (ido-1)/2
  ipp2 = ip+2
  ipph = (ip+1)/2

  if (ido < l1) go to 103

      do k=1,l1
         do i=1,ido
            ch(1,i,k,1) = cc(1,i,1,k)
         end do
      end do

      go to 106

  103 continue

  do i=1,ido
    do k=1,l1
      ch(1,i,k,1) = cc(1,i,1,k)
    end do
  end do

  106 do 108 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 107 k=1,l1
            ch(1,1,k,j) = cc(1,ido,j2-2,k)+cc(1,ido,j2-2,k)
            ch(1,1,k,jc) = cc(1,1,j2-1,k)+cc(1,1,j2-1,k)
 1007       continue
  107    continue
  108 continue
      if (ido == 1) go to 116
      if (nbd < l1) go to 112
      do 111 j=2,ipph
         jc = ipp2-j
         do 110 k=1,l1
            do 109 i=3,ido,2
               ic = idp2-i
               ch(1,i-1,k,j) = cc(1,i-1,2*j-1,k)+cc(1,ic-1,2*j-2,k)
               ch(1,i-1,k,jc) = cc(1,i-1,2*j-1,k)-cc(1,ic-1,2*j-2,k)
               ch(1,i,k,j) = cc(1,i,2*j-1,k)-cc(1,ic,2*j-2,k)
               ch(1,i,k,jc) = cc(1,i,2*j-1,k)+cc(1,ic,2*j-2,k)
  109       continue
  110    continue
  111 continue
      go to 116
  112 do 115 j=2,ipph
         jc = ipp2-j
         do 114 i=3,ido,2
            ic = idp2-i
            do 113 k=1,l1
               ch(1,i-1,k,j) = cc(1,i-1,2*j-1,k)+cc(1,ic-1,2*j-2,k)
               ch(1,i-1,k,jc) = cc(1,i-1,2*j-1,k)-cc(1,ic-1,2*j-2,k)
               ch(1,i,k,j) = cc(1,i,2*j-1,k)-cc(1,ic,2*j-2,k)
               ch(1,i,k,jc) = cc(1,i,2*j-1,k)+cc(1,ic,2*j-2,k)
  113       continue
  114    continue
  115 continue
  116 ar1 = 1.0E+00
      ai1 = 0.0E+00
      do 120 l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
         do 117 ik=1,idl1
            c2(1,ik,l) = ch2(1,ik,1)+ar1*ch2(1,ik,2)
            c2(1,ik,lc) = ai1*ch2(1,ik,ip)
  117    continue
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do 119 j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
            do 118 ik=1,idl1
               c2(1,ik,l) = c2(1,ik,l)+ar2*ch2(1,ik,j)
               c2(1,ik,lc) = c2(1,ik,lc)+ai2*ch2(1,ik,jc)
  118       continue
  119    continue
  120 continue
      do 122 j=2,ipph
         do 121 ik=1,idl1
            ch2(1,ik,1) = ch2(1,ik,1)+ch2(1,ik,j)
  121    continue
  122 continue
      do 124 j=2,ipph
         jc = ipp2-j
         do 123 k=1,l1
            ch(1,1,k,j) = c1(1,1,k,j)-c1(1,1,k,jc)
            ch(1,1,k,jc) = c1(1,1,k,j)+c1(1,1,k,jc)
  123    continue
  124 continue
      if (ido == 1) go to 132
      if (nbd < l1) go to 128
      do 127 j=2,ipph
         jc = ipp2-j
         do 126 k=1,l1
            do 125 i=3,ido,2
               ch(1,i-1,k,j) = c1(1,i-1,k,j)-c1(1,i,k,jc)
               ch(1,i-1,k,jc) = c1(1,i-1,k,j)+c1(1,i,k,jc)
               ch(1,i,k,j) = c1(1,i,k,j)+c1(1,i-1,k,jc)
               ch(1,i,k,jc) = c1(1,i,k,j)-c1(1,i-1,k,jc)
  125       continue
  126    continue
  127 continue
      go to 132
  128 do 131 j=2,ipph
         jc = ipp2-j
         do 130 i=3,ido,2
            do 129 k=1,l1
               ch(1,i-1,k,j) = c1(1,i-1,k,j)-c1(1,i,k,jc)
               ch(1,i-1,k,jc) = c1(1,i-1,k,j)+c1(1,i,k,jc)
               ch(1,i,k,j) = c1(1,i,k,j)+c1(1,i-1,k,jc)
               ch(1,i,k,jc) = c1(1,i,k,j)-c1(1,i-1,k,jc)
  129       continue
  130    continue
  131 continue
  132 continue
      if (ido == 1) return
      do 133 ik=1,idl1
         c2(1,ik,1) = ch2(1,ik,1)
  133 continue
      do 135 j=2,ip
         do 134 k=1,l1
            c1(1,1,k,j) = ch(1,1,k,j)
  134    continue
  135 continue
      if ( l1 < nbd ) go to 139
      is = -ido
      do 138 j=2,ip
         is = is+ido
         idij = is
         do 137 i=3,ido,2
            idij = idij+2
            do 136 k=1,l1
               c1(1,i-1,k,j) = wa(idij-1)*ch(1,i-1,k,j)-wa(idij)*ch(1,i,k,j)
               c1(1,i,k,j) = wa(idij-1)*ch(1,i,k,j)+wa(idij)*ch(1,i-1,k,j)
  136       continue
  137    continue
  138 continue
      go to 143
  139 is = -ido
      do 142 j=2,ip
         is = is+ido
         do 141 k=1,l1
            idij = is
            do 140 i=3,ido,2
               idij = idij+2
               c1(1,i-1,k,j) = wa(idij-1)*ch(1,i-1,k,j)-wa(idij)*ch(1,i,k,j)
               c1(1,i,k,j) = wa(idij-1)*ch(1,i,k,j)+wa(idij)*ch(1,i-1,k,j)
  140       continue
  141    continue
  142 continue
  143 continue

  return
end
subroutine r1fgkf (ido,ip,l1,idl1,cc,c1,c2,in1,ch,ch2,in2,wa)

!*****************************************************************************80
!
!! R1FGKF is an FFTPACK5.1 auxilliary function.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) l1

  real ( kind = 4 ) ai1
  real ( kind = 4 ) ai2
  real ( kind = 4 ) ar1
  real ( kind = 4 ) ar1h
  real ( kind = 4 ) ar2
  real ( kind = 4 ) ar2h
  real ( kind = 4 ) arg
  real ( kind = 4 ) c1(in1,ido,l1,ip)
  real ( kind = 4 ) c2(in1,idl1,ip)
  real ( kind = 4 ) cc(in1,ido,ip,l1)
  real ( kind = 4 ) ch(in2,ido,l1,ip)
  real ( kind = 4 ) ch2(in2,idl1,ip)
  real ( kind = 4 ) dc2
  real ( kind = 4 ) dcp
  real ( kind = 4 ) ds2
  real ( kind = 4 ) dsp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idij
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) ipp2
  integer ( kind = 4 ) ipph
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) nbd
  real ( kind = 4 ) tpi
  real ( kind = 4 ) wa(ido)

  tpi= 2.0E+00 * 4.0E+00 * atan( 1.0E+00 )
  arg = tpi/real ( ip, kind = 4 )
  dcp = cos(arg)
  dsp = sin(arg)
  ipph = (ip+1)/2
  ipp2 = ip+2
  idp2 = ido+2
  nbd = (ido-1)/2

      if (ido == 1) go to 119

      do ik=1,idl1
         ch2(1,ik,1) = c2(1,ik,1)
      end do

      do j=2,ip
         do k=1,l1
            ch(1,1,k,j) = c1(1,1,k,j)
         end do
      end do

      if ( l1 < nbd ) go to 107
      is = -ido
      do 106 j=2,ip
         is = is+ido
         idij = is
         do 105 i=3,ido,2
            idij = idij+2
            do 104 k=1,l1
               ch(1,i-1,k,j) = wa(idij-1)*c1(1,i-1,k,j)+wa(idij)*c1(1,i,k,j)
               ch(1,i,k,j) = wa(idij-1)*c1(1,i,k,j)-wa(idij)*c1(1,i-1,k,j)
  104       continue
  105    continue
  106 continue
      go to 111
  107 is = -ido
      do 110 j=2,ip
         is = is+ido
         do 109 k=1,l1
            idij = is
            do 108 i=3,ido,2
               idij = idij+2
               ch(1,i-1,k,j) = wa(idij-1)*c1(1,i-1,k,j)+wa(idij)*c1(1,i,k,j)
               ch(1,i,k,j) = wa(idij-1)*c1(1,i,k,j)-wa(idij)*c1(1,i-1,k,j)
  108       continue
  109    continue
  110 continue
  111 if (nbd < l1) go to 115
      do 114 j=2,ipph
         jc = ipp2-j
         do 113 k=1,l1
            do 112 i=3,ido,2
               c1(1,i-1,k,j) = ch(1,i-1,k,j)+ch(1,i-1,k,jc)
               c1(1,i-1,k,jc) = ch(1,i,k,j)-ch(1,i,k,jc)
               c1(1,i,k,j) = ch(1,i,k,j)+ch(1,i,k,jc)
               c1(1,i,k,jc) = ch(1,i-1,k,jc)-ch(1,i-1,k,j)
  112       continue
  113    continue
  114 continue
      go to 121
  115 do 118 j=2,ipph
         jc = ipp2-j
         do 117 i=3,ido,2
            do 116 k=1,l1
               c1(1,i-1,k,j) = ch(1,i-1,k,j)+ch(1,i-1,k,jc)
               c1(1,i-1,k,jc) = ch(1,i,k,j)-ch(1,i,k,jc)
               c1(1,i,k,j) = ch(1,i,k,j)+ch(1,i,k,jc)
               c1(1,i,k,jc) = ch(1,i-1,k,jc)-ch(1,i-1,k,j)
  116       continue
  117    continue
  118 continue
      go to 121
  119 do 120 ik=1,idl1
         c2(1,ik,1) = ch2(1,ik,1)
  120 continue
  121 do 123 j=2,ipph
         jc = ipp2-j
         do 122 k=1,l1
            c1(1,1,k,j) = ch(1,1,k,j)+ch(1,1,k,jc)
            c1(1,1,k,jc) = ch(1,1,k,jc)-ch(1,1,k,j)
  122    continue
  123 continue

      ar1 = 1.0E+00
      ai1 = 0.0E+00
      do 127 l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
         do 124 ik=1,idl1
            ch2(1,ik,l) = c2(1,ik,1)+ar1*c2(1,ik,2)
            ch2(1,ik,lc) = ai1*c2(1,ik,ip)
  124    continue
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do 126 j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
            do 125 ik=1,idl1
               ch2(1,ik,l) = ch2(1,ik,l)+ar2*c2(1,ik,j)
               ch2(1,ik,lc) = ch2(1,ik,lc)+ai2*c2(1,ik,jc)
  125       continue
  126    continue
  127 continue
      do 129 j=2,ipph
         do 128 ik=1,idl1
            ch2(1,ik,1) = ch2(1,ik,1)+c2(1,ik,j)
  128    continue
  129 continue

      if (ido < l1) go to 132
      do 131 k=1,l1
         do 130 i=1,ido
            cc(1,i,1,k) = ch(1,i,k,1)
  130    continue
  131 continue
      go to 135
  132 do 134 i=1,ido
         do 133 k=1,l1
            cc(1,i,1,k) = ch(1,i,k,1)
  133    continue
  134 continue
  135 do 137 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 136 k=1,l1
            cc(1,ido,j2-2,k) = ch(1,1,k,j)
            cc(1,1,j2-1,k) = ch(1,1,k,jc)
  136    continue
  137 continue
      if (ido == 1) return
      if (nbd < l1) go to 141
      do 140 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 139 k=1,l1
            do 138 i=3,ido,2
               ic = idp2-i
               cc(1,i-1,j2-1,k) = ch(1,i-1,k,j)+ch(1,i-1,k,jc)
               cc(1,ic-1,j2-2,k) = ch(1,i-1,k,j)-ch(1,i-1,k,jc)
               cc(1,i,j2-1,k) = ch(1,i,k,j)+ch(1,i,k,jc)
               cc(1,ic,j2-2,k) = ch(1,i,k,jc)-ch(1,i,k,j)
  138       continue
  139    continue
  140 continue
      return
  141 do 144 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 143 i=3,ido,2
            ic = idp2-i
            do 142 k=1,l1
               cc(1,i-1,j2-1,k) = ch(1,i-1,k,j)+ch(1,i-1,k,jc)
               cc(1,ic-1,j2-2,k) = ch(1,i-1,k,j)-ch(1,i-1,k,jc)
               cc(1,i,j2-1,k) = ch(1,i,k,j)+ch(1,i,k,jc)
               cc(1,ic,j2-2,k) = ch(1,i,k,jc)-ch(1,i,k,j)
  142       continue
  143    continue
  144 continue

  return
end
subroutine xerfft ( srname, info )

!*****************************************************************************80
!
!! XERFFT is an error handler for the FFTPACK routines.
!
!  Discussion:
!
!    XERFFT is an error handler for FFTPACK version 5.1 routines.
!    It is called by an FFTPACK 5.1 routine if an input parameter has an
!    invalid value.  A message is printed and execution stops.
!
!    Installers may consider modifying the stop statement in order to
!    call system-specific exception-handling facilities.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    31 July 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, character ( len = * ) SRNAME, the name of the calling routine.
!
!    Input, integer ( kind = 4 ) INFO, an error code.  When a single invalid 
!    parameter in the parameter list of the calling routine has been detected, 
!    INFO is the position of that parameter.  In the case when an illegal 
!    combination of LOT, JUMP, N, and INC has been detected, the calling 
!    subprogram calls XERFFT with INFO = -1.
!
  implicit none

  integer ( kind = 4 ) info
  character ( len = * ) srname

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'XERFFT - Fatal error!'

  if ( 1 <= info ) then
    write ( *, '(a,a,a,i3,a)') '  On entry to ', trim ( srname ), &
      ' parameter number ', info, ' had an illegal value.'
  else if ( info == -1 ) then
    write( *, '(a,a,a,a)') '  On entry to ', trim ( srname ), &
      ' parameters LOT, JUMP, N and INC are inconsistent.'
  else if ( info == -2 ) then
    write( *, '(a,a,a,a)') '  On entry to ', trim ( srname ), &
      ' parameter L is greater than LDIM.'
  else if ( info == -3 ) then
    write( *, '(a,a,a,a)') '  On entry to ', trim ( srname ), &
      ' parameter M is greater than MDIM.'
  else if ( info == -5 ) then
    write( *, '(a,a,a,a)') '  Within ', trim ( srname ), &
      ' input error returned by lower level routine.'
  else if ( info == -6 ) then
    write( *, '(a,a,a,a)') '  On entry to ', trim ( srname ), &
      ' parameter LDIM is less than 2*(L/2+1).'
  end if

  stop
end

