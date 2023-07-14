!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine emd_v3(indata, nvals, ncomp, rslt)
	  
! subroutine to decompose a time series using the empirical mode decomposition procedure
! where indata is the time series, nvals is rge length, ncomp is the number of imfs and
! rslt is the matrix containing successive IMFs and the residue
! RM181203 an update of the emd_v2 subroutine, added zero crossings stopping criteria
	  
  implicit none

! passed variables
  integer*4,   intent(in)  :: nvals, ncomp 
  real(kind=8),intent(in)  :: indata(nvals)      
  real(kind=8),intent(out) :: rslt(ncomp,nvals)	

! local variables
  integer*4                :: nmax, nmin, nimf, ki, k1, k2, ik, nzc, nsc
  integer*4                :: i, j, j1, im, j2, ii  
  real(kind=8)             :: Dv0, Dv1, Dv2
  real(kind=8)             :: maxess, vects, miness
  real(kind=8)             :: ispline
  real(kind=8), dimension(nvals) :: hmaxes, hmines, minenv, maxenv
  real(kind=8), dimension(nvals) :: maxees, mines, vect
  real(kind=8), dimension(nvals) :: diff, maxmin, maxes, mins
  real(kind=8), dimension(nvals) :: mean, ximf, b, rem, cc, dg
  logical                  :: cont_sift     
! make a copy of the input data
  do i=1,nvals
    ximf(i)=indata(i)
  enddo

! loop to decompose the input signal into ncomp successive IMFs
  nimf=0
  do im=1,ncomp
   
! at the beginning of the sifting process, h is the signal
  rem=ximf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! iterate sifting procedure until stopping criteria are satisfied 
  ii = 0
  nsc = 0
  cont_sift = .true.
  do while (cont_sift .eqv. .true.)
    ii = ii+1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i=2,nvals-1
        diff(i)=rem(i+1)-rem(i)
      end do
      diff(1)=0
     
! store extrema without distinction
      nmax=0
      nmin=0
      ik=0
      nzc=0
	  
! test for min and max points, store i in list	
      do i=2,nvals-2
        Dv0=diff(i-1)/dabs(diff(i-1))
        Dv1=diff(i)/dabs(diff(i))
        Dv2=diff(i+1)/dabs(diff(i+1))
        
        if(diff(i) .eq. 0) then ! we are on a zero
          if(dabs(Dv0-Dv1) >  0) then ! and it is a maximum
            ik=ik+1
            maxmin(ik)=i 
            nmax=nmax+1
          end if    
        elseif(dabs(Dv1-Dv2) > 0) then ! we are straddling a zero so
           ik=ik+1
           maxmin(ik)=i+1 ! define zero as at i+1 (not i)
           nmin=nmin+1
        end if
      end do

! find zero crossings RM181128
      do i=1,nvals-1
        if (rem(i) .lt. 0 .and. 0 .lt. rem(i+1) .or. rem(i) .gt. 0 .and. 0 .gt. rem(i+1)) then
          nzc = nzc + 1
        end if
      end do

! test whether we need to continue sifting      
	  if(nmax+nmin < 2) then ! then it is the residue, stop sifting
	    exit
	  endif

! divide maxmin into maxes and mins
	  k1=0
	  k2=0

	  if (maxmin(1) > maxmin(2)) then ! first one is a max not a min
	    do i=1,nmax+nmin
	      j1=1+(i-1)*2
	      j2=2+(i-1)*2
	      if(j1 .le. (nmax + nmin)) then
	        maxes(i)=maxmin(j1)
	        k1=k1+1
	      endif
	      if(j2 .le. (nmax + nmin)) then
	        mins(i)=maxmin(j2)
	        k2=k2+1
	      endif
	    enddo
	  else ! the other way around                               
	    do i=1,nmax+nmin
	      j1=1+(i-1)*2
	      j2=2+(i-1)*2
	      if(j1 .le. (nmax + nmin)) then
	        mins(i)=maxmin(j1)
	        k2=k2+1
	      endif
	      if(j2 .le. (nmax + nmin)) then
	        maxes(i)=maxmin(j2)
	        k1=k1+1
	      endif
	    enddo          
	  endif
	  
! make endpoints both maxes and mins
      do i=1,k1
        maxees(i+1)=maxes(i)
      enddo
      do i=1,k2
        mines(i+1)=mins(i)
      enddo      

      maxees(1)=1
      maxees(k1+2)=nvals
      mines(1)=1
      mines(k2+2)=nvals

! spline interpolate to get upper and lower envelopes and form imf
! list containing i for spline interpolation
	  do i=1,nvals
	    vect(i)=i
	  enddo

! maxenv
      ki=1     
      do i=1,nvals
        if(i .eq. maxees(ki)) then
          hmaxes(ki) = rem(i)
          ki=ki+1
        endif
      enddo
	  
      call spline(maxees, hmaxes, b, cc, dg, k1+2)

      do i=1,nvals
        vects=vect(i)
        maxess=ispline(vects, maxees, hmaxes, b, cc, dg, k1+2)
        maxenv(i)=maxess
      end do
	  
!minenv
      ki=1     
      do i=1,nvals
        if(i .eq. mines(ki)) then
          hmines(ki) = rem(i)
          ki=ki+1
        endif
      enddo

      call spline(mines, hmines, b, cc, dg, k2+2)

      do i=1,nvals
        vects=vect(i)
        miness=ispline(vects,mines, hmines, b, cc, dg, k2+2)
        minenv(i)=miness
      enddo

! print*, "no. of extrema: ", ik,"max:",k1,"min:",k2, "no. of zero crossings: ", nzc  
! calclate the mean of the upper and lower envelope and subtract from remaining signal
      mean=(maxenv+minenv)*0.5
  	  rem=rem-mean

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  	  
! check proto-mode function against stopping criteria    
    if (ik .eq. nzc .or. ik .eq. nzc-1 .or. ik .eq. nzc+1) then
      nsc = nsc+1
      if (nsc .ge. 5 .or. ii .ge. 10) then
        cont_sift = .false.
      end if
    else
      nsc = 0
      if (ii .ge. 10) then
        cont_sift = .false.
      end if
    end if
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

! store the extracted IMF in the result matric
    nimf=nimf+1 
    do i=1,nvals
      rslt(nimf,i)=rem(i)
    enddo
     
! stop criterion of the if we reach the end before n
    if(nmax+nmin < 2) then
      exit
    endif
	
! substract the extracted IMF from the signal
	ximf=ximf-rem
  enddo

  end subroutine emd_v3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine emd_v2(indata, nvals, ncomp, rslt)
	  
! subroutine to decompose a time series using the empirical mode decomposition procedure
! where indata is the time series, nvals is rge length, ncomp is the number of imfs and
! rslt is the matrix containing successive IMFs and the residue
! RM180913 this is an update of the emd subroutine
	  
  implicit none

! passed variables
  integer*4,   intent(in)  :: nvals, ncomp 
  real(kind=8),intent(in)  :: indata(nvals)      
  real(kind=8),intent(out) :: rslt(ncomp,nvals)	

! local variables
  integer*4                :: nmax, nmin, nimf, ki, k1, k2, ik
  integer*4                :: i, j, j1, im, j2, ii  
  real(kind=8)             :: Dv0, Dv1, Dv2
  real(kind=8)             :: maxess, vects, miness
  real(kind=8)             :: ispline
  real(kind=8), dimension(nvals) :: hmaxes, hmines, minenv, maxenv
  real(kind=8), dimension(nvals) :: maxees, mines, vect
  real(kind=8), dimension(nvals) :: diff, maxmin, maxes, mins
  real(kind=8), dimension(nvals) :: mean, ximf, b, rem, cc, dg
        
! make a copy of the input data
  do i=1,nvals
    ximf(i)=indata(i)
  enddo

! loop to decompose the input signal into ncomp successive IMFs
  nimf=0
  do im=1,ncomp
   
! at the beginning of the sifting process, h is the signal
  rem=ximf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! iterate sifting procedure for an s number of 5  
    do ii=1,5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i=2,nvals-1
        diff(i)=rem(i+1)-rem(i)
      end do
      diff(1)=0
     
! store extrema without distinction
      nmax=0
      nmin=0
      ik=0
	  
! test for min and max points, store i in list	
      do i=2,nvals-2
        Dv0=diff(i-1)/dabs(diff(i-1))
        Dv1=diff(i)/dabs(diff(i))
        Dv2=diff(i+1)/dabs(diff(i+1))
        
        if(diff(i) .eq. 0) then ! we are on a zero
          if(dabs(Dv0-Dv1) >  0) then ! and it is a maximum
            ik=ik+1
            maxmin(ik)=i 
            nmax=nmax+1
          end if    
        elseif(dabs(Dv1-Dv2) > 0) then ! we are straddling a zero so
           ik=ik+1
           maxmin(ik)=i+1 ! define zero as at i+1 (not i)
           nmin=nmin+1
        end if
      end do

! test whether we need to continue sifting      
	  if(nmax+nmin < 2) then ! then it is the residue, stop sifting
	    exit
	  endif

! divide maxmin into maxes and mins
	  k1=0
	  k2=0

	  if (maxmin(1) > maxmin(2)) then ! first one is a max not a min
	    do i=1,nmax+nmin
	      j1=1+(i-1)*2
	      j2=2+(i-1)*2
	      if(j1 .le. (nmax + nmin)) then
	        maxes(i)=maxmin(j1)
	        k1=k1+1
	      endif
	      if(j2 .le. (nmax + nmin)) then
	        mins(i)=maxmin(j2)
	        k2=k2+1
	      endif
	    enddo
	  else ! the other way around                               
	    do i=1,nmax+nmin
	      j1=1+(i-1)*2
	      j2=2+(i-1)*2
	      if(j1 .le. (nmax + nmin)) then
	        mins(i)=maxmin(j1)
	        k2=k2+1
	      endif
	      if(j2 .le. (nmax + nmin)) then
	        maxes(i)=maxmin(j2)
	        k1=k1+1
	      endif
	    enddo          
	  endif
	  
! make endpoints both maxes and mins
      do i=1,k1
        maxees(i+1)=maxes(i)
      enddo
      do i=1,k2
        mines(i+1)=mins(i)
      enddo      

      maxees(1)=1
      maxees(k1+2)=nvals
      mines(1)=1
      mines(k2+2)=nvals

! spline interpolate to get upper and lower envelopes and form imf
! list containing i for spline interpolation
	  do i=1,nvals
	    vect(i)=i
	  enddo

! maxenv
      ki=1     
      do i=1,nvals
        if(i .eq. maxees(ki)) then
          hmaxes(ki) = rem(i)
          ki=ki+1
        endif
      enddo
	  
      call spline(maxees, hmaxes, b, cc, dg, k1+2)

      do i=1,nvals
        vects=vect(i)
        maxess=ispline(vects, maxees, hmaxes, b, cc, dg, k1+2)
        maxenv(i)=maxess
      end do
	  
!minenv
      ki=1     
      do i=1,nvals
        if(i .eq. mines(ki)) then
          hmines(ki) = rem(i)
          ki=ki+1
        endif
      enddo

      call spline(mines, hmines, b, cc, dg, k2+2)

      do i=1,nvals
        vects=vect(i)
        miness=ispline(vects,mines, hmines, b, cc, dg, k2+2)
        minenv(i)=miness
      enddo
	  
! calclate the mean of the upper and lower envelope and subtract from remaining signal
      mean=(maxenv+minenv)*0.5
  	  rem=rem-mean

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  	  
	enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

! store the extracted IMF in the result matric
    nimf=nimf+1 
    do i=1,nvals
      rslt(nimf,i)=rem(i)
    enddo
     
! stop criterion of the if we reach the end before n
    if(nmax+nmin < 2) then
      exit
    endif
	
! substract the extracted IMF from the signal
	ximf=ximf-rem
  enddo

  end subroutine emd_v2

! ================================================
!module mysubs
!	implicit none 
! contains
 subroutine emd( xt, n, k, imf)

implicit none

        integer*4, intent(in) :: n, k 
        integer*4  :: kmax, kmin, ki, kimf, k1, k2, ik
	double precision, dimension(n), intent(in) :: xt(n)
        double precision :: eps, d11, d00
        double precision :: SD, d1, Dv0, Dv1, Dv2
        double precision  maxess, vects, miness
        double precision :: ispline
        integer*4 i, j, t, j2, count  
        double precision, dimension(n) ::   hmaxes, hmines, minenv
        double precision, dimension(n) :: maxees, mines, vect
        double precision, dimension(n) :: d, maxmin, maxes, mins
        double precision, dimension(n) :: m, c, b, h, cc, dg
        double precision, dimension(n) :: maxenv, prevh
        double precision, dimension(k,n), intent(out) :: imf(k,n)



! x : input time series
! n : length of input time series 
! k : number of Independent Mode Function (IMF)

! to test

!-------------------------------------------------------------------------
! loop to decompose the input signal into n successive IMFs
!imf :  Matrix which will contain the successive IMF, and the residue (last IMF)

  do i=1,n
    c(i) = xt(i)
  end do

  kimf=0
  do t=1,k

! loop on successive IMFs   
!   %-------------------------------------------------------------------------
!   % inner loop to find each imf
   
! at the beginning of the sifting process, h is the signal

    do i=1,n
       h(i) = c(i)
    end do
 
! Standard deviation which will be used to stop the sifting process
    SD = 1; 
    
!initialise arrays
    d=0
    maxmin=0
    maxes=0
    mins=0
    maxees=0
    mines=0
    m=0

! while the standard deviation is higher than 0.3 (typical value)
    count =0
    do  while (SD > 0.1) ! **
      
      do i=2,n-1
        d1 = h(i+1)-h(i)
        d(i) = d1
      end do
      d(1)=0 !** JP
     
 ! to store the optima (min and max without distinction so far)
      kmax = 0
      kmin = 0
      ik=0
      do i=2,n-2



        Dv0 = d(i-1)/dabs(d(i-1))
        Dv1 = d(i)/dabs(d(i))
        Dv2 = d(i+1)/dabs(d(i+1))
        
        if(d(i) .eq. 0) then
         ! we are on a zero
         ! it is a maximum
          if (dabs(Dv0-Dv1) >  0)  then 
            ik=ik+1
            maxmin(ik) = i 
            kmax = kmax+1
          end if
       ! we are straddling a zero so
            
        elseif( dabs(Dv1-Dv2) > 0 )  then
           ik=ik+1
           maxmin(ik) =  i+1      ! define zero as at i+1 (not i)
           kmin = kmin+1
        end if

      end do


      ! then it is the residue
      if(kmax+kmin < 2) then
!        print*,'exiting kmax,kmin:',kmax,kmin
        exit
      endif

    ! divide maxmin into maxes and mins
      k1=0
      k2=0

      if (maxmin(1) > maxmin(2)) then
        ! first one is a max not a min
         
        do i=1,kmax+kmin
          j=1+(i-1)*2
          j2 = 2+(i-1)*2
          if(j .le. (kmax + kmin)) then
            maxes(i) = maxmin(j)
            k1=k1+1
          end if
          if(j2 .le. (kmax + kmin)) then
            mins(i) = maxmin(j2)
            k2=k2+1
          end if

        end do
         !mins  = maxmin(2:2:length(maxmin));
       ! the other way around
      else                                
         !maxes = maxmin(2:2:length(maxmin));
         !mins  = maxmin(1:2:length(maxmin));
        do i=1,kmax+kmin
          j=1+(i-1)*2
          j2 = 2+(i-1)*2
          if(j .le. (kmax + kmin)) then
            mins(i) = maxmin(j)
            k2=k2+1
          end if
          if(j2 .le. (kmax + kmin)) then
            maxes(i) = maxmin(j2)
            k1=k1+1
          end if
        end do          

      end if

      ! make endpoints both maxes and mins
      
      do i=1,k1
        maxees(i+1) = maxes(i)
      end do
      do i=1,k2
        mines(i+1)  = mins(i)
      end do      

      maxees(1) = 1
      maxees(k1+2) = n
      mines(1) = 1
      mines(k2+2) = n

      do i=1,n
        vect(i) = i
      end do
      !-------------------------------------------------------------------------
      ! spline interpolate to get max and min envelopes; form imf
      ki=1     
      do i=1,n

        if(i .eq. maxees(ki)) then
          hmaxes(ki) = h(i)
          ki=ki+1
        end if

      end do

      ki=1     
      do i=1,n

        if(i .eq. mines(ki)) then
          hmines(ki) = h(i)
          ki=ki+1
        end if

      end do
      
      ! maxenv
      call spline (maxees, hmaxes, b, cc, dg, k1+2)

      do i=1,n
        vects = vect(i)
        maxess =  ispline(vects, maxees, hmaxes, b, cc, dg, k1+2)
        maxenv(i) = maxess
      end do

      !minenv
      
      call spline (mines, hmines, b, cc, dg, k2+2)

      do i=1,n
        vects = vect(i)
        miness = ispline(vects,mines, hmines, b, cc, dg, k2+2)
        minenv(i) = miness
      end do

      do i =1,n
        d11 = 0.0
         !  mean of max and min enveloppes
        d11 = maxenv(i) + minenv(i)
        d11 = d11/2
        m(i) = d11 
      end do

      do i = 1,n
        prevh(i) = h(i) ! copy of the previous value of h before modifying it
        h(i) = h(i) - m(i) ! substract mean to h
      end do

      ! calculate standard deviation
      eps = 0.0000001 ! to avoid zero values

      SD = 0
      do i = 1,n
        d1 =  (prevh(i) - h(i))
        d1 = d1*d1  
        d1 = d1/(prevh(i)*prevh(i) + eps)    
        SD = SD +  d1 
      end do
    
      count = count+1

       ! It is more judicious to use count >1000  **JP

      if(count>10000) then
        print*,'Error A- loop error - too many iterations (iter > 1000)'
        exit
      end if

    end do ! 1**

! store the extracted IMF in the matrix imf

    kimf=kimf+1 
    do i=1,n
      imf(kimf,i) = h(i)
    end do
   
! if size(maxmin,2)<2, then h is the residue
     
! stop criterion of the algo. if we reach the end before n
    if(kmax+kmin < 2) then
!size(maxmin,2) < 2
!      print*,' Error B - see kmax criterion'
      exit
    end if
   
    do i= 1,n
! substract the extracted IMF from the signal
      d00 = c(i)
      d11 = c(i) - h(i)
      c(i) = d11 
   ! write(30,*) d00, h(i), d11
    end do

  END DO

  end subroutine emd
! *************************************************************************



! *************************************************************************
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
  integer*4 n
  double precision :: x(n), y(n), b(n), c(n), d(n)
  integer*4 i, j, gap
  double precision :: h

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

!end module mysubs

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
  double precision :: ispline
  integer*4 n
  double precision ::  u, x(n), y(n), b(n), c(n), d(n)
  integer*4 i, j, k
  double precision :: dx

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




!_____________________________________________________________________


! =============
subroutine init_random_seed()
            implicit none
            integer, allocatable :: seed(:)
            integer :: i, n, un, istat, dt(8), pid, t(2), s
            integer(8) :: count, tms
          
            call random_seed(size = n)
            allocate(seed(n))
            ! First try if the OS provides a random number generator
            open(newunit=un, file="/dev/urandom", access="stream", &
                 form="unformatted", action="read", status="old", iostat=istat)
            if (istat == 0) then
               read(un) seed
               close(un)
            else
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
               call system_clock(count)
               if (count /= 0) then
                  t = transfer(count, t)
               else
                  call date_and_time(values=dt)
                  tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
                  t = transfer(tms, t)
               end if
               s = ieor(t(1), t(2))
               pid = getpid() + 1099279 ! Add a prime
               s = ieor(s, pid)
               if (n >= 3) then
                  seed(1) = t(1) + 36269
                  seed(2) = t(2) + 72551
                  seed(3) = pid
                  if (n > 3) then
                     seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
                  end if
               else
                  seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
               end if
            end if
            call random_seed(put=seed)
          end subroutine init_random_seed

