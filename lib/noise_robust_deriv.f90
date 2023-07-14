  subroutine noise_robust_deriv(data, deriv, delta, length, n)

! Tony Purcell April 12th 2013
! Subroutine to implement a noise-reducing numerical differentiation on
! a vector of function values with "length" number of entries.

! We have chosen a noise-robust technique that is smooth up to order 4 polynomial entries (x^4).
!
! The technique used is that described derived by Pavel Holoborodko
! http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/

    implicit none

! passed arguments    
    integer, intent(in)          :: length          ! number of entries in the input data vector
    integer, intent(in)          :: n               ! number of points to use in the filter. Must be an odd number
    real (kind=8), intent(in)    :: data(length)    ! input data to be time-differentiated
    real (kind=8), intent(out)   :: deriv(length)   ! output time-differentiated vector
    real (kind=8), intent(in)    :: delta           ! time interval between data points

! local counters
    integer :: i, j, n1, n2, n3

! coefficients of the filter
    real (kind=8), dimension(11) :: div
    real (kind=8), dimension(11, 5) :: c
    real (kind=8), dimension(5) :: d
    real (kind=8) :: div_d

! check the validity of the input parameters passed into the subroutine
    if ( ( n .gt. 11 ) .or. ( n .lt. 5 ) .or. ( mod(n,2) .ne. 1 ) ) then
      write(6, *) "ERROR: subroutine noise_robust_deriv - illegal value for N: ", n
      stop
    endif

    if ( length .lt. 1 ) then
      write(6, *) "ERROR: subroutine noise_robust_deriv - illegal value for length: ", length
      stop
    endif

    if ( delta .lt. 1.0d-16 ) then
      write(6, *) "ERROR: subroutine noise_robust_deriv - illegal value for delta: ", delta
      stop
    endif

! note that five point technique is only smooth up to x^2 terms
    c(5,1) = 2.0d0
    c(5,2) = 1.0d0
    div(5) = 8.0d0

! Here we use noise robust coefficients for the derivative series
    c(7,1) = 39.0d0
    c(7,2) = 12.0d0
    c(7,3) = -5.0d0
    div(7) = 96.0d0

    c(9,1) = 27.0d0
    c(9,2) = 16.0d0
    c(9,3) = -1.0d0
    c(9,4) = -2.0d0
    div(9) = 96.0d0
   
    c(11, 1) = 322.0d0
    c(11, 2) = 256.0d0
    c(11, 3) = 39.0d0
    c(11, 4) = -32.0d0
    c(11, 5) = -11.0d0
    div(11) = 1536.0d0

    d(1) = 1.0d0
    d(2) = 2.0d0
    d(3) = 0.0d0
    d(4) = -2.0d0
    d(5) = -1.0d0
    div_d = 8.0d0

! Calculate derivative over the interior of the input vector (i.e. away from the edges of the vector)
    n2 = (n - 1)/2
    do i = n2 + 1, length - n2 
      deriv(i) = 0.0d0
      do j = 1, n2
        deriv(i) = deriv(i) + c(n, j) * (data(i + j) - data(i - j)) 
      enddo
      deriv(i) = deriv(i)/(delta * div(n))
    enddo
 
! Use lower order formulations for points closer to the edges of the data set
    n1 = n
    n3 = n2
    do while ( n1 .gt. 5 )
      n1 = n1 - 2
      n3 = n3 - 1

      i=n3+1
      deriv(i) = 0.0d0
      do j = 1, n3
        deriv(i) = deriv(i) + c(n1, j) * (data(i + j) - data(i - j)) 
      enddo
      deriv(i) = deriv(i)/(delta * div(n1))

      i = length - n3
      deriv(i) = 0.0d0
      do j = 1, n3
        deriv(i) = deriv(i) + c(n1, j) * (data(i + j) - data(i - j)) 
      enddo
      deriv(i) = deriv(i)/(delta * div(n1))

    enddo

! Use one sided five point formulae for points 1, 2, n - 1, and n
    do i = 1, 2
      deriv(i) = 0.0d0
      do j = 1, 5
        deriv(i) = deriv(i) + d(j) * data(i + j - 1) 
      enddo
      deriv(i) = deriv(i)/(div_d * delta)
    enddo

    do i = length - 1, length
      deriv(i) = 0.0d0
      do j = 1, 5
        deriv(i) = deriv(i) - d(j) * data(i + 1 - j) 
      enddo
      deriv(i) = deriv(i)/(div_d * delta)
    enddo

    return
    end subroutine noise_robust_deriv
