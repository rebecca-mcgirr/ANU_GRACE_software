    subroutine pnm_calc(x, pnm, degree, max_deg, flag)
!    subroutine pnm_calc(x, pnm, dpnm, degree)

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! calculate the normalised associated Legendre function values at x up to degree
!
! Modified from Lambeck/Johnston/Zhang code that was part of the calsea program
!
! APP130919: Incorporated more extensive documentation and changed pnm to be two a dimensional array
!            Commented out code for calculating derivative terms
!            Detailed discussion of the normalisation used is given at the end of the code
!            The non-standard normalisation in turn affects the recursion relations between the
!            Legendre polynomials and derivations of these formula are also included.
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
  implicit none

!  integer*4, parameter:: max_deg = 512      ! Max degree of spherical harmonic expansion

  integer*4 degree, flag, max_deg

  real*8 x, pnm, dpnm

  dimension pnm(0:max_deg, 0:max_deg)
  dimension dpnm(0:max_deg, 0:max_deg)

  real*8 aa, as, s, sum, c, d

  integer i, j, l, m, n


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Initialise arrays, calculate P00 and P10, then apply recursion formula for Legendre polynomials (Spiegel - equ'n 25.20 with i = n + 1)
  pnm = 0.0d0
  dpnm = 0.0d0

  pnm(0, 0) = 1
  pnm(1, 0) = x
  do i = 2, degree
! PT/APP 160218: fixed indexing in this equation
    pnm(i, 0) = (dble(2*i - 1) * x * pnm(i-1, 0) - dble(i - 1) * pnm(i-2, 0))/dble(i)
  enddo                 ! end of i loop
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     do l=0,lmax-2    !First column of the matrix (m=0)
!        plm(l+2,0)=((2*(l+1)+1)*X*plm(l+1,0)-(l+1)*plm(l,0))/(l+2)
!     enddo

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Multiply by a factor that normalises to 4pi
  do i = 1, degree
    pnm(i, 0) = dsqrt(dble(2*i + 1)) * pnm(i, 0)
  enddo                 ! end of i loop
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  aa = 1.d0 - x*x
  as = dsqrt(aa)

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! If at either pole set all associated legendre polynomials to zero
  if ( aa .eq. 0.d0 ) then
    do n = 0, degree
      do m = 1, n
        pnm(n, m) = 0.d0
      enddo             ! end of m loop
    enddo               ! end of n loop
    return    
  endif
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Use formulae for degree 1 and 2 Associated legendre polynomials and apply normalisation (Spiegel - equ'ns 26.5-7 & pg 150 equ'n 26.15)
  pnm(1, 1) = as * dsqrt(3.D0)
  pnm(2, 1) = 3.d0 * as * x * dsqrt(5.D0/3.d0)
  pnm(2, 2) = 3.d0 * aa * dsqrt(5.D0/12.d0)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  do n = 3, degree

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Pnn(x) = (1 - x^2)^n/2 * 1/sqrt(2^n) * 1/sqrt(n!) * sqrt((2n - 1) * (2n - 3) * .... 1) * sqrt(4n + 2) (Derivation is given below)
    sum = 1.0d0
    do j = 1, n
      s = dble(2*j - 1)/dble(2 * j)
      sum = sum * dsqrt(aa * s) 
    enddo               ! end of j loop
    pnm(n, n) = sum * dsqrt(dble(4*n + 2))
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Pn(n-1)(x) = (1 - x^2)^(n-1)/2 * x * 1/sqrt(2^(n-1)) * 1/sqrt((n-1)!) * sqrt((2n - 1) * (2n - 3) * .... 1) * sqrt(4n + 2)
! Derivation is given below
    sum = 1.0d0
    do  j = 2, n
      s = dble(2*j - 1)/dble(2*j - 2)
      sum = sum * dsqrt(aa * s)
    enddo               ! end of j loop
    pnm(n, n-1) = x * dsqrt(dble(4*n + 2)) * sum
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Apply recursion formula for remaining associated Legendre polynomials (formula derived below)
    do  m = n-2, 1, -1
      c = 2.d0 * dble(m + 1)/dsqrt(dble((n + m + 1) * (n - m)))
      c = c * x/as
      d = dble((n + m + 2) * (n - m - 1))/dble((n + m + 1) * (n - m))
      d = dsqrt(d)
      pnm(n, m) = c * pnm(n, m + 1) - d * pnm(n, m + 2)  
    enddo               ! end of m loop
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  enddo                 ! end of n loop

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Calculate the derivative of the normalised spherical harmonic functiona, start with Pn0 and Pn1 then apply recursion
  if ( flag .eq. 1 ) then
    do n = 2, degree
      dpnm(n, 0) = 1/as * dsqrt(dble(n * (n + 1))/2.d0) * pnm(n, 1)
      dpnm(n, 1) = x/aa * pnm(n, 1) - dsqrt(dble(2 * n * (n + 1))/aa) * pnm(n, 0)
      do m = 2, m
        dpnm(n, m) = m * x/aa * pnm(n, m) - dsqrt(dble((m + n) * (n - m + 1))/aa) * pnm(n, m - 1)
      enddo               ! end of m loop
    enddo                 ! end of n loop
    pnm = dpnm
  endif
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  return
  end

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! NOTES ON NORMALISATION:
!
! There are a number of choices of normalisation that can be used for the Legendre polynomials. We will discuss the general principles below.
!
! We define:
!
! Snm = integral([-1,1]) Pnm(x) Pnm(x) dx) 
!
! If we use the unscaled forms for the Legendre polynomials from Spiegel - A Mathematical Handbook:
!
! Pn0(x) = 2^(-n)/(n!) d^n/dx^n (x^2 - 1)^n                 pg 146, equ'n 25.2
!
! Pnm(x) = (1-x^2)^(m/2) d^m/dx^m Pn0(x)                    pg 149, equ'n 26.2
!
! Then we have:
!
! Sn0 = 2/(2n + 1)                                          pg 147, equ'n 25.26
!
! Snm = 2/(2n + 1) (n + m)!/(n - m)!                        pg 150, equ'n 26.15
!
! These are the basic values for the Legendre polynomial, we are free to scale our polynomials to obtain values that are easier to work with:
!
! From basic trigonometric identities we have:
!
! Rm  = integral([0,2pi]) sin(my) sin(my) dy)
!     = integal([0,2pi]) 1/2 (1 - cos(2my)) dy)
!     = 1/2 (2pi - 0) - 1/2 1/2m (sin(4m pi) - sim(0))
!     = pi 
!
! Qm  = integral([0,2pi]) cos(my) cos(my) dy)
!     = integal([0,2pi]) 1/2 (1 + cos(2my)) dy)
!     = 1/2 (2pi - 0) + 1/2 1/2m (sin(4m pi) - sim(0))
!     = pi 
!
! where we have assumed m > 0. For m = 0 these identities become: R0 = 0, Q0 = 2 pi
!
! In this application we require, for reasons of mathematical convenience:
!
! 4 pi = integral([-1,1], integral([0,2pi], Pnm(x) Pnm(x) trig(my) trig(my) dy) dx) 
!      = integral([-1,1], Pnm(x) Pnm(x) dx) integral([0,2pi], trig(my) trig(my) dy)
!      = Snm pi
!
! In other words we want to normalise our associated Legendre polynomials so that:
!
! integral([-1,1]) Pnm(x) Pnm(x) dx) = 4            if m > 0
! integral([-1,1]) Pnm(x) Pnm(x) dx) = 2            if m = 0
!
! Letting Ynm(x) denote the unscaled form of the Legendre polynomials given above our forms for the Legendre polynomials become:
!
! Pn0(x) = Yn0(x) sqrt(2n + 1)
!
! Pnm(x) = Ynm(x) sqrt((4n + 2) (n - m)!/(n + m)!)
!
! NOTE: In previous applications of this code in the calsea program we required Snm = 2. The scaling factors used here therefore differ from
!       the scaling factors used there by a factor of sqrt(2)
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! SPECIFIC LEGENDRE POLYNOMIALS:
!
! Explicit forms for Legendre polynomials of low degree are given in Spiegel equ'ns 25.3-10 and 26.5-10
!
! Explicit forms for two special cases may also be obtained:
!
! ---------------------------------------------------------------------------------------------------------------------------------------------
!
! Pnn(x) = (1 - x*x)^(n/2) * 1/(2^n) * 1/n! * d^(2n)/dx^(2n) (x^2 - 1)^n * normalisation_term
!        = (1 - x*x)^(n/2) * 1/(2^n) * 1/n! * (2n)! * normalisation_term
!
! Applying normalisation factor of sqrt((4n + 2)/(2n)!) gives:
!
! Pnn(x) = (1 - x*x)^n/2 * 1/(2^n) * 1/n! * (2n)! * sqrt(4n + 2)/sqrt((2n)!)
!        = (1 - x*x)^n/2 * 1/(2^n) * 1/n! * sqrt((2n)!) * sqrt(4n + 2)
!        = (1 - x*x)^n/2 * 1/(2^n) * 1/n! * sqrt(2n * (2n - 2) * .... * 2) * sqrt((2n - 1) * (2n - 3) * .... 1) * sqrt(4n + 2)
!        = (1 - x*x)^n/2 * 1/sqrt(2^n) * 1/n! * sqrt(n * (n - 1) * .... * 1) * sqrt((2n - 1) * (2n - 3) * .... 1) * sqrt(4n + 2)
!        = (1 - x*x)^n/2 * 1/sqrt(2^n) * 1/n! * sqrt(n!) * sqrt((2n - 1) * (2n - 3) * .... 1) * sqrt(4n + 2)
!        = (1 - x*x)^n/2 * 1/sqrt(2^n) * 1/sqrt(n!) * sqrt((2n - 1) * (2n - 3) * .... 1) * sqrt(4n + 2)
!
! ---------------------------------------------------------------------------------------------------------------------------------------------
!
! Pn(n-1)(x) = (1 - x*x)^(n-1)/2 * 1/(2^n) * 1/n! * d^(2n-1)/dx^(2n-1) (x^2 - 1)^n * normalisation_term
!            = (1 - x*x)^(n-1)/2 * 1/(2^n) * 1/n! * (2n)! * x * normalisation_term
!
! Applying normalisation factor of sqrt((4n + 2)/(2n-1)!) gives:
!
! Pn(n-1)(x) = (1 - x*x)^(n-1)/2 * x * 1/(2^n) * 1/n! * (2n)! * sqrt(4n + 2)/sqrt((2n - 1)!)
!            = (1 - x*x)^(n-1)/2 * x * 1/(2^n) * 1/n! * (2n) * sqrt((2n - 1)!) * sqrt(4n + 2)
!            = (1 - x*x)^(n-1)/2 * x * 1/(2^(n-1)) * 1/(n-1)! * sqrt((2n - 2) * .... * 2) * sqrt((2n - 1) * (2n - 3) * .... 1) * sqrt(4n + 2)
!            = (1 - x*x)^(n-1)/2 * x * 1/sqrt(2^(n-1)) * 1/(n-1)! * sqrt((n - 1) * .... * 1) * sqrt((2n - 1) * (2n - 3) * .... 1) * sqrt(4n + 2)
!            = (1 - x*x)^(n-1)/2 * x * 1/sqrt(2^(n-1)) * 1/(n-1)! * sqrt((n-1)!) * sqrt((2n - 1) * (2n - 3) * .... 1) * sqrt(4n + 2)
!            = (1 - x*x)^(n-1)/2 * x * 1/sqrt(2^(n-1)) * 1/sqrt((n-1)!) * sqrt((2n - 1) * (2n - 3) * .... 1) * sqrt(4n + 2)
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! RECURSION FORMULAE:
!
! Spiegel gives explicit recursion formulae for Legendre polynomials:
!
! Yn0(x) = ((2n - 1) x Y(n-1)0(x) - (n - 1) Y(n-2)0(x)) / n                                     equ'n 25.20
!
! Ynm(x) = (2(m + 1) x/sqrt(1 - x^2) Yn(m+1)(x) - Yn(m+2)(x)) / ((n - m)(n + m + 1))            equ'n 26.13
!
! The first equation may be applied directly since we do the recursion before renormalising. For the associated polynomials that isn't the case.
!
! We can use the form for the normalisation scaling factors to relate them:
!
! sqrt((4n + 2) (n - m)!/(n + m)!) = sqrt((4n + 2) (n - m) (n + m + 1) (n - m - 1)!/(n + m + 1)!)
!                                  = sqrt((4n + 2) (n - m - 1)!/(n + m + 1)!) sqrt((n - m) (n + m + 1))
!                                  = sqrt((4n + 2) (n - m - 1) (n + m + 2) (n - m - 2)!/(n + m + 2)!) sqrt((n - m) (n + m + 1))
!                                  = sqrt((4n + 2) (n - m - 2)!/(n + m + 2)!) sqrt((n - m - 1) (n + m + 2)) sqrt((n - m) (n + m + 1))
!
! Setting Pnm(x) = Cnm Ynm(x) we can then write:
!
! Cnm = Cn(m+1) sqrt((n - m) (n + m + 1))
!     = Cn(m+2) sqrt((n - m) (n + m + 1)) sqrt((n - m - 1) (n + m + 2))
!
! Substituting back into the recursion formula given above yields:
!
! Pnm(x) = (2(m + 1) x/sqrt(1 - x^2) Pn(m+1)(x) - sqrt((n - m - 1) (n + m + 2)) Pn(m+2)(x)) / sqrt((n - m)(n + m + 1))
! 
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! RECURSION FORMULAE FOR DERIVATIVES:
!
! From the definition of the Associated Legendre Polynomials:
!
! Pn1(x) = sqrt(1 - x^2) sqrt(2) sqrt((n - 1)!/(n + 1)!) d/dx Pn0(x)
!        = sqrt(1 - x^2) sqrt(2) sqrt(1/(n(n+1))) d/dx Pn0(x)
!
! Which can be rearranged to give:
!
! d/dx Pn0(x) = 1/sqrt(1-x^2) * 1/sqrt(2) * sqrt(n*(n+1)) * Pn1(x)
!
! differentiating the definition of Pn1 gives:
!
! d/dx Pn1(x) = d/dx (sqrt(1 - x^2) sqrt(2) sqrt((n - 1)!/(n + 1)!) d/dx Pn0(x))
!             = 1/2 -2x/sqrt(1 - x^2) sqrt(2) sqrt(1/(n(n + 1))) d/dx Pn0(x)  
!                               + sqrt(1 - x^2) sqrt(2) sqrt(1/(n(n + 1))) d/dx d/dx Pn0(x)
!             = -x/(1 - x^2) sqrt(1 - x^2) sqrt(2) sqrt((n - 1)!/(n + 1)!) d/dx Pn0(x)) 
!                               + sqrt(1 - x^2) sqrt(2) sqrt(1/(n(n + 1)) (2 x d/dx Pn0(x) - n(n+1) Pn0(x))/(1 - x^2)
!             = -x/(1 - x^2) Pn1(x) + 2x/(1 - x^2) Pn1(x) - sqrt(2) sqrt(n(n + 1)) Pn0(x)/sqrt(1 - x^2)
!             = x/(1 - x^2) Pn1(x) - sqrt(2) sqrt(n(n + 1)) Pn0(x)/sqrt(1 - x^2)
!
! Where we have used the fact that Pn0(x) satisfies the Legendre Equation:
!
! d^2/dx^2 Pn0(x) = (2 x d/dx Pn0(x) - n(n+1) Pn0(x))/(1 - x^2)
!
! for m>1
!
! d/dx Ynm(x) = d/dx ((1 - x^2)^(m/2) d^m/dx^m Yn0(x))
!             = m/2 -2x (1 - x^2)^(m/2 - 1) d^m/dx^m Yn0(x) + (1 - x^2)^(m/2) d^(m+1)/dx^(m+1) Yn0(x)
!             = -mx/(1 - x^2) (1 - x^2)^(m/2) d^m/dx^m Pn0(x)) + 1/sqrt(1 - x^2) (1 - x^2)^((m+1)/2) d^(m+1)/dx^(m+1) Pn0(x)
!             = -mx/(1 - x^2) Ynm(x) + 1/sqrt(1 - x^2) Yn(m+1)(x)
!             = -mx/(1 - x^2) Ynm(x) + 2m x/(1 - x^2) Ynm(x) - (n + m)(n - m + 1)/sqrt(1 - x^2) Yn(m-1)(x)
!             = mx/(1 - x^2) Ynm(x) - (n + m)(n - m + 1) Yn(m-1)(x)/sqrt(1 - x^2)
!
! Where we have used the recursion formula (Spiegel equ'n 26.13)
!
! d/dx Pnm(x) = Cnm d/dx Ynm(x)
!             = Cnm (mx/(1 - x^2) Ynm(x) - (n + m)(n - m + 1) Yn(m-1)(x)/sqrt(1 - x^2))
!             = mx/(1 - x^2) Pnm(x) - Cnm/Cn(m-1) (n + m)(n - m + 1) Pn(m-1)(x)/sqrt(1 - x^2)
!             = mx/(1 - x^2) Pnm(x) - (n + m)(n - m + 1)/sqrt((n + m)(n - m + 1)) Pn(m-1)(x)/sqrt(1 - x^2)
!             = mx/(1 - x^2) Pnm(x) - sqrt((n + m)(n - m + 1)) Pn(m-1)(x)/sqrt(1 - x^2)
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
