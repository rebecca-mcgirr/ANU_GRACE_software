    subroutine legcalc_norm(X)

! calculate the normalised associated Legendre function values at x up to degree ndeg

! The normalisation is such that the integral of Pnm^2 from -1 to 1 is 2.
! This is slightly different from the normalisation used in the older
! sea-level programs.  In those programs, if m is not equal to 0, the
! normalisation had an extra factor sqrt(2), so that integral over the
! sphere of Pnm cos or sin m phi is 4 pi.  With this normalisation
! the integral over the sphere of Pnm exp(i m phi) is equal to 4 pi.
!
! Modified from Lambeck/Johnston/Zhang code that was part of the calsea program
!
! EKP091215: renormalised legendre polynomials (4pi) to be consistent with GRACE data
!            and added the calculation of the derivatives for Plm = PDlm
      
  use spherhar_mod
  implicit none

! input co-latitude in radians
  real*8 X

! local variables
  real*8 aa,as,s,sum,c,d
  integer i,j,l,m,l0,i0,i1,i2,i3,i4
  integer ntot

    Plm = 0.0d0
    PDlm = 0.0d0

    Plm(1)=1
    Plm(2)=X

    i1=1
    i2=2
! Evaluate the m=0 Legendre functions (un-normalised)
    do i=2,ndeg
      i0=i1
      i1=i2
      i2=i*(i+1)/2+1
      Plm(i2)=((2*i-1)*X*Plm(i1)-(i-1)*Plm(i0))/i 
       if(i2.eq.2)print*,'Plm(2)',Plm(2)
    enddo
! Multiply by a factor that normalises to 4pi
    do i=1,ndeg
      i2=i*(i+1)/2+1
      Plm(i2)=dsqrt(2.d0*i+1)*Plm(i2)
    enddo

    AA=1.d0-X*X
    IF(AA.eq.0.d0) then
      DO 200 l=0,ndeg
        l0=l*(l+1)/2+1
        do 200 m=1,l
200       Plm(m+l0)=0.d0
          RETURN    
    endif
    AS=dSQRT(AA)
    Plm(3)=AS*dSQRT(3.D0)
    Plm(5)=3.d0*AS*X*dSQRT(5.D0/3.d0)
    Plm(6)=3.d0*AA*dSQRT(5.D0/12.d0)


    DO l=3,ndeg
      SUM=1.d0
      DO  J=1,l
        S=dble(2*J-1)/(2.d0*J)
        SUM=SUM*dSQRT(AA*S) 
      enddo
      l0=l*(l+1)/2+1

      Plm(l+l0)=SUM*dSQRT(4.D0*l+2.d0)
      SUM=1.d0

      DO  J=2,l
        S=dble(2*J-1)/(2*J-2.d0)
        SUM=SUM*dSQRT(AA*S)
      enddo

      plm(l-1+l0)=X*dSQRT(4.D0*l+2.d0)*SUM
 
      DO  M=l-2,1,-1
        C=2.d0*(M+1)/dSQRT((l+M+1.D0)*(l-M))
        C=C*X/AS
        D=(l+M+2)*(l-M-1)/((l+M+1.D0)*(l-M))
        D=dSQRT(D)
        Plm(M+l0)=C*Plm(M+1+l0)-D*Plm(M+2+l0)  
       if(M+l0.eq.2)print*,'Plm(2)',Plm(2)
      enddo

    enddo
!do i=56,66
!print*, "Plm",i,Plm(i)
!enddo
!stop
! The following loops calculate the derivative of the normalised 
! associated legendre function based on the recursive relationships below:
! 
! for m=0
! Pn0' = 1/(sqrt(1-x^2)) * 1/sqrt(2) * sqrt(n*(n+1)) * Pn1
!
! for m=1
! Pn1' = x/(1-x^2)*Pn1 - sqrt(n*(1+n))/sqrt(1-x^2)*sqrt(2)*Pn0  
! 
! for m>1
! Pmn' = mx/(1-x^2)*Pnm -sqrt((m+n)*(n-m+1))/sqrt(1-x^2)*Pm-1,n

    DO i = 2,ndeg
      i2 = i*(i+1)/2+1
      PDlm(i2) = 1/AS*dsqrt(i*(i+1)/2.d0)*Plm(i2+1)
      PDlm(i2+1) = X/AA*Plm(i2+1)-dsqrt(2*i*(i+1)/AA)*Plm(i2)
     
      i4 = (i+1)*(i+2)/2
      m = 2


      DO j = i2+2 , i4

        PDlm(j) = m*X/AA*Plm(j) - dsqrt((m+i)*(i-m+1)/AA)*Plm(j-1)
        m=m+1
      ENDDO
    ENDDO

      return

      end

