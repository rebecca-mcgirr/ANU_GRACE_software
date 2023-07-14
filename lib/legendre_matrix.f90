!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine legendre_matrix(lmax,X,plm)
  implicit none

!This subroutine gives the normalized Legendre functions to the degree lmax  at the value X. The result is a (lmax x lmax)-size matrix with lines corresponding to l and columns to m

!Dummy variables
  double precision, intent(in)::X
  double precision, dimension(0:lmax,0:lmax), intent(out)::plm
  integer         , intent(in) :: lmax

!Local variables
  double precision:: aa,as,s,sum,c,d
  integer::i,j,l,m,l0,i0,i1,i2


  plm=0

  if (X>-1 .and. X<1) then    !Domain
     plm(0,0)=1
     plm(1,0)=X

     do l=0,lmax-2    !First column of the matrix (m=0)
        plm(l+2,0)=((2*(l+1)+1)*X*plm(l+1,0)-(l+1)*plm(l,0))/(l+2)
     enddo

     do l=1,lmax    !Normalization of the first column
        Plm(l,0)=dsqrt(2.d0*l+1)*Plm(l,0)
     enddo


!Definition for x=-1 or x=1
     AA=1.d0-X**2
     if(AA.eq.0) then
        do l=0,lmax
           do m=1,l
              Plm(l,m)=0
           enddo
        enddo
        RETURN
     endif


     AS=dSQRT(AA)
     Plm(1,1)=AS*dSQRT(1.5D0)
     Plm(2,1)=3*AS*X*dSQRT(2.5D0/3.d0)
     Plm(2,2)=3*AA*dSQRT(5.D0/24.d0)


     do l=3,lmax
        SUM=1.d0
        do  J=1,l
           S=dble(2*J-1)/(2.d0*J)
           SUM=SUM*dSQRT(AA*S)
        enddo
        Plm(l,l)=SUM*dSQRT(2.D0*l+1.d0)   !plm(l,l)
        SUM=1.d0

        do  J=2,l
           S=dble(2*J-1)/(2*J-2.d0)
           SUM=SUM*dSQRT(AA*S)
        enddo
        plm(l,l-1)=X*dSQRT(2.D0*l+1.d0)*SUM   !plm(l,l-1)

        do m=l-2,1,-1
           C=2.d0*(m+1)/dsqrt((l+m+1.D0)*(l-m))
           C=C*X/AS
           D=(l+m+2)*(l-m-1)/((l+m+1.D0)*(l-m))
           D=dsqrt(D)
           Plm(l,m)=C*Plm(l,m+1)-D*Plm(l,m+2)  !the rest
        enddo
     enddo


  else
!     write(*,*) 'domaine de definition! x=', X

  end if

  do l=0,lmax
     do m=1,l
        plm(l,m) = dsqrt(2.d0)*plm(l,m)
     enddo
  enddo


end subroutine legendre_matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

