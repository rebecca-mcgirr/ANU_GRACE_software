        subroutine matmult(M1,M2,M3,l,m,n)

!c  multiplies lxm matrix M1 by a mxn matrix M2 to give a lxn matrix M3
!c  written by Paul Tregoning 9th June 1992
!c  EKP130110 - updated for fortran 90 standard do loop construct
 
    implicit none
    integer :: i,j,k,l,m,n
    real*8 :: M1(l,m),M2(m,n),M3(l,n),temp
    temp = 0.0d0

    do k=1,n
      do i=1,l
        do j=1,m
          temp = M1(i,j) * M2(j,k) + temp
        enddo
        M3(i,k) = temp

        temp = 0.0d0
      enddo
    enddo
    return
    end

