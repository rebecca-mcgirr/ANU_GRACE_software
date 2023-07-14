      subroutine trnsp(A,B,IROW,ICOL)
!
! Transpose a Matrix

      implicit none

      integer*4 :: icol,irow,i,j
      real*8, dimension(irow, icol) :: a
      real*8, dimension(icol, irow) :: b

      do i=1,irow
      do j=1,icol
      b(j,i)=a(i,j)
      enddo 
      enddo

      return
      end
