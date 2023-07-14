! ---------------------------------------------
      subroutine invert(a,b,n)

!  subroutine stolen from someone else

      implicit none

      integer n,i,j,k
      real*8 a(n,n),b(n,n),z

      do  i=1,n
        do  j=1,n
          b(i,j) = 0.0
        enddo
        b(i,i) = 1.0
      enddo

      do k = 1,n
        do i = 1,n
	       if(i.ne.k)then
	         z = a(i,k)/a(k,k)
            do j=1,n
              a(i,j) = a(i,j)-a(k,j)*z
              b(i,j) = b(I,j) - b(k,j)*z
            enddo
	       endif
	     enddo
	     z = a(k,k)
	     do j=1,n
	       a(k,j) = a(k,j)/z
	       b(k,j) = b(k,j)/z
	     enddo
	   enddo

      return
      end


