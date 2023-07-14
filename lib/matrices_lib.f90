        subroutine transp(R,Rt,m,n)

!  transposes a mxn matrix
!  Written by Paul Tregoning
!  10/6/92

        implicit none
        integer i,j,m,n
        real*8 R(m,n),Rt(n,m)

!  set Rt to zero initially
        do 5 i=1,n
          do 8 j=1,m
             Rt(i,j)=0.0
8         continue
5       continue

!  perform the transpose
        do 10 i=1,m
        do 20 j=1,n
            Rt(j,i) = R(i,j)
20        continue
10      continue

        return
        end

! ======================================================================

