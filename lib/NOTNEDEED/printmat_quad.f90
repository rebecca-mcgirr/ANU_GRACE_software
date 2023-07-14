        subroutine printmat(mat,m,n,comment)

!  prints out a m x n matrix
!  P Tregoning 9th June 1992

        implicit none
        integer m,n,i,j,indx
        real*8 mat(m,n)
        character*80 comment

!  search for the first blank in the comment - only print up to it
        indx = index(comment,' ')
        write(*,'(a)')comment(1:indx)
        do 10 i=1,m
            write(*,*)(mat(i,j),j=1,n)
10      continue

        print*,' '
!        write(*,'(a)')' press return to continue : '
!        read(*,'(a)')ent

        return
        end

!!!!!!!*************!!!!!!!!!!!!!!****************
        subroutine printmat_int(mat,m,n,comment)

!  prints out a m x n matrix
!  P Tregoning 9th June 1992

        implicit none
        integer m,n,i,j,indx
        integer*4 mat(m,n)
        character*80 comment

!  search for the first blank in the comment - only print up to it
        indx = index(comment,' ')
        write(*,'(a)')comment(1:indx)
        do 10 i=1,m
            write(*,*)(mat(i,j),j=1,n)
10      continue

        print*,' '
!        write(*,'(a)')' press return to continue : '
!        read(*,'(a)')ent

        return
        end

