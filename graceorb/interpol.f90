subroutine interpol(x1, x2, y1, y2, x, y, n)

! Emma-Kate Potter, 6 Semptember 2010
! This subroutine does a linear interpolations between two data points (or set of data points
! in an array of data). 

! y1 and y2 are stored in an 1D array with n (n less than or equal to 4) elements

    implicit none
    integer :: i, n
    real(kind=8) ::  x1, x2, x
    real(kind=8) ::  grad 
    real(kind=8), dimension(n) ::  y1, y2, y 

    do i = 1, n
      grad = (y2(i)-y1(i))/(x2-x1)
      y(i) = y1(i)+grad*(x-x1)
    enddo 

    return
    end

