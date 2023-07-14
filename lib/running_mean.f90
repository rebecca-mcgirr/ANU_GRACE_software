  subroutine running_mean(nepochs,values,nsamp,values_smoothed)
  
! subroutine to calculate a running mean of a set of numbers, using a certain number of points.
!
! P. Tregoning
! 16 May 2013
!
! IN:
!     nepochs    : total number of epochs of the data to be smoothed
!     values     : column vector of data to be smoothed
!     nsamp      : number of points in the sliding window
!
! OUT:
!     values_smoothed : the smoothed values
  
  implicit none

  include '../includes/grace.param'
    
    
  integer*4 nepochs, nsamp, i, j
  double precision values(nepochs),values_smoothed(nepochs),sum,mean

  do i=(1+nsamp/2),nepochs-(1+nsamp/2)
    sum = 0.d0
    do j=i-(nsamp/2), i+(nsamp/2)
      sum = sum + values(j)
    enddo
! PT120906: now divide by the number of points as input
    values_smoothed(i) = sum / dble(nsamp+1)
  enddo

  return
  end

