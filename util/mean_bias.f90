  subroutine mean_bias(nepochs,accobs,bias)

! subroutine to determine the best possible mean bias estimate so that we can break the correlation between bias and scale in the least squares estimate.
! I will use the shadow epochs - unaffected by thrusts - to estimate a bias that will calibrate the observations properly.
!
! P. Tregoning
! 27 November 2015

  implicit none

  integer*4   , intent(in)     :: nepochs
  real(kind=8), intent(inout)  :: accobs(nepochs,5)
  real(kind=8), intent(inout)  :: bias(3)

! local variables
  integer*4     :: i,j,iepoch,nvals(3)
  real(kind=8)  :: sumvals(3)


  nvals = 0

! loop through all the epochs and sum the non-thrust observations during eclipse
  do iepoch = 1,nepochs
    if(accobs(iepoch,5) == 3)then   ! it is an eclipse observation

    ! loop over the three accelerometer axes
      do j=1,3
        nvals(j) = nvals(j) + 1
        sumvals(j) = sumvals(j) + accobs(iepoch,1+j)
      enddo
    endif
  enddo

! if there were no observations (i.e. never in eclipse) then just sum all the non-thrust observations
  if(nvals(2) == 0)then
    call status_update('STATUS','UTIL','mean_bias',' ','No eclipse obs - use mean of all obs to calculate bias',0)
    do iepoch = 1,nepochs
      if(accobs(iepoch,5) /= 1)then   ! it is not affected by thrusts

    ! loop over the three accelerometer axes
        do j=1,3
          nvals(j) = nvals(j) + 1
          sumvals(j) = sumvals(j) + accobs(iepoch,1+j)
        enddo
      endif
    enddo
  endif
   

! now, calculate the mean for each axis
  bias(:) = -1.0 *sumvals(:)/dble(nvals(:))*1.e9

! that's it
  return
  end
