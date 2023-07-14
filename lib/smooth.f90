  subroutine smooth(raw,smoothed,nvals,nsmooth)

! subroutine to run a smoothing filter over a time series. 
!
! P. Tregoning
! 25 June 2013

  implicit none

  character line*60,infile*60,arg*10
  integer maxval,ioerr,nvals,i,j,nsmooth
  double precision raw(nvals),smoothed(nvals),sum,mean

! ok, so now we start smoothing. 
  do i=(1+nsmooth/2),nvals-(1+nsmooth/2)
    sum = 0.d0
    do j=i-(nsmooth/2), i+(nsmooth/2)
      sum = sum + raw(j)
    enddo
! now divide by the number of points as input
    mean = sum / dble(nsmooth+1)
    smoothed(i) = mean
  enddo

  return
  end

