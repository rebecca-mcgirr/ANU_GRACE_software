  subroutine sort_longitude(nlon,longitudes)

! subroutine to reorganise a set of longitude values into increasing order, then return the pointer array to reflect the
! new order of the values
!
! P. Tregoning
! 17 September 2013

  implicit none

  integer*4,          intent(in)    :: nlon              ! number of longitude values
  double precision,   intent(inout) :: longitudes(nlon)  ! longitude values to be re-ordered

  integer*4 tmp_ptr,i,j
  double precision tmp_lon

! run a bubble sort to put them in ascending order 
  do i = nlon , 1, -1 
    do j=1,i
      if(longitudes(j) > longitudes(i) ) then
! swap the longitudes
        tmp_lon = longitudes(i)
        longitudes(i) = longitudes(j)
        longitudes(j) = tmp_lon
      endif
    enddo
  enddo

  return
  end

