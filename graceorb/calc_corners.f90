  subroutine calc_corners(rlat,rlon,nlat,nlon,corners)

! subroutine to calculate the four grid corners that encompass the input point
!
! P. Tregoning
! 10 September 2013
!
!      corners(1) = lower left
!      corners(2) = lower right
!      corners(3) = upper left
!      corners(4) = upper right

  implicit none

  integer*4, intent(out)       :: corners(4)                   ! four corners that surround the input point (LL, LR, UL, UR)
  integer*4, intent(in)        :: nlat,nlon                    ! number of lat, long grid points
  double precision, intent(in) :: rlat,rlon                    ! coordinates (in degrees) of required point
  double precision pi

  pi = 4.d0*datan(1.d0)

! special case when coordinates are 90N. (We don't need to handle South Pole because it is land)
  if(rlat == 90.d0)then
    corners(1) = int( (180.d0 - rlat*180.d0/pi)*181.d0/nlat )*nlon
    corners(2) = corners(1)
    corners(3) = 2
    corners(4) = 3
  else
    corners(1) = int( (180.d0 - rlat*180.d0/pi)*181.d0/nlat )*nlon + 360.d0/nlon * int(rlon*180.d0/pi)*nlat
    corners(2) = corners(1) + 1
    corners(3) = corners(1) + nlon
    corners(4) = corners(3) + 1
  endif

! check to see that we didn't wrap around the globe for corners(2) and corners(4)
  if ( corners(2) > mod(corners(3),nlon)*nlat ) then
    corners(2) = corners(2) - nlon
    corners(4) = corners(4) - nlon
  endif

  return
  end subroutine calc_corners









