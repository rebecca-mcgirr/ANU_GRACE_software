  subroutine which_primary_msc(luindx,first,lat,lon,msc_ternary,msc_primary)

! subroutine find in which primary and ternary mascon a lat/lon position resides. Code adapted from part of find_mascon.f90
!
! P. Tregoning
! 22 October 2014

  implicit none


! input variables
  integer*4    ::  luindx          ! index number of mascon_index   file
  logical      ::  first           ! logical to indicate whether the mascon index file should be read
  real(kind=8) ::  lat,lon         ! coordinates for which we seek the mascons

! output variables
  integer*4    ::  msc_ternary     ! ternary mascon number
  integer*4    ::  msc_primary     ! primary mascon number

! parameter declarations
  integer*4, parameter:: max_prime = 5000     ! Max number of primary cells
  integer*4, parameter:: max_tern  = 1500000  ! Max number of ternary cells
  integer*4, parameter:: ppd_tern = 6         ! Max number of ternary cells per degree
  integer*4, parameter:: max_tern_lat = 1 + (180 * ppd_tern) ! Max number of ternary latitude steps

! local variables
  integer*4 :: num_tern_lat        ! number of ternary mascon latitude bands
  integer*4 :: num_tern_long       ! number of ternary mascons in a particular latitude band
  integer*4 :: pts_per_deg         ! number of ternary mascons per degree of latitude
  integer*4 :: num_lat_cell        ! number of the latitude cell (?)
  integer*4 :: num_lon_cell        ! number of the longitude cell (?)
  integer*4 :: cell_lon            ! longitude cell number in the required latitude band
  integer*4 :: mcon_count          ! a counter
  integer*4 :: mcon_index(max_tern)! index of ternary mascons sequentially wrapping eastward and soutward from north pole
  real(kind=8) :: del_lat          ! size of latitude step
  real(kind=8) :: del_lon          ! size of lontigude step for a particular latitude band

  character :: message*250
  integer*4 :: i,j

! allocatable variables
  integer*4, allocatable :: lons_per_lat(:)     ! number of longitudes per each band of latitude

  save mcon_index, lons_per_lat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the entire mascon_index file only the first time
  if(first)then
! determine which mascon cell (at highest resolution), the satellite overflies
    read(luindx,*) num_tern_lat, pts_per_deg
    if ( num_tern_lat .gt. max_tern_lat ) then
      write(message,'(a,2i6)')"Too many ternary mascon latitudes ", num_tern_lat, max_tern_lat
      call status_update('FATAL','LIB','which_primary_msc',' ',message,0)
    endif
    if ( pts_per_deg .gt. 0.0d0 ) then
      del_lat = 1.0d0/dble(pts_per_deg)
    else
      del_lat = -dble(pts_per_deg)
    endif

! allocate the lons_per_lat variable
  allocate(lons_per_lat(num_tern_lat))

! set the latitude band to be the last one, so that the whole file gets read
   num_lat_cell = num_tern_lat

!   read through to this latitude band
    mcon_count = 0
    do i = 1, num_lat_cell
      read(luindx,*) num_tern_long
!     store it away for later use
      lons_per_lat(i) = num_tern_long
      do j = 1, num_tern_long
        mcon_count = mcon_count + 1
        if ( mcon_count .gt. max_tern ) then
          write(message,'(a,2i6)')"Too many ternary mascon cells ", mcon_count, max_tern            
          call status_update('FATAL','GRACEORB','graceorb',' ',message,0)
        endif
        read(luindx,*) mcon_index(mcon_count)
      enddo
    enddo 

    call status_update('STATUS','LIB','which_primary_msc',' ',"Have read the entire mascon index file",0)
    first = .false.
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now for each particular requested lat/lon

! which latitude band? 
  num_lat_cell = 1 + int((90.d0-lat)/del_lat + 0.5)

! read through to this latitude band
  mcon_count = 0
  do i = 1, num_lat_cell
    mcon_count = mcon_count + lons_per_lat(i)
    num_tern_long = lons_per_lat(i)
  enddo 

! determine which mascon cell (at highest resolution), the specified cell lies in
! which longitude?
  del_lon = 360.d0/dble(num_tern_long)
  num_lon_cell = 1 + int(lon/del_lon + 0.5)
  if ( num_lon_cell .gt. num_tern_long ) num_lon_cell = 1

! calculate the ternary mascon number
  msc_ternary = mcon_count - num_tern_long + num_lon_cell
! now extract the index for the primary mascon for this ternary mascon
  msc_primary = mcon_index(msc_ternary)


! rewind the file so that it might be read again
!  rewind(luindx)

  return
  end




