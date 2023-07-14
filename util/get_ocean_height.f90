  subroutine get_ocean_height(npoints,first,lat,lon,timearr,tidemod,ocean_ht)

! subroutine to interact with the two routines that interact with spotl .... phew !
!
! P. Tregoning
! 11 September 2013

  implicit none

  integer*4,        intent(in)  :: npoints                    ! number of points in the lat/lon vectors
  double precision, intent(in)  :: lat(npoints),lon(npoints)  ! coordinates of point(s) at which to estimate tidal height
  integer*4,        intent(in)  :: timearr(5)                 ! yr, doy, hr, min,sec of epoch to estimate tide
  character,        intent(in)  :: tidemod*10                 ! requested tide model
  logical,          intent(inout)  :: first                      ! flag for first time into subroutine

  double precision, intent(out) :: ocean_ht(npoints)          ! computed ocean heights

! local variables and things requred for spotl interaction
  integer*4, allocatable        :: cart_codes(:,:,:)          ! cartwright codes as read from the input spotl codes
  character, allocatable        :: darw_codes(:,:)*4          ! darwin codes
  real*4,    allocatable        :: amp(:,:)                   ! constituent amplitudes
  real*4,    allocatable        :: phase(:,:)                 ! constituent phases
  logical                       :: debug_flag                 ! whether to print to screen the phases and amplitudes
  integer*4                     :: ntmcon                     ! number of tidal constituents
  integer*4                     :: max_con                    ! maximum number of constituents
  integer*4 i,j
  character message*200

  save cart_codes,darw_codes,amp,phase,debug_flag

! now get the amplitudes and phases for this point
  if(first) then
! allocate the arrays
    max_con = 20
! set the debug flag to true so that we see the tidal consituents to be used
    debug_flag = .true.
    allocate(cart_codes(6,max_con,npoints))
    allocate(darw_codes(max_con,npoints))
    allocate(amp(max_con,npoints))
    allocate(phase(max_con,npoints))
    call status_update('STATUS','UTIL','get_ocean_height',' ','Computing amplitudes/phases',0)
    call spotl_point_setup_v2(lat,lon,npoints,tidemod,cart_codes,darw_codes,ntmcon,max_con,amp,phase,debug_flag) 
    debug_flag = .false.
  endif

! then compute the tidal height
  write(message,'(a,5i5,i8)')'Computing ocean heights for epoch',timearr,npoints
!  call status_update('STATUS','UTIL','get_ocean_height',' ',message,0)
!  do i=1,npoints
  
  do i=1,npoints
!   do i=1,50
!  if(first .and. i < 50)print*,'point,lat,lon',i,lat(i),lon(i),timearr
    if(amp(1,i) > 0.0 .and. amp(2,i) > 0.0)then     ! don't compute if the point is over land

!  do j=1,12
!    print*,'calling spotl i,amp,phase',i,amp(j,i),phase(j,i)
!  enddo
      call spotl_point_calc(timearr,ntmcon,cart_codes(:,:,i),darw_codes(:,i),amp(:,i),phase(:,i),ocean_ht(i))
!  print*,'ocean height computed in get_ocean_height: ',i,ocean_ht(i),timearr
    endif
  enddo

!  deallocate(cart_codes)
!  deallocate(darw_codes)
!  deallocate(amp)
!  deallocate(phase)

! set first to false
  if(first) first = .false.

! and we're done
  return
  end

