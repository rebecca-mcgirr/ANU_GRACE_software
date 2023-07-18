!------------------------------------------------------------------------------
! ANU/GRACE: Grace simulation Software
!------------------------------------------------------------------------------
!
! MODULE: tide_netcdf_mod
!
!> @author
!> Sebastien Allgeyer, ANU
!
! DESCRIPTION: 
!> This module deals with the tide input grid file (NetCDF)
!
! REVISION HISTORY:
! 2016-09-13 - Initial Version - SA
! 2016-09-13 - Add global variable to interface with GraceOrb - PT
! 2016-09-14 - Read the Global attribute of the netcdf file to compare with the
!              mascons identification - SA
! TODO       - Read the time vector to compare it with the time of the simulation
! 2017-03-21 - time in single precision (Don't know why I put it in double!) - SA)
!------------------------------------------------------------------------------
module tide_netcdf_mod
  use netcdf
  implicit none

  ! PT160913: define arrays to store the two epochs of ternary tide heights, and the tide epochs
  integer*4                :: tide_epoch1,tide_epoch_old
  real(kind=4),allocatable :: ternary_ht1(:),ternary_ht2(:),ternary_ht_current(:)
  real(kind=8)             :: dt_ocean_grid

  public T_open, T_close, T_read
  public T_get_info, T_get_nmascons
  public T_get_TS, T_get_msc_time, T_locate_mascons
  private check
  type netcdf_data
     private
     character(len=256)  ::  fname            !Name of the net cdf file
     integer             ::  ncid             !File ID for use in the API
     integer             ::  varid            !Netcdf Variable ID
     integer             :: lonvarid
     integer             :: latvarid
     integer             ::  n_mascons        !Netcdf Mascons numbers   
     integer             ::  n_time           !Number of time in the netcdf file    
     real, dimension(:), allocatable :: time_vect 
  end type netcdf_data

  type(netcdf_data)        :: tidedata              ! a structure to hold the netcdf data ?

  interface T_open
     module procedure open_tide_file
  end interface T_open

  interface T_close
     module procedure close_tide_file
  end interface T_close

  interface T_read
     module procedure read_tide_at_time
  end interface T_read

  interface T_get_TS
     module procedure read_fulltide_at_mascons
  end interface T_get_TS

  interface T_get_msc_time
     module procedure read_tide_msctime
  end interface T_get_msc_time

  interface T_get_nmascons
     module procedure get_nmascons
  end interface T_get_nmascons

  interface T_get_info
     module procedure get_info
  end interface T_get_info

  interface T_locate_mascons
     module procedure locate_mascons
  end interface T_locate_mascons

contains
  function open_tide_file(fname) result(self)
    character (len =*), intent(in):: fname
    type(netcdf_data):: self
    integer :: mascons_id
    integer :: time_id
    integer :: time_varid
    self%fname = fname
    call check(nf90_open(fname,NF90_NOWRITE, self%ncid))
    call check( nf90_inq_varid(self%ncid, "outtide", self%varid) )
    call check( nf90_inq_dimid(self%ncid, "msc3", mascons_id)) 
    call check( nf90_inq_dimid(self%ncid, "ntime", time_id)) 
    call check( nf90_inquire_dimension(self%ncid, mascons_id, len = self%n_mascons))
    call check( nf90_inquire_dimension(self%ncid, time_id, len = self%n_time))
    allocate(self%time_vect(self%n_time))
    call check( nf90_inq_varid(self%ncid, "time", time_varid) )
    call check( nf90_get_var(self%ncid, time_varid, self%time_vect))
    call check( nf90_inq_varid(self%ncid, "msc3_lon", self%lonvarid) )
    call check( nf90_inq_varid(self%ncid, "msc3_lat", self%latvarid) )

  end function  open_tide_file

  function get_nmascons(self) result(n_mascons)
    type(netcdf_data):: self
    integer :: n_mascons
    n_mascons = self%n_mascons
  end function get_nmascons

  function get_info(self, arg) result(string)
    type(netcdf_data), intent(in):: self
    character (len = 80) :: string
    character (len = *), intent(in) :: arg
    call check( nf90_get_att(self%ncid, nf90_global, arg, string))
  end function get_info


  subroutine close_tide_file(self)
    type(netcdf_data), intent(inout) :: self
    call check( nf90_close(self%ncid) )
  end subroutine close_tide_file


  subroutine  read_tide_at_time(self, idx_t, table)
    type(netcdf_data), intent(inout)         ::  self
    integer, intent(in):: idx_t
    real, intent(out):: table(:)
    integer ::  count_a(2)
    integer ::  start_a(2)
    start_a(1) = idx_t
    start_a(2) = 1
    count_a(1) = 1
    count_a(2) = self%n_mascons
    call check(nf90_get_var(self%ncid, self%varid, table, start=start_a, count=count_a))
  end subroutine read_tide_at_time

  subroutine read_fulltide_at_mascons(self, idx_msc, table)
    type(netcdf_data), intent(inout)         ::  self
    integer, intent(in):: idx_msc
    real, intent(out):: table(:)
    integer :: count_a(2)
    integer :: start_a(2)
    start_a(1) = 1
    start_a(2) = idx_msc
    count_a(1) = self%n_time
    count_a(2) = 1
    call check(nf90_get_var(self%ncid, self%varid, table, start=start_a, count=count_a))
  end subroutine read_fulltide_at_mascons

  subroutine locate_mascons(self, idx_msc, lon, lat )
    type(netcdf_data), intent(in) :: self
    integer, intent(in) :: idx_msc
    real, intent(out) :: lon, lat

    call check(nf90_get_var(self%ncid, self%lonvarid, lon, start = (/ idx_msc /) ) )
    call check( nf90_get_var(self%ncid, self%latvarid, lat, start = (/ idx_msc /) ) )

  end subroutine locate_mascons

  function read_tide_msctime(self, idx_msc, mjd, timevalue) result(tide)
    type(netcdf_data), intent(in) :: self
    integer, intent(in) :: idx_msc
    real, intent(in) :: timevalue
    real, intent(in) :: mjd
    integer  :: i
    integer :: count_a(2)
    integer :: start_a(2)
    real :: tidevalue(2)
    real :: tide
    real :: t

    t = timevalue + mjd
    !print*,  t

    if (t < self%time_vect(1) .or. t > self%time_vect(self%n_time)) then
       print *, "Request time : ", t, "Which is not included in the tide file", &
            self%time_vect(1), self%time_vect(self%n_time)
       stop "Request time error"
    endif


    do i=1,self%n_time
       if (t <= self%time_vect(i)) then
          exit
       endif
    enddo

    if (t == self%time_vect(i)) then
       start_a(1) = i
       start_a(2) = idx_msc
       count_a(1) = 1
       count_a(2) = 1
       call check(nf90_get_var(self%ncid, self%varid, tidevalue,  &
            start=start_a, count=count_a))
       !print *, t == self%time_vect(i) , tidevalue

       tide = tidevalue(1)
    else
       start_a(1) = i-1
       start_a(2) = idx_msc
       count_a(1) = 2
       count_a(2) = 1
       call check(nf90_get_var(self%ncid, self%varid, tidevalue, start=start_a,  &
            count=count_a))
       !! Linear interpolation
       !print *, t == self%time_vect(i) , tidevalue
       !print *, self%time_vect(i-1) , self%time_vect(i)
       tide = tidevalue(1) + (t - self%time_vect(i-1)) /  &
            (self%time_vect(i) - self%time_vect(i-1)) * (tidevalue(2) - tidevalue(1))
    end if
  end function read_tide_msctime


  subroutine check(status)
    integer, intent ( in) :: status
    if(status /= nf90_noerr) then 
       print *, trim(nf90_strerror(status))
       stop "Error in Netcdf"
    end if
  end subroutine check

end module tide_netcdf_mod
