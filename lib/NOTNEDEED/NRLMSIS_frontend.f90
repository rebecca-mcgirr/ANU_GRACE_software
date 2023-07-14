
  program NRLMSIS_frontend

! program to interface with the NRLMSIS subroutines to generate an estimate of the atmospheric density
! given location, day of year and several other controlling factors.
!
! P. Tregoning
! 23 October 2015

  use utils_constants
  use physics_msis

  implicit none

  character arg*50

! location varaibles
  real(kind=8) :: lat, lon, alt

! epoch variables
  integer*4    :: doy                ! day of year
  real(kind=8) :: sod                ! seconds of day
  real(kind=8) :: stl                ! local solar time
 
! solar flux values
  real(kind=8) :: F107A, F107

! magnetic index
  real(kind=8) :: ap(7)
  data ap/7*100./ 

! pressure
  real(kind=8) :: pressure

! don't know really what this one is, but it has a value of 48
  integer*4    :: mass

! output variables
  real(kind=8) :: d(9)
  real(kind=8) :: t(2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! Decode Command Line inputs !!!!!!!!!!!!!!!!!!!
  call getarg(1,arg)
  if (arg(1:1) == "")then
    print*,'NRLMSIS_frontend lat lon alt doy sec_of_day f107a f107 '
    print*,'e.g. NRLMSIS_frontend -35. 140.5 450.0 265 14400 180. 150. '
    stop
  endif
! coordinates
  read(arg,*)lat
  call getarg(2,arg)
  read(arg,*)lon
  call getarg(3,arg)
  read(arg,*)alt
! epoch
  call getarg(4,arg)
  read(arg,*)doy
  call getarg(5,arg)
  read(arg,*)sod
! solar flux
  call getarg(6,arg)
  read(arg,*)f107a
  call getarg(7,arg)
  read(arg,*)f107




! calculate local_solar_time
  stl = sod/3600. + lon/15.0


! call the NRMLSIS subroutine gtd7d to compute the total mass density at the input coords and altitude
  mass = 48 ! not sure why, but this seems to get a value of 48 in the example code
  call gtd7d(doy,sod,alt,lat,lon,stl,f107a,f107,ap,mass,d,t)  ! the total mass density was returned from the subroutine gtd7d in value d(6)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! Output the total mass density !!!!!!!!!!!!!!!!!
  write(*,'(3f12.3,i6,f10.1,2f6.1,e16.8)')lat,lon,alt,doy,sod,f107a,f107,d(6)
  end


