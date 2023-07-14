  subroutine setup_ocean_tides_v2(ocean_tide_ht_file,tin)

! subroutine to open the binary ocean tide height file, read the header and the coords of the global grid contained within and then
! the ocean height at each grid node.
!
! (much of this code does the same as the start of util/otm_latlon )
!
! P. Tregoning
! 17 September 2013
!
! MODS
! PT140603: the new grid file contains a full list of the coords of the grid points as part of the header. We therefore don't
!           need to calculate their coordinates
!

  use tide_mod

  implicit none
  
  integer*4 ioerr, irec,i,j,npoint


  character*150,     intent(in)  :: ocean_tide_ht_file       ! name of input binary ocean tide height grid file
  double precision,  intent(in)  :: tin                      ! start time of integration (seconds since 1 Jan 2000 1200UT

! variables read from the binary file header
  real*4                         :: dlat_ogrid               ! latitude spacing of grid
  character*10                   :: tidemod                  ! tide model contained in grid file
  integer*4                      :: iyr,imonth,iday          ! epoch of start of grid file (assumed to start at 00UT)
  real*4                         :: interval                 ! time interval of grid file 
  character*4                    :: units                    ! units of time interval

! time variables (to compute the seconds of day of the start of the integration (used to do the first read of the ocean grid file)
  double precision  :: jdate, mjdate,tstart
  integer*4         :: jd

! DEBUG variables
  integer*4 :: point_num
  real*8    :: tmpxyz(3)

! first, open the binary grid file
  lutide = 44
  open(lutide,file=ocean_tide_ht_file,status='old',access='direct',iostat=ioerr,recl=4)
  if(ioerr /= 0 )then
    call status_update('FATAL','GRACEORB','setup_ocean_tides',ocean_tide_ht_file,'Error opening binary tide grid file',0)
  endif

! read the header information
! PT140603: lat/lon coords of all grid nodes now read from the header
  call read_tide_grid_head(lutide,num_msc_grd,dlat_ogrid,tidemod,iyr,imonth,iday,interval,units,num_head_rec,otide_lat,otide_lon)

! convert the temporal spacing of the grid file into seconds
  if(units == "min ")then
    dt_ocean_grid = dble(interval)*60.d0
  else if (units == "sec ")then
    dt_ocean_grid = dble(interval)
  else if (units == "deg ")then
    dt_ocean_grid = dble(interval)/3600.d0
  endif

! calculate the seconds of day (converting "tin" into julian date, then into seconds of day)
  jdate = tin/86400.d0+2451545.d0   ! actual julian *date* (non-integer)
  mjdate = jdate-2400000.5d0        ! actual modified julian *date* (non-integer)
  jd = aint(jdate)                                 ! julian day (integer)
! PT140829: fixed bug in getting seconds of calendar day when hour > 12 hrs
  tstart = (jdate-aint(jdate-0.5))*86400.d0 -43200.d0  ! seconds of day for 0000-2400UT

! now, read the appropriate epochs of the ocean binary grid file, being before and after the epoch "tstart"
  ocean_grid_epoch = int(tstart/dt_ocean_grid)
  do i=1,num_msc_grd
    irec = num_head_rec + ocean_grid_epoch*num_msc_grd + i
    read(lutide,rec=irec)dh1(i)
    read(lutide,rec=irec+num_msc_grd)dh2(i)
  enddo
  ocean_grid_epoch_old = ocean_grid_epoch



  return

  end


