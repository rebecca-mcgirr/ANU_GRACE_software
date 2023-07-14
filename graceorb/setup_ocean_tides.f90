  subroutine setup_ocean_tides(ocean_tide_ht_file,tin)

! subroutine to open the binary ocean tide height file, read the header, calculate the coords of the global grid contained within and then
! derive the amplitudes and phases of each consituent at each grid node.
!
! (much of this code does the same as the start of util/otm_latlon )
!
! P. Tregoning
! 17 September 2013
!
!
! to determine the nearest grid node wrt any given lat/lon location we just need to know 
! a) how many latitude bands there are
! b) how many nodes on each latitude band
! c) the latitude band spacing
! d) the number of header records in the binary file
! e) the total number of grid nodes on the file (num_msc_grd)
!
! therefore, the following variables are all passed out of this routine via tide_mod.f90
!   deg_per_latnode, num_lat_bands, nmsc_lat_bands, num_head_rec, num_msc_grid , dt_ocean_grid

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

! set of latitude and longitudes of the grid nodes within the binary file (derived here)
  double precision               :: lat(165000),lon(165000)  ! (derived) grid coordinates of the nodes in the binary file

! time variables (to compute the seconds of day of the start of the integration (used to do the first read of the ocean grid file)
  double precision  :: jdate, mjdate,tstart
  integer*4         :: jd


! first, open the binary grid file
  lutide = 44
  open(lutide,file=ocean_tide_ht_file,status='old',access='direct',iostat=ioerr,recl=4)
  if(ioerr /= 0 )then
    call status_update('FATAL','GRACEORB','setup_ocean_tides',ocean_tide_ht_file,'Error opening binary tide grid file',0)
  endif

! read the header information
  call read_tide_grid_head(lutide,num_msc_grd,dlat_ogrid,tidemod,iyr,imonth,iday,interval,units,num_head_rec)

! derive the number of degrees of latitude for each dlat element
  deg_per_latnode = dlat_ogrid/60.d0

! convert the temporal spacing of the grid file into seconds
  if(units == "min ")then
    dt_ocean_grid = dble(interval)*60.d0
  else if (units == "sec ")then
    dt_ocean_grid = dble(interval)
  else if (units == "deg ")then
    dt_ocean_grid = dble(interval)/3600.d0
  endif

! now we need to construct a vector of latitudes and longitudes that define all the grid points. The next block of information in the
! binary grid file tells us how many latitude bands there are
  irec = num_head_rec + 1
  read(lutide,rec=irec)num_lat_bands

! read in how many secondary mascons lie in each band
  npoint = 0
  do i=1,num_lat_bands
    irec = irec+1
    read(lutide,rec=irec)nmsc_per_lat(i)
    do j=1,nmsc_per_lat(i)
! derive and store the coords
      npoint = npoint + 1
      lon(npoint) =  (j-1)*360.d0/dble(nmsc_per_lat(i))
      lat(npoint) = 90.d0-dble(i-1)*deg_per_latnode
    enddo
  enddo

  num_head_rec = irec

! calculate the seconds of day (converting "tin" into julian date, then into seconds of day)
  jdate = tin/86400.d0+2451545.d0   ! actual julian *date* (non-integer)
  mjdate = jdate-2400000.5d0        ! actual modified julian *date* (non-integer)
  jd = aint(jdate)                                 ! julian day (integer)
  tstart = (jdate-aint(jdate))*86400.d0 -43200.d0  ! seconds of day for jdate

! now, read the appropriate epochs of the ocean binary grid file, being before and after the epoch "tstart"
  ocean_grid_epoch = int(tstart/dt_ocean_grid)
  do i=1,num_msc_grd
    irec = num_head_rec + ocean_grid_epoch*num_msc_grd + i
    read(lutide,rec=irec)dh1(i)
    read(lutide,rec=irec+num_msc_grd)dh2(i)
! DEBUG:
!  print*,'setup_ocean_tide i,lat(i),lon(i),dh1(i),dh2(i)',i,lat(i),lon(i),dh1(i),dh2(i)
  enddo
  ocean_grid_epoch_old = ocean_grid_epoch



  return

  end


