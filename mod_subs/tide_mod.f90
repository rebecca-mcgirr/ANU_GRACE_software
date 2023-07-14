  module tide_mod

  implicit none

  integer(kind=4) :: tide_counter
  real(kind=8) :: tide_ti, tide_tf

! PT130917: unit number and name for tide grid file
  integer*4    :: lutide
  double precision               :: deg_per_latnode          ! number of grid nodes per degree of latitude (depends on spacing in grid file)
  integer*4                      :: num_lat_bands            ! number of latitude bands
  integer*4                      :: nmsc_per_lat(1000)       ! number of grid nodes per latitude band
  integer*4                      :: num_head_rec             ! record number of last line of header information
  integer*4                      :: num_msc_grd              ! number of grid nodes on the binary file
  double precision               :: dt_ocean_grid            ! temporal interval of binary grid file
  integer*4                      :: ocean_grid_epoch         ! epoch of ocean tide binary file that precedes the current epoch in graceorb
  integer*4                      :: ocean_grid_epoch_old     ! variable to save the value of ocean_tide_epoch to check when to update the read of the grid file
  real*4                         :: dh1(165000),dh2(165000)  ! ocean heights at epoch i and i+1, as read from the ocean binary grid file

! PT140605: lat/lon of ocean tide mascons (ie secondary mascon grid)
  real*4                         :: otide_lat(200000),otide_lon(200000) ! lat/lon of nodes on the ocean tide grid

! variables used in graceorb to define the tidal grid that is then converted to spherical harmonics
  integer*4                      :: nlat_fft,nlon_fft        ! grid for tidal heights to be converted to spherical harmonics via a FFT
  common/tide/deg_per_latnode,dh1,dh2,dt_ocean_grid,otide_lat,otide_lon  &
             ,num_msc_grd,lutide,num_lat_bands,nmsc_per_lat,num_head_rec &
             ,ocean_grid_epoch,ocean_grid_epoch_old,nlat_fft,nlon_fft


  end module tide_mod

