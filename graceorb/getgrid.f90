   subroutine getgrid(it, tideh)

! Emma-Kate Potter, 13 Oct, 2010
! calculate the grid of tide heights at time it
! by calling otmgrd2arr 

    use gauleg_mod
    use gm_mod
    use inmod_mod     ! PT130906: added to get the ocean tide model name into this routine

    implicit none

    real(kind=8), dimension(nlon, nlat) :: tideh
    real(kind=8), dimension(nlon*nlat) :: ocean_hts
    real(kind=8) :: time

    character*16 :: grd_file
    character*5  :: upperc

    integer :: j, ilat, ilon
    integer, dimension(5) :: it
    integer, dimension(3) :: tarray

! PT130906: use the ocean tide model name in GRACE.input to determine which grid file to open
    if(upperc(gt_oceantidemod) == "FES95")then
      grd_file = 'fes952.grid'
    else if (upperc(gt_oceantidemod) == "EOT11")then
      grd_file = 'eot11a.grid'
    else if (upperc(gt_oceantidemod) == "TPX70")then
      grd_file = 'tpxo70.grid'
    else
      write(message,'(a,a,a)')"Ocean tide model (",gt_oceantidemod,") not coded. Do it yourself please :-)"
      call status_update('FATAL','GRACEORB','getgrid',' ',message,0)
    endif

    call otmgrd2arr(it, grd_file, ocean_hts)

    j = 1 
    do ilat = 1, nlat
      do ilon = 1, nlon
        tideh(ilon, ilat) = ocean_hts(j)
        j = j + 1
      enddo
    enddo
      
    return
    end
