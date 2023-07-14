  program msc_to_grid

! program to take a temporal gravity field estimate (in mm EWH) and convert it into values on a regular grid. The gridded values
! can then be used to convert it to a spherical harmonic model, with coefficients of EWH (which can then be converted to dimensionless 
! or to anything else)
!
! P. Tregoning
! 10 May 2017

  use mascon_mod

  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! command line argument variables
  character*150   :: mascon_file            ! input mascon configuration file
  character*100   :: msc_soln_file          ! file containing estimated mascon EWH values
  character*100   :: output_grid            ! output file of lon,lat,EWH
  real(kind=8)    :: grid_int               ! grid step size (in decimal degrees)
  real(kind=8)    :: start_lat,start_lon    ! top left corner of the grid 
  character*3     :: fit_file               ! fit/vcv for reading a fit file rather than a vcv file

! local variables
  character*150   :: arg,message,line
  integer*4       :: lumsc,lusoln,lugrid
  integer*4       :: trimlen
  integer*4       :: ilat,ilon,nlat,nlon
  integer*4       :: imsc,nmsc
  real(kind=8)    :: tmplat,tmplon
  logical         :: debug
  integer*4       :: ioerr
  character*3     :: char3

! variables related to mascons
  integer*4       :: tern_number,prim_number
  real(kind=8)    :: tmpmsc,msc_ewh(100000)   ! make it bigger than it's ever likely to be
  real(kind=8)    :: lat_msc(100000)          ! all the EWH values on a particular latitude band
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer*4       :: pnum(100000)
   integer        :: grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! unit numbers
  lumsc  = 10
  lusoln = 11
  lugrid = 12
  grid = 13
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! decode command line
  call getarg(1,mascon_file)
  if(mascon_file(1:1) == "")then
    call status_update('FATAL','UTIL','msc_to_grid','' &
       ,"Runstring: msc_to_grid mascons_feb8_c addnorm.out output_grid 0.5 89.75 0.25 fit/vcv",0)
  endif
  call getarg(2,msc_soln_file)
  open(lusoln,file=msc_soln_file,status='old')
  call getarg(3,output_grid)
  open(lugrid,file=output_grid,status='unknown')
  open(grid, file="output_pnum.grd", status='unknown')
  call getarg(4,arg)
  read(arg,*)grid_int
  call getarg(5,arg)
  read(arg,*)start_lat
  call getarg(6,arg)
  read(arg,*)start_lon
  call getarg(7,fit_file)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! open and read the input mascon file
  call read_msc_hdr(lumsc,mascon_file,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                                ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
  call allocate_mascon_arrays
  call tern_lat_bands_ell(ternary_lat_spacing/60.d0)
  call read_mascon_file(lumsc,mascon_file)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the mascon EWH estimates from the solution file
  line = " "
  do while (line(7:9) /= " MC")
    read(lusoln,'(a)')line
  enddo
  backspace(lusoln)

print*,'fit or vcv?: ',fit_file
  ioerr = 0
  nmsc = 0
  do while (ioerr == 0)
    if(fit_file == "vcv")then
      read(lusoln,'(6x,a3,41x,f15.7)',iostat=ioerr,end=1001)char3,tmpmsc
    else
      read(lusoln,'(6x,a3,54x,f15.7)',iostat=ioerr,end=1001)char3,tmpmsc
    endif
    if(ioerr == 0 .and. char3 == " MC")then
      nmsc = nmsc+1
      msc_ewh(nmsc) = tmpmsc
    endif
  enddo  
1001 write(message,'(a,i5,a)')'Found ',nmsc,' mascons in file '
  call status_update('STATUS','UTIL','msc_to_grid',msc_soln_file,message,0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! how many latitude and longitude grid points will there be?
  nlat = nint(180.d0/grid_int)+1
  nlon = nint(360.d0/grid_int)+1

! loop over lat/lon grid 
  do ilat = 1,nlat
    tmplat = start_lat - (ilat-1)*grid_int

    do ilon = 1,nlon
      tmplon = start_lon + (ilon-1)*grid_int
!print*,'ilat,ilon,tmplat,tmplon',ilat,ilon,tmplat,tmplon

! for each point, find out which ternary (and hence which primary) the grid point resides
      debug = .false.
      call calc_which_ternary(debug,tmplat,tmplon,ternary_lat_spacing/60.d0,tern_number)
      prim_number = mcon_tern_ptr(tern_number,1)
      lat_msc(ilon) = msc_ewh(prim_number)
      pnum(ilon) = prim_number
!      print*,'point ',tmplat,tmplon,' is in primary ',prim_number

!! output the coords and the EWH
!      write(lugrid,'(2f10.3,f12.5,i6)')tmplat,tmplon,msc_ewh(prim_number),prim_number
    enddo
!   just output the EWH, in latitude bands
    write(lugrid,*)lat_msc(1:nlon)
    write(grid,*)pnum(1:nlon)
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call status_update('STATUS','UTIL','msc_to_grid',' ',"End of msc_to_grid",0)
  

  end 
