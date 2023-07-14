  program which_mascon

! a simple test program to interface with a new subroutine that, given a lat/lon of a point, will return the number of 
! the ternary mascon in which the point resides
!
! P. Tregoning
! 2 September 2016

  use mascon_mod

  implicit none

  real(kind=8) :: lat,lon
  character*50 :: arg
  integer*4    :: tern_number
  character*150 :: mascon_file

! other variables
  real(kind=8) :: pi
  real(kind=8) :: dlat      ! latitudinal spacing (in degrees)
  real(kind=8) :: lat_spacing
  character*150 :: message

! define pi
  pi = 4.d0*atan(1.d0)


! Get the command line information
  call getarg(1,arg)
  if(arg(1:1) == " ")then
    print*,"Runstring: which_mascon lat lon [mascon_file]"
    stop
  endif
  
  read(arg,*)lat
  call getarg(2,arg)
  read(arg,*)lon
  mascon_file = " "
  call getarg(3,mascon_file)

! set the ternary latitudinal spacing
  lat_spacing = 10.d0/60.d0     ! default is 10' ternary masons in latitude

! establish the array that contains how many ternary latitude bands and how many ternary mascons per band
  call tern_lat_bands_ell(lat_spacing)

! find which ternary it is in
  call calc_which_ternary(.true.,lat,lon,lat_spacing,tern_number)

! get the primary coordinates if a mascon file was entered
  if(mascon_file(1:1) /= " ")then
! allocate the array sizes for the mascons, then read in the mascon file
  call read_msc_hdr(20,mascon_file,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                                ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
  call allocate_mascon_arrays
  call read_mascon_file(20,mascon_file) 

! PT160920: check that there is a valid ternary mascon for this location in the mascon file
  if(mcon_tern_ptr(tern_number,1) == 0)then
    write(message,'(a,2f11.4)')"Error: there is no ternary mascon in the input mascon file at location:",lat,lon
    call status_update('FATAL','UTIL','which_mascon',' ',message,0)
  endif

  write(*,'(a,2f12.6,a,i8,a,i5,a,2f12.6,a)')'Point (',lat,lon,') is in ternary mascon ' &
     ,tern_number,'. Primary mascon: ',mcon_tern_ptr(tern_number,1),',     Primary coords: (' &
     ,mcon_prim(mcon_tern_ptr(tern_number,1),1:2)*180.d0/pi,')'
  else
    print*,'Point ',lat,lon,' is in ternary number ',tern_number
  endif

  stop
  end

