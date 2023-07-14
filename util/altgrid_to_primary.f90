  program altgrid_to_primary
  
! program to convert the lat/lon of the altimetry gridded data into the primary mascon geometry of our mascons. That is,
! identify in which primary mascon each altimetry grid point resides. This makes a lookup table for use in python.
!
! P. Tregoning
! 29 March 2022

  use mascon_mod      ! defines all the variables needed for the mascons
  
  
  implicit none
  
! command line arguments
  character*150  :: grid_file, mascon_file, output_file                ! self explanatory
  
! file unit numbers
  integer*4      :: lualt=10,lumsc=11,luout=12

! mascon variables
  real(kind=8)   :: lat_spacing,tmplat,tmplon,tmp_sec,tmp_sec_err
  integer*4      :: tern_number,prim_number
  integer*4      :: tmp_ilat,tmp_ilon
! local variables
  character*100  :: arg,message
  integer*4      :: ioerr,npoints
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! decode the runstring
  grid_file = " "
  call getarg(1,grid_file)
  if(grid_file(1:1) == " ")then
    message = "Runstring: altgrid_to_primary altgrid.dat mascons_stage5_V006_200km altgrid.lookup"
    call report_stat('FATAL','UTIL','altgrid_to_primary',' ',message,0)
  else
    open(lualt,file=trim(grid_file))
  endif
  call getarg(2,mascon_file)
  call getarg(3,output_file)
  open(luout,file=output_file,status='unknown')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set the ternary latitudinal spacing
  lat_spacing = 10.d0/60.d0     ! default is 10' ternary masons in latitude

! establish the array that contains how many ternary latitude bands and how many ternary mascons per band
  call tern_lat_bands_ell(lat_spacing)

! get the primary coordinates if a mascon file was entered
  if(mascon_file(1:1) /= " ")then
! allocate the array sizes for the mascons, then read in the mascon file
    call read_msc_hdr(lumsc,mascon_file,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                                ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
    call allocate_mascon_arrays
    call read_mascon_file(lumsc,mascon_file) 
  else
    call report_stat('FATAL','UTIL','altgrid_to_primary',mascon_file,"no mascon file name entered",0)
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now, read through until the end of the file is reached
  ! skip the first and second line in the file
  read(lualt,'(a)')arg
  !write(*,'(a,a)')"first  line:",arg
  read(lualt,'(a)')arg
  !write(*,'(a,a)')"second line:",arg
  
  npoints = 0
  ioerr = 0
  do while (ioerr == 0)
    read(lualt,*,end=1000,iostat=ioerr)tmp_ilat,tmp_ilon,tmplat,tmplon,tmp_sec,tmp_sec_err
    if(tmplon < 0.d0)tmplon = tmplon + 360.d0
    
    if(ioerr == 0)then
      npoints = npoints + 1
      ! find which ternary it is in
      call calc_which_ternary(.false.,tmplat,tmplon,lat_spacing,tern_number) 
! DEBUG:
!print*,tern_number,mcon_tern_ptr(tern_number-10:tern_number+10,1)

      ! find which primary mascon this ternary mascon belongs to
      prim_number = mcon_tern_ptr(tern_number,1)
      !write out the input information plus the primary mascon in which it resides
      !write(*,'(2i10,2f10.4,i8,a)')tmp_ilat,tmp_ilon,tmplat,tmplon,prim_number," altimetry grid vs primary mascon number"
      write(luout,'(2i10,4f10.4,i8)')tmp_ilat,tmp_ilon,tmplat,tmplon,tmp_sec,tmp_sec_err,prim_number
    endif
  enddo
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
1000  continue
      write(message,'(a,i6,a)')"End of file reached, but first found ",npoints," grid points in file"
      call report_stat('STATUS','UTIL','altgrid_to_primary',grid_file,message,0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
     
