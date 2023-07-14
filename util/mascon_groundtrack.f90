  program mascon_groundtrack

! program to get the mascon groundtrack (from GNV1B file) and mascon EWH (from vcv file) and 
! output the information for whatever use (plotting, histograms etc)
!
! P. Tregoning
! 26 June 2017

  use mascon_mod

  implicit none

! command line arguments
  character  :: gnv1b_file*100         ! GNV1B file from which we get GRACE A coords
  character  :: mascon_file*150        ! mascon prim/sec/tern file
  character  :: mcon_ewh_file*100      ! file containing mascon EWH estimates

! coords
  integer*4    :: gracesec
  real(kind=8) :: xyz(3),sph(3),lat,lon,radius       ! GRACE A coordinates

! mascons
  real(kind=8),allocatable :: mcon_EWH_sigma(:)
  integer*4    :: prim_number                    ! pointer to the primary mascon
  integer*4    :: tern_number                    ! pointer to the ternary mascon
  integer*4    :: n_msc                             ! number of mascons in the vcv file

! unit numbers
  integer*4,parameter    :: lu_gnv=10,lu_msc=11,lu_ewh=12

! other stuff
  integer*4     :: ioerr,i,j
  character     :: sat*1,blah*1
  real(kind=8)  :: pi

  character*100 :: line
  character*200 :: message

  pi = 4.d0*datan(1.d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! decode command line
! GNV1B file
  call getarg(1,gnv1b_file)
  if(gnv1b_file(1:1) == "")then
    write(message,'(a)')"Runstring: mascon_groundtrack GNV1B_file mascon_stage4 mascon_vcv "
    call status_update('FATAL','UTIL','mascon_groundtrack',' ',message,0)
  else
    open(lu_gnv,file=gnv1b_file,status='old',iostat=ioerr)
    if(ioerr /= 0)then
      write(message,'(a)')"Error opening file. Does it exist? "
      call status_update('FATAL','UTIL','mascon_groundtrack',gnv1b_file,message,0)   
    endif
  endif
! mascon coord file
  call getarg(2,mascon_file)

! mascon soluton file
  call getarg(3,mcon_ewh_file)
  open(lu_ewh,file=mcon_ewh_file,status='old',iostat=ioerr)
  if(ioerr /= 0)then
    write(message,'(a)')"Error opening vcv file. Does it exist? "
    call status_update('FATAL','UTIL','mascon_groundtrack',mcon_ewh_file,message,0)   
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read in the mascon primary/secondary/ternary coords
  call read_msc_hdr(lu_msc,mascon_file,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                                ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)

! establish the information to identify each ternary mascon from lat/lon 
  call tern_lat_bands_ell(ternary_lat_spacing/60.d0)

! allocate the array sizes for the mascons
  call allocate_mascon_arrays

! read in the mascon information
  call read_mascon_file(lu_msc,mascon_file)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read in the mascon EWH estimates
! first, get the number of mascons. This is on the 5th line
! PT170626: don't do this ... just use max_prim and assume all is well.
!  do i=1,4
!    read(lu_ewh,'(a)')line
!  enddo
!  read(lu_ewh,'(40x,i7)')n_msc
  n_msc = max_prim

! check that there are the right number of mascons
  if (n_msc /= max_prim)then
    write(message,'(a,i7,a,i7,a)')"Error: number of mascons in vcv file (",n_msc &
            ,") does not match number in mascon file (",max_prim,")"
    call status_update('FATAL','UTIL','mascon_groundtrack',' ',message,0)
  endif

  allocate(mcon_EWH_sigma(n_msc))

! now skip down to the mascons
  line=" "
  do while (line(7:9) /= " MC")
    read(lu_ewh,'(a)')line
  enddo
  backspace(lu_ewh)

! now read them in
  do i=1,n_msc
    read(lu_ewh,'(47x,2f17.7)')mcon_EWH(i),mcon_EWH_sigma(i)
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! skip the header of the GNV1B file
  line=" "
  do while (line(1:13) /= "END OF HEADER")
    read(lu_gnv,'(a)')line
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!  MAIN LOOP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! loop through the entire GNV1B file, determine in which mascon the satellite resides, match 
! up the coords with the mascon EWH value and output

  ioerr = 0
  do while (ioerr == 0)
! get xyz coords
    read(lu_gnv,*,end=1000)gracesec,sat,blah,xyz(1:3)
! convert to lat/lon
    call cart_to_sph(xyz,sph) 
    lat = sph(1)*180.d0/pi
    lon = sph(2)*180.d0/pi

!print*,xyz,sph,lat,lon
! find in which ternary the satellite resides
    call calc_which_ternary(.false.,lat,lon,10.d0/60.d0,tern_number)

! which primary mascon is this
    prim_number = mcon_tern_ptr(tern_number,1)

! output
    write(*,'(4f15.4,i15,2i8,a)')lat,lon,mcon_EWH(prim_number),mcon_EWH_sigma(prim_number),gracesec,prim_number,tern_number &
                                 ," mascon EWH"
  enddo

1000 call status_update('STATUS','UTIL','mascon_groundtrack',gnv1b_file,'Reached end of GNV1B file',0)
  end



