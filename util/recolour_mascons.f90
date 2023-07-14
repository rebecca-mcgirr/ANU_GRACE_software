  program recolour_mascons

! program to change the colour code of particular specified mascons. This is a diagnostic program and is 
! not meant to be very smart .....
!
!  P. Tregoning  
!  20 March 2019

  use mascon_mod

  implicit none

! input/output files
  character :: msc_input*150       ! name of input  primary mascon file containing primary mascon coords
  character :: msc_output*150      ! name of output primary mascon file containing the recomputed CoM of the adjusted primary mascons
  character :: recolour_file*100   ! file listing primary mascon numbers and the colour to which they should be changed

! file unit numbers
  integer*4 :: lumsc_in, lumsc_out, lucolour

! variables to write out the mascon information
  real(kind=8) :: msc_crds(3)
  character*6  :: msc_code
  integer*4    :: pointers(2),colours(2)
  character*15 :: region

! variables for reading the recolour file
  integer*4    :: tmp_prim,tmp_colour,tmp_tern

! variables for writing date and time into header of output file
  character(8)  :: date
  character(10) :: time
  character(5)  :: timezone
  character*12  :: hashcode          ! temporary variable used to calculate the new character code for the output mascon file
 
! counters
  integer*4 :: i,j,k,imsc,itern
  integer*4 :: imsc_prim       ! number of input primary mascon when reading the input ternary mascon file
  integer*4 :: ntern           ! number of input ternary mascons per input primary mascon

! other stuff
  integer*4 :: ioerr
  character :: message*250
  character :: arg*100
!  real(kind=8) :: tmprad,tmparea,tmpdensity,tmpdepth,tmpdummy
!  integer*4    :: trimlen
  real(kind=8) :: pi





  pi = 4.d0*datan(1.d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      decode command line
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call getarg(1,msc_input)
  if (msc_input(1:1) == " ")then
    write(message,*)'Runstring: recolour_mascons msc_infile msc_outfile prim_colour_file'
    call status_update('FATAL','UTIL','recolour_mascons',' ',message,0)
  endif
  call getarg(2,msc_output)
  call getarg(3,recolour_file)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      open files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  lumsc_in   = 10
  lumsc_out  = 11
  lucolour   = 12

! output mascon file
  open(lumsc_out,file=msc_output,status='unknown',iostat=ioerr)
  if (ioerr /= 0)then
    call status_update('FATAL','UTIL','recolour_mascons',msc_output,'Error opening output mascon file',0)
  endif

  call status_update('STATUS','UTIL','recolour_mascons',msc_output,'Opened output mascon file',0)

! input recolour file
  open(lucolour,file=recolour_file,status='old',iostat=ioerr)
  if (ioerr /= 0)then
    call status_update('FATAL','UTIL','recolour_mascons',recolour_file,'Error opening input recolour mascon file',0)
  endif

  call status_update('STATUS','UTIL','recolour_mascons',recolour_file,'Opened input recolour mascon file',0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      read the mascon information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call read_msc_hdr(lumsc_in,msc_input,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                                ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
  call allocate_mascon_arrays
  call read_mascon_file(lumsc_in,msc_input)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   read the file containing primary mascon numbers and their new colours
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ioerr = 0
  do while (ioerr == 0)
    read(lucolour,*,iostat=ioerr,end=20)tmp_prim,tmp_colour
    if(ioerr == 0)then
! just change it here, directly
      do itern=1,nint(mcon_prim(tmp_prim,8))
        tmp_tern = mcon_tern(tmp_prim,itern,7)
        tern_colours(tmp_tern,1:2) = tmp_colour
      enddo
      write(message,'(a,i7,a,i5)')"Change colour of primary mascon",tmp_prim," to ",tmp_colour
      call status_update('STATUS','UTIL','recolour_mascons',recolour_file,message,0)
    endif      
  enddo
20 continue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!                                    OUTPUT
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call status_update('STATUS','UTIL','recolour_mascons',msc_output,"Outputting re-coloured mascon file",0)
! add a new line to the header to represent what has been done in running this program
  call date_and_time(date,time,timezone)
  write(msc_hdr_lines(n_hdr_lines+1),'(a)')"# modified by recolour_mascons to change colour of primary mascons."
! generate new unique code for the output program
  call hash_header(msc_hdr_lines(2:n_hdr_lines+1),n_hdr_lines+1,hashcode)
! write the new code and mascon numbers to the first line of the header
  write(msc_hdr_lines(1),'(a,6i8)')hashcode,max_prim,max_sec,max_tern &
                                ,max_sec_per_prim,max_tern_per_prim,max_tern_per_sec
! and write out all header lines to the output file
  do i=1,n_hdr_lines+1
    write(lumsc_out,'(a)')msc_hdr_lines(i)
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now the mascons

  do imsc=1,max_prim
! write out the primary mascon
      msc_code = prim_flags(imsc)
      msc_crds(1:2) = mcon_prim(imsc,1:2)*180.d0/pi
      msc_crds(3) = mcon_prim(imsc,3)
      pointers = 0
      colours = 0
      region = mcon_region(imsc)
      call write_mascon_record(lumsc_out,msc_code,imsc,1,nint(mcon_prim(imsc,8)),msc_crds &
            ,mcon_prim(imsc,4),mcon_prim(imsc,5),mcon_prim(imsc,12),mcon_prim(imsc,6),mcon_prim(imsc,11) &
            ,mcon_prim(imsc,10) , region,pointers,colours)
!print*,"wrote out ",msc_code,imsc,new_prim_total,msc_prim_out(imsc,6),msc_prim_out(imsc,7)*100.d0

! PT161030: also write out a secondary mascon line
      pointers(1) = imsc
      pointers(2) = 0
      colours = sec_colour(imsc)
      msc_code = "S"//prim_flags(imsc)(2:6)
      call write_mascon_record(lumsc_out,msc_code,imsc,1,nint(mcon_prim(imsc,8)),msc_crds &
            ,mcon_prim(imsc,4),mcon_prim(imsc,5),mcon_prim(imsc,12),mcon_prim(imsc,6),mcon_prim(imsc,11) &
            ,mcon_prim(imsc,10) , region,pointers,colours)

! and write out the ternary lines
      do itern=1,nint(mcon_prim(imsc,8))
        tmp_tern = mcon_tern(imsc,itern,7)
        msc_code = tern_flag(tmp_tern)
        msc_crds(1:2) = mcon_tern(imsc,itern,1:2)*180.d0/pi
        msc_crds(3) = mcon_tern(imsc,itern,3)
        pointers = imsc
        colours(1) = tern_colours(tmp_tern,1)
        colours(2) = tern_colours(tmp_tern,1)
        call write_mascon_record(lumsc_out,msc_code,nint(mcon_tern(imsc,itern,7)),1,nint(mcon_tern(imsc,itern,8)),msc_crds &
            ,mcon_tern(imsc,itern,4),mcon_tern(imsc,itern,5),mcon_tern(imsc,itern,8),mcon_tern(imsc,itern,6) &
            ,mcon_prim(imsc,11),mcon_prim(imsc,10), region,pointers,colours)
      enddo

  enddo  ! end of loop over writing out primary mascons.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  call status_update('STATUS','UTIL','recolour_mascons',' ',"End of recolour_mascons",0)



  end
