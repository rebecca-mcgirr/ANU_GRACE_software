  program extract_ocean_mascons_JP20190807

! program to read a mascon file and extract out the primary/secondary/ternary mascons that are ocean to make the ocean mascon file. We
! can then reconfigure the ocean mascon file if we so desire ... add/change tidal amplitude flags, change primary ocean mascon shapes/sizes etc
!
! P. Tregoning
! 05 December 2016
! 22 July 2019 modified by J. Pfeffer to conserve several ocean primaries (not all ternaries are merged in a single primary as previously)
! 07 August 2019 corrected by J. Pfeffer to have the same number of primaries and secondaries (same conditions on water density and inclusion of at least one ternary)
! 08 August 2019 corrected by J. Pfeffer to have the correct number of ternaries in ocean primaries (do not count islands) 
  use mascon_mod

  implicit none

  character   :: infile*150,outfile*150,line*150
  integer*4   :: lumsc_in, lumsc_out

  integer*4   :: trimlen

! variables related to writing out the mascons
  real(kind=8) :: msc_crds(3)
  real(kind=8) :: pi
  integer*4    :: pointers(2),colours(2)
  integer*4    :: itern,isec,nmsc_current,imsc,nsec,ntern
  REAL , DIMENSION(:), allocatable :: ntern_ocean
  character    :: date*8,time*10,timezone*5,hashcode*12

! variables to count the ocean mascons
  integer*4    :: new_n_prim                       ! total number of ocean primaries mascons
  integer*4    :: new_n_sec                        ! total number of ocean secondary mascons
  integer*4    :: n_water_tern                            ! number of water ternarys within a mixed/ocean primary mascon
  integer*4    :: total_ocean_ternarys                    ! total number of ocean ternary mascons
  integer*4    :: new_max_tern_per_prim                   ! maximum number of ternaries in an ocean mascon

! local variables
  character*250  :: message

  pi = 4.d0*datan(1.d0)

! define unit numbers
  lumsc_in  = 10
  lumsc_out = 11

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! decode runstring and open files
  call getarg(1,infile)
  if(infile(1:1) == " ")then
    write(message,'(a)')"Runstring: extract_ocean_msc mascon_infile mascon_outfile"
    call status_update('FATAL','UTIL','extract_ocean_msc',' ',message,0)
  endif
  write(message,'(a)')" Running new version of extract ocean mascons modified by J. Pfeffer on August 07, 2019"
  call status_update('STATUS','UTIL','extract_ocean_msc','',message,0)
  call getarg(2,outfile)
  open(lumsc_out,file=outfile,status='unknown')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read input mascon file header, allocate mascon arrays and read in all arrays
  call read_msc_hdr(lumsc_in,infile,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                                ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
  call allocate_mascon_arrays
  allocate(ntern_ocean(max_prim))
  call read_mascon_file(lumsc_in,infile)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! loop through all the primary mascons to identify how many primary mascons
! contain ocean ternarys and how many ocean ternarys there are in total
!  total_ocean_prim=0
!  total_ocean_sec=0
  total_ocean_ternarys = 0
  new_max_tern_per_prim = 0
  do imsc = 1,total_prim
    if(mcon_prim(imsc,11) < 100.)then
      if(mcon_tern(imsc,1,6) > 1000.) then
        new_n_prim=new_n_prim+1 
        do isec=1,nint(mcon_prim(imsc,7))
          new_n_sec=new_n_sec+1
          n_water_tern=0
          do itern = 1,nint(mcon_sec(imsc,1,8))
            if(mcon_tern(imsc,itern,6) > 1000.) then
              total_ocean_ternarys  = total_ocean_ternarys  + 1  
              n_water_tern = n_water_tern +1
            endif    
          enddo
          ntern_ocean(imsc)=n_water_tern
      	enddo   
        if(nint(mcon_prim(imsc,8)) > new_max_tern_per_prim)then
          new_max_tern_per_prim = mcon_prim(imsc,8)
          print *, 'Maximum number of ternaries in primaries',new_max_tern_per_prim
          print *, 'Region : ',mcon_region(imsc)
        endif  
      endif  
      
    endif
  enddo
  write(message,'(i7,a,i7,a)')new_n_prim," of a total of ",max_prim," primary mascons are ocean mascons."
  call status_update('STATUS','UTIL','extract_ocean_mascons',' ',message,0)
  write(message,'(i7,a,i7,a)')new_n_sec," of a total of ",max_sec," secondary mascons are ocean mascons."
  call status_update('STATUS','UTIL','extract_ocean_mascons',' ',message,0)
  write(message,'(i7,a,i7,a)')total_ocean_ternarys," of a total of ",max_tern," ternary mascons are ocean mascons."
  call status_update('STATUS','UTIL','extract_ocean_mascons',' ',message,0)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                          O U T P U T    A L L    O C E A N    M A S C O N S
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! update the numbers in the first line of the header for the mascon file
!  write(msc_hdr_lines(1)(13:57),'(2i6,i9,i6,2i9)')1,1,total_ocean_ternarys,1,total_ocean_ternarys,total_ocean_ternarys
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! add a new line to the header to represent what has been done in running this program
  call date_and_time(date,time,timezone)
  write(msc_hdr_lines(n_hdr_lines+1),'(a,a,a,1x,a8,1x,a10,1x,a5)') &
       "# extract_ocean_mascons from file ",infile(1:trimlen(infile)) &
       ," Timetag: ",date,time,timezone
! generate new unique code for the output program
  call hash_header(msc_hdr_lines(2:n_hdr_lines+1),n_hdr_lines+1,hashcode)

! write the new code and mascon numbers to the first line of the header
  write(msc_hdr_lines(1)(1:81),'(a,i6,5i8)')hashcode(1:8),new_n_prim,new_n_sec,max_tern,max_sec_per_prim &
                                           ,new_max_tern_per_prim,new_max_tern_per_prim
! write the header lines out to the mascon output file
  write(lumsc_out,'(a)')msc_hdr_lines(1:n_hdr_lines+1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write out the primary and secondary records for the ocean mascon file
! PT170531: added a geoid height of "0.d0" to this primary record. Changed the number of ternarys to be only water ternarys
  write(message,'(a,i6,a)')"Writing out ",new_n_prim," ocean primary mascons"
  call status_update('STATUS','UTIL','extract_ocean_mascons',outfile,message,0)
  nmsc_current = 0
  ! now, loop through and output only the water ones
  do imsc = 1,total_prim 
    if(mcon_prim(imsc,11) < 100.)then	
      if(mcon_tern(imsc,1,6) > 1000.) then
        !print *, "Writing out ",prim_flags(imsc)," ocean primary mascons"
        nmsc_current = nmsc_current+ 1
	    msc_crds(1:2) = mcon_prim(imsc,1:2)*180.d0/pi
	    msc_crds(3)   = mcon_prim(imsc,3)
	    call write_mascon_record(lumsc_out,prim_flags(imsc),nmsc_current,1,int(ntern_ocean(imsc)),msc_crds     &
	   					  ,mcon_prim(imsc,4),mcon_prim(imsc,5),0.d0,mcon_prim(imsc,6),mcon_prim(imsc,11) &
						  ,mcon_prim(imsc,10),mcon_region(imsc),pointers,colours)
!     secondary record for this primary. We still can only have one secondary per primary at this point.
	    pointers(1) = nmsc_current
	    pointers(2) = 0
	! PT161121: we don't want this, we want the secondary mascons to be sequential along with the primary mascons
	!      nsec = mcon_sec(imsc,1,7)           ! this is the unique secondary mascon number as read from the input mascon file
	    colours = 0
	    nsec = 1
	    call write_mascon_record(lumsc_out,sec_flags(imsc),nmsc_current,nmsc_current,int(ntern_ocean(imsc)),msc_crds &
							  ,mcon_prim(imsc,4),mcon_prim(imsc,5),0.d0,mcon_prim(imsc,6),mcon_prim(imsc,11) &
							  ,mcon_prim(imsc,10),mcon_region(imsc),pointers,colours)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! ternary records for this secondary
        do itern = 1,nint(mcon_prim(imsc,8))
	      pointers(1) = nmsc_current
	      pointers(2) = nsec
	  ! PT170531: only water ternarys
		  if(mcon_tern(imsc,itern,6) > 1000.) then
		    ntern = nint(mcon_tern(imsc,itern,7))
		    msc_crds(1:2) = mcon_tern(imsc,itern,1:2)*180.d0/pi
		    msc_crds(3)   = mcon_tern(imsc,itern,3)
		    call write_mascon_record(lumsc_out,tern_flag(ntern),ntern,nsec,int(ntern_ocean(imsc)),msc_crds &
							  ,mcon_tern(imsc,itern,4),mcon_tern(imsc,itern,5),0.d0,mcon_tern(imsc,itern,6),mcon_prim(imsc,11) &
							  ,mcon_prim(imsc,10),mcon_region(imsc),pointers,tern_colours(ntern,:))
           endif
        enddo
      endif  
    endif
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(message,'(a,i7,a)')"Written ",nmsc_current," primary ocean mascons out to file"
  call status_update('STATUS','UTIL','extract_ocean_mascons',outfile,message,0)
  end
  


