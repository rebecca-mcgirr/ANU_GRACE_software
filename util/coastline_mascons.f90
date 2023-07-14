  program coastline_mascons

! program to take the mascons_stage3 mascon file and merge the small segments of mascons adjacent to coastlines
! into neighbouring mascons of similar densities. The intention is to maintain as much as possible the N/S boundaries
! of the east and west edges of the mascons, since this seems to affect the accuracy with which the primary mascon
! signals can be estimated.
!
! P. Tregoning (+A. Purcell, R McGirr, S. Allgeyer, H. McQueen)
! 1 September 2021

  use mascon_mod   ! provides all the mascon arrays


  implicit none

! command line arguments
  character*150 :: mascon_file_in,mascon_file_out      ! input and output mascon file names
  integer*4     :: min_tern_per_prim                   ! minimum number of ternaries in a primary before which it is merged

! local variables
  integer*4, parameter :: lumsc_in=10,lumsc_out=11
  integer*4            :: imsc
  character*100        :: arg
  character*250        :: message
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! decode runstring
  call getarg(1,mascon_file_in)
  if(mascon_file_in(1:1) == "")then
    call status_update('FATAL','UTIL','coastline_mascons',' ' &
      ,"Runstring: coastline_mascons mascons_stage3 mascons_stage3_coast min_tern_per_prim",0)
    stop
  endif
  call getarg(2,mascon_file_out)
  call getarg(3,arg)
  read(arg,*)min_tern_per_prim
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      read the input mascon information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call read_msc_hdr(lumsc_in,mascon_file_in,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                                ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
  call allocate_mascon_arrays
  call read_mascon_file(lumsc_in,mascon_file_in)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      loop through and identify small primary mascons that need to be merged
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do imsc=1,max_prim
    if(mcon_prim(imsc,8) < 0.9*nint(n_tern_per_prim))then
      write(message,'(3(a,i6))')'Primary',imsc,' has too few ternary mascons (',nint(mcon_prim(imsc,8)) &
          ,'). Minimum required is: ',nint(0.9*min_tern_per_prim)
      call status_update('STATUS','UTIL','coastline_mascons',' ',message,0)
    endif
  enddo

  end 


