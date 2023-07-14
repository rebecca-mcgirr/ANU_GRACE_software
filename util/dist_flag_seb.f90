      Program dist_flag

! Program to compute the minimum distance between any two ternary mascons in any two primary mascons.
! This program creates the dist_flag file that we use in GRACEORB to know whether to break into 
! secondary or ternary mascons.
!
! Program reads the new mascon file that contains all info on primary, secondary and ternary mascons.
!
! Based on code "distance.f" written by Tony Purcell (so blame him!!)
!
! P. Tregoning
! 7 September 2016
!
! MODS:
! PT181031: output the 20 closest primary mascon pairs
! PT101102: provide the option of just outputting information regarding the nearest primary mascon (land and water) of each primary
! PT190321: get minimum number of ternarys and minimum centre of mass separation values from command line

  use mascon_mod

  implicit none

  character*150 :: mascon_file, dist_flag_file      ! input and output files
  integer*4     :: lumsc                            ! unit number of mascon file

! variables to set thresholds for identifying problems with primary mascons
  integer*4     :: min_tern                         ! minimum number of ternarys permitted in a primary mascon
  real(kind=8)  :: min_dist_sep                     ! minimum permitted separation of centres of mass of neighbouring primary mascons

! variables related to outputting mascon info rather than the dist_flag file
  character*4   :: info_flag                        ! flag to indicate whether to output the dist_flag file or just primary mascon info
  integer*4     :: closest_msc(3)                   ! closest mascons: land, ocean, same_type
  real(kind=8)  :: closest_dist(3)                  ! associated distances to the closest mascons of the three types
  integer*4     :: max_do_loop                      ! variable to set how many mascons we should loop through (different for "info" and full run)

! counters and pointers
  integer*4     :: iprim1,iprim2,itern1,itern2

! coordinate variables
  real(kind=8),dimension(3) :: xyz_1,xyz_2,xyz_1a,xyz_2a
  real(kind=8)              :: base_dist,dist,min_dist
  real(kind=8)              :: sitelat,sitelon

! variables for determining the closest mascon pairs
  real(kind=8), allocatable :: msc_to_msc_dist(:,:)        ! distances between mascons and the two primary mascon numbers

! local variables
  real(kind=8)  :: threshold(3)
  real(kind=8)  :: pi,rad_fact
  logical       :: debug
  character*250 :: message,arg
  real(kind=8)  :: amag3

  pi = 4.0d0 * atan(1.0d0)
  rad_fact = pi/180.0d0
  debug = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT181102: set thresholds by which decisions are made as to whether to use primary/secondary/ternary mascons
  threshold(1) = sqrt(2.0 * (1.0 - cos(25. * rad_fact)))
  threshold(1) = threshold(1) * 6371200. + sqrt(2.0) * 330000.
  threshold(2) = sqrt(2.0 * (1.0 - cos(20. * rad_fact)))
  threshold(2) = threshold(2) * 6371200. + sqrt(2.0) * 330000.
! PT161130: change ternary/secondary threshold from 5 deg up to 19 deg while we don't have any secondary mascons in the mascon file
  threshold(3) = sqrt(2.0 * (1.0 - cos(19.0d0 * rad_fact)))
  threshold(3) = threshold(3) * 6371200. + sqrt(2.0) * 330000.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! decode runstring
  call getarg(1,mascon_file)
  if(mascon_file(1:1) == " ")then
    write(message,'(a)')"Runstring: dist_flag input_mascon_file output_flag_file [info]"
    call status_update('FATAL','UTIL','dist_flag',mascon_file,message,0)
  endif

  call getarg(2,dist_flag_file)
  call getarg(3,info_flag)
  if(info_flag == "info")then
    call status_update('STATUS','UTIL','dist_flag',mascon_file,"Will output only information about nearest primary mascons",0)
    call getarg(4,arg)
    read(arg,*)min_tern
    call getarg(5,arg)
    read(arg,*)min_dist_sep
  else
    min_tern = 100
    min_dist_sep = 150.d0
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read through primary mascon file to determine number of mascon cells and
! convert coords to XYZ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  lumsc=21
  call status_update('STATUS','UTIL','dist_flag',mascon_file,"Read in the mascons",0)
  call read_msc_hdr(lumsc,mascon_file,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                                ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
  call allocate_mascon_arrays
  call read_mascon_file(lumsc,mascon_file)
  call calc_mascon_xyz("N")
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize some variables
  allocate(msc_to_msc_dist(total_prim,total_prim))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! open the output file if required
  if (info_flag /= "info")then
    open(9,file=dist_flag_file)
    write(9,'(a)')msc_hdr_lines(1)        ! write out the first header line of the mascon file, which contains the unique code for the mascon file
    write(9,*) total_prim                 ! the second line is the total number of primary mascons
    write(9,'(i1)') 3                     ! not sure what this one is ....  is it the flag for primary1 wrt primary1 ?
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! loop over all primary mascons, and all ternary mascons within, outputting flags
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call status_update('STATUS','UTIL','dist_flag',mascon_file,"Start loop over all the mascons",0)
! PT190406: start at mascon 1, so we get information for them all
  do iprim1 = 1,total_prim              ! total_prim assigned value when mascon file is read in read_mascon_file
    closest_dist = 1.e9      ! a ridiculously large distance
    closest_msc  = 2.e6      ! a number higher than the number of ternary mascons

    xyz_1 = mcon_prim(iprim1,1:3)
    if(mod(iprim1,300) == 0 .and. info_flag /= "info")then
      write(message,'(a,i7,a,i7,a)')'Mascon ',iprim1,' of',total_prim,' in file '
      call status_update('STATUS','UTIL','dist_flag',mascon_file,message,0)
    endif

    ! set the max counter value for the do loop over mascons 
    if(info_flag == "info")then
      max_do_loop = total_prim
    else
      max_do_loop = iprim1-1
    endif

    do iprim2 = 1,max_do_loop
      xyz_2 = mcon_prim(iprim2,1:3)
!     calculate the distance between them
      base_dist = sqrt((xyz_1(1)-xyz_2(1))**2+(xyz_1(2)-xyz_2(2))**2+(xyz_1(3)-xyz_2(3))**2)
! PT181031: store this distance so that we can identify the closest pairs of primary mascons
      msc_to_msc_dist(iprim1,iprim2) = base_dist
! PT181031: identify the closest pairs of primary mascons. Let's set a limit of less than 150km, then output all pairs within it
      if(msc_to_msc_dist(iprim1,iprim2) < min_dist_sep*1.e3 .and. iprim2 /= iprim1)then
        write(message,'(a,i6,i6,f15.3,a,2f8.1)')"Primary mascon centres of mass closer than 150 km together: " &
                                         ,iprim1,iprim2,base_dist/1.e3," km separation" &
                                         ,mcon_prim(iprim1,6),mcon_prim(iprim2,6)
        call status_update('STATUS','UTIL','dist_flag',mascon_file,message,0)
      endif
! PT181102: output also if there are fewer than 100 ternary mascons in the primary
      sitelat = dasin(mcon_prim(iprim1,3)/amag3(mcon_prim(iprim1,1:3)))*180.d0/pi
      sitelon = datan(mcon_prim(iprim1,2)/mcon_prim(iprim1,1))*180.d0/pi
      if(mcon_prim(iprim1,2) >=0.d0 .and. mcon_prim(iprim1,1) < 0.d0)then
        sitelon = 180.d0+sitelon
      elseif(mcon_prim(iprim1,2)  < 0.d0 .and. mcon_prim(iprim1,1) < 0.d0)then
        sitelon = 180.d0 + sitelon  
      elseif(mcon_prim(iprim1,2) < 0.d0 .and. mcon_prim(iprim1,1) > 0.d0)then
        sitelon = 360.d0 + sitelon
      endif
      if(mcon_prim(iprim1,8) < min_tern .and. iprim1 == iprim2)then ! only print this out once per primary mascon
        write(message,'(a,i7,a,i7,a,i3,a,f10.3,a,f10.3,a)')"Primary mascon ",iprim1," has fewer than",min_tern &
            ," ternary mascons: ", nint(mcon_prim(iprim1,8)),' (',sitelon,' , ',sitelat,')'
        call status_update('STATUS','UTIL','dist_flag',mascon_file,message,0)
      endif


!--------------------------------------------------------------------
      if(iprim1 /= iprim2)then
!       is this the closest mascon of this type
        if(mcon_prim(iprim2,6) < 1010.d0)then   ! it is a land mascon
          if(base_dist < closest_dist(1) )then       !!!.and. prim_flags(iprim2)(1:5) /= "PCasp")then
            closest_msc(1) = iprim2
            closest_dist(1) = base_dist
          endif
        endif
        if(mcon_prim(iprim2,6) >= 1010.d0)then   ! it is an ocean mascon
          if(base_dist < closest_dist(2) )then
            closest_msc(2) = iprim2
            closest_dist(2) = base_dist
          endif
        endif
      endif
!--------------------------------------------------------------------



! PT181102: only do this if we want the full running of dist_flag
      if(info_flag /= "info")then

        if ( base_dist .gt. threshold(1) ) then
          write(9,*) iprim1, iprim2, base_dist, 1 !,base_dist,' primary to primary exceeds threshold 1:',threshold(1)
          if(debug)write(*,*) 1 !,base_dist,' primary to primary exceeds threshold 1:',threshold
        else
          min_dist = base_dist

          do itern1 = 1,int(mcon_prim(iprim1,8))       ! loop over the number of ternary mascons in the first of the primary mascon pair
            xyz_1a = mcon_tern(iprim1,itern1,1:3)
            do itern2 = 1, int(mcon_prim(iprim2,8))    ! loop over the number of ternary mascons in the second of the primary mascon pair
              xyz_2a = mcon_tern(iprim2,itern2,1:3)
              dist = sqrt((xyz_1a(1)-xyz_2a(1))**2+(xyz_1a(2)-xyz_2a(2))**2+(xyz_1a(3)-xyz_2a(3))**2)
               min_dist = min(min_dist, dist)
               if ( min_dist .lt. threshold(3)) goto 100
            enddo
          enddo

100       if ( min_dist .gt. threshold(2) ) then
            write(9,*) iprim1, iprim2, min_dist, 1 !,min_dist,' min_dist exceeds threshold(2)',threshold(2)
            if(debug)write(*,*) 1 !,min_dist,' min_dist exceeds threshold(2)',threshold(2)
          elseif ( min_dist .gt. threshold(3) ) then
            write(9,*)  iprim1, iprim2, min_dist, 2 !,min_dist,' a secondary? ',threshold(3)
            if(debug)write(*,*) 2 !,min_dist,' a secondary? ',threshold(3)
          else
            write(9,*) iprim1, iprim2, min_dist, 3 !,min_dist,' must be a ternary'
            if(debug)write(*,*) 3 !,min_dist,' must be a ternary'
          endif
        endif
      endif      ! end of if statement overjust outputting info

    enddo        ! end of loop over second mascon
    if(info_flag /= "info")then
      write(9,*)  iprim1, iprim2,  0, 3 !,min_dist,' this one is from the end'

    else
!   output info on the closest mascons to this primary
      write(message,'(a,i7,a,f7.1,2(a,i7,f10.1,a))')"Mascon ",iprim1, "  density", mcon_prim(iprim1,6) &
                                          ,"      Closest  land primary: ",closest_msc(1),closest_dist(1)/1.e3," km." &
                                          ,"      Closest ocean primary: ",closest_msc(2),closest_dist(2)/1.e3," km."
      call status_update('STATUS','UTIL','dist_flag',mascon_file,message,0)
    endif

  enddo          ! end of loop over first mascon

  if(info_flag /= "info")close(9)


  stop
  end
