  program voronoi_mascons

! program to take a set of primary mascons that define a particular region (e.g. a continent, an island or a drainage basin) and 
! reshape them to be roughly equal area and equal shape. This is done by reassigning the ternary mascons of the oversize primary
! mascons to the undersize primary mascons, until they are all roughly the same.
!
! Many of the subroutines and some of the logic replicates what is done in util/configure_mascons
!
! P. Tregoning
! 5 November 2014
!
! PT160826 : The script that creates the command file for this program is hitting a limit of the length of arguments that it can echo to 
!            the file. I will get around this by writing the primary mascon numbers out one line at a time instead of all on the same line.
!            This means that I have to change the way that they are read in within this program.

  implicit none

! command line arguments
  character    :: cmdfile*100                   ! command file containing a) name of ternary file b) # and list of primary mascons to adjust
  real(kind=8) :: req_area                      ! required area of reshaped primary mascons (will dictate the size of the output mascons)
  integer*4    :: ioerr
  logical      :: info_only                     ! flag to indicate whether to stop after outputting info on how many ternary and primary mascons

! command file variables
  character :: msc_tern_file*100                ! name of input ternary mascon file
  integer*4 :: nmsc_prim_adj                    ! number of input primary mascons to be adjusted
  real(kind=8), allocatable :: msc_prim_in(:,:) ! number to adjust, lat/lon/hgt/area/density/%land/#ternaries/tide_adjust?

! output primary and ternary mascon information
  integer*4                 :: nmsc_prim_out      ! number of primary mascons in the recomputed primary mascon configuration. 
  real(kind=8), allocatable :: msc_prim_out(:,:)  ! number of adjusted primary mascon, lat/lon/hgt/area/density/tide_adjust/depth
  real(kind=8), allocatable :: msc_tern_in(:,:,:) ! primary mascon, ternary number, lat/lon/radius/area/density/depth (or height)
  real(kind=8), allocatable :: msc_tern_out(:,:,:)! primary mascon, ternary number, lat/lon/radius/area/density/depth (or height)
  real(kind=8)              :: sum_area           ! total area of the input primary mascons, calcuated from the input ternary mascons
  integer*4                 :: old_prim_number    ! number of output old primary mascon (the reconfigured ones go to the bottom of the file)
  character                 :: msc_tern_outfile*100   ! name of output ternary mascon file
  integer*4                 :: random_colour(10000)   ! array of random numbers used to assign colours when we plot the primary mascons
  character                 :: prim_flag*5            ! type of primary mascon ("Deep ", "Land ", or "Shelf")
  logical                   :: shelf_mascon           ! logical to indicate whether we are treating continental shelf mascons or not

! temp parameters to read the ternary mascon info
  integer*4    :: tmp_pri, tmp_n_tern
  real(kind=8) :: tmp_land, tmplat, tmplon       ! variables to read the input primary mascon information
  integer*4    :: ntern                          ! the sum of all ternary mascons in all the primary mascons to be reshaped
  real(kind=8) :: tmparea,tmpdensity,tmptide     ! variables to read the input primary mascon information
  real(kind=8) :: tmpdepth
  character    :: char1*1, char7*7
  real(kind=8) :: density,tides                  ! based on land/water primary mascon, write out mascon density and tide estimate status
 
! parameters related to looping over the transfer process
  integer*4 :: iter
  integer*4 :: max_iter
  integer*4 :: n_terns(int(1e6))

! parameters concerned with transferring mascons
  integer*4    :: i_big, closest_node
  real(kind=8) :: tmp_big
  real(kind=8), allocatable :: dist_pr_2_ter(:,:) ! distances between primary CoM and all ternarys within the primary. Second col is pointer to the ternary #
  integer*4    :: itern
  real(kind=8), allocatable :: tern_crds(:,:)    ! array to store all ternary mascons, without info on which primary they originally came from
  real(kind=8) :: rand_number
  real(kind=8) :: tern_dist, min_dist
  real(kind=8) :: tmp_tern_in(6)                 ! array to transfer unchanged ternary mascons from input to output file

! unit numbers
  integer*4, parameter      :: lucmd  = 10        ! unit number of command file
  integer*4, parameter      :: lutern = 11        ! unit number of input  ternary mascon file
  integer*4, parameter      :: lutern_out = 12    ! unit number of output ternary mascon file
  integer*4, parameter      :: lutern_region = 13 ! unit number of output file containing only mascons for the requested region

! certain parameters to dimension arrays etc
  integer*4, parameter :: max_tern = 500000      ! maximum number of ternary mascons in the total region
  integer*4, parameter :: max_prim = 5000        ! maximum number of primary mascons that can be reconfigured
  real(kind=8)         :: pi

! counters and logicals
  integer*4 :: i,j,k
  integer*4 :: nfound                            ! count how many of the requested primary mascons have been found in the ternary file
  logical   :: use_primary     
  integer*4 :: counter

! other stuff
  integer*4 :: indx, trimlen
  character :: line*200, message*250,arg*100
  logical   :: debug



  pi = 4.d0*atan(1.d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  get the command line and cmdfile information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call getarg(1,cmdfile)
  if(cmdfile(1:1) == " ")then
    call status_update('WARNING','UTIL','voronoi_mascons',' ','Runstring: voronoi_mascons cmdfile max_area',0)
    stop
  endif
  open(lucmd,file=cmdfile,status='old',iostat=ioerr)
  if(ioerr /= 0)then
    call status_update('FATAL','UTIL','voronoi_mascons',cmdfile,'Error opening command file. Does it exist?',0)
  endif
  call getarg(2,arg)
  read(arg,*)req_area
  call getarg(3,arg)
  read(arg,*)max_iter
! flag on whether to stop after outputting info on how many ternary and primary mascons
  call getarg(4,arg)
  if (arg(1:9) == "info_only")then
    info_only = .true.
  else
    info_only = .false.
  endif

! read the info from the command file
  read(lucmd,'(a)')msc_tern_file  
  read(lucmd,*)nmsc_prim_adj
! allocate the array of input primary mascons to be adjusted
  allocate(msc_prim_in(nmsc_prim_adj,8))
! read the primary mascon numbers of the primary mascons to be adjusted
! PT160826: change this to read one primary mascon number per line (instead of all from the same line)
  do i=1,nmsc_prim_adj
    read(lucmd,*)msc_prim_in(i,1)
  enddo

! generate the name of the output ternary mascon file
  indx = index(msc_tern_file," ")
  msc_tern_outfile = msc_tern_file(1:indx-1)//".reconfigured"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  open the input and output ternary mascon files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  indx = index(msc_tern_file," ")
  open(lutern,file=msc_tern_file(1:index(msc_tern_file," ")),status='old',iostat=ioerr)
  if(ioerr /= 0)then
    call status_update('FATAL','UTIL','voronoi_mascons',msc_tern_file,'Error opening ternary mascon file. Does it exist?',0)
  else
    call status_update('STATUS','UTIL','voronoi_mascons',msc_tern_file(1:index(msc_tern_file," ")) &
                         ,'Opened input ternary mascon file',0)
  endif
! now the output file
  open(lutern_out,file=msc_tern_outfile(1:index(msc_tern_outfile," ")),status='unknown',iostat=ioerr)

! now the output file for mascons of just this region
  open(lutern_region,file="region_mascons",status='unknown',iostat=ioerr)

! allocate the ternary mascon array
  allocate(msc_tern_in(nmsc_prim_adj,max_tern,7))    ! primary #, ternary #, lat/lon/radius/area/density/depth (or height)
  allocate(dist_pr_2_ter(max_tern,2))                  ! distances between the primary CoM and all ternaries within the primary
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   read through the input ternary file to extract out the required primary/ternary information for !
!   the primary mascons to be reshaped                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ioerr = 0
  nfound = 0
  ntern = 0
  old_prim_number = 0
  shelf_mascon = .false.
  do while (ioerr == 0 )
    read(lutern,*,iostat=ioerr,end=1000)tmp_pri,char7,tmp_n_tern,tmp_land,tmplat,tmplon,tmparea,tmpdepth,tmpdensity,tmptide
!print*,tmp_pri,char7,tmp_n_tern,tmp_land,tmplat,tmplon,tmparea,tmpdepth,tmpdensity,tmptide
    if(ioerr == 0 )then
!     check whether we want to adjust this primary mascon
      use_primary = .false.
      do j=1,nmsc_prim_adj
        if(tmp_pri == msc_prim_in(j,1))then ! it is one of the primary mascons that we want to adjust. Save the information
          use_primary = .true.
! PT141117: is it the continental shelf mascon that needs to be divided?
          if(char7 == "PShelf ")then
            print*,"It is the continental shelf mascon"
            shelf_mascon = .true.
          endif
          nfound = nfound + 1
          msc_prim_in(j,2) = tmplat
          msc_prim_in(j,3) = tmplon
          msc_prim_in(j,6) = tmp_land
          msc_prim_in(j,7) = tmp_n_tern
          msc_prim_in(j,4) = 0.0   ! set the primary mascon area to zero
          random_colour(j) = nint(rand()*dble(nmsc_prim_adj)*2)
          do k=1,tmp_n_tern
            ntern = ntern + 1
            read(lutern,*)msc_tern_in(j,k,1:6)
! PT160825: write the ternarys to the screen if info_only
            if(info_only)print*,msc_tern_in(j,k,1:2),random_colour(j)," included ternary"
          enddo
        endif
      enddo
      if(.not.use_primary)then
!       skip over all the ternaries for this primary mascon - we don't need them - but write them out
        old_prim_number = old_prim_number + 1
        write(lutern_out,*)old_prim_number," ",char7,tmp_n_tern,tmp_land,tmplat,tmplon,tmparea,tmpdepth,tmpdensity,tmptide
        do i=1,tmp_n_tern
          read(lutern,'(a)')line
          write(lutern_out,'(a)')line
        enddo
      endif

    endif
  enddo

1000 continue

! we now know how many ternary mascons need to be divided up. Allocate an array and then transfer them
! all into one array
  allocate(tern_crds(ntern,7))
  counter = 0
  do j=1,nmsc_prim_adj
    do k=1,int(msc_prim_in(j,7))
      counter = counter + 1
      tern_crds(counter,1:6) = msc_tern_in(j,k,1:6)
    enddo
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   now, based on the areas of each requested mascon and the input area, how many output primary    !
!   mascons do we need                                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  sum_area = sum(tern_crds(1:ntern,4))
  nmsc_prim_out = nint(sum_area/req_area)
print*,'sum_area,req_area,nmsc_prim_out',sum_area,req_area,nmsc_prim_out
  write(message,'(a,i8,a,e12.6)')'There  are    a total of ',ntern,' ternary mascons. Total area (km^2) = ',sum_area/1.d6
  call status_update('STATUS','UTIL','voronoi_mascons',' ',message,0)
  write(message,'(a,i8,a,i6,a)')"There will be a total of ",nmsc_prim_out,' reshaped primary mascons (from ',nmsc_prim_adj &
                                   ,' input primary mascons)'
  call status_update('STATUS','UTIL','voronoi_mascons',' ',message,0)
! stop if that's all we wanted to know
  if(info_only)stop

! allocate the output primary mascon array
  allocate(msc_prim_out(nmsc_prim_out,7))
  allocate(msc_tern_out(nmsc_prim_out,max_tern,7))    ! primary #, ternary #, lat/lon/radius/area/density/depth (or height)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   generate the appropriate number of seed points to start the nucleation process      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(nmsc_prim_out == nmsc_prim_adj)then
!   same number of input and output primary mascons, use the input primary mascon coords as the seed points.
    msc_prim_out(:,1:2) = msc_prim_in(:,2:3)
  else if (nmsc_prim_out < nmsc_prim_adj) then
!   fewer output than input primary mascons, use the first X input primary mascon coords as the seed points.
    msc_prim_out(1:nmsc_prim_out,1:2) = msc_prim_in(1:nmsc_prim_out,2:3)
  else if (nmsc_prim_out > nmsc_prim_adj) then
!   more output than input primary mascons, we need all input primary mascon coords plus some extras.
    msc_prim_out(1:nmsc_prim_adj,1:2) = msc_prim_in(:,2:3)
!    write(message,'(a,2i6,a)')'Need more primary mascons than there were for this region (',nmsc_prim_out , nmsc_prim_adj,')'
!    call status_update('STATUS','UTIL','voronoi_mascons',' ',message,0)
    do i=nmsc_prim_adj+1,nmsc_prim_out
      call srand(i)
! PT141107: I don't know why I need to do this twice, but doing it only once results in a value of zero!!!
      rand_number = rand()
      rand_number = rand()
      msc_prim_out(i,1:2) = tern_crds(int(rand_number*ntern),1:2)
    enddo
  endif

! PT141113: assign random numbers to be used for mascon colours when plotting
  do i=1,nmsc_prim_out
    random_colour(i) = nint(rand()*dble(nmsc_prim_out))
    random_colour(i) = nint(rand()*dble(nmsc_prim_out)*2)
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   loop through all the ternary mascons and assign each to the closest primary coordinate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! reset to zero the output ternary mascon information
!  msc_tern_out = 0.0

  do iter = 1,max_iter 
    n_terns = 0
    do i=1,ntern
! calculate the distance from this mascon to each of the primary seed coordinates
      min_dist = 1.d10
      do j=1,nmsc_prim_out
        call MASCON_msc_to_msc(tern_crds(i,1),tern_crds(i,2),msc_prim_out(j,1),msc_prim_out(j,2),tern_dist)
        if(tern_dist < min_dist)then
          closest_node = j
          min_dist = tern_dist
        else
        endif
      enddo
!   assign this ternary to the closest seed coordinate
      n_terns(closest_node) = n_terns(closest_node) + 1
      tern_crds(i,7) = closest_node
      msc_tern_out(closest_node,n_terns(closest_node),:) = tern_crds(i,:)
! DEBUG:
! output all the ternaries and the code of their new primary mascon
!      if(iter == max_iter)print*,tern_crds(i,1:2),tern_crds(i,7),' is the new primary. Iteration',iter &
!                           ,' Colour ',random_colour(int(tern_crds(i,7)) ),i
    enddo

! calculate the new centres of mass of the primary mascons
    call status_update('STATUS','UTIL','voronoi_mascons',' ',"Calculating new centres of mass",0)

! !!$OMP PARALLEL DO  shared(msc_tern_out,msc_prim_out, n_terns)
    do j=1,nmsc_prim_out
!print*,'calculating centre of mass of reconfigured mascon ',j
! assign how many ternary mascons there are in each primary
      if(mod(j,10) == 0)then
        write(message,'(a,i6,a,i6,a,i6,a,i6)')'Calculating new CoM for reconfigured mascon',j,' of',nmsc_prim_out &
                                     ,' reconfigured mascons. Iteration ',iter,' of',max_iter
        call status_update('STATUS','UTIL','voronoi_mascons',' ',message,0)
      endif
      msc_prim_out(j,4) = n_terns(j)
      call MASCON_calc_CoM(max_tern,msc_tern_out(j,:,:),msc_prim_out(j,:))
    enddo
! !!$OMP END PARALLEL DO 

  enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   write out the information on the newly configured primary mascons if more than the original
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  debug = .false.
    do i=1,nmsc_prim_out
      if(mod(i,10) == 0)then
        write(message,'(a,i6,a,i6,a)')'Writing out reconfigured mascon',i,' of',nmsc_prim_out,' reconfigured mascons'
        call status_update('STATUS','UTIL','voronoi_mascons',' ',message,0)
      endif
      if(int(msc_prim_out(i,4)) > 0)call MASCON_percent_land(debug,int(msc_prim_out(i,4)),max_tern,msc_tern_out(i,:,:) &
                               ,msc_prim_out(i,3))
! assign a density based on whether it is land or water
      if(msc_prim_out(i,3) < 0.5)then
        density = 1029
        tides = 0
        prim_flag = "Deep "
      else
        density = 1000
! based on mean depth, decide whether to turn on the estimation of tidal constituents. (#&^% NOT IMPLEMENTED ON YET
        tides = 0
        prim_flag = "Land "
      endif
! if it is the continental shelf mascon being broken up, set values accordingly
      if(shelf_mascon)then
        prim_flag = "Shelf"
        tides = 31
        density = 1029
      endif
      write(lutern_out,*)old_prim_number+i," P",prim_flag," ",int(msc_prim_out(i,4)),msc_prim_out(i,3)*100.d0 &
                            ,msc_prim_out(i,1:2),msc_prim_out(i,5),msc_prim_out(i,7),density,tides
      write(*,*)old_prim_number+i," P",prim_flag," ",int(msc_prim_out(i,4)),msc_prim_out(i,3)*100.d0 &
                            ,msc_prim_out(i,1:2),msc_prim_out(i,5),msc_prim_out(i,7),density,tides
! PT160829: also write it out to a separate file that contains only the elements for the requested region
      write(lutern_region,*)old_prim_number+i," P",prim_flag," ",int(msc_prim_out(i,4)),msc_prim_out(i,3)*100.d0 &
                            ,msc_prim_out(i,1:2),msc_prim_out(i,5),msc_prim_out(i,7),density,tides
      do itern = 1,int(msc_prim_out(i,4))
        write(lutern_out,*)msc_tern_out(i,itern,:)
        write(lutern_region,*)msc_tern_out(i,itern,:)
      enddo
    enddo      

  call status_update('STATUS','UTIL','voronoi_mascons',' ',"End of voronoi_mascons",0)

  end





