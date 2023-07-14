!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subroutines for mascon manipulations
!
! P. Tregoning
! 27 October 2014
!
! Contains:
!             MASCON_msc_to_msc          :  calculates the distance between two mascons
!             MASCON_transfer_ternary    :  transfers a ternary mascon from one primary to another
!             MASCON_percent_land        :  calculates the percentage of land in a primary mascon, from the ternary mascon information
!             MASCON_nearest_mscs        :  finds the 10 nearest primary mascons to a particular primary mascon
!             MASCON_read_seed           :  reads the input seed file for new primary mascons
!             MASCON_primary_perimeter   :  identifies which of the ternary mascons lie on the perimeter of the primary mascon
!             voronoi_reshape            :  reshape a selection of primary mascons into a particular size
!
!             make_tern_array            :  combine a subset of ternarys (and their info) into one single array
!             assign_ternarys            :  assign available ternary mascons to their nearest primary CoM
!             jiggle_mascons             :  evens out primary mascon sizes using subsets of mascons
!             voronoi_reshape_v3         :  cut-down version of voronoi_reshape to just do some of the tasks
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine voronoi_reshape_v3(nprim_to_make,ntern_to_reshape,nmsc_to_use,fix_prims,sum_terns,tern_crds &
                               ,msc_reshaped_prims,msc_tern_reshaped)

! subroutine that will:
! 1. select "nmsc_to_use" random ternarys to act as the a priori CoM values
! 2. run a voronoi cell algorithm to reshape this subset of primary mascons into equal area shapes
! 3  compute the new centres of mass, percent land etc and apply all the attributes 
! 4. return the newly configured primary mascon information in msc_reshaped_prims
!
! I will try to use a different approach to calculating the new CoM and percent land to speed up the 
! calculations ... the existing subroutines are far to slow for the simple computations that are being done!
!
! P. Tregoning
! 14 February 2020

  use mascon_mod

  implicit none

! passed variables
  integer*4    ,intent(in)   :: nprim_to_make         ! total number of primaries to be reshaped
  integer*4    ,intent(in)   :: ntern_to_reshape      ! used to dimension the number of ???
  integer*4    ,intent(in)   :: nmsc_to_use           ! number of subset primaries to be reshaped
  integer*4    ,intent(in)   :: fix_prims(nmsc_to_use)! the numbers of the primary mascons to be averaged
  integer*4    ,intent(in)   :: sum_terns             ! number of ternarys in the subset of primaries
  real(kind=8) ,intent(in)   :: tern_crds(max_tern,nvar_tern) ! single array containing subset ternarys
  real(kind=8) ,intent(inout):: msc_reshaped_prims(nprim_to_make,nvar_prim)  ! all primary mascon information
  real(kind=8),intent(inout) :: msc_tern_reshaped(nprim_to_make,max_tern_per_prim,nvar_tern)   ! the new allocation of ternarys to primarys. Mimics mcon_tern array.

! local variables
  integer*4   :: iprim,itern,extra_prim
  real(kind=8),allocatable :: new_prims(:,:)          ! temporary array for the newly shaped subset primaries
  real(kind=8),allocatable :: new_terns(:,:,:)        ! temporary array for the newly shaped subset ternarys
  logical,     allocatable :: rand_numbers(:)         ! array of flags to say whether a particular random number has been used or not
  real(kind=8)             :: shortest_dist(2)        ! array to sort ternary-to-CoM distances in ascending order
  real(kind=8)             :: RMS,tmp_dist,pi
  integer*4                :: iter                    ! number of iterations to converge the reshaping
  character*200            :: message
  logical                  :: new_CoM
  real(kind=8)             :: old_RMS,diff_RMS

  allocate(new_prims(nmsc_to_use,nvar_prim))
  allocate(new_terns(nmsc_to_use,max_tern_per_prim,nvar_tern))
  allocate(rand_numbers(sum_terns))

  pi = 4.d0*datan(1.d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  new_CoM = .false.
  if(new_CoM)then
  ! select "nmsc_to_use" random ternarys to act as a priori CoMs
    call srand(sum_terns)
    rand_numbers(1:sum_terns) = .false.
    do iprim = 1,nmsc_to_use
      call get_random_value(rand_numbers,sum_terns,extra_prim)
      new_prims(iprim,1:3) = tern_crds(extra_prim,1:3)   ! transfer the coords of the ternary to the primary
      new_prims(iprim,4:12) = 0.d0                       ! set all the rest to zero
    enddo
  else
  ! just use the CoMs that were transferred in
    do iprim=1,nmsc_to_use
      new_prims(iprim,1:3) = msc_reshaped_prims(fix_prims(iprim),1:3)
    enddo
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now the loop to convergence, assigning ternarys to CoMs, then recalculating etc
! calculate the distance of each ternary to these 
  RMS = 1.d5
  old_RMS = 0.d0
  diff_RMS = 1.e6
  iter = 1

  write(message,*)"Fix_prims: ",fix_prims(1:nmsc_to_use)
  call status_update('STATUS','UTIL','voronoi_reshape_v3',' ',message,0)
  write(message,*)"Original number of ternarys:",(nint(msc_reshaped_prims(fix_prims(iprim),8)),iprim=1,nmsc_to_use)
  call status_update('STATUS','UTIL','voronoi_reshape_v3',' ',message,0)


  do while ( iter < 5 .or. (iter < 20 .and. dabs(diff_RMS) > 0.5d0) )
    new_prims(:,8) = 0.d0

    do itern = 1,sum_terns
      shortest_dist = 1.d6
      do iprim=1,nmsc_to_use
        call MASCON_msc_to_msc(new_prims(iprim,1)*180.d0/pi,new_prims(iprim,2)*180.d0/pi &
                            ,tern_crds(itern,1)*180.d0/pi,tern_crds(itern,2)*180.d0/pi, tmp_dist )
!print*,'itern,iprim,tmp_dist',itern,iprim,tmp_dist
        if(iprim == 1 .or. tmp_dist < shortest_dist(2))then
          shortest_dist(2) = tmp_dist   ! this is now the distance to the closest CoM
          shortest_dist(1) = iprim      ! this is now the closest CoM
        endif
      enddo
      ! assign the ternary to the closest CoM
      new_prims(nint(shortest_dist(1)),8) = new_prims(nint(shortest_dist(1)),8) + 1
      new_terns(nint(shortest_dist(1)),nint(new_prims(nint(shortest_dist(1)),8)),:) = tern_crds(itern,:)
!print*,'itern,nint(shortest_dist(1)),nint(shortest_dist(1)),8)',itern,nint(shortest_dist(1)),new_prims(nint(shortest_dist(1)),8)
!print*,'tern info:',new_terns(nint(shortest_dist(1)),nint(new_prims(nint(shortest_dist(1)),8)),:)
!print*,'tern_crds(itern,:)',tern_crds(itern,:)
    enddo
!   print a comment of the number of ternarys in each mascon sent to this subroutine
    write(message,*)"Original number of ternarys:",nint(new_prims(1:nmsc_to_use,8)) &
                     ," Mean: ",nint(sum(new_prims(1:nmsc_to_use,8))/nmsc_to_use)
    !if(iter==1)call status_update('STATUS','UTIL','voronoi_reshape_v3',' ',message,0)

    ! we now need to recalculate the CoM of each of the new primary mascons and update the information into the "new_prims" array
    do iprim = 1,nmsc_to_use
      !write(message,'(a,i7,a,3i7)')"Calculating CoM for mascon",iprim," iteration ",iter,nint(new_prims(iprim,8)) &
      !                             , nint(new_prims(iprim,8) - sum_terns/nmsc_to_use)
      !call status_update('STATUS','UTIL','voronoi_reshape_v3',' ',message,0)
      call MASCON_calc_CoM(max_tern_per_prim,nvar_prim,nvar_tern,new_terns(iprim,:,:),new_prims(iprim,:))
      RMS = RMS + (new_prims(iprim,8) - sum_terns/nmsc_to_use)**2
    enddo

    ! calculate the metric - perhaps just use the RMS of the difference of actual ternarys per primary compared to the ideal number (i.e. the average)
    RMS = dsqrt(RMS/dble(nmsc_to_use))
    write(message,'(a,i6,a,i4,a,f15.1)')'Iteration: ',iter,'.   RMS of ',nmsc_to_use,' ternary differences from nominal: ',RMS   
    if(fix_prims(1) ==242)call status_update('STATUS','UTIL','voronoi_reshape_v3',' ',message,0)
    iter = iter+1    
    diff_RMS = RMS-old_RMS
    old_RMS = RMS
  enddo

  ! transfer from new_prims to msc_reshaped_prims
  do iprim=1,nmsc_to_use
!print*,'old msc_reshaped_prims:',iprim,fix_prims(iprim),msc_reshaped_prims(fix_prims(iprim),:)
    msc_reshaped_prims(fix_prims(iprim),:) = new_prims(iprim,:)
!print*,'new msc_reshaped_prims:',iprim,fix_prims(iprim),msc_reshaped_prims(fix_prims(iprim),:)
    msc_tern_reshaped(fix_prims(iprim),1:nint(new_prims(iprim,8)),:) = new_terns(iprim,1:nint(new_prims(iprim,8)),:)
  enddo
!  call status_update('STATUS','UTIL','voronoi_reshape_v3',' ',message,0)
  write(message,*)"Modified number of ternarys:",nint(new_prims(1:nmsc_to_use,8))
  call status_update('STATUS','UTIL','voronoi_reshape_v3',' ',message,0)

  return
  end subroutine voronoi_reshape_v3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine jiggle_mascons(region,rms_limit,cross_0E,ntern_to_reshape,nprim_to_make &
                           ,msc_reshaped_prims,msc_tern_reshaped)

! subroutine to identify mascons with too many/too few ternarys, select a subset of primary mascons that
! surround the offending mascon and then reshape them all. This has the effect of averaging the number 
! of ternarys in each of the mascons in the subset, thereby bringing the overall pattern towards each
! primary having the same number of ternary mascons.
!
! P. Tregoning
! 13 February 2020

  use mascon_mod

  implicit none

! passed variables
  character*15,intent(in)    :: region                                  ! name of region being reshaped
  real(kind=8),intent(in)    :: rms_limit                               ! convergence criterion for voronoi_reshape
  logical     ,intent(in)    :: cross_0E                                ! does the region cross the Greenwich meridian?
  integer*4,   intent(in)    :: ntern_to_reshape                        ! total number of ternarys to be reassigned into new primaries
  integer*4,   intent(in)    :: nprim_to_make                           ! total number of new primary mascons to be made
  !real(kind=8),intent(inout) :: tern_crds(ntern_to_reshape,nvar_tern)   ! the complete list of ternarys to be reassigned
  real(kind=8),intent(inout) :: msc_tern_reshaped(nprim_to_make,max_tern_per_prim,nvar_tern) ! the new allocation of ternarys to primarys. Mimics mcon_tern array.
  real(kind=8),intent(inout) :: msc_reshaped_prims(nprim_to_make,nvar_prim)  ! information on reshaped primary mascons. Contains the CoM info to begin with.

! local variables
  integer*4             :: nloops            ! the number of times the loop has been executed to fix a subset of mascons
  integer*4,parameter   :: max_loops = 500   ! the maximum number of loops
  integer*4             :: iprim,iprim2,itern
  character*200         :: message
  integer*4             :: fix_prim          ! the primary mascon identified to be fixed
  real(kind=8),allocatable :: nterns(:,:)       ! array of pointers and numbers of ternarys per new primary mascons
  real(kind=8)          :: pi

! variables for working out which nearby primaries to include
  real(kind=8),allocatable :: fix_dist(:,:)
  integer*4                :: nmsc_to_use
  integer*4                :: fix_prims(1000)  ! we should never want to use as many as 1000 neighbouring mascons!
  real(kind=8)             :: near_mscs        ! the maximum distance away a neighbouring mascon can be and still be used in the subset
  logical                  :: debug, debug_rms
  integer*4                :: sum_terns        ! the sum of all the ternarys in the subset of primary mascons to be reshaped
  real(kind=8),allocatable :: tern_crds(:,:)
  real*8                   :: rand_number

  allocate(nterns(nprim_to_make,4))
  allocate(fix_dist(nprim_to_make,4))
  allocate(tern_crds(max_tern,nvar_tern))

  pi = 4.d0*datan(1.d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first, assign some metadata for the primary mascons passed in !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call status_update('STATUS','UTIL','jiggle_mascons',' ',"Assign metadata for new primary mascons",0)
!!$OMP PARALLEL DO private (iprim) &
!!$OMP             ,shared(nprim_to_make,ntern_to_reshape,msc_tern_reshaped,msc_reshaped_prims)
  do iprim = 1,nprim_to_make

    !call MASCON_percent_land(.false.,nvar_tern,nint(msc_reshaped_prims(iprim,8)),max_tern_per_prim &
    !            ,msc_tern_reshaped(iprim,:,:),msc_reshaped_prims(iprim,11))

    ! assign a density, based on the percentage of land
    if (msc_reshaped_prims(iprim,11) >= 0.5)then
      msc_reshaped_prims(iprim,6) = 1000.d0
    else
      msc_reshaped_prims(iprim,6) = 1029.d0
    endif
    ! identify how many secondary mascons there are in this primary. Hardwire it to 1
    msc_reshaped_prims(iprim,7) = 1.d0
    ! set the tidal amplitude to zero
    msc_reshaped_prims(iprim,10) = 0.d0

    ! make an array of primary numbers
    nterns(iprim,1) = iprim
  enddo
!!$OMP END PARALLEL DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the main loop starts here
  call status_update('STATUS','UTIL','jiggle_mascons',' ',"Starting main loop of jiggle_mascons",0)
  nloops = 1

  do while (nloops < max_loops)

  ! first, make a list of the number of ternarys per primary, in ascending order
    nterns(:,2) = msc_reshaped_prims(:,8)
!print*,'calling bubble_sort nterns(:,2)',nint(nterns(:,2))
    call bubble_sort_R8(nprim_to_make,nterns(:,1),nterns(:,2),nterns(:,3),nterns(:,4)) 
    ! select a primary mascon to repair
    if(mod(dble(nloops),10.d0) == 0.d0)then
      ! pick a random primary mascon to work on
      call srand(nprim_to_make)
      call random_number(rand_number) 
      fix_prim = int(rand_number*nprim_to_make+0.5)
      near_mscs = 1000.d3    ! set this to 1000 km. May need to adjust as we see how the algorithm performs.

      write(message,'(i5,a,i4,a,i5,a,i6,a)')nloops,'.  Random primary ',fix_prim,' has     ' &
                        ,nint(msc_reshaped_prims(fix_prim,8)), " ternarys (nominal is",nint(dble(ntern_to_reshape)/dble(nprim_to_make)),")"
    else if(mod(dble(nloops),2.d0) == 0.d0)then
      ! for even loops, choose to fix the primary with the most ternarys (if it is over the maximum allowed)
      fix_prim = nint(nterns(nprim_to_make,3))
      if(nloops < 20)then
        near_mscs = 1000.d4    ! set this to 1000 km. May need to adjust as we see how the algorithm performs.
      else
        near_mscs = 1000.d3    ! set this to 1000 km. May need to adjust as we see how the algorithm performs.
      endif

      write(message,'(i5,a,i4,a,i5,a,i6,a)')nloops,'.  Primary ',fix_prim,' has     ' &
                        ,nint(msc_reshaped_prims(fix_prim,8)), " ternarys (nominal is",nint(dble(ntern_to_reshape)/dble(nprim_to_make)),")"
    else if (mod(dble(nloops),2.d0) > 0.d0)then
    ! for odd loops, choose to fix the primary with the fewest ternarys (if it is under the maximum allowed)
      fix_prim = nint(nterns(1,3))
      near_mscs = 500.d3    ! set this to 1000 km. May need to adjust as we see how the algorithm performs.
      write(message,'(i5,a,i4,a,i5,a,i6,a)')nloops,'.  Primary ',fix_prim,' has only' &
                        ,nint(msc_reshaped_prims(fix_prim,8)), " ternarys (nominal is",nint(dble(ntern_to_reshape)/dble(nprim_to_make)),")"
    endif
    call status_update('STATUS','UTIL','assign_ternarys',' ',message,0)


    ! identify the other primary mascons within XXX km distance of the selected primary. For this, we need to calculate all the distances 
    ! between the CoMs and then bubble sort them
    do iprim2 = 1,nprim_to_make
      fix_dist(iprim2,1) = iprim2
      call MASCON_msc_to_msc(msc_reshaped_prims(fix_prim,1)*180.d0/pi,msc_reshaped_prims(fix_prim,2)*180.d0/pi &
                            ,msc_reshaped_prims(iprim2,1)*180.d0/pi,msc_reshaped_prims(iprim2,2)*180.d0/pi, fix_dist(iprim2,2) )
    enddo
    call bubble_sort_R8(nprim_to_make,fix_dist(:,1),fix_dist(:,2),fix_dist(:,3),fix_dist(:,4)) 
    
    ! carve out the ternarys that lie within this subset of primary mascons. We will use primary mascons within near_mscs km
    nmsc_to_use = 1
    fix_prims(1) = fix_prim
    do iprim=2,nprim_to_make   ! start at two because the closest mascon will be itself at zero distance
      if(fix_dist(iprim,4) <= near_mscs .and. nmsc_to_use < 7)then
         nmsc_to_use = nmsc_to_use + 1
         fix_prims(nmsc_to_use) = fix_dist(iprim,3)
      endif
    enddo
!    print*,nmsc_to_use,"nearest mascons (to msc ",fix_prims(1),") :",fix_prims(2:nmsc_to_use)
    
    ! make an array of all the ternarys in the subset of mascons
    tern_crds = 0.d0
    sum_terns = 0
    do iprim=1,nmsc_to_use
      tern_crds(1+sum_terns:sum_terns+nint(msc_reshaped_prims(fix_prims(iprim),8)),:) =  &
                                            msc_tern_reshaped(fix_prims(iprim),1:nint(msc_reshaped_prims(fix_prims(iprim),8)),:)
      sum_terns = sum_terns+msc_reshaped_prims(fix_prims(iprim),8)
    enddo

    ! now call a subroutine that will:
    ! 1. select "nmsc_to_use" random ternarys to act as the a priori CoM values
    ! 2. run a voronoi cell algorithm to reshape this subset of primary mascons into equal area shapes
    ! 3  compute the new centres of mass, percent land etc and apply all the attributes 
    ! 3. return the newly configured primary mascon information in msc_reshaped_prims
    call voronoi_reshape_v3(nprim_to_make,ntern_to_reshape,nmsc_to_use,fix_prims,sum_terns,tern_crds &
                           ,msc_reshaped_prims,msc_tern_reshaped) 

    ! that's it - from the single run through this loop we assume that each mascon in the subset of primary mascons
    ! now has the same number of ternarys. That is, this process will have averaged the ternarys across the subset 
    ! of primary mascons.

    ! increment the loop number and do it again.
    nloops = nloops + 1
  enddo


  return

  end subroutine jiggle_mascons
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine assign_ternarys(ntern_to_reshape,nprim_to_make,tern_crds,msc_reshaped_prims,msc_tern_reshaped)

! assign each of a list of ternarys to the closest primary CoM
!
! P. Tregoning
! 12 February 2020

  use mascon_mod

  implicit none

! passed variables
  integer*4,   intent(in)    :: ntern_to_reshape                        ! total number of ternarys to be reassigned into new primaries
  integer*4,   intent(in)    :: nprim_to_make                           ! total number of new primary mascons to be made
  real(kind=8),intent(inout) :: tern_crds(ntern_to_reshape,nvar_tern)   ! the complete list of ternarys to be reassigned
  real(kind=8),intent(inout) :: msc_tern_reshaped(nprim_to_make,max_tern_per_prim,nvar_tern) ! the new allocation of ternarys to primarys. Mimics mcon_tern array.
  real(kind=8),intent(inout) :: msc_reshaped_prims(nprim_to_make,nvar_prim)  ! information on reshaped primary mascons. Contains the CoM info to begin with.

! local variables
  integer*4    :: iprim,itern                ! do loop counters
  real(kind=8) :: pi
  real(kind=8),allocatable :: tern_to_CoM(:,:,:)
  integer*4,allocatable    :: nterns(:)
  integer*4                :: closest_prim
  real(kind=8)             :: shortest_dist
  logical                  :: debug

  debug = .false.

  pi = 4.d0*datan(1.d0)
  allocate(tern_to_CoM(ntern_to_reshape,nprim_to_make,4))   ! ternary number, primary mascon number, dist to primary, 
  allocate(nterns(nprim_to_make))

! calculate the distance from each ternary to each primary mascon centre of mass
  call status_update('STATUS','UTIL','assign_ternarys',' ',"Calculating distances between ternarys and CoMs",0)

!$OMP PARALLEL DO private (iprim,itern) &
!$OMP            ,shared(tern_crds,pi,msc_reshaped_prims ,nprim_to_make,ntern_to_reshape,tern_to_CoM)
  do itern = 1,ntern_to_reshape
    do iprim = 1, nprim_to_make
        tern_to_CoM(itern,iprim,1) = iprim     ! assign the primary mascon number to the ternary/primary distance to be calculated
        ! calculates the distance from this ternary to this new primary CoM
        call MASCON_msc_to_msc(tern_crds(itern,1)*180.d0/pi,tern_crds(itern,2)*180.d0/pi &
           ,msc_reshaped_prims(iprim,1)*180.d0/pi,msc_reshaped_prims(iprim,2)*180.d0/pi,tern_to_CoM(itern,iprim,2) )
      enddo
    enddo
!$OMP END PARALLEL DO

! find the closest CoM for each ternary
  call status_update('STATUS','UTIL','reshape_mascons',' ',"Finding the closest CoM for each ternary",0)
  nterns = 0
  do itern = 1,ntern_to_reshape
    shortest_dist = 1.d10
    closest_prim = 0
    do iprim = 1,nprim_to_make
!print*,ntern_to_reshape,nprim_to_make,itern,iprim
      if(tern_to_CoM(itern,iprim,2) < shortest_dist)then
        closest_prim = nint(tern_to_CoM(itern,iprim,1))
        shortest_dist = tern_to_CoM(itern,iprim,2)
if(debug)then
  print*,'new closest prim and dist = ',itern,iprim,closest_prim,shortest_dist
endif

      endif
    enddo

    ! add this ternary to the relevant primary mascon information
    tern_crds(itern,9) = nint(tern_to_CoM(itern,closest_prim,1))
    nterns(closest_prim) = nterns(closest_prim) + 1
    msc_tern_reshaped(nint(tern_to_CoM(itern,closest_prim,1)),nterns(closest_prim),1:8) = tern_crds(itern,1:8)
    msc_reshaped_prims(nint(tern_to_CoM(itern,closest_prim,1)),8) = nterns(closest_prim)   ! save the number of ternarys assigned to this new primary

  enddo

! DEBUG
debug = .false.
if(debug)then
  do iprim=1,nprim_to_make
    print*,'primary and number of ternarys:',iprim,nterns(iprim)
  enddo
  print*,'total assigned ternarys = ',sum(nterns(:))
endif


  return
  end subroutine assign_ternarys
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine make_tern_array(nprim_to_reshape,ntern_to_reshape,ocean_mascons,region_prims,cross_0E,tern_crds)

! given a subset of primary mascons, make a single array containing all the ternary information from the 
! sum of the primary mascons
!
! P. Tregoning
! 12 February 2020

  use mascon_mod

  implicit none

! passed variables
  integer*4, intent(in)      :: nprim_to_reshape                        ! number of primary mascons from which to take the ternary mascons
  integer*4, intent(in)      :: ntern_to_reshape                        ! total number of ternary mascons to be reshaped into new primaries
  logical  , intent(in)      :: ocean_mascons                           ! flag as to whether we are handling an ocean mascon file or not
  integer*4, intent(in)      :: region_prims(nprim_to_reshape)          ! vector of primary mascon numbers to be reshaped
  logical,      intent(in)   :: cross_0E                                ! true if region polygon and mascons span the 0East longitude (causes problems with 359E and 1E etc)
  real(kind=8),intent(inout) :: tern_crds(ntern_to_reshape,nvar_tern)   ! array of all the available ternary mascons (including all their attributes

! local variables
  integer*4    :: iprim,itern,counter,do_counter
  logical      :: debug
  real(kind=8) :: pi

  pi = 4.d0*datan(1.d0)
  debug = .false.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! transfer all the ternary mascons into one array
  counter = 0
  do iprim=1,nprim_to_reshape
    if(debug)print*,'incorporating ',nint(mcon_prim(region_prims(iprim),8)) &
                   ,' ternarys from original primary mascon:',region_prims(iprim)
! PT190523: we need a different do-loop counter here, depending on whether we are doing an ocean_mascon file or a normal_mascon file
    if(ocean_mascons)then
      do_counter = ntern_to_reshape
    else
      do_counter = nint(mcon_prim(region_prims(iprim),8))
    endif
    do itern=1,do_counter

      counter = counter+1
      ! PT190523: distinguish between ternarys from an ocean mascon file or from a normal mascon file
      if(ocean_mascons)then
        tern_crds(counter,1:8) = mcon_ocean_tern(itern,1:8)
      else
        tern_crds(counter,1:8) = mcon_tern(region_prims(iprim),itern,1:8)
      endif
      if(debug .and. itern < 10)print*,'tern_crds:',iprim,itern,tern_crds(counter,1:2)*180./pi,tern_crds(counter,3:7)

    enddo
  enddo

  return
  end subroutine make_tern_array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine voronoi_reshape_v2(ocean_mascons,region,rms_limit,cross_0E,n_to_reshape,region_prims,msc_reshaped_prims,n_reshaped &
                             ,msc_tern_reshaped,ntern_reshape,debug,debug_rms)

! subroutine to reshape the primary mascons within a particular region to make them roughly equal area and equal shape.
! The philosophy of the algorithm is based on voronoi cells and placing each mascon in the cell closest to it. The 
! centroids of the cells are then recalculated, and the process repeated until the exchange of ternary mascons from
! one primary to another is completed.
!
! Adapted from the program voronoi_mascons.f90 and turned into this subroutine.
!
! P. Tregoning
! 14 November 2016
!
! PT190408: modified to limit the area of primary mascons, so that the algorithm doesn't create primary mascons that are too small!
! PT200212: modified to loop through the new primarys and fill them, rather than loop through the ternarys and assign to the closest
!           primary. We have the listing of the 1st, 2nd, 3rd-closest etc primaries for each ternary. Therefore, if a new primary 
!           doesn't have enough ternarys in it from 1st choice, start using 2nd-closest etc until the primary has enough ternarys.


  use mascon_mod 
  use omp_lib

  implicit none

! passed variables
  logical     , intent(in)   :: ocean_mascons                      ! indicate whether to use ternarys from the ocean mascon file or the normal file (T/F)
  character*15, intent(in)   :: region
  real(kind=8), intent(in)   :: rms_limit                          ! lower value of rms of centroid adjustments, below which we stop the adjustment algorithm
  logical,      intent(in)   :: cross_0E                           ! true if region polygon and mascons span the 0East longitude (causes problems with 359E and 1E etc)
  integer*4,    intent(in)   :: n_to_reshape                       ! number of input primary mascons to be reshaped
  integer*4,    intent(in)   :: n_reshaped                         ! required number of primary mascons of the required size for the region provided
  integer*4,    intent(in)   :: ntern_reshape                      ! total number of ternarys in the "n_to_reshape" primary mascons
  integer*4,    intent(in)   :: region_prims(n_to_reshape)         ! vector of numbers of primary mascons to be reshaped
  real(kind=8),intent(inout) :: msc_reshaped_prims(n_reshaped,nvar_prim)                         ! reshaped primary mascons (lat/lon/rad/area/depth/density/#sec/#tern/#1st sec/bit-mapped tide/%land
  real(kind=8),intent(inout) :: msc_tern_reshaped(n_reshaped, nint(ntern_reshape/1.0),nvar_tern) ! the ternary mascons to be reshaped (not sure about this for the moment ...)
  logical     ,intent(inout) :: debug,debug_rms

! variables related to deciding which is the closest primary mascon
  real(kind=8)  :: min_dist,tern_dist
  integer*4     :: closest_node
  real(kind=8),allocatable :: tern_crds(:,:)         ! array of all the available ternary mascons (including all their attributes
  integer*4   ,allocatable :: n_terns(:)             ! the number of ternary mascons in each reshaped primary mascon
  logical,     allocatable :: rand_numbers(:)        ! array of flags to say whether a particular random number has been used or not

! variables to compute the RMS of the primary mascon adjustments
  real(kind=8),allocatable :: prim_crd_old(:,:)
  real(kind=8) :: RMS,dist,dlon
  integer*4    :: niter
  real(kind=8),parameter :: earthrad = 6378136.46d0

! local variables
  integer*4    :: extra_prim                         ! temporary random ternary number to act as seed coords for extra primary mascons
  integer*4    :: random_tern                        ! a random ternary number
  integer*4    :: itern,iprim,itern2,k,counter                      ! counters
  character*200 :: message
  real(kind=8) :: pi
  real(kind=8) :: rand_number
  character*150:: message2
  real(kind=8) :: tmp_calc                           ! temporary variable in calculating spherical distance between two mascons

! PT190523: variables specific for ocean_mascon file analysis
  integer*4    :: do_counter

! PT190408: variables related to limiting the area of primary mascons
  integer*4    :: permitted_terns_per_prim                   ! limit the number of permitted ternarys in a new primary
  real(kind=8),allocatable :: tern_to_CoM(:,:,:)             ! array used to order tern-CoM distances in increasing size
  logical      :: found_closest
  integer*4    :: itmp
  real(kind=8),allocatable :: tern_index(:,:)                ! array used to order ternary mascons by descending latitude
  integer*4    :: start_counter,end_counter,d_counter
  integer*4    :: nclosest                                   ! count of the number of closest ternarys to include in a new primary
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! declare and allocate certain variables
  pi = 4.d0*datan(1.d0)

  allocate(tern_crds(ntern_reshape,nvar_tern))
  allocate(n_terns(n_reshaped))
  allocate(prim_crd_old(n_reshaped,2))
  allocate(rand_numbers(ntern_reshape))
! PT190408: allocate the array of distances of all ternarys to all CoMs
  allocate(tern_to_CoM(ntern_reshape,n_reshaped,4))   ! ternary number, primary mascon number, dist to primary, 
  allocate(tern_index(ntern_reshape,4))

  n_terns = 0
! PT190408: permitted ternarys per primary. Make it 110%
  permitted_terns_per_prim = nint(1.05d0*dble(ntern_reshape)/dble(n_reshaped))
! debug
  write(message,'(a,i7,a,i10,a)')'Will limit primary mascons to a maximum of ',permitted_terns_per_prim &
                          ,' ternarys (there are',ntern_reshape,' ternarys in total)'
  call status_update('STATUS','UTIL','voronoi_reshape',' ',message,0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT161118: shift all longitudes 90 deg east if the mascons span
!           the 0 deg longitude line
  if(cross_0E)then
    write(message,'(a)')"Shifting longitudes by 90 deg for region "
    if(debug)call status_update('STATUS','UTIL','voronoi_reshape',' ',message,0)
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! transfer all the ternary mascons into one array
  counter = 0
  do iprim=1,n_to_reshape

! PT190523: we need a different do-loop counter here, depending on whether we are doing an ocean_mascon file or a normal_mascon file
    if(ocean_mascons)then
      do_counter = ntern_reshape
    else
      do_counter = nint(mcon_prim(region_prims(iprim),8))
    endif
    do itern=1,do_counter
      counter = counter+1
! PT190523: distinguish between ternarys from an ocean mascon file or from a normal mascon file
      if(ocean_mascons)then
        tern_crds(counter,1:8) = mcon_ocean_tern(itern,1:8)
      else
        tern_crds(counter,1:8) = mcon_tern(region_prims(iprim),itern,1:8)
      endif
      if(debug)print*,'tern_crds:',iprim,itern,tern_crds(counter,1:2)*180./pi,tern_crds(counter,3:7)
! PT161118: shift all longitudes 90 deg east if the mascons span the 0 deg longitude line
      if(cross_0E)then
        tern_crds(counter,2) = tern_crds(counter,2) + pi/2.d0
        if(tern_crds(counter,2) > 2.d0*pi)tern_crds(counter,2) = tern_crds(counter,2) - 2.d0*pi
      endif
    enddo
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! set up our array of starting primary mascons
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do iprim=1,n_to_reshape
! Pt161117: only do this if the mascon number is less than or equal to the total number we want.
    if(iprim <= n_reshaped)then
      msc_reshaped_prims(iprim,1:nvar_prim) = mcon_prim(region_prims(iprim),1:nvar_prim)    ! lat/lon/rad/area/depth/density/.... what are the others?
! PT161118: shift all longitudes 90 deg east if the mascons span the 0 deg longitude line
      if(cross_0E)then
        msc_reshaped_prims(iprim,2) = msc_reshaped_prims(iprim,2) + pi/2.d0
        if(msc_reshaped_prims(iprim,2) > 2.d0*pi)msc_reshaped_prims(iprim,2) = msc_reshaped_prims(iprim,2)-2.d0*pi
      endif
    endif
  enddo

! DEBUG
if(debug)print*,'msc_reshaped_prims(10):',msc_reshaped_prims(10,1:2)*180./pi,msc_reshaped_prims(10,8)

! if we need more primary mascons than we found in the region, choose some random ternary mascons and use them as seed points
  if(n_reshaped > n_to_reshape )then
    write(message,'(a,i6,a)')'Need',n_reshaped - n_to_reshape,' seed primary mascons'
    if(debug)call status_update('STATUS','UTIL','voronoi_reshape',' ',message,0)
    call srand(ntern_reshape)
    rand_numbers = .false.
    do iprim=1,n_reshaped - n_to_reshape + 1   ! PT161125: added 1 so that entries with only one point are all randomly selected
      call get_random_value(rand_numbers,ntern_reshape,extra_prim)

! and assign the coord information etc
      msc_reshaped_prims(iprim+n_to_reshape-1,1:3) = tern_crds(extra_prim,1:3)
    enddo
  endif

! if there were too many primary mascons for the way we want to break up the region, do something different
  if(n_reshaped < n_to_reshape)then
    write(message,'(a,i6,a)')'Randomly seed',n_reshaped,'  primary mascons'
    call status_update('STATUS','UTIL','voronoi_reshape',' ',message,0)
!   randomly select n_reshaped ternary mascons from the entire list and use them as the starting coords
    call srand(ntern_reshape)
    do iprim=1,n_reshaped
      random_tern = nint(rand(0)*ntern_reshape)
      msc_reshaped_prims(iprim,1:3) = tern_crds(random_tern,1:3)
if(debug)print*,'random coords',iprim,msc_reshaped_prims(iprim,1:2)*180./pi
    enddo
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! make a copy of the starting coords for the primary mascons
  prim_crd_old(:,1:2) = msc_reshaped_prims(:,1:2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! assign ternarys to their closest primary mascon node

  niter = 0
  RMS = 1.d9
  do while (RMS > rms_limit)
    n_terns = 0
    tern_to_CoM = 0.d0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PT190408: calculate and store the distance from each ternary to each !!
  !           current CoM

    do iprim=1,n_reshaped
      tern_to_CoM(:,iprim,1) = iprim
    enddo

!$OMP PARALLEL DO private (iprim,itern) &
!$OMP            ,shared(tern_crds,pi,msc_reshaped_prims ,n_reshaped,ntern_reshape,tern_to_CoM)
    do itern=1,ntern_reshape     ! ntern_reshape is the total number of ternarys to be placed into new primary mascon shapes
      do iprim=1,n_reshaped      ! n_reshaped is the total number of new primary mascon shapes
! DEBUG
if(debug)print*,'tern crds',itern,tern_crds(itern,1)*180.d0/pi,tern_crds(itern,2)*180.d0/pi
        ! calculates the distance from this ternary to this new primary CoM
        call MASCON_msc_to_msc(tern_crds(itern,1)*180.d0/pi,tern_crds(itern,2)*180.d0/pi &
           ,msc_reshaped_prims(iprim,1)*180.d0/pi,msc_reshaped_prims(iprim,2)*180.d0/pi,tern_to_CoM(itern,iprim,2) )
      enddo
    enddo
!$OMP END PARALLEL DO
    write(message,'(a,i6)')"have calculated ternarys to primary CoMs for iteration",int(niter)
    if(debug)call status_update('STATUS','UTIL','voronoi_reshape',' ',message,0)

  ! PT190408: need to make a list of the order of primary mascons, in increasing distance
!$OMP PARALLEL DO private (itern) &
!$OMP            ,shared(ntern_reshape,tern_to_CoM)
    do itern=1,ntern_reshape
      call bubble_sort_R8(n_reshaped,tern_to_CoM(itern,:,1),tern_to_CoM(itern,:,2),tern_to_CoM(itern,:,3),tern_to_CoM(itern,:,4))
    enddo
!$OMP END PARALLEL DO
    write(message,'(a,i6)')"have done the bubble sort for iteration",int(niter)
    if(debug)call status_update('STATUS','UTIL','voronoi_reshape',' ',message,0)
    
  ! PT190408: we now have the distance from each ternary to each CoM. Now need to loop through the ternarys and assign
  !           to primarys, ensuring that we don't overfill any primary. 

  ! PT190408: can I loop through the ternarys in increasing latitude order? I think this might help convergence ....
    do itern=1,ntern_reshape
      tern_index(itern,1) = itern
    enddo

    start_counter = 1
    end_counter = ntern_reshape
    d_counter = 1
    tern_index(:,3) = tern_index(:,1)
    permitted_terns_per_prim = nint(1.d0*dble(ntern_reshape)/dble(n_reshaped))


    do itern=1,ntern_reshape
      ! which primary was closest
      found_closest = .false.
      itern2 = tern_index(itern,3)

      do iprim=1,n_reshaped
        if(.not. found_closest .and. n_terns(nint(tern_to_CoM(itern2,iprim,3))) < permitted_terns_per_prim)then   ! we can add this ternary to this primary
          found_closest = .true.
          n_terns(nint(tern_to_CoM(itern2,iprim,3))) = n_terns(nint(tern_to_CoM(itern2,iprim,3))) + 1
          tern_crds(itern,9) = nint(tern_to_CoM(itern2,iprim,3))
          do k=1,8
            msc_tern_reshaped(nint(tern_to_CoM(itern2,iprim,3)),n_terns(nint(tern_to_CoM(itern2,iprim,3))),k) = tern_crds(itern2,k)
          enddo
        endif
      enddo
    enddo
    write(message,'(a,i6)')"have found the closest CoM for iteration",int(niter)
    if(debug)call status_update('STATUS','UTIL','voronoi_reshape',' ',message,0)
       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now calculate the centroids of the modified set of primary mascons and the percent land
    RMS  = 0.d0

!$OMP PARALLEL DO private (iprim,dlon,dist,tmp_calc) &
!$OMP             ,shared(n_reshaped,ntern_reshape,msc_tern_reshaped,msc_reshaped_prims,prim_crd_old) &
!$OMP             ,reduction(+:RMS)
    do iprim = 1,n_reshaped
      msc_reshaped_prims(iprim,8) = n_terns(iprim)


      if(msc_reshaped_prims(iprim,8) > 0)then
        call MASCON_calc_CoM(ntern_reshape ,nvar_prim,nvar_tern,msc_tern_reshaped(iprim,:,:),msc_reshaped_prims(iprim,:) )
        call MASCON_percent_land(.false.,nvar_tern,nint(msc_reshaped_prims(iprim,8)),ntern_reshape &
                ,msc_tern_reshaped(iprim,:,:),msc_reshaped_prims(iprim,11))
      endif

! reset the density based on the percent land
      if(msc_reshaped_prims(iprim,11) > 0.5)then
        msc_reshaped_prims(iprim,6) = 1000.d0
      else if(msc_reshaped_prims(iprim,11) >= 0)then
        msc_reshaped_prims(iprim,6) = 1029.d0
      else ! RM190522: set GIA mascons
        msc_reshaped_prims(iprim,6) = 3300.d0
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! and now the rms of the coordinate movements. Do this by calculating the great circle distance
! between the old and new lat/lon, multiplied by radius ofthe Earth.
      dlon = msc_reshaped_prims(iprim,2) - prim_crd_old(iprim,2)
! PT170608: stop it wrapping right around the Earth for longitude differences > 180 degrees
      if(dlon > pi)dlon = 2.d0*pi - dlon
      if(dlon == 0.d0 .and. dabs(prim_crd_old(iprim,1)-msc_reshaped_prims(iprim,1)) < 1.d-7)then
        dist = 0.d0
      else
        tmp_calc = dcos(prim_crd_old(iprim,1))*dcos(msc_reshaped_prims(iprim,1))  &
             +dsin(prim_crd_old(iprim,1))*dsin(msc_reshaped_prims(iprim,1))*dcos(dlon)
! PT170608: try to mitigate NaNs coming from rounding errors pushing this slightly over 1.0
        if(tmp_calc > 1.d0)tmp_calc = 1.d0
        dist = dacos( tmp_calc )
      endif
      dist = dist*earthrad
      RMS = RMS + dist**2
    enddo
!$OMP END PARALLEL DO

    RMS = dsqrt(RMS/dble(n_reshaped))
    if(debug_rms)then
!       print*,"Iteration",niter," RMS of coord adjustments: ",RMS," Cutoff RMS level:",rms_limit &
!              ," Region: ",region
       write(message,'(a,i5,a,f15.3,a,f15.3,a,a)')"Iteration",int(niter)," RMS of coord adjustments: ",dble(RMS) &
              ," Cutoff RMS level:",dble(rms_limit)," Region: ",region
       call status_update('STATUS','UTIL','voronoi_reshape',' ',message,0)
    endif

    if(RMS > rms_limit)then
      niter = niter + 1
      prim_crd_old(:,1:2) = msc_reshaped_prims(:,1:2)
    endif
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(cross_0E)then
    write(message,'(a)')"Shifting longitudes back by 90 deg for region "
    if(debug)call status_update('STATUS','UTIL','voronoi_reshape',' ',message,0)

! first the primary mascons
    do iprim=1,n_reshaped
      msc_reshaped_prims(iprim,2) = msc_reshaped_prims(iprim,2) - pi/2.d0
      if(msc_reshaped_prims(iprim,2) < 0.d0)msc_reshaped_prims(iprim,2) = msc_reshaped_prims(iprim,2) + 2.d0*pi
    enddo

! now the ternary mascons
    do iprim=1,n_reshaped
      msc_tern_reshaped(iprim,:,2) = msc_tern_reshaped(iprim,:,2) - pi/2.d0
      do itern=1,nint(msc_reshaped_prims(iprim,8))
        if(msc_tern_reshaped(iprim,itern,2) < 0.d0)msc_tern_reshaped(iprim,itern,2) = msc_tern_reshaped(iprim,itern,2) +2.d0*pi
      enddo
    enddo

  endif


  return
  end subroutine voronoi_reshape_v2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine voronoi_reshape(ocean_mascons,region,rms_limit,cross_0E,n_to_reshape,region_prims,msc_reshaped_prims,n_reshaped &
                             ,msc_tern_reshaped,ntern_reshape,debug,debug_rms)

! subroutine to reshape the primary mascons within a particular region to make them roughly equal area and equal shape.
! The philosophy of the algorithm is based on voronoi cells and placing each mascon in the cell closest to it. The 
! centroids of the cells are then recalculated, and the process repeated until the exchange of ternary mascons from
! one primary to another is completed.
!
! Adapted from the program voronoi_mascons.f90 and turned into this subroutine.
!
! P. Tregoning
! 14 November 2016
!
! PT190408: modified to limit the area of primary mascons, so that the algorithm doesn't create primary mascons that are too small!

  use mascon_mod 
  use omp_lib

  implicit none

! passed variables
  logical     , intent(in)   :: ocean_mascons                      ! indicate whether to use ternarys from the ocean mascon file or the normal file (T/F)
  character*15, intent(in)   :: region
  real(kind=8), intent(in)   :: rms_limit                          ! lower value of rms of centroid adjustments, below which we stop the adjustment algorithm
  logical,      intent(in)   :: cross_0E                           ! true if region polygon and mascons span the 0East longitude (causes problems with 359E and 1E etc)
  integer*4,    intent(in)   :: n_to_reshape                       ! number of input primary mascons to be reshaped
  integer*4,    intent(in)   :: n_reshaped                         ! required number of primary mascons of the required size for the region provided
  integer*4,    intent(in)   :: ntern_reshape                      ! total number of ternarys in the "n_to_reshape" primary mascons
  integer*4,    intent(in)   :: region_prims(n_to_reshape)         ! vector of numbers of primary mascons to be reshaped
  real(kind=8),intent(inout) :: msc_reshaped_prims(n_reshaped,nvar_prim)                         ! reshaped primary mascons (lat/lon/rad/area/depth/density/#sec/#tern/#1st sec/bit-mapped tide/%land
  real(kind=8),intent(inout) :: msc_tern_reshaped(n_reshaped, nint(ntern_reshape/1.0),nvar_tern) ! the ternary mascons to be reshaped (not sure about this for the moment ...)
  logical     ,intent(inout) :: debug,debug_rms

! variables related to deciding which is the closest primary mascon
  real(kind=8)  :: min_dist,tern_dist
  integer*4     :: closest_node
  real(kind=8),allocatable :: tern_crds(:,:)         ! array of all the available ternary mascons (including all their attributes
  integer*4   ,allocatable :: n_terns(:)             ! the number of ternary mascons in each reshaped primary mascon
  logical,     allocatable :: rand_numbers(:)        ! array of flags to say whether a particular random number has been used or not

! variables to compute the RMS of the primary mascon adjustments
  real(kind=8),allocatable :: prim_crd_old(:,:)
  real(kind=8) :: RMS,dist,dlon
  integer*4    :: niter
  real(kind=8),parameter :: earthrad = 6378136.46d0

! local variables
  integer*4    :: extra_prim                         ! temporary random ternary number to act as seed coords for extra primary mascons
  integer*4    :: random_tern                        ! a random ternary number
  integer*4    :: i,j,k,counter                      ! counters
  character*200 :: message
  real(kind=8) :: pi
  real(kind=8) :: rand_number
  character*150:: message2
  real(kind=8) :: tmp_calc                           ! temporary variable in calculating spherical distance between two mascons

! PT190523: variables specific for ocean_mascon file analysis
  integer*4    :: do_counter

! PT190408: variables related to limiting the area of primary mascons
  integer*4    :: permitted_terns_per_prim                   ! limit the number of permitted ternarys in a new primary
  real(kind=8),allocatable :: tern_to_CoM(:,:,:)             ! array used to order tern-CoM distances in increasing size
  logical      :: found_closest
  integer*4    :: itmp,itern
  real(kind=8),allocatable :: tern_index(:,:)                ! array used to order ternary mascons by descending latitude
  integer*4    :: start_counter,end_counter,d_counter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! declare and allocate certain variables
  pi = 4.d0*datan(1.d0)

  allocate(tern_crds(ntern_reshape,nvar_tern))
  allocate(n_terns(n_reshaped))
  allocate(prim_crd_old(n_reshaped,2))
  allocate(rand_numbers(ntern_reshape))
! PT190408: allocate the array of distances of all ternarys to all CoMs
  allocate(tern_to_CoM(ntern_reshape,n_reshaped,4))
  allocate(tern_index(ntern_reshape,4))

  n_terns = 0
! PT190408: permitted ternarys per primary. Make it 110%
  permitted_terns_per_prim = nint(1.05d0*dble(ntern_reshape)/dble(n_reshaped))
! debug
  write(message,'(a,i7,a,i10,a)')'Will limit primary mascons to a maximum of ',permitted_terns_per_prim &
                          ,' ternarys (there are',ntern_reshape,' ternarys in total)'
  call status_update('STATUS','UTIL','voronoi_reshape',' ',message,0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT161118: shift all longitudes 90 deg east if the mascons span
!           the 0 deg longitude line
  if(cross_0E)then
    write(message,'(a)')"Shifting longitudes by 90 deg for region "
    if(debug)call status_update('STATUS','UTIL','voronoi_reshape',' ',message,0)
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! transfer all the ternary mascons into one array
  counter = 0
  do i=1,n_to_reshape

! PT190523: we need a different do-loop counter here, depending on whether we are doing an ocean_mascon file or a normal_mascon file
      if(ocean_mascons)then
        do_counter = ntern_reshape
      else
        do_counter = nint(mcon_prim(region_prims(i),8))
      endif
      do j=1,do_counter
        counter = counter+1
! PT190523: distinguish between ternarys from an ocean mascon file or from a normal mascon file
        if(ocean_mascons)then
          tern_crds(counter,1:8) = mcon_ocean_tern(j,1:8)
        else
          tern_crds(counter,1:8) = mcon_tern(region_prims(i),j,1:8)
        endif
if(debug)print*,'tern_crds:',i,j,tern_crds(counter,1:2)*180./pi,tern_crds(counter,3:7)
! PT161118: shift all longitudes 90 deg east if the mascons span the 0 deg longitude line
        if(cross_0E)then
          tern_crds(counter,2) = tern_crds(counter,2) + pi/2.d0
          if(tern_crds(counter,2) > 2.d0*pi)tern_crds(counter,2) = tern_crds(counter,2) - 2.d0*pi
        endif
      enddo
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! set up our array of starting primary mascons
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do i=1,n_to_reshape
! Pt161117: only do this if the mascon number is less than or equal to the total number we want.
    if(i <= n_reshaped)then
      msc_reshaped_prims(i,1:nvar_prim) = mcon_prim(region_prims(i),1:nvar_prim)    ! lat/lon/rad/area/depth/density/.... what are the others?
! PT161118: shift all longitudes 90 deg east if the mascons span the 0 deg longitude line
      if(cross_0E)then
        msc_reshaped_prims(i,2) = msc_reshaped_prims(i,2) + pi/2.d0
        if(msc_reshaped_prims(i,2) > 2.d0*pi)msc_reshaped_prims(i,2) = msc_reshaped_prims(i,2)-2.d0*pi
      endif
    endif
  enddo

! DEBUG
if(debug)print*,'msc_reshaped_prims(10):',msc_reshaped_prims(10,1:2)*180./pi,msc_reshaped_prims(10,8)

! if we need more primary mascons than we found in the region, choose some random ternary mascons and use them as seed points
  if(n_reshaped > n_to_reshape )then
    write(message,'(a,i6,a)')'Need',n_reshaped - n_to_reshape,' seed primary mascons'
    if(debug)call status_update('STATUS','UTIL','voronoi_reshape',' ',message,0)
    call srand(ntern_reshape)
    rand_numbers = .false.
    do i=1,n_reshaped - n_to_reshape + 1   ! PT161125: added 1 so that entries with only one point are all randomly selected
      call get_random_value(rand_numbers,ntern_reshape,extra_prim)

! and assign the coord information etc
      msc_reshaped_prims(i+n_to_reshape-1,1:3) = tern_crds(extra_prim,1:3)
    enddo
  endif

! if there were too many primary mascons for the way we want to break up the region, do something different
  if(n_reshaped < n_to_reshape)then
    write(message,'(a,i6,a)')'Randomly seed',n_reshaped,'  primary mascons'
    call status_update('STATUS','UTIL','voronoi_reshape',' ',message,0)
!   randomly select n_reshaped ternary mascons from the entire list and use them as the starting coords
    call srand(ntern_reshape)
    do i=1,n_reshaped
      random_tern = nint(rand(0)*ntern_reshape)
      msc_reshaped_prims(i,1:3) = tern_crds(random_tern,1:3)
if(debug)print*,'random coords',i,msc_reshaped_prims(i,1:2)*180./pi
    enddo
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! make a copy of the starting coords for the primary mascons
  prim_crd_old(:,1:2) = msc_reshaped_prims(:,1:2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! assign ternarys to their closest primary mascon node

  niter = 0
  RMS = 1.d9
  do while (RMS > rms_limit)
    n_terns = 0
    tern_to_CoM = 0.d0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PT190408: calculate and store the distance from each ternary to each !!
  !           current CoM

    do i=1,n_reshaped
      tern_to_CoM(:,i,1) = i
    enddo

!$OMP PARALLEL DO private (i,j) &
!$OMP            ,shared(tern_crds,pi,msc_reshaped_prims ,n_reshaped,ntern_reshape,tern_to_CoM)
    do i=1,ntern_reshape     ! ntern_reshape is the total number of ternarys to be placed into new primary mascon shapes
      do j=1,n_reshaped      ! n_reshaped is the total number of new primary mascon shapes
! DEBUG
if(debug)print*,'tern crds',i,tern_crds(i,1)*180.d0/pi,tern_crds(i,2)*180.d0/pi
        ! calculates the distance from this ternary to this new primary CoM
        call MASCON_msc_to_msc(tern_crds(i,1)*180.d0/pi,tern_crds(i,2)*180.d0/pi &
           ,msc_reshaped_prims(j,1)*180.d0/pi,msc_reshaped_prims(j,2)*180.d0/pi,tern_to_CoM(i,j,2) )
      enddo
    enddo
!$OMP END PARALLEL DO
    write(message,'(a,i6)')"have calculated ternarys to primary CoMs for iteration",int(niter)
    if(debug)call status_update('STATUS','UTIL','voronoi_reshape',' ',message,0)

  ! PT190408: need to make a list of the order of primary mascons, in increasing distance
!$OMP PARALLEL DO private (i) &
!$OMP            ,shared(ntern_reshape,tern_to_CoM)
    do i=1,ntern_reshape
      call bubble_sort_R8(n_reshaped,tern_to_CoM(i,:,1),tern_to_CoM(i,:,2),tern_to_CoM(i,:,3),tern_to_CoM(i,:,4))
    enddo
!$OMP END PARALLEL DO
    write(message,'(a,i6)')"have done the bubble sort for iteration",int(niter)
    if(debug)call status_update('STATUS','UTIL','voronoi_reshape',' ',message,0)
    
  ! PT190408: we now have the distance from each ternary to each CoM. Now need to loop through the ternarys and assign
  !           to primarys, ensuring that we don't overfill any primary. 

  ! PT190408: can I loop through the ternarys in increasing latitude order? I think this might help convergence ....
    do i=1,ntern_reshape
      tern_index(i,1) = i
    enddo
    ! PT190408: for even iterations, sort on latitude. For odd iterations, sort on longitude
    if(region(1:6) == "Arctic" .or. region(1:5) == "SthnO" .or. region(1:4) == "Pine" .or. region(1:5) == "Asian" &
        .or. region(1:5) == "Green" .or. 1 == 1 )then
      start_counter = 1
      end_counter = ntern_reshape
      d_counter = 1
      tern_index(:,3) = tern_index(:,1)
      permitted_terns_per_prim = nint(10.1d0*dble(ntern_reshape)/dble(n_reshaped))
! PT200129: put this back to doing every case (by doing it if 1 = 1)
    elseif(region(1:9) == "Antarctic" ) then !.or. 1 == 1)then
      start_counter = 1
      end_counter = ntern_reshape
      d_counter = 1
      tern_index(:,3) = tern_index(:,1)
      permitted_terns_per_prim = nint(1.d0*dble(ntern_reshape)/dble(n_reshaped))
! PT190716: set these mod tests to -1.d0 from 0.d0, so that they will fail.
    elseif(mod(dble(niter),2.0) == -1.d0)then
      call bubble_sort_R8(ntern_reshape,tern_index(:,1),tern_crds(:,1),tern_index(:,3),tern_index(:,4))
      if(mod(dble(niter),4.0) == 0.d0)then
        start_counter = 1
        end_counter = ntern_reshape
        d_counter = 1
      else
        start_counter = ntern_reshape
        end_counter = 1
        d_counter = -1
      endif
    else
! PT200129: change this to sort on latitude (from longitude)
      call bubble_sort_R8(ntern_reshape,tern_index(:,1),tern_crds(:,1),tern_index(:,3),tern_index(:,4))
      if(mod(dble(niter+1),4.0) == -1.d0)then
        start_counter = 1
        end_counter = ntern_reshape
        d_counter = 1
      else
        start_counter = ntern_reshape
        end_counter = 1
        d_counter = -1
      tern_index(:,3) = tern_index(:,1)
      endif
    endif


print*,'permitted_terns_per_prim = ',permitted_terns_per_prim
print*,'set permitted_terns_per_prim to 252'
permitted_terns_per_prim = 252
    do i=start_counter,end_counter,d_counter
      ! which primary was closest
      found_closest = .false.
      itern = tern_index(i,3)
!      itern = i
      do j=1,n_reshaped
        if(.not. found_closest .and. n_terns(nint(tern_to_CoM(itern,j,3))) < permitted_terns_per_prim)then   ! we can add this ternary to this primary
          found_closest = .true.
          n_terns(nint(tern_to_CoM(itern,j,3))) = n_terns(nint(tern_to_CoM(itern,j,3))) + 1
          tern_crds(itern,9) = nint(tern_to_CoM(itern,j,3))
          do k=1,8
            msc_tern_reshaped(nint(tern_to_CoM(itern,j,3)),n_terns(nint(tern_to_CoM(itern,j,3))),k) = tern_crds(itern,k)
          enddo
        endif
      enddo
    enddo
 write(message,'(a,i6)')"have found the closest CoM for iteration",int(niter)
if(debug)call status_update('STATUS','UTIL','voronoi_reshape',' ',message,0)
       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now calculate the centroids of the modified set of primary mascons and the percent land
    RMS  = 0.d0

!$OMP PARALLEL DO private (i,dlon,dist,tmp_calc) &
!$OMP             ,shared(n_reshaped,ntern_reshape,msc_tern_reshaped,msc_reshaped_prims,prim_crd_old) &
!$OMP             ,reduction(+:RMS)
    do i = 1,n_reshaped
      msc_reshaped_prims(i,8) = n_terns(i)


      if(msc_reshaped_prims(i,8) > 0)then
        call MASCON_calc_CoM(ntern_reshape ,nvar_prim,nvar_tern,msc_tern_reshaped(i,:,:),msc_reshaped_prims(i,:) )
        call MASCON_percent_land(.false.,nvar_tern,nint(msc_reshaped_prims(i,8)),ntern_reshape &
                ,msc_tern_reshaped(i,:,:),msc_reshaped_prims(i,11))
      endif

! reset the density based on the percent land
      if(msc_reshaped_prims(i,11) > 0.5)then
        msc_reshaped_prims(i,6) = 1000.d0
      else if(msc_reshaped_prims(i,11) >= 0)then
        msc_reshaped_prims(i,6) = 1029.d0
      else ! RM190522: set GIA mascons
        msc_reshaped_prims(i,6) = 3300.d0
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! and now the rms of the coordinate movements. Do this by calculating the great circle distance
! between the old and new lat/lon, multiplied by radius ofthe Earth.
      dlon = msc_reshaped_prims(i,2) - prim_crd_old(i,2)
! PT170608: stop it wrapping right around the Earth for longitude differences > 180 degrees
      if(dlon > pi)dlon = 2.d0*pi - dlon
      if(dlon == 0.d0 .and. dabs(prim_crd_old(i,1)-msc_reshaped_prims(i,1)) < 1.d-7)then
        dist = 0.d0
      else
        tmp_calc = dcos(prim_crd_old(i,1))*dcos(msc_reshaped_prims(i,1))  &
             +dsin(prim_crd_old(i,1))*dsin(msc_reshaped_prims(i,1))*dcos(dlon)
! PT170608: try to mitigate NaNs coming from rounding errors pushing this slightly over 1.0
        if(tmp_calc > 1.d0)tmp_calc = 1.d0
        dist = dacos( tmp_calc )
      endif
      dist = dist*earthrad
      RMS = RMS + dist**2
    enddo
!$OMP END PARALLEL DO

    RMS = dsqrt(RMS/dble(n_reshaped))
    if(debug_rms)then
!       print*,"Iteration",niter," RMS of coord adjustments: ",RMS," Cutoff RMS level:",rms_limit &
!              ," Region: ",region
       write(message,'(a,i5,a,f15.3,a,f15.3,a,a)')"Iteration",int(niter)," RMS of coord adjustments: ",dble(RMS) &
              ," Cutoff RMS level:",dble(rms_limit)," Region: ",region
       call status_update('STATUS','UTIL','voronoi_reshape',' ',message,0)
    endif

    if(RMS > rms_limit)then
      niter = niter + 1
      prim_crd_old(:,1:2) = msc_reshaped_prims(:,1:2)
    endif
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(cross_0E)then
    write(message,'(a)')"Shifting longitudes back by 90 deg for region "
    if(debug)call status_update('STATUS','UTIL','voronoi_reshape',' ',message,0)

! first the primary mascons
    do i=1,n_reshaped
      msc_reshaped_prims(i,2) = msc_reshaped_prims(i,2) - pi/2.d0
      if(msc_reshaped_prims(i,2) < 0.d0)msc_reshaped_prims(i,2) = msc_reshaped_prims(i,2) + 2.d0*pi
    enddo

! now the ternary mascons
    do i=1,n_reshaped
      msc_tern_reshaped(i,:,2) = msc_tern_reshaped(i,:,2) - pi/2.d0
      do j=1,nint(msc_reshaped_prims(i,8))
        if(msc_tern_reshaped(i,j,2) < 0.d0)msc_tern_reshaped(i,j,2) = msc_tern_reshaped(i,j,2) +2.d0*pi
      enddo
    enddo

  endif


  return
  end subroutine voronoi_reshape

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine MASCON_msc_to_msc(lat1,lon1,lat2,lon2,dist)

  implicit none

! IN 
  real(kind=8) :: lat1,lon1,lat2,lon2          ! coordinates of two mascons

! OUT
  real(kind=8) :: dist

! LOCAL
  real(kind=8) :: xyz(2,3)
  real(kind=8), parameter :: earthrad = 6378136.45
  real(kind=8) :: pi, tmp_dist
  integer*4 :: j

  pi = 4.d0*atan(1.d0)

! convert lat/lon of mascon 1 into XYZ (assume spherical for now)
  xyz(1,1) = earthrad*cos(lat1*pi/180.d0)*cos(lon1*pi/180.d0)
  xyz(1,2) = earthrad*cos(lat1*pi/180.d0)*sin(lon1*pi/180.d0)
  xyz(1,3) = earthrad*sin(lat1*pi/180.d0)

! convert lat/lon of mascon 2 into XYZ (assume spherical for now)
  xyz(2,1) = earthrad*cos(lat2*pi/180.d0)*cos(lon2*pi/180.d0)
  xyz(2,2) = earthrad*cos(lat2*pi/180.d0)*sin(lon2*pi/180.d0)
  xyz(2,3) = earthrad*sin(lat2*pi/180.d0)


! calculate the distance
  dist = sum(sqrt( (xyz(1,:)-xyz(2,:))**2))

  return
  end subroutine MASCON_msc_to_msc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine MASCON_transfer_ternary(debug,take_from,give_to,nmsc_prim_in,max_itern,nvar_prim,nvar_tern,itern,msc_ternary1 &
                                     ,msc_ternary2,msc_ternary3,msc_prim1, msc_prim2, num_tern_old)

! subroutine to transfer a ternary mascon from one primary to another. It involves a few steps:
!
! 1. changing the number of mascons in each primary
! 2. recalculating the area of each primary
! 3. (eventually) recalculating the coords of centre of mass of each primary   !@ don't know how to do this yet!!
!
! P. Tregoning
! 29 October 2014

  implicit none

  logical  ,    intent(inout) :: debug
  integer*4,    intent(in)    :: take_from,give_to     ! primary mascon numbers that are losing and gaining a ternary
  integer*4,    intent(in)    :: num_tern_old          ! original number of ternary mascons in the mixed primary mascon
  integer*4,    intent(in)    :: nmsc_prim_in          ! number of input primary mascons
  integer*4,    intent(in)    :: max_itern             ! maximum number of ternary mascons
  integer*4,    intent(in)    :: itern                 ! number of the ternary mascon (in the primary mascon) to be transferred 
  integer*4,    intent(in)    :: nvar_prim,nvar_tern   ! numbers of attributes (columns) of primary and ternary mascon arrays
! ternary mascon information is (3rd column)
! 1. lat, 2. lon, 3. radius, 4. area, 5. density, 6. depth, 7. geoid/ellipsoid separation. 
! PT170531: it is dimensioned as 8 in the main program!!
! PT161104: increase by a factor of 3 the maximum number of ternary mascons (to allow for reshaping of them beyond the input sizes)
  real(kind=8), intent(inout) :: msc_ternary1(max_itern,nvar_tern)   ! information on the original set of ternary mascons to lose a ternary ( lat/lon/rad/area/depth/densitygeoid )
  real(kind=8), intent(inout) :: msc_ternary2(max_itern,nvar_tern)   ! information on the modified set of ternary mascons to gain a ternary
  real(kind=8), intent(inout) :: msc_ternary3(max_itern,nvar_tern)   ! information on the set of ternary mascons to gain a ternary
! primary mascon information is:
! 1. lat, 2. lon,  3. % land, 3. number of ternary mascons
! PT170531: these should be "10", not "9"
  real(kind=8), intent(inout) :: msc_prim1(nvar_prim)          ! number of ternary mascons in primary mascon to lose a ternary
  real(kind=8), intent(inout) :: msc_prim2(nvar_prim)          ! number of ternary mascons in primary mascon to gain a ternary

! local variables
  integer*4  :: i,dtern
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!DEBUG
if(debug)then
  print*,'input variables in mascon_transfer_ternary:'
  print*,'take_from,give_to',take_from,give_to
  print*,'num_tern_old',num_tern_old
  print*,'nmsc_prim_in',nmsc_prim_in
  print*,'max_itern,itern',max_itern,itern
  print*,'nvar_prim,nvar_tern',nvar_prim,nvar_tern
endif


  if(debug)print*,'taking from primary',take_from,' and giving to primary',give_to
  if(debug)print*,'there are currently ',num_tern_old-nint(msc_prim1(8)),' fewer ternaries than to begin with',num_tern_old,msc_prim1(8) !PT210902: fixed but. Index was 5 but should be 8
  dtern = num_tern_old-int(msc_prim1(8))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first, give the ternary away to the other primary mascon
  msc_prim2(8) = msc_prim2(8) + 1   ! add one to the number of ternary mascons in the receiving primary mascon
  msc_ternary3(int( msc_prim2(8) ),:) = msc_ternary1(itern,:)
  if(debug)print*,'new line in receiving primary:',msc_ternary3(int( msc_prim2(8) ),:),msc_prim2(8) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now, remove the ternary from the list by shuffling all the others over the top of it
  msc_prim1(8) = msc_prim1(8) - 1   ! reduce by one the number of ternary mascons in the donating primary mascon
! DEBUG
  if(debug)print*,'    remove the entire row:',int(msc_prim1(8)),itern,msc_ternary1(itern,:)
  do i=itern,max_itern-1
!    print*,'replacing row ',i-dtern,i,dtern
    msc_ternary2(i-dtern,:) = msc_ternary1(i+1,:)
  enddo
  msc_ternary2(max_itern,:) = 0.d0
  if(debug)print*,'    the replaced row is  :',int(msc_prim1(8)),itern,msc_ternary2(itern-dtern,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! we're done
  if(debug)print*,"we're done. Return to calling code"
  return
  end subroutine MASCON_transfer_ternary


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine MASCON_percent_land(debug,nvar_tern,ntern,max_tern,msc_ternary,percent_land)

! subroutine to compute the percentage of land in a primary mascon, based on the ternary mascons
!
! P. Tregoning
! 30 October 2014

  implicit none

  logical     , intent(in ) :: debug
  integer*4   , intent(in ) :: max_tern                        ! maximum number of ternary mascons per primary mascon
  integer*4   , intent(in ) :: ntern                           ! number of ternary mascons in the primary mascon
  integer*4   , intent(in)  :: nvar_tern                       ! number of attributes (ie columns) of ternary mascon arrays
  real(kind=8), intent(in)  :: msc_ternary(max_tern,nvar_tern) ! the ternary mascons for a particular primary mascon
  real(kind=8), intent(out) :: percent_land                    ! percentage (between 0 and 1) of land or -1 for gia mascons

! local variables
  integer*4 :: nland, nocean, ngia, i, j
  real*8    :: pi

  pi = 4.d0*atan(1.d0)

! trap to catch poor input data
  if (ntern == 0 ) then
    call status_update('FATAL','UTIL','mascon_percent_land',"Input number of ternary mascons is zero",0)
  endif

! loop over the ternary mascons and count how many are land
  nland = 0
  nocean = 0
  ngia = 0
  do i=1,ntern
    !print*,'percent_land: i,msc_ternary(i,5)',i,msc_ternary(i,:)
    if (msc_ternary(i,5) >= 0.d0 .and. msc_ternary(i,6) < 3000.d0)then
      nland = nland + 1
    elseif (msc_ternary(i,5) < 0.d0 .and. msc_ternary(i,6) < 3000.d0)then
      !if(debug)print*,'ocean: ',msc_ternary(i,1:2)*180./pi,msc_ternary(i,3:6)
      nocean = nocean + 1
    else
      ngia = ngia + 1
    endif
  enddo

! calculate the percentage of land
  if (ngia == 0) then
    percent_land = dble(nland)/dble(ntern)
  else ! RM190522: set GIA mascons to -1
    percent_land = -0.01
  endif
  if(debug)print*,'ntern,nland,nocean,ngia,percent_land',ntern,nland,nocean,ngia,percent_land

  return
  end subroutine MASCON_percent_land
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine MASCON_nearest_mscs(debug,imsc,nsurround,total_prim,msc_tern_mult,max_tern_per_prim,nvar_tern &
                                ,nvar_prim,msc_prim_in,msc_tern_out &
                                ,max_distance,msc_surrounding,required_density,nposs)
!
! subroutine to find the nearest "nsurround" primary mascons to a particular primary mascons "imsc"
!
! P. Tregoning
! 31 October 2014

  implicit none 

  logical,      intent(inout)    :: debug
  integer*4,    intent(in )   :: imsc                 ! the particular primary mascon
  integer*4,    intent(in )   :: nsurround            ! the number of surrounding mascons that we want to find
  integer*4,    intent(in )   :: total_prim           ! number of primary mascons
  integer*4,    intent(in )   :: msc_tern_mult,max_tern_per_prim,nvar_tern ! needed to dimension the msc_tern_in array
  integer*4,    intent(in )   :: nvar_prim            ! number of primary mascon variables
  real(kind=8), intent(in )   :: msc_prim_in(total_prim,nvar_prim)  ! primary mascon coords, % land and # ternary mascons, area, blank, av depth.
  real(kind=8), intent(in )   :: msc_tern_out(total_prim*msc_tern_mult,max_tern_per_prim*msc_tern_mult,nvar_tern)  ! Ternary mascon coords, % land and # ternary mascons, area, blank, av depth.
  real(kind=8), intent(in )   :: max_distance         ! maximum distance away that each mascon can be considered a nearest neighbour
  real(kind=8), intent(in )   :: required_density     ! density required to match the density of ternaries being given away
  real(kind=8), intent(out)   :: msc_surrounding(nsurround,4) ! info on surrounding primary mascons

! local variables
  integer*4    :: i,j,itern
  real(kind=8),allocatable :: dist(:)                 ! distance from each mascon to the particular mascon
  real(kind=8),allocatable :: xyz(:,:)                ! distance from each mascon to the particular mascon
  real(kind=8), parameter  :: earthrad = 6378136.46 
  real(kind=8)             :: pi
  integer*4                :: nposs                   ! number of possible mascons within X km of the requested mascon
  real(kind=8)             :: check_dists(1000,2)      ! distance info on mascons within X km of the requested mascon. Values are 1. distance, 2. mascon number.
  real(kind=8)             :: tmp_dist(2)
  real(kind=8)             :: dist_tern

  pi = 4.0*atan(1.d0)


! allocate variables
  allocate(dist(total_prim ))
  allocate( xyz(total_prim ,3))

! calculate the coordinates of the particular mascon
  xyz(imsc,1) = earthrad*cos(msc_prim_in(imsc,1))*cos(msc_prim_in(imsc,2))
  xyz(imsc,2) = earthrad*cos(msc_prim_in(imsc,1))*sin(msc_prim_in(imsc,2))
  xyz(imsc,3) = earthrad*sin(msc_prim_in(imsc,1))

! loop through all primary mascons and calculate the distance from each to the particular mascon
  nposs = 0
  do i=1,total_prim

! PT210903: I don't want the land primaries containing Iceland to be used by Greenland primaries
! for 300 km mascons_stage2 it is mascon 222
! for 200 km mascons_stage2_200km it is 479 and 557
    if(total_prim == 4582 .and. (imsc == 221 .or. imsc ==222 .or. imsc==223 .or. imsc==275 .or. imsc==276))then      ! it is likely that this is mascons_stage2_300km
      cycle
    elseif(total_prim == 10314 .and. ((imsc>=478.and.imsc<=480).or.(imsc>=557.and.imsc<=559)) )then ! it is likey that this is mascons_stage2_200km
      cycle
    endif

!if(i == 61)debug = .true.
    if(debug)print*,'i,imsc,nint(msc_prim_in(imsc,8)',i,imsc,nint(msc_prim_in(imsc,8))
! PT161103: add the condition that there are at least some ternary mascons in the primary under consideration
    if (i /= imsc .and. nint(msc_prim_in(i,8)) > 0 .and. dabs(msc_prim_in(i,6) - required_density) < 1.d-3 ) then
      xyz(i,1) = earthrad*dcos(msc_prim_in(i,1))*dcos(msc_prim_in(i,2))
      xyz(i,2) = earthrad*dcos(msc_prim_in(i,1))*dsin(msc_prim_in(i,2))
      xyz(i,3) = earthrad*dsin(msc_prim_in(i,1))
      dist(i) = sum(dsqrt( (xyz(i,:)-xyz(imsc,:))**2))
      if(debug)print*,'xyz of mascon and of mascon i',i,xyz(i,:),msc_prim_in(i,1),msc_prim_in(i,2),earthrad
      if(debug)print*,imsc,i,dist(i),max_distance
      if(dist(i) < max_distance)then
      if(debug)print*,'close mascon',imsc,i,dist(i),max_distance,msc_prim_in(i,1)*180./pi,msc_prim_in(i,2)*180./pi &
         ,msc_prim_in(imsc,1)*180./pi &
         ,msc_prim_in(imsc,2)*180./pi
        nposs = nposs + 1

        ! PT210903: calculate the shortest distance to any of the ternary mascons within this primary. Use it as 
        !           "dist(i)" rather than the distance to the primary CoM
        do itern = 1,nint(msc_prim_in(i,8))  ! use the ternaries in the reshaped mascons here
          if( dabs(msc_tern_out(i,itern,6) - required_density) < 1.d-3 )then
            xyz(i,1) = earthrad*dcos(msc_tern_out(i,itern,1))*dcos(msc_tern_out(i,itern,2))
            xyz(i,2) = earthrad*dcos(msc_tern_out(i,itern,1))*dsin(msc_tern_out(i,itern,2))
            xyz(i,3) = earthrad*dsin(msc_tern_out(i,itern,1))
            dist_tern = sum(dsqrt( (xyz(i,:)-xyz(imsc,:))**2))
if(debug)print*,'imsc,i,itern,dist,dist_tern',imsc,i,itern,dist(i),dist_tern
            if(dist_tern < dist(i))dist(i) = dist_tern 
          endif           
        enddo
        check_dists(nposs,1) = dist(i)
        check_dists(nposs,2) = dble(i)
        if(debug)print*,'near   mascon',imsc,i,dist(i),max_distance,msc_prim_in(i,1)*180./pi,msc_prim_in(i,2)*180./pi &
         ,msc_prim_in(imsc,1)*180./pi &
         ,msc_prim_in(imsc,2)*180./pi
      else
        if(debug)then
          if(dabs(msc_prim_in(i,6) - required_density) < 1.d-3 )then
            print*,'far   mascon',imsc,i,dist(i),max_distance,msc_prim_in(i,1)*180./pi,msc_prim_in(i,2)*180./pi &
            ,msc_prim_in(imsc,1)*180./pi &
            ,msc_prim_in(imsc,2)*180./pi
          else
            print*,'wrong density  mascon',imsc,i,dist(i),max_distance,msc_prim_in(i,1)*180./pi,msc_prim_in(i,2)*180./pi &
            ,msc_prim_in(imsc,1)*180./pi &
            ,msc_prim_in(imsc,2)*180./pi
          endif
        endif

      endif
    endif
!if(i==61)debug = .false.
  enddo

! we have now identified our closest mascons.
! run a bubble sort to put them in ascending order 
  do i = nposs , 1, -1 
    do j=1,i
      if(check_dists(j,1) > check_dists(i,1) ) then
! swap the distances
        tmp_dist(1) = check_dists(i,1)
        check_dists(i,1) = check_dists(j,1)
        check_dists(j,1) = tmp_dist(1)
! and the pointers
        tmp_dist(2) = check_dists(i,2)
        check_dists(i,2) = check_dists(j,2)
        check_dists(j,2) = tmp_dist(2)
      endif
    enddo
  enddo

! ok, now we want the first "nsurround" of them
  do i=1,nsurround
    if(i <= nposs)then
      msc_surrounding(i,1:2) = msc_prim_in(nint(check_dists(i,2)),1:2)
      msc_surrounding(i,3)   = check_dists(i,2)
      msc_surrounding(i,4)   = msc_prim_in(nint(check_dists(i,2)),7) ! PT161028: changed to "7" from "3" to match new order of variables in msc_prim_in array
    else
      msc_surrounding(i,:) = 0.d0
    endif

! DEBUG
if(imsc == 36)debug = .true.
  if(debug .and. check_dists(i,1) > 0.d0)print*,'surround mascon ',int(check_dists(i,2)),' dist',check_dists(i,1),' info: ',msc_prim_in(int(check_dists(i,2)),:)
  enddo
if(imsc == 36)debug = .false.

! deallocate variables
  deallocate(dist)
  deallocate(xyz)

  return
  end subroutine MASCON_nearest_mscs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine MASCON_read_seed(luseed,max_seed,nseed,msc_seed)

! subroutine to read a file of "dummy" primary mascons that contain zero ternary mascons. These coords
! then act as seed primary mascons to draw transfers of ternary mascons into regions (typically small
! land regions) in places in the world where the mixed mascons can't transfer to neighbouring regions.
!
! P. Tregoning
! 4 November 2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  implicit none

  integer*4,   intent(in)   :: luseed                ! unit number of input seed primary mascon file
  integer*4,   intent(in)   :: max_seed              ! maximum number of input seed primary mascon file
  real(kind=8),intent(out)  :: msc_seed(max_seed,6)  ! coords, density, max dist over which to search for nearest surrounding mascons
  integer*4,   intent(out)  :: nseed                 ! number of seed primary mascons found in input file

! local variables
  integer*4    :: ioerr,i,j
  real(kind=8) :: tmpvals(6)       ! temp storage of values read from input file (coatdeg,colatmin,londeg,lonmin,density, dist)
  character    :: message*200

! read in all the values
  ioerr = 0
  do while (ioerr == 0)
    read(luseed,*,iostat=ioerr,end=1000)(tmpvals(j),j=1,6)
    if(ioerr == 0)then
      nseed = nseed + 1
      if(nseed > max_seed) then
        write(message,'(a,i5,a,i5,a)')"Number of seed primary mascons (",nseed,") exceeds max number (",max_seed,")"
        call status_update('FATAL','UTIL','mascon_read_seed',' ',message,0)
      endif
      msc_seed(nseed,:) = tmpvals(:)
    endif
  enddo

1000 return
  end subroutine MASCON_read_seed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine MASCON_calc_CoM(max_tern,nvar_prim,nvar_tern,msc_tern,msc_prim_out)

! subroutine to calculate the new centre of mass (i.e. the "coordinates") of the new primary mascons after
! the redistribution of ternary mascons.
!
! P. Tregoning
! 4 November 2014
!
! PT161119: add new logic to handle the cases where the ternary mascons span 0E longitude (where 
!           we must average , e.g. 359E and 2E)

  implicit none

  integer*4       :: max_tern              ! maximum number of ternary mascons per primary mascon
  integer*4       :: max_tern2
  integer*4       :: nvar_prim             ! number of primary mascon attributes
  integer*4       :: nvar_tern             ! number of ternary mascon attributes
!  parameter (nvar_prim = 9)
  real(kind=8)    :: msc_tern(max_tern,nvar_tern)  ! list of ternary mascons of new primary mascon  (coords in radians?)
  real(kind=8)    :: msc_prim_out(nvar_prim)       ! coords, %land, # ternary mascons, area

! local variables
  real(kind=8) :: tmpsum,tmparea,tmpdepth,tmpgeoid,pi
  integer*4    :: i,j,ntern
  real(kind=8) :: minlon,maxlon
  real(kind=8),allocatable :: longitudes(:) 
  logical      :: debug
  character*200 :: message

  pi = 4.d0*datan(1.d0)

  allocate(longitudes(1000000))

  debug = .false.

! ntern is the number of ternary mascons
  ntern = nint(msc_prim_out(8))

! PT141117: trap for when there are no ternary mascons at all in this primary
  if (ntern == 0)then
    write(message,'(a,i6,a)')"There are ",ntern," ternary mascons. Cannot compute a centroid"
    call status_update('WARNING','UTIL','calc_CoM',' ',message,0)
    msc_prim_out(1:2) = 0.d0
  endif


! Find the centre of mass of longitude
! check whether the longitudes cross the 0E line. If so, make the westward ones negative
  minlon = 2.d0*pi
  maxlon = 0
  do i=1,ntern
    if(msc_tern(i,2) > maxlon)maxlon = msc_tern(i,2)
    if(msc_tern(i,2) < minlon)minlon = msc_tern(i,2)
  enddo
if(debug)then
  print*,'diff in longitude',maxlon,minlon,(maxlon - minlon)*180.d0/pi
endif


  do i=1,ntern
    longitudes(i) = msc_tern(i,2) 
! PT181123: change this from 340 deg to 330 deg (Spitzbergen and Greenland are separated by 337 deg of longitude!)
    if ((maxlon - minlon)*180.d0/pi > 300.d0)then  ! our longitudes have crossed 0E in that they are spearated by more than 330 degrees
!if(debug)print*,'shifting longitudes',(maxlon - minlon)*180.d0/pi
! PT181124: this code is wrong! We need to add 2pi to the low longitudes (ie 0-30deg)
!      longitudes(i) = longitudes(i) + pi/2.d0
!      if(longitudes(i) > 2.d0*pi)then
!if(debug)print*,'shifting back this longitude',longitudes(i)*180./pi
!        longitudes(i) = longitudes(i) - 2.d0*pi
!      endif

      if(longitudes(i) < pi/2.d0)then
if(debug)print*,'orig longitude, pi/2',longitudes(i)*180./pi,pi/2.d0
        longitudes(i) = longitudes(i) + 2.d0*pi
if(debug)print*,'shifting by 2*pi this longitude',i,longitudes(i)*180./pi
      endif

    endif
  enddo

  tmpsum   = 0.d0
  tmparea  = 0.d0
  tmpdepth = 0.d0
  tmpgeoid = 0.d0
!$OMP PARALLEL DO private (i) &
!$OMP             ,shared(ntern,longitudes,msc_tern) &
!$OMP             ,reduction(+:tmpsum,tmparea,tmpdepth,tmpgeoid)
  do i=1,ntern
    tmpsum   = tmpsum + longitudes(i)*msc_tern(i,4)/1.d6
    tmparea  = tmparea + msc_tern(i,4)/1.d6
    tmpdepth = tmpdepth + msc_tern(i,5)
    tmpgeoid = tmpgeoid + msc_tern(i,8)
  enddo
!$OMP END PARALLEL DO
  msc_prim_out(2) = tmpsum/tmparea
! PT161123: undo, if necessary, the longitude shifting that made the longitudes monotonic
! PT181123: change this from 340 deg to 330 deg (Spitzbergen and Greenland are separated by 337 deg of longitude!)
  if(msc_prim_out(2) > 2.d0*pi)msc_prim_out(2) = msc_prim_out(2) - 2.d0*pi

! save the average depth and geoid/ellipsoid separation
  msc_prim_out(5) = tmpdepth / dble(ntern)
  msc_prim_out(12) = tmpgeoid / dble(ntern)

! now for latitude
  tmpsum  = 0.d0
  tmparea = 0.d0
!$OMP PARALLEL DO private (i) &
!$OMP             ,shared(ntern,msc_tern) &
!$OMP             ,reduction(+:tmpsum,tmparea)
  do i=1,ntern
    tmpsum = tmpsum + msc_tern(i,1)*msc_tern(i,4)/1.d6
    tmparea = tmparea + msc_tern(i,4)/1.d6
  enddo
!$OMP END PARALLEL DO
  msc_prim_out(1) = tmpsum/tmparea

! PT141111: save the area of the primary mascon
  msc_prim_out(4) = tmparea*1.d6

  debug = .false.
  if(debug)then
    print*,'New centroid coords                           ',msc_prim_out(1:2)*180./pi,' area',msc_prim_out(4)
  endif

! that's all
  return
  end subroutine MASCON_calc_CoM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


