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



