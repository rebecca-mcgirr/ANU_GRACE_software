!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Subroutines related to mascon computations done inside graceorb
!!
!! calc_mascon_acc_part    ::  calculates the gravitational acceleration and partials for mascons
!! generate_mascon_vector  :: generates the list of mascons (combination of primary, secondary, ternary) for each epoch
!! read_ternary_EWH        :: calculates mean EWH values for primary  and secondary mascons, given EWH values on ternary mascons
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calc_mascon_acc_part(efpos,t,mcon_efacc)

! subroutine to handle the computation of accelerations from gravity mascons as well as
! their partial derivatives
! Code moved here from gravcalc.f90
!
! P. Tregoning
! 9 September 2016

  use mascon_mod       ! all the required mascon arrays
  use gm_mod           ! this gives us Gconst
  use timing_mod
  use inmod_mod        ! this gives us gt_mcon(1), whether to compute mascon partials or not

  implicit none

  real(kind=8),intent(in)  :: efpos(3)              ! satellite XYZ (earth-fixed)
  real(kind=8),intent(in)  :: t                     ! time in seconds
  real(kind=8),intent(out) :: mcon_efacc(3)         ! sum of accelerations from the a priori mascon EWH values

! local variables
  integer*4      :: imsc,itern,k,m,j
  real(kind=8)   :: mcon_dist,tern_mcon_dist
  real(kind=8)   :: tmp_part(3)
  real(kind=8)   :: def_fact,mcon_fact,tmp_mcon_efacc          ! temporal storage of values used in repeated computations
  logical        :: on_ternaries

! debug variables
  real(kind=8)   :: msc_sph(3),sat_sph(3)

! initialise values to zero
  mcon_efacc = 0.d0
! PT140428: initialize the mascon partials
  if(gt_mcon(1) == "Y")mcon_efacc_part = 0.d0

  on_ternaries = .true.
  if(on_ternaries)then

    mcon_efacc = 0.d0
!$OMP PARALLEL DO private (mcon_dist,tern_mcon_dist,def_fact,mcon_fact,tmp_part,j,imsc,itern) &
!$OMP         shared(mcon_tern,efpos,GConst,gt_mcon) &
!$OMP         reduction(+:mcon_efacc_part,mcon_efacc)
    do imsc = 1, total_prim

! PT211231: compute the distance satellite to primary mascon. If > 8000 km, compute using point source
      mcon_dist = dsqrt((mcon_prim(imsc,1)-efpos(1))**2 + (mcon_prim(imsc,2)-efpos(2))**2 + (mcon_prim(imsc,3)-efpos(3))**2)

      ! PT221129: reduce from 20,000 km to 1000 km
      if(mcon_dist < 1000.d3)then
        ! PT211231: compute on the ternaries and then sum it up for the primary
        do itern = 1,nint(mcon_prim(imsc,8))
          tern_mcon_dist = dsqrt( (mcon_tern(imsc,itern,1)-efpos(1))**2 + (mcon_tern(imsc,itern,2)-efpos(2))**2 &
                + (mcon_tern(imsc,itern,3)-efpos(3))**2)

!DEBUG:
!if(imsc == 292 .and. itern == 1)print*,"tern_mcon_dist = ",tern_mcon_dist
          def_fact = dble(mcon_tern(imsc,itern,6) * mcon_tern(imsc,itern,4))  ! itern,6 is density, itern,4 is area
          ! PT220102: use the sat/tern distance in km, then correct later for 1e-9
          mcon_fact = Gconst * def_fact/(tern_mcon_dist*1.e-3)**3

          do j=1,3
            tmp_part(j) = (mcon_tern(imsc,itern,j)-efpos(j)) * mcon_fact ! PT160909: the deformation components are no longer computed + def_fact * up_def_vector(j,m) + def_fact * tang_def_vector(j,m)
            ! PT210806: only save the partial if we want partials for the mascons
            if(gt_mcon(1) == "Y")mcon_efacc_part(imsc,j)   = mcon_efacc_part(imsc,j)   + tmp_part(j)
            mcon_efacc(j) = mcon_efacc(j) + tmp_part(j)*1.e-9 * mcon_tern_EWH(nint(mcon_tern(imsc,itern,7))) 
          enddo
! DEBUG
! PT220104: print out the lat/lon of the mascon and the satellite
          call cart_to_sph(mcon_tern(imsc,itern,1:3),msc_sph)
          call cart_to_sph(efpos(1:3),sat_sph)
!          print*,msc_sph(1:2)*180.d0/pi,sat_sph(1:2)*180.d0/pi,t,' ternary mascon and sat spherical'

        enddo  ! end of loop over ternary mascons within a primary
        
        ! PT220102: now rescale the partial by 1.e-9
        if(gt_mcon(1) == "Y")mcon_efacc_part(imsc,:) = mcon_efacc_part(imsc,:) * 1.e-9

      
      else
        ! PT211231: treat the primary as a point source
        def_fact = dble(mcon_prim(imsc,6) * mcon_prim(imsc,4))  ! PT211231: imcs,6 is density, imsc,4 is area
        mcon_fact = Gconst * def_fact/mcon_dist**3
        
        do j=1,3
          tmp_part(j) = (mcon_prim(imsc,j)-efpos(j)) * mcon_fact ! PT160909: the deformation components are no longer computed + def_fact * up_def_vector(j,m) + def_fact * tang_def_vector(j,m)
          ! PT210806: only save the partial if we want partials for the mascons
          if(gt_mcon(1) == "Y")mcon_efacc_part(imsc,j)   = mcon_efacc_part(imsc,j)   + tmp_part(j)
          mcon_efacc(j) = mcon_efacc(j) + tmp_part(j) * mcon_prim_EWH(imsc) 
        enddo
              
      endif
! DEBUG
! PT220104: print out the lat/lon of the mascon and the satellite
        call cart_to_sph(mcon_prim(imsc,1:3),msc_sph)
        call cart_to_sph(efpos(1:3),sat_sph)
!        print*,msc_sph(1:2)*180.d0/pi,sat_sph(1:2)*180.d0/pi,t,' primary mascon and sat spherical'

!! DEBUG: PT211231
!if(imsc == 908)then
!  mcon_dist = dsqrt((mcon_prim(imsc,1)-efpos(1))**2 + (mcon_prim(imsc,2)-efpos(2))**2 + (mcon_prim(imsc,3)-efpos(3))**2)
!  print*,t,mcon_efacc_part(imsc,:),mcon_dist,k,mcon_num(imsc),'  scaled,t,mcon_efacc_part,mcon_dist'
!!  print*,mcon_prim(imsc,1:3),efpos(:)
!endif

    enddo  ! end of loop over primary mascons
!$OMP END PARALLEL DO

    !stop ' stopped after one epoch of 2700 km limited'

  else
    ! first, decide whether each primary mascon should be broken into secondary/ternary and make a vector of them
    call generate_mascon_vector(efpos,t)

    ! now, compute the required stuff for the gravitational effect of the deformation of the load
    ! PT160907: comment this out for now .... it takes a lot of time to compute and we think it is insignificant anyway
    !    call defgravcalc(efpos)

    ! PT160914: uncomment this if you want to time the acceleration computation
    !   call system_clock ( clock_count1, clock_rate, clock_max )

    ! loop over all the mascons in the vector of mascons
    do imsc = 1,total_prim
      tmp_mcon_efacc = 0.d0

      do k = 1, mcon_num(imsc)
        ! PT130904: derive the value of m explicitly so that it will work in openMP
        m = sum(mcon_num(1:imsc-1))+k
        mcon_dist = dsqrt((mcon_xyz(1,m)-efpos(1))**2 + (mcon_xyz(2,m)-efpos(2))**2 + (mcon_xyz(3,m)-efpos(3))**2)

        def_fact = dble(mcon_rho(m) * mcon_area(m)) 
        mcon_fact = Gconst * def_fact/mcon_dist**3
        do j=1,3
          tmp_part(j) = (mcon_xyz(j,m)-efpos(j)) * mcon_fact ! PT160909: the deformation components are no longer computed  "+ def_fact * up_def_vector(j,m) + def_fact * tang_def_vector(j,m) "
          ! PT210806: only save the partial if we want partials for the mascons
          if(gt_mcon(1) == "Y")mcon_efacc_part(imsc,j)   = mcon_efacc_part(imsc,j)   + tmp_part(j)
          mcon_efacc(j) = mcon_efacc(j) + tmp_part(j) * mcon_EWH_vector(m) 

        enddo
!if(imsc==10315)print*,'debug for 10315',t,tmp_part(:) * mcon_EWH_vector(m)

! DEBUG
! PT220104: print out the lat/lon of the mascon and the satellite
        call cart_to_sph(mcon_xyz(:,m),msc_sph)
        call cart_to_sph(efpos(1:3),sat_sph)
!        print*,msc_sph(1:2)*180.d0/pi,sat_sph(1:2)*180.d0/pi,t,imsc,k,m,' mascon and sat spherical'
      enddo     ! end of (sub)mascon loop
    end do    !  end of mascon loop
    ! PT160914: uncomment to time and output the tide acceleration computation
    !  call system_clock ( clock_count2, clock_rate, clock_max )
    !  print '("mascon acceleration loop Time = ",f10.3," seconds. Mascon efacc=",3f10.3," seconds. t=",f10.2)' &
    !          ,real(clock_count2-clock_count1)/real(clock_rate) &
    !          ,mcon_efacc*1.e9,t
    
    !stop 'stopped after one epoch of dist_flag approach'
    
  endif



  return
  end subroutine calc_mascon_acc_part
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calc_ocean_acc_part(sat,rlat,rlong,efpos,jd,t,ocean_efacc,output_tide)

! subroutine to handle the computation of accelerations from ocean mascons as well as
! their partial derivatives
! Code moved here from gravcalc.f90
!
! P. Tregoning
! 9 September 2016
!
! MODS:
! PT170125: added the computations of the partial derivatives to allow tide amplitudes to be estimated.

  use mascon_mod       ! all the required mascon arrays
  use gm_mod           ! this gives us Gconst
!  use tide_mod         ! this gives us all the tidal grid variables
  use tide_netcdf_mod  ! definitions/functions for interfacing with the tide grid file
  use timing_mod
  use accel_mod


  implicit none

  character*1 ,intent(in)  :: sat                   ! which GRACE satellite
  real(kind=8),intent(in)  :: rlat,rlong            ! satellite lat/lon (in radians)
  real(kind=8),intent(in)  :: efpos(3)              ! satellite XYZ (earth-fixed)
  integer*4   ,intent(in)  :: jd                    ! julian day of orbit
  real(kind=8),intent(in)  :: t                     ! time in seconds
  real(kind=8),intent(out) :: ocean_efacc(3)        ! sum of accelerations from the ocean tide height values
  logical,intent(in)       :: output_tide           ! logical as to whether to output satellite/tide info to the screen

! ocean tide height variables
  real(kind=8) ::  ocean_ht         ! value of ocean height passed back from SPOTL routines
  real(kind=8) ::  ocean_tide_corr  ! computation of tidal corrections using a priori mascon tide amplitudes
  real(kind=8) ::  dt_tides         ! seconds since 2000/01/01 12:00:00
  real(kind=8) ::  time             ! time of epoch in decimal julian day
  integer*4    ::  timearr(5)       ! yr,doy,hr,min,sec to calculate ocean tide model
  integer*4    ::  date(5)          ! yr,mo,day,hr,min  (used to get from "time" to "timearr")
  real(kind=8) ::  sec              ! seconds of hour
  character*20 ::  grd_file         ! name of tidal grid file
  character*5  ::  upperc           ! upper case function
  real(kind=8) ::  secs_of_day      ! variable to output the seconds of day since 00UT
  real(kind=8) ::  dtime

! variables for the computation of the ocean tide acceleration
  real(kind=8) :: def_fact,ocean_mcon_dist,mcon_fact,tmp_part(3)
  real(kind=8) :: ocean_efacc_noampl(3), ocean_efacc_ampl(3)     


! local variables
  integer*4     :: i,j,l,itern,imsc
  integer*4     :: irec
  character*250 :: message
  integer*4     :: tern_number
  real(kind=8)  :: tide_ht,tide_range(2)
  logical       :: debug
  logical       :: bitmap
  integer*4     :: ocean_tern

! debug variables
  real*8 :: tmplat,tmplon,xyz(3),xyz_ocean(3),amag3

  rad_fact = 4.d0*atan(1.d0)/180.d0      ! conversion from degrees to radians


  ocean_efacc_noampl = 0.d0
  ocean_efacc_ampl   = 0.d0
  ocean_efacc        = 0.d0
  ocean_eftid_part   = 0.d0

! calculate the time since 0 GRACE seconds (to use in the tidal computations)
  dt_tides = (dble(jd) - 2451545.d0) * 86400.d0 + t   ! seconds since 2000/01/01 12:00:00


  if(t < 43200.d0)then
! PT140829: fix bug of wrap-around when 12hr < t < 24 hr
   tide_epoch1 = int((t+43200.d0)/dt_ocean_grid) + 1
  else
    tide_epoch1 = int((t-43200.d0)/dt_ocean_grid) + 1
  endif

! check whether to update the ocean tide grid to a new temporal epoch
  if(tide_epoch1 > tide_epoch_old) then
    ternary_ht1(:) = ternary_ht2(:)
    tide_epoch_old = tide_epoch1
    call read_tide_at_time(tidedata,tide_epoch1+1,ternary_ht2)
    ternary_ht2 = ternary_ht2 / 100.d0
  endif

! we need to interpolate the two epochs of the tide grids to the current epoch
  dtime = mod(t,dt_ocean_grid)/dt_ocean_grid
  ternary_ht_current = ternary_ht1 + (ternary_ht2-ternary_ht1) * dtime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now we need code to compute the accelerations of tide heights 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     
! PT160914: uncomment this if you want to time the tide acceleration computation
!   call system_clock ( clock_count1, clock_rate, clock_max )
 

! PT170130: now loop through all non-tidal-amplitude ternarys, using OMP to speed it up
!
!$OMP PARALLEL DO private (ocean_mcon_dist,def_fact,mcon_fact,tmp_part,i,j,l,imsc,itern) &
!$OMP         shared(mcon_ocean_prim,mcon_ocean_tern,efpos,GConst,ocean_eftid_part,msc_tide,mcon_tide_period) &
!$OMP         reduction(+:ocean_efacc,ocean_efacc_noampl)
  do i = 1,n_tern_no_ampl
    itern = ocean_ternary_no_ampl(i)
    ocean_mcon_dist = dsqrt((mcon_ocean_tern(itern,1)-efpos(1))**2 + (mcon_ocean_tern(itern,2)-efpos(2))**2 &
                + (mcon_ocean_tern(itern,3)-efpos(3))**2)

    def_fact = dble(mcon_ocean_tern(itern,6) * mcon_ocean_tern(itern,4))  ! itern,6 is density, itern,4 is area
    mcon_fact = Gconst * def_fact/ocean_mcon_dist**3

    do j=1,3
      tmp_part(j) = (mcon_ocean_tern(itern,j)-efpos(j)) * mcon_fact ! PT160909: the deformation components are no longer computed + def_fact * up_def_vector(j,m) + def_fact * tang_def_vector(j,m)
!      mcon_efacc_part(imsc,j)   = mcon_efacc_part(imsc,j)   + tmp_part(j)
      ocean_efacc_noampl(j) = ocean_efacc_noampl(j) + tmp_part(j) * ternary_ht_current(itern) 


    enddo
  enddo
!$OMP END PARALLEL DO

! PT170130: now loop through all tidal-amplitude ternarys, including computations of the contribution of a priori
!           tidal amplitudes to the tide height. Also calculate the partials of the tidal amplitude parameters

  do i=1,n_tern_ampl
    itern = ocean_ternary_ampl(i)
    imsc = nint(mcon_ocean_tern(itern,7))   ! PT161205: this is the pointer to the primary mascon in which this ternary resides

    ocean_mcon_dist = dsqrt((mcon_ocean_tern(itern,1)-efpos(1))**2 + (mcon_ocean_tern(itern,2)-efpos(2))**2 &
                + (mcon_ocean_tern(itern,3)-efpos(3))**2)

    def_fact = dble(mcon_ocean_tern(itern,6) * mcon_ocean_tern(itern,4))  ! itern,6 is density, itern,4 is area
    mcon_fact = Gconst * def_fact/ocean_mcon_dist**3

    do j=1,3
      tmp_part(j) = (mcon_ocean_tern(itern,j)-efpos(j)) * mcon_fact ! PT160909: the deformation components are no longer computed + def_fact * up_def_vector(j,m) + def_fact * tang_def_vector(j,m)
!      mcon_efacc_part(imsc,j)   = mcon_efacc_part(imsc,j)   + tmp_part(j)
      ocean_efacc_ampl(j) = ocean_efacc_ampl(j) + tmp_part(j) * ternary_ht_current(itern) 

!! and the tidal partials
      do l=1,max_msc_tides
!          ! add up the contribution of this ternary to the partial derivative for the tidal amplitude parameters
          ocean_eftid_part(j,l,1,imsc) = ocean_eftid_part(j,l,1,imsc) + tmp_part(j)*dsin(mcon_tide_period(l)*dt_tides)  !   sine component
          ocean_eftid_part(j,l,2,imsc) = ocean_eftid_part(j,l,2,imsc) + tmp_part(j)*dcos(mcon_tide_period(l)*dt_tides)  ! cosine component
!
!          ! also add on the contribution of this ternary to any tidal height from a priori tidal amplitude values
          ! PT170302: the sin and cos terms were the wrong way around in this computation. Fixed it to match the partials above.
          ocean_efacc_ampl(j) = ocean_efacc_ampl(j) + msc_tide(l,1,imsc) * tmp_part(j)*dsin(mcon_tide_period(l)*dt_tides) &
                                          + msc_tide(l,2,imsc) * tmp_part(j)*dcos(mcon_tide_period(l)*dt_tides)

! DEBUG
!if(j == 1 .and. imsc == 2080 )then !.and. l == 1)then
!  print*,'imsc,j,l,ocean_eftid_part(j,l,1,imsc)',imsc,itern,j,l,ocean_eftid_part(j,l,2,imsc)
!endif

      enddo

    enddo

! DEBUG
!    print*,'n_tern_ampl,i,itern,imsc,jd,t:',n_tern_ampl,i,itern,imsc,jd,t


  enddo

! add together the contributions of ternarys with and without tidal amplitudes
  ocean_efacc = ocean_efacc_noampl + ocean_efacc_ampl


! PT160914: uncomment to time and output the tide acceleration computation
!  call system_clock ( clock_count2, clock_rate, clock_max )
!  print '("ternary ocean acceleration loop Time = ",f10.3," seconds. Ocean efacc=",3f10.3," seconds. t=",f10.2)' &
!          ,real(clock_count2-clock_count1)/real(clock_rate) &
!          ,ocean_efacc*1.e9,t






! derive the ocean tide height from the global grid at the location beneath the satellites
!SA 200130 Preprocessor flags to remove this logging information in production mode
#ifndef _NDEBUG
  if( (output_tide .and. mod(t,300.) < 1.e-4)  )then
!  if( (output_tide .and. mod(t,300.) == 0.)  )then
    if(t >= 43200.d0)secs_of_day = t-43200.d0
    if(t <  43200.d0)secs_of_day = t+43200.d0
    debug = .false.
    call calc_which_ternary(debug,rlat*180.d0/pi,rlong*180.d0/pi,10.d0/60.d0,tern_number)

! PT170130: replaced the code to output the tide height, now that it is done inside a subroutine that deals only
!           with ocean ternary mascons
    ocean_tern = mcon_ocean_tern_ptr(tern_number,2)
    if(ocean_tern == 0 )then
      tide_ht = 0.d0          ! the ternary was not in the ocean ternary file, and is therefore over land
    else
      tide_ht = ternary_ht_current(mcon_ocean_tern_ptr(tern_number,2))
    endif
    write(message,'(a,a,a,i6,a,f8.4,a,2f11.5,a,3f10.3,a)')'    GRACE ',sat,' Seconds of day:',int(secs_of_day)&
             ,'      Tide ht (from netcdf ),',tide_ht  &
             ,' m. Coords (',rlat*180.d0/pi,rlong*180.d0/pi, ' )',ocean_efacc*1.e9,' nm/sec^2'
        call status_update('STATUS','GRACEORB','calc_ocean_acc_part',' ',message,0)

  endif
#endif
  return
  end subroutine calc_ocean_acc_part
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine generate_mascon_vector(efpos,t)

! subroutine to generate the list of mascons (primary, secondary or ternary) to be used to calculate
! the accelerations and partials. This replaces Tony's original mascon_calc.f90 code, but is based on it.
!
! P. Tregoning
! 7 September 2016

  use mascon_mod
  use accel_mod

  implicit none

  real(kind=8), intent(in) :: efpos(3)          ! satellite XYZ coordinates (earth-fixed)
  real(kind=8), intent(in) :: t                 ! time in seconds of current epoch

! satellite variables
!  real(kind=8)     :: sat_lat,sat_lon,sat_rad   ! lat, lon, radius of the satellite location


  integer*4        :: num_cell_sat_old          ! value of previous primary cell over which the satellite is located

! local variables
  real(kind=8)     :: pi            
  integer*4        :: i,j,k,iprim,isec,itern
  character*250    :: message
  logical          :: debug

! debug
  real*8           :: msc_sph(3)
  
  pi = 4.d0*atan(1.d0)
  rad_fact = pi/180.d0

! PT210817: if we don't want to use the flag file then just return a vector of the primary mascons
  if(.not. use_flag_file)then
    do iprim = 1,total_prim
      mcon_num(iprim) = 1
      mcon_xyz(:,iprim) = mcon_prim(iprim,1:3)
      mcon_area(iprim)  = mcon_prim(iprim,4)
      mcon_rho(iprim)   = mcon_prim(iprim,6)
    enddo

  else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! determine latitude and longitude occupied by the satellite (in degrees)
    sat_rad = sqrt(efpos(1)**2 + efpos(2)**2 + efpos(3)**2)
    ! PT161026: sat_lat is meant to be "latitude", not "colatitude" !!!
    sat_lat = 90.d0 - acos(efpos(3)/sat_rad) /rad_fact
    sat_lon = atan2(efpos(2), efpos(1)) / rad_fact
    ! PT121221: put longitude between 0 and 360 degrees
    if ( sat_lon .lt. 0.0d0 ) sat_lon = sat_lon + 360.0d0

    ! determine which ternary mascon cell (at highest resolution), the satellite overflies
    debug = .false.
    call calc_which_ternary(debug,sat_lat,sat_lon,ternary_lat_spacing/60.d0,num_cell_sat)

    ! now work out which primary this corresponds to
    num_cell_sat = mcon_tern_ptr(num_cell_sat,1)

    ! PT160920: catch circumstances where the ternary doesn't exist in the input ternary file
    if(num_cell_sat == 0)then
      write(message,'(a,2f10.3)')"Error: there is no ternary mascon in the input mascon file at location:",sat_lat,sat_lon
      call status_update('WARNING','GRACEORB','generate_mascon_vector',' ',message,0)
      num_cell_sat = 1  ! jsut set it to 1 so that the code can continue to work ....
    endif


    ! PT130918: check whether the satellites are over the same or a different primary mascon
   
    if(num_cell_sat == num_cell_sat_old .and. t /= 43200. ) then     ! we don't need to do the computations below
      return
    else                                          ! update the num_cell_sat_old and continue with the computations
      num_cell_sat_old = num_cell_sat
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                    !!
!!                    Loop through primary mascon cells                               !!
!!                                                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    tmp_mcon = 0

! PT/AP220104: this code is for a lower-triangular without diagonal elements    
    do iprim = 1, num_cell_sat - 1
      dist_flag(iprim) = mcon_flag(1+int(dble(iprim - 1) + (dble(num_cell_sat-1) * (dble(num_cell_sat - 2)))/2.d0))
    end do
    
    dist_flag(num_cell_sat) = 3
      
    do iprim = num_cell_sat + 1, total_prim
      dist_flag(iprim) = mcon_flag(2+int(dble(num_cell_sat - 1) + (dble(iprim - 1) * (dble(iprim-2)))/2.d0))
    end do

!! PT220106: lower-triangular WITH diagonal elements (which is what dist_flag output up to end of 2021)
!    do iprim = 1, num_cell_sat - 1
!      dist_flag(iprim) = mcon_flag(int(dble(iprim ) + (dble(num_cell_sat-1) * (dble(num_cell_sat - 2)))/2.d0))
!    end do
!    
!    dist_flag(num_cell_sat) = 3
!      
!    do iprim = num_cell_sat + 1, total_prim
!      dist_flag(iprim) = mcon_flag(int(dble(num_cell_sat ) + (dble(iprim ) * (dble(iprim-1)))/2.d0))
!    end do


!! DEBUG
!print*,'mcon_flag',mcon_flag(67808825:67808845)

    do iprim = 1, total_prim


      ! extract distance flag appropriate for the distance between the current cell and the cell containing the satellite
      j = min(iprim, num_cell_sat)
      k = max(iprim, num_cell_sat) - 1
      !    dist_flag(iprim) = mcon_flag(j + (k * (k + 1))/2)
      ! RM190503: fixed integer multiplication for large values of j and k
      !dist_flag(iprim) = mcon_flag(int(dble(j) + (dble(k) * (dble(k - 1)))/2.d0))

! PT220104: what if the satellite is in mascon "iprim"? Set it to break into ternaries
      !if(num_cell_sat == iprim)dist_flag(iprim) = 3
      

! debug
! PT220104: print out the dist flag value
call cart_to_sph(mcon_prim(iprim,1:3),msc_sph)
!print*,'lopmsc', iprim,' dist_flag',dist_flag(iprim),1+int(dble(num_cell_sat - 1) + (dble(iprim - 1) * (dble(iprim-2)))/2.d0)
!,msc_sph(1:2)*180.d0/pi,num_cell_sat,sat_lat,sat_lon 

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                       !
!     P R I M A R Y    M A S C O N      !
!                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! If distance is greater than the coarsest level (20 degrees) use point-mass approximation
      if ( dist_flag(iprim) .eq. 1 ) then
        mcon_num(iprim) = 1
        tmp_mcon = tmp_mcon + 1
        if ( tmp_mcon .gt. 80 * max_mcon1 ) then
          write(message,'(a,2i8)') "ERROR: combined primary mascons exceed array limit: ", tmp_mcon, 80 * max_mcon1
          call status_update('FATAL','GRACEORB','generate_mascon_vector',' ',message,0)
        endif
        mcon_xyz(:,tmp_mcon) = mcon_prim(iprim,1:3)
        mcon_area(tmp_mcon)  = mcon_prim(iprim,4)
        mcon_rho(tmp_mcon)   = mcon_prim(iprim,6)
        mcon_EWH_vector(tmp_mcon)   = mcon_prim_EWH(iprim)

!DEBUG
!PT220104: output distance to primary to be used a sprimary
!print*,'mascon',iprim,' dist to satellite:',dsqrt( (mcon_prim(iprim,1)-efpos(1))**2 &
!            +(mcon_prim(iprim,2)-efpos(2))**2+(mcon_prim(iprim,3)-efpos(3))**2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  S E C O N D A R Y    M A S C O N 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! If distance is in intermediate range (<= 20, >= 5 degreed) use secondary resolution
      else if ( dist_flag(iprim) .eq. 2 ) then
        mcon_num(iprim) = mcon_prim(iprim,7)
        do isec = 1, int(mcon_prim(iprim,7))
          tmp_mcon = tmp_mcon + 1
          if ( tmp_mcon .gt. 80 * max_mcon1 ) then
            write(message,'(a,2i8)') "ERROR: combined secondary mascons exceed array limit: ", tmp_mcon, 80 * max_mcon1
            call status_update('FATAL','GRACEORB','generate_mascon_vector',' ',message,0)
          endif
          mcon_xyz(:,tmp_mcon) = mcon_sec(iprim,isec,1:3)
          mcon_area(tmp_mcon) = mcon_sec(iprim,isec,4)
          mcon_rho(tmp_mcon) = mcon_sec(iprim,isec,6)
          mcon_EWH_vector(tmp_mcon)   = mcon_sec_EWH(iprim)

        enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  T E R N A R Y    M A S C O N 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! If distance is in near-field range (<= 5 degrees) use ternary resolution
      else
        mcon_num(iprim) = mcon_prim(iprim,8)
!print*,'mcon_prim(iprim,8)=',iprim,dist_flag(iprim),mcon_prim(iprim,8)
!print*,'mascon',iprim,' terndist to satel:',dsqrt( (mcon_prim(iprim,1)-efpos(1))**2 &
!            +(mcon_prim(iprim,2)-efpos(2))**2+(mcon_prim(iprim,3)-efpos(3))**2)

        do itern = 1, int(mcon_prim(iprim,8))
          tmp_mcon = tmp_mcon + 1
          if ( tmp_mcon .gt. 800 * max_mcon1 ) then
            write(message,'(a,2i8)') "ERROR: combined ternary mascons exceed array limit: ", tmp_mcon, 800 * max_mcon1
            call status_update('FATAL','GRACEORB','generate_mascon_vector',' ',message,0)
          endif
          mcon_xyz(:,tmp_mcon) = mcon_tern(iprim,itern,1:3)
          mcon_area(tmp_mcon)  = mcon_tern(iprim,itern,4)
          mcon_rho(tmp_mcon)   = mcon_tern(iprim,itern,6)
          mcon_EWH_vector(tmp_mcon)   = mcon_tern_EWH(nint(mcon_tern(iprim,itern,7)))
        enddo

      endif
      if(dabs(mcon_EWH_vector(tmp_mcon)) > 100.d0)then
        write(message,'(a,i8,f15.5)')'erroneously large EWH value:',tmp_mcon,mcon_EWH_vector(tmp_mcon)
        call status_update('FATAL','graceorb','generate_mascon_vector',' ',message,0)
      endif

    enddo  ! end of loop over all primary mascons

  endif  ! end of if statement on whether or not to use the input flag file for distances
  return
  end subroutine generate_mascon_vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ternary_tidal_ampl

! subroutine to make two separate lists of ternary mascons: those involved in tidal amplitude
! modelling and estimation and those not involved.
!
! P. Tregoning
! 30 January 2017

  use mascon_mod   ! provides the declarations of the arrays and counters used in this subroutine


  implicit none

  integer*4      :: itern, imsc
  character*150  :: message

  allocate(ocean_ternary_no_ampl(total_ocean_tern))
  allocate(ocean_ternary_ampl(total_ocean_tern))

!  first, separate ternarys into those for which we want to estimate tidal amplitudes 
!  and those we don't. This includes modelling the a priori the tidal height from the tidal amplitudes
  n_tern_no_ampl = 0
  n_tern_ampl = 0
  do itern = 1,total_ocean_tern
    imsc = nint(mcon_ocean_tern(itern,7))      ! PT161205: this is the pointer to the primary mascon in which this ternary resides
    if (mcon_ocean_prim(imsc,2) == 0)then      ! no tidal amplitude modelling or estimation involved with this primary/ternary
      n_tern_no_ampl = n_tern_no_ampl + 1
      ocean_ternary_no_ampl(n_tern_no_ampl) = itern
    else                                             ! we do have tidal amplitude modelling and estimation
      n_tern_ampl = n_tern_ampl + 1
      ocean_ternary_ampl(n_tern_ampl) = itern
    endif
  enddo

  write(message,'(a,i8,a,i7,a)')"There are",n_tern_no_ampl,' ocean ternarys without tidal amplitude data and' &
                                 ,n_tern_ampl," in primary mascons to have tidal amplitudes modelled/estimated"
  call status_update('STATUS','GRACEORB','ternary_tidal_ampl',' ',message,0)

  return
  end subroutine ternary_tidal_ampl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_ternary_EWH(luin)

! read a file - in msc_to_plot.dat format - containing a priori EWH values on each ternary mascon.
! The values are stored in vectors mcon_tern_EWH and mcon_prim_EWH
!
! P. Tregoning
! 18 August 2021

  use mascon_mod

  implicit none

! passed variables
  integer*4, intent(in) :: luin

! local variables
  real(kind=8) :: lat,lon,junk(4),tmpval,tmpsum
  integer*4    :: ioerr,iprim,itern,tmptern


! loop through all ternary values in the file
  do itern = 1,max_tern
    read(luin,*,iostat=ioerr)lat,lon,tmpval,junk,tmptern
    if(ioerr ==0)then
      mcon_tern_EWH(tmptern) = tmpval/1.d3
    endif
  enddo

! now, calculate the average value for each primary mascon. Sum up the ternaries as a volume, then divide by area of the primary
  do iprim = 1, total_prim
    mcon_prim_EWH(iprim) = 0.d0
    do itern = 1,nint(mcon_prim(iprim,8))
      mcon_prim_EWH(iprim) = mcon_prim_EWH(iprim) + mcon_tern_EWH(nint(mcon_tern(iprim,itern,7)))*mcon_tern(iprim,itern,4)
!print*,'iprim,itern,mcon_tern(iprim,itern,7)',iprim,itern,mcon_tern(iprim,itern,7) &
!       ,mcon_tern_EWH(nint(mcon_tern(iprim,itern,7))),mcon_tern(iprim,itern,4),mcon_prim_EWH(iprim)
    enddo
    mcon_prim_EWH(iprim) = mcon_prim_EWH(iprim)/mcon_prim(iprim,4)
    !PT210820: kluge here to make mcon_sec_EWH. This is wrong if the number of secondaries in a primary is NOT 1
    if(mcon_prim(iprim,7) /= 1)then
      call status_update('FATAL','GRACEORB','read_ternary_EWH',' ','code fails at present if nsec per prim > 1',0)
    else
      mcon_sec_EWH(iprim) = mcon_prim_EWH(iprim)
    endif
  enddo

! PT/RM 210819: write out a file of EWH on primary mascons.
  open(50,file='mcon_prim_EWH.dat',status='unknown')
  write(50,*)"Primary Mascon    EWH (m)"
  do iprim = 1,total_prim
    write(50,'(i5". MC",i5.5," (m)            ",f17.7,f17.7,f17.7)')iprim,iprim,0.d0,mcon_prim_EWH(iprim),mcon_prim(iprim,8)
  enddo
  close(50)
  call status_update('STATUS','GRACEORB','read_ternary_EWH','mcon_prim_EWH.dat',"Output average EWH values on primary mascons",0)

! PT210826: calculate and output the RMS of the ternary anomalies about the mean of each primary
  open(51,file='mcon_prim_EWH.rms',status='unknown')
  write(51,*)"V2 GRACEORB   VCV "
  write(51,*)" PARAMETER                     ZEROES             RMStern            EWHprim"
  do iprim = 1,total_prim
    tmpsum = 0.d0
    do itern = 1,nint(mcon_prim(iprim,8))
!print*,'ternary residuals',iprim,itern,mcon_prim_EWH(iprim) , mcon_tern_EWH(nint(mcon_tern(iprim,itern,7)))
      tmpsum = tmpsum + (mcon_prim_EWH(iprim) - mcon_tern_EWH(nint(mcon_tern(iprim,itern,7))) )**2
    enddo
   
    write(51,'(i5". MC",i5.5," (m)            ",f17.7,f17.7,f17.7)')iprim,iprim,0.d0,dsqrt(tmpsum/mcon_prim(iprim,8)) &
            ,mcon_prim_EWH(iprim)
  enddo
  close(51)
  call status_update('STATUS','GRACEORB','read_ternary_EWH','mcon_prim_EWH.rms',"Output RMS of ternary EWH values per primary mascon",0)

  return

  end subroutine read_ternary_EWH
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











 
