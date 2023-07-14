   subroutine tidegrid(jd, tin, otidecoef)  

! Emma-Kate Potter, 13 Oct, 2010
! This subroutine calls for and returns the grid of tidal heights to oceanfield 
!
! tide_ti  : epoch of tide grid that precedes the epoch of observation
! tide_tf  : epoch of tide grid that follows  the epoch of observation

    use gm_mod
    use gauleg_mod
    use tide_mod
    use inmod_mod     ! PT130920: added to get the ocean tide model name into this routine

    implicit none

    integer*4 :: i, j
    integer*4 :: jd 
    integer*4, dimension(5) :: date, datei, datef
    integer:: tarray(3)
    integer*4, dimension(5) :: it, iti, itf

    real*8 :: tin, time
    real*8 :: temp 
    real*8, allocatable, dimension(:,:) :: tideoh, tideohi, tideohf
    real*8 :: sec, seconds, epoch, gradient
    real*8, dimension(maxnml) :: otidecoefi            ! tide coefficients that precede epoch of observation
    real*8, dimension(maxnml) :: otidecoeff            ! tide coefficients that follow the epoch of observation
    real*8, dimension(maxnml) :: otidecoef             ! tide coefficients linearly interpolated to epoch of observation

    character tidemod*10                                            ! ocean tide model name (e.g. eot11a, tpxo70)
    double precision, allocatable          :: ocean_ht(:)           ! vector of ocean heights for each tide grid point
    integer*4                              :: counter
    double precision, allocatable          :: lat_fft(:),lon_fft(:) ! vectors for latitudes/longitudes for the tidal grid for the FFT
    double precision, allocatable          :: lat_weight(:)         ! gauleg returns a weight value, so we dimension it here as well
    double precision, allocatable          :: tmp_lat(:),tmp_lon(:)
    character*5 upperc
    logical first_flag, debug_flag

! APP: 130218 - Changed data declaration so that initial and end times for the current 10 
!               minute interval are stored in a common block
!    integer :: counter = 1
!    real*8 :: ti = 0.d0, tf = 0.d0

    save :: tideohi, tideohf 
    save :: datei, datef
    save :: first_flag
    save :: lat_fft, lon_fft,ocean_ht
! PT131202: we have to save the values here for subsequent use. Why did the code never do that????
     save :: otidecoefi,otidecoeff

!  print*,'top of tidegrid tide_counter = ',tide_counter
!  if(tide_counter > 0) then
!    print*,'top of tidegrid tide_counter > 1 in tidegrid'
!      print*,'tidegrid: ocean_ht(3830:3860)=',ocean_ht(3830:3860)
!    print*,'are we still alive?'
!  endif

    pi = 4.d0*datan(1.d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    time = dble(jd)+tin/86400.d0
    call JD_to_YMDHMS(time, date, sec) ! return the YMDHMS for this JD
    call ymdhms_to_ydoyhms(date, sec, it)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    tide_counter = tide_counter+1
! size of grid
    nlat_fft = 256
    nlon_fft = 512
    if(tide_counter == 1)then
      allocate(tideoh(nlon_fft,nlat_fft))
      allocate(tideohi(nlon_fft,nlat_fft))
      allocate(tideohf(nlon_fft,nlat_fft))
      allocate(ocean_ht(nlat_fft*nlon_fft))
      allocate( tmp_lon(nlon_fft))
      allocate( tmp_lat(nlat_fft))
      allocate( lat_fft(nlat_fft*nlon_fft))
      allocate( lon_fft(nlat_fft*nlon_fft))
      allocate( lat_weight(nlat_fft))

! for the latitudes, we need to find the abscissas for degree 256
      tmp_lat = 0.d0
      lat_weight = 0.d0
      call gauleg(-1.d0,1.d0,tmp_lat,lat_weight,nlat_fft)    ! returns dcos(colat)

      do i=1,nlat_fft
        tmp_lat(i) = 90.d0-dacos(tmp_lat(i))*180.d0/pi
        do j=1,nlon_fft
          lat_fft((i-1)*nlon_fft + j) = tmp_lat(i)
          lon_fft((i-1)*nlon_fft + j) = dble(j-1)*360.d0/dble(nlon_fft)
        enddo
      enddo
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE the initial and final ocean tide height coefficients for this 10 min interval
    if (tide_counter.eq.1) then
      date(5) = 10 * int(date(5)/10)         ! find the closest (lower) 10 min timepoint
      sec = 0.d0
      call YMDHMS_to_JD(date, sec, epoch)    ! determine the jd for this 10min timepoint
      call ymdhms_to_ydoyhms(date, sec, iti) ! determine ydoyhms for this 10min timepoint
      tide_ti = epoch

! get the right name for the tide model
      if(upperc(gt_oceantidemod) == "FES95")then
        tidemod = 'fes952    '
      else if (upperc(gt_oceantidemod) == "EOT11")then
        tidemod = 'eot11a    '
      else if (upperc(gt_oceantidemod) == "TPX70")then
        tidemod = 'tpxo70    '
      else
        write(message,'(a,a,a)')"Ocean tide model (",gt_oceantidemod,") not coded. Do it yourself please :-)"
        call status_update('FATAL','GRACEORB','tidegrid',' ',message,0)
      endif
  
! interact with the spotl routines to get the amplitudes and phases
      first_flag = .true.
      debug_flag = .false.
!  print*,'calling spotl 1st time iti=',iti
      call spotl_get_ocean_height(nlat_fft*nlon_fft,first_flag,lat_fft,lon_fft,iti,tidemod,ocean_ht,debug_flag)
! DEBUG
!  call spotl_get_ocean_height(131001,first_flag,lat_fft(1:131000),lon_fft(1:131000),iti,tidemod,ocean_ht(1:131000),debug_flag)
!  do i=1,131000 
!    print*,'tidegrid 131000 points: ',i,lat_fft(i),lon_fft(i),ocean_ht(i)
!  enddo
!  stop    
      counter = 0
      do i=1,nlat
        do j=1,nlon
          counter = counter + 1
          tideohi(j,i) = ocean_ht(counter)
!  print*,'tideohi lat,lon, ocean_ht',lat_fft((i-1)*nlon + j),lon_fft((i-1)*nlon + j),tideohi(j,i),ocean_ht(counter),nlat,nlon,i,j
        enddo
      enddo
         
      call sphharm(tideohi,otidecoefi)       ! calculate C&S coeff for this timepoint
  ! second grid
      tide_tf = tide_ti+10.d0/60.d0/24.d0    ! adds 10 min to find next 10 min timepoint             
      call JD_to_YMDHMS(tide_tf, date, sec)  ! return the YMDHMS for this JD
      call ymdhms_to_ydoyhms(date, sec, itf)
!  print*,'calling spotl 2nd time itf=',itf
      call spotl_get_ocean_height(nlat_fft*nlon_fft,.false.,lat_fft,lon_fft,itf,tidemod,ocean_ht,debug_flag)
      counter = 0
      do i=1,nlat
        do j=1,nlon
          counter = counter + 1
          tideohf(j,i) = ocean_ht(counter)
! DEBUG
!  print*,'tideohf lat,lon, ocean_ht',lat_fft((i-1)*nlon + j),lon_fft((i-1)*nlon + j),tideohf(j,i),ocean_ht(counter),nlat,nlon,i,j
        enddo
      enddo
!  stop 'stopped after second time'
! DEBUG: print out the entire tide grid in format for analyse_hs
!      do i=nlat,1,-1
!        print*,(tideohf(j,i),j=1+nlon/2,nlon),(tideohf(j,i),j=1,nlon/2)," tide grid for analyse_hs"
!      enddo
!      stop ' stopped after printing out tide grid'
      call sphharm(tideohf,otidecoeff)       ! calculate C&S coeff for this timepoint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  subsequent calls. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else ! if counter>1
      debug_flag = .false.

      if (time.gt.tide_tf) then
        tide_ti = tide_tf
        tideohi = tideohf
        otidecoefi = otidecoeff               !  transfers the second set of coefficients to being the first set.
        sec = 0.0d0
        tide_tf = tide_ti+10.d0/60.d0/24.d0   ! adds 10 minutes to get the time of the new second grid
        call JD_to_YMDHMS(tide_tf, date, sec) ! return the YMDHMS for this JD
        call ymdhms_to_ydoyhms(date, sec, itf)
        call spotl_get_ocean_height(nlat_fft*nlon_fft,first_flag,lat_fft,lon_fft,itf,tidemod,ocean_ht,debug_flag)
        counter = 0
        do i=1,nlat
          do j=1,nlon
            counter = counter + 1
            tideohf(j,i) = ocean_ht(counter)
! DEBUG
!  print*,'tideohf lat,lon, ocean_ht',lat_fft((i-1)*nlon_fft + j),lon_fft((i-1)*nlon_fft + j),tideohf(j,i),nlat,nlon,i,j
          enddo
        enddo
        call sphharm(tideohf,otidecoeff)      ! calculate C&S coeff for this timepoint
      endif 
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INTERPOLATE (linear) for each grid point to determine ocean height for current epoch
     do i = 1,maxnml 
       gradient = (otidecoeff(i)-otidecoefi(i))/(tide_tf-tide_ti)
       otidecoef(i) = otidecoefi(i)+gradient*(time-tide_ti)
!       if(i < 100) print*,'i,otidecoefi(i),otidecoeff(i),otidecoef(i)',i,otidecoefi(i),otidecoeff(i),otidecoef(i),tin &
!                   , tide_ti,time-tide_ti
     enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   return
   end
