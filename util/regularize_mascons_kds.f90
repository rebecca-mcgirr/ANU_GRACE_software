    Program regularize_mascons

! Program based on:	correl_apr
! Author:	Anthony Purcell
! Date:		November 6, 2012
! Description:	Program to calculate the distance between mascons to determine
!		the apriori correlation matrix weighted by distance.
!               Zero correlation applied between land-based and ocean-based
!               mascons (as determined from input density values).
! Modified:     APP130826 - to read in different latitudinal and longitudinal
!               distance scales with different values over land and ocean
! PT140922: add different variances for water and land. Original default was 1.0 m for all. More appropriate is perhaps 0.1 m for land, 0.05 m for ocean.
! PT140923: store the whole matrix, then invert it in this program to output the constraint matrix required in gracefit/gracesim
!
! PT140923: rather than generating a covariance matrix that then needs to be inverted, just build directly the matrix that needs to be added
!           onto the AtWA inside gracefit/gracesim/addnorm. This is the essential difference between this program and correl_apr
! PT140930: use the hydrological seasonal amplitude as the sigma for land mascons between +-60N/S
! PT141020: restrict the hydrological seasonal amplitude to 100 mm max
! PT161128: renamed make_msc_constraint_nov2014 to regularize_mascons and updated to read new mascon file
! PT170503: add detailed regularisation for Hudson Bay, Baltic Sea, Sumatra-Anderman earthquake region etc
!           (only it seems that I never did this on 3 May 2017!!)
!
! PT170704: changed the regularisation for Hudson Bay.
! PT170705: added different length scales for continents, Antarctica, Oceans and Greenland. This changes the runstring.

    use mascon_mod


    implicit none

    integer max_num_prime

    parameter (max_num_prime = 5000)
      
    real*8 a, a_prime, area, b, b_prime, c, c_prime, cos_a, d_prime, dist
    real*8 scale_land, scale_ocean, scale_Ant, scale_Grn
    real*8 pi, r, r_prime, rad_fact,  scale, scale_factor, sigma_ocean, sigma_land, sigma_Ant, sigma_grn, variance, zero
    real*8 test_1, test_2
    real*8 junk

    real(kind=8), allocatable :: whole_matrix(:,:),msc_const(:,:)
    real(kind=8), allocatable :: lon_prime(:),lat_prime(:),rho(:),density(:)
    real(kind=8), allocatable :: sigma_msc(:)


    integer i, ioerr, j, k, l, m, n, count, num_cells
    integer londeg, lonmin, latdeg, latmin,  ht

! PT161128: unit numbers
    integer*4,parameter :: lumsc_in=21,luout=22

! PT150722: variables for the length scale calculation
    integer*4    :: length_scl_pwr
    real(kind=8) :: average_colat,dist_ratio
    real(kind=8),allocatable :: colat_prime(:)
    logical :: zero_correl

! PT170706: command line option to set off-diagonal elements to zero
    logical :: zero_offdiag

    character mascon_file*200, correl_file*80, arg*50, message*256

    integer*4    :: trimlen


! get the name of the mascon file from the command line
    call getarg(1,mascon_file)

    if ( mascon_file == " " ) then
      write(6,10)
10 format(/,'##############################################################################################' &
   ,/,'                                  Program regularize_mascons                    ' &
   ,/,'Purpose: Produce mascon correlation matrix using different length scales for land/ocean & lat/long' &
   ,//,'Usage: regularize_mascons mascon-file output-file scl(land) scl(ocean) scl(Ant) scl(Grn) ' &
   ,  ' sigma_land, sigma_ocean sigma_Ant sigma_Grn (sigmas in metres) [zero_offdiag] ' &
   ,/,'##############################################################################################')
      stop
    endif

    open(unit=lumsc_in,file=mascon_file,status='old',iostat=ioerr)
    if(ioerr /= 0)then
      write(message,'(a,a,a)')"Error opening input file ",mascon_file,". Does it exist?"
      call status_update('FATAL','regularize_mascons ','regularize_mascons',' ',message,ioerr)
    endif

    call getarg(2,correl_file)
    if ( correl_file == " " ) then
      write(message,'(a)')"Too few command line parameters."
      call status_update('FATAL','regularize_mascons ','regularize_mascons',' ',message,ioerr)
    else
      call status_update('STATUS','UTIL','regularize_mascons',correl_file,"Output msc constraint file",0)
    endif

! read in distance scale factors for latitude & longitude on land and over ocean
    call getarg(3,arg)
    read(arg,*) scale_land
    call getarg(4,arg)
    read(arg,*) scale_ocean
    call getarg(5,arg)
    read(arg,*) scale_Ant
    call getarg(6,arg)
    read(arg,*) scale_Grn
    call getarg(7,arg)
    read(arg,*) sigma_land
    call getarg(8,arg)
    read(arg,*) sigma_ocean
    call getarg(9,arg)
    read(arg,*) sigma_Ant
    call getarg(10,arg)
    read(arg,*) sigma_Grn
    call getarg(11,arg)
    if(arg(1:12) == "zero_offdiag")then
      zero_offdiag = .true.
      call status_update('STATUS','UTIL','regularize_mascons',' ',"Will zero out off-diagonal elements",0)
    else
      zero_offdiag = .false.
      call status_update('STATUS','UTIL','regularize_mascons',' ',"Will write out non-zero off-diagonal elements",0)
    endif

! PT170705: don't worry about testing this anymore
!    test_1 = min(lat_scale_land, lon_scale_land, lat_scale_ocean, lon_scale_ocean)
!    test_2 = max(lat_scale_land, lon_scale_land, lat_scale_ocean, lon_scale_ocean)
!    if ( ( test_1 .lt. 1.0d4 ) .or. ( test_2 .gt. 2.0d7 ) ) then
!      write(message,'(a)')"Length scale parameters should lie between 10,000 and 20,000,000 m."
!      call status_update('FATAL','regularize_mascons ','regularize_mascons',' ',message,0)
!    else
      write(message,'(a,i9,i9,i9,i9,a)') "Scale parameters (lat/lon land, lat/lon ocean, lat/lon Antarctica, lat/lon Greenland: " &
                        ,int(scale_land),int(scale_ocean), int(scale_Ant), int(scale_Grn)," metres"
      call status_update('STATUS','regularize_mascons ','regularize_mascons',' ',message,0)
!    endif
    write(message,'(a,e11.3,a,e11.3,a)')'Sigma between mascons on land: ',sigma_land,' Sigma over oceans: ',sigma_ocean,' m.'
    call status_update('STATUS','UTIL','regularize_mascons',' ',message,0)
      write(message,'(a,e11.3,a,e11.3,a)') "Sigma on Antarctic/Greenland mascons : ",sigma_Ant," /",sigma_Grn," metres"
      call status_update('STATUS','regularize_mascons ','regularize_mascons',' ',message,0)

    pi = 4.0d0 * datan(1.0d0)
    rad_fact = pi/180.0d0
    r_prime = 6.3712d6
    zero = 0.0d0
    variance = 1.0d0

! PT161128: use the new subroutines to read the mascon file
  call read_msc_hdr(lumsc_in,mascon_file,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                                ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
  call allocate_mascon_arrays
  call read_mascon_file(lumsc_in,mascon_file(1:trimlen(mascon_file)))

!!! Dimension the entire matrix
    allocate(whole_matrix(total_prim,total_prim))
    allocate(msc_const(total_prim,total_prim))
    allocate(lat_prime(total_prim))
    allocate(colat_prime(total_prim))
    allocate(lon_prime(total_prim))
    allocate(density(total_prim))
    allocate(rho(total_prim))
    allocate(sigma_msc(total_prim))
    msc_const = 0.d0

! PT161128: transfer the lat/lon values into the lat_prime/lon_prime values used in this program
    do i = 1, total_prim
      lat_prime(i)   = mcon_prim(i,1)/rad_fact
      colat_prime(i) = (90.d0 - lat_prime(i))*rad_fact
      lon_prime(i)   = mcon_prim(i,2)
    enddo
    close (21)

! PT170703: set up the uncertainties, including special cases
    do i=1,total_prim

! general cases
      if(mcon_prim(i,6) > 1000.0)then
          sigma_msc(i) = sigma_ocean
      else 
          sigma_msc(i) = sigma_land
      endif

! PT170703: special case of Hudson Bay. Change the density to treat it as "land" so that the GIA signal 
!           is not damped out as though it is ocean
      if(      mcon_prim(i,1)*180.d0/pi > 51.d0   .and. mcon_prim(i,1)*180.d0/pi < 70.5  &
               .and. mcon_prim(i,2)*180.d0/pi > 265.d0  .and. mcon_prim(i,2)*180.d0/pi < 296.d0)then
        sigma_msc(i) = sigma_land
        mcon_prim(i,6) = 1000.d0
      endif

! PT170709: set the uncertainty for Greenland
      if(lat_prime(i) > 60.0 .and. lon_prime(i)*180.0/pi > 300.0 .and. lon_prime(i)*180.0/pi < 338.0) then   ! we are in Greenland
        sigma_msc(i) = sigma_Grn
      endif

!! PT141020: add a different constraint for Antarctica
! PT170517: only do this for land Antarctic mascons
      if(lat_prime(i) < -60.0 .and. mcon_prim(i,6) < 1020.d0)then
        sigma_msc(i) = sigma_Ant
      endif
    enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! M A I N    L O O P  !!!!!!!
    do k = 1, total_prim
      if(mod(k,500) == 0 )then
        write(message,'(a,i5)')'Correlations for mascon: ',k
        call status_update('STATUS','UTIL','regularize_mascons',' ',message,0)
      endif
      b = colat_prime(k)
      b_prime = dabs(dcos(b))
      do n = 1, total_prim
        if ( mcon_prim(k,6) .ne. mcon_prim(n,6) ) then
          rho(n) = zero
        elseif ( k .eq. n ) then
          rho(n) = 1.d0/sigma_msc(n)**2
        else
          c = colat_prime(n)
          c_prime = dabs(dcos(c))

! check whether we should be using the land or ocean scaling factors
          if ( mcon_prim(n,6) > 1000.d0 ) then
            scale_factor = scale_ocean/scale_ocean
            scale = scale_ocean
            variance = sigma_msc(k)*sigma_msc(n)
            zero_correl = .false.
          else
! PT170705: replace this with a subroutine that determines geographically varying values
!            scale_factor = scale_land/scale_land
!            scale = scale_land
            call geogr_scl_factors(mcon_prim(k,1:2),mcon_prim(n,1:2),scale_land,scale_Ant,scale_Grn,scale_factor,scale,zero_correl)
            variance = sigma_msc(k)*sigma_msc(n)
          endif

          d_prime = max(b_prime, c_prime)
          if ( d_prime .gt. 0.99997d0 ) then

! Treat as a special case the instance where either point lies at a pole)
            dist = r_prime * dabs(b - c)
          else

! we scale the longitudinal distance to reflect the influence of latitude
! and the difference in the latitudinal and longitudinal scale factors
            a_prime = dabs(lon_prime(k) - lon_prime(n))
            if ( a_prime .gt. pi ) a_prime = 2.0d0 * pi - a_prime
! APP/PT 140923: multiply this by sin(colat) to account for the reduction in longitude separation closer to the poles. Use mean co-latitude
            a_prime = a_prime * scale_factor  

! Now we use the cosine rule from spherical trigonometry to obtain the angular
! distance between the two mascons and scale by the radius of the sphere
            cos_a = dcos(b) * dcos(c) + dsin(b) * dsin(c) * dcos(a_prime)
            if(cos_a .gt. 1.d0)cos_a = 1.d0
            if(cos_a .lt. -1.d0)cos_a = -1.d0

            a = dabs(dacos(cos_a))
            dist = a * r_prime

! use the distance between the two mascons to calculate correlation
! enter a zero value if the correlation is below some threshold value

          endif
          length_scl_pwr = 1
          average_colat = (colat_prime(k)+colat_prime(n))/2.0
          dist_ratio = dist/(scale/sin(average_colat))
          dist_ratio = dist/scale
! PT150722: undo the latitude dependence for polar regions
!          if(lat_prime(k) < -60 .or. lat_prime(k) > 70)dist_ratio = dist/scale
          rho(n) = 1.d0/variance * exp(-(dist_ratio**length_scl_pwr))


! PT170718: can we try setting the contribution to zero if the distance exceeds some value?
          ! if(dist > 1500.d3)rho(n) = 0.d0   ! anything greater than 1500 km won't contribute

        endif

!if(k==8)print*,"debug:",k,n,msc_const(k,k),rho(n)
! PT140923: now we simply add the values into the final matrix
        msc_const(k,k) = msc_const(k,k) +  1.d0 * rho(n)
        msc_const(n,n) = msc_const(n,n) +  1.d0 * rho(n)
        msc_const(n,k) = msc_const(n,k) + (-1.d0) * rho(n)

! ************************************************************************************************************************************************
! **********     R E M O V E    U N W A N T E D    C O R R E L A T I O N S      ******************************************************************
! ************************************************************************************************************************************************
        if(k .ne. n) then
! PT150501: limit the land correlations of Greenland to only Greenland (i.e. don't let it correlate with Canada)
          if(mcon_prim(k,6) == 1000.0 .and. abs(msc_const(n,k)) > 0.0)then
!!!!!! Greenland !!!!!!!!!
! covariances between a Greenland mascon and Canada
            if(lat_prime(k) > 60.0 .and. lon_prime(k)*180.0/pi > 300.0 .and. lon_prime(k)*180.0/pi < 338.0) then   ! we are in Greenland
! PT150506: set the Greenland sigma to 0.4 m
! PT150723: no, set it to the sigma_polar value
!              print*,'Greenland mascon sigma set to ',sigma_Grn,' m. lat/lon',lat_prime(k),lon_prime(k)/rad_fact
! PT170709: we shouldn't have to set this here - it has been done above
!              msc_const(k,k) = 1.0/(sigma_Grn**2)

              if (lon_prime(n)*180.0/pi < 300.0) then  ! other point is in or west of Canada
!                print*,'k,n, setting Grn/CAN covariance to zero',k,n,msc_const(n,k),lat_prime(k),lon_prime(k)*180.0/pi
                msc_const(n,k) = 0.0
              endif

! covariances between a Greenland mascon and Europe
              if (lon_prime(n)*180.0/pi > 338.0 .or. lon_prime(n)*180.0/pi < 50.0) then
!                print*,'k,n, setting Grn/EUR covariance to zero',k,n,msc_const(n,k),lat_prime(k),lon_prime(k)*180.0/pi
                msc_const(n,k) = 0.0
              endif
            endif
! covariances between a Canadian mascon and Greenland
            if(lat_prime(k) > 60.0 .and. lon_prime(k)*180.0/pi < 300.0 .and. lon_prime(k)*180.0/pi > 200.0) then   ! we are in North America
              if(lat_prime(n) > 60.0 .and. lon_prime(n)*180.0/pi > 300.0 .and. lon_prime(n)*180.0/pi < 338.0) then   ! we are in Greenland
!                print*,'k,n, setting CAN/Grn covariance to zero',k,n,msc_const(n,k),lat_prime(k),lon_prime(k)*180.0/pi &
!                                                                                       ,lat_prime(n),lon_prime(n)*180.0/pi
                msc_const(n,k) = 0.0
              endif
            endif
! covariances between a European mascon and Greenland
            if(lat_prime(k) > 60.0 .and. (lon_prime(k)*180.0/pi > 338.0 .or. lon_prime(k)*180.0/pi < 50.0) ) then   ! we are in Europe
              if(lat_prime(n) > 60.0 .and. lon_prime(n)*180.0/pi > 300.0 .and. lon_prime(n)*180.0/pi < 338.0) then   ! we are in Greenland
!                print*,'k,n, setting EUR/Grn covariance to zero',k,n,msc_const(n,k),lat_prime(k),lon_prime(k)*180.0/pi &
!                                                                                       ,lat_prime(n),lon_prime(n)*180.0/pi
                msc_const(n,k) = 0.0
              endif
            endif


          endif

          if(mcon_prim(k,6) > 1000.0 .and. abs(msc_const(n,k)) > 0.0)then
! PT150501: limit the correlation of the Pacific to the Atlantic across the southern part of South America
            if ( lat_prime(k) < -36.0 .and. lat_prime(k) > -66.0 .and. lon_prime(k)*180.0/pi > 292.0 &
                                                                      .and. lon_prime(k)*180.0/pi < 305.0 ) then   ! point k is in the Atlantic off South America
              if(lat_prime(n) < -36.0 .and. lat_prime(n) > -66.0 .and. lon_prime(n)*180.0/pi < 291.0 &
                                                                      .and. lon_prime(n)*180.0/pi > 260.0) then ! point n is in the Pacific  off South America 
!                print*,'k,n, setting ATL/pac covariance to zero',k,n,msc_const(n,k),lat_prime(k),lon_prime(k)*180.0/pi &
!                                                                                       ,lat_prime(n),lon_prime(n)*180.0/pi
                msc_const(n,k) = 0.0
              endif
            endif
            if(lat_prime(k) < -36.0 .and. lat_prime(k) > -66.0 .and. lon_prime(k)*180.0/pi < 291.0 &
                                                                      .and. lon_prime(k)*180.0/pi > 260.0) then ! point k is in the Pacific  off South America 
              if(lat_prime(n) < -36.0 .and. lat_prime(n) > -66.0 .and. lon_prime(n)*180.0/pi > 292.0 &
                                                                      .and. lon_prime(n)*180.0/pi < 305.0) then   ! point k is in the Atlantic off South America
!                print*,'k,n, setting atl/PAC covariance to zero',k,n,msc_const(n,k),lat_prime(k),lon_prime(k)*180.0/pi &
!                                                                                       ,lat_prime(n),lon_prime(n)*180.0/pi
                msc_const(n,k) = 0.0
              endif
            endif
          endif

! PT170705: if we decided in geogr_scale_factor that we didn't want correlations, set to zero
          if(zero_correl)then
            msc_const(n,k) = 0.d0
          endif


! PT150501: mirror the values into the upper triangular matrix
          msc_const(k,n) =  msc_const(n,k)
        endif

! ************************************************************************************************************************************************

      enddo   ! end of "n" loop
    enddo   ! end of "k" loop

! now write it out
    call status_update('STATUS','UTIL','regularize_mascons',' ',"Writing out mascon constraint matrix",0)
    open(luout,file=correl_file,status='unknown')
! write header lines containing the mascon file code plus info on how this constraint file was created
    write(luout,'(a)')msc_hdr_lines(1)
    write(luout,'(a,4f10.2,a,4e11.3,a)')'# scale dists land, ocean, Ant, Grn:',scale_land/1.d3,scale_ocean/1.d3 &
                  ,scale_Ant/1.d3,scale_Grn/1.d3 &
                  ,' km. Sigmas land/ocean/Ant/Grn: ',sigma_land,sigma_ocean,sigma_Ant,sigma_Grn,' m.'
    do i=1,total_prim
      if(mod(i,500) == 0 )then
        write(message,'(a,i5)')'Writing out correlations for mascon: ',i
        call status_update('STATUS','UTIL','regularize_mascons',' ',message,0)
      endif
!! PT170216: set off-diagonal terms to zero
! PT170706: now control this with a command line option
      if(zero_offdiag)then
        do j=1,total_prim
          if(j /= i)msc_const(i,j) = 0.d0
        enddo
      endif
      write(luout,*)(msc_const(i,j),j=1,total_prim)
    enddo
    close(luout)
    call status_update('STATUS','UTIL','regularize_mascons',correl_file,"Written out mascon constraint matrix",0)

    stop
    end
! ************************************************************************************************************************************************



! ************************************************************************************************************************************************
  subroutine geogr_scl_factors(latlon1,latlon2,scale_land,scale_Ant,scale_Grn,scale_factor,scale,zero_correl)

! subroutine to determine which scale distance to use (land, Ant or Grn), based on the locations of the two sites

  implicit none

  real(kind=8),  intent(in ) :: latlon1(2),latlon2(2)           ! coords (in radians) of the two mascons in question
  real(kind=8),  intent(in ) :: scale_land,scale_Ant,scale_Grn  ! user-defined length scales for continents, Antarctica and Greenland
  real(kind=8),  intent(out) :: scale_factor,scale              ! the appropriate length scale to use for this pair of mascons
  logical     ,  intent(out) :: zero_correl

  real(kind=8) :: pi
  logical :: debug

  pi = 4.d0*datan(1.d0)
  debug = .false.

if(debug)print*,'site coords:',latlon1*180/pi,latlon2*180/pi
! is the first site in Antarctica?
  if(latlon1(1)*180.d0/pi < -60.d0)then   ! yes, it is!
if(debug)print*,'first site in Antarctica'    
  ! is the second site also in Antarctica ?
    if(latlon2(1)*180.d0/pi < -60.d0)then   ! yes, it is!
if(debug)print*,'second site in Antarctica as well'
    ! use the Antarctica length scale
      scale = scale_Ant
      scale_factor = 1.d0
      zero_correl = .false. ! turn off while debugging for Greenland

! PT170718: try treating differently east and west Antarctica. Make the scale factor for West Antarctica half that of East Ant. Also,
!           set the off-diagonal terms to zero if not in the same hemisphere

      if( (latlon1(2)*180.d0/pi < 180.0 .and. latlon2(2)*180.d0/pi .ge. 180.d0)  &
                .or. (latlon1(2)*180.d0/pi > 180.0 .and. latlon2(2)*180.d0/pi .le. 180.d0) ) then  ! in different hemispheres
!         zero_correl = .true.
      else if ( latlon1(2)*180.d0/pi > 180.0 ) then  ! West Antarctica
!         scale = scale / 2.d0
      endif

      return

    else  ! the second site is not in Antarctica. Set the zero  correlation flag so that the off-diagonaly term will be computed as zero
      scale = scale_Ant
      scale_factor = 1.d0  
      zero_correl = .true.
      return
    endif

! is the first site in Greenland ?
  else if (latlon1(1)*180.d0/pi > 60.0 .and. latlon1(2)*180.0/pi > 300.0 .and. latlon1(2)*180.0/pi < 338.0) then  ! yes, we are in Greenland
if(debug)print*,'first site is in greenland'  
  ! is the second site also in Greenland ?
    if (latlon2(1)*180.d0/pi > 60.0 .and. latlon2(2)*180.0/pi > 300.0 .and. latlon2(2)*180.0/pi < 338.0) then  ! yes, we are in Greenland
      scale = scale_Grn
      scale_factor = 1.d0
if(debug)print*,'both sites in greenland',latlon1*180/pi,latlon2*180/pi
      zero_correl = .false.
      return

    else    ! the second site is not in Greenland. Set the scale factor to zero so that the correlation will be computed as zero
      scale = scale_Grn
      scale_factor = 1.d0  ! we can't set this to zero .. it makes the maths wrong in the main program
      zero_correl = .true.
if(debug)print*,'first site in greenland, second not',latlon1*180/pi,latlon2*180/pi
      return
    endif

  else
! we must be on some other continent. Use the scale_land
if(debug)print*,'first site on land somewhere else'
    scale = scale_land
    scale_factor = 1.d0
    ! check whether the second site is in Antarctica or Greenland. If so, set correl to zero
    if(latlon2(1)*180.d0/pi < -60.d0)then 
      zero_correl = .true.
if(debug)print*,'second site in Antarctica'
    elseif (latlon2(1)*180.d0/pi > 60.0 .and. latlon2(2)*180.0/pi > 300.0 .and. latlon2(2)*180.0/pi < 338.0) then
      zero_correl = .true.
if(debug)print*,'second site in Greenland'
    else
if(debug)print*,'second site on land somewhere as well'
      zero_correl = .false.
    endif

    return
  endif

    



  end subroutine geogr_scl_factors





