!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Subroutines related to tides (SPOTL stuff moved out of here to an unused subroutine spotl_tide_lib.f90)
!
! write_tide_head  :   write a header on the binary file that contains the ocean tide heights (every 10 mins) for a 24-hr period
! which_tide_mascon :  return which secondary tide mascon a particular XYZ location is in
! compute_sph_tide_acc : computes the accelerations at a satellite, given a set of spherical harmonic 
!                        coefficients and a satellite location
! primary_tide_vol  : calculate an integrated ocean height from the secondary tide heights within a primary mascon
! check_ternary_tides :  determine which ternary mascons are oceanic but within a secondary mascon of both oceanic and continental ternary mascons
! setup_netcdf_tides  : subroutine to open and read header information from the netcdf tide grid file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!************************!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_tide_grid_head(luout,num_msc,dlat,tidemod,iyr,imonth,iday,interval,units,irec,mcon_lat,mcon_lon)

! subroutine to write a tide grid file binary header
!
! P. Tregoning
! 11 September 2013

  implicit none

  integer*4,   intent(in)    ::  luout              ! unit number of binary output file
  integer*4,   intent(in)    ::  num_msc            ! number of secondary mascon tide elements
  real*4,      intent(in)    ::  dlat               ! latitudinal spacing of mascon bands
  character*10 ,intent(in)   ::  tidemod            ! ocean tide model being written out
  integer*4,   intent(in)    :: iyr,imonth,iday     ! year, month, day of the ocean tide height data
  real*4,      intent(in)    :: interval            ! grid time interval within file
  character*4, intent(in)    :: units               ! units of grid time interval
  real*4   , intent(in)      :: mcon_lat(num_msc)   ! latitudes  of all the mascons
  real*4   , intent(in)      :: mcon_lon(num_msc)   ! longitudes of all the mascons

  integer*4, intent(out)     ::  irec               ! record number of last line written out
 
  character*4 tidemod_short                         ! 4-byte version of tide model name
  integer*4 trimlen
  character*200 message
  integer*4  i

! first line of the file is the version number of grid file
  irec = 1
  write(luout,rec=irec)"1.00"

! date of the 24 hours of tidal information. year, month, day - each on a separate line
  irec=irec+1
  write(luout,rec=irec)iyr
  irec=irec+1
  write(luout,rec=irec)imonth
  irec=irec+1
  write(luout,rec=irec)iday

! ocean tide model
  irec=irec+1
  if(tidemod(1:6) == "eot11a")then
    tidemod_short = "et11"
  else if (tidemod(1:3) == "got")then
    tidemod_short = "got4"
  else if (tidemod(1:3) == "fes")then
    tidemod_short = "fs12"
  else
    write(message,'(a,a,a)')'tide model (',tidemod(1:trimlen(tidemod)),') not coded'
    call status_update('FATAL','LIB','write_tide_grid_head',' ',message,0)
  endif
  write(luout,rec=irec)tidemod_short

! number of mascon tide elements per epoch
  irec=irec+1
  write(luout,rec=irec)num_msc

! epoch interval of each grid entry in the file
  irec=irec+1
  write(luout,rec=irec)interval

! units of grid time interval
  irec=irec+1
  write(luout,rec=irec)units

! latitude spacing of each band of elements (in minutes)        
  irec=irec+1
  write(luout,rec=irec)dlat

! PT150806: to be compatible with the fes2012 grids, we must now include the latitude coords of all elements
  do i=1,num_msc
    irec=irec+1
    write(luout,rec=irec)mcon_lat(i)
  enddo

! PT150806: to be compatible with the fes2012 grids, we must now include the longitude coords of all elements
  do i=1,num_msc
    irec=irec+1
    write(luout,rec=irec)mcon_lon(i)
  enddo
!  print*,'last header line was record',irec
!stop

  return
  end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!************************!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_tide_grid_head(luin,num_msc,dlat,tidemod,iyr,imonth,iday,interval,units,irec,lat,lon)

! subroutine to read a tide grid file binary header
!
! P. Tregoning
! 11 September 2013

  implicit none

  integer*4,   intent(in)     ::  luin               ! unit number of binary output file
  integer*4,   intent(out)    ::  num_msc            ! number of secondary mascon tide elements
  real*4,      intent(out)    ::  dlat               ! latitudinal spacing of mascon bands
  character*10 ,intent(out)   ::  tidemod            ! ocean tide model being read in
  integer*4,   intent(out)    ::  iyr,imonth,iday    ! year, month, day of the ocean tide height data
  real*4,      intent(out)    ::  interval           ! grid time interval within file
  character*4, intent(out)    ::  units              ! units of grid time interval
  integer*4, intent(out)      ::  irec               ! record number of last line of header
  real*4   , intent(out)      ::  lat(200000),lon(200000)  ! lat/lon values of each grid node (hardwired declaration > 164285 secondary mscs)

  character*4 tidemod_short                         ! 4-byte version of tide model name
  integer*4 trimlen
  character*200 message
  character*4 versn                                 ! version number of the ocean tide height grid file

  integer*4 i

! first line of the file is the version number of grid file
  irec = 1
  read(luin,rec=irec)versn

! date of the 24 hours of tidal information. year, month, day - each on a separate line
  irec=irec+1
  read(luin,rec=irec)iyr
  irec=irec+1
  read(luin,rec=irec)imonth
  irec=irec+1
  read(luin,rec=irec)iday

! ocean tide model
  irec=irec+1
  read(luin,rec=irec)tidemod_short
  if(tidemod_short == "et11") then
    tidemod(1:6) = "eot11a"
    tidemod_short = "et11"
    write(message,'(a,a4,a,3i4,a,a)')'Grid version: ',versn,'   Epoch: ',iyr,imonth,iday,'   Tide model: ',tidemod(1:6)
    call status_update('STATUS','LIB','read_tide_grid_head',' ',message,0)
  else if (tidemod_short == "got4")then
    tidemod(1:6) = "got4p7"
    write(message,'(a,a4,a,3i4,a,a)')'Grid version: ',versn,'   Epoch: ',iyr,imonth,iday,'   Tide model: ',tidemod(1:6)
    call status_update('STATUS','LIB','read_tide_grid_head',' ',message,0)
  else if (tidemod_short == "fs12")then
    tidemod(1:6) = "fes12"
    write(message,'(a,a4,a,3i4,a,a)')'Grid version: ',versn,'   Epoch: ',iyr,imonth,iday,'   Tide model: ',tidemod(1:6)
    call status_update('STATUS','LIB','read_tide_grid_head',' ',message,0)
  else
    write(message,'(a,a,a)')'tide model found in grid file (',tidemod_short,') not coded'
    call status_update('FATAL','LIB','read_tide_grid_head',' ',message,0)
  endif

! number of mascon tide elements per epoch
  irec=irec+1
  num_msc = 0
  read(luin,rec=irec)num_msc

! epoch interval of each grid entry in the file
  irec=irec+1
  read(luin,rec=irec)interval

! units of grid time interval
  irec=irec+1
  read(luin,rec=irec)units

! latitude spacing of each band of elements (in minutes)        
  irec=irec+1
  read(luin,rec=irec)dlat

! PT140603: read the latitude coords of all elements
  do i=1,num_msc
    irec=irec+1
    read(luin,rec=irec)lat(i)
  enddo

! PT140603: read the longitude coords of all elements
  do i=1,num_msc
    irec=irec+1
    read(luin,rec=irec)lon(i)
!print*,'i,lat(i),lon(i)',i,lat(i),lon(i)
  enddo
  write(message,'(a,i7,a,f6.2,a1,a4,a,f6.2,a)')'Number of points: ',num_msc,'   Temporal interval: ',interval," ",units &
                                             ,'   Lat spacing: ',dlat,"'"
  call status_update('STATUS','LIB','read_tide_grid_head',' ',message,0)

! increment for the blank line (for who knows what later on .....)
!  irec=irec+1
!  print*,'last header line was record',irec
!stop

  return
  end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine integrate_msc_accel(satxyz,msc_XYZ,msc_area,msc_rho,msc_ht,acc,partial)

! subroutine to compute the XYZ accelerations caused by a height of water on a single mascon 
!
! P. Tregoning
! 27 November 2013

  implicit none

  real(kind=8), intent(in)   :: satxyz(3)        ! cartesian coordinates of satellite location
  real(kind=8), intent(in)   :: msc_XYZ(3)       ! mascon cartesian coordinates
  real(kind=8), intent(in)   :: msc_area         ! mascon area    (m^2)
  real(kind=8), intent(in)   :: msc_rho          ! mascon density (kg/m^3)
  real(kind=4), intent(in)   :: msc_ht           ! mascon equivalent water height (m)
  real(kind=8), intent(out)  :: acc(3)           ! XYZ accelerations acting on the satellite (m/s^2)
  real(kind=8), intent(out)  :: partial(3)       ! XYZ acceleration partials

! local variables
  integer*4     :: j
  real(kind=8)  :: tmp_msc_acc(3),d_crd(3),msc_mass,msc_acc_fact,msc_partial_fact,dist
  real(kind=8)  :: Gconst

  Gconst = 6.674e-11 ! Gravitational constant

! compute the distance between the satellite and the mascon
  d_crd =  msc_XYZ - satxyz
  dist = dsqrt( d_crd(1)**2 + d_crd(2)**2 + d_crd(3)**2 )

! compute various factors
  msc_mass         = msc_rho * msc_area * msc_ht
  msc_acc_fact     = Gconst * msc_mass/dist**3
  msc_partial_fact = Gconst * msc_rho * msc_area / dist**3

  tmp_msc_acc = 0.d0
  acc(j) = 0.d0
  do j=1,3
! compute accelerations
    tmp_msc_acc(j) = d_crd(j) * msc_acc_fact
! PT131127: these are insignificant (1e-14), so just ignore them here ....
!    tmp_mcon_acc(j) = tmp_mcon_acc(j) + msc_mass *   up_def_vector(j,m)
!    tmp_mcon_acc(j) = tmp_mcon_acc(j) + msc_mass * tang_def_vector(j,m)
    acc(j) = acc(j) + tmp_msc_acc(j)
! partial derivative
    partial(j) = partial(j) + d_crd(j) * msc_partial_fact
! PT131127: these are insignificant (1e-14), so just ignore them here ....
!    partial(j) = mcon_partial(j) + mcon_rho_area *   up_def_vector(j,m)
!    partial(j) = mcon_partial(j) + mcon_rho_area * tang_def_vector(j,m)
  enddo


! that's it, we're done !!
  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine which_tide_mascon(xyz,point_num)

! subroutine to compute which tide (i.e. secondary) mascon a particular XYZ location is in.
!
! Coordinates of ocean tide mascons (ie secondary mascon grid) are passed in through the common in tide_mod.f90 (otide_lat, otide_lon) 
!
! P. Tregoning
! 5 June 2014

  use tide_mod      ! needed for coords of tide mascons          (otide_lat, otide_lon)
  use accel_mod     ! needed for the number of secondary mascons (num_mcon2)

  implicit none

  real(kind=8),   intent(in )  :: xyz(3)
  integer*4,      intent(out)  :: point_num

  real(kind=8)                 :: tol_dlat                ! XYZ needs to be within tol_dlat of a mascon latitude 
  real(kind=8)                 :: tmplat,tmplon
  real(kind=8)                 :: pi
  integer*4                    :: i,j
  real(kind=8)                 :: dist,tmp_tmp_lat,tmp_tmplon(1000),dlon

! debug variables
  logical debug



  pi = 4.d0*datan(1.d0)

! compute the lat/lon from the input XYZ
  dist = dsqrt(xyz(1)**2+xyz(2)**2+xyz(3)**2)
  tmplat = dasin(xyz(3)/dist)*180.d0/pi
  if(dabs(xyz(1)) < 0.01d0 .and. dabs(xyz(2)) < 0.01d0)then
    if(xyz(3) > 0)point_num = 1
    if(xyz(3) < 0)point_num = 164838
    return
  else if(xyz(1) /= 0.d0)then
    tmplon = datan2(xyz(2),xyz(1))*180.d0/pi
  else
    if(xyz(2) >= 0.d0) tmplon = 90.d0
    if(xyz(2) <  0.d0) tmplon = 270.d0
  endif

  if(tmplon < 0 ) tmplon = tmplon + 360.d0

!print*,'which_tide_mascon tmplat,tmplon',tmplat,tmplon
  tol_dlat = 0.25d0
! loop through the latitudes until we find the right band of the requested point
  i = 1
! PT160919: doesn't the sign here need to be ">" rather than "<" ?
!           NO. The otide_lat are ordered from -90 to +90 in the tide grid file, so it should be "<"
  do while (tmplat-otide_lat(i) < tol_dlat .and. i < 164838) !.and. dabs(tmplat-otide_lat(i)) < tol_dlat)
!    print*,'i,otide_lat(i),tmplat,dtol_lat',i,otide_lat(i),tmplat,tmplat-otide_lat(i),tol_dlat
    i = i + 1
    if(dabs(tmplat-otide_lat(i)) <= tol_dlat)then
      tmp_tmp_lat = otide_lat(i)
! make a list of all the longitudes with this latitude
      j=1
      dlon = 1.d8
      do while (otide_lat(i) ==  tmp_tmp_lat)
        tmp_tmplon(j) =  otide_lon(i)
! check whether the point is closer than earlier ones. If so, update closest point
        if(dabs(tmp_tmplon(j) - tmplon) < dlon)then  
          dlon = dabs(tmp_tmplon(j) - tmplon)
          point_num = i
        endif
        j = j+1
        i = i+1  
      enddo
    endif        
  enddo

! one final case: if i = num_msc then the returned secondary mascon is the last one
  if (i == 164838) point_num = i

  return

  end  subroutine which_tide_mascon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine compute_sph_acc(rlat,rlong,rx,C,S,maxdeg,Plm,PDlm,flag,acc)

! subroutine to compute the accelerations acting on a satellite, given a set of spherical harmonic coefficients
! and a satellite location. The coefficients need to be Stokes' coefficients for this to work, so that it is a
! generic subroutine that can be applied to any of the components (dealiasing, tides, static gravity etc)
! The code is based on the lower part of sphharmfield.f90
!
! P. Tregoning
! 29 July 2016

  implicit none

  real(kind=8), intent(in)  :: rlat,rlong,rx       ! satellite location (in radians and metres)
  integer*4,    intent(in)  :: maxdeg             ! dimension of the spherical harmonic matrices (actually set in coeff_mod.f90)
  real(kind=8), intent(in)  :: C(0:maxdeg,0:maxdeg),S(0:maxdeg,0:maxdeg)  ! the spherical harmonic coefficients
  real(kind=8), intent(out) :: acc(3)              ! the computed accelerations
  character*6, intent(in)  :: flag               ! flag to indicate scaling factor to convert the spherical harmonic coefficients to different units (if needed)
  

! local variables
  real(kind=8)  ::rlong_tmp       ! a temporary longitude variable so that we can convert it to -180/180 range
  real(kind=8)  :: Plm((maxdeg+1)*(maxdeg+2)/2),Pdlm((maxdeg+1)*(maxdeg+2)/2) ! associated Legendre functions and derivatives
  real(kind=8)  :: sumrt,sumlatt,sumlont,multl,Ae=6378136.46,sumlat,sumlon,sumr
  real(kind=8)  :: pi,grav(3),g_cart(3)
  real(kind=8)  :: gm                                        ! we need GM of the Earth
  integer*4     :: nm,ideg,iord
  real(kind=8)  :: ocean_ht
  real(kind=8)  :: rho_w=1000., rho_av = 5515.d0
  logical       :: factor_flg
  real(kind=8)  :: tmpfactor
  character*100 :: message

! define some values
  gm = 398600.4418D9
  pi = 4.d0*atan(1.d0)


! initialise the summation variables to zero
  sumrt = 0.d0
  sumlatt = 0.d0
  sumlont = 0.d0 
  ocean_ht = 0.d0

! calculate the accelerations
  nm = 1         ! start at degree 1
  do ideg = 0,  maxdeg 
    multl = (Ae/rx)**ideg

    if(flag == "EWH   ") then
      tmpfactor = 3.*rho_w/(Ae*rho_av)/(2.*dble(ideg)+1.)
    else if (flag == "press ")then
      tmpfactor = 1.0 ! code this correctly later when Iknow how
    else if (flag == "stokes")then
      tmpfactor = 1.0
    else
      write(message,'(a,a,a)')"Unknown spherical harmonic type: ",flag," Code it yourself."
      call status_update('FATAL','GRACEORB','compute_sph_acc',' ',message,0)
    endif
!    factor = 3.*rho_w/(Ae*rho_av)*1./(2.*dble(ideg)+1.)

    sumrt = sumrt+Plm(nm)*C(ideg,0)*dble(ideg+1)*multl*tmpfactor   ! * cos(m*lambda) 
    sumlatt = sumlatt+dcos(rlat)*PDlm(nm)*C(ideg,0)*multl*tmpfactor
    sumlont = sumlont+0.d0
    nm = nm + 1

   do iord = 1, ideg    
      sumr =  Plm(nm)*multl*dble(ideg+1)*(C(ideg,iord)* &
              dcos(iord*rlong)+S(ideg,iord)*dsin(iord*rlong) )*tmpfactor
      sumlat = dcos(rlat)*( PDlm(nm)*multl*(C(ideg,iord)* &
               dcos(iord*rlong)+S(ideg,iord)*dsin(iord*rlong) ))*tmpfactor
      sumlon = (dble(iord)*Plm(nm)*multl*(-C(ideg,iord)* &
               dsin(iord*rlong)+S(ideg,iord)*dcos(iord*rlong) ))/dcos(rlat)*tmpfactor
      sumrt = sumrt + sumr 
      sumlatt = sumlatt + sumlat 
      sumlont = sumlont + sumlon 

      nm = nm + 1
    enddo
  enddo

! SCALE to give gravitational potential at earth's surface by multiplying 
! by GM/R values of gravitational constants and earth radius 
! grav() array gives r, theta, phi elements of gravity vector in spherical coordinates         
  grav(1) = -gm/((rx)**2)*(sumrt) !  PT160928: exclude the "+1" when calculating tidal accelerations    +1.d0) 
  grav(2) = gm/((rx)**2)*(sumlatt) 
  grav(3) = gm/((rx)**2)*(sumlont) 

! set longitude range to -180/180
  rlong_tmp = rlong
  if (rlong.gt.pi)rlong_tmp = rlong-2.d0*pi

! convert calculated gravity into cartesian coordinates:
  call cart_convert(rlat,rlong_tmp,grav,acc)

  return
  end subroutine compute_sph_acc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine primary_tide_vol(num_prim,t,dt_tides,ocean_ht)

! subroutine to calculate the integrated tidal volume over a primary mascon by using the secondary mascons
!
! P. Tregoning
! 9 August 2016

  use tide_mod   ! added for binary ocean tide height grid variables
  use accel_mod  ! added for primary and secondary mascon arrays
  
  implicit none

  integer*4 , intent(in)   :: num_prim              ! primary mascon number
  real(kind=4), intent(in) ::  t,dt_tides           ! time values for tide interpolation and time argument for periodic tides
  real(kind=8), intent(out):: ocean_ht              ! the integrated tide volume over the ocean part of the primary mascon

! local variables
  integer*4  :: isecondary, point_num, num_sec
  real(kind=8) :: primary_volume,tmp_ht ,tmp_area

  primary_volume = 0.0
  tmp_area = 0.

! for this primary mascon, loop over all the secondary mascons
  num_sec = num_mcon2(num_prim)
  do isecondary = 1, num_sec
! get the point number for the secondary mascon 
    point_num = mcon_sec2sec(isecondary,num_prim)

! get the tide height for this point for this epoch
    tmp_ht = dh1(point_num) + (dh2(point_num)-dh1(point_num)) * mod(t,dt_ocean_grid)/dt_ocean_grid
    tmp_ht = tmp_ht / 1.d2  ! (convert from FES2012 cm to metres)
    tmp_area = tmp_area + mcon_area2(isecondary,num_prim)

! if the secondary mascon is part ocean, part land then ..... what do we do?
    if(tmp_ht == 0. .and. mcon_fract_ocean2(mcon_sec2sec(isecondary,num_prim)) > 0.)then
      print*,"***** Problem: secondary mascon tide height is zero but there is ocean area ...",num_prim,point_num,tmp_ht &
         , mcon_fract_ocean2(mcon_sec2sec(isecondary,num_prim))
    endif
! multiply by the area times the fraction of the secondary mascon that is ocean
    primary_volume = primary_volume + tmp_ht*mcon_area2(isecondary,num_prim)*mcon_fract_ocean2(mcon_sec2sec(isecondary,num_prim))

! DEBUG
!if(num_prim == 18)then
!   print*,'primary #: ',num_prim,' secondary #',isecondary,' tide height: ',tmp_ht &
!          ,dh1(point_num)/100.,dh2(point_num)/100.,point_num &
!          ,' % ocean: ',mcon_fract_ocean2(mcon_sec2sec(isecondary,num_prim)),mcon_sec2sec(isecondary,num_prim)
!endif
  enddo

! now, convert back to just a height spread over the entire primary mascon
  ocean_ht = primary_volume / mcon_area1(num_prim)

! DEBUG
!  print*,num_prim,primary_volume,ocean_ht,mcon_area1(num_prim),tmp_area,' mascon, tide volume, ocean ht'

  end subroutine primary_tide_vol
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine check_ternary_tides(tin)

! subroutine to work out which are the ocean ternary mascons within a secondary mascon of mixed ocean/continent
! ternary mascons
!
! P. Tregoning
! 10 August 2016

  use accel_mod          ! provides the matrices of the primary, secondary and ternary mascons

  implicit none

  real*8, intent(in)  :: tin         ! time of the start of the integration

! local variables
  integer*4 iprim,isec,itern,n_coastal_tern,n_tern
  real*8 :: pi, jd,jd_secs
  integer*4 :: date(5)
  character message*200


print*,'inside check_ternary_tides'
print*,'just returning without doing anything'
return

  pi = 4.d0*atan(1.d0)
  n_coastal_tern = 0
  n_tern = 0

! convert the start time of integration (in GRACE seconds) to yymmdd
  jd = 51544.5d0 + (tin + 0.05)/86400.d0
!print*,'calling jd_to_ymdhms'
  call jd_to_ymdhms(jd, date, jd_secs)

! loop through all the primary mascons
  do iprim = 1, num_mcon
  write(message,'(a,i6,a,i8)')'primary number ',iprim, '  coastal ternary number ',n_coastal_tern
  call status_update('STATUS','TIDE_LIB','check_ternary_tides',' ',message,0)

! loop through all ternary mascons for this primary
    do itern = 1,num_mcon3(iprim)
      n_tern = n_tern + 1

! is the ternary oceanic AND in a secondary mascon comprised of ocean and continent ternary mascons?
      isec = mcon_tern2sec(itern,iprim)
      if(mcon_fract_ocean2(isec) /= 0. .and. mcon_fract_ocean2(isec) /= 1.0)then           ! it is a mix of oceanic and continental secondary mascons
        if(mcon_rho3(itern,iprim) > 1000.) then            ! it is an oceanic ternary mascon
          n_coastal_tern = n_coastal_tern + 1              ! increment the number of coastal ternary mascons
          coastal_ternary(n_coastal_tern,1,1) = n_tern     ! store the pointer to which ternary mascon

!          print*,iprim,isec,itern,' is a coastal ternary mascon in a mixed secondary mascon',n_coastal_tern  &
!               , 90. - mcon_colat3(itern,iprim) *180./pi, mcon_lng3(itern,iprim)*180./pi

! now calculate the tide heights on this particular ternary mascon for each 10' tide epoch, from 0h to 24h
          call calc_ternary_tide(n_coastal_tern,90.-mcon_colat3(itern,iprim)*180./pi,mcon_lng3(itern,iprim)*180./pi &
                                 ,date(1), date(2), date(3),145, coastal_ternary(n_coastal_tern,2,:))
        endif
      endif
    enddo
  enddo

  stop

  end subroutine check_ternary_tides
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calc_ternary_tide(n_coastal_tern,lat,lon,yy,mm,dd,nvals,ternary_tide)

! subroutine to calculate the 10' tide heights on a particular ternary mascon and return it in a vector
!
! P. Tregoning
! 10 August 2016

  implicit none

  integer*4   , intent(in)  :: n_coastal_tern  ! the number of the coastal ternary being evaluated (for debug)
  integer*4   , intent(in)  :: yy,mm,dd        ! year, month,day for which we want to compute the tide heights
  real(kind=8), intent(in)  :: lat,lon         ! coordinates of the ternary mascon
  integer*4,    intent(in)  :: nvals           ! number of tidal values to compute (should always be 145 for 24 hours)
  real(kind=8), intent(out) :: ternary_tide(nvals)  ! computed tide height every 10' from 0h to 24h
  character :: message*200
  integer*4 :: date(5),i
  real*8    :: secs,mjd


! for now, just run the fes_slev code, write results to a file and then read it in from the file
! we need to convert YMD into CNES MJD (JD0 is 1/1/1950 = 33282)
  date(1) = yy
  date(2) = mm
  date(3) = dd
  date(4:5) = 0
  secs = 0.
  call ymdhms_to_mjd(date,secs,mjd)
  mjd = mjd - 33282     ! converts it to CNES JD 

  write (message,*)"fes_slev_10min ",mjd,lat,lon,date(1:3)," > tide.tmp"

  call system(message)

! open and read in the tide heights (in cm)
  open(513,file="tide.tmp",status='old')
  read(513,'(a)')message ! read a header line
  do i=1,145
    read(513,'(a68,f10.3)',end=1000)message(1:68),ternary_tide(i)
  enddo
  close(513)
  return

1000 continue
! if we are still here it is because the tide program thinks that there is no tide for this ternary mascon. Set the amplitudes to zero and then return
  ternary_tide(:) = 0.
  print*,'coastal ternary mascon',n_coastal_tern,' has no tide', lat,lon
  close(513)
  return

  end subroutine calc_ternary_tide
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine setup_netcdf_tides(tide_grid_file,tin)

! subroutine to open and read the header information of the netcdf tide grid file.
! 
! P. Tregoning
! 13 September 2016

  use netcdf              ! generic netcdf stuff
  use tide_netcdf_mod     ! definitions/functions for interfacing with the tide grid file

  implicit none

  character*150, intent(in) :: tide_grid_file        ! name of the tide grid netcdf file
  real(kind=8),intent(in)   :: tin                   ! start time of integration (seconds since 1 Jan 2000 1200UT)

! local variables
  character*250 :: message
  real(kind=8)  :: jdate,tstart
  integer*4     :: jd,i,n_mascons
  real(kind=8)  :: tide_range(2)

! first, open the file
  tidedata = T_open(tide_grid_file)

! now get the tide file header info
! need to update to new format of Seb's subroutine
!  call T_get_info(tidedata)

! get the number of mascons
  n_mascons = T_get_nmascons(tidedata)

! now allocate the tide height array
  allocate(ternary_ht1(n_mascons))
  allocate(ternary_ht2(n_mascons))
  allocate(ternary_ht_current(n_mascons))

  write(message,'(a,i8,a)')"There are ",n_mascons," mascons found in tide grid file"
  call status_update('STATUS','GRACEORB','setup_netcdf_tides',tide_grid_file,message,0)

! PT160913: set some variables that will, eventually, be settable from the header information ... 
!           code taken from setup_ocean_tides_v2.f90 
  dt_ocean_grid = 10.d0*60.d0   ! 10 min (temporal) epochs converted into seconds.
! calculate the seconds of day (converting "tin" into julian date, then into seconds of day)
  jdate = tin/86400.d0+2451545.d0   ! actual julian *date* (non-integer)
  jd = aint(jdate)                                 ! julian day (integer)
!! PT140829: fixed bug in getting seconds of calendar day when hour > 12 hrs
  tstart = (jdate-aint(jdate-0.5))*86400.d0 -43200.d0  ! seconds of day for 0000-2400UT

! now calculate the tide epochs required. The first epoch of the day (ie seconds-of-day = 0) is epoch "1".
  tide_epoch1 = int(tstart/dt_ocean_grid)+1
  tide_epoch_old = tide_epoch1

! now, extract out the tide heights for these two epochs
   call T_read(tidedata,tide_epoch1,ternary_ht1)
   call T_read(tidedata,tide_epoch1+1,ternary_ht2)

! convert from FES2012 cm tide heights to metres
  ternary_ht1 = ternary_ht1/100.
  ternary_ht2 = ternary_ht2/100.

  tide_range(1) = minval(ternary_ht1)
  tide_range(2) = maxval(ternary_ht1)
  write(message,'(a,i7,a,2(f8.3," m"))')'Tide extremes for epoch',tide_epoch1,': ',tide_range
  call status_update('STATUS','GRACEORB','setup_netcdf_tides',tide_grid_file,message,0)
  tide_range(1) = minval(ternary_ht2)
  tide_range(2) = maxval(ternary_ht2)
  write(message,'(a,i7,a,2(f8.3," m"))')'Tide extremes for epoch',tide_epoch1+1,': ',tide_range
  call status_update('STATUS','GRACEORB','setup_netcdf_tides',tide_grid_file,message,0)

  return
  end subroutine setup_netcdf_tides
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






