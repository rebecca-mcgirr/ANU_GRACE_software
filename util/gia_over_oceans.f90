  program gia_over_oceans

! program to calculate the integrated GIA correction to be applied to the integrated ocean EWH change.
! The program will read in the name of the spherical harmonic GIA file, the mascon geometry file 
! plus the GRACE fit/vcv files for which the ! GIA correction is required, then output a file of: 
!
! decima_year   GIA_correction
!
! [I've already written this program once before but can't find it, so I am doing it again ..... ]
!
! P. Tregoning
! 3 August 2020

  use mascon_mod

  implicit none

! command line variables
  character*100 :: GIA_stokes_file        ! spherical harmonic file of stoke's coefficients for the GIA model
  character*150 :: mascon_file            ! mascon geometry file that matches the information used in the GRACE solutions
  character*100 :: grace_solns_file       ! list of GRACE solutions for which GIA corrections are required
  integer*4,parameter :: luGIA=10,lumsc=11,lugrace=12   ! unit numbers for input files

! spherical harmonic variables
  real(kind=8), allocatable :: C_gia(:,:),S_gia(:,:)   ! C/S coefficients of GIA model
  integer*4                 :: maxdeg                  ! maximum degree of GIA spherical harmonic model
  integer*4                 :: ideg,iord               ! counters
  real(kind=8)              :: tmpC,tmpS
  integer*4                 :: tmpdeg,tmpord
  real(kind=8),allocatable  :: Plm(:,:)                ! normalized Legendre functions

! constants
  real(kind=8),parameter :: earthrad=6378136.5d0,p_w=1000.d0,p_av=5515.d0,G = 6.67428e-11

! Load Love number variables
  real(kind=8) :: love_h(0:256),love_l(0:256),love_k(0:256)

! calculation variables
  real(kind=8) :: year, month              ! year and month of input grace file
  real(kind=8) :: dec_yr(1000)             ! array of decimal_year epochs of the GRACE solutions
  integer*4    :: nepochs                  ! number of GRACE epochs
  real(kind=8) :: dt                       ! time since first epoch (used to calculate the GIA rate correction)
  real(kind=8) :: tmpsum                   ! temporary summation of spherical harmonic value
  real(kind=8) :: sum_area,sum_volume      ! calculate area and volume of ocean changes
  real(kind=8) :: GIA_rate                 ! GIA rate integrated over the oceans

! other variables
  character*100 :: line,message
  integer*4     :: ioerr,itrash,i,imsc
  integer*4     :: nmsc_ocean              ! number of ocean mascons
  real(kind=8)  :: colat,rlong,pi,vfact

  pi = 4.d0*datan(1.d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! decode runstring
  call getarg(1,GIA_stokes_file)
  if(GIA_stokes_file(1:1) == "")then
    call status_update('FATAL','UTIL','gia_over_oceans',' ','Runstring: gia_over_oceans GIA.sph mascons_stage5 grace.list',0)
  endif
  call getarg(2,mascon_file)
  call getarg(3,grace_solns_file)


! open the input files
  open(luGIA,file=GIA_stokes_file,status='old',iostat=ioerr)
  if(ioerr /=0)then
    call status_update('FATAL','UTIL','gia_over_oceans',GIA_stokes_file,"Error opening input GIA spherical harmonic file",0)
  else
    call status_update('STATUS','UTIL','gia_over_oceans',GIA_stokes_file,"Opened input GIA spherical harmonic file",0)
  endif

  open(lumsc,file=mascon_file,status='old',iostat=ioerr)
  if(ioerr /=0)then
    call status_update('FATAL','UTIL','gia_over_oceans',mascon_file,"Error opening input mascon file",0)
  else
    call status_update('STATUS','UTIL','gia_over_oceans',mascon_file,"Opened input mascon file",0)
  endif

  open(lugrace,file=grace_solns_file,status='old',iostat=ioerr)
  if(ioerr /=0)then
    call status_update('FATAL','UTIL','gia_over_oceans',grace_solns_file,"Error opening input GRACE solutions list file",0)
  else
    call status_update('STATUS','UTIL','gia_over_oceans',grace_solns_file,"Opened input GRACE solutions list file",0)
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! open the Load Love number file
  open(unit=11,file='Load_Love2_CM.dat',status='old')
! read 14 header lines (including the degree 0 line)
  do i=1,14
    read(11,'(a)')line
  enddo
! read in the coefficients
  do i=1,256
    read(11,*)itrash,love_h(i),love_l(i),love_k(i)
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the spherical harmonic coefficients
  maxdeg = 0
  ioerr = 0
  do while (ioerr == 0)
    read(luGIA,*,iostat=ioerr,end=1000)tmpdeg,tmpord,tmpC,tmpS
    if(ioerr == 0 .and.tmpdeg > maxdeg)maxdeg = tmpdeg
  enddo

1000 write(message,'(a,i6)')'Maximum degree of GIA file: ',maxdeg
  call status_update('STATUS','UTIL','gia_over_oceans',GIA_stokes_file,message,0)

! dimension the arrays
  allocate(C_gia(0:maxdeg,0:maxdeg))
  allocate(S_gia(0:maxdeg,0:maxdeg))
  allocate(Plm(0:maxdeg,0:maxdeg))

! read them in
  rewind(luGIA)
  ioerr = 0
  do while (ioerr == 0)
    read(luGIA,*,iostat=ioerr,end=1001)tmpdeg,tmpord,tmpC,tmpS
!print*,tmpdeg,tmpord,tmpC,tmpS
    if(ioerr == 0)then
      C_gia(tmpdeg,tmpord) = tmpC
      S_gia(tmpdeg,tmpord) = tmpS
    endif
  enddo
1001 call status_update('STATUS','UTIL','gia_over_oceans',GIA_stokes_file,"Have read GIA spherical harmonic coefficients",0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get the mascon information
  call read_msc_hdr(lumsc,mascon_file,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                                ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
  call allocate_mascon_arrays
  call read_mascon_file(lumsc,mascon_file)  
  close(lumsc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! integrate the GIA signal - in terms of EWH - over the oceans. Because it is linear over decades, we need only do
! this calculation once, then use it as a rate correction

  ! loop over the mascons
  sum_volume = 0.d0
  sum_area = 0.d0
  nmsc_ocean = 0
  call status_update('STATUS','UTIL','gia_over_oceans',' ',"Calculating GIA rate integrated over ocean mascons",0)
  do imsc = 1,max_prim
!print*,'imsc,mcon_prim(imsc,1:2)',imsc,mcon_prim(imsc,1:2)*180.d0/pi
    if(mcon_prim(imsc,6) > 1010.d0)then    ! it is an ocean mascon
      nmsc_ocean = nmsc_ocean + 1
      sum_area = sum_area + mcon_prim(imsc,4)
      colat = (pi/2.d0-mcon_prim(imsc,1))  ! latitude  is in radians already
      rlong = mcon_prim(imsc,2)            ! longitude is in radians already

      ! generate the legendre polynomials
      call legendre_matrix(maxdeg,dcos(colat),Plm)
    
      tmpsum = 0.d0
      do ideg=2,maxdeg
        vfact = earthrad*p_av/(3.d0*p_w)*(2.d0*dble(ideg)+1.d0)/(1.d0+love_k(ideg))*1.d3
        tmpsum = tmpsum + vfact*Plm(ideg,0)*C_gia(ideg,0)

        do iord = 1,ideg
          tmpsum = tmpsum + vfact*Plm(ideg,iord)*(C_gia(ideg,iord)*dcos(iord*rlong)+S_gia(ideg,iord)*dsin(iord*rlong) )
        enddo
      enddo
      sum_volume = sum_volume + tmpsum * mcon_prim(imsc,4)

    endif

  enddo
  ! convert to an EWH (in mm) over the oceans
  GIA_rate = sum_volume/sum_area 
  write(message,'(a,f8.3,a)')'GIA rate integrated over the oceans = ',GIA_rate,' mm/yr'
  call status_update('STATUS','UTIL','gia_over_oceans',GIA_stokes_file,message,0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! loop through all the GRACE epochs and calculate the GIA correction integrated over all ocean mascons
  ioerr = 0
  nepochs = 1
  do while (ioerr == 0)
    read(lugrace,'(a)',iostat=ioerr,end=1002)line
    if(ioerr == 0)then
      ! get the epoch from the file name
      read(line(9:12),*)year
      read(line(14:15),*)month
      dec_yr(nepochs) = year + month/12.d0

      ! time difference from first requested epoch
      dt = dec_yr(nepochs) - dec_yr(1)
                
      ! output the information
      write(*,'(f20.6,f15.4,5x,i4,3x,i2.2)')dec_yr(nepochs),-1.d0*GIA_rate*dt,int(year),int(month)
      nepochs = nepochs + 1
    endif

  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

1002 call status_update('STATUS','UTIL','gia_over_oceans',GIA_stokes_file,"End of Program",0)

  end












  
