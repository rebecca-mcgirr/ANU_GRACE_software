    Program regularise_mascons

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
! PT180220: for a test, turn OFF the zero correlation between land and ocean. This is to test simulated data with signals that cross coasts
! PT180322: tightly constrain any primary mascon with < 30 ternary mascons in it .... why do these primarys even exist?
! PT180427: add a different sigma uncertainty for Alaska glacier mascons
! PT180831: read mascon values and use to adjust off-diagonal correlations when mascons are not similar
! PT180910: no longer use off-diagonal terms. Use input vcv file to apply larger uncertainties when (mascon EWH) > 1 m
! PT190507: add a separate sigma for each of Antarctic and Arctic ocean regions
! RM190521: added a seperate sigma for GIA regions
! PT190719: relax the regularisation in ocean areas related to the Argentine Gyre and the two large earthquakes
!           (Sumatra-Andaman and Tohoku-Oki)
! PT190904: changed to read information from a command file (rather than so much information on the command line)
! PT190920: loosen the constraints on W_Ant mascons where altimetry shows a large mass loss
! PT200110: add different constraints for Australia and Amazon, via the .cmd file
! PT200309: add different constraint for Patagonia, via the .cmd file
! PT200310: also added Baffin Island, via .cmd file
! PT200908: add a different constraint to the coastal mascons in Greenland in the V004c solution, hardwired by mascon number!

    use mascon_mod


    implicit none

    integer max_num_prime

    parameter (max_num_prime = 5000)
      
    real*8 a, a_prime, area, b, b_prime, c, c_prime, cos_a, d_prime, dist
    real*8 scale_land, scale_ocean, scale_Ant, scale_Grn
    real*8 pi, r, r_prime, rad_fact,  scale, scale_factor, variance, zero
    real*8 test_1, test_2
    real*8 junk

! variables to store the different sigmas to use in the regularisation
    real(kind=8) :: sigma_land, sigma_ocean          ! default values for land and ocean uncertainty
    real(kind=8) :: sigma_Ant_ext,sigma_Ant_int      ! sigmas for the outer ring of the Antarctic continent, and the interior of Antarctica
    real(kind=8) :: sigma_W_Ant                      ! for W Ant mass loss regions
    real(kind=8) :: sigma_Grn                        ! for Greenland mascons
    real(kind=8) :: sigma_alaska                     ! for Alaskan mascons
    real(kind=8) :: sigma_Caspian                    ! Caspian Sea mascons
    real(kind=8) :: sigma_Ant_O,sigma_Arctic_O       ! for the oceans around Antarctica and in the Arctic region
    real(kind=8) :: sigma_Sumatra                    ! region of the Sumatra-Andaman earthquake
    real(kind=8) :: sigma_Tohoku                     ! region of the Tohoku-Oki earthquake
    real(kind=8) :: sigma_gyre                       ! Atlantic gyre near the Faulkland Islands
    real(kind=8) :: sigma_Ant_GIA                    ! Antarctic GIA mascons
    real(kind=8) :: sigma_Ant_Alt                    ! Antarctic mascons where altimetry shows large mass loss
    real(kind=8) :: sigma_Aus                        ! Australian mascons
    real(kind=8) :: sigma_Amazon                     ! Amazon mascons
    real(kind=8) :: sigma_Patagonia                  ! Patagonia glacier mascons (just southern SAmerica, in fact)
    real(kind=8) :: sigma_Baffin                     ! Baffin Is mascons (just southern SAmerica, in fact)
    real(kind=8) :: sigma_Africa
    real(kind=8) :: sigma_NAmerica
    real(kind=8) :: sigma_SAmerica
    real(kind=8) :: sigma_Europe
    real(kind=8) :: sigma_Asia
    real(kind=8) :: sigma_Svalbard
    real(kind=8) :: sigma_Iceland
    real(kind=8) :: min_tern                         ! RM200312: variable to store threshold for tightly constraining small primarys

    real(kind=8), allocatable :: whole_matrix(:,:),msc_const(:,:)
    real(kind=8), allocatable :: lon_prime(:),lat_prime(:),rho(:),density(:)
    real(kind=8), allocatable :: sigma_msc(:)
    character*5,  allocatable :: msc_label(:)      ! character tag for the sigma region in which the mascon resides

    integer i, ioerr, j, k, l, m, n, count, num_cells
    integer londeg, lonmin, latdeg, latmin,  ht

! PT161128: unit numbers
    integer*4,parameter :: lumsc_in=21,luout=22,lusoln_in=23

! PT150722: variables for the length scale calculation
    integer*4    :: length_scl_pwr
    real(kind=8) :: average_colat,dist_ratio
    real(kind=8),allocatable :: colat_prime(:)
    logical :: zero_correl

! PT170706: command line option to set off-diagonal elements to zero
    logical :: zero_offdiag

! PT180831: variables for using mascon EWH values to scale off-diagonal correlations
    logical      :: msc_soln
    real(kind=8) :: d_EWH,max_EWH
    real(kind=8) :: ewh_scl
    real(kind=8) :: msc_EWH(50000)    ! we should never have more than 50000 mascons !
    integer*4    :: n_msc
    character*150:: line

    character mascon_file*150, correl_file*150, arg*50, message*256, soln_file*150

    integer*4    :: trimlen

! PT190904: new command file from which to take the regularisation information.
    character*100 :: command_file
    integer*4, parameter :: lucmd = 24

!PT200908: variables for constraining Greenland coastal mascons
    real(kind=8) :: coastal_Greenland_sigma
    logical      :: coastal_Greenland

! get the name of the mascon file from the command line
    call getarg(1,command_file)

    if ( command_file == " " ) then
      write(6,10)
10 format(/,'##############################################################################################' &
   ,/,'                                  Program regularise_mascons                    ' &
   ,/,'Purpose: Produce mascon correlation matrix using different length scales for land/ocean & lat/long' &
   ,//,'Usage: regularise_mascons command_file mascon-file output-file vcv_file (use "none" if you do not want one) ' &
   ,/,'##############################################################################################')
      stop
    else
      open(unit=lucmd,file=command_file,status='old',iostat=ioerr)
      if(ioerr /= 0)then
        write(message,'(a,a,a)')"Error opening input command file ",command_file,". Does it exist?"
        call status_update('FATAL','regularise_mascons ','regularise_mascons',' ',message,ioerr)
      endif
    endif

    call getarg(2, mascon_file)
    open(unit=lumsc_in,file=mascon_file,status='old',iostat=ioerr)
    if(ioerr /= 0)then
      write(message,'(a,a,a)')"Error opening input file ",mascon_file,". Does it exist?"
      call status_update('FATAL','regularise_mascons ','regularise_mascons',' ',message,ioerr)
    endif

    call getarg(3,correl_file)
    if ( correl_file == " " ) then
      write(message,'(a)')"Too few command line parameters."
      call status_update('FATAL','regularise_mascons ','regularise_mascons',' ',message,ioerr)
    else
      call status_update('STATUS','UTIL','regularise_mascons',correl_file,"Output msc regularisation file",0)
    endif

    call getarg(4,soln_file)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!! extract constraint information from the command file
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first, the sigmas for the typical mascon regions
  call command_storecommand(lucmd,'sigma_land'           ,sigma_land     ,1,1,"Land               sigma",1)
  call command_storecommand(lucmd,'sigma_ocean'          ,sigma_ocean    ,1,1,"Ocean              sigma",1)
  call command_storecommand(lucmd,'sigma_W_Antarctica'   ,sigma_W_Ant    ,1,1,"West Antarctica    sigma",1)
  call command_storecommand(lucmd,'sigma_Antarctic_ext'  ,sigma_Ant_ext  ,1,1,"Antarctic exterior sigma",1)
  call command_storecommand(lucmd,'sigma_Antarctic_int'  ,sigma_Ant_int  ,1,1,"Antarctic interior sigma",1)
  call command_storecommand(lucmd,'sigma_Greenland'      ,sigma_grn      ,1,1,"Greenland          sigma",1)
  call command_storecommand(lucmd,'sigma_Alaska'         ,sigma_alaska   ,1,1,"Alaska             sigma",1)

  call command_storecommand(lucmd,'sigma_Caspian'        ,sigma_Caspian  ,1,1,"Caspian Sea        sigma",1)
  call command_storecommand(lucmd,'sigma_Ant_ocean'      ,sigma_Ant_O    ,1,1,"Antarctic ocean    sigma",1)
  call command_storecommand(lucmd,'sigma_Arctic_ocean'   ,sigma_Arctic_O ,1,1,"Arctic ocean       sigma",1)

! Values for different regions/continents (self-explanatory). These three are optional (above are mandatory)
  sigma_Aus       = sigma_land
  sigma_Amazon    = sigma_land
  sigma_Patagonia = sigma_land
  sigma_Baffin    = sigma_land
  sigma_Africa    = sigma_land
  sigma_NAmerica  = sigma_land
  sigma_SAmerica  = sigma_land
  sigma_Europe    = sigma_land
  sigma_Asia      = sigma_land
  sigma_Svalbard  = sigma_land
  sigma_Iceland   = sigma_land
  call command_storecommand(lucmd,'sigma_Australia'      ,sigma_Aus      ,1,1,"Australia          sigma",0)
  call command_storecommand(lucmd,'sigma_Amazon'         ,sigma_Amazon   ,1,1,"Amazon             sigma",0)
  call command_storecommand(lucmd,'sigma_Patagonia'      ,sigma_Patagonia,1,1,"Patagonia          sigma",0)
  call command_storecommand(lucmd,'sigma_Baffin'         ,sigma_Baffin   ,1,1,"Baffin Is          sigma",0)
  call command_storecommand(lucmd,'sigma_Africa'         ,sigma_Africa   ,1,1,"Africa             sigma",0)
  call command_storecommand(lucmd,'sigma_NAmerica'       ,sigma_NAmerica ,1,1,"N America          sigma",0)
  call command_storecommand(lucmd,'sigma_SAmerica'       ,sigma_sAmerica ,1,1,"S America          sigma",0)
  call command_storecommand(lucmd,'sigma_Europe'         ,sigma_Europe   ,1,1,"Europe             sigma",0)
  call command_storecommand(lucmd,'sigma_Asia'           ,sigma_Asia     ,1,1,"Asia               sigma",0)
  call command_storecommand(lucmd,'sigma_Svalbard'       ,sigma_Svalbard ,1,1,"Svalbard           sigma",0)
  call command_storecommand(lucmd,'sigma_Iceland'        ,sigma_Iceland  ,1,1,"Iceland            sigma",0)
  if(sigma_Aus /= sigma_land)then
     write(message,'(a,a,e25.10,a)') "Australia          sigma",' found and stored (',sigma_Aus,')'
     call status_update('STATUS','GRACEFIT','command_storeCommand',' ',message,0)
  endif
  if(sigma_Amazon /= sigma_land)then
     write(message,'(a,a,e25.10,a)') "Amazon             sigma",' found and stored (',sigma_Amazon,')'
     call status_update('STATUS','GRACEFIT','command_storeCommand',' ',message,0)
  endif
  if(sigma_Patagonia /= sigma_land)then
     write(message,'(a,a,e25.10,a)') "Patagonia          sigma",' found and stored (',sigma_Patagonia,')'
     call status_update('STATUS','GRACEFIT','command_storeCommand',' ',message,0)
  endif
  if(sigma_Baffin /= sigma_land)then
     write(message,'(a,a,e25.10,a)') "Baffin Island      sigma",' found and stored (',sigma_Baffin,')'
     call status_update('STATUS','GRACEFIT','command_storeCommand',' ',message,0)
  endif
  if(sigma_Africa /= sigma_land)then
     write(message,'(a,a,e25.10,a)') "Africa             sigma",' found and stored (',sigma_Africa,')'
     call status_update('STATUS','GRACEFIT','command_storeCommand',' ',message,0)
  endif
  if(sigma_NAmerica /= sigma_land)then
     write(message,'(a,a,e25.10,a)') "North America      sigma",' found and stored (',sigma_NAmerica,')'
     call status_update('STATUS','GRACEFIT','command_storeCommand',' ',message,0)
  endif
  if(sigma_SAmerica /= sigma_land)then
     write(message,'(a,a,e25.10,a)') "South America      sigma",' found and stored (',sigma_SAmerica,')'
     call status_update('STATUS','GRACEFIT','command_storeCommand',' ',message,0)
  endif
  if(sigma_Europe /= sigma_land)then
     write(message,'(a,a,e25.10,a)') "Europe             sigma",' found and stored (',sigma_Europe,')'
     call status_update('STATUS','GRACEFIT','command_storeCommand',' ',message,0)
  endif
  if(sigma_Asia /= sigma_land)then
     write(message,'(a,a,e25.10,a)') "Asia               sigma",' found and stored (',sigma_Asia,')'
     call status_update('STATUS','GRACEFIT','command_storeCommand',' ',message,0)
  endif
  if(sigma_Svalbard /= sigma_land)then
     write(message,'(a,a,e25.10,a)') "Svalbard           sigma",' found and stored (',sigma_Svalbard,')'
     call status_update('STATUS','GRACEFIT','command_storeCommand',' ',message,0)
  endif
  if(sigma_Iceland /= sigma_land)then
     write(message,'(a,a,e25.10,a)') "Iceland            sigma",' found and stored (',sigma_Iceland,')'
     call status_update('STATUS','GRACEFIT','command_storeCommand',' ',message,0)
  endif

! now some exception regions due to specific geophysical signals
  call command_storecommand(lucmd,'sigma_Sumatra_Andaman',sigma_Sumatra  ,1,1,"Sumatra-Andaman    sigma",1)
  call command_storecommand(lucmd,'sigma_Tohoku_Oki'     ,sigma_Tohoku   ,1,1,"Tohoku-Oki         sigma",1)
  call command_storecommand(lucmd,'sigma_gyre'           ,sigma_gyre     ,1,1,"Atlantic gyre      sigma",1)

! GIA mascon sigmas
  call command_storecommand(lucmd,'sigma_Ant_GIA'        ,sigma_Ant_GIA  ,1,1,"Antarctic GIA      sigma",1)

! Antarctic Altimetry mascons
  call command_storecommand(lucmd,'sigma_Ant_Alt'        ,sigma_Ant_Alt  ,1,1,"Antarctic Alt      sigma",1)

! RM200312: variable to store threshold for tightly constraining small primarys (optional)
  min_tern =  100
  call command_storecommand(lucmd,'min_tern'      ,min_tern      ,1,1,"Minimum ternaries          ntern",0)
  if(min_tern /= 100)then
     write(message,'(a,a,i6,a)') "Minimum ternaries  ntern",' found and stored (',int(min_tern),')'
  else
     write(message,'(a,a,i6,a)') "Minimum ternaries  ntern",' set to default (',int(min_tern),')'
  endif
  call status_update('STATUS','GRACEFIT','command_storeCommand',' ',message,0)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT180831: read in a solution file for mascon EWH values. Can be either a .fit or a .vcv (gracefit/gracesim/addnorm)
  if(soln_file(1:4) /= "none" .and. soln_file(1:1) /= " ")then
    msc_soln = .true.
    open(lusoln_in,file=soln_file,status='old',iostat=ioerr)
    if(ioerr /= 0)then
      write(message,'(a,a,a)')"Error opening input soln file ",soln_file,". Does it exist?"
      call status_update('FATAL','regularise_mascons ','regularise_mascons',' ',message,ioerr)
    else
      call status_update('STATUS','regularise_mascons ','regularise_mascons',soln_file,"Use EWH info to regularise mascons",ioerr)
    endif
  else
    msc_soln = .false.
    call status_update('STATUS','regularise_mascons ','regularise_mascons',' ',"Will not use EWH info to regularise mascons",ioerr)
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output requested information to the screen
!    write(message,'(a,i9,i9,i9,i9,a)') "Scale parameters (lat/lon land, lat/lon ocean, lat/lon Antarctica, lat/lon Greenland: " &
!                        ,int(scale_land),int(scale_ocean), int(scale_Ant), int(scale_Grn)," metres"
!    call status_update('STATUS','regularise_mascons ','regularise_mascons',' ',message,0)

  write(message,'(a,e11.3,a,e11.3,a)')'Sigma between mascons on land: ',sigma_land,' Sigma over oceans: ',sigma_ocean,' m.'
  call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
  write(message,'(a,e11.3,a,e11.3,a)') "Sigma on Antarctic exterior/interior mascons : ",sigma_Ant_ext," /",sigma_Ant_int," metres"
  call status_update('STATUS','regularise_mascons ','regularise_mascons',' ',message,0)
  write(message,'(a,e11.3,a)') "Sigma on W_Ant     mascons       : ",sigma_W_Ant," metres"
  call status_update('STATUS','regularise_mascons ','regularise_mascons',' ',message,0)
  write(message,'(a,e11.3,a)') "Sigma on Greenland mascons       : ",sigma_Grn," metres"
  call status_update('STATUS','regularise_mascons ','regularise_mascons',' ',message,0)
  write(message,'(a,e11.3,a)') "Sigma on Alaska mascons          : ",sigma_alaska," metres"
  call status_update('STATUS','regularise_mascons ','regularise_mascons',' ',message,0)

  write(message,'(a,e11.3,a,e11.3,a)') "Sigma on Antarctic/Arctic ocean mascons : ",sigma_Ant_o," /",sigma_Arctic_o," metres"
  call status_update('STATUS','regularise_mascons ','regularise_mascons',' ',message,0)

  write(message,'(a,e11.3,a)') "Sigma on Sumatra-Andaman mascons : ",sigma_Sumatra," metres"
  call status_update('STATUS','regularise_mascons ','regularise_mascons',' ',message,0)
  write(message,'(a,e11.3,a)') "Sigma on Tohoku-Oki mascons      : ",sigma_tohoku," metres"
  call status_update('STATUS','regularise_mascons ','regularise_mascons',' ',message,0)
  write(message,'(a,e11.3,a)') "Sigma on Atlantic gyre mascons   : ",sigma_gyre," metres"
  call status_update('STATUS','regularise_mascons ','regularise_mascons',' ',message,0)

  write(message,'(a,e11.3,a)') "Sigma on Antarctic GIA mascons   : ",sigma_Ant_GIA," metres"
  call status_update('STATUS','regularise_mascons ','regularise_mascons',' ',message,0)

  write(message,'(a,e11.3,a)') "Sigma on Australian mascons      : ",sigma_Aus," metres"
  call status_update('STATUS','regularise_mascons ','regularise_mascons',' ',message,0)
  write(message,'(a,e11.3,a)') "Sigma on Amazon mascons          : ",sigma_Amazon," metres"
  call status_update('STATUS','regularise_mascons ','regularise_mascons',' ',message,0)

  pi = 4.0d0 * datan(1.0d0)
  rad_fact = pi/180.0d0
  r_prime = 6.3712d6
  zero = 0.0d0
  variance = 1.0d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT180831: if required, read in the mascon EWH values from the solution file
  if(msc_soln)then

! read the first line to find out what sort of file it is
    read(lusoln_in,'(a)')line
    if(line(1:15) == "V2 GRACEFIT  FIT" .or. line(1:16) == "V2 GRACESIM  FIT")then            ! it is a gracefit/gracesim .fit file. 

    else if (line(1:16) == "V2 GRACEFIT  VCV" .or. line(1:16) == "V2 GRACESIM  VCV" &
             .or. line(1:16) == "V2 ADDNORM   VCV")then   ! it is a gracefit/gracesim/addnorm .vcv file
      ! read through until we find a line with the mascon label
      do while (line(7:9) /= " MC")
        read(lusoln_in,'(a)')line
      enddo
      backspace(lusoln_in)
      ! now read until there are no more mascon lines
      do while(line(7:9) == " MC")
        read(lusoln_in,'(a)',iostat=ioerr,end=2001)line
        if(line(7:9) == " MC")then
          n_msc = n_msc + 1
          read(line(50:64),*)msc_ewh(n_msc)
        endif
      enddo
2001  write(message,'(a,i6,a)')"Have read",n_msc," mascon EWH values"
      call status_update('STATUS','regularise_mascons ','regularise_mascons',soln_file,message,0)

                            
    else if (line(1:16) == "V2 ADDNORM    VCV")then   ! it is an addnorm .vcv file
      call status_update('STATUS','ADDNORM','addnorm',soln_file,"reading addnorm vcv solution file",0)
      ! read through until we find a line with the mascon label
      do while (line(7:9) /= " MC")
        read(lusoln_in,'(a)')line
      enddo
      backspace(lusoln_in)
      ! now read until there are no more mascon lines
      do while(line(7:9) == " MC")
        read(lusoln_in,'(a)',iostat=ioerr,end=2002)line
        if(line(7:9) == " MC")then
          n_msc = n_msc + 1
          read(line(53:65),*)msc_ewh(n_msc)
        endif
      enddo
2002  write(message,'(a,i6,a)')"Have read",n_msc," mascon EWH values"
      call status_update('STATUS','regularise_mascons ','regularise_mascons',' ',message,0)

    else if (line(1:16) == "V2 ADDNORM   FIT")then   ! it is an addnorm .fit file
      ! read through until we find a line with the mascon label
      do while (line(7:9) /= " MC")
        read(lusoln_in,'(a)')line
      enddo
      backspace(lusoln_in)
      ! now read until there are no more mascon lines
      do while(line(7:9) == " MC")
        read(lusoln_in,'(a)',iostat=ioerr,end=2003)line
        if(line(7:9) == " MC")then
          n_msc = n_msc + 1
          read(line(63:78),*)msc_ewh(n_msc)
        endif
      enddo
2003  write(message,'(a,i6,a)')"Have read",n_msc," mascon EWH values"
      call status_update('STATUS','regularise_mascons ','regularise_mascons',' ',message,0)
    endif

  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT161128: use the new subroutines to read the mascon file
  call read_msc_hdr(lumsc_in,mascon_file,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                                ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
  call allocate_mascon_arrays
! PT180220: using trimlen here doesn't work - the subroutine is expecting a C*150, so pass the whole thing through
!  call read_mascon_file(lumsc_in,mascon_file(1:trimlen(mascon_file)))
  call read_mascon_file(lumsc_in,mascon_file)

!!! Dimension the entire matrix
  allocate(whole_matrix(total_prim,total_prim))
  allocate(msc_const(total_prim,total_prim))
  allocate(lat_prime(total_prim))
  allocate(colat_prime(total_prim))
  allocate(lon_prime(total_prim))
  allocate(density(total_prim))
  allocate(rho(total_prim))
  allocate(sigma_msc(total_prim))
  allocate(msc_label(total_prim))
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
! RM190521: added GIA case
    if(mcon_prim(i,6) > 3000.0)then
      sigma_msc(i) = sigma_Ant_GIA
      msc_label(i) = "GIA  "
    else if(mcon_prim(i,6) > 1000.0)then
      sigma_msc(i) = sigma_ocean
      msc_label(i) = "Ocean"
!print*,'ocean mascon ',i,' sigma ',sigma_ocean
    else 
      sigma_msc(i) = sigma_land
      msc_label(i) = "Land "
    endif

! PT200312: add a bunch of continent constraints
      ! SAmerica
      if(lon_prime(i)*180.0/pi > 258.d0 .and. lon_prime(i)*180.0/pi < 327.d0 &
        .and. lat_prime(i) > -58.d0 .and. lat_prime(i) < 14.d0  &
        .and. mcon_prim(i,6) < 1020.d0)then
        write(message,'(a,i7,2f9.4,a,f10.3)')'Setting  mascon',i,lon_prime(i)*180.0/pi,lat_prime(i) &
                                            ,' to SAmerica uncertainty of',sigma_SAmerica
        !call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
        sigma_msc(i) = sigma_SAmerica
        msc_label(i) = "SAmerica"
      endif
      ! NAmerica
      if(lon_prime(i)*180.0/pi > 190.d0 .and. lon_prime(i)*180.0/pi < 300.d0 &
        .and. lat_prime(i) > 13.d0 .and. lat_prime(i) < 70.d0  &
        .and. mcon_prim(i,6) < 1020.d0)then
        write(message,'(a,i7,2f9.4,a,f10.3)')'Setting  mascon',i,lon_prime(i)*180.0/pi,lat_prime(i) &
                                            ,' to NAmerica uncertainty of',sigma_NAmerica
        !call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
        sigma_msc(i) = sigma_NAmerica
        msc_label(i) = "NAmerica"
      endif
      ! Australia
      if(lon_prime(i)*180.0/pi > 112.d0 .and. lon_prime(i)*180.0/pi < 155.d0 &
        .and. lat_prime(i) > -45.d0 .and. lat_prime(i) < -10.d0  &
        .and. mcon_prim(i,6) < 1020.d0)then
        write(message,'(a,i7,2f9.4,a,f10.3)')'Setting  mascon',i,lon_prime(i)*180.0/pi,lat_prime(i) &
                                            ,' to Australia uncertainty of',sigma_Aus
        !call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
        sigma_msc(i) = sigma_Aus
        msc_label(i) = "Australia"
      endif
      ! Africa. Use a single rectangle, which will not encompass correctly the African continent
      if( (     (lon_prime(i)*180.0/pi > 0.d0 .and. lon_prime(i)*180.0/pi < 52.d0) &
           .or. lon_prime(i)*180.0/pi > 320.d0 ) &
        .and. lat_prime(i) > -38.d0 .and. lat_prime(i) < 35.d0  &
        .and. mcon_prim(i,6) < 1020.d0)then
        write(message,'(a,i7,2f9.4,a,f10.3)')'Setting  mascon',i,lon_prime(i)*180.0/pi,lat_prime(i) &
                                            ,' to Africa uncertainty of',sigma_Africa
        !call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
        sigma_msc(i) = sigma_Africa
        msc_label(i) = "Africa"
      endif

      ! Asia. Use a single rectangle, which will not encompass correctly the Asian continent
      if(lon_prime(i)*180.0/pi > 52.d0 .and. lon_prime(i)*180.0/pi < 180.d0 &
        .and. lat_prime(i) > -10.d0 .and. lat_prime(i) < 85.d0  &
        .and. mcon_prim(i,6) < 1020.d0)then
        write(message,'(a,i7,2f9.4,a,f10.3)')'Setting  mascon',i,lon_prime(i)*180.0/pi,lat_prime(i) &
                                            ,' to Asia uncertainty of',sigma_Asia
        !call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
        sigma_msc(i) = sigma_Asia
        msc_label(i) = "Asia"
      endif

      ! Europe. Use a single rectangle, which will not encompass correctly the Asian continent
      if( (  (lon_prime(i)*180.0/pi > 0.d0 .and. lon_prime(i)*180.0/pi < 52.d0) .or. lon_prime(i)*180.0/pi > 340.d0) &
        .and. lat_prime(i) > 35.d0 .and. lat_prime(i) < 85.d0  &
        .and. mcon_prim(i,6) < 1020.d0)then
        write(message,'(a,i7,2f9.4,a,f10.3)')'Setting  mascon',i,lon_prime(i)*180.0/pi,lat_prime(i) &
                                            ,' to Europe uncertainty of',sigma_Europe
        !call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
        sigma_msc(i) = sigma_Europe
        msc_label(i) = "Europe"
      endif

! PT170709: set the uncertainty for Greenland
! PT180808: apply the density test to ensure that the uncertainty is applied only to land mascons
    if(lat_prime(i) > 60.d0 .and. lat_prime(i) < 74.d0 .and. lon_prime(i)*180.0/pi > 300.0 .and. lon_prime(i)*180.0/pi < 338.5d0 &
            .and. mcon_prim(i,6) < 1020.d0 ) then   ! we are in Greenland
      sigma_msc(i) = sigma_Grn
      msc_label(i) = "Grnld"
    else if(lat_prime(i) > 74.d0 .and. lon_prime(i)*180.0/pi > 285.0 .and. lon_prime(i)*180.0/pi < 355.d0 &
            .and. mcon_prim(i,6) < 1020.d0 ) then   ! we are in northern Greenland
      sigma_msc(i) = sigma_Grn
      msc_label(i) = "Grnld"
    endif

!! PT141020: add a different constraint for Antarctica
! PT170517: only do this for land Antarctic mascons
! PT190507: now apply a different sigma for Antarctic oceans
      if(lat_prime(i) < -60.0 .and. mcon_prim(i,6) < 1020.d0)then
        ! PT190904: we now differentiate between mascons in the Antarctic interior and those around the periphery (+ west Antarctica)
        if(mcon_prim(i,2)*180.d0/pi > 200.d0 .and. mcon_prim(i,2)*180.d0/pi < 290.d0 .and. lat_prime(i) > -78.0)then
          sigma_msc(i) = sigma_W_Ant
        else if(mcon_prim(i,2)*180.d0/pi > 200.d0 .and. mcon_prim(i,2)*180.d0/pi < 260.d0 .and. lat_prime(i) > -88.0)then
          sigma_msc(i) = sigma_Ant_ext
        else if (mcon_prim(i,2)*180.d0/pi > 260.d0 .and. mcon_prim(i,2)*180.d0/pi < 340.d0 .and. lat_prime(i) > -84.0)then
          sigma_msc(i) = sigma_Ant_ext
        else if (lat_prime(i) > -73.0)then
          sigma_msc(i) = sigma_Ant_ext
        else
          sigma_msc(i) = sigma_Ant_int
        endif
        msc_label(i) = "Antar"
      else if(lat_prime(i) < -60.0 .and. mcon_prim(i,6) >= 1020.d0)then
        sigma_msc(i) = sigma_Ant_o
        msc_label(i) = "Ant_o"
      endif

!! PT190507: Arctic Ocean mascons now have a different sigma. Define as > 58 N, but > 60N around Alaska
!      if(      ( lon_prime(i)*180.d0/pi > 265.d0 .and. lon_prime(i)*180.d0/pi < 360.d0 .and. lat_prime(i) > 58.0 ) &
!          .or. ( lon_prime(i)*180.d0/pi < 265.d0 .and. lat_prime(i) > 60.0 )   &
!          .and. mcon_prim(i,6) >= 1020.d0)then
!        sigma_msc(i) = sigma_Arctic_o
!        msc_label(i) = "Arc_o"
!      endif

! PT200319: brute force testing to find error in the Arctic ...
      if( lat_prime(i) > 58.0 .and.      lon_prime(i)*180.d0/pi > -1.d0 .and. lon_prime(i)*180.d0/pi < 361.d0    &
          .and. mcon_prim(i,6) >= 1020.d0)then
!print*,"constrain Arctic mascon",i,' by ',sigma_Arctic_o, lat_prime(i),lon_prime(i)*180.d0/pi
        sigma_msc(i) = sigma_Arctic_o
        msc_label(i) = "Arc_o"
      else if (lat_prime(i) > 60.0 .and. mcon_prim(i,6) >= 1020.d0)then
        sigma_msc(i) = 1.d3
!if(mcon_prim(i,6) >= 1020.d0)print*,'not touching sigma for mascon',i,' which is ' &
!     ,sigma_msc(i),lat_prime(i),lon_prime(i)*180.d0/pi,mcon_prim(i,6)
      endif


! PT180427: use a different constraint for Alaskan mascons located in areas of significant mass loss
! RM200312: broadened area defined of Alaska
      if(lon_prime(i)*180.0/pi > 202.d0 .and. lon_prime(i)*180.0/pi < 224.d0 &
        .and. lat_prime(i) > 58.d0 .and. lat_prime(i) < 64.5d0 .and. mcon_prim(i,6) < 1020.d0)then
        write(message,'(a,i7,2f9.4,a,f10.3)')'Setting  mascon',i,360.d0-lon_prime(i)*180.0/pi,lat_prime(i) &
                                            ,' to Alaskan uncertainty of',sigma_alaska
        call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
        sigma_msc(i) = sigma_alaska
        msc_label(i) = "Alask"
      endif
      if(lon_prime(i)*180.0/pi > 220.d0 .and. lon_prime(i)*180.0/pi < 230.d0 &
        .and. lat_prime(i) > 57.d0 .and. lat_prime(i) < 62.d0 .and. mcon_prim(i,6) < 1020.d0)then
        write(message,'(a,i7,2f9.4,a,f10.3)')'Setting  mascon',i,360.d0-lon_prime(i)*180.0/pi,lat_prime(i) &
                                            ,' to Alaskan uncertainty of',sigma_alaska
        call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
        sigma_msc(i) = sigma_alaska
        msc_label(i) = "Alask"
     endif
      if(lon_prime(i)*180.0/pi > 224.d0 .and. lon_prime(i)*180.0/pi < 234.d0 &
        .and. lat_prime(i) > 51.d0 .and. lat_prime(i) < 58.0d0 .and. mcon_prim(i,6) < 1020.d0)then
        write(message,'(a,i7,2f9.4,a,f10.3)')'Setting  mascon',i,360.d0-lon_prime(i)*180.0/pi,lat_prime(i) &
                                            ,' to Alaskan uncertainty of',sigma_alaska
        call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
        sigma_msc(i) = sigma_alaska
        msc_label(i) = "Alask"
      endif

! RM200312: use different constraint for Svalbard mascons
      if(lon_prime(i)*180.0/pi > 5.d0 .and. lon_prime(i)*180.0/pi < 35.d0 &
        .and. lat_prime(i) > 76.d0 .and. lat_prime(i) < 81.d0 .and. mcon_prim(i,6) < 1020.d0)then
        write(message,'(a,i7,2f9.4,a,f10.3)')'Setting  mascon',i,360.d0-lon_prime(i)*180.0/pi,lat_prime(i) &
                                            ,' to Svalbard uncertainty of',sigma_Svalbard
        call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
        sigma_msc(i) = sigma_Svalbard
        msc_label(i) = "Svalb"
      endif

! RM200312: use different constraint for Iceland mascons
      if(lon_prime(i)*180.0/pi > 332.d0 .and. lon_prime(i)*180.0/pi < 349.d0 &
        .and. lat_prime(i) > 52.d0 .and. lat_prime(i) < 67.d0 .and. mcon_prim(i,6) < 1020.d0)then
        write(message,'(a,i7,2f9.4,a,f10.3)')'Setting  mascon',i,360.d0-lon_prime(i)*180.0/pi,lat_prime(i) &
                                            ,' to Svalbard uncertainty of',sigma_Iceland
        call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
        sigma_msc(i) = sigma_Iceland
        msc_label(i) = "Icela"
      endif

! PT190719: set a different ocean constraint around the Argentine Gyre
      if(lon_prime(i)*180.0/pi > 300.d0 .and. lon_prime(i)*180.0/pi < 330.d0 &
        .and. lat_prime(i) > -50.d0 .and. lat_prime(i) < -30.0d0 .and. mcon_prim(i,6) > 1020.d0)then
        write(message,'(a,i7,2f9.4,a,f10.3)')'Setting  mascon',i,lon_prime(i)*180.0/pi,lat_prime(i) &
                                            ,' to Argentine Gyre uncertainty of',sigma_gyre
        !call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
        sigma_msc(i) = sigma_gyre
        msc_label(i) = "Argentine Gyre"
      endif

! PT190719: set a 10 cm constraint for the ocean around the Sumatra-Andaman earthquake region
      if(lon_prime(i)*180.0/pi > 85.d0 .and. lon_prime(i)*180.0/pi < 110.d0 &
        .and. lat_prime(i) > -5.d0 .and. lat_prime(i) < 20.0d0 .and. mcon_prim(i,6) > 1020.d0)then
        write(message,'(a,i7,2f9.4,a,f10.3)')'Setting  mascon',i,lon_prime(i)*180.0/pi,lat_prime(i) &
                                            ,' to Sumatra-Andaman earthquake uncertainty of',sigma_Sumatra
        !call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
        sigma_msc(i) = sigma_Sumatra
        msc_label(i) = "Sumatra_Andaman_Ocean"
      endif

! PT190719: set a 10 cm constraint for the ocean around the Tohoku-Oki earthquake region
      if(lon_prime(i)*180.0/pi > 125.d0 .and. lon_prime(i)*180.0/pi < 160.d0 &
        .and. lat_prime(i) > 20d0 .and. lat_prime(i) < 50.0d0 .and. mcon_prim(i,6) > 1020.d0)then
        write(message,'(a,i7,2f9.4,a,f10.3)')'Setting  mascon',i,lon_prime(i)*180.0/pi,lat_prime(i) &
                                            ,' to Tohoku-Oki earthquake uncertainty of',sigma_tohoku
        !call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
        sigma_msc(i) = sigma_tohoku
        msc_label(i) = "Tohoku-Oki_Ocean"
      endif

! PT190904: Caspian Sea
!      if(prim_flags(i)(1:5) == "PCasp")then
      if(mcon_region(i)(1:4) == "Casp")then
        write(message,'(a,i7,2f9.4,a,f10.3)')'Setting  mascon',i,lon_prime(i)*180.0/pi,lat_prime(i) &
                                            ,' to Caspian Sea uncertainty of',sigma_Caspian
        call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
        sigma_msc(i) = sigma_Caspian
        msc_label(i) = "Caspian"
      endif

! PT190905: increase the sigma of mascons seaward of Pine Island glacier, where there is 8 m EWH increase over GRACE era
        if (mcon_prim(i,2)*180.d0/pi > 258.d0 .and. mcon_prim(i,2)*180.d0/pi < 275.d0 .and. lat_prime(i) > -75.0 &
             .and. lat_prime(i) < -70.0  .and. mcon_prim(i,6) < 1020.d0)then
        write(message,'(a,i7,2f9.4,a,f10.3)')'Setting  mascon',i,lon_prime(i)*180.0/pi,lat_prime(i) &
                                            ,' to 1xW_Ant uncertainty:',sigma_W_Ant*1.d0
        call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
          sigma_msc(i) = sigma_W_Ant * 1.d0
        endif

! PT200110: Australian mascons
      if(prim_flags(i)(1:4) == "PAUS")then
        write(message,'(a,i7,2f9.4,a,f10.3)')'Setting  mascon',i,lon_prime(i)*180.0/pi,lat_prime(i) &
                                            ,' to Australian uncertainty of',sigma_Aus
        call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
        sigma_msc(i) = sigma_Aus
        msc_label(i) = "AUS"//prim_flags(i)(5:6)
      endif

! PT200110: Amazon basin mascons
! PT200806: mascons_stage5_V004 no longer uses PAmaz code. Identify by mcon_region(i) now instead
!      if(prim_flags(i)(1:5) == "PAmaz")then
      if(mcon_region(i)(1:6) == "AMAZON")then
        write(message,'(a,i7,2f9.4,a,f10.3)')'Setting  mascon',i,lon_prime(i)*180.0/pi,lat_prime(i) &
                                            ,' to Amazon uncertainty of',sigma_Amazon
        call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
        sigma_msc(i) = sigma_Amazon
        msc_label(i) = "Amazon"
      endif

! PT200110: Patagonia mascons
      if(lon_prime(i)*180.0/pi > 260.d0 .and. lon_prime(i)*180.0/pi < 320.d0 &
        .and. lat_prime(i) > -57d0 .and. lat_prime(i) < -45.0d0  &
        .and. mcon_prim(i,6) < 1020.d0)then
        write(message,'(a,i7,2f9.4,a,f10.3)')'Setting  mascon',i,lon_prime(i)*180.0/pi,lat_prime(i) &
                                            ,' to Patagonia uncertainty of',sigma_Patagonia
        call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
        sigma_msc(i) = sigma_Patagonia
        msc_label(i) = "Patagonia"
      endif

! PT200110: Baffin Island mascons. Define by one rectangle
      if(lon_prime(i)*180.0/pi > 275.d0 .and. lon_prime(i)*180.0/pi < 298.8d0 &
        .and. lat_prime(i) > 64.7d0 .and. lat_prime(i) < 73.7d0  &
        .and. mcon_prim(i,6) < 1020.d0)then
        write(message,'(a,i7,2f9.4,a,f10.3)')'Setting  mascon',i,lon_prime(i)*180.0/pi,lat_prime(i) &
                                            ,' to Baffin Island uncertainty of',sigma_Baffin
        call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
        sigma_msc(i) = sigma_Baffin
        msc_label(i) = "Baffin"
      endif


! PT170703: special case of Hudson Bay. Change the density to treat it as "land" so that the GIA signal 
!           is not damped out as though it is ocean
    if(      mcon_prim(i,1)*180.d0/pi > 51.d0   .and. mcon_prim(i,1)*180.d0/pi < 70.5  &
             .and. mcon_prim(i,2)*180.d0/pi > 265.d0  .and. mcon_prim(i,2)*180.d0/pi < 296.d0)then
      sigma_msc(i) = sigma_land
      mcon_prim(i,6) = 1000.d0
      msc_label(i) = "Hudson"
    endif


! PT180322: try tighty constraining any primary mascon with < 10 ternary mascons in it
! PT200224: make it < 100
! RM200312: make it < min_tern (default = 100)
      if(mcon_prim(i,8) < int(min_tern))then
        sigma_msc(i) = 0.001
        write(message,'(a,i6,a,i4,a)')'mascon ',i,' has only ',int(mcon_prim(i,8)),' ternary mascons. Constrain to 1 mm'
        call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)

      endif




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  PT180910: scale the uncertainty if the input EWH < -0.5m   (ie significant mass loss 
      if(msc_soln)then

! PT190826: change this to relax the constraints on any mascon with |adjustment| > 0.5
!        if(msc_ewh(i) < -0.5d0 .and. msc_label(i) /= "Land " .and. msc_label(i) /= "Ocean")then
        if(dabs(msc_ewh(i)) > 0.2d0 )then
!        if(dabs(msc_ewh(i))  > 0.3d0 .and. msc_label(i) /= "Land " .and. msc_label(i) /= "Ocean")then
          write(message,'(a,i7,a,f8.3,a,f8.3,a,f8.3,3a,2f12.4)')'increase uncertainty on mascon',i,' from',sigma_msc(i),' to ' &
                                      ,sigma_msc(i)*3.d0," EWH value is ",msc_ewh(i)," m (",msc_label(i),")" &
                                      ,lon_prime(i)*180.d0/pi,lat_prime(i)
          call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
          sigma_msc(i) = sigma_msc(i)*3.d0
        else if(dabs(msc_ewh(i)) > 0.1d0 )then 
          write(message,'(a,i7,a,f8.3,a,f8.3,a,f8.3,3a,2f12.4)')'increase uncertainty on mascon',i,' from',sigma_msc(i),' to ' &
                                      ,sigma_msc(i)*2.d0," EWH value is ",msc_ewh(i)," m (",msc_label(i),")" &
                                      ,lon_prime(i)*180.d0/pi,lat_prime(i)
          call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
          sigma_msc(i) = sigma_msc(i)*2.d0
        endif
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    enddo


! PT200908: different constraint around the Greenland coastal mascons, by mascon number in the V004c mascon file
    coastal_Greenland = .true.
    if(coastal_Greenland)then
      coastal_Greenland_sigma = 0.5d0 
      sigma_msc(3233) = coastal_Greenland_sigma
      sigma_msc(3231) = coastal_Greenland_sigma
      sigma_msc(3244) = coastal_Greenland_sigma
      sigma_msc(3238) = coastal_Greenland_sigma
      sigma_msc(3258) = coastal_Greenland_sigma
      sigma_msc(3225) = coastal_Greenland_sigma
      sigma_msc(3215) = coastal_Greenland_sigma
      sigma_msc(3256) = coastal_Greenland_sigma
      sigma_msc(3210) = coastal_Greenland_sigma
      sigma_msc(3224) = coastal_Greenland_sigma
      sigma_msc(3252) = coastal_Greenland_sigma
      sigma_msc(3226) = coastal_Greenland_sigma
      sigma_msc(3250) = coastal_Greenland_sigma
      sigma_msc(3216) = coastal_Greenland_sigma
      sigma_msc(3223) = coastal_Greenland_sigma
      sigma_msc(3221) = coastal_Greenland_sigma
      sigma_msc(3245) = coastal_Greenland_sigma
    endif





! DEBUG
!PT190509: output the region for each mascon. Need to know this for EMSC3032_2019 assignment 4
!do k=1,total_prim
!  write(*,'(i7,a,a,f10.4,a)'),k," ",msc_label(k),sigma_msc(k)," region."
!enddo   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! M A I N    L O O P  !!!!!!!
    do k = 1, total_prim
      if(mod(k,500) == 0 )then
        write(message,'(a,i5)')'Correlations for mascon: ',k
        call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
      endif
      b = colat_prime(k)
      b_prime = dabs(dcos(b))
      do n = 1, total_prim
! PT190220: re-order this so that it is easy to turn off setting the coast/ocean combination to zero
        if ( k .eq. n ) then
          rho(n) = 1.d0/sigma_msc(n)**2
          msc_const(k,k) = rho(n)
        else if ( mcon_prim(k,6) .ne. mcon_prim(n,6) ) then
          rho(n) = zero
          msc_const(n,k) = 0.d0
        else
          c = colat_prime(n)
          c_prime = dabs(dcos(c))
          variance = sigma_msc(k)*sigma_msc(n)

! check whether we should be using the land or ocean scaling factors
          if ( mcon_prim(n,6) > 1000.d0 ) then
! PT200918: scale_ocean doesn't seem to have a value. Set it here to 1000 km
            scale_ocean = 600.d3
            scale_factor = scale_ocean/scale_ocean
            scale = scale_ocean
            zero_correl = .false.   ! use this to turn on/off correlation between ocean mascons
          else
            scale_land = 300.d3  ! PT191203: set this here since it is no longer read from the command line
            scale_Ant  = 300.d3
            scale_Grn  = 600.d3
            call geogr_scl_factors(mcon_prim(k,1:2),mcon_prim(n,1:2),scale_land,scale_Ant,scale_Grn,scale_Grn &
                                   ,mcon_prim(k,6), scale_factor,scale,zero_correl)

! PT180907: for a test, set off-diagonal correlations to zero for all land
            zero_correl = .true.

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
!          dist_ratio = dist/(scale/sin(average_colat))
          dist_ratio = dabs(dist/scale)
          rho(n) = 1.d0/variance * exp(-(dist_ratio**length_scl_pwr))


! PT170718: can we try setting the contribution to zero if the distance exceeds some value?
          if(dist > 1500.d3)rho(n) = 0.d0   ! anything greater than 1500 km won't contribute


! PT180831: scale the off-diagonal term by the relative difference in EWH values of each mascon, if requested to use a mascon solution
          if(msc_soln)then
            d_EWH = dabs(msc_EWH(k)-msc_EWH(n))
            max_EWH = max(dabs(msc_EWH(k)),dabs(msc_EWH(n)))

            if( max_EWH > 0.1d0)then              ! only calculate a variable scale value if at least one EWH value is > 0.1m
              ewh_scl = 1.0 * (1.d0 - 1.d0*d_EWH/max_EWH)
            else
              ewh_scl = 1.d0                      ! set it to 1. for small EWH values
            endif
          ! set the bounds to be zero or 1
            if(ewh_scl < 0.5d0)then
              ewh_scl = 0.d0
            else if(ewh_scl > 1.d0)then
              ewh_scl = 1.d0
            else
              ewh_scl = ewh_scl/1.d0
            endif

          else
            ewh_scl = 1.d0
          endif

! PT180904: try this again, putting just the rho(n) on the off-diagonal (the k == n case puts the value on the diagonal)
          if(.not. zero_correl) then
! PT200918: try permitting off-diagonal terms for ocean but not land
            if(mcon_prim(n,6) > 1010.d0) then   ! oean mascon
              msc_const(n,k) = rho(n)  !* ewh_scl
!print*,'n,k,mcon_prim(n,6),mcon_prim(k,6),rho(n)',n,k,mcon_prim(n,6),mcon_prim(k,6),rho(n)
            else
              msc_const(n,k) = 0.d0
            endif
          else
            msc_const(n,k) = 0.d0
          endif

          ! msc_const(n,k) = 0.d0
        endif

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

              if (lon_prime(n)*180.0/pi < 300.0) then  ! other point is in or west of Canada. This includes Siberia, Alaska etc
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
    call status_update('STATUS','UTIL','regularise_mascons',' ',"Writing out mascon constraint matrix",0)
    open(luout,file=correl_file,status='unknown')
! write header lines containing the mascon file code plus info on how this constraint file was created
    write(luout,'(a)')msc_hdr_lines(1)
    write(luout,'(a,4f10.2,a,11e11.3,a)')'# scale dists land, ocean, Ant, Grn:',scale_land/1.d3,scale_ocean/1.d3 &
                  ,scale_Ant/1.d3,scale_Grn/1.d3 &
                  ,' km. Sigmas land/ocean/W_Ant/Ant_ext/Ant_int/Grn/Alaska/Ant_alt/Casp/Patagonia/Baffin: ' &
                  ,sigma_land,sigma_ocean,sigma_W_Ant,sigma_Ant_ext,sigma_Ant_int,sigma_Grn,sigma_alaska,sigma_Ant_alt,sigma_caspian &
                  ,sigma_Patagonia,sigma_Baffin,' m.'
    do i=1,total_prim
      if(mod(i,500) == 0 )then
        write(message,'(a,i5)')'Writing out correlations for mascon: ',i
        call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
      endif
!! PT170216: set off-diagonal terms to zero
! PT170706: now control this with a command line option
! PT181202: to be sure, hardwire it here so that I know I'm not getting off-diagonal terms!
! PT200918: modify this to test what happens when we have off-diagonal correlations on only the ocean
      zero_offdiag = .true.
      if(zero_offdiag)then
        do j=1,total_prim
          if(j /= i .and. mcon_prim(j,6) < 1010.d0)msc_const(i,j) = 0.d0
        enddo
      endif
      write(luout,*)(msc_const(i,j),j=1,total_prim)
! DEBUG
!print*,mcon_prim(i,1:2)*180.d0/pi,1.d0/dsqrt(msc_const(i,i)),i,' regularisation sigmas' 
    enddo
    close(luout)
    call status_update('STATUS','UTIL','regularise_mascons',correl_file,"Written out mascon constraint matrix",0)

    stop
    end
! ************************************************************************************************************************************************



! ************************************************************************************************************************************************
  subroutine geogr_scl_factors(latlon1,latlon2,scale_land,scale_Ant,scale_Grn,scale_alaska,density,scale_factor,scale,zero_correl)

! subroutine to determine which scale distance to use (land, Ant or Grn), based on the locations of the two sites

  implicit none

  real(kind=8),  intent(in ) :: latlon1(2),latlon2(2)           ! coords (in radians) of the two mascons in question
  real(kind=8),  intent(in ) :: scale_land,scale_Ant,scale_Grn  ! user-defined length scales for continents, Antarctica and Greenland
  real(kind=8),  intent(in)  :: scale_alaska                    ! currently passed in as scale_Grn
  real(kind=8),  intent(in)  :: density                         ! density of the kth mascon
  real(kind=8),  intent(out) :: scale_factor,scale              ! the appropriate length scale to use for this pair of mascons
  logical     ,  intent(out) :: zero_correl

  real(kind=8) :: pi
  logical :: debug

  pi = 4.d0*datan(1.d0)
  debug = .false.

if(debug)print*,'site coords:',latlon1*180/pi,latlon2*180/pi

! is the first site in Antarctica?
  if(latlon1(1)*180.d0/pi < -60.d0)then   ! yes, it is!
    ! use the Antarctica length scale
    scale = scale_Ant

    if(debug)print*,'first site in Antarctica'    

    ! is the second site also in Antarctica ?
    if(latlon2(1)*180.d0/pi < -60.d0)then   ! yes, it is!
      if(debug)print*,'second site in Antarctica as well'
      scale_factor = 1.d0
      zero_correl = .true. ! PT180905: set true for a test (should really be false!)

!! PT170718: try treating differently east and west Antarctica. Make the scale factor for West Antarctica half that of East Ant. Also,
!!           set the off-diagonal terms to zero if not in the same hemisphere
!
!      if( (latlon1(2)*180.d0/pi < 180.0 .and. latlon2(2)*180.d0/pi .ge. 180.d0)  &
!                .or. (latlon1(2)*180.d0/pi > 180.0 .and. latlon2(2)*180.d0/pi .le. 180.d0) ) then  ! in different hemispheres
!!         zero_correl = .true.
!      else if ( latlon1(2)*180.d0/pi > 180.0 ) then  ! West Antarctica
!!         scale = scale / 2.d0
!      endif

      return

    else  ! the second site is not in Antarctica. Set the zero  correlation flag so that the off-diagonaly term will be computed as zero
      scale_factor = 1.d0  
      zero_correl = .true.
      return
    endif

! is the first site in Greenland ?
  else if (latlon1(2)*180.d0/pi > 60.0 .and. latlon1(2)*180.0/pi > 300.0 .and. latlon1(2)*180.0/pi < 338.0) then  ! yes, we are in Greenland
    if(debug)print*,'first site is in greenland'  
  ! is the second site also in Greenland ?
      if (latlon2(1)*180.d0/pi > 60.0 .and. latlon2(2)*180.0/pi > 300.0 .and. latlon2(2)*180.0/pi < 338.0) then  ! yes, we are in Greenland
      scale = scale_Grn
      scale_factor = 1.d0
      if(debug)print*,'both sites in greenland',latlon1*180/pi,latlon2*180/pi
      zero_correl = .true. ! PT180905: set true for a test (should really be false!)
      return

    else    ! the second site is not in Greenland. Set the scale factor to zero so that the correlation will be computed as zero
      scale = scale_Grn
      scale_factor = 1.d0  ! we can't set this to zero .. it makes the maths wrong in the main program
      zero_correl = .true.
      if(debug)print*,'first site in greenland, second not',latlon1*180/pi,latlon2*180/pi
      return
    endif

! PT180427: use a different constraint for Alaskan mascons located in areas of significant mass loss
  else if( (latlon1(2)*180.0/pi > 208.d0 .and. latlon1(2)*180.0/pi < 224.d0 &
        .and. latlon1(2) > 58.d0 .and. latlon1(2) < 63.5d0 .and. density < 1020.d0) .or. &
        (latlon1(2)*180.0/pi > 220.d0 .and. latlon1(2)*180.0/pi < 230.d0 &
        .and. latlon1(2) > 57.d0 .and. latlon1(2) < 59.5d0 .and. density < 1020.d0) .or. &   
        (latlon1(2)*180.0/pi > 224.d0 .and. latlon1(2)*180.0/pi < 233.d0 &
        .and. latlon1(2) > 52.d0 .and. latlon1(2) < 58.0d0 .and. density < 1020.d0) )then

        scale = scale_alaska
        scale_factor = 1.d0
        zero_correl = .true.


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
      zero_correl = .true.  ! PT180906: set to .true. for a test (should be .false. I think ... !)
    endif

    return
  endif

    



  end subroutine geogr_scl_factors





