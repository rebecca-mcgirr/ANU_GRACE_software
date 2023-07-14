    Program regularise_mascons_eclipse

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
!
! PT180914: stripped out all the off-diagonal code. Handle Greenland and Antarctica differently during times of eclipse. This is to
!           deal with the times at the end of the GRACE mission when there are no KBRR observations when the satellites are in eclipse (ie
!           over Greenland in NH winter, Antarctica in SH winter)

    use mascon_mod


    implicit none

    integer max_num_prime

    parameter (max_num_prime = 5000)
      
    real*8 a, a_prime, area, b, b_prime, c, c_prime, cos_a, d_prime, dist
    real*8 scale_land, scale_ocean, scale_Ant, scale_Grn
    real*8 pi, r, r_prime, rad_fact,  scale, scale_factor, sigma_ocean, sigma_land, sigma_Ant, sigma_grn, variance, zero
    real*8 sigma_alaska
    real*8 test_1, test_2
    real*8 junk

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

    character mascon_file*150, correl_file*80, arg*50, message*256, soln_file*150

    integer*4    :: trimlen

! PT180914: variables to deal with part of the world being in eclipse
    character*3  :: eclipse

! get the name of the mascon file from the command line
    call getarg(1,mascon_file)

    if ( mascon_file == " " ) then
      write(6,10)
10 format(/,'##############################################################################################' &
   ,/,'                                  Program regularise_mascons                    ' &
   ,/,'Purpose: Produce mascon correlation matrix using different length scales for land/ocean & lat/long' &
   ,//,'Usage: regularise_mascons_eclipse mascon-file output-file  ' &
   ,  ' sigma_land, sigma_ocean sigma_Ant sigma_Grn (sigmas in metres) sigma_Alaska fit_file/vcv_file Ant/Grn/Nil' &
   ,/,'##############################################################################################')
      stop
    endif

    open(unit=lumsc_in,file=mascon_file,status='old',iostat=ioerr)
    if(ioerr /= 0)then
      write(message,'(a,a,a)')"Error opening input file ",mascon_file,". Does it exist?"
      call status_update('FATAL','regularise_mascons ','regularise_mascons',' ',message,ioerr)
    endif

    call getarg(2,correl_file)
    if ( correl_file == " " ) then
      write(message,'(a)')"Too few command line parameters."
      call status_update('FATAL','regularise_mascons ','regularise_mascons',' ',message,ioerr)
    else
      call status_update('STATUS','UTIL','regularise_mascons',correl_file,"Output msc constraint file",0)
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read in distance scale factors for latitude & longitude on land and over ocean
    call getarg(3,arg)
    read(arg,*) sigma_land
    call getarg(4,arg)
    read(arg,*) sigma_ocean
    call getarg(5,arg)
    read(arg,*) sigma_Ant
    call getarg(6,arg)
    read(arg,*) sigma_Grn
    call getarg(7,arg)
    read(arg,*) sigma_alaska
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT180831: read in a solution file for mascon EWH values. Can be either a .fit or a .vcv (gracefit/gracesim/addnorm)
    call getarg(8,soln_file)
    if(soln_file(1:4) /= "none")then
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
! PT180914: read in whether it is Antarctica or Greenland without KBR observations
    call getarg(9,eclipse)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output requested information to the screen
    write(message,'(a,e11.3,a,e11.3,a)')'Sigma for mascons on land: ',sigma_land,' Sigma over oceans: ',sigma_ocean,' m.'
    call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
    write(message,'(a,e11.3,a,e11.3,a,e11.3,a)') "Sigma on Antarctic/Greenland/Alaska mascons : ",sigma_Ant," /" &
                                         ,sigma_Grn," /",sigma_alaska," metres"
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
      call status_update('STATUS','regularise_mascons ','regularise_mascons',' ',message,0)

!    else if (line(1:16) == "V2 ADDNORM   VCV")then   ! it is an addnorm .vcv file
!      call status_update('STATUS','ADDNORM','addnorm',soln_file,"reading addnorm vcv solution file",0)

    else if (line(1:16) == "V2 ADDNORM   FIT")then   ! it is an addnorm .fit file
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
          read(line(63:78),*)msc_ewh(n_msc)
        endif
      enddo
2002  write(message,'(a,i6,a)')"Have read",n_msc," mascon EWH values"
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

! DEBUG
open(932,file='ant_data.dat',status='unknown')

! PT170703: set up the uncertainties, including special cases
    do i=1,total_prim

! general cases
      if(mcon_prim(i,6) > 1000.0)then
          sigma_msc(i) = sigma_ocean
          msc_label(i) = "Ocean"
      else 
          sigma_msc(i) = sigma_land
          msc_label(i) = "Land "
      endif

! PT170703: special case of Hudson Bay. Change the density to treat it as "land" so that the GIA signal 
!           is not damped out as though it is ocean
      if(      mcon_prim(i,1)*180.d0/pi > 51.d0   .and. mcon_prim(i,1)*180.d0/pi < 70.5  &
               .and. mcon_prim(i,2)*180.d0/pi > 265.d0  .and. mcon_prim(i,2)*180.d0/pi < 296.d0)then
        sigma_msc(i) = sigma_land
        mcon_prim(i,6) = 1000.d0
        msc_label(i) = "Hudson"
      endif

! PT170709: set the uncertainty for Greenland
! PT180808: apply the density test to ensure that the uncertainty is applied only to land mascons
      if(lat_prime(i) > 60.0 .and. lon_prime(i)*180.0/pi > 300.0 .and. lon_prime(i)*180.0/pi < 338.0 &
            .and. mcon_prim(i,6) < 1020.d0 ) then   ! we are in Greenland
        sigma_msc(i) = sigma_Grn
        msc_label(i) = "Grnld"
      endif

!! PT141020: add a different constraint for Antarctica
! PT170517: only do this for land Antarctic mascons
      if(lat_prime(i) < -60.0 .and. mcon_prim(i,6) < 1020.d0)then
        sigma_msc(i) = sigma_Ant
        msc_label(i) = "Antar"
      endif

! PT180427: use a different constraint for Alaskan mascons located in areas of significant mass loss
      if(lon_prime(i)*180.0/pi > 206.d0 .and. lon_prime(i)*180.0/pi < 224.d0 &
        .and. lat_prime(i) > 58.d0 .and. lat_prime(i) < 63.5d0 .and. mcon_prim(i,6) < 1020.d0)then
        write(message,'(a,i7,2f9.4,a,f8.3)')'Setting  mascon',i,360.d0-lon_prime(i)*180.0/pi,lat_prime(i) &
                                            ,' to Alaskan uncertainty of',sigma_alaska
        call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
        sigma_msc(i) = sigma_alaska
        msc_label(i) = "Alask"
      endif
      if(lon_prime(i)*180.0/pi > 220.d0 .and. lon_prime(i)*180.0/pi < 230.d0 &
        .and. lat_prime(i) > 57.d0 .and. lat_prime(i) < 59.5d0 .and. mcon_prim(i,6) < 1020.d0)then
        write(message,'(a,i7,2f9.4,a,f8.3)')'Setting  mascon',i,360.d0-lon_prime(i)*180.0/pi,lat_prime(i) &
                                            ,' to Alaskan uncertainty of',sigma_alaska
        call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
        sigma_msc(i) = sigma_alaska
        msc_label(i) = "Alask"
     endif
      if(lon_prime(i)*180.0/pi > 224.d0 .and. lon_prime(i)*180.0/pi < 233.d0 &
        .and. lat_prime(i) > 52.d0 .and. lat_prime(i) < 58.0d0 .and. mcon_prim(i,6) < 1020.d0)then
        write(message,'(a,i7,2f9.4,a,f8.3)')'Setting  mascon',i,360.d0-lon_prime(i)*180.0/pi,lat_prime(i) &
                                            ,' to Alaskan uncertainty of',sigma_alaska
        call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
        sigma_msc(i) = sigma_alaska
        msc_label(i) = "Alask"
      endif

! PT180322: try tighty constraining any primary mascon with < 10 ternary mascons in it
      if(mcon_prim(i,8) < 10)then
        sigma_msc(i) = 0.001
        write(message,'(a,i6,a,i2,a)')'mascon ',i,' has only ',int(mcon_prim(i,8)),' ternary mascons. Constrain to 1 mm'
        call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)

      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  PT180910: scale the uncertainty if the input EWH < -0.5m   (ie significant mass loss 
      if(msc_soln)then

! when Antarctica is always in eclipse
        if(eclipse == "Ant" .and. msc_label(i) == "Antar" .and. mcon_prim(i,1) < -60.d0*pi/180.d0)then
          if(dabs(msc_ewh(i)) > 0.01d0 )then
            write(message,'(a,i7,a,f8.3,a,f8.3,a,f8.3,4a)')'increase uncertainty on Antarctic mascon',i,' from',sigma_msc(i) &
                                      ,' to ',dabs(msc_ewh(i))," EWH value is ",msc_ewh(i)," m (",msc_label(i),")"
            call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
            sigma_msc(i) = dabs(msc_ewh(i))
          else
            sigma_msc(i) = 0.001d0
          endif
! DEBUG
write(932,*)mcon_prim(i,1:2)*180.d0/pi,"MC",i,msc_ewh(i),sigma_msc(i)

! when Greenland is always in eclipse
        else if (eclipse == "Grn")then
          call status_update('STATUS','UTIL','regularise_mascons',' ',"Weight Greenland mascons in a special way during eclipse",0)
          call status_update('FATAL','UTIL','regularise_mascons',' ',"Code not yet written for Greenland mascons",0)

! when we don't want to handle eclipses in a special way.
        else
          if(msc_ewh(i) < -0.5d0 .and. msc_label(i) /= "Land " .and. msc_label(i) /= "Ocean")then
            write(message,'(a,i7,a,f8.3,a,f8.3,a,f8.3,4a)')'increase uncertainty on mascon',i,' from',sigma_msc(i),' to ' &
                                      ,sigma_msc(i)*2.d0," EWH value is ",msc_ewh(i)," m (",msc_label(i),")"
            call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
            sigma_msc(i) = sigma_msc(i)*2.d0
          endif
        endif
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    enddo

    close(932)
   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! M A I N    L O O P  !!!!!!!
    do k = 1, total_prim
      msc_const(k,k) = 1.d0/sigma_msc(k)**2
    enddo   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now write it out
    call status_update('STATUS','UTIL','regularise_mascons',' ',"Writing out mascon constraint matrix",0)
    open(luout,file=correl_file,status='unknown')
! write header lines containing the mascon file code plus info on how this constraint file was created
    write(luout,'(a)')msc_hdr_lines(1)
    write(luout,'(a,4f10.2,a,5e11.3,a)')'# scale dists land, ocean, Ant, Grn:',scale_land/1.d3,scale_ocean/1.d3 &
                  ,scale_Ant/1.d3,scale_Grn/1.d3 &
                  ,' km. Sigmas land/ocean/Ant/Grn/Alaska: ',sigma_land,sigma_ocean,sigma_Ant,sigma_Grn,sigma_alaska,' m.'
    do i=1,total_prim
      if(mod(i,500) == 0 )then
        write(message,'(a,i5)')'Writing out mascon: ',i
        call status_update('STATUS','UTIL','regularise_mascons',' ',message,0)
      endif

      write(luout,*)(msc_const(i,j),j=1,total_prim)
    enddo


    close(luout)
    call status_update('STATUS','UTIL','regularise_mascons',correl_file,"Written out mascon constraint matrix",0)

    stop
    end
! ************************************************************************************************************************************************



