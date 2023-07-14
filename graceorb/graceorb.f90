    program graceorb

! E-K Potter 24th Jan 2011
!   This program will read in the tim, ICS and integration time steps for each satellite

! List of variables:
! tin		integration start time
! efic		6 position and velocity ICs in earth fixed coordinates
! tim 		step size (in seconds)
! step		number of integration steps
! sat		GRACE satellite (A or B)

! MODIFIED: APP 121004
! To allow user defined data file names read from standard input
!
! PT130414: file "in_file" prepared by sh_graceorb will now be a command line argument, opened and read by graceorb rather than being
!           piped in as a default source
! PT130820: similarly, file "GRACE.input" will now be a command line argument (rather than reading it from in_file_X). Make it the first one.
! PT130905: bit-map which of the satellite ICs to update a priori
! PT140820: add reading of twice-per-rev parameters from input vcv file
! PT160902: revamp the reading of the mascon information, plus some general cleanup of code and moving things into subroutines.
!
! PT170609: resstructure of declaring the record length of the GTORB files. This info now defined in mod_subs/gtorb_mod
! PT170624: output the mascon number and apriori mcon_EWH (for debug)
! PT190917: permit mascon a priori values to be read from a different file from the ICs

    use sat_mod
    use bsscl_mod
    use inmod_mod
    use GPSant_mod    
    use accel_mod
    use dealias_mod
    use soln_mod          ! PT190312: changed from write_soln_mod
    use rotation_mod
    use gm_mod
    use mascon_mod        ! the new definition of the mascon arrays
! PT140328: pass the roll/pitch/yaw adjustment angles through to generate_inert_efixed
    use rpy_mod
! PT160913: for the new netcdf tide grid files
    use netcdf              ! generic netcdf stuff
    use tide_netcdf_mod     ! definitions/functions for interfacing with the tide grid file
    use gtorb_mod           ! PT170609: this declares the record length information for the GTORB file
    use accred_mod          ! PT180620: this declares the arrays for the accelerometer and star camera data
 
    implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer*4 :: step
    integer*4 :: i, j, k, n, num_read, ioerr, ideg, iord, icount
    integer*4 :: mcon_tot
    integer*4 :: param_offset
    integer*4 :: luin
    
    integer, dimension(3) :: tarray 
! PT140820: increased this to 10 parameters by adding once- and twice-per-rev along-track accelerations
    real(kind=8), dimension(10) :: efic
    real(kind=8), dimension(10) :: efadj
    real(kind=8) :: tin, tim

    logical bitmap
    character*2 version

! PT130414: declare the name of the file on the command line 
    character in_file*50,arg*50,line*150
    character frmt*20

! PT131107: array to store the input vcv file parameter descriptors
    character*30, allocatable :: prm_input(:) 

! PT130902: record number in binary GTORB file. Passed back to gracorb so that the mascon info can be written to the end of the file
! PT170609: now declared in mod_subs/gtorb_mod
!    integer*4 :: nrec
    character message_msc*36504

! PT140410: add another three command line values, so that we have 6 IC adjustments for pos/vel
    real(kind=8) :: pos_adj(3),vel_adj(3)

! PT140605: declare an array to store temporarily XYZ coords related to mascons
    real(kind=8) :: xyz(3)
! PT140612: pointers to the first mascon and tidal mascon parameters in the a priori list read from the input VCV file
    integer*4 :: imascon, imsc_tide,num_secondary_mscs
    real(kind=8) :: sum_land_area, sum_ocean_area

! PT190917: variables to permit reading mascon a priori values from a non-vcv file
    logical   :: found_msc_apriori
    integer*4 :: lumsc_apr


! PT181204: declare a variable of the total number of parameters. Required to use read_soln_v3 to read a priori values
!           from a vcv file
!    integer*4 :: total_params
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT140227: declare the value of pi
    pi = 4.d0*datan(1.d0)
    rad_fact = pi/180.d0


! PT160205: declare the min distance of mascon-to-satellite. This gets reduced as calcs are made in mascon_calc
    mindist_msc2sat = 1.d9

! allocate  dynamic variables here 
! PT181204: hardwire the total number of parameters to be big enough to get me out of trouble here
!           just before going to AGU!
!    total_params = 20*2+20000  ! 20 IC params per satellite plus 20000 mascons 
! RM190312: increasing the number of params to get me out of trouble
!    total_params = 20*2+48000  ! 20 IC params per satellite plus 20000 mascons 
!    allocate(prm_input(total_params))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!                READ COMMAND LINE VALUES                        !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call orbversn

! get the name of the control file (i.e. GRACE.input) from the command line 
    call getarg(1,input_file)

! PT210615: print runstring if no entry from command line
    if(input_file(1:1) == " ")then
      call status_update('STATUS','GRACEORB','graceorb',' ',"Runstring: graceorb GRACE.input in_file_A/B/C/D",0)
      stop
    endif

! get and open the command line file
    call getarg(2,in_file)
    open(unit=43,file=in_file,status='old',iostat=ioerr)
    if(ioerr /= 0)then
      call status_update('FATAL','GRACEORB','graceorb',in_file,"Error opening command line file. Does it exist?",ioerr)
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     rollerr = 0.0
     pitcherr = 0.0
     yawerr = 0.0

! PT140328: store these in rot_ang_errors, defined in rpy_mod.f90. Convert from milliradians to radian
! PT190405: I don't think this works - rpy_perturb is reset to zero in gravcalc, then passed to accelerom
!           which then recalculates rot_ang_errors. So setting it here has no impact!
    rot_ang_errors(1) = rollerr   * 1.d-3
    rot_ang_errors(2) = pitcherr  * 1.d-3
    rot_ang_errors(3) = yawerr    * 1.d-3

! PT140410: interpret these values as XYZ IC adjustments instead, and also read another three for the velocity adjustments
    pos_adj = 0.d0
    vel_adj = 0.d0
  
! APP130405: inserting angular correction terms for each axis. We are passing roll/pitch/yaw errors in here in radians.
!    call rpy2rotmat(rot_ang_errors(1),rot_ang_errors(2),rot_ang_errors(3) , rotcorr, rotcorr_deriv_roll, rotcorr_deriv_pitch &
!                    , rotcorr_deriv_yaw)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                             !!
! now read the information from the GRACEORB input file (replacing the read from unit=5)
!!                                                                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call read_graceorb_infile(43)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                             !!
!   read in ICs and integration step size/duration information
!!                                                                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    luin = 20
    call read_graceorb_IC(luin,efic,step, tin, GPSant_magL1, GPSant_magL2, GPSant_cosL1, GPSant_cosL2)
    call read_l1b_files()
    write(message,'(a,a,a,a,a,i5,a,f20.3)')"Integrating satellite GRACE ",sat," RL",RL,". Number of epochs: ",step &
                                          ,". Start time: ",tin
    call status_update('STATUS','GRACEORB','graceorb',' ',message,0)
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                             !!
! APP/PT130523:  read all the dealiasing information into one array
!!                                                                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   read in the dealiasing coefficients for each epoch in the dealiasing file
    luin = 50
    call read_dealias(luin,dealias_file)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                             !!
!   Read the GRACE.input file and print information to screen 
!!                                                                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   open GRACE.input file, which specifies which components of orbit integration to include, models to use, and
!   accelerometer bias and scale parameters
    lugt(1) = 22
    open (unit=lugt(1), file=input_file, status='old')
    call read_GRACEinput(1) ! reads in components of orbit integration to include and which models to use

!   print components of integration to the screen
    call inred_print

! PT170624: temporarily open a file to dump out mascon number and mcon_EWH beneath satellite
!    open(574,file="mcon_groundtrack.ewh",status='unknown')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! PT180620: moved read of accelerometer data and star camera data to here from ORBITCALC
! PT171003: only do this if we want to use the accelerometer data
! RM190703: call gmset to define GM for the Earth
    call gmset
    if (gt_acc(1).eq."Y")then
      call status_update('STATUS','GRACEORB','graceorb',' ' ,'     reading star camera and accelerometer data',0)
! PT180619: use new routine that handles either GRACE or GRACE FO data
!      call accred(tin)
      call read_acc_sca_thr(tin,efic(1:3))
    else
      call status_update('STATUS','GRACEORB','graceorb',' ' ,' NOT reading star camera and accelerometer data',0)
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                             !!
!!                   S E T U P   T H E   O C E A N   T I D E   G R I D         !!
!!                                                                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (gt_oceantide(1) .eq. 'G')then
!     call setup_ocean_tides_v2(ocean_tide_ht_file,tin)
! PT160913: update to reading the new netcdf version of the ocean tide grid
      call setup_netcdf_tides(ocean_tide_ht_file,tin)
    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!                   R E A D    T H E    M A S C O N S
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! PT170207: change to read the mascon header file irrespective of whether we want to model mascon partials. We need the dimensions
!           so that we can configure the ocean mascon arrays
!  if (gt_mcon(1).eq."Y") then
! read the mascon header
    call read_msc_hdr(21,combined_mascon_file,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                                ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
 ! endif
! open the ocean mascon file and read the header (all done in the subroutine)
  if (gt_oceantide(1).eq.'G') then
    call read_msc_hdr(22,ocean_mascon_file,msc_ocn_hdr_lines,max_ocn_hdr_lines,n_ocn_hdr_lines,msc_ocn_hdr_code &
                                ,max_ocean_prim,max_ocean_sec,max_ocean_tern &
                                ,max_ocean_tern_per_prim,max_ocean_tern_per_sec,max_ocean_sec_per_prim)
  endif

! if either, need to allocate the arrays
! PT210805: gt_mcon can have the following values:
!    N: we don't want anything to do with mascons
!    Y: we want to model the apriori mascon accelerations and compute partials
!    O: we want just mascon accelerations, no mascon partials
  if (gt_oceantide(1).eq.'G' .or. gt_mcon(1).eq."Y" .or. gt_mcon(1).eq."O" ) then
    ! allocate the array sizes for the mascons
    call allocate_mascon_arrays
  endif

! establish the information to identify each ternary mascon from lat/lon 
! PT/HMcQ 161201: updated to the ellipsoidal geometry version of the ternary macson pattern
  call tern_lat_bands_ell(ternary_lat_spacing/60.d0)


! PT160905: write a new routine to read the mascon file that contains all info for the primary, secondary and ternary mascons
  num_mcon = 0

  if (gt_mcon(1).eq."Y" .or. gt_mcon(1).eq."O" ) then
    if(gt_mcon(1).eq."Y" )write(message,'(a,a)')"   Will compute mascon accelerations and partials       : ", combined_mascon_file
    if(gt_mcon(1).eq."O" )write(message,'(a,a)')"   Will compute just mascon accelerations (no partials) : ", combined_mascon_file
    call status_update('STATUS','GRACEORB','graceorb',' ',message,0)

    ! read in the mascon information
    call read_mascon_file(21,combined_mascon_file)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Now we read in the deformational component of gravity. this is calculated
    ! for a point load as a function of angular distance and altitude
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! PT130213: use the actual deformation file name rather than hardwiring it
    call read_def_grav(26,def_grav_file)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Now we read in the primary mascon distance flags
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    luin = 26
    ! PT220912: turn off reading this file. We no longer use the information contained in it and it takes time to read
    !call read_prim_dist_flags(luin,mascon_flag_file)

    ! PT181023: add the geoid/ell separation and/or topography to mascon radii if necessary.
    call adjust_mascon_surface("GRACEORB",mascon_surface)

!! DEBUG
!do i=1,total_prim
!  print*,mcon_prim(i,1:2)*180.d0/pi,mcon_prim(i,3),i," primary lat, lon, rad"
!enddo
!stop
  else
    call status_update('STATUS','GRACEORB','graceorb',' ',"NOT Including MASCONS",0)
  endif

!!! and/or read ocean mascons
  if (gt_oceantide(1).eq.'G') then
    write(message,'(a,a)')"    Including OCEAN MASCON FILE       : ", ocean_mascon_file
    call status_update('STATUS','GRACEORB','graceorb',' ',message,0)

    call read_ocean_mascons(22,ocean_mascon_file)
    ! separate the ternary mascons into those involving tidal amplitudes and those not
    call ternary_tidal_ampl
  endif

! convert the mascon lat/lon/rad coords into cartesian xyz (quicker than doing it later at each epoch ...)
  call calc_mascon_xyz(gt_oceantide(1))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT190403: define the size of the GTORB record length required, based on the 
!           number of parameters. This is done here now (rather than hardwiring a
!           value) because we now know the number of mascon parameters, so can
!           dimension the record length to be only as large as necessary. This
!           helps to control the size of the output GTORB files.
!
! Line contains:  epoch (I4), pos/vel (6xR*8), quat(4xR*8), IC partials (36xR*8), bias (18xR*8), scale (18xR*8), 1/rev and 2/rev (4*R8), GPSant (3x6xR*8)
!                 , roll/pitch/yaw (3*6*R*8), 6x num_mascxR*8, 5x2 tidal amplitudes x 1000 tidal_masconsxR*8
! PT210805: define record length to include partials only if gt_mcon = "Y"
   if(gt_mcon(1).eq."Y") then
     GTORB_recl = 4 + (6+4+36+18+18)*8 +      (2+2+3+3)*6*8        + (max_prim)*6*8 + (1000*5*2)*6*8
     write(message,'(a,i6,a,i9)')"    GTORB record length (including ",max_prim," primary mascons) set to: ",GTORB_recl
     call status_update("STATUS","GRACEORB","graceorb",' ',message,0)
   else
     ! PT210812: remove the tidal mascon parameter partials if we are not wanting partials
     GTORB_recl = 4 + (6+4+36+18+18)*8 +      (2+2+3+3)*6*8        +          0*6*8 + (0*5*2)*6*8
     write(message,'(a,i9)')"    GTORB record length (no mascons) set to: ",GTORB_recl
     call status_update("STATUS","GRACEORB","graceorb",' ',message,0)
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!                   R E A D    T H E    A   P R I O R I    F I L E
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT130513: only open and read the apriori_file if information is required from it.
! PT190927: modify this to read only the ICs from the apriori_file. The mascons are now read from a different file.
  if (apriori_file.ne."".and.(gt_use_apriori_ics(1)> 0 ))then
    call INPUT_read_apriori(31,tin,apriori_file,imascon,imsc_tide,"N",0,found_msc_apriori)
  endif

! PT210823: if we want a priori mascons, get them from the msc_apriori_file
  if(gt_use_apriori_mcs(1) .eq. "Y" .and. .not. found_msc_apriori)then
    lumsc_apr = 24
    open(lumsc_apr,file=msc_apriori_file,status='old',iostat=ioerr)
    if(ioerr /= 0)then
      call status_update('FATAL','GRACEORB','graceorb',msc_apriori_file,'Error opening file. Does it exist?',ioerr)
    else
      call status_update('STATUS','GRACEORB','graceorb',msc_apriori_file,'Reading a priori mascon values',ioerr)
      ! HM190919: should only reach here if ICs were read from a list file by read_IC_list_v1
      !			  - therefore start storing mascon apriori values at index 19*2+1=39
      imascon = 39
      call input_read_msc_apriori(lumsc_apr)    
    endif
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! form a complete a priori and solution vector, being the combination of the ICs
! and the mascons. This is a bit convoluted, but should work .....
! PT210811: add the case where we don't want ICs but may want mascons. Here, the mascon a priori values go into array elements starting at row 39.
  if(gt_use_apriori_ics(1)> 0 .or. gt_use_apriori_mcs(1) .eq. "Y")then

! first, save a copy of the IC information
    allocate(apriori_IC(n_ICs))
    allocate(soln_IC(n_ICs))

    apriori_IC(1:n_ICs) = apriori(1:n_ICs)
    soln_IC(1:n_ICs)    = soln(1:n_ICs)

! now, deallocate the apriori and solution arrays
    if(gt_use_apriori_ics(1)> 0)then
      deallocate(soln)
      deallocate(apriori)
    endif

! now reallocate with the size of n_IC + total_prim
    allocate(apriori(n_ICs+total_prim))
    allocate(soln(n_ICs+total_prim))

! transfer the IC information into the new arrays
    apriori(1:n_ICs) = apriori_IC(1:n_ICs)
    soln(1:n_ICs) = soln_IC(1:n_ICs)

! now transfer the mascon information
    if(gt_use_apriori_mcs(1) == "Y")then
      imascon = n_ICs + 1
      apriori(n_ICs+1:n_ICs+total_prim) = mcon_prim_EWH(1:total_prim)
      soln(n_ICs+1:n_ICs+total_prim)    = mcon_prim_EWH(1:total_prim) 
    endif 
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!          U P D A T E    T H E    A   P R I O R I   P A R A M E T E R S
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(gt_use_apriori_ics(1) > 0)then
    call INPUT_update_apriori_params(apriori_file,gt_use_apriori_mcs(1),gt_use_apriori_ics(1) &
                                    ,imascon,total_prim,mcon_prim_EWH,efic)
  endif

  if((gt_mcon(1) == "Y" .or. gt_mcon(1) == "O") .and. gt_use_apriori_mcs(1) == "Y")then
    ! use the mascon apriori file for mascons instead of the one used for the ICs
    call INPUT_update_apriori_params(msc_apriori_file,gt_use_apriori_mcs(1),gt_use_apriori_ics(1) &
                                    ,imascon,total_prim,mcon_prim_EWH,efic)
  else if((gt_mcon(1) == "Y" .or. gt_mcon(1) == "O") .and. gt_use_apriori_mcs(1) == "N")then
    ! PT210823: set the mascon EWH arrays to zero
    mcon_prim_EWH = 0.d0
    mcon_sec_EWH  = 0.d0
    mcon_tern_EWH = 0.d0
    call status_update('STATUS','GRACEORB','graceorb',' ',"Setting all mascon apriori values to zero",0)
  endif

  if(gt_use_apriori_ics(1) == 0 .and. gt_mcon(1) == "N")then
    write(message,'(a,a)') "NOT updating any apriori parameter values for mascons or orbital ICs"
    call status_update('STATUS','GRACEORB','graceorb',' ',message,0)
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! A  P R I O R I    T I D A L    M A S C O N    P A R A M E T E R S
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT170302: initialize the values to zero before perhaps updating with a priori information
! PT170703: but only if we have at least one ocean primary mascon (ie tide modelling is on)
  if(max_ocean_prim > 0)msc_tide = 0.d0

  if( bitmap(gt_use_apriori_ics(1),8) )then   ! (value 128)
    gt_use_apriori_msc_tide(1) = "Y"
    call INPUT_update_msc_tide()
    write(message,'(a,a)') "    Using gracefit format apriori file for mascon tidal amplitudes"
  else
    write(message,'(a,a)') "NOT Using gracefit format apriori file for mascon tidal amplitudes" 
 gt_use_apriori_msc_tide(1) = "Y"  ! PT170303: set to "Y" otherwise the partials for tidal amplitudes are wrong in ICpart_calc ....
  endif
 call status_update('STATUS','GRACEORB','graceorb',' ',message,0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT140722: find out how many tidal constituent amplitudes there are to estimate
! PT170127: updated this code to use the ocean mascons, not the "gravity" mascons
    num_mcon_tide = 0
    do i=1,total_ocean_prim
      do j=1,max_msc_tides
        if( bitmap(mcon_ocean_prim(i,2),j) )then
          num_mcon_tide = num_mcon_tide + 1
        endif
      enddo
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    read (gt_intstpsz(1),*) tim
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Begin orbit calculation
    call etides_IERS2010    ! PT211013: defines the solid body tide coefficients, Doodson numbers etc
    call orbitcalc(tin, efic, tim, step) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                              !!
!!  O U T P U T    M A S C O N    A N D    T I D A L     I N F O R M A T I O N  !!
!!                                                                              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT130725: write out here the a priori mascon values to the end of the GTORB file
    if (.not. ish5) then
    nrec=nrec+1
    write(15,rec=nrec)"END OF PARTIALS"
    message = " "
    write(message,'(a18,i6)')"NUMBER OF MASCONS:",total_prim
    nrec=nrec+1
    write(15,rec=nrec)message
    write(message,'(a,i6,a,i8,a)')"Writing ",total_prim,' as number of mascon values to GTORB file (record ',nrec,')'
    call status_update('STATUS','GRACEORB','graceorb',' ',message,0)
    endif
    if(total_prim > 0)then
      nrec=nrec+1
! write out mascon a priori values
      if (ish5) then
        call h5_write_apr(orb_file, msc = mcon_prim_EWH )
      else
        write(15,rec=nrec)(mcon_prim_EWH(i),i=1,total_prim)
      endif
    endif

!!!! ************ &&&&&&&&&&&&&& ************** !!!!!!!!!!!!!
! write out the a priori amplitudes for sine/cos for every tide for every mascon


    nrec = nrec + 1
    if(.not. ish5) then 
      write(message,'(a,i7)')'NUMBER OF MASCON TIDAL AMPLITUDES:',num_mcon_tide*2
      write(15,rec=nrec)message
    endif 
    write(message,'(a,i8,a,i8,a)')"Written # mascon tidal amplitudes to GTORB file (record ",nrec,')'
    call status_update('STATUS','GRACEORB','graceorb',' ',message,0)
    if(num_mcon_tide > 0)then
      write(message,'(a,i8,a,i8,a)')"Writing ",num_mcon_tide*2 &
               ," a priori mascon tidal amplitudes to GTORB file (starting record ",nrec+1,')'
      call status_update('STATUS','GRACEORB','graceorb',' ',message,0)
! now write out the amplitude (sine then cosine) for each constituent for all mascons (most will be zero). One component per constituent per record.
      if (ish5) then
        call h5_write_apr(orb_file, msc_tide = msc_tide )
      else
        do j=1,max_msc_tides
          do k=1,2
            nrec = nrec + 1
            write(15,rec=nrec)(msc_tide(j,k,i),i =1,max_ocean_prim)
          enddo
        enddo
      endif
    endif
    if (ish5) then
      call orb_close(orb_file)
    else 
      close(15)
    endif

! PT160205: write out the minimum distance of each mascon to the satellite throughout the entire orbit integration
    open(15,file='msc_sat_mindist',status='unknown')
    do i=1,num_mcon
      write(15,'(2f10.5,e18.10,i5)')mcon_lng1(i)*180./pi,90.0 - mcon_colat1(i)*180./pi,mindist_msc2sat(i)/1.d3,mcon_rho1(i)
    enddo



    write(message,'(a,a)')"End of GRACEORB for satellite ",sat
    call status_update('STATUS','GRACEORB','graceorb',' ',message,0)
    
    end
