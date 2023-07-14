
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!    Input subroutines for GRACEORB
!!
!!  read_graceorb_infile  : reads the information from the in_file_A/B_YYYY-MM-DD file
!!  read_graceorb_IC      : reads the input IC file to get pos/vel etc
!!  inred_print           : output analysis options after reading them in 
!!  input_read_msc_apriori : read the apriori values for the mascons
!!  read_prim_dist_flag
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_graceorb_infile(luin)

! subroutine to read in the relevant file names for nearly all the input files required by graceorb
! Code moved from graceorb.f90 
!
! P. Tregoning
! 2 September 2016
!
! PT180606: read additional lines in the file, being information on how to model non-linear behaviour in L1B accelerometer data

  use inmod_mod

  implicit none

  integer*4     ,intent(in) :: luin         ! unit number of file (already opened by graceorb.f90)
  character*250             :: line
  integer   :: length,i,ioerr
  character*50:: dummy 
! read in all the file names (variables are declared and passed via inmod_mod.f90
    read (luin,31) line
    read (luin,31) output_file
    read (luin,31) IC_file
    read (luin,31) combined_mascon_file
    read (luin,31) ocean_mascon_file
    read (luin,31) mascon_flag_file
    read (luin,31) def_grav_file
! APP130130 : modified dealiasing calculation to only use one file
    read (luin,31) dealias_file
    read (luin,31) dummy
    read (luin,31) dummy
    read (luin,31) static_grav_coeff_file
    read (luin,31) pole_coeff_file
! PT130917: read the binary ocean tide height file name
    read(luin,31)  ocean_tide_ht_file         ! unit number 44 assigned in tide_mod.f90
! PT180606: read info on how to model non-linear accelerometer data for each axis
    do i=1,3
      read(luin,*)acc_fit_model(i),acc_fit_start(i),acc_fit_end(i)
31    format(a150)
    enddo
! PT181023: read whether to model mascons on the ellipsoid, on the geoid or on the topography
    mascon_surface = " "
    read(luin,31,iostat=ioerr,end=1000) mascon_surface
1000  if(ioerr /= 0 )then  ! old file that didn't contain this info
        mascon_surface = "ellipsoid"
      endif

! PT181023: until I turn this on, just hardwire the value. Three options: ellipsoid, geoid, topography
!  mascon_surface = "ellipsoid"

  length = len(trim(output_file))
  if( output_file(length-2:length) == '.h5') then
    ish5 = .true.
  else
    ish5 = .false.
  endif


  return
  end subroutine read_graceorb_infile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_graceorb_IC(luin,efic,step,tin, GPSant_magL1, GPSant_magL2, GPSant_cosL1, GPSant_cosL2)

! subroutine to read the IC and GPS antenna offset information from an input file
! Code moved from graceorb.f90
!
! P. Tregoning
! 2 September 2016

  use bsscl_mod        ! contains the declaration of the "scl(3)" array
  use inmod_mod
  use sat_mod          ! defines sat and RL

  implicit none

  integer*4   , intent(in) :: luin               ! input unit number for IC file
  real(kind=8),intent(out) :: efic(10)           ! earth-fixed position and velocity
  integer*4   ,intent(out) :: step               ! number of steps in integration
  real(kind=8),intent(out) :: tin                ! start time (in GRACE seconds) of orbit integration
  real(kind=8),intent(out) :: GPSant_magL1,GPSant_magL2,GPSant_cosL1(3),GPSant_cosL2(3)      ! GPS antenna vectors (mag and dirn cosines) 

  integer*4                :: ioerr


! open the IC file
  open (unit=luin, file=IC_file, status='old',iostat=ioerr)
  if(ioerr /= 0)then
    write(message,'(a,a,a)')"Error opening IC file ",IC_file,". Does it exist?"
    call status_update('FATAL','GRACEORB','graceorb',' ',message,ioerr)
  endif

! read the info from the IC file
  read (luin,*) tin      ! integration start time 
  write(message,'(a,f18.7)')"Integration start time: ",tin
  call status_update('STATUS','GRACEORB','graceorb',' ',message,0)
  read (luin,*) efic(1)  ! x IC 
  read (luin,*) efic(2)  ! y IC
  read (luin,*) efic(3)  ! z IC
  read (luin,*) efic(4)  ! xdot IC
  read (luin,*) efic(5)  ! ydot IC
  read (luin,*) efic(6)  ! zdot IC
  read (luin,*) step     ! number of steps in integration
  read (luin,*) sat      ! satellite name (A or B)
  read (luin,*) RL       ! RL
  read (luin,*) c0x(1), c0y(1), c0z(1) !PT18069: don't read any further, c1x(1), c1y(1), c1z(1), c2x(1), c2y(1), c2z(1)! acc bias parameters pre mjd=52705
  read (luin,*) c0x(2), c0y(2), c0z(2) !PT18069: don't read any further, c1x(2), c1y(2), c1z(2), c2x(2), c2y(2), c2z(2)! acc bias parameters post mjd=52705
  read (luin,*) scl(1), scl(2), scl(3) ! accelerometer scale parameters
  read (luin,*) GPSant_magL1, GPSant_cosL1(1), GPSant_cosL1(2), GPSant_cosL1(3) ! GPS antenna L1 vector (mag and dirn cosines)
  read (luin,*) GPSant_magL2, GPSant_cosL2(1), GPSant_cosL2(2), GPSant_cosL2(3) ! GPS antenna L2 vector (mag and dirn cosines)
! PT140820: by default, set the once- and twice-per-rev amplitudes for the along-track accelerations are zero - unless overwritten when reading 
!           in a vcv file
  efic(7:10) = 0.d0
! close the input file
  close(luin)

! PT180622: set the other bias coefficients to zero (we should remove them!)
  c1x = 0.d0
  c1y = 0.d0
  c1z = 0.d0
  c2x = 0.d0
  c2y = 0.d0
  c2z = 0.d0

  return
  end subroutine read_graceorb_IC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_l1b_files()

! subroutine written by Sebastien Allgeyer?

  use inmod_mod
  use sat_mod  
  implicit none
  character(len=150), dimension(2) :: GNV, THR, ACC, SCA
  character (len=150) :: KBR, LRI
  integer :: uin 
  integer*4 :: ioerr

  namelist /l1b_files/ GNV, THR, ACC, SCA, KBR, LRI
  uin = 0
  open(newunit=uin, file='L1B_files.txt',status='old',iostat=ioerr)
  if(ioerr /= 0)then
    call status_update('FATAL','GRACEORB','read_l1b_files','L1B_files.txt',"Error opening file. Does it exist?",0)
  endif

  read(uin, nml=l1b_files) 
  if ((sat .eq. "A") ) then
    SCA_file = SCA(1)
    ACC_file = ACC(1)
    THR_file = THR(1)
  else if (sat .eq. "B") then
    SCA_file = SCA(2)
    ACC_file = ACC(2)
    THR_file = THR(2)
  else if (sat .eq. "C") then 
    SCA_file = SCA(1)
    ACC_file = ACC(1)
    THR_file = THR(1)
  else if (sat .eq. "D") then
    SCA_file = SCA(2)
    ACC_file = ACC(2)
    THR_file = THR(2)
  end if

end subroutine read_l1b_files



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine inred_print

! subroutine to print out the analysis strategy options 
! Code moved from graceorb.f90
!
! P. Tregoning
! 2 September 2016

    use inmod_mod

    implicit none

    if (gt_planet(1).eq.'Y') then
      write(message,'(a)') "    Including PLANETARY POINT MASS EFFECTS"
    else
      write(message,'(a)') "NOT Including PLANETARY POINT MASS EFFECTS"
    endif
    call status_update('STATUS','GRACEORB','graceorb',' ',message,0)
  
    if (gt_sunmoon(1).eq.'Y') then
      write(message,'(a)') "    Including SUN AND MOON "
    else
      write(message,'(a)') "NOT Including SUN AND MOON "
    endif
    call status_update('STATUS','GRACEORB','graceorb',' ',message,0)
  
    if (gt_solidtide(1).eq.'Y') then
      write(message,'(a)') "    Including SOLID EARTH TIDE "
    else
      write(message,'(a)') "NOT Including SOLID EARTH TIDE "
    endif
    call status_update('STATUS','GRACEORB','graceorb',' ',message,0)
  
    if (gt_poletide(1).eq.'Y') then
      poltype = "both"
      write(message,'(a)') "    Including POLE TIDES (hard-wired both solid and ocean)"
    else
      write(message,'(a)') "NOT Including POLE TIDES"
    endif
    call status_update('STATUS','GRACEORB','graceorb',' ',message,0)
  
    if (gt_genrel(1).eq.'Y') then
      write(message,'(a)') "    Including GENERAL RELATIVITY"
    else
      write(message,'(a)') "NOT Including GENERAL RELATIVITY"
    endif
    call status_update('STATUS','GRACEORB','graceorb',' ',message,0)
  
    if (gt_acc(1).eq.'Y') then
      write(message,'(a)') "    Including ACCELEROMETERS"
    else
      write(message,'(a)') "NOT Including ACCELEROMETERS"
    endif
    call status_update('STATUS','GRACEORB','graceorb',' ',message,0)
  
    if (gt_statfield(1).eq.'Y') then
      write(message,'(a,a)') "    Including STATIC FIELD              : ", gt_statfieldmod(1)
    else
      write(message,'(a)') "NOT Including STATIC FIELD"
    endif
    call status_update('STATUS','GRACEORB','graceorb',' ',message,0)

! PT140204: allow this flag to be either 'Y' for ocean tides from spherical harmonics or 'G' for ocean tides from a grid  
    if (gt_oceantide(1).eq.'Y' .or. gt_oceantide(1) .eq. 'G') then
      write(message,'(a,a)') "    Including OCEAN TIDES               : ", gt_oceantidemod(1)
    else
      write(message,'(a)') "NOT Including OCEAN TIDES"
    endif
    call status_update('STATUS','GRACEORB','graceorb',' ',message,0)
  
    if (gt_atmtide(1).eq.'Y') then
      write(message,'(a,a)') "    Including ATMOSPHERIC TIDES         : ", gt_atmtidemod(1)
    else
      write(message,'(a)') "NOT Including ATMOSPHERIC TIDES"
    endif
    call status_update('STATUS','GRACEORB','graceorb',' ',message,0)
  
    if (gt_dealias(1).eq.'Y') then
      write(message,'(a,a)') "    Including NON_TIDAL OCEAN/ATMOSPHERE: ", gt_dealiasmod(1)
    else
      write(message,'(a)') "NOT Including NON_TIDAL OCEAN/ATMOSPHERE"
    endif
    call status_update('STATUS','GRACEORB','graceorb',' ',message,0)

    return
    end subroutine inred_print
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_dealias(luin,dealias_file)

! subroutine to read the dealias information. Code moved from graceorb.f90 to here
!
! P. Tregoning
! 5 September 2016
!
! MODS
! PT161018: modified to read the new format dealiasing file, which has a header line (containing max degree) then  deg  ord  C  S

  use dealias_mod

  implicit none

  integer*4, intent(in)    :: luin                   ! unit number of input dealias file (=50)
  character, intent(in)    :: dealias_file*150       ! name of dealiasing file

! local variables
  integer*4                :: ioerr
  integer*4                :: ideg, iord, icount,maxdeg_aod1b
  character*250            :: message
  integer*4                :: tmpdeg,tmpord
  real(kind=8)             :: tmpC,tmpS

  open (unit=luin, file=dealias_file, status='old', iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(message,'(a,a,a)')"Error opening dealiasing file ",dealias_file,". Does it exist?"
    call status_update('FATAL','GRACEORB','read_dealias',' ',message,ioerr)
  endif

  num_dealias_epochs = 0
  ioerr = 0

! PT161018: read the maximum degree of the dealiasing file
  read(luin,*)maxdeg_aod1b

  do while ( ioerr .eq. 0 )
    num_dealias_epochs = num_dealias_epochs + 1
    read(luin, *, end = 100, iostat = ioerr) epoch_time(num_dealias_epochs)
    icount = 0
    do ideg = 0 , maxdeg_aod1b
      do iord = 0 , ideg
        read(luin, *, iostat = ioerr) tmpdeg,tmpord,tmpC,tmpS
! PT161019: the new dealiasing file starts at degree zero but the graceorb code expects it
!           to start at degree 2 in the vector. Therefore, ignore the degree 0 and 1 terms
        if(tmpdeg > 1)then
          icount = icount + 1
          CT(num_dealias_epochs,icount) = tmpC
          ST(num_dealias_epochs,icount) = tmpS
        endif
      enddo
    enddo
  enddo
  if ( ioerr /= 0 ) then
    write(message,'(a,a,a,i1,a,i4,a,i4)')"Error reading dealiasing file ", dealias_file,". At epoch: ", num_dealias_epochs, &
                                         " Degree: ",ideg,", Order: ",iord
    call status_update('FATAL','GRACEORB','graceorb',' ',message,ioerr)
  endif
100 close(luin)
  num_dealias_epochs = num_dealias_epochs - 1
!  write(6,*) num_dealias_epochs


  end subroutine read_dealias
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_def_grav(luin,def_grav_file)

! subroutine to read in the gravitational effect of the elastic deformation due 
! to the load. This is probably insignificant, but we do it anyway.
! 
! Code moved to this subroutine from graceorb.f90
!
! P. Tregoning
! 6 September 2016

  use accel_mod       ! this pretty much defines all the variables that we need here

  implicit none

  integer*4, intent(in)    :: luin
  character*150,intent(in) :: def_grav_file

! local variables
  integer*4     :: i,j
  character*250 :: message

  open (unit=luin, file=def_grav_file, status='old')
  read(luin,*) theta0_def, dtheta_def, ntheta_def
  if ( ntheta_def .gt. max_def_colat ) then
    write(message,'(a,a,2i6)')"Too many deformation co-latitude values in ",def_grav_file,ntheta_def, max_def_colat
    call status_update('FATAL','GRACEORB','graceorb',' ',message,0)
  endif
  read(luin,*) rad0_def, drad_def, nrad_def
  if ( nrad_def .gt. max_def_rad ) then
    write(message,'(a,a,2i6)')"Too many deformation radius values in ",def_grav_file,nrad_def, max_def_rad
    call status_update('FATAL','GRACEORB','graceorb',' ',message,0)
  endif
  write(message,'(a,a)')"    Reading deformation values from file ",def_grav_file
  call status_update('STATUS','GRACEORB','graceorb',' ',message,0)
  do i = 1, ntheta_def
    do j = 1, nrad_def
      read(luin,*) def_tang(j,i), def_up(j,i)
    enddo
  enddo
  close(luin)


  end subroutine read_def_grav
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_prim_dist_flags(luin,mascon_flag_file)

! subroutine to read in the gravitational effect of the elastic deformation due 
! to the load. This is probably insignificant, but we do it anyway.
! 
! Code moved to this subroutine from graceorb.f90
!
! P. Tregoning
! 6 September 2016
!
! PT210817: of mascon_flag_file is "none" then just return

!  use accel_mod       ! this pretty much defines all the variables that we need here
  use mascon_mod

  implicit none

  integer*4, intent(in)    :: luin
  character*150,intent(in) :: mascon_flag_file

! local variables
  integer*4     :: ioerr,mcon_tot,i,j
  character*250 :: message
  character*100 :: msc_flag_hdr

! PT210817: don't do anything if "none" is the mascon flag file name
  if(trim(mascon_flag_file) == "none")then
    call status_update('STATUS','GRACEORB','read_prim_dist_flags',mascon_flag_file,"Not reading the mascon flag file",0)
    use_flag_file = .false.
    return
  else
    use_flag_file = .true.
  endif


! allocate the array for the primary mascon distance flags
!print*,"total_prim:",total_prim,"total_prim+1:",total_prim+1,"total_prim+2",total_prim+2,"num:",int((dble(total_prim+1) * dble(total_prim+2))/2.d0)
  allocate(mcon_flag(int((dble(total_prim+1) * dble(total_prim+2))/2.d0))) 
!  allocate(mcon_flag(((total_prim+1) * (total_prim+2))/2))
 
  open (unit=luin, file=mascon_flag_file, status='old',iostat=ioerr)
  if(ioerr .ne. 0)then
    write(message,'(a,a)')"Error opening mascon flag file: ",mascon_flag_file
    call status_update('FATAL','GRACEORB','read_prim_dist_flags',' ',message,0)
  endif

! PT161124: read the header line, then the second line has the number of primary mascons
  read(luin,'(a)')msc_flag_hdr
  write(message,'(a,a)')"    Mascon Flag header: ",msc_flag_hdr
  call status_update('STATUS','GRACEORB','read_prim_dist_flags',' ',message,0)

! check that it's compatible with the mascon input file
  if(msc_flag_hdr(1:8) /= msc_hdr_code)then
    write(message,'(a,a,a,a)')"Incompatible mascon and mascon flag file codes: ",msc_hdr_code," vs ",msc_flag_hdr(1:8)
    call status_update('FATAL','GRACEORB','read_prim_dist_flags',' ',message,0)
  endif

  ! skip the second line of the file
  read(luin,'(a)')message
  
  mcon_tot = 0
  write(message,'(a,a)')"    Now reading in mascon flags from file: ",mascon_flag_file 
  call status_update('STATUS','GRACEORB','read_prim_dist_flags',' ',message,0)
  do i = 1, total_prim
    do j = 1, i-1   ! PT/AP220104: bug here, it was "i" but needs to be "i-1"
      mcon_tot = mcon_tot + 1
      read(luin,*) mcon_flag(mcon_tot)
    enddo
  enddo 
  close (luin)

  return
  end subroutine read_prim_dist_flags
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine INPUT_read_apriori(luin,tin,apriori_file,imascon,imsc_tide,gt_use_apriori_mcs,total_prim,found_msc_apriori)

! open and read the input vcv file that contains a priori values for parameters
! Code moved to subroutine from graceorb.f90
!
! P. Tregoning
! 6 September 2016
!
! HM190502: modified read_IC_list_v1 to use soln_mod and allocate apriori arrays
! PT190917: set here a flag to indicate whether the a priori mascon values have been read or not

  use soln_mod      ! defines the prmnam array, soln and apriori vectors 

  implicit none

  integer*4, intent(in)       :: luin                 ! unit number of a priori VCV file
  character*80,intent(in)     :: apriori_file         ! name of a priori VCV file
  integer*4,intent(out)       :: imascon              ! pointer to the first mascon parameter in the apriori solution vector
  integer*4,intent(out)       :: imsc_tide            ! pointer to the first mascon tide parameter
  character*1,intent(in)      :: gt_use_apriori_mcs   ! need to know if we need to dimension apriori array to include mascons
  integer*4,intent(in)        :: total_prim           ! total number of primary mascons
  integer*4,intent(in)        :: tin                  ! GRACE seconds of requested ICs
  logical  , intent(out)      :: found_msc_apriori    ! true if a vcv file was read for a priori values, false if just ICs read

! local variables
  integer*4     :: ioerr,i
  character*250 :: message
  character*6   :: version
  integer*4     :: extra_params                       ! how many mascon params to add to dimensioning of "apriori" and "soln" (0 or total_prim)

! variables for reading the IC list file
!  real(kind=8)  :: ICs(19,2)    ! pos,ve,scl,bs,1pr,2pr,rpy = 19 parameters per satellite x 2 satellites

! PT181204: variables to pass to/from read_soln_v3
  integer*4     :: sTid,sScl,sBias,sEmp,sMasc,n_emp,nScl,nBias

  open(luin, file = apriori_file, status = "old",iostat=ioerr)

  if(ioerr /= 0)then
    write(message,'(a,a)')"Error opening gracefit format apriori file  : ", apriori_file
    call status_update('FATAL','GRACEORB','graceorb',' ',message,0)
  endif
  version = " "
  read(luin,'(a)')version(1:6)

  if ( version(1:2) .eq. "V2" ) then
    call read_soln_v2('GRACEORB',luin, 0,imascon,imsc_tide)
    if(gt_use_apriori_mcs == "Y")found_msc_apriori = .true.
    n_ICs = maxparm

  else if ( version(1:2) .eq. "V3" ) then
! PT181203: V3 has 5-digit mascon numbers
    call read_soln_v3('GRACEORB', luin, .false., sTid, sScl, sBias, sEmp, sMasc, n_emp, nScl, nBias)
    imascon = sMasc
    imsc_tide = sTid
    if(gt_use_apriori_mcs == "Y")found_msc_apriori = .true.
    n_ICs = maxparm

  else if(version(1:6) == "# YYYY")then
    n_ICs = 38
! search an IC list for the target date
!    if(gt_use_apriori_mcs == "Y")then
! PT190927: changed to be zero, since we only want to dimension for the ICs
!      !extra_params = total_prim
!      extra_params = 0
!    else
!      extra_params = 0
!    endif
    extra_params = 0
    call read_IC_list_v1(luin,tin,extra_params)
    imascon = 0
    imsc_tide = 0
    found_msc_apriori = .false.
  else
    write(message,'(a)')"Error: input apriori file not recognised. Cannot proceed"
    call status_update('FATAL','GRACEORB','graceorb',apriori_file,message,0)
  endif
  close(luin)


  return

  end subroutine INPUT_read_apriori
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine INPUT_read_msc_apriori(luin)

! subroutine to read only mascon a priori values from a file. The subroutine will read from any .fit or .vcv format (gracefit/gracesim/addnorm)
! This subroutine is required because we need a way of reading in a priori mascon values in cases where ICs are read from the file of only IC 
! values (i.e. not from a VCV file)
!
! P. Tregoning
! 17 September 2019

  use mascon_mod   ! provides the number of primary mascons 
  use soln_mod     ! provides the arrays for "a priori" and "soln" values for the mascons
  use inmod_mod    ! provides the values gt_use_apriori_mcs, gt_mcon

  implicit none

  integer*4    , intent(in)  :: luin              ! unit number of input a priori file

! local variables
  integer*4     :: imsc,ioerr,itern
  character*100 :: line,first_line
  character*3   :: apm_type

! first, check whether to impose zero apriori mascon values
  if( (gt_mcon(1) == "O" .or. gt_mcon(1) == "Y") .and. gt_use_apriori_mcs(1) == "N")then
    mcon_prim_EWH = 0.d0
    mcon_sec_EWH  = 0.d0
    mcon_tern_EWH = 0.d0
    call status_update('STATUS','GRACEORB','INPUT_read_msc_apriori',' ',"Setting apriori mascon values to zero",0)
    return
  endif

! read the first line of the input a priori file
  read(luin,'(a)')first_line

! PT210817: check whether this is a TERNARY EWH file or a fit/vcv file
  if(first_line(1:14) == "V3 TERNARY EWH")then
    call read_ternary_EWH(luin)
    return
  else

    ! read through to the first mascon file line
    do while (line(8:9) /= "MC"  .and. line(10:11) /= "MC")
      read(luin,'(a)')line
    enddo
    backspace(luin)
    ! find out how many mascons there are in the file
    ioerr = 0 
    do while (ioerr == 0)
      read(luin,'(a)',iostat=ioerr,end=1000)message
      if(ioerr == 0)total_msc = total_msc + 1
    enddo
1000  write(message,'(a,i8,a)')"A total of ",total_msc," mascons found in file"

    !PT210820: check that the number of mascons matches those in the mascon file
    if(total_msc /= total_prim)then
      write(message,'(a,i8,a,i8,a)')"Number of mascons in apriori file (",total_msc,") does not match number in mascon file (" &
                                         ,total_prim,")"
      call status_update('FATAL','GRACEORB','graceorb',' ',message,0)
    endif

    ! place the file at the line of the first mascon entry
    rewind(luin)
    line = " "
    do while (line(8:9) /= "MC" .and. line(10:11) /= "MC")
      read(luin,'(a)')line
    enddo
    backspace(luin)

    ! now, dimension the apriori mascon array
    allocate(apriori_msc(total_prim,2))

    ! 190919_HM :  stop if the apriori mascon file type can't be resolved 
    !
    ! what sort of file is it?
    if(first_line(15:17) == "FIT" .or. line(14:16) == "FIT")then
      apm_type="FIT"
    else if (first_line(15:17) == "VCV" .or. line(14:16) == "VCV" )then
      apm_type="VCV"
    else if (first_line(1:3) == "V2 " )then
      apm_type="VCV"
    else
      write(message,'(a)')"Error: input mascon apriori file not recognised. Cannot proceed"
      call status_update('FATAL','GRACEORB','graceorb',' ',message,0)
    endif

    line = " "
    if(apm_type == "FIT")then
      imsc = 0
      do while (imsc < total_prim)
        read(luin,'(a)')line
        if(line(8:9) == "MC")then
          imsc = imsc + 1
          read(line,'(35x,f12.5,12x,f19.5)')apriori_msc(imsc,1),apriori_msc(imsc,2)
        endif
      enddo

    else if (apm_type == "VCV")then
      do imsc=1,total_prim
        read(luin,'(a)')line
        read(line,'(30x,f17.7,f19.9)')apriori_msc(imsc,1),apriori_msc(imsc,2)
      enddo
    else
      print*,'apm_type is neither vcv nor fit - how did I get here?'
    endif

! PT210820: now need to fill in the mcon_prim_EWH and mcon_tern_EWH arrays
! RM/PT211011: fixed bug here (removed the sum of apriori_msc(:,2) and apriori_msc(:,1)
   mcon_prim_EWH = apriori_msc(:,2)
    mcon_sec_EWH = apriori_msc(:,2)
    do imsc = 1,total_prim
      do itern=1,nint(mcon_prim(imsc,8))
        mcon_tern_EWH(nint(mcon_tern(imsc,itern,7))) = mcon_prim_EWH(imsc)
      enddo
    enddo
 
  endif

  end subroutine INPUT_read_msc_apriori

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine INPUT_update_apriori_params(apriori_file,use_apriori_msc,use_apriori_ics,imascon, num_mcon,mcon_prim_EWH,efic)

! subroutine to update the apriori parameter values if so requested. Code moved
! to subroutine from graceorb.f90
!
! P. Tregoning
! 6 September 2016
!
! PT190312: updated to use mod_subs/soln_mod

  use soln_mod              ! defines the soln and apriori vectors
  use bsscl_mod             ! defines the bias and scale vectors
  use GPSant_mod            ! defines the GPS antenna vectors
  use sat_mod               ! defines the satellite A/B

  implicit none

  character*80,intent(in)     :: apriori_file        ! name of apriori VCV file
  character*1                 :: use_apriori_msc     ! Y/N
  integer*4                   :: use_apriori_ics     ! bit-mapped integer
  integer*4,    intent(in)    :: imascon             ! pointer to parameter number for first primary mascon
  integer*4, intent(in)       :: num_mcon            ! number of primary mascons
!  character*30,intent(in)     :: prm_input(maxparm)  ! character labels of input parameters
  real(kind=8), intent(out)   :: mcon_prim_EWH(num_mcon)  ! a priori EWH values for each primary mascon
  real(kind=8), intent(out)   :: efic(10)            ! earth-fixed a priori ICs

! local variables
  character*250  :: message
  integer*4      :: i
  integer*4      :: param_offset 

  logical bitmap


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first, check whether to update the mascon apriori values
  if (use_apriori_msc.eq."Y") then
    write(message,'(a,a)') "    Using mascon apriori values from gracefit apriori file : ", apriori_file
    call status_update('STATUS','GRACEORB','INPUT_update_apriori_params',' ',message,0)
    do i = 1, num_mcon
! PT140612  use pointer to first mascon as passed out of read_soln_v2
!        print*,'imascon,num_mcon,i',imascon,num_mcon,i,soln(imascon+i)
! PT160209: pointer "imascon" is to the first mascon, so need to subtract 1 here when pointing into the soln() vector
      mcon_prim_EWH(i) = soln(imascon + i -1)
    enddo
  else
    mcon_prim_EWH = 0.d0
    write(message,'(a,a)') "NOT Using mascon apriori values from gracefit apriori file : ", apriori_file
    call status_update('STATUS','GRACEORB','graceorb',' ',message,0)
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! DEBUG
!print*,'here is prmnam'
!print*,prmnam
!print*,'there it was'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! update ICs ?
    if (use_apriori_ics > 0) then
! extracting ICs for orbital parameters and accelerometer parameters depends on which satellite
! we are modelling
! copy across ICs for position  (value 1)
      if( bitmap(use_apriori_ics,1) )then
        call which_ICparam("X0  ",prmnam,sat,param_offset)
        if(param_offset > 0 )then
          call status_update('STATUS','GRACEORB','graceorb',' ',"    Updating orbital IC:  position",0)
          do i = 1, 3
! PT140410: also add on the IC adjustment from the command line (for the forward modelling test)
            efic(i) = soln(param_offset-1 + i) 
          enddo
        else
          call status_update('STATUS','GRACEORB','graceorb',' ',"NOT Updating orbital IC:  position (not found in file)",0)
        endif
      endif

! copy across ICs for velocity  (value 2)
      if( bitmap(use_apriori_ics,2) )then
        call which_ICparam("XV0 ",prmnam,sat,param_offset)
        if(param_offset > 0 )then
          call status_update('STATUS','GRACEORB','graceorb',' ',"    Updating orbital IC:  velocity",0)
          do i = 1, 3
! PT140410: also add on the IC adjustment from the command line (for the forward modelling test)
            efic(i+3) = soln(param_offset-1 + i) 
          enddo
        else
          call status_update('STATUS','GRACEORB','graceorb',' ',"NOT Updating orbital IC:  velocity (not found in file)",0)
        endif
      endif

! copy across accelerometer scales  (value 4)
      if( bitmap(use_apriori_ics,3) )then
        call which_ICparam("sclx",prmnam,sat,param_offset)
        if(param_offset > 0 )then
          call status_update('STATUS','GRACEORB','graceorb',' ',"    Updating orbital IC:  accelerometer scale value",0)
          do i = 1, 3
            scl(i) = soln(param_offset-1 + i)
          enddo
        else
          call status_update('STATUS','GRACEORB','graceorb',' '  &
                             ,"NOT Updating orbital IC:  accelerometer scale value (not found in file)",0)
        endif
      endif

! copy across quadratic coefficients for accelerometer biases  (value 8)
! We use a trivial quadratic formula where the coefficients of t and t^2 are set to 0
! Two identical sets of coefficients are used so the derived accelerometer biases will be independent of time
      if( bitmap(use_apriori_ics,4) )then
        call which_ICparam("bsx ",prmnam,sat,param_offset)
        if(param_offset > 0 )then
          call status_update('STATUS','GRACEORB','graceorb',' ',"    Updating orbital IC:  accelerometer bias value",0)
          do i = 1, 2
            c0x(i) = soln(param_offset)
            c0y(i) = soln(param_offset + 1)
            c0z(i) = soln(param_offset + 2)
            c1x(i) = 0.0d0
            c1y(i) = 0.0d0
            c1z(i) = 0.0d0
            c2x(i) = 0.0d0
            c2y(i) = 0.0d0
            c2z(i) = 0.0d0
          enddo
! PT140528: set a flag that the bias values are from an iterated file 
          vcv_bias = .true.
        else
          call status_update('STATUS','GRACEORB','graceorb',' ' &
                                          ,"NOT Updating orbital IC:  accelerometer bias value (not found in file)",0)
        endif
      else
        vcv_bias = .false.
      endif

! PT140820: copy across the once-per-rev sine/cos amplitudes  (value 16)
      if( bitmap(use_apriori_ics,5) )then
        call which_ICparam("1prS",prmnam,sat,param_offset)
        if(param_offset > 0 )then
          call status_update('STATUS','GRACEORB','graceorb',' ',"    Updating orbital IC:  once-per-rev along-track acc value",0)
          do i = 1, 2
            efic(6+i) = soln(param_offset-1 + i)
          enddo
        else
          call status_update('STATUS','GRACEORB','graceorb',' '  &
                             ,"NOT Updating orbital IC:  once-per-rev along-track acc value (not found in file)",0)
        endif
      endif

! PT140820: copy across the twice-per-rev sine/cos amplitudes (value 32)
      if( bitmap(use_apriori_ics,6) )then
        call which_ICparam("2prS",prmnam,sat,param_offset)
        if(param_offset > 0 )then
          call status_update('STATUS','GRACEORB','graceorb',' ',"    Updating orbital IC:  twice-per-rev along-track acc value",0)
          do i = 1, 2
            efic(6+i+2) = soln(param_offset-1 + i)
          enddo
        else
          call status_update('STATUS','GRACEORB','graceorb',' '  &
                             ,"NOT Updating orbital IC:  twice-per-rev along-track acc value (not found in file)",0)
        endif
      endif



! PT130528: copy across the GPS antenna offsets (value 64)
      if( bitmap(use_apriori_ics,7) )then
        call which_ICparam("GanX",prmnam,sat,param_offset)
        if(param_offset > 0 )then
          call status_update('STATUS','GRACEORB','graceorb',' ',"    Updating orbital IC:  GPS antenna vector",0)
          do i = 1, 3
            gpsantoff(i) = soln(param_offset-1 + i)
          enddo
        else
          call status_update('STATUS','GRACEORB','graceorb',' '  &
                                                ,"NOT Updating orbital IC:  GPS antenna vector (not found in file)",0)
        endif
      endif
    else
      write(message,'(a,a)') "NOT Using orbital ICs from gracefit apriori file : ", apriori_file
      call status_update('STATUS','GRACEORB','graceorb',' ',message,0)
    endif

  return
  end subroutine INPUT_update_apriori_params

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine INPUT_update_msc_tide()

! subroutine to work out which mascon tidal amplitude parameters require an
! update of the a priori amplitude values. Code moved to subroutine from graceorb.f90
!
! P. Tregoning
! 6 September 2016
!
! PT161205: updated to use the new mascon file arrays

  use accel_mod         ! defines tide arrays
  use soln_mod          ! defines the soln and apriori arrays.
  use mascon_mod

  implicit none

! local variables
  character*250 :: message
  real(kind=8)  :: pi
  integer*4     :: i,j
  logical       :: bitmap
  integer*4     :: param_offset

  pi = 4.d0*datan(1.d0)

! define here which tides we will permit to be adjusted. Variable "mcon_tide_name" defined in accel_mod.f90
  mcon_tide_name(1) = "M2"
  mcon_tide_name(2) = "O1"
  mcon_tide_name(3) = "S2"
  mcon_tide_name(4) = "K1"
  mcon_tide_name(5) = "K2"
  msc_tide = 0.d0                 ! set the a priori values to zero to begin with
! set the periods of these tidal constituents
  mcon_tide_period(1) = 2.d0*pi/(12.42d0*3600.d0)
  mcon_tide_period(2) = 2.d0*pi/(25.82d0*3600.d0)
  mcon_tide_period(3) = 2.d0*pi/(12.00d0*3600.d0)
  mcon_tide_period(4) = 2.d0*pi/(23.93d0*3600.d0)
  mcon_tide_period(5) = 2.d0*pi/(11.97d0*3600.d0)

  do i=1,total_ocean_prim
    if( mcon_ocean_prim(i,2) > 0) then    ! we want to estimate at least one set of tidal amplitudes for this mascon
      do j=1,max_msc_tides
        if( bitmap(mcon_ocean_prim(i,2),j))then
          call which_MSCtide_param(nparam,mcon_tide_name(j),prmnam,i,param_offset)
          if(param_offset > 0)then
            msc_tide(j,1,i) = soln(param_offset)
            msc_tide(j,2,i) = soln(param_offset+1)
            write(message,'(a,a,a,i5,a,2f8.4,a)')"    Updating ",mcon_tide_name(j)," tide for Mascon:",i &
                                            ,' (Ampl sin/cos: ',msc_tide(j,1,i),msc_tide(j,2,i),')'
            call status_update('STATUS','GRACEORB','graceorb',' ',message,0)
          endif
        endif
      enddo
    endif
  enddo


  return
  end subroutine INPUT_update_msc_tide
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




