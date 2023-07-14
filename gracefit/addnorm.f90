  program addnorm 
! Program to read in a series of normal equations and their RHS, stack them and then invert for the solution of all the parameters
!
! The normal equations will be stacked with the mascon parameters in the top left, then each set of ICs making blocks down the diagonal.
! There will be, of course, values relating each set of ICs to the mascons, but there will be no values relating each set of ICs to any other 
! set of ICs. 
!
! So the stacked normal equations will have:
! 1) the top nmascon_t rows full, 
! 2) the first nmascon_t columns full,  
! 3) blocks of ICs down the diagonals (with potentially different numbers of ICs per block - this is not important)
! 4) zeroes everywhere else
!
! P. Tregoning
! 28 October 2013
!
! MODS
! PT140605: add (temporary?) constraints to the scale factor estimates so that they don't adjust from their a priori values
! PT140709: add the handling of mascon tidal amplitudes (new parameters in graceorb/gracefit)
! PT161130: updated to read and ignore the header lines in the regularisation file
! PT200108: get addnorm output filename stem from the command line
! PT200112: read the appropriate mascon file so that the CoM constraint can be applied taking the different areas of the primary
!           mascons into account. May also make it that sum(land) = -sum(ocean)
!           Mascon file name to be read from the command line after the gracesim/gracefit argument.
! PT/JP200122 : read mascon time series information from a file containing annual ampl, RMS values, low-frequency signals etc
! RM210208: Add option to run SVD on AtWA 6th character flag passed through .cmd
! PT210216: modified to:
!             a) write out a binary file of stacked AtWA and AtWb, and
!             b) be able to read in said file so that making monthly solutions becomes quicker
!
! PT220527: add option to scale by the number of files found in a snorm
! PT220603: scale by ndays/2 for GRACE-FO solutions (2 .norm files per day with 12-hour inversions)
! PT220607: include the correlations between mascons into the off-diagonal regularisation terms
! PT220620: add adaptive regularisation option

  use mascon_mod   ! all the matrices for reading the mascon attributes
  use addnorm_mod  ! all the variables needed for reading the addnorm command file

  implicit none


  character     :: arg*100             ! variable to read the command line arguments
  character     :: outfile_stem*100    ! read from command line and used to create output file names
  character*300 :: message             ! string to use in status_update calls
  character*100 :: cmdfile             ! controlling command file name
  character*100 :: gracesim_cmdfile    ! KS160701 gracesim command file name
  character*150 :: mascon_file         ! mascon file name
  character*100 :: line                ! variable to read first line of addnorm command file
  character*1   :: char1,char2,char3   ! variables to read the analysis strategy from the input command file
  character*1   :: char4,char5,char6   ! variables to read the analysis strategy from the input command file

  integer*4,parameter :: lucmd=200    ! 200 is the unit number of command file
  integer*4,parameter :: lufile=201   ! 201 is the unit number of each normal equation file (binary)
  integer*4,parameter :: luWRMS=991   ! 991 is the unit number of the file containing a vector of mascon constraints

!  integer*4     :: nfiles             ! number of normal equation files to be added
  integer*4     :: output_nfiles      ! number of normal equation files, output to fit/vcv files
  integer*4     :: nparam_t,nmascon_t ! number of parameters and mascons for the actual inversion
  integer*4     :: nIC_t              ! number of IC parameters for inversion
  integer*4     :: nmsc_tidal_t       ! number of mascon tidal amplitude parameters for inversion

  integer*4     :: max_msc_tide = 5   ! maximum number of mascon tidal constituents

! pointers to rows and columns
  integer*4     :: iIC                ! first row/column for the first IC parameter in the combined normal equations
  integer*4     :: iArctic_msc        ! first row/column for the first of the time-varying mascons in the Arctic region
  integer*4     :: iMascons           ! first row/column for the first static mascon
  integer*4     :: iparm              ! parameter counter for outputting the solution
  integer*4     :: iMsc_tidal         ! first row/column for the first mascon tidal amplitude

! regularisation file mascon code
  character     :: mascon_code*6      ! mascon code as read from first line of header of regularisation file
  character     :: diag_only*9        ! RM201125 flag for diagonal only reg file

  logical IC_only
  logical mascon_only
  logical mascon_plus_ICs
  logical tidal_msc

!  character, allocatable     :: normeq_files(:)*100   ! names of normal equation files to be processed
  real(kind=8 ), allocatable :: tmpvals(:)            ! line-by-line read of the binary normal equations
  integer*4, allocatable     :: rec_len(:)            ! record length of each binary file to be read
  integer*4, allocatable     :: nparams(:)            ! number of parameters in each binary file
  integer*4, allocatable     :: nmascons(:)           ! number of mascons in each binary file
  integer*4, allocatable     :: nmsc_tidal(:)         ! number of mascons in each binary file
  integer*4, allocatable     :: nICs(:)               ! number of IC parameter in each binary file
  real(kind=8 ), allocatable :: versions(:)           ! version numbers of binary normal equation files
  real(kind=8 ), allocatable :: normeq_all(:,:)       ! compiled normal equations to be inverted
  real(kind=8 ), allocatable :: msc_const(:,:)        ! length-dependent constraints on mascons
  real(kind=8 ), allocatable :: VCV(:,:)              ! compiled normal equations to be inverted
  real(kind=8 ), allocatable :: AtWb(:)               ! stacked RHS to be used in the soluton
  real(kind=8 ), allocatable :: adjust(:)             ! LS solution of (AtWA)^-1 * AtWb
  real(kind=8 ), allocatable :: apr_tmp(:)            ! a priori parameter values (read from 2nd line of input binary files)
  real(kind=8 ), allocatable :: apr_prm(:)            ! a priori msc and/or IC parameter values 
  character*30,  allocatable :: prmnam(:)             ! character labels for each parameter to be estimated (used in output file)
  character*30,  allocatable :: prmnam_allICs(:)      ! character labels for each parameter to be estimated, including all ICs for all days (used in output file)
  character*30,  allocatable :: prmnam_all(:)         ! character labels for each parameter to be estimated, including all ICs for all days (used in output file)
  character*30,  allocatable :: prm_tidal(:)          ! character labels for each tidal amplitude parameter to be estimated (used in output file)
  integer*4                  :: msc_offset            ! amount by which mascon parameter codes need to be shifted forward to put mascons before ICs
  integer*4                  :: IC_offset             ! amount by which IC parameter codes need to be shifted to put ICs after mascons

! date parameters
  integer*4, allocatable     :: year(:),month(:),day(:),hr(:),minute(:),sec(:)      ! epoch of start of orbit integration to which the normal equations refer
  real(kind=8 ), allocatable :: duration(:)                                         ! duration of the orbit integration  
  integer*4,  allocatable     :: param_type(:)                                       ! (coded) parameter types on the binary norm file
  character*1, allocatable   :: sat1(:),sat2(:)                                     ! sat code for IC parameter descriptions on binary norm file
! ocean/land mascon flags
  logical,     allocatable   :: mcon_ocean(:)          ! ocean/land logical flag for each mascon (T for ocean, F for land)
  integer*4,   allocatable   :: mcon_tides(:,:)        ! bit-mapped value of mascon tidal amplitudes
  integer*4,   allocatable   :: tidal_lookup(:,:,:)    ! total compiled list of the mascons and tide amplitudes to be estimated
  integer*4,   allocatable   :: tide_list_epoch(:,:,:) ! list of mascon/tidal constituents estimated per set of normal equations
  integer*4                  :: k_tidal                ! counter for looking up the tidal lookup table
! constraint values for all parameters
  real(kind=8), allocatable  :: param_const(:)         ! diagonal constraints to apply to each of the parameters
  real(kind=8)               :: IC_const               ! (temporary) single value for the constraint on all IC parameters
  real(kind=8)               :: tmsc_const             ! (temporary) value for the constraint on tidal mascon parameters

! regularisation variables
  character*150              :: msc_WRMS_file          ! file of constraint values (in mm) to apply to mascons
  character*20               :: msc_const_model        ! type of model based on msc metrics
  integer*4                  :: ymd(3),doy             ! variables to compute decimal year
  real(kind=8)               :: dec_yr                 ! middle of the month of the first file in the input file list
  real(kind=8)               :: loose_const,tight_const! loosest and tightest constraints permitted on mascons

! scale parameter sigmas to read from the input command file
  real(kind=8)               :: scale_sigmas(3,2)
  real(kind=8)               :: bias_sigmas(3,2)
  character*1                :: satnam(2)
  integer*4                  :: isat
! variables to aid in the mass conservation constraints
  real(kind=8)  :: total_area                          ! total area of all primary mascons. Used to normalise the weighting by area
  real(kind=8)  :: mass_offset_corrn                   ! corrn to apply because the true answer does not conserve mass 
  real(kind=8)  :: scale_area   ,scale_mass            ! scale factor for area and mass to improve numerical stability
  real(kind=8)  :: sum_EWH

  integer*4     :: i,j,k,l,irec,irow,icol              ! counters
  integer*4     :: tide_row(2)          ! the row in the total tidal amplitude block that a particular tidal parameter relates to
  real(kind=8 ) :: tmp
  integer*4 trimlen
  integer*4     :: iflag                ! flag to  convert from ascii labels of parameters to parameter codes (1) or vice versa (2)
  integer*4     :: LUADD                ! unit number of ascii output file in VCV format
  integer*4     :: LUAtWA               ! unit number of ascii output file in AtWA format
  integer*4     :: LUAtWb               ! unit number of ascii output file in AtWb format

! RM210208: SVD related variables
  real(kind=8), allocatable  :: S(:)            ! contains the singular values
  real(kind=8), allocatable  :: U(:,:)          ! contains the left singular vectors
  real(kind=8), allocatable  :: VT(:,:)         ! contains the right singular values
  real(kind=8), allocatable  :: work(:)         ! contains the unconverged superdiagonal elements
  real(kind=8), allocatable  :: iwork(:)         ! contains the unconverged superdiagonal elements
  real(kind=8), allocatable  :: dummy(:)        ! dummy work array
  real(kind=8), allocatable  :: normeq_tmp(:,:) ! contains a copy of the AtWA
  real(kind=8)               :: CN         ! condition number
  integer*4                  :: lwork      ! dimensions work array
  integer*4                  :: info       ! info flag
  integer*4                  :: LUSVD      ! unit number of ascii output file containing S

! PT210217: variables to read the dates of stacked files from the snorm file. Dimension for 50 - there will never be that many days in a month!
! PT220509: with 2x12hr solutions there can be more than 50 days. Increase to 100
  integer*4    :: ndates,idate
  integer*4,parameter    :: max_epochs=10000
  integer*4    :: year_s(max_epochs),month_s(max_epochs),day_s(max_epochs),hr_s(max_epochs),min_s(max_epochs),sec_s(max_epochs)
  character*100:: normeq_files_s(max_epochs)                                              ! C*100, has to match the size of normeq_files defined above
  real(kind=8) :: duration_s(max_epochs)
  integer*4    :: i_IC,counter
  character*1  :: tmpsat

! PT140507: dummy variable to read unwanted AtWb values when doing a mascon-only inversion
! PT150828: increased to 30+4556 to account for 1/rev and 2/rev if needed
  real(kind=8) :: junk(30+4556)

! PT181216: declare variables for the correlation matrix
  real(kind=8),allocatable :: corr(:,:)

! PT150722: error flag
  integer*4 :: ioerr

! DEBUG SVD
  real(kind=8), allocatable  :: US(:,:),S_diag(:,:)
  real(kind=8), allocatable  :: SVD(:,:)
  real(kind=8), allocatable  :: X(:), Y(:)


! PT211206: variables to change how the command file is read
  real(kind=8) :: tmp_R8
  integer*4    :: lcmd
  
! PT220608: variables for computing correlation-dependent regularisation
  real(kind=8) :: ave_sigma
  integer*4    :: nmasc_correlated
  
! PT220620: variables related to adaptive regularisation
  real(kind=8 ), allocatable :: adaptive_adjust(:)    ! LS solution of (AtWA)^-1 * AtWb for first adaptive reg solution (2 cm Tikhonov)
  real(kind=8) :: delta_adj,adapt_scale_fact
  
! ****************************************************


! ****************************************************
! read the command line arguments
!
! command file
  call getarg(1,cmdfile)
  if(cmdfile(1:1) == ' ')then
     write(message,'(a,a)')'Runstring: addnorm addnorm.cmd mascons_stage5_V006_200km output_stem '&
                         ,' [msc_WRMS_file msc_const_model]'
     call status_update('FATAL','ADDNORM','addnorm',' ',message,0)
  endif
! PT211206: removed gracesim command file as the second argument. Read all info from the addnorm command file.

! mascon file name
  call getarg(2,mascon_file)
! output file name stem
  call getarg(3,outfile_stem)

  msc_WRMS_file = ""
  msc_const_model = ""

!! PT211029: conservation of mass constraint
!  call getarg(5, arg)
!  read(arg,*)cons_mass_sigma
!
!! mascon constraint?
!  call getarg(6,msc_WRMS_file)
!
!! mascon constraint model?
!  if(msc_WRMS_file(1:1) /= " ")then
!    call getarg(7,msc_const_model)
!  endif

! what else?
! ****************************************************



! ****************************************************
! open the command file
  call output_openFILE(lucmd,cmdfile,'old','ADDNORM','FATAL','Error opening command file','sequential',0)
! read all the command file information
  call read_addnorm_cmd(lucmd,IC_solve,msc_solve,tmsc_solve,output_snorm,output_SVD,output_ascii,msc_scl_factor,cons_mass_sigma &
                       ,mass_offset,msc_param_const,msc_constraint_file, tikhonov,lambda,nfiles,normeq_files )
! ****************************************************


! ****************************************************
! initialize the mascon_code to be XXXXXX
  mascon_code = "XXXXXX"

! allocate the arrays of the number of files and their record lengths, and numbers of parameters
!  allocate(normeq_files(nfiles))    ! PT211206: now allocated and read in read_addnorm_cmd
  allocate(rec_len(nfiles))
  allocate(versions(nfiles))
  allocate(nparams(nfiles))
  allocate(nICs(nfiles))
  allocate(nmascons(nfiles))
  allocate(nmsc_tidal(nfiles))
! allocate the variables for date/duration of orbit integrations
  allocate(year(nfiles))
  allocate(month(nfiles))
  allocate(day(nfiles))
  allocate(hr(nfiles))
  allocate(minute(nfiles))
  allocate(sec(nfiles))
  allocate(duration(nfiles))
! variables that indicate which satellite is first in the binary input files
  allocate(sat1(nfiles))
  allocate(sat2(nfiles))


! ****************************************************

! ****************************************************
!! KS160701 open gracesim command file
!!  call input_openCommand(gracesim_cmdfile)
!!  call command_readSetup()

! ********************************************************
!
!        r e a d   m a s c o n    f i l e
!
! ********************************************************
! read the mascon header
    call read_msc_hdr(21,mascon_file,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                                ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
! PT/HMcQ 161201: updated to the ellipsoidal geometry version of the ternary macson pattern
  call tern_lat_bands_ell(ternary_lat_spacing/60.d0)
  call allocate_mascon_arrays

! read in the mascon information
!print*,'turned off reading mascon file'
  call read_mascon_file(21,mascon_file)
  close(21)
! calculate the sum of the area of all mascons
  total_area = sum(mcon_prim(:,4))
! ********************************************************



! ********************************************************
!
! R E A D   R E C O R D   1    O F   B I N A R Y   F I L E
!
! ********************************************************
! open and determine the record length of each binary file
  do i=1,nfiles
     ! we don't yet know the record length - it depends how many parameters are in the file. So we need to open it and read the first 8-byte value to be able to work out the number of parameters. Get the number of mascons while we are at it.
     call output_openFILE(lufile,trim(normeq_files(i)),'old','ADDNORM','FATAL','Error opening normal equation file','direct',24)
     read(lufile,rec=1)versions(i),nparams(i),nmascons(i),nICs(i),nmsc_tidal(i)
     ! PT140717: if we don't want to estimate tidal amplitudes then set this to zero. It's used later on to dimension the number of parameters etc
     if(.not. tmsc_solve)nmsc_tidal(i) = 0
     !  print*,"set the param numbers right *#&^$*Q. PT150721: set to only 24 ICs (we don't estimate 1/rev anymore)"
     !    nICs(i) = 28
     !    nmascons(i) = 4556
     ! print*,'msc_tidal', nmsc_tidal(i)
     ! check that it is version 1
     if(versions(i) /= 1.0 .and. versions(i) /= 2.0)then
        write(message,'(a,f6.2)')'Code not written for binary file version',versions(i)
        call status_update('FATAL','ADDNORM','addnorm',normeq_files(i),message,0)
     endif
     ! check that the number of mascons matches, so that they can be stacked properly
     if(i == 1)then
        nmascon_t = nmascons(1)
     else
        if(nmascons(i) /= nmascon_t)then
           write(message,'(a,a,a,i6,a,a,a,i6,a)')'number mascons on file ',normeq_files(i),' (',nmascons(i) & 
                ,') not = number in first file ',normeq_files(1),' (',nmascon_t,')'
           call status_update('FATAL','ADDNORM','addnorm',' ',message,0)
        endif
     endif

     ! PT140709: we need to read the bit-mapped mascon tidal amplitudes (record 5), so that we can work out how many parameters there are in the whole solution
     if(i==1)then
        allocate(mcon_tides(nmascon_t,nfiles))
        allocate(tidal_lookup(nmascon_t,max_msc_tide,2))
        allocate(prm_tidal(nmascon_t*max_msc_tide*2))
        allocate(tide_list_epoch(nmascon_t*max_msc_tide*2,2,nfiles))
     endif
     ! reopen 
     close(lufile)
     call output_openFILE(lufile,normeq_files(i),'old','ADDNORM','FATAL','Error opening normal equation file','direct',nparams(i)*8)
     read(lufile,rec=5)(mcon_tides(j,i),j=1,nmascons(i))
     !   print *,  mcon_tides
     ! create/add to a list of mascon tidal amplitudes that need to be estimated across all solutions (this gives us nmsc_tidal_t that we need to dimension the inversion)
     call use_tidal_msc(nmascon_t,max_msc_tide,nmsc_tidal_t,mcon_tides(:,i),tidal_lookup,prm_tidal,tide_list_epoch(:,:,i))
     !   print *, tidal_lookup,prm_tidal,tide_list_epoch(:,:,i)
     write(message,'(a,i6)')"Number of tidal amplitudes to estimate (sin and cos): ",nmsc_tidal_t*2
     call status_update('STATUS','ADDNORM','addnorm',normeq_files(i),message,0)
     close(lufile)
  enddo
! We have now read the first line of each file, so can now dimension the required problem.
! ****************************************************




! ****************************************************
! we can now allocate the space for our entire normal equations, based on how many mascon and IC parameters there are
!
! ** Mascon only solution **
  if(msc_solve .and. .not. IC_solve .and. .not. tmsc_solve) then
     nparam_t = nmascon_t
     iMascons = 1
     nIC_t = 0
     nmsc_tidal_t = 0
  endif

! ** IC only **
  if(.not. msc_solve .and. IC_solve .and. .not. tmsc_solve)then
     ! We want only ICs. Add up the number of IC parameters
     nparam_t = 0
     do i=1,nfiles
        nparam_t = nparam_t + nICs(i)
        print*,'i,nparam_t,nICs(i)',i,nparam_t,nICs(i),versions(i),nparams(i),nmascons(i),nmsc_tidal(i),nIcs(i)
     enddo
     nIC_t = nparam_t
     iIC = 1                ! ICs start at the beginning if we don't want to estimate mascons
     nmascon_t    = 0         ! we don't want to estimate mascons if it is an IC_only solution
     iMascons     = 0
     nmsc_tidal_t = 0
     iMsc_tidal   = 0
  endif

! ** Mascon + IC **
  if(msc_solve .and. IC_solve .and. .not. tmsc_solve) then
     ! We want (at least) ICs. Add up the number of IC parameters
     nparam_t = 0
     do i=1,nfiles
        nparam_t = nparam_t + nICs(i)
     enddo
     nIC_t = nparam_t
     nparam_t = nparam_t + nmascon_t
     iMascons     = 1
     iMsc_tidal   = 0
     nmsc_tidal_t = 0
     iIC = nmascon_t + 1                     ! ICs start after mascons if we estimate mascons but not msc_tidal
  endif

! ** Mascon + IC + tidal amplitudes **
  if(msc_solve .and. IC_solve .and. tmsc_solve)then
     ! we want everything
     nparam_t = 0
     do i=1,nfiles
        nparam_t = nparam_t + nICs(i)
     enddo
     nIC_t = nparam_t
     nparam_t = nparam_t + nmascon_t + nmsc_tidal_t*2  ! multiply by two because there is a sine and cosine parameter for each constituent
     iMascons = 1
     iMsc_tidal = nmascon_t + 1
     iIC = nmascon_t + nmsc_tidal_t*2 + 1        ! ICs start after mascons and tidal amplitudes
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! now we can allocate the least squares arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(normeq_all(nparam_t,nparam_t))
  allocate(AtWb(nparam_t))
  allocate(adjust(nparam_t))
  allocate(apr_prm(nparam_t))
  allocate(VCV(nparam_t,nparam_t))
  allocate(prmnam(nparam_t))
  write(message,'(a,i8,a)')'allocate prmnam_all with ',nparam_t,' elements'
  call status_update('STATUS','ADDNORM','addnorm',' ',message,0)
  allocate(prmnam_all(nparam_t))
  allocate(msc_const(nparam_t,nparam_t))
  ! ocean/land flag
  allocate(mcon_ocean(nmascon_t))
  ! parameter constraints
  allocate(param_const(nparam_t))

! ***********************************************************
! initialise some variables
  prmnam = " "
  normeq_all = 0.d0
  AtWb = 0.d0
  param_const = 0.d0
! ****************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT140714: transfer the tidal mascon labels into prmnam (if necessary)
  if(tmsc_solve)then
     !    print*,'transferring tmsc labels into prmnam'
     do i=1,nmsc_tidal_t*2
        prmnam(iMsc_tidal-1+i) = prm_tidal(i)
     enddo
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! *******************************************************************************************************************************************
! output some diagnostic information
  if(msc_solve)   write (message(1:18),'(a18)')'Estimate mascons. '
  if(IC_solve)    write (message(19:32),'(a14)')'Estimate ICs. '
  if(tmsc_solve)  write (message(33:67),'(a35)')'Estimate mascons tidal amplitudes. '
  write(message(69:204),'(4(a,i6))')'# parameters to estimate: ',nparam_t,'. # mascons to estimate: ',nmascon_t &
       ,' # ICs to estimate: ',nIC_t,' # mascon tidal amplitudes to estimate: ',nmsc_tidal_t*2
  call status_update('STATUS','ADDNORM','addnorm',' ',message,0)
! *******************************************************************************************************************************************


! ***********************************************************************************************************************************************



! now open each binary file properly and read the information
  do i=1,nfiles
    call output_openFILE(lufile,normeq_files(i),'old','ADDNORM','FATAL','Error opening normal equation file','direct',nparams(i)*8)

    ! allocate the temporary variable that reads in the normal equations
    allocate(tmpvals(nparams(i)))
    allocate(apr_tmp(nparams(i)))
    allocate(param_type(nparams(i)))

    ! ****************************************************
    !
    !  H E A D E R   R E A D   O F   B I N A R Y   F I L E
    !
    ! ****************************************************
    ! line 1 contains version, # parameters, date and duration of orbit integration
    read(lufile,rec=1)versions(i),nparams(i),nmascons(i),nICs(i),nmsc_tidal(i),year(i),month(i),day(i),hr(i),minute(i),sec(i) &
          ,duration(i)
    write(message,'(a,f3.1,4(a,i6),a,6i5,a,f7.3 )') &
            "File version: ",versions(i)," Num params: ",nparams(i)," Num mascons: ",nmascons(i)," Num ICs: "&
            ,nICs(i), " Num tidal mascons: ",nmsc_tidal(i)," Date: ",year(i),month(i),day(i),hr(i),minute(i),sec(i) &
          ," Duration: ",duration(i)
    call status_update('STATUS','ADDNORM','addnorm',' ',message,0)

    ! line 2 of version 1.0 and version 2.0 files are the a priori values of the parameters
    read(lufile,rec=2)(apr_tmp(j),j=1,nparams(i))

    if(versions(i) == 1.0)then
      if(msc_solve) apr_prm(iMascons:nmascon_t)      = apr_tmp(nICs(i)+1:nICs(i)+nmascons(i))  ! store the a priori mascon values
      if(IC_solve ) apr_prm(iIC:iIC+nICs(i)-1)       = apr_tmp(1:nICs(i))             ! store the a priori IC values
      if(tmsc_solve)apr_prm(iMsc_tidal:iMsc_tidal+nmsc_tidal(i)-1) = apr_tmp(nICs(i)+nmascons(i)+1:+nICs(i)+nmascons(i)+nmsc_tidal(i)) 
    else if (versions(i) == 2.0)then
      if(msc_solve) apr_prm(1:nmascons(i))           = apr_tmp(1:nmascons(i))  ! store the a priori mascon values
      if(IC_solve ) apr_prm(nmascons(i):nparam_t)   = apr_tmp(nmascons(i):nparam_t)             ! store the a priori IC values
      if(tmsc_solve)stop 'code not written to estimate tidal mascons from snorm files' 
    endif

    ! line 3 of version 1.0 files contains (coded) descriptors for the parameters
    irec=3
    if(nmascons(i) == 0) then                ! the parameters will be only ICs
       read(lufile,rec=irec)sat1(i),(param_type(j),j=1,nICs(i)/2),sat2(i),(param_type(j),j=nICs(i)/2+1,nICs(i))
    else if (nICs(i) == 0) then      ! the parameters will be only mascons and possibly tidal amplitudes (ie no satellite indicator first).
       read(lufile,rec=irec)(param_type(j),j=1,nmascons(i)+nmsc_tidal(i))
    else
       ! PT210216: the read is different for norm and snorm files
       if(versions(i) == 1.0)then
         read(lufile,rec=irec)sat1(i),(param_type(j),j=1,nICs(i)/2),sat2(i),(param_type(j),j=nICs(i)/2+1,nICs(i)) &
            ,(param_type(nICs(i)+j),j=1,nmascons(i)+nmsc_tidal(i))
         ! set the mascon offset to be the number of ICs. This will shift the parameter codes correctly to before the ICs 
         msc_offset = -1 * nICs(i)
         IC_offset = nmascons(i)
       else if (versions(i) == 2.0)then
         read(lufile,rec=irec)param_type(1:nmascons(i)),param_type(nmascons(i)+1:nparam_t),sat1(i),sat2(i)
         ! no offset required for moving mascons before ICs - they are already first in the snorm files
         msc_offset = 0
         IC_offset = 0
       endif
    endif

    ! convert the parameter codes into character labels to be written out in gracefit vcv format
    iflag = 2    ! flag to convert from parameter codes to ascii labels for output file
    call output_codes2labels(iflag,i,IC_solve,msc_solve,tmsc_solve,nparams(i),nmascons(i),nmsc_tidal_t,nICs(i) &
          ,param_type,sat1(i),sat2(i),IC_offset,msc_offset,prmnam)
   
! store them in the prmnam_all vector
    prmnam_all(1:nmascons(i)) = prmnam(1:nmascons(i))
    prmnam_all(nmascons(i)+sum(nICs(1:i-1))+1:nmascons(i)+sum(nICs(1:i))) = prmnam(nmascons(i)+1:nmascons(i)+nICs(i))
    

    ! PT190906: store the satellite names, taking them from the values in the first file read.
    if(i == 1)then
      satnam(1) = sat1(i)
      satnam(2) = sat2(i)
    endif


    ! PT210217: record 4 are the dates of the files used in a stacked norm file
    !    irec=4
    if(versions(i) == 2.0)then
      call status_update('STATUS','ADDNORM','addnorm',' ','Reading dates of stacked normal equations',0)
      irec = 4
      read(lufile,rec=irec)ndates,(year_s(idate),month_s(idate),day_s(idate),hr_s(idate),min_s(idate),sec_s(idate) &
                          ,duration_s(idate), normeq_files_s(idate),idate=1,ndates)
! DEBUG: check the file names in the .snorm file
do idate=1,ndates
  print*,idate," ",normeq_files_s(idate)
enddo

    else
      ! PT220901: set ndates to nfiles if it wasn't a snorm
      ndates = nfiles
    endif


    ! PT220527: if a scale factor of -99 is entered and a snorm file was read in then change the msc_scale_factor to be
    !           the number of files that made up the snorm file
    !if(versions(i) == 2.0 .and. msc_scl_factor == -99.d0)then
    !  write(message,'(a,i4,a)')'Setting scale factor to number of days in snorm file (',ndates,')'
    !  call status_update('STATUS','ADDNORM','addnorm',' ',message,0)
    !  ! PT220603: divide by 2 if it is GRACE-FO since 2x12-hr norm files used
    !  ! PT220901: don't do this. Herb's gadi solutions are 24hr
    !  if(year_s(1) > 2017.d0)then
    !    msc_scl_factor = dble(ndates)/1.d0
    !  else
    !    msc_scl_factor = dble(ndates)
    !  endif
    !endif

    ! PT220829: if a negative scale factor has been entered then we need to change the scale factor to be a multiple of the
    !           number of days times the absolute value of the scale factor value entered
    if(msc_scl_factor < 0.d0)then
      write(message,'(a,i4,a,f7.2)')'Setting scale factor to number of days (',ndates,') times value entered' &
                                    ,dabs(msc_scl_factor)
      call status_update('STATUS','ADDNORM','addnorm',' ',message,0)
      msc_scl_factor = dble(ndates)*dabs(msc_scl_factor)
    else
      msc_scl_factor = 1.d0
    endif

    ! PT210218: each set of ICs in a snorm file will have been given the same day number, being "1". Need to number them sequentially
    if(versions(i) == 2.0)then
      do idate = 1,ndates
        do j = 1,24
          write(prmnam_all(nmascons(i)+(idate-1)*24+j)(1:3),'(i3)') idate
        enddo
      enddo
    endif


    ! PT140710: record 5 is the bit-mapped tidal amplitudes. We've already read and stored this, so just skip it here.
    irec = 5


    ! ****************************************************
    !
    !  N O R M E Q   R E A D   O F   B I N A R Y   F I L E
    !
    ! ****************************************************
! PT210217: we now have two different reads: the ICs-then-MSC files of GRACEFIT/GRACESIM or the MCs-then-ICs of stacked normal equations from ADDNORM
  if(versions(i) == 2.0)then
    ! just read it in as it is
    irec = 6
    do k=1,nparams(i)
      read(lufile,rec=irec)normeq_all(k,:)
      irec = irec + 1
    enddo



  else if (versions(i) == 1.0)then
    ! now read in all the normal equations for this file. They are found stored in the file as ICs first, then mascons.
    irec=6
    call status_update('STATUS','ADDNORM','addnorm',normeq_files(i),'Reading input normal equations and AtWb',0)
    do k=1,nparams(i)
      read(lufile,rec=irec)(tmpvals(j),j=1,nparams(i))


      ! store as required in the big normal equation matrix
      !
      ! ** Mascon-only solution **
      if(msc_solve .and. .not. IC_solve .and. .not. tmsc_solve) then
        if (k > nICs(i) .and. k <= nICs(i)+nmascons(i))then   ! the mascon/IC elements refer to IC parameter, k, whilever k <= nICs(i)
          do j=1,nmascon_t
            normeq_all(iMascons-1+k-nICs(i),j) = normeq_all(iMascons-1+k-nICs(i),j) + tmpvals(nICs(i)+j)
          enddo
        endif
      endif

      ! ** IC-only solution **
      if (IC_solve .and. .not. msc_solve .and. .not. tmsc_solve) then
        if(k <= nICs(i)) then
          do j=1,nICs(i)
            normeq_all(iIC-1+k,iIC-1+j) = tmpvals(j)
          enddo
        endif
      endif

      ! ** Mascon + IC (with or without tidal mascon) **
      if (msc_solve .and. IC_solve ) then
        if (k <= nICs(i))then  
          ! the IC elements 
          do j=1,nICs(i)
            normeq_all(iIC-1+k,iIC-1+j) = tmpvals(j)
          enddo

          ! additionally, the mascon/IC elements that refer to IC parameter, k, whilever k <= nICs(i)
          do j=1,nmascon_t
            normeq_all(iMascons-1+j,iIC-1+k) = tmpvals(nICs(i)+j)
            ! PT/RM 210521: don't need this in the lower-triangular anymore
            ! normeq_all(iIC-1+k,iMascons-1+j) = tmpvals(nICs(i)+j)
          enddo
        else                    ! the mascon/mascon elements come below the bottom-right corner of the input IC normal equations
          if( k < nICs(i)+nmascons(i)+1)then
            do j=1,nmascon_t
              normeq_all(iMascons-1+k-nICs(i),j) = normeq_all(iMascons-1+k-nICs(i),j) + tmpvals(nICs(i)+j)
            enddo
          endif
        endif
      endif

      ! ** add the tidal amplitudes? **
      if(tmsc_solve) then
        if (k > nICs(i)+nmascons(i))then  ! we have a tidal amplitude parameter
          k_tidal = nint( (dble(k)+0.5-nICs(i)-nmascons(i))/2.d0 )     ! k_tidal is the nth tidal constituent (not double counting for sine and cosine)
          !
          !
          ! First, work our which row in the tidal amplitude block this parameter relates to
          tide_row(1) = tidal_lookup( tide_list_epoch(k_tidal,1,i),tide_list_epoch(k_tidal,2,i),2)*2-1
          ! distinguish between the sine and cosine parameter rows
          ! PT170301: fix this so that it no longer depends on there being an even number of mascons ....
          if(mod(dble(k-(nICs(i)+nmascons(i))),2.d0) == 0)then   ! it is the second row, therefore the cosine so add one row to the index
            tide_row(1) = tide_row(1) + 1
          endif

          !print*,'row',k,' is tidal amplitude ',tide_row(1),k_tidal

          ! store the block of ICs wrt this tidal amplitude.
          normeq_all(iMsc_tidal-1+tide_row(1),iIC:iIC+nICs(i)-1) = tmpvals(1:nICs(i))
          normeq_all(iIC:iIC+nICs(i)-1,iMsc_tidal-1+tide_row(1)) = tmpvals(1:nICs(i))

          ! stack the block of MCs wrt this tidal amplitude
          normeq_all(iMsc_tidal-1+tide_row(1),iMascons:iMascons+nmascon_t-1) = &
               normeq_all(iMsc_tidal-1+tide_row(1),iMascons:iMascons+nmascon_t-1) + tmpvals(nICs(i)+1:nICs(i)+nmascons(i))
          normeq_all(iMascons:iMascons+nmascon_t-1,iMsc_tidal-1+tide_row(1)) = &    
               normeq_all(iMascons:iMascons+nmascon_t-1,iMsc_tidal-1+tide_row(1)) + tmpvals(nICs(i)+1:nICs(i)+nmascons(i))

          ! now, stack the entries of this mascon/tidal amplitude combination against all other tidal amplitudes in this set of normal equations.
          ! This is just a matter of looking up the row for the second tidal amplitude, then stacking away the values
          do l = 1,nmsc_tidal(i)/2
            tide_row(2) = tidal_lookup( tide_list_epoch(l,1,i),tide_list_epoch(l,2,i),2)*2 - 1
            normeq_all(iMsc_tidal-1+tide_row(1),iMsc_tidal-1+tide_row(2))   =  &
                      normeq_all(iMsc_tidal-1+tide_row(1),iMsc_tidal-1+tide_row(2)) + tmpvals(nICs(i)+nmascon_t+tide_row(2))
            normeq_all(iMsc_tidal-1+tide_row(1),iMsc_tidal-1+tide_row(2)+1) = &
                      normeq_all(iMsc_tidal-1+tide_row(1),iMsc_tidal-1+tide_row(2)+1) + tmpvals(nICs(i)+nmascon_t+tide_row(2)+1)
          enddo
        endif
      endif

      ! increment the record number for the binary file
      irec = irec + 1

    enddo    ! end of reading all the lines of this particular input normal equations file

  endif ! end of if to read version 1.0 or version 2.0 input files (.norm or .snorm)
! **********************************************************************************************************************************


! ***********
! 
!     RHS
!
! ***********
! now read in the RHS from this file. It simply follows on from the normal equations in the binary file - all on a single line.
  if(versions(i) == 2.0)then
    read(lufile,rec=irec)AtWb(:)

  else if (versions(i) == 1.0)then
    if(IC_solve)then
      ! read the IC values
      read(lufile,rec=irec)(AtWb(iIC-1+j),j=1,nICs(i) )
    endif
    if(msc_solve)then
      ! the mascon AtWb comes after the IC values
      read(lufile,rec=irec)(junk(j),j=1,nICs(i)),(tmpvals(j),j=1,nmascons(i))
      do j=1,nmascon_t
        AtWb(iMascons-1+j) = AtWb(iMascons-1+j) + tmpvals(j)
      enddo
    endif
    if(tmsc_solve)then
      read(lufile,rec=irec)(junk(j),j=1,nICs(i)+nmascons(i)),(tmpvals(j),j=1,nmsc_tidal(i))  ! PT170301: needed to multiply this by two ?
      ! PT140914: we need some logic here to ensure that the right value is put into the right row
      !  tide_list_epoch: array to give us mascson and tide for a given nmsc_tidal entry
      !  tidal_lookup   : array to give us row number in the entire list, given mascon and tide
      do j=1,nmsc_tidal(i)
        if(mod(j,2) > 0)then
          tide_row(1) = tidal_lookup(tide_list_epoch(nint((j+0.5)/2),1,i),tide_list_epoch(nint((j+0.5)/2),2,i),2)*2 - 1
          AtWb(iMsc_tidal-1+tide_row(1)) = AtWb(iMsc_tidal-1+tide_row(1)) + tmpvals(j)
          AtWb(iMsc_tidal-1+tide_row(1)+1) = AtWb(iMsc_tidal-1+tide_row(1)+1) + tmpvals(j+1)
        endif
      enddo
    endif

  endif
    ! ***********

    ! now we increment the pointer for the location of the ICs in the combined normal equations
    iIC = iIC + nICs(i)

    ! deallocate the temporary arrays used to read the input binary files
    deallocate(tmpvals)
    deallocate(apr_tmp)
    deallocate(param_type)

  enddo   ! end of loop over input binary files
! ****************************************************

! ****************************************************
! PT210216: write out the stacked normal equations - mascons first - into a binary file
  if(output_snorm)then
    write(message,'(a,a)')'Writing out binary version of stacked AtWA to ' &
         ,trim(outfile_stem)//".snorm ('snorm' = 'stacked norm' .... I didn't choose the name!)"
    call status_update('STATUS','ADDNORM','addnorm',' ',message,0)
    LUAtWA = 21
    rec_len = nparam_t*8
    open(LUAtWA,file=trim(outfile_stem)//'.snorm',status='unknown',access='direct',form='unformatted',recl=rec_len(1))

    ! make these snorm files version 2.0
    irec = 1
    write(LUAtWA,rec=irec)2.d0,nparam_t,nmascon_t,nparam_t-nmascon_t-nmsc_tidal_t*2,nmsc_tidal_t*2 &
       ,year(1),month(1),day(1),hr(1),minute(1),sec(1) ,duration(1) 

    ! a priori parameter values. Mascons first, then ICs
    irec=2
    write(LUAtWA,rec=irec)(apr_prm(i),i=1,nparam_t)

    ! Record 3: coded info related to parameter types. 
    irec=3
    ! convert the parameter codes into character labels to be written out in gracefit vcv format
    IC_solve   = .false.
    msc_solve  = .false.
    tmsc_solve = .false.
    if(nIC_t > 0  ) IC_solve   = .true.
    if(nmascon_t > 0 ) msc_solve  = .true.
    if(nmsc_tidal_t > 0) tmsc_solve = .true.

    ! PT210216: we need a new character array to store the parameter codes for all the mascons + ICs for each day
    allocate(prmnam_allICs(nparam_t))
    allocate(param_type(nparam_t))
    prmnam_allICs = " "
    prmnam_allICs(1:nmascon_t)(1:29) = prmnam(1:nmascon_t)(2:30)  ! mascon labels
    do j=1,nfiles
      prmnam_allICs(nmascon_t+(j-1)*24+1:nmascon_t+(j)*24)(1:27) = prmnam(nmascon_t+1:nmascon_t+24)(4:30)
    enddo
    iflag = 1    ! flag to convert from ascii labels to parameter codes for output file
    call output_codes2labels(iflag,1,IC_solve,msc_solve,tmsc_solve,nparam_t,nmascon_t,nmsc_tidal_t,nparam_t-nmascon_t &
       ,param_type,sat1,sat2,1,0,prmnam_allICs)
    write(LUAtWA,rec=irec)param_type(1:nparam_t),sat1(1),sat2(1)

    ! record 4 is now going to be a list of the dates of all the stacked files. Number of files first, then yr/mo/day/hr/min for each
    ! PT210706: also write out the file names of the daily .norm files that have been stacked.
    irec = 4
    write(LUAtWA,rec=irec)nfiles,  ( year(i),month(i),day(i),hr(i),minute(i),sec(i) ,duration(i), normeq_files(i) ,i=1,nfiles)

    ! record 5: bit-mapped tidal amplitudes that were estimated
    ! PT170611: we write out the number of ocean mascons, not nmascons_t
    irec = 5
    write(LUAtWA,rec=irec)(mcon_tides(i,1),i=1,total_ocean_prim_ampl)

    !************* WRITE normal equations to normeq FILE ***********
    do i= 1, nparam_t
      irec=irec+1
      write(LUAtWA,rec=irec)(normeq_all(i,j),j=1,nparam_t)
    enddo
    !****************************************************************

    !******************* WRITE RHS to normeq FILE *******************
    write(message,'(a)')'Writing out stacked (AtWb)^t ' 
    call status_update('STATUS','ADDNORM','addnorm',trim(outfile_stem)//".snorm",message,0)
    irec=irec+1
    write(LUAtWA,rec=irec)(AtWb(i),i=1,nparam_t)
    !****************************************************************
    close (LUAtWA)

  endif  ! end of if statement to output binary stacked normal equations file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  O U T P U T    T H E    S T A C K E D    N O R M    E Q S
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RM201020: if output_norm is true write out the stacked norm equations in ascii format
  if(output_ascii)then
    write(message,'(a,a)')'Writing out ascii version of combined normal equations to ',trim(outfile_stem)//'.AtWA'
    call status_update('STATUS','ADDNORM','addnorm',' ',message,0)
    LUAtWA = 21
    open(LUAtWA,file=trim(outfile_stem)//'.AtWA',status='unknown')
    ! PT180828: hardwire some values into the first line
    write(LUAtWA,*)1.0,nparam_t,nmascon_t,nparam_t-nmascon_t
    do i=1,nparam_t
      write(LUAtWA,*)(normeq_all(i,j),j=1,nparam_t)
    enddo
    close(LUAtWA)
  endif
! RM201119: write out AtWb too
  if(output_ascii)then
    write(message,'(a,a)')'Writing out ascii version of combined AtWb to ',trim(outfile_stem)//'.AtWb'
    call status_update('STATUS','ADDNORM','addnorm',' ',message,0)
    LUAtWb = 22
    open(LUAtWb,file=trim(outfile_stem)//'.AtWb',status='unknown')
    ! PT180828: hardwire some values into the first line
    write(LUAtWb,*)1.0,nparam_t,nmascon_t,nparam_t-nmascon_t
    do i=1,nparam_t
      write(LUAtWb,*)(AtWb(i))
    enddo
    close(LUAtWb)
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ********************************
!  Unregularised correlation matrix
!  ********************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(use_correl_in_reg)then
    call status_update('STATUS','ADDNORM','addnorm',normeq_files(1),"Compute unregularised VCV",0)
    allocate(normeq_tmp(nparam_t,nparam_t))
    normeq_tmp = normeq_all
    call Chol_normSolve(nparam_t,normeq_tmp,AtWb,adjust)
    VCV = normeq_tmp
    deallocate(normeq_tmp)
    
    call status_update('STATUS','ADDNORM','addnorm',' ',"Compute unregularised correlation matrix",0)
    allocate(corr(nparam_t,nparam_t))
    do irow = 1, nparam_t
       do icol = 1, nparam_t
         corr(irow,icol) = VCV(irow,icol)/( dsqrt(VCV(irow,irow))*dsqrt(VCV(icol,icol)) ) 
       enddo
    enddo
! DEBUG
!print*,"corr(4145:4155): ",corr(4150,4145:4155)
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ********************************
!  Adaptive regularisation : first inversion
!  ********************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(use_adaptive_reg)then
    write(message,'(a,f8.4,a)')"Compute 1st adaptive reg soln with sigma",adapt_sigma,' m'
    adapt_sigma = 0.03d0
    call status_update('STATUS','ADDNORM','addnorm',normeq_files(1),message,0)
    allocate(normeq_tmp(nparam_t,nparam_t))
    normeq_tmp = normeq_all
    ! add Tikhonov regularisation of n/2cm^^2
    do irow=1,nmascon_t
      normeq_tmp(irow,irow) = normeq_tmp(irow,irow) + ndates/adapt_sigma**2
    enddo
    
    ! save off the first adjustment here
    allocate(adaptive_adjust(nparam_t))
    call Chol_normSolve(nparam_t,normeq_tmp,AtWb,adaptive_adjust)
    deallocate(normeq_tmp)

! debug
    open(506,file=trim(outfile_stem)//'.adaptive_1st',status='unknown')   ! 1 header line, then param number+full correlation matrix
    write(506,*)'1st adaptive regularisation Tikhonov soln (2cm sigma)'

    do i=1,nmascon_t
      write(506,*)mcon_prim(i,1:2)/(4.0*atan(1.d0)/180.d0),adaptive_adjust(i)
    enddo
    close(506)    


  endif



    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  A D D    C O N S T R A I N T S
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialise the constraint matrix to zero
  msc_const = 0.d0

! ****************************************************
! KS160701 get the mascon constraint from the gracesim command file
! PT211206: move this to read_addnorm_cmd
!  call command_storeCommand(10,'apr_mascons',msc_param_const,1,1,'Mascon Constraint',0)
  if(nmascon_t > 0 .and. msc_param_const /= 0.d0)then
     write(message,'(a,e15.6,a)')'Adding',msc_param_const,' m constraint to mascon parameters'
     call status_update('STATUS','ADDNORM','addnorm',' ',message,0)
     do i=1,nmascon_t
        param_const(i) = 1.d0/msc_param_const**2
     enddo


  else
     call status_update('STATUS','ADDNORM','addnorm',' ',"NOT adding mascon constraints",0) 
  endif

! satellite pos/vel
  IC_const = 1.d-0
  if(nIC_t > 0 .and. IC_const /= 0.d0)then
     write(message,'(a,e15.6,a)')'Adding',IC_const,' m constraint to IC parameters'
     call status_update('STATUS','ADDNORM','addnorm',' ',message,0)
     !iIC = nmascon_t+nmsc_tidal_t*2+1
     do i=1,nfiles
        iIC = nmascon_t+nmsc_tidal_t*2 + (i-1)*nICs(i)+1
        do j=1,nICs(i)
           param_const(iIC -1 + j) = 1.d0/IC_const**2
        enddo
     enddo
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! constrain the scales and biases
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT190906: get the along-track, cross-track and radial scale constraints from the gracefit/gracesim command file
  do isat = 1, 2
     do i = 1,3
        call command_storeCommand(lucmd,'apr_sat_scale_'//satnam(isat),scale_sigmas(i,isat),i, &
             3,'GRACE '//satnam(isat)//' acceleration scale a priori constraints',1)
        call command_storeCommand(lucmd,'apr_sat_bias_'//satnam(isat),bias_sigmas(i,isat),i, &
             3,'GRACE '//satnam(isat)//' acceleration bias a priori constraints',1)
     enddo
  enddo

  if(nIC_t > 0 )then
    write(message,'(a,a1,a,3e15.4)')'Adding scale constraints satellite:',satnam(1),':',scale_sigmas(:,1)
    call status_update('STATUS','ADDNORM','addnorm',' ',message,0)
    write(message,'(a,a1,a,3e15.4)')'Adding scale constraints satellite:',satnam(2),':',scale_sigmas(:,2)
    call status_update('STATUS','ADDNORM','addnorm',' ',message,0)
    write(message,'(a,a1,a,3e15.4)')'Adding bias  constraints satellite:',satnam(1),':',bias_sigmas(:,1)
    call status_update('STATUS','ADDNORM','addnorm',' ',message,0)
    write(message,'(a,a1,a,3e15.4)')'Adding bias  constraints satellite:',satnam(2),':',bias_sigmas(:,2)
    call status_update('STATUS','ADDNORM','addnorm',' ',message,0)
    iIC = nmascon_t+nmsc_tidal_t*2+1

    do i=1,nfiles

! PT210222: redo the IC constraint logic so that it works for snorm files as well as multiple norm files
      do j=iIC,iIC+nICs(i)-1
  
        ! do a character check on the prmnam_all code to identify the scale and bias parameters

        ! along-track scale
        if(prmnam_all(j)(16:19) == "sclx" .and. prmnam_all(j)(12:12) == sat1(i))param_const(j)     = 1.d0/scale_sigmas(1,1)**2
        if(prmnam_all(j)(16:19) == "sclx" .and. prmnam_all(j)(12:12) == sat2(i))param_const(j)     = 1.d0/scale_sigmas(1,2)**2
        ! along-track bias
        if(prmnam_all(j)(16:18) == "bsx" .and. prmnam_all(j)(12:12) == sat1(i))param_const(j)     = 1.d0/bias_sigmas(1,1)**2
        if(prmnam_all(j)(16:18) == "bsx" .and. prmnam_all(j)(12:12) == sat2(i))param_const(j)     = 1.d0/bias_sigmas(1,2)**2

        ! cross-track scale
        if(prmnam_all(j)(16:19) == "scly" .and. prmnam_all(j)(12:12) == sat1(i))param_const(j)     = 1.d0/scale_sigmas(2,1)**2
        if(prmnam_all(j)(16:19) == "scly" .and. prmnam_all(j)(12:12) == sat2(i))param_const(j)     = 1.d0/scale_sigmas(2,2)**2
        ! cross-track bias
        if(prmnam_all(j)(16:18) == "bsy" .and. prmnam_all(j)(12:12) == sat1(i))param_const(j)     = 1.d0/bias_sigmas(2,1)**2
        if(prmnam_all(j)(16:18) == "bsy" .and. prmnam_all(j)(12:12) == sat2(i))param_const(j)     = 1.d0/bias_sigmas(2,2)**2

        ! radial scale
        if(prmnam_all(j)(16:19) == "sclz" .and. prmnam_all(j)(12:12) == sat1(i))param_const(j)     = 1.d0/scale_sigmas(3,1)**2
        if(prmnam_all(j)(16:19) == "sclz" .and. prmnam_all(j)(12:12) == sat2(i))param_const(j)     = 1.d0/scale_sigmas(3,2)**2
        ! radial bias
        if(prmnam_all(j)(16:18) == "bsz" .and. prmnam_all(j)(12:12) == sat1(i))param_const(j)     = 1.d0/bias_sigmas(3,1)**2
        if(prmnam_all(j)(16:18) == "bsz" .and. prmnam_all(j)(12:12) == sat2(i))param_const(j)     = 1.d0/bias_sigmas(3,2)**2

      enddo
      iIC = iIC + nICs(i)

    enddo
 
  endif

! tidal mascon amplitudes
  tmsc_const = 1.d-6
  if(nmsc_tidal_t > 0 .and. tmsc_const /= 0.d0)then
     do i=1,nmsc_tidal_t*2
        if(    prmnam(nmascon_t+i)(16:18) == "M2s")then
           call command_storeCommand(10,'apr_M2_ampl',tmsc_const,1,2,'M2 sin tide amplitude constraint',0)
           write(message,'(a,e10.4)')"Constraint on M2 sin: ",tmsc_const
           call status_update('STATUS','ADDNORM','addnorm',' ',message,0)
           param_const(nmascon_t+i) = 1.0/tmsc_const**2
        elseif(prmnam(nmascon_t+i)(16:18) == "M2c")then
           call command_storeCommand(10,'apr_M2_ampl',tmsc_const,2,2,'M2 cos tide amplitude constraint',0)
           write(message,'(a,e10.4)')"Constraint on M2 cos: ",tmsc_const
           call status_update('STATUS','ADDNORM','addnorm',' ',message,0)
           param_const(nmascon_t+i) = 1.0/tmsc_const**2
        elseif(prmnam(nmascon_t+i)(16:18) == "O1s")then
           call command_storeCommand(10,'apr_O1_ampl',tmsc_const,1,2,'O1 sin tide amplitude constraint',0)
           write(message,'(a,e10.4)')"Constraint on O1 sin: ",tmsc_const
           call status_update('STATUS','ADDNORM','addnorm',' ',message,0)
           param_const(nmascon_t+i) = 1.0/tmsc_const**2
        elseif(prmnam(nmascon_t+i)(16:18) == "O1c")then
           call command_storeCommand(10,'apr_O1_ampl',tmsc_const,2,2,'O1 cos tide amplitude constraint',0)
           write(message,'(a,e10.4)')"Constraint on O1 cos: ",tmsc_const
           call status_update('STATUS','ADDNORM','addnorm',' ',message,0)
           param_const(nmascon_t+i) = 1.0/tmsc_const**2
        elseif(prmnam(nmascon_t+i)(16:18) == "S2s")then
           call command_storeCommand(10,'apr_S2_ampl',tmsc_const,1,2,'S2 sin tide amplitude constraint',0)
           write(message,'(a,e10.4)')"Constraint on S2 sin: ",tmsc_const
           call status_update('STATUS','ADDNORM','addnorm',' ',message,0)
           param_const(nmascon_t+i) = 1.0/tmsc_const**2
        elseif(prmnam(nmascon_t+i)(16:18) == "S2c")then
           call command_storeCommand(10,'apr_S2_ampl',tmsc_const,2,2,'S2 cos tide amplitude constraint',0)
           write(message,'(a,e10.4)')"Constraint on S2 cos: ",tmsc_const
           call status_update('STATUS','ADDNORM','addnorm',' ',message,0)
           param_const(nmascon_t+i) = 1.0/tmsc_const**2
        elseif(prmnam(nmascon_t+i)(16:18) == "K1s")then
           call command_storeCommand(10,'apr_K1_ampl',tmsc_const,1,2,'K1 sin tide amplitude constraint',0)
           write(message,'(a,e10.4)')"Constraint on K1 sin: ",tmsc_const
           call status_update('STATUS','ADDNORM','addnorm',' ',message,0)
           param_const(nmascon_t+i) = 1.0/tmsc_const**2
        elseif(prmnam(nmascon_t+i)(16:18) == "K1c")then
           call command_storeCommand(10,'apr_K1_ampl',tmsc_const,2,2,'K1 cos tide amplitude constraint',0)
           write(message,'(a,e10.4)')"Constraint on K1 cos: ",tmsc_const
           call status_update('STATUS','ADDNORM','addnorm',' ',message,0)
           param_const(nmascon_t+i) = 1.0/tmsc_const**2
        elseif(prmnam(nmascon_t+i)(16:18) == "K2s")then
           call command_storeCommand(10,'apr_K2_ampl',tmsc_const,1,2,'K2 sin tide amplitude constraint',0)
           write(message,'(a,e10.4)')"Constraint on K2 sin: ",tmsc_const
           call status_update('STATUS','ADDNORM','addnorm',' ',message,0)
           param_const(nmascon_t+i) = 1.0/tmsc_const**2
        elseif(prmnam(nmascon_t+i)(16:18) == "K2c")then
           call command_storeCommand(10,'apr_K2_ampl',tmsc_const,2,2,'K2 cos tide amplitude constraint',0)
           write(message,'(a,e10.4)')"Constraint on K2 cos: ",tmsc_const
           call status_update('STATUS','ADDNORM','addnorm',' ',message,0)
           param_const(nmascon_t+i) = 1.0/tmsc_const**2
        endif
     enddo

  else
     call status_update('STATUS','ADDNORM','addnorm',' ',"Not estimating tidal amplitudes",0)
  endif

!****************************************************************


! ****************************************************************
! Add mascon length-dependent constraint
!    call status_update('STATUS','ADDNORM','addnorm',' ',"**NOT ** Adding length-dependent constraints to mascons",0)
  if(.not.tikhonov .and. msc_WRMS_file(1:1) == "")then
    write(message,'(a,f6.2)')"Adding spatial constraints to mascons. Scale factor used: ",msc_scl_factor
    call status_update('STATUS','ADDNORM','addnorm',msc_constraint_file,message,0)
    ioerr = 0
    open(990,file=trim(msc_constraint_file),status='old',iostat=ioerr)
    if(ioerr /= 0)then
       call status_update('FATAL','ADDNORM','addnorm',msc_constraint_file,'Error opening file. Does it exist?',0)
    endif

    ! PT171031: extract the mascon code from the first header line
    ! RM201125: read entire line to determine if file is diag_only
    read(990,'(a)')line
    read(line,'(1x,a6,53x,a9)')mascon_code,diag_only

    ! PT161130: loop over all header lines
    message(1:1) = "#"
    do while (message(1:1) == "#")
       ! read a header line
       read(990,'(a)')message
       ! PT150721: check that it really was a header line
       if(message(1:1) == "#")then
          call status_update('STATUS','ADDNORM','addnorm',' ',message,0)
       else
          backspace(990)
       endif
    enddo

    ! RM201125: read in entire reg file if not diag_only
    if(diag_only /= "diag_only")then
      do i=1,nmascon_t
         read(990,*)(msc_const(i,j),j=1,nmascon_t)
      enddo
    else
    ! RM201125: else just fill in diag elements from vector reg file
      do i=1,nmascon_t
         read(990,*)msc_const(i,i)
      enddo
    endif
    close(990)
  else if (msc_WRMS_file(1:1) /= "")then

! PT200122: we now have multiple options for what to read/use, so put it in a subroutine.
!           mascon constraint is added to diagonal elements within the subroutine
!           Models coded: landWRMS, WRMS, LOWFRQ_ONLY, LOWFRQ+ANNUAL
!
! PT200131: the model value is now set through a command line argument
    ymd(1) = year(1)
    ymd(2) = month(1)
    ymd(3) = 15
    call ymd_to_doy(ymd,doy)
    dec_yr = dble(year(1)) + dble(doy)/365.d0
    call read_msc_ts_wrms(luWRMS,msc_WRMS_file,msc_const_model,dec_yr,nparam_t,nmascon_t,msc_const)

! PT200203: set minimum and maximum values for the macson constraints
    loose_const = 0.3d0  ! 0.3 m is the loosest constraint we allow
    tight_const = 0.03d0 ! 3 cm  is the tightest constraint we allow
    do i=1,nmascon_t
      if(msc_const(i,i) < 1.d0/loose_const**2 )then
        msc_const(i,i) = 1.d0/loose_const**2
      else if (msc_const(i,i) > 1.d0/tight_const**2 )then
        msc_const(i,i) = 1.d0/tight_const**2
      endif
    enddo

  else
    ! PT180618: set up the Tikhonov matrix, with the input lambda value down the diagonal
    write(message,'(a,f10.2)')"Adding Tikhonov constraints to mascons. Diagonal used: ",lambda
    call status_update('STATUS','ADDNORM','addnorm',msc_constraint_file,message,0)
    do i=1,nmascon_t
       msc_const(i,i) = lambda
    enddo
  endif

  ! PT220620: upscale some of the diagonal elements if we want adaptive regularisation
  if(use_adaptive_reg)then
    do i=1,nmascon_t
      ! PT220620: upscale the constraints if we want adaptive regularisation, but only for continental mascons
      if( mcon_prim(i,6) < 1010.d0)then
        ! PT220621: we want a scale factor that changes linearly from 1.0 for zero adjustments to a maximum value 
        !           of 3.0 for abs(adjustments) > 0.5
        delta_adj = dabs(adaptive_adjust(i))
        if(delta_adj > 0.5d0)delta_adj = 0.5d0
        adapt_scale_fact = 1.d0 + (delta_adj/0.5d0) * 2.d0
        !print*,"Mascon ",i," with adjustment ",adaptive_adjust(i)," should upscale by ",adapt_scale_fact
        msc_const(i,i) = msc_const(i,i)/adapt_scale_fact
      endif
    enddo
  endif

  ! add in the parameter constraints
  ! PT220608: add them into the msc_const matrix rather than directly into the normal equations
  do i=1,nparam_t
     normeq_all(i,i) = normeq_all(i,i) + param_const(i)
     !msc_const(i,i) = msc_const(i,i) + param_const(i)
  enddo

  ! PT220608: if requested, include the off-diagonal terms based on mascon correlations from the unregularised inversion
  if(use_correl_in_reg)then
    call status_update('STATUS','ADDNORM','addnorm',' ',"Including off-diag terms in regularisation, based on mascon correlations",0)
    do i=imascons,imascons+nmascon_t-1
! DEBUG
!j=1
!print*,"i,j,dsqrt(1/msc_const(i,i))",i,j,1.d0/dsqrt(msc_const(i-imascons+1,i-imascons+1)) &
!      ,1.d0/dsqrt(msc_const(j-imascons+1,j-imascons+1))

      do j=imascons,imascons+nmascon_t-1
        if(dabs(corr(i-imascons+1,j-imascons+1)) > 0.7 .and. i /= j .and. mcon_prim(i,6) == mcon_prim(j,6))then
          nmasc_correlated = nmasc_correlated + 1    
          ! find the average sigma constraint of the two mascons
          ave_sigma = ( 1.d0/dsqrt(msc_const(i-imascons+1,i-imascons+1)) + 1.d0/dsqrt(msc_const(j-imascons+1,j-imascons+1)) ) /2.d0
! DEBUG
!print*,"i,j,dsqrt(1/msc_const(i,i))",i,j,1.d0/dsqrt(msc_const(i-imascons+1,i-imascons+1)) &
!      ,1.d0/dsqrt(msc_const(j-imascons+1,j-imascons+1)),ave_sigma


          msc_const(i-imascons+1,j-imascons+1) = msc_const(i-imascons+1,j-imascons+1) &
             + (-0.21d0)* corr(i-imascons+1,j-imascons+1) *1.d0/ave_sigma**2 !     *400.d0 !* (msc_const(i-imascons+1,i-imascons+1) + msc_const(j-imascons+1,j-imascons+1))/2.d0


        endif
      enddo
    enddo 
    write(message,'(a,i6,a)')"A total of ",nmasc_correlated," mascon pairs had high correlations"
    call status_update('STATUS','ADDNORM','addnorm',normeq_files(1),message,0)

    deallocate(corr)
  endif
  

  if( (diag_only == "diag_only" .or. tikhonov) .and. .not. use_correl_in_reg )then
    do i=imascons,imascons+nmascon_t-1
      normeq_all(i,i) = normeq_all(i,i) + msc_scl_factor * msc_const(i-imascons+1,i-imascons+1)
    enddo
  else
    do i=imascons,imascons+nmascon_t-1
      do j=imascons,imascons+nmascon_t-1
        normeq_all(i,j) = normeq_all(i,j) + msc_scl_factor * msc_const(i-imascons+1,j-imascons+1)
      enddo
    enddo 
  endif

!****************************************************************
! PT150722: Add mass conserving constraints on mascons
! PT150722: this code works ONLY if the a priori mascon values are zero. If not, we MUST add stuff to AtWb as well !!!!
  !cons_mass_sigma = 1.5d+14
  ! PT211029: just for testing, set msc_scl_factor to 1.0 for adding the CoM constraint
  msc_scl_factor = 1.d0
  mass_offset = 1.11652531846531015625e+3   ! PT211101: this is the amount (in kg) by which Bec's simulations do not coserve mass when using a zero ocean.
  mass_offset = 0.d0
  scale_area = 1.d-8
  scale_mass = 1.d-11 
! sum_EWH = sum(apr_prm(imascons:imascons+nmascon_t-1))
  sum_EWH = 330.8

  if(cons_mass_sigma > 0.0) then
     write(message,'(a,e15.6,a)')"Applying mass conservation with sigma of ",cons_mass_sigma," m."
     call status_update('STATUS','ADDNORM','addnorm',' ',message,0)

     if( dabs(sum_EWH) >0.05d0)then
      mass_offset_corrn = mass_offset/sum_EWH
     else
       mass_offset_corrn = 0.d0
     endif

     do i=imascons,imascons+nmascon_t-1

        do j=i,imascons+nmascon_t-1
           ! PT200112: we now scale this by the areas of both mascons (in km^2). Also, once per file (so multiply by "nfiles").
           ! PT210222: instead of number of files, multiply by msc_scl_factor
           ! PT211029: as written originally, this is a volume constraint. We need to use the densities as well.

!if(i<10 .and. j < 10)write(*,'(2i7,9e15.6)')i,j &
!             ,(mcon_prim(i,6) * mcon_prim(i,4) - mass_offset_corrn) &
!                               * (mcon_prim(j,6) * mcon_prim(j,4) - mass_offset_corrn)/ cons_mass_sigma**2 &
!             , normeq_all(i,j) &
!             , mass_offset_corrn &
!             , mcon_prim(i,6) , mcon_prim(i,4),mass_offset_corrn &
!             , (mcon_prim(i,6) * mcon_prim(i,4) - mass_offset_corrn) &
!             , (mcon_prim(j,6) * mcon_prim(j,4) - mass_offset_corrn) &
!             ,  mcon_prim(i,6) * mcon_prim(i,4)

 !           normeq_all(i,j) = normeq_all(i,j) + (scale_mass*(mcon_prim(i,6) * mcon_prim(i,4) - mass_offset_corrn)) &
 !                              * (scale_mass*(mcon_prim(j,6) * mcon_prim(j,4) - mass_offset_corrn))/ cons_mass_sigma**2 

           normeq_all(i,j) = normeq_all(i,j) + &
                msc_scl_factor*( mcon_prim(i,6) * mcon_prim(j,6)) &
                              *(mcon_prim(i,4)*scale_mass)*(mcon_prim(j,4)*scale_mass) / cons_mass_sigma**2  
!
        enddo
!        ! PT211030: need to include in the AtWb as well if there is an offset mass term to be considered 
!        !           (e.g. for simulations that don't conserve mass)
        AtWb(i) = AtWb(i) + msc_scl_factor*mass_offset * mcon_prim(i,6)*(mcon_prim(i,4)*scale_mass)/cons_mass_sigma**2

     enddo
  else
     call status_update('STATUS','ADDNORM','addnorm',' ',"NOT applying conservation of mass condition",0)
  endif


!****************************************************************
! RM210505: Moved this up so don't have to duplicate normeq_all
! solve the combined normal equations
! PT/SA190820: replace LS_normSolve with Chol_normsolve
! DEBUG
! copy nparam_all to tmp variable
  allocate(normeq_tmp(nparam_t,nparam_t))
  normeq_tmp = normeq_all
  call status_update('STATUS','ADDNORM','addnorm',outfile_stem,"Calling Chol_normSolve to invert the normal equations",0)
  call Chol_normSolve(nparam_t,normeq_all,AtWb,adjust)
  VCV = normeq_all
!****************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  O U T P U T    T H E    S I N G U L A R    V A L U E S
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RM210208: if output_SVD is true run SVD and output singular values
  if(output_SVD)then
    write(message,'(a)')'Calculating SVD of the AtWA matrix using DGESDD... '
    call status_update('STATUS','ADDNORM','addnorm',' ',message,0)

    ! allocate necessary arrays
    allocate(S(nparam_t))
    allocate(U(nparam_t,nparam_t))
    allocate(VT(nparam_t,nparam_t))
    allocate(iwork(8*nparam_t))
    allocate(dummy(1))

    ! mirror the normal equations
    do irow=1,nparam_t
      normeq_tmp(:,irow) = normeq_tmp(irow,:)
    enddo

    ! first call to dgesvd to calculate optimal size of work array 
    lwork = -1
    call DGESDD('S', nparam_t, nparam_t, normeq_tmp, nparam_t, S, U, nparam_t, VT, nparam_t, dummy, lwork, iwork, info)
    write(message,'(a,i10)')'Calculated optimal size of work array = ',int(dummy(1))
    call status_update('STATUS','ADDNORM','addnorm',' ',message,0)

    allocate(work(int(dummy(1))))

    ! Calculate SVD
    call DGESDD('S', nparam_t, nparam_t, normeq_tmp, nparam_t, S, U, nparam_t, VT, nparam_t, work, lwork, iwork, info)

    ! Calculate condition number
    CN = S(1)/S(nparam_t)    
    write(message,'(a,e15.6)')'SVD condition number = ',CN
    call status_update('STATUS','ADDNORM','addnorm',' ',message,0)

    ! Write out the Sigmas (Don't delete!)
    !write(message,'(a,a)')'Writing out ascii version of SVD sigmas ',trim(outfile_stem)//'.sigma'
    !call status_update('STATUS','ADDNORM','addnorm',' ',message,0)
    !LUSIG = 23
    !open(LUSIG,file=trim(outfile_stem)//'.sigma',status='unknown')
    !write(LUSIG,*)1.0,nmascon_t,CN
    !do i=1,nmascon_t
    !  write(LUSIG,*)(S(i))
    !enddo
    !close(LUSIG)

    ! RM210422: Write out the U, Vt, diag S and AtWb to a binary file (to be used as input to ~/gt/util/svd_invert)
    write(message,'(a,a)')'Writing out binary version of U and Vt matrices ',trim(outfile_stem)//'.SVD'
    call status_update('STATUS','ADDNORM','addnorm',' ',message,0)
    ! write out U matrix followed by Vt matrix
    rec_len = nparam_t*8
    open(LUSVD,file=trim(outfile_stem)//'.SVD',status='unknown',access='direct',form='unformatted',recl=rec_len(1))
    ! write header info
    irec = 1
    if(versions(1) == 2.0)then
      write(LUSVD,rec=irec)ndates,nparam_t,nmascon_t,nparam_t-nmascon_t
    else if(versions(1) == 1.0)then
      write(LUSVD,rec=irec)nfiles,nparam_t,nmascon_t,nparam_t-nmascon_t
    endif
    irec = 2
    write(LUSVD,rec=irec)msc_scl_factor,lambda,msc_param_const,trim(msc_constraint_file)
    if(versions(1) == 2.0)then
      do i=1,ndates
        irec=irec+1
        write(LUSVD,rec=irec)i,year_s(i),month_s(i),day_s(i),hr_s(i),min_s(i),sec_s(i),duration_s(i),normeq_files(1)
      enddo
    else if(versions(1) == 1.0)then
      do i=1,nfiles
        irec=irec+1
        write(LUSVD,rec=irec)i,year(i),month(i),day(i),hr(i),minute(i),sec(i),duration(i),normeq_files(i)
      enddo
    endif
    !************************* WRITE U matrix ***********************
    do i= 1, nparam_t
      irec=irec+1
      write(LUSVD,rec=irec)(U(i,j),j=1,nparam_t)
    enddo
    !****************************************************************
    !************************* WRITE Vt matrix **********************
    do i= 1, nparam_t
      irec=irec+1
      write(LUSVD,rec=irec)(VT(i,j),j=1,nparam_t)
    enddo
    !****************************************************************
    !*********************** WRITE singular vals ********************
    irec=irec+1
    write(LUSVD,rec=irec)(S(i),i=1,nparam_t)
    !****************************************************************
    !*************************** WRITE AtWb *************************
    irec=irec+1
    write(LUSVD,rec=irec)(AtWb(i),i=1,nparam_t)
    !****************************************************************
    close(LUSVD)
  endif

! ****************************************************
!
!         O U T P U T     S O L U T I O N
!
! ****************************************************
! write out the solution
! DEBUG: write to file addnorm.out
  open(30,file=trim(outfile_stem)//'.fit',status='unknown')      ! in a sort of .fit file version
  write(30, '(a)') 'V2 ADDNORM   FIT solution:  '
  LUADD = 31
  open(LUADD,file=trim(outfile_stem)//'.vcv',status='unknown')   ! in a sort of .vcv file version
! PT190301: added an extra blank so that "VCV" is in columns 15-17 to match GRACEFIT/GRACESIM
  write(LUADD, '(a)') 'V2 ADDNORM    VCV solution:  '

  if(versions(1) == 2.0)then
    output_nfiles = ndates
  else
    output_nfiles = nfiles
  endif
  write(message,'(a,f15.4,a,i6)')'Write out solution with mascon constraints of (m) : ',msc_param_const &
       ,' Number of files: ',output_nfiles
  call status_update('STATUS','ADDNORM','addnorm',' ',message,0)
  write(30,'(a,f15.4,a,i6,a,i6)')'Solution with mascon constraints of (m) : ',msc_param_const,' Number of files: ' &
                            ,output_nfiles,"  Number of mascons: ",nmascon_t
  write(LUADD,'(a,f15.4,a,i6,a,i6)')'Solution with mascon constraints of (m) : ',msc_param_const,' Number of files: ' &
                            ,output_nfiles,"  Number of mascons: ",nmascon_t
! record the name of the regularization file used
! PT190301: fixed bug here (was missing the last character variable in the format statement)
  write(30,'(a,f6.3,a,a)')    "Regularisation scale factor and constraint file: ",msc_scl_factor,"   ",msc_constraint_file
  write(LUADD,'(a,f6.3,a,a)') "Regularisation scale factor and constraint file: ",msc_scl_factor,"   ",msc_constraint_file

! write the input file names, their orbit integration start times and durations
! PT210217: if a snorm was read then "ndates" contains the number of files used in the stacking. Transfer to "nfiles"
  if(versions(1) == 2.0)then
    do i=1,ndates
       write(30,'(a,i3,1x,a,6i5,a,f10.5, a)')"File: ",i &
          ," Start: ",year_s(i),month_s(i),day_s(i),hr_s(i),min_s(i),sec_s(i) &
          ," Duration: ",duration_s(i), trim(normeq_files_s(i)) 
       write(LUADD,'(a,i3,1x,a,6i5,a,f10.5,a )')"File: ",i  & 
          ," Start: ",year_s(i),month_s(i),day_s(i),hr_s(i),min_s(i),sec_s(i) &
          ," Duration: ",duration_s(i), trim(normeq_files_s(i))
    enddo
  else if (versions(1) == 1.0)then
    do i=1,nfiles
      ! PT170601: re-ordered to put the filename (of variable length) at the end of the line
       write(30,'(a,i3,1x,a,6i5,a,f10.5,1x,a)')"File: ",i &
          ," Start: ",year(i),month(i),day(i),hr(i),minute(i),sec(i) &
          ," Duration: ",duration(i) &
          ,trim(normeq_files(i)) 
       write(LUADD,'(a,i3,1x,a,6i5,a,f10.5,1x,a)')"File: ",i  & 
          ," Start: ",year(i),month(i),day(i),hr(i),minute(i),sec(i) &
          ," Duration: ",duration(i) &
          ,trim(normeq_files(i)) 
    enddo
  endif


  write(LUADD,'(/,a)') ' SOLUTION A PRIORI AND VECTOR:'
  write(LUADD,'(a)') ' PARAMETER                     A PRIORI             VECTOR            SIGMA'   
  write(30,'(a)') ' PARAMETER                     A PRIORI       '   

!****************************************************************

  if ( IC_solve )then
     iparm = nmascon_t + nmsc_tidal_t*2
     ! write out the IC adjustments (for each file) first
     do j=1,nfiles
        do k=1,nICs(j)
           iparm = iparm + 1
           write(30 ,'(a30,f17.5,f10.2,f21.7,f12.5)') prmnam_all(iparm),apr_prm(iparm)*1.d3,adjust(iparm)*1.d3 &
                ,apr_prm(iparm)*1.d3+adjust(iparm)*1.d3,dsqrt(VCV(iparm,iparm))*1.d3
           write(LUADD,'(a30,f17.7,f17.7,f17.7)') prmnam_all(iparm),apr_prm(iparm),apr_prm(iparm)+adjust(iparm),dsqrt(VCV(iparm,iparm))
        enddo
        write(30,*)" "
     enddo
  endif

! now write out the mascon estimates
  if(msc_solve) then
     do j=1,nmascon_t
        write(30 ,'(a29,f18.5,f12.5,f19.5,f12.5)') prmnam_all(j)(2:30),apr_prm(j),adjust(j) &
             ,apr_prm(j)+adjust(j),dsqrt(VCV(j,j))
        write(LUADD,'(a29,f18.7,f17.7,f17.7)') prmnam_all(j)(2:30),apr_prm(j),apr_prm(j)+adjust(j),dsqrt(VCV(j,j))
     enddo
  endif

! now tidal mascon amplitudes
  if(tmsc_solve)then
     iparm = nmascon_t
     do j=1,nmsc_tidal_t*2
        iparm = iparm + 1
        write(30,'(a30,f17.5,f12.5,f19.5,f12.5)') prmnam_all(iparm),apr_prm(iparm),adjust(iparm) &
             ,apr_prm(iparm)+adjust(iparm),dsqrt(VCV(iparm,iparm))
        write(LUADD,'(a30,f17.7,f17.7,f17.7)') prmnam_all(iparm),apr_prm(iparm),apr_prm(iparm)+adjust(iparm),dsqrt(VCV(iparm,iparm))
     enddo
  endif

! PT200203: output the inverse square root of the diagonal constraint matrix for the mascons
  if(msc_solve)then
    write(30,'(a)')"MASCON CONSTRAINT SIGMAS"
    do iparm=1,nmascon_t
      write(30,'(i7,f15.6,a)')iparm,dsqrt(1.d0/msc_const(iparm,iparm) ),"  mascon constraint (m)"
    enddo
  endif


  close(30)

! Write the vcv solution   
  write(LUADD,*)' VCV SOLUTION '
! PT190821: to save space, don't write out the actual VCV matrix ....
  call status_update('STATUS','ADDNORM','addnorm',' ',"As of 2019-08-21, VCV matrix is NOT written out to the .vcv file",0)
! DO NOT DELETE THESE THREE LINES!
!  call status_update('STATUS','ADDNORM','addnorm',' ',"VCV matrix IS CURRENTLY written out to the .vcv file",0)
! PT191029: loop over rows, then columns (not the other way, as it was). Change variable names to irow,icol
!  do irow = 1, nparam_t
!     write(LUADD,*) ((VCV(irow, icol)), icol = 1, nparam_t)   
!  enddo
  close(LUADD)

! PT220608: compute and output the correlation matrix, if requested
  if(output_correl)then
    call status_update('STATUS','ADDNORM','addnorm',' ',"Compute the correlation matrix",0)

    allocate(corr(nparam_t,nparam_t))
    do irow = 1, nparam_t
       do icol = 1, nparam_t
         corr(irow,icol) = VCV(irow,icol)/( dsqrt(VCV(irow,irow))*dsqrt(VCV(icol,icol)) ) 
       enddo
    enddo


    ! write out the matrix
    call status_update('STATUS','ADDNORM','addnorm',trim(outfile_stem)//'.correl' &
             ,"Write out the correlation matrix",0)
    open(352,file=trim(outfile_stem)//'.correl',status='unknown')   ! 1 header line, then param number+full correlation matrix
    write(352,*)'Correlation matrix'
    do i = 1, nparam_t
       write(352,*)i,(corr(i,j),j = 1, nparam_t)
    enddo
  endif
  
  call status_update('STATUS','ADDNORM','addnorm',outfile_stem,'End of ADDNORM.',0)

end program



