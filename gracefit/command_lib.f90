!********************************************************************************************************************************
!  File: command_lib.f90
!
!  Purpose: Set of subroutines used to open and read in data
!
!  Author: Thomas Greenspan
!          (some subroutines taken or inspired from subroutines by Simon McClusky or Paul Tregoning)
!
!  API:
!       command_storeCommand     : Store an indexed value of specified command
!       command_readSetup        : Read in setup as directed by command file
!       command_readDataWeights  : Read in data weights on observationsfrom command file
!       command_readAprConst     : Read constraints on parameters from command file
!       command_printSetup       : Prints out setup if requested
!       command_close            : Closes command file
!
!  August 1, 2013
!
!********************************************************************************************************************************
!********************************************************************************************************************************
! command_storeCommand: Read and store an indexed command value (ex: 3rd value associated with the given
!                       command), given the command file unit number, the command to be searched for, the
!                       variable in which to store value, the index of the command value to be stored, the
!                       total number of values that should be associated with given command and a message
!                       concerning the type of data being read in.
!                       NOTE: If total number of values indicated is 0, no message is printed
!                             There must also be at least the number of values associated with the command
!                             as indicated or, in case the indicated number is 0, at least 1 value.
! Author: Thomas Greenspan
!********************************************************************************************************************************

subroutine command_storeCommand(file_unit_num,command_name,value,ivalue,end_num,data_message,print_message)

  implicit none

  !********************  Variable declarations ********************

  integer*4,    intent(in) :: file_unit_num     ! File unit number of command file to be read from
  integer*4,    intent(in) :: ivalue            ! Index number of value to be read in from command
  integer*4,    intent(in) :: end_num           ! Total number of values that should be associated with given command
  integer*4,    intent(in) :: print_message     ! Indicator whether to print out message or not 1 = yes, 0 = no
  character(*), intent(in) :: command_name      ! Name of command to be searched for
  character(*), intent(in) :: data_message      ! Message concerning type of data being read in
  double precision, intent(out) :: value        ! Vairable in which command value is to be stored

  integer*4      :: lcmd       ! Variable that records whether command was found and if there was a value associated with it
  !                                  -1 if command not found, 0 if command found but no value, > 0 if command found with value
  character(256) :: wcmd       ! Variable in which value of command stored
  character(128) :: message                  ! Message printed out
  double precision :: dummy_int(end_num-1)   ! Dummy variable used to store unwanted data from wcmd
  !****************************************************************

  ! Assure that the values given makes sense (one or greater for index and total values, 0 or 1 for print message indicator)
  if(ivalue < 1)  call status_update('FATAL','GRACEFIT','command_storeCommand',' ','Index of value must be one or greater',0)
  if(end_num < 1) call status_update('FATAL','GRACEFIT','command_storeCommand',command_name, &
       'Total number of values must be one or greater',0)
  if(print_message /= 0 .and. print_message /= 1) &
       call status_update('FATAL','GRACEFIT','command_storeCommand',' ','print indicator must be 0 or 1',0)

  !************************* STORE COMMAND ************************
  lcmd = 0
  call get_keyword(file_unit_num,command_name,wcmd,lcmd,1)
  ! Report failure to find command or value associated with it
  if(lcmd <= 0 .and. print_message == 1)then
     write(message,'(a,a,a)') 'Command file not open or command "',command_name,'" not found or no value associated with it'
     call status_update('FATAL','GRACEFIT','command_storeCommand',' ',message,lcmd)
  endif
  ! If command exists and has values, store wanted value
  if(lcmd > 0)then
     read(wcmd,*) dummy_int(1:ivalue-1),value,dummy_int(ivalue+1:end_num-1)  ! Last part just there as a safety check
  endif
  !****************************************************************

  if(print_message == 1 .and. ivalue == end_num)then
     write(message,'(a,a,e25.10,a)') data_message,' found and stored (',value,')'
     call status_update('STATUS','GRACEFIT','command_storeCommand',' ',message,0)
  endif
  !****************************************************************

  return
end subroutine command_storeCommand

!********************************************************************************************************************************
!********************************************************************************************************************************
! command_readSetup: Read in what observations and parameters are to be used and set
!                    common variables (in gracefit.h90) accordingly
!                    NOTE: Actual number for mascons is read in from GTORB header. nmascons_t is
!                          temporarily used as an indicator. nparam_t and prmnam are also updated
!                          when real number of mascons is read in. (done in header_readGTORB)
! Author: Thomas Greenspan
!********************************************************************************************************************************

subroutine command_readSetup()
  use gracefit_mod
  use soln_mod       ! PT190905: added to get declaration of prmnam
  implicit none


  include 'input.h90'

  !********************  Variable declarations ********************

  integer*4 :: isatt
  integer*4 :: i
  double precision :: satellites(2)  ! Indicators on which satellites to consider
  double precision :: gobs(3)        ! Indicators on which gps observations to consider
  double precision :: antenna_obs(3) ! Indicators on which antenna conditions to use
  double precision :: shadow         ! Indicator  on whether to use shadow conditions
  double precision :: kband(3)       ! Indicators on which kband observations to use
  double precision :: mean_accel     ! Indicator  on whether to use mean acceleration conditions
  double precision :: gprm(2)        ! Indicators on which gps parameters to use
  double precision :: scale_ind      ! Indicators on whether to use parameters for scale
  double precision :: bias_ind       ! Indicators on whether to use parameters for bias
  double precision :: twopr_ind      ! Indicator  on whether to use parameters for twice-per-rev along-track accelerations
  double precision :: onepr_ind      ! Indicator  on whether to use parameters for  once-per-rev along-track accelerations
  double precision :: antenna_offset ! Indicator on whether to use parameters for antenna offsets
  double precision :: mascons        ! Number of mascons to be used for parameters
  double precision :: tidal_mascons  ! Number of mascons to be used for parameters
  double precision :: gnv_int        ! Interval at which GNV1B data is read in

  double precision :: tmp_start_neq  ! epoch at which to start stacking normal equations
  double precision :: tmp_end_neq    ! epoch at which to stop  stacking normal equations
  double precision :: est_rpy        ! indicator on whether to estimate rpy orientation angle corrections
  double precision :: tmp_gracefit_step  ! step at which to stop running gracefit
  double precision :: tmp_mask       ! use/not use outlier.mask file


  character(1), parameter, dimension(3) :: UCOORD = ['X','Y','Z']
  character(1), parameter, dimension(3) :: LCOORD = ['x','y','z']
  character(1), parameter, dimension(2) :: RevCRD = ['S','C'    ]
  character*1000 message

! temporary real value
  real(kind=8) :: tmp_R8

  !****************************************************************

  !********************* SATELLITES TO BE USED ********************

! PT180625: read the mission (GRACE: 0, GRACE FO: 1, GRACE II: 2)
  call command_storeCommand(LUCMD,'mission',tmp_R8,1,1,' ',0)
  mission = nint(tmp_R8)

  do isatt = 1, 2
     call command_storeCommand(LUCMD,'satellites',satellites(isatt),isatt,2,' ',0)
     if(int(satellites(isatt)) /= 0 .and. int(satellites(isatt)) /= 1) &
          call status_update('FATAL','GRACEFIT','command_readSetup',' ','satellite indicators must be 0 or 1',0)
  enddo
  SAT_1 = int(satellites(1))
  SAT_2 = int(satellites(2))
  ! Total number of satellites
  nsat_t = SAT_1 + SAT_2
  if(nsat_t == 0) call status_update('FATAL','GRACEFIT','command_readSetup',' ',&
       'No satellites to be considered? Then we are done!',0)
  ! Names of satellites
! PT180625: this is now mission dependent
  satnam = '~'
  if(mission == 0)then
    if(SAT_1 == 1) satnam(1) = 'A'
    if(SAT_2 == 1) satnam(nsat_t) = 'B'
  else if (mission == 1)then
    if(SAT_1 == 1) satnam(1) = 'C'
    if(SAT_2 == 1) satnam(nsat_t) = 'D'
  else if (mission == 2)then
    if(SAT_1 == 1) satnam(1) = 'E'
    if(SAT_2 == 1) satnam(nsat_t) = 'F'
  endif    
  !****************************************************************

  !**************************** EPOCHS ****************************

  ! Read in interval at which GNV1B data should be read in
  call command_storeCommand(LUCMD,'gnv_int',gnv_int,1,1,' ',0)
  if(int(gnv_int) < 0) call status_update('FATAL','GRACEFIT','command_readSetup',' ','GNV interval must be positive',0)
  if(mod(int(gnv_int),int(epoch_interval)) /= 0) call status_update('FATAL','GRACEFIT','command_readSetup',' ',&
       'GNV interval must be a multiple of the epoch interval',0)
  gnv_interval = int(gnv_int/epoch_interval)

  ! read epochs for which normal wquations should be stacked
  call command_storeCommand(LUCMD,'start_neq',tmp_start_neq,1,1,' ',0)
  if(tmp_start_neq <= 0)then
     start_neq = 1
  else
     start_neq = int(tmp_start_neq)
  endif
  !PT220125: initialize this to zero
  tmp_end_neq = 0
  call command_storeCommand(LUCMD,'end_neq',tmp_end_neq,1,1,' ',0) ! set to nepochs_t in gracefit if not given a value here
  end_neq = int(tmp_end_neq)

! PT191126: read to what step in the process gracefit should be run
! 0: all of it; 1: just to calculating prefit residuals; 2: just to stacking normal equations
  gracefit_step = 0  ! set a default value
  call command_storeCommand(LUCMD,'gracefit_step',tmp_gracefit_step,1,1,' ',0)
  gracefit_step = int(tmp_gracefit_step)



  !****************************************************************

  !************************** PARAMETERS **************************

  ! GPS parameters
  do i = 1, 2
     call command_storeCommand(LUCMD,'gprm',gprm(i),i,2,' ',0)
     if(int(gprm(i)) /= 0 .and. int(gprm(i)) /= 1) &
          call status_update('FATAL','GRACEFIT','command_readSetup',' ','gprm indicators must be 0 or 1',0)
  enddo
  if(int(gprm(2)) == 1 .and. int(gprm(1)) == 0) &
       call status_update('FATAL','GRACEFIT','command_readSetup',' ', &
       'Cannot have velocity parameters without position parameters',0)
  ngprm_t = sum(int(gprm(1:2))*3)
  ! Scale parameters
  call command_storeCommand(LUCMD,'scale',scale_ind,1,1,' ',0)
  if(int(scale_ind) /= 0 .and. int(scale_ind) /= 1) &
       call status_update('FATAL','GRACEFIT','command_readSetup',' ','Scale indicator must be 0 or 1',0)
  nScale_t = int(scale_ind)*3
  ! Bias parameters
  call command_storeCommand(LUCMD,'bias',bias_ind,1,1,' ',0)
  if(int(bias_ind) /= 0 .and. int(bias_ind) /= 1) &
       call status_update('FATAL','GRACEFIT','command_readSetup',' ','Bias indicator must be 0 or 1',0)
  nBias_t = int(bias_ind)*3
  ! Once-per-rev parameters
  call command_storeCommand(LUCMD,'once-per-rev',onepr_ind,1,1,' ',0)
  if( int(onepr_ind) /= 0 .and. int(onepr_ind)  /= 1 ) &
       call status_update('FATAL','GRACEFIT','command_readSetup',' ','Once-per-rev indicators must be 0 or 1',0)
  if(int(onepr_ind) == 1)n_onepr_t = 2
  ! Twice-per-rev parameters
  call command_storeCommand(LUCMD,'twice-per-rev',twopr_ind,1,1,' ',0)
  if( int(twopr_ind) /= 0 .and. int(twopr_ind)  /= 1 ) &
       call status_update('FATAL','GRACEFIT','command_readSetup',' ','Twice-per-rev indicators must be 0 or 1',0)
  if(int(twopr_ind) == 1)ntwopr_t = 2
  ! Antenna offset parameters
  call command_storeCommand(LUCMD,'antenna_offset',antenna_offset,1,1,' ',0)
  if(int(antenna_offset) /= 0 .and. int(antenna_offset) /= 1) &
       call status_update('FATAL','GRACEFIT','command_readSetup',' ','antenna offset parameter indicator must be 0 or 1',0)
  nantoff_t = int(antenna_offset)*3
  ! Mascon parameters
  call command_storeCommand(LUCMD,'mascons',mascons,1,1,' ',0)
  if(int(mascons) /= 0 .and. int(mascons) /= 1) &
       call status_update('FATAL','GRACEFIT','command_readSetup',' ','mascons indicator must be 0 or 1',0)
  nmascons_t = int(mascons) ! Set as indicator for now until actual number is read in from GTORB header
  ! Mascon tidal amplitude parameters
  call command_storeCommand(LUCMD,'mascon_tidal',tidal_mascons,1,1,' ',0)
  if(int(tidal_mascons) /= 0 .and. int(tidal_mascons) /= 1) &
       call status_update('FATAL','GRACEFIT','command_readSetup',' ','tidal mascons indicator must be 0 or 1',0)
  est_msc_tides = int(tidal_mascons)
  ! PT170405: roll/pitch/yaw orientation offset parameters
  call command_storeCommand(LUCMD,'est_rpy',est_rpy,1,1,' ',0)
  if(int(est_rpy) /= 0 .and. int(est_rpy) /= 1) &
       call status_update('FATAL','GRACEFIT','command_readSetup',' ','rpy orientation indicator must be 0 or 1',0)
  nrpy_t = int(est_rpy)*3

  ! Total number of orbital parameters
  norbprm_t = ngprm_t + nScale_t + nBias_t + n_onepr_t + ntwopr_t + nrpy_t
  ! Total number of satellite-specific parameters
  nsatprm_t = norbprm_t + nantoff_t
  ! Total number of parameters (NOTE: Mascons and tidal mascons not included yet. Included in header_readGTORB)
  nparam_t = nsatprm_t*nsat_t
  !****************************************************************

  !************************* OBSERVATIONS *************************

  ! GPS-observations
  do i = 1, 2
     call command_storeCommand(LUCMD,'gobs',gobs(i),i,2,' ',0)
     if(int(gobs(i)) /= 0 .and. int(gobs(i)) /= 1) &
          call status_update('FATAL','GRACEFIT','command_readSetup',' ','gobs indicators must be 0 or 1',0)
  enddo
  if(int(gobs(2)) == 1 .and. int(gobs(1)) == 0) &
       call status_update('FATAL','GRACEFIT','command_readSetup',' ', &
       'Cannot have velocity observations without position observations',0)
  ngobs_t = sum(int(gobs(1:2))*3)

  ! Antenna offset conditions
  do i = 1, 2
     call command_storeCommand(LUCMD,'antenna_obs',antenna_obs(i),i,2,' ',0)
     if(int(antenna_obs(i)) /= 0 .and. int(antenna_obs(i)) /= 1) &
          call status_update('FATAL','GRACEFIT','command_readSetup',' ','antenna_obs indicators must be 0 or 1',0)
  enddo
  if(int(antenna_obs(2)) == 1 .and. int(antenna_obs(1)) == 0) &
       call status_update('FATAL','GRACEFIT','command_readSetup',' ', &
       'Cannot have velocity antenna conditions without position antenna conditions',0)
  if(int(antenna_obs(2)) == 1 .and. ngobs_t < 6)then
     call status_update('WARNING','GRACEFIT','command_readSetup',' ', &
        'Cannot have velocity antenna conditions without velocity observations. Velocity antenna offsets will not be considered',0)
     antenna_obs(2) = 0
  endif
  nCoM_t = sum(int(antenna_obs(1:2))*3)

  ! Shadow condition
  call command_storeCommand(LUCMD,'shadow',shadow,1,1,' ',0)
  if(int(shadow) /= 0 .and. int(shadow) /= 1) &
       call status_update('FATAL','GRACEFIT','command_readSetup',' ','shadow condition indicator must be 0 or 1',0)
  ! PT140919: fix bug in setting ncond_t if bias or scale are not estimated. There are two condition equations per satellite for
  !           each of scale and bias. That is, 8 condition equations if both scale and bias are being estimated.
  if(int(shadow) == 1 )then
     ! PT150609: increased from 2 to 4, since we now add another condition observation 1/2 revolution away from the shadow obs
     ! PT160420: with the new acc_Z constraint (the mean of lead/trail satellite obs during eclipse = -5 nm/s^2, which accounts
     !           for the atmospheric drag component that occurs in the SRF Z direction due to tilting of the satellites onto the
     !           LOS vector), we now need only 1 condition equation for acc_Z for both satellites. However, we also add a new
     !           condition equation that makes the
     ncond_t = 4
  endif
  if(int(shadow) == 1 .and. ncond_t == 0)then
     call status_update('WARNING','GRACEFIT','command_readSetup',' ',&
          'Shadow condition useless if neither scale or bias are being estimated. Condition will not be applied',0)
     shadow = 0
  endif

  ! Kband observations
  do i = iRANGE, iRACC
     call command_storeCommand(LUCMD,'kband',kband(i),i,3,' ',0)
     if(int(kband(i)) /= 0 .and. int(kband(i)) /= 1) &
          call status_update('FATAL','GRACEFIT','command_readSetup',' ','kband indicators must be 0 or 1',0)
     if(nsat_t < 2 .and. int(kband(i)) == 1)then
        call status_update('WARNING','GRACEFIT','command_readSetup',' ',&
             'kband observations require both satellites to be considered. kband will not be considered',0)
        kband(i) = 0
     endif
  enddo
  nkbr_t = int(kband(iRANGE))
  nkbrr_t = int(kband(iRRATE))
  nkbra_t = int(kband(iRACC))
! PT181003: this statement is now wrong. We now use range acceleration!
!  if (nkbra_t > 0) call status_update('WARNING','GRACEFIT','command_readSetup',' ', &
!       'Program currently does not support looking at Range Acceleration',0)
  !   Total kband observations
  nkobs_t = nkbr_t + nkbrr_t + nkbra_t

  ! Conditions on mean acceleration
  call command_storeCommand(LUCMD,'mean_accel',mean_accel,1,1,' ',0)
  if(int(mean_accel) /= 0 .and. int(mean_accel) /= 1) &
       call status_update('FATAL','GRACEFIT','command_readSetup',' ','mean acceleration condition indicator must be 0 or 1',0)
  if(nsat_t < 2 .and. int(mean_accel) == 1)then
     call status_update('WARNING','GRACEFIT','command_readSetup',' ',&
          'Mean acceleration condition requires both satellites to be considered. Condition will not be applied',0)
     mean_accel = 0
  endif
  if(int(mean_accel) == 1 .and. nScale_t + nBias_t == 0)then
     call status_update('WARNING','GRACEFIT','command_readSetup',' ',&
          'Mean acceleration condition useless if neither scale or bias are being estimated. Condition will not be applied',0)
     mean_accel = 0
  endif
  naobs_t = int(mean_accel)*3

  ! PT150527: conditions on relation between accelerometer scales
  call command_storeCommand(LUCMD,'use_accel_scl_const',use_accel_scale_const,1,1,' ',0)

  ! PT180726: name of model(s) to use to linearise the accelerometer obs
  acc_obs_model = 0.d0
  do i = 1, 3
    call command_storeCommand(LUCMD,'accel_obs_model',acc_obs_model(i),i,3,'accelerometer linearisation model',0)
  enddo

  ! Total number of sat specific observations
  nsatobs_t = ngobs_t + nCoM_t + ncond_t*nsat_t  !@# Not used right now
  ! Total number of observations used
  nobs_t = ngobs_t*nsat_t + nCoM_t*nsat_t + ncond_t*nsat_t + nkobs_t + naobs_t
  write(message,'(a,i7)')'Total observations per epoch (including conditional eqtns): ',nobs_t
  call status_update('STATUS','GRACEFIT','command_readSetup',' ',message,0)
  !****************************************************************

  !**************************** INDICES ***************************

  ! Indices for observations
  !@# TEMPORARY ordering...
  igobs  = 1
  ikbr   = ngobs_t*nsat_t + 1
  ikbrr  = ngobs_t*nsat_t + nkbr_t + 1
  ikbra  = ngobs_t*nsat_t + nkbr_t + nkbrr_t + 1
  iCoM   = ngobs_t*nsat_t + nkbr_t + nkbrr_t + nkbra_t + 1
  icond  = ngobs_t*nsat_t + nkbr_t + nkbrr_t + nkbra_t + nCoM_t*nsat_t + 1
  iacobs = ngobs_t*nsat_t + nkbr_t + nkbrr_t + nkbra_t + nCoM_t*nsat_t + ncond_t*nsat_t + 1
  ! if (nkbra_t == 0) ikbra = -9999
  ! if (nkbrr_t == 0 ) ikbrr =-9999
  ! if (nkbr_t == 0) ikbr = -9999

  ! Indices for parameters
  iICparm  = 1
  iScale   = ngprm_t + 1
  iBias    = ngprm_t + nScale_t + 1
  ionepr   = ngprm_t + nScale_t + nBias_t + 1
  itwopr   = ngprm_t + nScale_t + nBias_t + n_onepr_t + 1
  iantoff  = ngprm_t + nScale_t + nBias_t + n_onepr_t + ntwopr_t + 1
  irpy     = ngprm_t + nScale_t + nBias_t + n_onepr_t + ntwopr_t + nantoff_t + 1
  imascons = nsatprm_t*nsat_t + 1
  !****************************************************************

  !*********************** OBSERVATION NAMES **********************

  obsnam = ' '
  !   Satellite specific observations/conditions  !@# CHANGE ORDERING THIS IS RIDICULOUS
  do isatt = 1, nsat_t
     do i = igobs, igobs+ngobs_t-1
        if(ngobs_t > 0 .and. i <= 3)then        ! position
           write(obsnam(i+(isatt-1)*ngobs_t),'(i3,a,a2,a,a1,a)') i+(isatt-1)*ngobs_t,'. SAT',satnam(isatt),&
                ':  ',UCOORD(i-igobs+1),'p     (m)'
        else if(ngobs_t > 3 .and. i <= 6)then   ! velocity
           write(obsnam(i+(isatt-1)*ngobs_t),'(i3,a,a2,a,a1,a)') i+(isatt-1)*ngobs_t,'. SAT',satnam(isatt),&
                ':  ',UCOORD(i-igobs+1-3),'v     (m/s)'
        endif
     enddo  ! End of observation loop
  enddo  ! End of satellite loop

  if(nkbr_t /= 0) write(obsnam(ikbr),'(1x,i2,a)')   ikbr, '.   R            (m)'      ! Range
  if(nkbrr_t /= 0) write(obsnam(ikbrr),'(1x,i2,a)') ikbrr,'.   RR           (m/s)'    ! Range Rate
  if(nkbra_t /= 0) write(obsnam(ikbra),'(1x,i2,a)') ikbra,'.   RA           (m/s^2)'  ! Range Acceleration

  do isatt = 1, nsat_t
     do i = iCoM, iCoM+nCoM_t-1
        if(nCoM_t > 0 .and. i <= iCoM-1+3)then        ! position
           write(obsnam(i+(isatt-1)*nCoM_t),'(i3,a,a2,a,a1,a)') i+(isatt-1)*nCoM_t,'. SAT',satnam(isatt),&
                ': ant',UCOORD(i-iCoM+1),  'p   (m)'
        else if(nCoM_t > 3 .and. i <= iCoM-1+6)then   ! velocity
           write(obsnam(i+(isatt-1)*nCoM_t),'(i3,a,a2,a,a1,a)') i+(isatt-1)*nCoM_t,'. SAT',satnam(isatt),&
                ': ant',UCOORD(i-iCoM+1-3),'v   (m/s)'
        else if(nCoM_t > 6 .and. i <= iCoM-1+9)then   ! acceleration
           write(obsnam(i+(isatt-1)*nCoM_t),'(i3,a,a2,a,a1,a)') i+(isatt-1)*nCoM_t,'. SAT',satnam(isatt),&
                ': ant',UCOORD(i-iCoM+1-6),'a   (m/s^2)'
        endif
     enddo  ! End of observation loop
  enddo  ! End of satellite loop

  do isatt = 1, nsat_t
     do i = icond, icond+ncond_t-3  ! PT150609: reduce this from -1 to -3 because we now have an additional 2 conditional obs equations
        write(obsnam(i+(isatt-1)*ncond_t),'(i3,a,a2,a,a1,a)') i+(isatt-1)*ncond_t,'. SAT',satnam(isatt),&
             ': cond',UCOORD(i+1-icond+1),'a  (um/s^2)'
     enddo  ! End of observation loop
  enddo  ! End of satellite loop

  do i = iacobs, iacobs +naobs_t-1
     write(obsnam(i),'(i3,a,a1,a)') i,'.  Mean',UCOORD(i-iacobs+1),'a        (um/s^2)'
  enddo
  !************************ PARAMETER NAMES ***********************

  prmnam = ' '
  !   Satellite specific parameters
  do isatt = 1, nsat_t
     do i = 1, nsatprm_t
        if(ngprm_t > 0 .and. i <= 3)then            ! position
           write(prmnam(i+(isatt-1)*nsatprm_t),'(i3,a,a2,a,a1,a)') i+(isatt-1)*nsatprm_t,'. SAT',satnam(isatt), &
                ': ',UCOORD(i),'0   (m)'
        else if(ngprm_t > 3 .and. i < iScale)then   ! velocity
           write(prmnam(i+(isatt-1)*nsatprm_t),'(i3,a,a2,a,a1,a)') i+(isatt-1)*nsatprm_t,'. SAT',satnam(isatt), &
                ': ',UCOORD(i-3),'V0  (m/s)'
        else if(i < iBias)then                      ! Scale
           write(prmnam(i+(isatt-1)*nsatprm_t),'(i3,a,a2,a,a1,a)') i+(isatt-1)*nsatprm_t,'. SAT',satnam(isatt), &
                ': scl',LCOORD(i-iScale+1),' (n/a)'
        else if(i < ionepr)then                     ! Bias
           write(prmnam(i+(isatt-1)*nsatprm_t),'(i3,a,a2,a,a1,a)') i+(isatt-1)*nsatprm_t,'. SAT',satnam(isatt), &
                ': bs',LCOORD(i-iBias+1),'  (um/s^2)'
        else if(i < itwopr)then                    ! once-per-rev
           write(prmnam(i+(isatt-1)*nsatprm_t),'(i3,a,a2,a,a1,a)') i+(isatt-1)*nsatprm_t,'. SAT',satnam(isatt), &
                ': 1pr',RevCRD(i-ionepr+1),'  (um/s^2)'
        else if(i < iantoff)then                    ! two-per-rev
           write(prmnam(i+(isatt-1)*nsatprm_t),'(i3,a,a2,a,a1,a)') i+(isatt-1)*nsatprm_t,'. SAT',satnam(isatt), &
                ': 2pr',RevCRD(i-itwopr+1),'  (um/s^2)'
        else if (i < irpy) then                     ! Antenna offsets
           write(prmnam(i+(isatt-1)*nsatprm_t),'(i3,a,a2,a,a1,a)') i+(isatt-1)*nsatprm_t,'. SAT',satnam(isatt), &
                ': Gan',UCOORD(i-iantoff+1),' (m)'
        else                                        ! rpy offsets
           write(prmnam(i+(isatt-1)*nsatprm_t),'(i3,a,a2,a,a1,a)') i+(isatt-1)*nsatprm_t,'. SAT',satnam(isatt), &
                ': rpy',UCOORD(i-iantoff+1),' (mrad)'
        endif
     enddo  ! End of nsatprm loop
  enddo  ! End of satellite loop

  ! NOTE: Mascons not taken care of yet. Dealt with in header_readGTORB, called in input_readGTORB

  !****************************************************************

  !*********************** KBRR OMC DECOMPOSITION *****************
  call command_storeCommand(LUCMD,'kbrr_decomp_components',use_kbrr_decomp(1),1,2,' ',1)
  ! PT170911: the EMD component at which to start summing them to reconstruct the signals
  call command_storeCommand(LUCMD,'kbrr_decomp_components',use_kbrr_decomp(2),2,2,' ',1)

  call command_storeCommand(LUCMD,'kb_NR', use_kb_nr(1),1,2," ",1)
  call command_storeCommand(LUCMD,'kb_NR', use_kb_nr(2),2,2," ",1)

! PT/RmCG190722: add the fft options for kb data
  call command_storeCommand(LUCMD,'kb_cos_fft',use_kb_cos_fft(1),1,3,' ',1)
  call command_storeCommand(LUCMD,'kb_cos_fft',use_kb_cos_fft(2),2,3,' ',1)
  call command_storeCommand(LUCMD,'kb_cos_fft',use_kb_cos_fft(3),3,3,' ',1)

! PT/RMcG190722: over-ride the decomp setting with the cosine fft if both are set 
  if(use_kb_cos_fft(1) > 0)use_kbrr_decomp(1) = 0.d0

  !**************** MASCON length-dependent constraint ************
  call command_storeCommand(LUCMD,'msc_len_const',use_msc_len_const,1,1,' ',1)

  !**PT150721 ***** MASCON mass conservation constraint ************
  call command_storeCommand(LUCMD,'msc_conserve_mass',use_msc_conserve_mass,1,1,' ',1)


  !****************************************************************
  ! PT211206: do we want to use the outlier.mask file or not?
  call command_storeCommand(LUCMD,'use_outlier_mask',tmp_mask,1,1,'Do we use outlier.mask?',0)
  if(int(tmp_mask) == 1)then
    use_outlier_mask = .true.
  else
    use_outlier_mask = .false.
  endif

  !****************************************************************

  return
end subroutine command_readSetup

!********************************************************************************************************************************
!********************************************************************************************************************************
! command_readDataWeights: Reads in uncertainties from command file and places values into
!                          data weight array accordingly
! Author: Thomas Greenspan
!         (taken in large part from sections of read_input.f90 by R. King itself based on subroutine of same name by S. McClusky)
!********************************************************************************************************************************

subroutine command_readDataWeights(apr_wght)

! MODS
! PT201014: increased dimensioning of apr_wght from nobs_t to nobs_t*nepochs_t to make it a full vector to match ATransp matrix

  use gracefit_mod
  implicit none


  include 'input.h90'

!********************  Variable declarations ********************

  double precision, intent(out) :: apr_wght(nobs_t*nepochs_t)   ! Weights to be applied to data !@# change to nobs_t

  integer*4  :: isat                  ! Counter variable that runs through satellites
  integer*4  :: i                     ! Counter variable
  double precision :: print_out       ! Indicator on whether the data weights and uncertainties should be printed out
  double precision :: sigmas(nobs_t)  ! Uncertainties for the observations
  integer :: iepoch
!****************************************************************

!****************** DATA WEIGHTS TO BE APPLIED ******************

!   GPS position data weights
  do isat = 1, nsat_t
     do i = igobs,igobs+ngobs_t-1
        call command_storeCommand(LUCMD,'gps_sig_'//satnam(isat),sigmas(i+(isat-1)*ngobs_t),i-igobs+1,ngobs_t,&
             'GRACE '//satnam(isat)//' GPS a priori data weights',1)
     enddo
  enddo
!   GPS antenna CoM conditional equations on position
  do isat = 1, nsat_t
     do i = iCoM,iCoM+nCoM_t-1
        call command_storeCommand(LUCMD,'apr_CoM_gps_off_sig_'//satnam(isat),sigmas(i+(isat-1)*nCoM_t),i-iCoM+1,nCoM_t,&
             'GRACE '//satnam(isat)//' GPS Ant CoM a priori data weights',1)
     enddo
  enddo
!   Data weights for conditions on acceleration during eclipse
  do isat = 1, nsat_t
     do i = icond,icond+ncond_t-1
        call command_storeCommand(LUCMD,'apr_acc_sig_'//satnam(isat),sigmas(i+(isat-1)*ncond_t),i-icond+1,ncond_t,&
             'GRACE '//satnam(isat)//' shadow conditions on accel a priori data weights',1)
     enddo
  enddo

!   Data weights for conditions on accelerometer scales entering/leaving eclipse
  if(use_accel_scale_const /= 0) call command_storeCommand(LUCMD,'acc_scl_cond_sig',sigma_acc,1,1,&
       'Accelerometer scale condition a priori data weights',1)

!   K-Band range data weights
  if(nkbr_t /= 0) call command_storeCommand(LUCMD,'kbr_sig',sigmas(ikbr),1,nkbr_t,'KB range a priori data weight',1)
!   K-Band range rate data weights
  if(nkbrr_t /= 0) call command_storeCommand(LUCMD,'kbrr_sig',sigmas(ikbrr),1,nkbrr_t,&
       'KB range rate a priori data weight',1)
!   K-Band range acceleration data weights
  if(nkbra_t /= 0) call command_storeCommand(LUCMD,'kbra_sig',sigmas(ikbra),1,nkbra_t,&
       'KB range acceleration a priori data weight',1)
!   Mean accelerometer difference data weights
  do i = iacobs,iacobs+naobs_t-1
     call command_storeCommand(LUCMD,'mean_dacc_sig',sigmas(i),i-iacobs+1,naobs_t,&
          'Mean accelerometer difference a priori data weights',1)
  enddo
!****************************************************************

! PT201014: Set data weights for first epoch
  apr_wght = 1.d0
  do i = 1,nobs_t
     apr_wght(i) = 1/sigmas(i)**2
  enddo
! PT201014: now duplicate for all other epochs
  do iepoch = 1, nepochs_t-1
  !Original of PT apr_wght(nobs_t*iepoch+1:(nobs_t+1)*iepoch)
    apr_wght(nobs_t*iepoch+1:nobs_t*(iepoch+1)) = apr_wght(1:nobs_t)
  enddo
!****************************************************************

! Print out apr_wght and uncertainties if requested
  call command_storeCommand(LUCMD,'wght_print_out',print_out,1,1,' ',0)
  if(int(print_out) /= 0 .and. int(print_out) /= 1) &
       call status_update('FATAL','GRACEFIT','command_readDataWeights',' ','print out indicator must be 0 or 1',0)
  if(int(print_out) == 1)then
     write(*,'(a)') ' Uncertainties and data weights on observations:'
     write(*,'(a,a)') ' Observations       units           Uncertainty        ',&
          'Data weight (units the inverse squared of given ones)'
     do i = 1, nobs_t
        write(*,'(a30,f18.9,1x,f23.7)') obsnam(i), sigmas(i), apr_wght(i)
     enddo
  endif
!****************************************************************

  return
end subroutine command_readDataWeights

!********************************************************************************************************************************
!********************************************************************************************************************************
! command_readAprConst: Reads in uncertainties on parameters from command file and places values into
!                       a priori constraints array accordingly
! Author: Thomas Greenspan
!         (taken in large part from sections of read_input.f90 by R. King itself based on subroutine of same name by S. McClusky)
!********************************************************************************************************************************

subroutine command_readAprConst(apr_const,apr_tide_ampl_const)
  use gracefit_mod
  use soln_mod
  implicit none


  include 'input.h90'

  !********************  Variable declarations ********************

  double precision, intent(out) :: apr_const(nparam_t)                    ! A priori constraint values
  double precision, intent(out) :: apr_tide_ampl_const(max_mcon_tides,2)  ! A priori constraint values on each tidal constituent

  integer*4  :: isatt                    ! Counter variable that goes through satellites
  integer*4  :: i                       ! Counter variable
  double precision :: print_out         ! Indicator on whether the constraints and uncertainties should be printed out
  double precision :: apr_mascons       ! A priori mascon constraint value
  double precision :: apr_tide_mascons  ! A priori mascon tidal amplitude constraint value (same for each constituent at this stage)
  double precision :: sigmas(nparam_t)  ! sigmas for observations
  character*200    :: message
  !****************************************************************

  !******************* CONSTRAINTS TO BE APPLIED ******************
  !   IC a priori constraints
  do isatt = 1, nsat_t
     do i = iICparm, iICparm+ngprm_t-1
        call command_storeCommand(LUCMD,'apr_sat_'//satnam(isatt),sigmas(i+(isatt-1)*nsatprm_t),i-iICparm+1,ngprm_t, &
             'GRACE '//satnam(isatt)//' IC a priori constraints',1)
     enddo
  enddo
  !   acceleration scale a priori constraints
  do isatt = 1, nsat_t
     do i = iScale,iScale+nScale_t-1
        call command_storeCommand(LUCMD,'apr_sat_scale_'//satnam(isatt),sigmas(i+(isatt-1)*nsatprm_t),i-iScale+1, &
             nScale_t,'GRACE '//satnam(isatt)//' acceleration scale a priori constraints',1)
     enddo
  enddo
  !   acceleration bias a priori constraints
  do isatt = 1, nsat_t
     do i = iBias,iBias+nBias_t-1
        call command_storeCommand(LUCMD,'apr_sat_bias_'//satnam(isatt),sigmas(i+(isatt-1)*nsatprm_t),i-iBias+1, &
             nBias_t,'GRACE '//satnam(isatt)//' acceleration bias a priori constraints',1)
     enddo
  enddo
  !   once-per-rev a priori constraints
  do isatt = 1, nsat_t
     do i = ionepr,ionepr+n_onepr_t-1
        call command_storeCommand(LUCMD,'apr_1pr_'//satnam(isatt),sigmas(i+(isatt-1)*nsatprm_t),i-ionepr+1, &
             n_onepr_t,'GRACE '//satnam(isatt)//' once-per-rev a priori constraints',1)
     enddo
  enddo
  !   twice-per-rev a priori constraints
  do isatt = 1, nsat_t
     do i = itwopr,itwopr+ntwopr_t-1
        call command_storeCommand(LUCMD,'apr_2pr_'//satnam(isatt),sigmas(i+(isatt-1)*nsatprm_t),i-itwopr+1, &
             ntwopr_t,'GRACE '//satnam(isatt)//' twice-per-rev a priori constraints',1)
     enddo
  enddo
  !   antenna offset a priori constraints
  do isatt = 1, nsat_t
     do i = iantoff,iantoff+nantoff_t-1
        call command_storeCommand(LUCMD,'apr_sat_gps_off_'//satnam(isatt),sigmas(i+(isatt-1)*nsatprm_t),i-iantoff+1, &
             nantoff_t,'GRACE '//satnam(isatt)//' antenna offset a priori constraints',1)
     enddo
  enddo
  !   roll/pitch/yaw offset a priori constraints
  do isatt = 1, nsat_t
     do i = irpy,irpy+nrpy_t-1
        call command_storeCommand(LUCMD,'apr_sat_rpy_off_'//satnam(isatt),sigmas(i+(isatt-1)*nsatprm_t),i-irpy+1, &
             nrpy_t,'GRACE '//satnam(isatt)//' roll/pitch/yaw offset a priori constraints',1)
     enddo
  enddo
  !   mascon a priori constraint
  if(nmascons_t /= 0) then
     call command_storeCommand(LUCMD,'apr_mascons',apr_mascons,1,1,'Mascon Constraint',1)
     sigmas(imascons:imascons+nmascons_t-1) = apr_mascons
  endif
  !SA ADD "then" "endif"
  ! PT140618: mascon tidal amplitude constraint - read a default value for the same constraint on each consituent
  if(est_msc_tides /= 0) then
     call command_storeCommand(LUCMD,'apr_tide_mascons',apr_tide_mascons,1,1,'Mcon Tidal Ampl. Constraint',1)
     print *, 'READ APR TIDE' ,apr_tide_mascons, size(sigmas), imsctide, imsctide+nmsc_tid_constit_t*2-1
     sigmas(imsctide:imsctide+nmsc_tid_constit_t*2-1) = apr_tide_mascons
  endif
  ! PT150826: now also read in specific constraints for each constituent
  if(est_msc_tides /= 0) then
     !     M2
     do i=1,2
        call command_storeCommand(LUCMD,'apr_M2_ampl',apr_tide_ampl_const(1,i),i,2,' ',0)
     enddo
     !     O1
     do i=1,2
        call command_storeCommand(LUCMD,'apr_O1_ampl',apr_tide_ampl_const(2,i),i,2,' ',0)
     enddo
     !     S2
     do i=1,2
        call command_storeCommand(LUCMD,'apr_S2_ampl',apr_tide_ampl_const(3,i),i,2,' ',0)
     enddo
     !     K1
     do i=1,2
        call command_storeCommand(LUCMD,'apr_K1_ampl',apr_tide_ampl_const(4,i),i,2,' ',0)
     enddo
     !     K2
     do i=1,2
        call command_storeCommand(LUCMD,'apr_K2_ampl',apr_tide_ampl_const(5,i),i,2,' ',0)
     enddo
     ! Output information on constraints
     write(message,'(a,2e12.3,a)')'M2 apr constraints (sin, cos):',apr_tide_ampl_const(1,:)," m"
     call status_update('STATUS','GRACEFIT','command_readDAprConst',' ',message,0)
     write(message,'(a,2e12.3,a)')'O1 apr constraints (sin, cos):',apr_tide_ampl_const(2,:)," m"
     call status_update('STATUS','GRACEFIT','command_readDAprConst',' ',message,0)
     write(message,'(a,2e12.3,a)')'S2 apr constraints (sin, cos):',apr_tide_ampl_const(3,:)," m"
     call status_update('STATUS','GRACEFIT','command_readDAprConst',' ',message,0)
     write(message,'(a,2e12.3,a)')'K1 apr constraints (sin, cos):',apr_tide_ampl_const(4,:)," m"
     call status_update('STATUS','GRACEFIT','command_readDAprConst',' ',message,0)
     write(message,'(a,2e12.3,a)')'K2 apr constraints (sin, cos):',apr_tide_ampl_const(5,:)," m"
     call status_update('STATUS','GRACEFIT','command_readDAprConst',' ',message,0)
  endif


  !**************************************************************

  ! Set a priori constraints
  do i = 1,nparam_t
     apr_const(i) = 1/sigmas(i)**2
  enddo
  !**************************************************************

  ! Print out apr_const and uncertainties if requested
  call command_storeCommand(LUCMD,'const_print_out',print_out,1,1,' ',0)
  if(int(print_out) /= 0 .and. int(print_out) /= 1) &
       call status_update('FATAL','GRACEFIT','command_readDAprConst',' ','print out indicator must be 0 or 1',0)
  if(int(print_out) == 1)then
     write(*,'(a)') ' Uncertainties and constraints on parameters:'
     write(*,'(a)') ' Parameters     units       Uncertainty         Constraints (units the inverse squared of given ones)'
     do i = 1, nparam_t
        write(*,'(a24,f16.12,3x,f20.7)') prmnam(i), sigmas(i), apr_const(i)
     enddo
  endif
  !****************************************************************

  return
end subroutine command_readAprConst

!********************************************************************************************************************************
!********************************************************************************************************************************
! command_readTol: Read in and store orbit misfit and kbrr prefit tolerances
! Author: Thomas Greenspan
!         (taken in large part from sections of read_input.f90 by R. King itself based on subroutine of same name by S. McClusky)
!********************************************************************************************************************************

subroutine command_readTol(kbrr_prefit_tol)
  use gracefit_mod
  implicit none


  include 'input.h90'

  !********************  Variable declarations ********************

  double precision, intent(out) :: kbrr_prefit_tol ! Maximal kbrr prefit tolerance
  character*100 message

  !****************************************************************

  if (nkbrr_t /= 0) then
     call command_storeCommand(LUCMD,'kbrr_prefit_tol',kbrr_prefit_tol,1,1,'KBRR prefit tolerance',1)
     write(message,'(a,1x,f18.2,1x,a)') 'KBRR prefit tolerance is: ',kbrr_prefit_tol*m_um,' (um)'
     call status_update('STATUS','GRACEFIT','command_readTol',' ',message,0)
  endif
  !**************************************************************

  return
end subroutine command_readTol

!********************************************************************************************************************************
!********************************************************************************************************************************
! command_printSetup: Prints out setup if requested in command file
! Author: Thomas Greenspan
!********************************************************************************************************************************

subroutine command_printSetup()
  use gracefit_mod
  use soln_mod
  implicit none


  include 'input.h90'

  !********************  Variable declarations ********************

  integer*4        :: i            ! Counter variable
  double precision :: run_program  ! Indicator as to whether to run the program or print out the setup
  character(16)    :: kbnam        ! Names of kband observations
  !****************************************************************

  ! Get command from file
  call command_storeCommand(LUCMD,'run_program',run_program,1,1,' ',0)
  if(int(run_program) /= 0 .and. int(run_program) /= 1) &
       call status_update('FATAL','GRACEFIT','command_printSetup',' ','run_program indicator must be 0 or 1',0)
  !**************************************************************

  ! If running the program requested, return without printing anything
  if(int(run_program) == 1) return

  write(*,'(a)') ' '
  write(*,'(a)') ' Printing setup for GRACEFIT '
  write(*,'(a)') ' '
  write(*,'(a)') ' '
  write(*,'(a,i6,a,f4.1,a)') ' Total number of epochs: ',nepochs_t,'       (hours spanned = ',spanhrs,')'
  write(*,'(a,i12,a)')  ' Starting epoch: ',int(starting_epoch),' (grace-seconds)'
  write(*,'(a,i4,a)') ' Epoch interval: ',int(epoch_interval),' (s)'
  write(*,'(a,i6,a,i4,a)') ' Maximum number of GPS observations used: ',nGPSobs_t(1),&
       '    (GPS interval : ',int(gnv_interval*epoch_interval), ' s)'
  write(*,'(a,i6)') ' Minimum number of GPS observations allowed: ',minGPSobs
  write(*,'(a,i6)') ' Maximum number of Kband observations used: ',nKBobs_t
  write(*,'(a,i6)') ' Minimum number of Kband observations allowed: ',minKBobs
  write(*,'(a)') ' '
  write(*,'(a,i6)') ' Total number of observations: ',nobs_t
  if(ngobs_t > 0) write(*,'(a,i4)')  '    G-observations per satellite: ',ngobs_t
  if(nCoM_t > 0)  write(*,'(a,i4)')  '    Antenna offset observations per satellite: ',nCoM_t
  if(ncond_t > 0)  write(*,'(a,i4)') '    Shadow condition observations per satellite: ',ncond_t
  kbnam = ' '
  if(nkbr_t > 0) kbnam = ' R'
  if(nkbrr_t > 0) kbnam = kbnam(1:nkbr_t*2)//' RR'
  if(nkbra_t > 0) kbnam = kbnam(1:nkbr_t*2+nkbrr_t*3)//' RA'
  kbnam = kbnam(1:nkbr_t*2+nkbrr_t*3+nkbra_t*3)//')'
  kbnam = adjustl(kbnam)
  if(nkobs_t > 0) write(*,'(a,i4,a,a)')  '    Kband observations: ',nkobs_t,'   (',kbnam
  if(naobs_t > 0)  write(*,'(a,i4)') '    Mean acceleration condition observations: ',naobs_t
  write(*,'(a,i4)') '    Total satellite-specific observations: ',nsatobs_t
  write(*,'(a)') ' '
  write(*,'(a,i4)') ' Total number of parameters:   ',nparam_t
  write(*,'(a,i4)') '    IC parameters:   ',ngprm_t
  write(*,'(a,i4)') '    Scale parameters:   ',nScale_t
  write(*,'(a,i4)') '    Bias parameters:   ',nBias_t
  write(*,'(a,i4)') '    Antenna offset parameters:   ',nantoff_t
  write(*,'(a,i4)') '    Total orbital parameters:   ',norbprm_t
  write(*,'(a,i4)') '    Total satellite-specific parameters:   ',nsatprm_t
  write(*,'(a,i4)') '    Mascon parameters:   ',nmascons_t
  write(*,'(a)') ' '
  write(*,'(a)') ' Observations:'
  do i = 1, nobs_t
     write(*,'(a)')obsnam(i)
  enddo
  write(*,'(a)') ' '
  write(*,'(a)') ' Parameters:'
  do i = 1, nsatprm_t*nsat_t
     write(*,'(a)')prmnam(i)
  enddo
  write(*,'(a,i6)') ' Number of mascons: ',nmascons_t
  write(*,'(a)') ' '
  write(*,'(a)') ' Normal end of program '

  ! Close command file
  close(LUCMD)
  stop
  !**************************************************************

end subroutine command_printSetup

!********************************************************************************************************************************
!********************************************************************************************************************************
! command_close: close command file
! Author: Thomas Greenspan
!********************************************************************************************************************************

subroutine command_close()

  implicit none

  include 'input.h90'

  !****************************************************************

  ! close file
  close(LUCMD)
  !**************************************************************

end subroutine command_close

!********************************************************************************************************************************
