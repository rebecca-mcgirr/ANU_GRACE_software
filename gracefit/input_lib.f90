!********************************************************************************************************************************
!  File: input_lib.f90
!
!  Purpose: Set of subroutines used to open and read in data
!
!  Author: Thomas Greenspan
!          (some subroutines taken or inspired from subroutines by Simon McClusky or Paul Tregoning)
!
!  API:
!       input_openFile      : Open a file
!       input_openCommand   : Opens command file
!       input_openh5GTORB   : Open hdf5 GTORB file(s)
!       input_openGTORB     : Open binary GTORB file(s)
!       input_openAllFiles  : Open all files needed as input for GRACEFIT
!       input_skipHeader    : Skip header of a file
!       input_readGTORBh5   : Read in data from hdf5 GTORB file(s)
!       input_readGTORB     : Read in data from binary GTORB file(s)
!       input_readACC1B     : Read in data from ACC1B file(s)
!       input_readTHR1B     : Read in data from THR1B file(s)
!       input_readGNV1B     : Read in data from GNV1B file(s)
!       input_readKBR1B     : Read in data from KBR1B file
!       input_isValid       : Validate that the data read in (makes sure satellites have the same epoch
!                             and that it corresponds correctly to the index)
!       input_readGTORBh5_v1: Read in data from v1 hdf5 GTORB file(s)
!       input_readGTORBh5_Amat: Read in data from v1 hdf5 GTORB file(s)
!
!  July 11, 2013
!
!********************************************************************************************************************************
!********************************************************************************************************************************
! input_openFile: Opens a file and writes a status report given the file unit number, name and status as
!                 well as the name of the program and subroutine, the type of error if there is one and
!                 the message written in case of an error.
! Author: Thomas Greenspan
! MODS
! PT131104: added access_tpye and rec_len so that binary, direct access files can be opened as well.
!********************************************************************************************************************************

subroutine input_openFile(file_unit_num,file_name,file_status,program_name,subroutine_name,access_type,rec_len &
     ,error_type,error_message)

  implicit none

  !********************  Variable declarations ********************
  integer*4,    intent(in) :: file_unit_num   ! File unit number to be opened
  character(*), intent(in) :: file_name       ! Name of file to be opened
  character(*), intent(in) :: file_status     ! Status of file to be opened
  character(*), intent(in) :: program_name    ! Program name
  character(*), intent(in) :: subroutine_name ! Subroutine name
  character(*), intent(in) :: access_type     ! "sequential" or "direct    "
  integer*4,    intent(in) :: rec_len         ! record length for binary files (zero otherwise, and not used)
  character(*), intent(in) :: error_type      ! Type of error if any occurs
  character(*), intent(in) :: error_message   ! Message in case of error

  integer*4 :: ioerr      ! Standard error variable
  !****************************************************************


  ! Open file and report stat
  if(access_type == "sequential")then
     open(unit=file_unit_num,file=file_name,status=file_status,access='sequential',iostat=ioerr)
  else
     open(unit=file_unit_num,file=file_name,status=file_status,access='direct',recl=rec_len,iostat=ioerr)
  endif
  if (ioerr /= 0) call status_update(error_type,program_name,subroutine_name,file_name,error_message,ioerr)
  call status_update('STATUS',program_name,subroutine_name,trim(file_name),'Opened file: ',ioerr)
  !****************************************************************

  return
end subroutine input_openFile

!********************************************************************************************************************************
!********************************************************************************************************************************
! input_openCommand: Open the command file
! Author: Thomas Greenspan
!********************************************************************************************************************************

subroutine input_openCommand(cmdnam)

  implicit none

  include 'input.h90'

  !********************  Variable declarations ********************

  character(*), intent(in) :: cmdnam  ! Name of command file
  !****************************************************************

  ! Open the command file
  call input_openFile(LUCMD,cmdnam,'old','GRACEFIT','input_openCommand',"sequential",0,'FATAL','Error opening the command file: ')

  !****************************************************************

  return
end subroutine input_openCommand

!********************************************************************************************************************************
!********************************************************************************************************************************
! input_openCommand: Open the command file
! Author: Thomas Greenspan
!********************************************************************************************************************************

subroutine input_openh5GTORB(gt_fname, luoffset)
  use gracefit_mod
  implicit none
   
  include 'input.h90'
  character(len=*), intent(inout) :: gt_fname(maxsat)
  integer :: luoffset
  if (SAT_1 == 1) then
     call orb_open(gt_fname(1), orb_files(1+luoffset))
  else
     gt_fname(1) = 'N/A'
  endif
  if (SAT_1 == 1) then
     call orb_open(gt_fname(2), orb_files(2+luoffset))
  else
     gt_fname(2) = 'N/A'
  endif
end subroutine input_openh5GTORB



subroutine input_openGTORB(gt_fnam,lu_offset,access_type,rec_len)

  use gracefit_mod
  implicit none
   
  include 'input.h90'

  !********************  Variable declarations ********************

  character(*), intent(inout)  :: gt_fnam(maxsat)    ! Name of GTORB files
  integer*4,    intent(in)     :: lu_offset          ! offset to indicate whether theoretical or truth files (affects unit numbers used)
  character(*),intent(in)      :: access_type        ! either "sequential" or "direct    "
  integer*4  ,intent(out)      :: rec_len(2)         ! record length for binary, direct access files (= 219348)

  !****************************************************************

  ! PT131104: define the record length for direct access files
  if (access_type == "sequential") then
     rec_len = 0
  endif

  ! open the GRACE GTORB file(s)
  if(SAT_1 == 1)then
     if(access_type(1:6) == "direct")call GTORB_record_length(LUGTORB(1+lu_offset),calling_prog,gt_fnam(1),rec_len(1))
     call input_openFile(LUGTORB(1+lu_offset),gt_fnam(1),'old',calling_prog,'input_openGTORB',access_type,rec_len(1),'FATAL',&
          'Error opening the GRACE A GTORB file: ')
  else
     gt_fnam(1) = 'N/A'
  endif
  if(SAT_2 == 1)then
     if(access_type(1:6) == "direct")call GTORB_record_length(LUGTORB(nsat_t+lu_offset),calling_prog,gt_fnam(2),rec_len(2))
     call input_openFile(LUGTORB(nsat_t+lu_offset),gt_fnam(2),'old',calling_prog,'input_openGTORB',access_type,rec_len(2),'FATAL',&
          'Error opening the GRACE B GTORB file: ')
  else
     gt_fnam(2) = 'N/A'
  endif
  !****************************************************************

  return
end subroutine input_openGTORB

!********************************************************************************************************************************
!********************************************************************************************************************************
! input_openAllFiles: Opens all the files needed for the program gracefit given the date string, the name
!                     of the command file and the name of the GTORB files given as input.
! Author: Thomas Greenspan
!         (taken in part from subroutine open_orbfiles.f90)
!********************************************************************************************************************************

subroutine input_openAllFiles(date_str)
  use gracefit_mod         
  implicit none
   
  include 'input.h90'

  !********************  Variable declarations ********************

  character(10), intent(in) :: date_str     ! String containing date that was input

! local variables
  integer*4     :: isat                     ! Counter that goes through satellites
  character(64) :: ut1fil,polfil,nutfil     ! Names of ut1, pole and nutation tables respectively
  character*4   :: suffix                   ! end of the L1B 
  character*1   :: CH1                      ! character in KBR1B file name (X: GRACE; Y: GRACE FO; Z? GRACE II)
  character*5   :: acc_code                 ! name of accelerometer file (ACC1B for GRACE, ACT1B for GRACE FO)
  !****************************************************************

   character(len=150), dimension(2) :: GNV,SCA,ACC,THR 
   character (len=150) :: KBR, LRI
   integer :: uin ,ioerr
   namelist /l1b_files/ GNV,SCA,ACC,THR, KBR, LRI
   uin = 0
! PT190903: check if file exists and stop fatally if not
   open(newunit=uin, file='L1B_files.txt',status='old',iostat=ioerr)
   if(ioerr /= 0)then
     call status_update('FATAL',calling_prog,'input_openallfiles','L1B_files.txt',"Error opening file. Does it exist?",0)
   endif
   read(uin, nml=l1b_files)
   gn_fnam(1)=gnv(1)
   gn_fnam(2)=gnv(2)
   acc_fnam(1) = acc(1)
   acc_fnam(2) = acc(2)
   thr_fnam(1) = thr(1)
   thr_fnam(2) = thr(2)
   kb_fnam = kbr
   lri_fnam = lri
!   ! set the suffix for the file name based on the mission
!   if(mission == 0)then
!     suffix = ".asc"
!     RL_NUM = 02
!     CH1 = "X"
!     acc_code = "ACC1B"
!   else if (mission == 1) then
!     suffix = ".txt"
! ! PT190526: GRACEFO data have been released as RL04
!     RL_NUM = 04
!     CH1 = "Y"
!     acc_code = "ACT1B"
!   else if (mission == 2) then
!     suffix = ".txt"
!     RL_NUM = 99
!     CH1 = "Z"
!     acc_code = "ACC1B"
!   endif

!   !   Form the file names for GNV1B and KBR1B
!   do isat = 1, nsat_t
!      write(gn_fnam(isat), '(a,a10,a,a1,a1,i2.2,a)') "GNV1B_",date_str,"_",satnam(isat),"_",RL_NUM,suffix
!      write(acc_fnam(isat),'(a5,a1,a10,a,a1,a1,i2.2,a)') acc_code,"_",date_str,"_",satnam(isat),"_",RL_NUM,suffix
!      write(thr_fnam(isat),'(a,a10,a,a1,a1,i2.2,a)') "THR1B_",date_str,"_",satnam(isat),"_",RL_NUM,suffix
!   enddo
!   write(kb_fnam,'(a,a10,a,a1,a1,i2.2,a)') "KBR1B_",date_str,"_",CH1,"_",RL_NUM,suffix

!   write(lri_fnam,'(a,a10,a,a1,a1,i2.2,a)')"LRI1B_",date_str,"_",CH1,"_",RL_NUM,suffix
  !  Create file names for the ut1, pole, and nutation files
! PT200706: delete all references to these files!
!  ut1fil = 'ut1.'
!  polfil = 'pole.'
!  nutfil = 'nutabl.'
  !************************** OPEN FILES **************************

  do isat = 1, nsat_t
     !     open level 1B GPS position file(s)
     call input_openFile(LUGN(isat),gn_fnam(isat),'old','GRACEFIT','input_openAllFiles',"sequential",0,'FATAL',&
          'Error opening the GRACE '//satnam(isat)//' GNV1B file: ')
     !     open level 1B GRACE ACC file(s)
     call input_openFile(LUAC(isat),acc_fnam(isat),'old','GRACEFIT','input_openAllFiles',"sequential",0,&
          'FATAL','Error opening the GRACE '//satnam(isat)//' ACC1B file: ')
     !     open level 1B GRACE ACC file(s)
     call input_openFile(LUTH(isat),thr_fnam(isat),'old','GRACEFIT','input_openAllFiles',"sequential",0,&
          'FATAL','Error opening the GRACE '//satnam(isat)//' THR1B file: ')
  enddo
  !   open level 1B GRACE KBR file
  call input_openFile(LUKB,kb_fnam,'old','GRACEFIT','input_openAllFiles',"sequential",0,'FATAL','Error opening the KBR1B file: ')

! PT180625: if GRACE FO, open the LRI file
  if(mission == 1)then
print*,"Do not (yet) open a LRI file for GRACE FO - we don't have one!!!"
    !call input_openFile(LULRI,lri_fnam,'old','GRACEFIT','input_openAllFiles',"sequential",0,'FATAL','Error opening the LRI1B file: ')
  endif

! PT200706: delete these
!  ! Open the ut1, pole, and nutation files (error FATAL because they are needed for computation of sun position)
!  !   open ut1 file
!  call input_openFile(iUT1,ut1fil,'old','GRACEFIT','input_openAllFiles',"sequential",0,'FATAL','Error opening ut1 table: ')
!  !   open pole file
!  call input_openFile(iPOLE,polfil,'old','GRACEFIT','input_openAllFiles',"sequential",0,'FATAL','Error opening pole table: ')
!  !   open nutation file
!  call input_openFile(iNUT,nutfil,'old','GRACEFIT','input_openAllFiles',"sequential",0,'FATAL','Error opening nutation table: ')
  !****************************************************************

  return
end subroutine input_openAllFiles

!********************************************************************************************************************************
!********************************************************************************************************************************
! input_skipHeader: Skips the header of a file given its unit number and its name.
! Author: Thomas Greenspan
!********************************************************************************************************************************

subroutine input_skipHeader(file_unit_num,file_name)

  !********************  Variable declarations ********************

  integer*4, intent(in)    :: file_unit_num   ! Unit number of file whose header is going to be skipped
  character(*), intent(in) :: file_name       ! Name of file whose header is going to be skipped

  character(150) :: line         ! Variable where each line is stored until the END OF HEADER is found
  character(150) :: message      ! Message to be written in case of failure
  !****************************************************************

  !   read through header
  line = ' '
  do while (line(1:13) /= 'END OF HEADER')
     read(file_unit_num,'(a)',end=13)line(1:13)
  enddo
  !****************************************************************

  return   ! return if successful

13 write(message,'(a,a,a)')'No end of header found reading ',file_name,'file'
  call status_update('FATAL','GRACEFIT','input_lib',' ',message,0)

  return
end subroutine input_skipHeader



subroutine input_readGTORBh5(satics_t,apr_ScaleBias,GPSantoff,rvec,apr_prm &
     ,sciframe_quat,srf2trf_rotmat,srf2trf_deriv_rotmat,lu_offset, mcon_ocean &
     ,msc_tide_amp, mcon_tides)

  use gracefit_mod
  implicit none

   

  !********************  Variable declarations ********************

  double precision, intent(in)  :: satics_t(maxrvecprm+1,nsat_t)              ! A priori values for IC's (read from header)
  double precision, intent(in)  :: apr_ScaleBias(max_SB,nsat_t)               ! A priori values for bias and scale (read from header)
  double precision, intent(in)  :: GPSantoff(7,maxsat)                        ! Quaternion and offsets for antenna
  double precision, intent(out) :: rvec(nrvec_t,nsat_t,nepochs_t)             ! Pos, vel and partials for each satellite at every epoch
  double precision, intent(out) :: apr_prm(nparam_t)                          ! A priori values for parameters
  double precision, intent(out) :: sciframe_quat(4,nsat_t,nepochs_t)          ! Quaternions for each satellites at every epoch
  double precision, intent(out) :: srf2trf_rotmat(3,3,nepochs_t,nsat_t)       ! SRF to TRF rotation matrix for GRACE A and B
  double precision, intent(out) :: srf2trf_deriv_rotmat(3,3,nepochs_t,nsat_t) ! Differentiated SRF to TRF rotation matrix for A and B
  integer*4,        intent(in)  :: lu_offset                                  ! offset to indicate whether theoretical or truth files (affects unit numbers used)
  logical,          intent(out) :: mcon_ocean(nmascons_t)                     ! logical array to indicate whether an ocean or land mascon
  real(kind=8)    , intent(out) :: msc_tide_amp(max_mcon_tides,2,nmascons_t)  ! a priori mascon tidal amplitudes for each component/constituent/mascon
  integer*4        ,intent(in)  :: mcon_tides(4556,2)                         ! bit-mapped indicator of which tidal amplitudes to estimate per mascon


  integer*4        :: isat                    ! Counter for satellites
  integer*4        :: iepoch                  ! Counter for epochs
  integer*4        :: i,j,k                   ! Counter variables
  integer*4        :: mascon_num              ! Number of mascons
  integer*4        :: tmp_epoch               ! temporary epoch variable to then be converted to R*8
  double precision :: epoch(nsat_t,nepochs_t) ! Epochs recorded in GTORB file (only used to validate data)
  character(150)   :: line                    ! Line to be read 
  character(150)   :: message                 ! Message printed when done reading file
  integer*4        :: ioerr                   ! Standard error variable
  integer*4        :: irec_eoh                ! storing the input irec value for use with the second GTORB file
  integer*4        :: mascon_tidal_amp_num    ! value read off the bottom of the GTORB file
  integer*4        :: count_mcon_tides        ! counter used to transfer apriori mcon tidal amplitudes to apr_prm variable
  logical          :: bitmap

  ! PT131024: variables to segment and store the parials according to the parameters being estimated
  integer*4        :: pntr_rvec               ! pointer to which row we are up to in rvec, based on which parameters are estimated
  double precision, allocatable :: tmprvec(:,:) ! temporary storage of all partials of a single epoch, before transferring reqired values to rvec

  integer :: local_offset

  local_offset = lu_offset


  ! Header already read in beforehand
  ! Set IC a priori values based on data taken from header
  allocate(tmprvec(36+18+18+24+24+nmascons_t*6 + nmsc_tid_constit_t*2*6, nepochs_t))
  apr_prm = 0.d0
  do isat = 1, nsat_t
     do i = iICparm, iICparm+ngprm_t-1
        apr_prm(i+(isat-1)*nsatprm_t) = satics_t(i-iICparm+1+1,isat)
     enddo
     do i = iScale, iScale+nScale_t-1
        apr_prm(i+(isat-1)*nsatprm_t) = apr_ScaleBias(i-iScale+1,isat)
     enddo
     do i = iBias, iBias+nBias_t-1
        apr_prm(i+(isat-1)*nsatprm_t) = apr_ScaleBias(i-iBias+3+1,isat)
     enddo
     ! PT140821: added code for a priori values for once- and twice-per-rev along-track acceleration
     do i = ionepr, ionepr+n_onepr_t-1
        apr_prm(i+(isat-1)*nsatprm_t) = satics_t(i-ionepr+iICparm+ngprm_t+1,isat)
     enddo
     do i = itwopr, itwopr+ntwopr_t-1
        apr_prm(i+(isat-1)*nsatprm_t) = satics_t(i-itwopr+iICparm+ngprm_t+3,isat)
     enddo
     do i = iantoff, iantoff+nantoff_t-1
        apr_prm(i+(isat-1)*nsatprm_t) = GPSantoff(i-iantoff+1+4,isat)
     enddo
  enddo

  ! store the original pointer of the input record number (= the "END OF HEADER" record)

  !*********************** READ GTORB FILES ***********************

  ! Read in positions, velocities and partials as well as quaternions for each satellite every epoch
  ioerr = 0
  do isat = 1, nsat_t
        call status_update('STATUS','GRACEFIT','input_readGTORBh5',' ','reading orbit coords and partials',0)
     call orb_read(orb_files(isat+local_offset), time =epoch(isat,:), coord = rvec(1:nrvecprm_t,isat,:),  &
          sciframe_quat = sciframe_quat(:,isat,:), partials = tmprvec(1:114,:) )
     if (nmascons_t > 0 .and. lu_offset == 0) then
        call status_update('STATUS','GRACEFIT','input_readGTORBh5',' ','reading mascon partials',0)
        call orb_read(orb_files(isat+local_offset), msc_part = tmprvec(115:114+nmascons_t*6,:) )
     endif
     if (nmsc_tid_constit_t > 0) then
        call status_update('STATUS','GRACEFIT','input_readGTORBh5',' ','reading tidal amplitude partials',0)
        call orb_read(orb_files(isat+local_offset), &
             msc_part = tmprvec(114+nmascons_t*6+1:114+nmascons_t*6+nmsc_tid_constit_t*12,:))
     endif

     ! now, transfer into rvec the partials of the parameters to be estimated
     ! they are in tmprvec in order: 36x partials wrt ICs, 18x partials wrt acc bias, 18x partials wrt scls, 6x partials wrt each mascon, 6x2 partials wrt each tidal constituent amplitude

     if(ngprm_t == 3 ) then  ! GPS positions to be estimated but not velocity
        rvec(nrvecprm_t+1:nrvecprm_t+3 ,isat,:) = tmprvec(1:3,:)      ! X wrt XYZ                     
        rvec(nrvecprm_t+4:nrvecprm_t+6,isat,:) = tmprvec(7:9,:)       ! Y wrt XYZ   
        rvec(nrvecprm_t+7:nrvecprm_t+9,isat,:) = tmprvec(13:15,:)     ! Z wrt XYZ 
        pntr_rvec = nrvecprm_t+9
     else if(ngprm_t == 6) then  ! GPS positions and velocities to be estimated
        rvec(nrvecprm_t+1:nrvecprm_t+36 ,isat,:) = tmprvec(1:36,:)     ! all pos/vel wrt all pos/vel   
        pntr_rvec = nrvecprm_t+36
     endif
     ! do we need partials for acc biases?
     if(nBias_t == 3) then
        rvec(pntr_rvec+1:pntr_rvec+18,isat,:) = tmprvec(37:54,:)
        pntr_rvec = pntr_rvec + 18
     endif

     ! do we need partials for acc scales?
     if(nScale_t == 3) then
        rvec(pntr_rvec+1:pntr_rvec+18,isat,:) = tmprvec(55:72,:)
        pntr_rvec = pntr_rvec + 18
     endif

     ! PT140821: do we ned partials for once-per-rev (SINE and COSINE)?
     if(n_onepr_t == 2) then
        rvec(pntr_rvec+1:pntr_rvec+12,isat,:) = tmprvec(73:84,:)
        pntr_rvec = pntr_rvec + 12
     endif

     ! PT140821: do we ned partials for twice-per-rev (SINE and COSINE)?
     if(ntwopr_t == 2) then
        rvec(pntr_rvec+1:pntr_rvec+12,isat,:) = tmprvec(85:96,:)
        pntr_rvec = pntr_rvec + 12
     endif

     ! PT170410: do we need partials for roll/pitch/yaw ?
     if(nrpy_t == 3) then
        rvec(pntr_rvec+1:pntr_rvec+18,isat,:) = tmprvec(97:114,:)
        pntr_rvec = pntr_rvec + 18
     endif

     ! now store the partials of mascon wrt pos (and maybe vel)
     if(nmascons_t > 0 .and. lu_offset == 0)then
       call status_update('STATUS','GRACEFIT','input_readGTORBh5',' ','now store the partials of mascon wrt pos (and maybe vel)',0)
        if(ngprm_t == 3)then
           do j=1,nmascons_t
              rvec(pntr_rvec+1+(j-1)*3:pntr_rvec+3+(j-1)*3,isat,:) = tmprvec(97+(j-1)*3:99+(j-1)*3,:)    ! XYZ wrt mascons (don't transfer vel partials wrt mascons)
              pntr_rvec = pntr_rvec + 3
           enddo
        elseif(ngprm_t == 6)then
           ! PT140904: changed the index into tmprvec from "73+" to "97+" to account for the presence of the 1/rev and 2/rev partials now in the GTORB files
           ! PT170410: changed again from "97+" tp "115+" to account for roll/pitch/yaw
           rvec(pntr_rvec+1:pntr_rvec+nmascons_t*6,isat,:) = tmprvec(115:115+nmascons_t*6 - 1,:)  ! both position and velocity mascon partials
           pntr_rvec = pntr_rvec + nmascons_t*6
        endif
     endif

     call status_update('STATUS','GRACEFIT','input_readGTORBh5',' ','now tidal amplitudes',0)
     ! PT140617: now store the partials of mascon tide amplitudes wrt pos (and maybe vel)
     if(est_msc_tides > 0)then
        if(ngprm_t == 3)then
           do j=1,nmsc_tid_constit_t*2
              ! PT150825: the indices into tmprvec were completely wrong here ... they didn't account for the 1/rev and 2/rev as per the mascon code above, nor that it should be nmascons_t*6
              rvec(pntr_rvec+1+(j-1)*3:pntr_rvec+3+(j-1)*3,isat,:) = &
                   tmprvec(115+nmascons_t*3+(j-1)*3:117+nmascons_t*3+(j-1)*3,:)    ! XYZ wrt mascons (don't transfer vel partials wrt mascons)
              pntr_rvec = pntr_rvec + 3
           enddo
        elseif(ngprm_t == 6)then
           ! PT150825: changed the index here from "73+" to "97+" as well ... should have done that a year ago when I did it for the mascons above!!
           ! PT170410: updated for the roll/pitch/yaw (so now starts at 115, not 97)
           rvec(pntr_rvec+1:pntr_rvec+nmsc_tid_constit_t*2*6,isat,:) = &
                tmprvec(115+nmascons_t*6:117+nmascons_t*6+nmsc_tid_constit_t*2*6 - 1,:)  ! both position and velocity mascon partials
           ! DEBUG
           !  print*,'rvec(M2s:M2c,isat,iepoch)',rvec(pntr_rvec+1:pntr_rvec+2,isat,iepoch)
        endif
     endif
     !  print*,'after mascon tides pntr_rvec = ',pntr_rvec+nmsc_tid_constit_t*2*6

 !    write(*,*) epoch(1,1:10), epoch(2,1:10)
 !    write(*,*) epoch(1,nepochs_t-9:nepochs_t), epoch(2,nepochs_t-9:nepochs_t)
     call status_update('STATUS','GRACEFIT','input_readGTORBh5',' ','now quaternians',0)
     do iepoch = 1 , nepochs_t
        call rotation_quat2mat_3d (sciframe_quat(:,isat,iepoch),srf2trf_rotmat(:,:,iepoch,isat))
        ! Validate recorded data once both satellites have been read in
        if(isat == nsat_t) call input_isValid(iepoch, epoch(1,iepoch), epoch(nsat_t,iepoch),'input_readGTORBh5','FATAL')
     enddo


!!! Apriori transfer to be done here! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! ******** Mascon  a priori values  ************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Read in mascon a priori values if necessary
     if(nmascons_t > 0)then
       call status_update('STATUS','GRACEFIT','input_readGTORBh5',' ','now msc_apriori values',0)
        call orb_read(orb_files(isat+local_offset), msc_apr =  apr_prm(imascons:imascons+nmascons_t-1) )
        do i=1,nmascons_t
           if(apr_prm(imascons+i-1) > 900.d0)then   ! it is an ocean mascon
              apr_prm(imascons+i-1) = apr_prm(imascons+i-1) - 1000.d0
              mcon_ocean(i) = .true.
           else
              mcon_ocean(i) = .false.
           endif
        enddo
     endif

     if(est_msc_tides > 0) then


        call orb_read(orb_files(isat+local_offset), tmsc_apr =  msc_tide_amp )

        call status_update('STATUS','GRACEFIT','input_readGTORBh5',' ',' Read mascon tidal amplitudes',0)

        do i=1,nmascons_t
           do j=1,max_mcon_tides
              if(bitmap(mcon_tides(i,1),j))then
                 do k=1,2
                    count_mcon_tides = count_mcon_tides + 1
                    apr_prm(imsctide+count_mcon_tides -1) = msc_tide_amp(j,k,i)
                 enddo
              endif
           enddo
        enddo
        write(message,'(a,i7,a)')' Transferred ',count_mcon_tides,' a priori mascon tidal amplitudes to apr_prm'
        call status_update('STATUS','GRACEFIT','input_readGTORBh5',' ',message,0)
     endif
  end do

!Numerically differentiate the SRF to TRF rotation matrix elements
  call status_update('STATUS','GRACEFIT','input_readGTORBh5',' ','differentiate SRF to TRF matrices',0)
  do isat = 1, nsat_t
     do i = 1, 3
        do j = 1, 3
           call noise_robust_deriv(srf2trf_rotmat(i,j,:,isat),srf2trf_deriv_rotmat(i,j,:,isat) &
                ,epoch_interval,nepochs_t,5)
        enddo
     enddo
  enddo


end subroutine input_readGTORBh5












!********************************************************************************************************************************
!********************************************************************************************************************************
! input_readGTORB: Reads in all needed data from GTORB files and calculates both the rotation matrix
!                  from SRF to TRF and its numerical derivative. Input is the vector in which positions,
!                  velocities and partials will be stored as well as the array of quaternions to be
!                  taken from the files and the two arrays mentioned above.
! Author: Thomas Greenspan
!         (modified version of part of preread_files.f90 by Paul Tregoning)
!********************************************************************************************************************************

subroutine input_readGTORB(access_type,irec,satics_t,apr_ScaleBias,GPSantoff,rvec,apr_prm &
     ,sciframe_quat,srf2trf_rotmat,srf2trf_deriv_rotmat,lu_offset, mcon_ocean &
     ,msc_tide_amp, mcon_tides)

  use gracefit_mod

  implicit none

   
  include 'input.h90'

  !********************  Variable declarations ********************

  character(*),     intent(in)  :: access_type                                ! "sequential" or "direct    "
  integer*4,      intent(inout) :: irec                                       ! record number of last record read from binary file (zero otherwise)
  double precision, intent(in)  :: satics_t(maxrvecprm+1,nsat_t)              ! A priori values for IC's (read from header)
  double precision, intent(in)  :: apr_ScaleBias(max_SB,nsat_t)               ! A priori values for bias and scale (read from header)
  double precision, intent(in)  :: GPSantoff(7,maxsat)                        ! Quaternion and offsets for antenna
  double precision, intent(out) :: rvec(nrvec_t,nsat_t,nepochs_t)             ! Pos, vel and partials for each satellite at every epoch
  double precision, intent(out) :: apr_prm(nparam_t)                          ! A priori values for parameters
  double precision, intent(out) :: sciframe_quat(4,nsat_t,nepochs_t)          ! Quaternions for each satellites at every epoch
  double precision, intent(out) :: srf2trf_rotmat(3,3,nepochs_t,nsat_t)       ! SRF to TRF rotation matrix for GRACE A and B
  double precision, intent(out) :: srf2trf_deriv_rotmat(3,3,nepochs_t,nsat_t) ! Differentiated SRF to TRF rotation matrix for A and B
  integer*4,        intent(in)  :: lu_offset                                  ! offset to indicate whether theoretical or truth files (affects unit numbers used)
  logical,          intent(out) :: mcon_ocean(nmascons_t)                     ! logical array to indicate whether an ocean or land mascon
  real(kind=8)    , intent(out) :: msc_tide_amp(max_mcon_tides,2,nmascons_t)  ! a priori mascon tidal amplitudes for each component/constituent/mascon
  integer*4        ,intent(in)  :: mcon_tides(4556,2)                         ! bit-mapped indicator of which tidal amplitudes to estimate per mascon


  integer*4        :: isat                    ! Counter for satellites
  integer*4        :: iepoch                  ! Counter for epochs
  integer*4        :: i,j,k                   ! Counter variables
  integer*4        :: mascon_num              ! Number of mascons
  integer*4        :: tmp_epoch               ! temporary epoch variable to then be converted to R*8
  double precision :: epoch(nsat_t,nepochs_t) ! Epochs recorded in GTORB file (only used to validate data)
  character(150)   :: line                    ! Line to be read 
  character(150)   :: message                 ! Message printed when done reading file
  integer*4        :: ioerr                   ! Standard error variable
  integer*4        :: irec_eoh                ! storing the input irec value for use with the second GTORB file
  integer*4        :: mascon_tidal_amp_num    ! value read off the bottom of the GTORB file
  integer*4        :: count_mcon_tides        ! counter used to transfer apriori mcon tidal amplitudes to apr_prm variable
  logical          :: bitmap

  ! PT131024: variables to segment and store the parials according to the parameters being estimated
  integer*4        :: pntr_rvec               ! pointer to which row we are up to in rvec, based on which parameters are estimated
  double precision, allocatable :: tmprvec(:) ! temporary storage of all partials of a single epoch, before transferring reqired values to rvec


  !****************************************************************

  ! define the size of the temporary partials storage vector
  ! PT140821: add another 24 to this for the once- and twice-per-rev along-track parameters
  allocate(tmprvec(36+18+18+24+24+nmascons_t*6 + nmsc_tid_constit_t*2*6))

  ! Header already read in beforehand
  ! Set IC a priori values based on data taken from header
  apr_prm = 0.d0
  do isat = 1, nsat_t
     do i = iICparm, iICparm+ngprm_t-1
        apr_prm(i+(isat-1)*nsatprm_t) = satics_t(i-iICparm+1+1,isat)
     enddo
     do i = iScale, iScale+nScale_t-1
        apr_prm(i+(isat-1)*nsatprm_t) = apr_ScaleBias(i-iScale+1,isat)
     enddo
     do i = iBias, iBias+nBias_t-1
        apr_prm(i+(isat-1)*nsatprm_t) = apr_ScaleBias(i-iBias+3+1,isat)
     enddo
     ! PT140821: added code for a priori values for once- and twice-per-rev along-track acceleration
     do i = ionepr, ionepr+n_onepr_t-1
        apr_prm(i+(isat-1)*nsatprm_t) = satics_t(i-ionepr+iICparm+ngprm_t+1,isat)
     enddo
     do i = itwopr, itwopr+ntwopr_t-1
        apr_prm(i+(isat-1)*nsatprm_t) = satics_t(i-itwopr+iICparm+ngprm_t+3,isat)
     enddo
     do i = iantoff, iantoff+nantoff_t-1
        apr_prm(i+(isat-1)*nsatprm_t) = GPSantoff(i-iantoff+1+4,isat)
     enddo
  enddo

  ! store the original pointer of the input record number (= the "END OF HEADER" record)
  irec_eoh = irec

  !*********************** READ GTORB FILES ***********************

  ! Read in positions, velocities and partials as well as quaternions for each satellite every epoch
  ioerr = 0

  do isat = 1, nsat_t
     irec = irec_eoh             ! this needs to be here to reset the record number for the second satellite
     do iepoch = 1,nepochs_t
        if(access_type == "sequential")then
           ! PT140821: no longer allow anything other than direct access files
           call status_update('FATAL','input_readGTORB',' ',' ',"Only direct access binary files are handled since 20 Aug 2014",0)
           !          read(LUGTORB(isat+lu_offset),*,iostat=ioerr,end=101)&
           !            epoch(isat,iepoch),(rvec(i,isat,iepoch),i=1,nrvecprm_t),(sciframe_quat(i,isat,iepoch),i=1,4),&
           !            (tmprvec(i),i = 1,(12+nmascons_t)*6)
        else
           irec=irec+1
           read(LUGTORB(isat+lu_offset),rec=irec,iostat=ioerr)&
                tmp_epoch,(rvec(i,isat,iepoch),i=1,nrvecprm_t),(sciframe_quat(i,isat,iepoch),i=1,4),&
                ! PT140821: replaced the hardwired value of "12" by norbprm_t
                ! PT140904: no, this is wrong. We need to read here as many values as have been written out in the GTORB file, even if we
                !           don't want to use them all. At present, this is 3 pos, 3 vel, 3 bias, 3 scale, 2 1/rev, 2 2/rev = 16.
                ! PT170410: we now also have roll/pitch/yaw so it is now 19
                (tmprvec(i),i = 1,(19+nmascons_t+nmsc_tid_constit_t*2)*6)
           epoch(isat,iepoch)=dble(tmp_epoch)        ! this is to account for GRACEORB writing out as integer but gracefit reading as R*8
        endif


        ! now, transfer into rvec the partials of the parameters to be estimated
        ! they are in tmprvec in order: 36x partials wrt ICs, 18x partials wrt acc bias, 18x partials wrt scls, 6x partials wrt each mascon, 6x2 partials wrt each tidal constituent amplitude

        if(ngprm_t == 3 ) then  ! GPS positions to be estimated but not velocity
           rvec(nrvecprm_t+1:nrvecprm_t+3 ,isat,iepoch) = tmprvec(1:3)      ! X wrt XYZ                     
           rvec(nrvecprm_t+4:nrvecprm_t+6,isat,iepoch) = tmprvec(7:9)       ! Y wrt XYZ   
           rvec(nrvecprm_t+7:nrvecprm_t+9,isat,iepoch) = tmprvec(13:15)     ! Z wrt XYZ 
           pntr_rvec = nrvecprm_t+9
        else if(ngprm_t == 6) then  ! GPS positions and velocities to be estimated
           rvec(nrvecprm_t+1:nrvecprm_t+36 ,isat,iepoch) = tmprvec(1:36)     ! all pos/vel wrt all pos/vel   
           pntr_rvec = nrvecprm_t+36
        endif

        ! do we need partials for acc biases?
        if(nBias_t == 3) then
           rvec(pntr_rvec+1:pntr_rvec+18,isat,iepoch) = tmprvec(37:54)
           pntr_rvec = pntr_rvec + 18
        endif

        ! do we need partials for acc scales?
        if(nScale_t == 3) then
           rvec(pntr_rvec+1:pntr_rvec+18,isat,iepoch) = tmprvec(55:72)
           pntr_rvec = pntr_rvec + 18
        endif

        ! PT140821: do we ned partials for once-per-rev (SINE and COSINE)?
        if(n_onepr_t == 2) then
           rvec(pntr_rvec+1:pntr_rvec+12,isat,iepoch) = tmprvec(73:84)
           pntr_rvec = pntr_rvec + 12
        endif

        ! PT140821: do we ned partials for twice-per-rev (SINE and COSINE)?
        if(ntwopr_t == 2) then
           rvec(pntr_rvec+1:pntr_rvec+12,isat,iepoch) = tmprvec(85:96)
           pntr_rvec = pntr_rvec + 12
        endif

        ! PT170410: do we need partials for roll/pitch/yaw ?
        if(nrpy_t == 3) then
           rvec(pntr_rvec+1:pntr_rvec+18,isat,iepoch) = tmprvec(97:114)
           pntr_rvec = pntr_rvec + 18
        endif


        ! now store the partials of mascon wrt pos (and maybe vel)
        if(nmascons_t > 0)then
           if(ngprm_t == 3)then
              do j=1,nmascons_t
                 rvec(pntr_rvec+1+(j-1)*3:pntr_rvec+3+(j-1)*3,isat,iepoch) = tmprvec(97+(j-1)*3:99+(j-1)*3)    ! XYZ wrt mascons (don't transfer vel partials wrt mascons)
                 pntr_rvec = pntr_rvec + 3
              enddo
           elseif(ngprm_t == 6)then
              ! PT140904: changed the index into tmprvec from "73+" to "97+" to account for the presence of the 1/rev and 2/rev partials now in the GTORB files
              ! PT170410: changed again from "97+" tp "115+" to account for roll/pitch/yaw
              rvec(pntr_rvec+1:pntr_rvec+nmascons_t*6,isat,iepoch) = tmprvec(115:115+nmascons_t*6 - 1)  ! both position and velocity mascon partials
              pntr_rvec = pntr_rvec + nmascons_t*6
           endif
        endif

        ! PT140617: now store the partials of mascon tide amplitudes wrt pos (and maybe vel)
        if(est_msc_tides > 0)then
           if(ngprm_t == 3)then
              do j=1,nmsc_tid_constit_t*2
                 ! PT150825: the indices into tmprvec were completely wrong here ... they didn't account for the 1/rev and 2/rev as per the mascon code above, nor that it should be nmascons_t*6
                 rvec(pntr_rvec+1+(j-1)*3:pntr_rvec+3+(j-1)*3,isat,iepoch) = &
                      tmprvec(115+nmascons_t*3+(j-1)*3:117+nmascons_t*3+(j-1)*3)    ! XYZ wrt mascons (don't transfer vel partials wrt mascons)
                 pntr_rvec = pntr_rvec + 3
              enddo
           elseif(ngprm_t == 6)then
              ! PT150825: changed the index here from "73+" to "97+" as well ... should have done that a year ago when I did it for the mascons above!!
              ! PT170410: updated for the roll/pitch/yaw (so now starts at 115, not 97)
              rvec(pntr_rvec+1:pntr_rvec+nmsc_tid_constit_t*2*6,isat,iepoch) = &
                   tmprvec(115+nmascons_t*6:117+nmascons_t*6+nmsc_tid_constit_t*2*6 - 1)  ! both position and velocity mascon partials
              ! DEBUG
              !  print*,'rvec(M2s:M2c,isat,iepoch)',rvec(pntr_rvec+1:pntr_rvec+2,isat,iepoch)
           endif
        endif
        !  print*,'after mascon tides pntr_rvec = ',pntr_rvec+nmsc_tid_constit_t*2*6

        call rotation_quat2mat_3d (sciframe_quat(:,isat,iepoch),srf2trf_rotmat(:,:,iepoch,isat))

        ! Assert there was no error reading file
        if(ioerr /= 0) call status_update('FATAL','GRACEFIT','input_readGTORB',' ','Error reading GTORB input file',ioerr)
        ! Validate recorded data once both satellites have been read in
        if(isat == nsat_t) call input_isValid(iepoch, epoch(1,iepoch), epoch(nsat_t,iepoch),'input_readGTORB','FATAL')
     enddo ! end of epoch loop

     write(message,'(a,i6,a,a)')"Read    ",iepoch-1,' epochs from GTORB_',SATNAM(isat)
     call status_update('STATUS','GRACEFIT','input_readGTORB',' ',message,0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! ******** Mascon  a priori values  ************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Read in mascon a priori values if necessary
     if(nmascons_t > 0)then
        if(access_type == "sequential")then
           read(LUGTORB(isat+lu_offset),'(a)',iostat=ioerr) line  ! this is the "END OF PARTIALS" line
           read(LUGTORB(isat+lu_offset),'(a)',iostat=ioerr) line  
        else
           irec=irec+2
           read(LUGTORB(isat+lu_offset),rec=irec,iostat=ioerr) line  
        endif
        if(line(1:17) /= 'NUMBER OF MASCONS') call status_update('FATAL','GRACEFIT','input_readGTORB',' ', &
             'expecting mascons at end of GTORB file',0)
        read(line(19:),*,iostat=ioerr)mascon_num
        if(mascon_num /= nmascons_t) then
           call status_update('WARNING','GRACEFIT','input_readGTORB',' ', &
                'Number of mascons in GTORB file not the same as indicated in header',0)
        endif
        write(message,'(a,i6,a)')'Reading ',nmascons_t,' a priori mascon values'
        call status_update('STATUS','GRACEFIT','input_readGTORB',' ',message,0)
        if(access_type == "sequential")then
           read(LUGTORB(isat+lu_offset),*,iostat=ioerr) apr_prm(imascons:imascons+nmascons_t-1)
        else
           irec = irec+1
           read(LUGTORB(isat+lu_offset),rec=irec,iostat=ioerr) apr_prm(imascons:imascons+nmascons_t-1)
           ! PT140527: I've changed graceorb to add 1000 m to the a priori value of any ocean mascon, so that I can distinguish between
           !           ocean and land mascons within gracefit. Subtract it off here, and make an array of ocean/land identifiers for mascons
           do i=1,nmascons_t
              if(apr_prm(imascons+i-1) > 900.d0)then   ! it is an ocean mascon
                 apr_prm(imascons+i-1) = apr_prm(imascons+i-1) - 1000.d0
                 mcon_ocean(i) = .true.
              else
                 mcon_ocean(i) = .false.
              endif
           enddo
        endif
     endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! ******** Mascon tidal amplitude a priori values  ***************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! PT140619: now read in the a priori mascon tidal amplitudes for sine/consine components of each tidal consitutent (most will be zero)
     if(est_msc_tides > 0) then
        irec = irec + 1
        read(LUGTORB(isat+lu_offset),rec=irec,iostat=ioerr) line  
        if(line(1:33) /= 'NUMBER OF MASCON TIDAL AMPLITUDES') call status_update('FATAL','GRACEFIT','input_readGTORB',' ', &
             'expecting mascon tidal amplitudes at end of GTORB file',0)
        read(line(35:),*,iostat=ioerr)mascon_tidal_amp_num
        ! Amplitudes are written out as 4556 mascons/record, with sine then cosine in alternating records. M2 O1 S2 K1 K2 constituents at present, making 10 records in total.
        count_mcon_tides = 0
        do j=1,max_mcon_tides
           do k=1,2
              irec = irec + 1
              read(LUGTORB(isat+lu_offset),rec=irec)(msc_tide_amp(j,k,i),i = 1,nmascons_t)
              ! DEBUG
              !  print*,'amplitudes for 4535 4536',j,k,msc_tide_amp(j,k,4535:4536)
           enddo
        enddo
        call status_update('STATUS','GRACEFIT','input_readGTORB',' ',' Read mascon tidal amplitudes',0)
     endif

     ! transfer to apr_prm the amplitudes that are being estimated
     if(est_msc_tides > 0) then
        do i=1,nmascons_t
           do j=1,max_mcon_tides
              if(bitmap(mcon_tides(i,1),j))then
                 do k=1,2
                    count_mcon_tides = count_mcon_tides + 1
                    apr_prm(imsctide+count_mcon_tides -1) = msc_tide_amp(j,k,i)
                 enddo
              endif
           enddo
        enddo
        write(message,'(a,i7,a)')' Transferred ',count_mcon_tides,' a priori mascon tidal amplitudes to apr_prm'
        call status_update('STATUS','GRACEFIT','input_readGTORB',' ',message,0)
     endif

  enddo ! end of satellite loop
  !****************************************************************

  ! Numerically differentiate the SRF to TRF rotation matrix elements
  do isat = 1, nsat_t
     do i = 1, 3
        do j = 1, 3
           call noise_robust_deriv(srf2trf_rotmat(i,j,:,isat),srf2trf_deriv_rotmat(i,j,:,isat) &
                ,epoch_interval,nepochs_t,5)
        enddo
     enddo
  enddo
  !****************************************************************

  ! Close files, they should no longer be necessary
  do isat = 1, nsat_t
     close(LUGTORB(isat+lu_offset))
  enddo

  ! PT131104: reset irec back to the last line of the header file (required so that gracesim can now read the real GTORB files)
  irec = irec_eoh

  return   ! return if successful

101 call status_update('FATAL','GRACEFIT','input_readGTORB',' ',&
       'Number of epochs in GTORB file do not match number indicated in header',0)

  return
end subroutine input_readGTORB

!********************************************************************************************************************************
!********************************************************************************************************************************
! input_readACC1B: Reads in all needed data from ACC1B files and stores it in the array acc that it is
!                  given. Validates input (concerning epochs) as it is recorded.
! Author: Thomas Greenspan
!         (modified version of part of preread_files.f90 by Paul Tregoning)
!
! PT180625: modified to use the subroutines in lib/acc_lib.f90 so that it can read either GRACE or GRACE FO L1B ascii obs
!********************************************************************************************************************************
subroutine input_readACC1B(thrusters_off)
  use gracefit_mod
  use accred_mod    ! arrays for the accelerometer observations
  use inmod_mod     ! provides model options for removing non-linear effects in accelerometer obs

  implicit none

   
  include 'input.h90'

  !********************  Variable declarations ********************
! passed variables
  logical, intent(in)   :: thrusters_off(nepochs_t,nsat_t)   ! indicates whether each obs is affected by thrusts

! local variables
  integer*4      :: isat          ! Counter for satellites
  integer*4      :: iepoch        ! Counter for epochs
  integer*4      :: i             ! Counter variable
  character(2)   :: dummy_char    ! Dummy variable used to store unwanted data from ACC1B
!  character(150) :: message       ! Message printed when done reading file
  integer*4      :: ioerr         ! Standard error variable
  integer*4      :: acc_span(2,2) ! start/end time (in gracesec) for accelerometer data
  integer*4      :: nvals_acc(2)  ! temp counter of the number of accelerometer values in a ACC1B file
  real(kind=8)   :: params(3)     ! parameters of the model fit to the accelerometer obs
  real(kind=8),allocatable :: tmp_thr_obs(:,:)

  integer*4      :: start_ep,end_ep ! start and end epochs for fitting linearising model to accelerometers
  !****************************************************************

  ! PT180625: read the header of the ACC1B file(s). This will read either GRACE or GRACE FO format
  do isat = 1, nsat_t
     call acc_read_hdr(LUAC(isat),calling_prog,acc_fnam(isat),mission,acc_span(isat,:) )
     nvals_acc(isat) = acc_span(isat,2) - acc_span(isat,1) + 1
  enddo

! dimension the accelerometer observation array, using whichever satellite has the most observations
  if(nvals_acc(2) > nvals_acc(1))then
    allocate(acc_obs(2,nvals_acc(2),5))
    ACCn = nvals_acc(2)
  else
    allocate(acc_obs(2,nvals_acc(1),5))
    ACCn = nvals_acc(1)
  endif

! now read in the ACC1B observations   
  do isat = 1, nsat_t
     call acc_read_data(LUAC(isat),calling_prog,acc_fnam(isat),mission,nvals_acc,3,acc_span(isat,:),acc_obs(isat,:,:) )
! convert to um/s^2
     acc_obs(isat,:,2:4) = acc_obs(isat,:,2:4)*m_um
  enddo

! make a temporary integer array of thrust on/off information (rather than the logical as stored in gracefit)
! PT180726: need to make this tmp_thr_obs array to have the same temporal separation as the ACC1B data (ie 1 second) 
!           I've hardwired this to account for a 5-sec to 1-sec difference ..... not sure how to do this generically !!!!
  allocate(tmp_thr_obs(ACCn,nsat_t))
  tmp_thr_obs = 0.d0
  do isat = 1,nsat_t
    do iepoch =1,nepochs_t
      if(thrusters_off(iepoch,isat))then
        tmp_thr_obs(iepoch:iepoch+4,isat) = 1.d0
!        print*,'thrust unaffected epoch',iepoch,isat,thrusters_off(iepoch,isat)
      else
!        print*,'is it or not?',iepoch,isat
      endif
    enddo
  enddo

! PT180702: now apply the relevant model to remove non-linear behaviour in the accelerometer obs
  do isat = 1,nsat_t
    do i=1,3
      if( nint(acc_obs_model(i)) == 0)then       ! we don't want to model anything - leave obs as they are
        write(message,'(a,i3,a,i3)')"Removing    no     model from ACC1B obs for axis: ",i," satellite ",isat

      else if (nint(acc_obs_model(i)) == 2)then  ! remove a quadratic and rate (but leave offset)
        start_ep = acc_span(isat,1)-nint(starting_epoch)+1
        end_ep   = acc_span(isat,2)-nint(starting_epoch)+1

        call acc_fit_quadratic(.false.,calling_prog,ACCn,acc_obs(isat,:,i+1),tmp_thr_obs(:,isat),start_ep,end_ep,params) 
        write(message,'(2(a,i3),a,3e13.4,a)')"Fitted   quadratic model to   ACC1B obs for axis: ",i," satellite ",isat &
             ," (",params,") for accel, rate, offset. nm/s^2"

     else
        write(message,'(a,a,a)')"Model `",acc_obs_model(i)," not coded. Please feel free to code it yourself :-) "
      endif

      call status_update('STATUS',calling_prog,"input_readACC1B",' ',message,0)
    enddo
  enddo

  ! Close files, they should no longer be necessary
  do isat = 1, nsat_t
     close(LUAC(isat))
  enddo

  return   ! return if successful

  end subroutine input_readACC1B

!********************************************************************************************************************************
!********************************************************************************************************************************
! input_readTHR1B: Reads in all needed data from THR1B files, stores it in the array thr and sets logical
!                  array thrusters_off accordingly (true if there are no thrusts affecting acceleration
!                  data at given iepoch and false otherwise). Validates input (concerning epochs) as it is
!                  recorded.Input is the array, thr, of thruster data and the logical array, thrusters_off.
! Author: Thomas Greenspan
!         (modified version of part of preread_files.f90 by Paul Tregoning and myself)
!********************************************************************************************************************************

subroutine input_readTHR1B(thrusters_off)
  use gracefit_mod

  implicit none

   
  include 'input.h90'

!********************  Variable declarations ********************

  logical,          intent(out) :: thrusters_off(nepochs_t,nsat_t) ! Data on whether thrusters are having any effect on accelerations

  integer*4      :: isat           ! Counter for satellites
  integer*4      :: ithrust        ! counter for thrusts
  integer*4      :: iepoch, kepoch ! Counters for epochs
  integer*4      :: i              ! Counter variable
  integer*4      :: dummy_ints(14) ! Dummy variable used to store unwanted data from THR1B
  character(2)   :: dummy_char     ! Dummy variable used to store unwanted data from THR1B
  character(150) :: message        ! Message printed when done reading file
  integer*4      :: ioerr          ! Standard error variable

  integer*4      :: n_thrusts(2)   ! number of thrusts in the THR1B file for each satellite
  integer*4      :: max_thrusts    ! the max number of thrusts of the two satellites
  integer*4, allocatable      :: thr_obs(:,:,:)   ! thrust observations for each satellite

  integer*4      :: thr_start,thr_end
!****************************************************************


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the headers of the two thruster files
  do isat = 1,nsat_t
    call thr_read_hdr(LUTH(isat),calling_prog,thr_fnam(isat),mission,n_thrusts(isat) )
  enddo

! allocate the thrust observation array.
  if (n_thrusts(1) >= n_thrusts(2))then
    max_thrusts = n_thrusts(1)
  else
    max_thrusts = n_thrusts(2)
  endif
  allocate(thr_obs(2,max_thrusts,8))
  thr_obs = 0

! read the thrust observations from the THR1B files
  do isat = 1,nsat_t
    call thr_read_data(LUTH(isat),calling_prog,thr_fnam(isat),mission,n_thrusts(isat),max_thrusts,thr_obs(isat,:,:) )
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set the flags to indicate the epochs affected by 


  thrusters_off = .true.
! Set thrusters_off according to data in the THR1B files
  do isat= 1, nsat_t
    do ithrust = 1,n_thrusts(isat)
      if(thr_obs(isat,ithrust,1) > nint(starting_epoch) .and.  &
        (thr_obs(isat,ithrust,1) - nint(starting_epoch))/epoch_interval < nepochs_t )then
        ! PT180725: we want to remove 70 seconds of data either side of the thrust epoch
!print*,'in if statement',thr_obs(isat,ithrust,1)-70 , nint(starting_epoch)
        if(thr_obs(isat,ithrust,1)-70 < nint(starting_epoch))then
          thr_start = 1
        else
          thr_start = nint( (thr_obs(isat,ithrust,1)-nint(starting_epoch)-70)/epoch_interval)
! PTRMcG 180927: prevent zero value if a thrust occurs in the first few seconds of a day
          if(thr_start == 0)thr_start = 1
       endif
        if(thr_obs(isat,ithrust,1)-nint(starting_epoch)+70 > nepochs_t*nint(epoch_interval) )then
          thr_end = nepochs_t
        else
          thr_end = nint( (thr_obs(isat,ithrust,1)-nint(starting_epoch)+70)/epoch_interval)
        endif
      
      ! now, change the flag to indicate that the epochs are affected by thrusts
!DEBUG
!print*,'isat,ithrust',isat,ithrust,thr_start,thr_end,thr_obs(isat,ithrust,1),nint(starting_epoch) &
!         ,(thr_obs(isat,ithrust,1) - nint(starting_epoch))/epoch_interval,nepochs_t

        thrusters_off(thr_start:thr_end,isat) = .false.
      endif
    enddo
  enddo


! Close files, they should no longer be necessary
  do isat = 1, nsat_t
     close(LUTH(isat))
  enddo

! deallocate variables
  deallocate(thr_obs)

  return
end subroutine input_readTHR1B
! *****************************************************************************************





!********************************************************************************************************************************
!********************************************************************************************************************************
! input_readGNV1B: Reads in all needed data from GNV1B files and stores it in the arrays gvec and uvec
!                  that it takes as input. Validates input (concerning epochs) as it is recorded.
!                  NOTE: The purpose of uvec is to hold data from the the GNV1B files that is
!                        currently not needed. It has been left as a comment to ease the process
!                        of adding it in if it becomes necessary in the future. (130717)
! Author: Thomas Greenspan
!
! PT180627: modified to use the lib subroutines to read GRACE and GRACE FO data
!********************************************************************************************************************************

subroutine input_readGNV1B(gvec,num_epochs_used,n_gnv)
  use gracefit_mod

  implicit none

   
  include 'input.h90'

!********************  Variable declarations ********************
  double precision, intent(out)      :: gvec(maxgprm,nsat_t,nepochs_t)     ! Pos and vel for each satellite at every epoch
  integer*4, intent(out)             :: num_epochs_used(nsat_t)            ! Number of GPS epochs that will be used in gracefit

! PT180628: variables for interacting with the lib routines for reading L1B data
  real(kind=8)          :: gnv_step(2)
  integer*4             :: gnv_span(2)        ! span (in seconds) from first to last GNV1B observation
  integer*4             :: n_gnv

! local variables
  integer*4        :: isat                    ! Counter for satellites
  integer*4        :: iepoch                  ! Counter for epochs
  character(256)   :: message                 ! Message printed when done reading file
  integer*4        :: ioerr                   ! Standard error variable
!****************************************************************

! PT180627: read the headers of the files and read the data 
  do isat = 1,nsat_t
    call gnv_read_hdr(LUGN(isat),calling_prog ,gn_fnam(isat),mission,num_epochs_used(isat),gnv_step(isat),gnv_span(isat) )
  enddo

! allocate the array of epochs for gvec (unused in this subroutine but required when calling gnv_read_data)
  if(gnv_span(2) > gnv_span(1))n_gnv = int(gnv_span(2) / gnv_step(2))
  if(gnv_span(1) >= gnv_span(2))n_gnv = int(gnv_span(1) / gnv_step(1))
  allocate( gvec_epochs(n_gnv,2))

! read the GNV1B data
  do isat = 1,2
    call gnv_read_data(LUGN(isat),calling_prog ,gn_fnam(isat),mission,maxgprm,nsat_t,nepochs_t,isat,starting_epoch &
                       ,epoch_interval,num_epochs_used(isat),gnv_step(isat),gvec_epochs(:,isat),n_gnv,gvec)
  enddo

! Close files, they should no longer be necessary
  do isat = 1, nsat_t
     close(LUGN(isat))
  enddo

  return
end subroutine input_readGNV1B

!********************************************************************************************************************************
!********************************************************************************************************************************
! input_readKBR1B: Reads in all needed data from KBR1B file and stores it in the arrays kbrange and
!                  kblt and kbant that it takes as input. Validates input (concerning epochs) as it
!                  is recorded.
!                  NOTE: The purpose of kbion_cor and kbsnr is to hold data from the the KBR1B file 
!                        that is currently not needed. It has been left as a comment to ease the
!                        process of adding it in if it becomes necessary in the future. (130717)
! Author: Thomas Greenspan
!
! MODS
! PT201014: added use_obs to the argument list and filled out the values within this subroutine
!********************************************************************************************************************************

subroutine input_readKBR1B(kbrange,kblt,kbant)
  use gracefit_mod

  implicit none

   
  include 'input.h90'

  !********************  Variable declarations ********************

  double precision, intent(out) :: kbrange(3,nepochs_t)      ! Range value (biased), rate, and acceleration
  double precision, intent(out) :: kblt(3,nepochs_t)         ! Light time corrections for range, rate and acceleration
  double precision, intent(out) :: kbant(3,nepochs_t)        ! Antenna phase center corrections for range, rate and acceleration
 
  integer*4        :: iepoch           ! Counter for epochs
  integer*4        :: kb_iepoch        ! Index as to where the epoch read in from the KBR1B file should go
  integer*4        :: i                ! Counter variable
  double precision :: epoch            ! Epochs recorded in KBR1B file (only used to validate data)
  integer*4        :: num_epochs       ! Number of epochs read in from file for each satellite
  integer*4        :: missing_interval ! Number of seconds missing (at each occasion)
  double precision :: dummy_var        ! Dummy variable used to store unwanted data from KBR1B
  character(150)   :: message          ! Message printed when done reading file
  integer*4        :: ioerr            ! Standard error variable

  integer*4        :: n_kbr

! PT180916: variables to enable reading of the error flags in the KBR1B files
  integer*4 :: dummy_values(4)
  integer*4 :: error_flags
  !****************************************************************

! Read in header (make necessary checks)
  call kbr_read_hdr(LUKB,calling_prog,kb_fnam,mission,n_kbr)

!*********************** READ KBR1B FILES ***********************
! PT181119: move all the reading of the KBR data into subroutine kbr_read_data
  call kbr_read_data(LUKB,calling_prog,kb_fnam,mission,nepochs_t,starting_epoch,epoch_interval,kbrange,kblt,kbant,num_epochs)

! PT201014: loop through epochs and set use_obs to true for KBR,KBRR,KBRA if kbr obs exist
! SA201028: Seems that this doesn't work!!! 
!  do iepoch=1,nepochs_t
!    if(kbrange(iepoch,1) /= 0.d0)use_obs(:,iepoch) = .true.
!  enddo

! wrap it up with some output information
  nKBobs_t = num_epochs
  write(message,'(a,i6,a,a)') 'read ',num_epochs,' epochs from KBR1B'
  call status_update('STATUS',calling_prog,'input_readKBR1B',' ',message,0)
  ! Verify there are enough observations
  write(message,'(a,i6,a)')'Not enough Kband observations (minimum: ',minKBobs,')'
  if(nKBobs_t < minKBobs) call status_update('FATAL',calling_prog,'input_readKBR1B',' ',message,0)


  return
end subroutine input_readKBR1B

!********************************************************************************************************************************
!********************************************************************************************************************************
! input_isValid: Checks that the epochs read in from a file are correct given the index, iepoch, they are read
!                in at, the epoch recorded in grace seconds for GRACE A and B, the subroutine the files are
!                read in from and the error type if there is one.
! Author: Thomas Greenspan
!********************************************************************************************************************************

subroutine input_isValid(iepoch, epoch_A, epoch_B, subroutine_name, error_type)
  use gracefit_mod

  implicit none

   

  !********************  Variable declarations ********************

  integer*4, intent(in)        :: iepoch          ! Index from which the data is being read
  double precision, intent(in) :: epoch_A         ! Epoch read in from GRACE A
  double precision, intent(in) :: epoch_B         ! Epoch read in from GRACE B
  character(*), intent(in)     :: subroutine_name ! Subroutine name
  character(*), intent(in)     :: error_type      ! Type of error if any occurs
  !****************************************************************

  ! Assert epochs are the same for GRACE A and B
  if(epoch_A /= epoch_B) call status_update(error_type,'GRACEFIT',subroutine_name,' ',&
       'Error reading in files: epochs not the same for GRACE A and B',0)
  ! Assert epochs in grace seconds correspond correctly to their index, iepoch
  if(epoch_A /= starting_epoch+epoch_interval*(iepoch-1)) then
     print*,'starting_epoch,epoch_interval,iepoch-1',starting_epoch,epoch_interval,iepoch-1,epoch_A,epoch_B &
          , starting_epoch+epoch_interval*(iepoch-1)
     call status_update(error_type,'GRACEFIT',subroutine_name,' ',&
          'Error reading in files: epochs do not correspond correctly to iepoch',0)
  endif
  !****************************************************************

  return
end subroutine input_isValid

!********************************************************************************************************************************





subroutine input_readGTORBh5_v1(satics_t,apr_ScaleBias,GPSantoff,rvec, part, apr_prm &
   ,sciframe_quat,srf2trf_rotmat,srf2trf_deriv_rotmat,lu_offset, mcon_ocean &
   ,msc_tide_amp, mcon_tides)

use gracefit_mod
implicit none

 

!********************  Variable declarations ********************

double precision, intent(in)  :: satics_t(maxrvecprm+1,nsat_t)              ! A priori values for IC's (read from header)
double precision, intent(in)  :: apr_ScaleBias(max_SB,nsat_t)               ! A priori values for bias and scale (read from header)
double precision, intent(in)  :: GPSantoff(7,maxsat)                        ! Quaternion and offsets for antenna
double precision, intent(out) :: rvec(6 ,nsat_t,nepochs_t)                  ! Pos, vel and partials for each satellite at every epoch
double precision, intent(out) :: part(nobs_t,nparam_t,nepochs_t)           ! Pos, vel and partials for each satellite at every epoch
double precision, intent(out) :: apr_prm(nparam_t)                          ! A priori values for parameters
double precision, intent(out) :: sciframe_quat(4,nsat_t,nepochs_t)          ! Quaternions for each satellites at every epoch
double precision, intent(out) :: srf2trf_rotmat(3,3,nepochs_t,nsat_t)       ! SRF to TRF rotation matrix for GRACE A and B
double precision, intent(out) :: srf2trf_deriv_rotmat(3,3,nepochs_t,nsat_t) ! Differentiated SRF to TRF rotation matrix for A and B
integer*4,        intent(in)  :: lu_offset                                  ! offset to indicate whether theoretical or truth files (affects unit numbers used)
logical,          intent(out) :: mcon_ocean(nmascons_t)                     ! logical array to indicate whether an ocean or land mascon
real(kind=8)    , intent(out) :: msc_tide_amp(max_mcon_tides,2,nmascons_t)  ! a priori mascon tidal amplitudes for each component/constituent/mascon
integer*4        ,intent(in)  :: mcon_tides(4556,2)                         ! bit-mapped indicator of which tidal amplitudes to estimate per mascon


integer*4        :: isat                    ! Counter for satellites
integer*4        :: iepoch                  ! Counter for epochs
integer*4        :: i,j,k                   ! Counter variables
integer*4        :: mascon_num              ! Number of mascons
integer*4        :: tmp_epoch               ! temporary epoch variable to then be converted to R*8
double precision :: epoch(nsat_t,nepochs_t) ! Epochs recorded in GTORB file (only used to validate data)
character(150)   :: line                    ! Line to be read 
character(150)   :: message                 ! Message printed when done reading file
integer*4        :: ioerr                   ! Standard error variable
integer*4        :: irec_eoh                ! storing the input irec value for use with the second GTORB file
integer*4        :: mascon_tidal_amp_num    ! value read off the bottom of the GTORB file
integer*4        :: count_mcon_tides        ! counter used to transfer apriori mcon tidal amplitudes to apr_prm variable
logical          :: bitmap

! PT131024: variables to segment and store the parials according to the parameters being estimated
integer*4        :: pntr_rvec               ! pointer to which row we are up to in rvec, based on which parameters are estimated
! double precision, allocatable :: tmprvec(:,:) ! temporary storage of all partials of a single epoch, before transferring reqired values to rvec

integer :: local_offset

local_offset = lu_offset


! Header already read in beforehand
! Set IC a priori values based on data taken from header
! allocate(tmprvec(36+18+18+24+24+nmascons_t*6 + nmsc_tid_constit_t*2*6, nepochs_t))
apr_prm = 0.d0
do isat = 1, nsat_t
   do i = iICparm, iICparm+ngprm_t-1
      apr_prm(i+(isat-1)*nsatprm_t) = satics_t(i-iICparm+1+1,isat)
   enddo
   do i = iScale, iScale+nScale_t-1
      apr_prm(i+(isat-1)*nsatprm_t) = apr_ScaleBias(i-iScale+1,isat)
   enddo
   do i = iBias, iBias+nBias_t-1
      apr_prm(i+(isat-1)*nsatprm_t) = apr_ScaleBias(i-iBias+3+1,isat)
   enddo
   ! PT140821: added code for a priori values for once- and twice-per-rev along-track acceleration
   do i = ionepr, ionepr+n_onepr_t-1
      apr_prm(i+(isat-1)*nsatprm_t) = satics_t(i-ionepr+iICparm+ngprm_t+1,isat)
   enddo
   do i = itwopr, itwopr+ntwopr_t-1
      apr_prm(i+(isat-1)*nsatprm_t) = satics_t(i-itwopr+iICparm+ngprm_t+3,isat)
   enddo
   do i = iantoff, iantoff+nantoff_t-1
      apr_prm(i+(isat-1)*nsatprm_t) = GPSantoff(i-iantoff+1+4,isat)
   enddo
enddo

! store the original pointer of the input record number (= the "END OF HEADER" record)

!*********************** READ GTORB FILES ***********************

! Read in positions, velocities and partials as well as quaternions for each satellite every epoch
ioerr = 0
do isat = 1, nsat_t
      call status_update('STATUS','GRACEFIT','input_readGTORBh5_v1',' ','reading orbit coords and partials',0)
      call orb_read(orb_files(isat+local_offset), time =epoch(isat,:), coord = rvec(:,isat,:),  &
        sciframe_quat = sciframe_quat(:,isat,:))
      call orb_read_in_part(orb_files(isat+local_offset), (/6*(isat-1),12*(isat-1),0/), (/6,12,nepochs_t/), partials = part)

   if (nmascons_t > 0 .and. lu_offset == 0) then
      call status_update('STATUS','GRACEFIT','input_readGTORBh5_v1',' ','reading mascon partials',0)
      call orb_read_in_part(orb_files(isat+local_offset), (/6*(isat-1),24,0/), (/6,nmascons_t,nepochs_t/), msc_part = part)


   endif
   if (nmsc_tid_constit_t > 0) then
      call status_update('WARNING','GRACEFIT','input_readGTORBh5_v1',' ','reading tidal amplitude partials not working yet',0)
 !     call orb_read(orb_files(isat+local_offset), &
 !          msc_part = tmprvec(114+nmascons_t*6+1:114+nmascons_t*6+nmsc_tid_constit_t*12,:))
   endif

   ! now, transfer into rvec the partials of the parameters to be estimated
   ! they are in tmprvec in order: 36x partials wrt ICs, 18x partials wrt acc bias, 18x partials wrt scls, 6x partials wrt each mascon, 6x2 partials wrt each tidal constituent amplitude
   ! print *, nparam_t/2, size(part)
   ! SA Change, we assume we need to solve for pos, vel, baises, scale
   ! if(ngprm_t == 6) then  ! GPS positions and velocities to be estimated
   !    do i = 1 , 12
   !       part( (isat-1)*6+1:(isat)*6, i+12*(isat-1), : ) = tmprvec((i-1)*6+1:i*6 ,:)     ! all pos/vel wrt all pos/vel   
   !    enddo
   !    !pntr_rvec = nrvecprm_t+36
   ! endif
   ! ! do we need partials for acc biases?
   ! if(nBias_t == 3) then
   !    rvec(pntr_rvec+1:pntr_rvec+18,isat,:) = tmprvec(37:54,:)
   !    pntr_rvec = pntr_rvec + 18
   ! endif

   ! ! do we need partials for acc scales?
   ! if(nScale_t == 3) then
   !    rvec(pntr_rvec+1:pntr_rvec+18,isat,:) = tmprvec(55:72,:)
   !    pntr_rvec = pntr_rvec + 18
   ! endif

   ! ! PT140821: do we ned partials for once-per-rev (SINE and COSINE)?
   ! if(n_onepr_t == 2) then
   !    rvec(pntr_rvec+1:pntr_rvec+12,isat,:) = tmprvec(73:84,:)
   !    pntr_rvec = pntr_rvec + 12
   ! endif

   ! ! PT140821: do we ned partials for twice-per-rev (SINE and COSINE)?
   ! if(ntwopr_t == 2) then
   !    rvec(pntr_rvec+1:pntr_rvec+12,isat,:) = tmprvec(85:96,:)
   !    pntr_rvec = pntr_rvec + 12
   ! endif

   ! ! PT170410: do we need partials for roll/pitch/yaw ?
   ! if(nrpy_t == 3) then
   !    rvec(pntr_rvec+1:pntr_rvec+18,isat,:) = tmprvec(97:114,:)
   !    pntr_rvec = pntr_rvec + 18
   ! ! endif
   ! print *, '****', nparam_t

   ! ! now store the partials of mascon wrt pos (and maybe vel)
   ! if(nmascons_t > 0 .and. lu_offset == 0)then
   !   call status_update('STATUS','GRACEFIT','input_readGTORB',' ','now store the partials of mascon wrt pos (and maybe vel)',0)
   !    if(ngprm_t == 3)then
   !       do j=1,nmascons_t
   !          rvec(pntr_rvec+1+(j-1)*3:pntr_rvec+3+(j-1)*3,isat,:) = tmprvec(97+(j-1)*3:99+(j-1)*3,:)    ! XYZ wrt mascons (don't transfer vel partials wrt mascons)
   !          pntr_rvec = pntr_rvec + 3
   !       enddo
   !    elseif(ngprm_t == 6)then
   !       ! PT140904: changed the index into tmprvec from "73+" to "97+" to account for the presence of the 1/rev and 2/rev partials now in the GTORB files
   !       ! PT170410: changed again from "97+" tp "115+" to account for roll/pitch/yaw
   !       ! rvec(pntr_rvec+1:pntr_rvec+nmascons_t*6,isat,:) = tmprvec(115:115+nmascons_t*6 - 1,:)  ! both position and velocity mascon partials
   !       print *, '****', nparam_t
   !       ! do j=1,nmascons_t
   !       !    part( (isat-1)*6+1:(isat)*6, 12*2 + j   , : )  = tmprvec(115+(j-1)*6:114+j*6,:) 
   !       ! enddo
   !       !pntr_rvec = pntr_rvec + nmascons_t*6
   !    endif
   ! endif

   ! call status_update('STATUS','GRACEFIT','input_readGTORB',' ','now tidal amplitudes',0)
   ! ! PT140617: now store the partials of mascon tide amplitudes wrt pos (and maybe vel)
   ! if(est_msc_tides > 0)then
   !    if(ngprm_t == 3)then
   !       do j=1,nmsc_tid_constit_t*2
   !          ! PT150825: the indices into tmprvec were completely wrong here ... they didn't account for the 1/rev and 2/rev as per the mascon code above, nor that it should be nmascons_t*6
   !          rvec(pntr_rvec+1+(j-1)*3:pntr_rvec+3+(j-1)*3,isat,:) = &
   !               tmprvec(115+nmascons_t*3+(j-1)*3:117+nmascons_t*3+(j-1)*3,:)    ! XYZ wrt mascons (don't transfer vel partials wrt mascons)
   !          pntr_rvec = pntr_rvec + 3
   !       enddo
   !    elseif(ngprm_t == 6)then
   !       ! PT150825: changed the index here from "73+" to "97+" as well ... should have done that a year ago when I did it for the mascons above!!
   !       ! PT170410: updated for the roll/pitch/yaw (so now starts at 115, not 97)
   !       rvec(pntr_rvec+1:pntr_rvec+nmsc_tid_constit_t*2*6,isat,:) = &
   !            tmprvec(115+nmascons_t*6:117+nmascons_t*6+nmsc_tid_constit_t*2*6 - 1,:)  ! both position and velocity mascon partials
   !       ! DEBUG
   !       !  print*,'rvec(M2s:M2c,isat,iepoch)',rvec(pntr_rvec+1:pntr_rvec+2,isat,iepoch)
   !    endif
   ! endif
   !  print*,'after mascon tides pntr_rvec = ',pntr_rvec+nmsc_tid_constit_t*2*6

!    write(*,*) epoch(1,1:10), epoch(2,1:10)
!    write(*,*) epoch(1,nepochs_t-9:nepochs_t), epoch(2,nepochs_t-9:nepochs_t)
   call status_update('STATUS','GRACEFIT','input_readGTORBh5_v1',' ','now quaternians',0)
   do iepoch = 1 , nepochs_t
      call rotation_quat2mat_3d (sciframe_quat(:,isat,iepoch),srf2trf_rotmat(:,:,iepoch,isat))
      ! Validate recorded data once both satellites have been read in
      if(isat == nsat_t) call input_isValid(iepoch, epoch(1,iepoch), epoch(nsat_t,iepoch),'input_readGTORBh5_v1','FATAL')
   enddo


!!! Apriori transfer to be done here! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! ******** Mascon  a priori values  ************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Read in mascon a priori values if necessary
   if(nmascons_t > 0)then
     call status_update('STATUS','GRACEFIT','input_readGTORBh5_v1',' ','now msc_apriori values',0)
      call orb_read(orb_files(isat+local_offset), msc_apr =  apr_prm(imascons:imascons+nmascons_t-1) )
      do i=1,nmascons_t
         if(apr_prm(imascons+i-1) > 900.d0)then   ! it is an ocean mascon
            apr_prm(imascons+i-1) = apr_prm(imascons+i-1) - 1000.d0
            mcon_ocean(i) = .true.
         else
            mcon_ocean(i) = .false.
         endif
      enddo
   endif

   if(est_msc_tides > 0) then


      ! call orb_read(orb_files(isat+local_offset), tmsc_apr =  msc_tide_amp )

      ! call status_update('STATUS','GRACEFIT','input_readGTORBh5_v1',' ',' Read mascon tidal amplitudes',0)

      ! do i=1,nmascons_t
      !    do j=1,max_mcon_tides
      !       if(bitmap(mcon_tides(i,1),j))then
      !          do k=1,2
      !             count_mcon_tides = count_mcon_tides + 1
      !             apr_prm(imsctide+count_mcon_tides -1) = msc_tide_amp(j,k,i)
      !          enddo
      !       endif
      !    enddo
      ! enddo
      ! write(message,'(a,i7,a)')' Transferred ',count_mcon_tides,' a priori mascon tidal amplitudes to apr_prm'
      ! call status_update('STATUS','GRACEFIT','input_readGTORBh5_v1',' ',message,0)
      write(message,'(a)')' Tidal mascons not implemented yet'
      call status_update('WARNING','GRACEFIT','input_readGTORBh5_v1',' ',message,0)
   endif
end do

!Numerically differentiate the SRF to TRF rotation matrix elements
call status_update('STATUS','GRACEFIT','input_readGTORBh5_v1',' ','differentiate SRF to TRF matrices',0)
do isat = 1, nsat_t
   do i = 1, 3
      do j = 1, 3
         call noise_robust_deriv(srf2trf_rotmat(i,j,:,isat),srf2trf_deriv_rotmat(i,j,:,isat) &
              ,epoch_interval,nepochs_t,5)
      enddo
   enddo
enddo

! deallocate(tmprvec)
end subroutine input_readGTORBh5_v1





subroutine input_readGTORBh5_Amat(satics_t,apr_ScaleBias,GPSantoff,rvec, part, apr_prm &
   ,sciframe_quat,srf2trf_rotmat,srf2trf_deriv_rotmat,lu_offset, mcon_ocean &
   ,msc_tide_amp, mcon_tides)

use gracefit_mod
implicit none

 
!********************  Variable declarations ********************

double precision, intent(in)  :: satics_t(maxrvecprm+1,nsat_t)              ! A priori values for IC's (read from header)
double precision, intent(in)  :: apr_ScaleBias(max_SB,nsat_t)               ! A priori values for bias and scale (read from header)
double precision, intent(in)  :: GPSantoff(7,maxsat)                        ! Quaternion and offsets for antenna
double precision, intent(out) :: rvec(6 ,nsat_t,nepochs_t)                  ! Pos, vel and partials for each satellite at every epoch
double precision, intent(out) :: part(nobs_t*nepochs_t,nparam_t)            ! Pos, vel and partials for each satellite at every epoch
double precision, intent(out) :: apr_prm(nparam_t)                          ! A priori values for parameters
double precision, intent(out) :: sciframe_quat(4,nsat_t,nepochs_t)          ! Quaternions for each satellites at every epoch
double precision, intent(out) :: srf2trf_rotmat(3,3,nepochs_t,nsat_t)       ! SRF to TRF rotation matrix for GRACE A and B
double precision, intent(out) :: srf2trf_deriv_rotmat(3,3,nepochs_t,nsat_t) ! Differentiated SRF to TRF rotation matrix for A and B
integer*4,        intent(in)  :: lu_offset                                  ! offset to indicate whether theoretical or truth files (affects unit numbers used)
logical,          intent(out) :: mcon_ocean(nmascons_t)                     ! logical array to indicate whether an ocean or land mascon
real(kind=8)    , intent(out) :: msc_tide_amp(max_mcon_tides,2,nmascons_t)  ! a priori mascon tidal amplitudes for each component/constituent/mascon
integer*4        ,intent(in)  :: mcon_tides(4556,2)                         ! bit-mapped indicator of which tidal amplitudes to estimate per mascon


integer*4        :: isat                    ! Counter for satellites
integer*4        :: iepoch                  ! Counter for epochs
integer*4        :: i,j,k                   ! Counter variables
integer*4        :: mascon_num              ! Number of mascons
integer*4        :: tmp_epoch               ! temporary epoch variable to then be converted to R*8
double precision :: epoch(nsat_t,nepochs_t) ! Epochs recorded in GTORB file (only used to validate data)
character(150)   :: line                    ! Line to be read 
character(150)   :: message                 ! Message printed when done reading file
integer*4        :: ioerr                   ! Standard error variable
integer*4        :: irec_eoh                ! storing the input irec value for use with the second GTORB file
integer*4        :: mascon_tidal_amp_num    ! value read off the bottom of the GTORB file
integer*4        :: count_mcon_tides        ! counter used to transfer apriori mcon tidal amplitudes to apr_prm variable
logical          :: bitmap

! PT131024: variables to segment and store the parials according to the parameters being estimated
integer*4        :: pntr_rvec               ! pointer to which row we are up to in rvec, based on which parameters are estimated
! double precision, allocatable :: tmprvec(:,:) ! temporary storage of all partials of a single epoch, before transferring reqired values to rvec

integer :: local_offset

local_offset = lu_offset


! Header already read in beforehand
! Set IC a priori values based on data taken from header
! allocate(tmprvec(36+18+18+24+24+nmascons_t*6 + nmsc_tid_constit_t*2*6, nepochs_t))
apr_prm = 0.d0
do isat = 1, nsat_t
   do i = iICparm, iICparm+ngprm_t-1
      apr_prm(i+(isat-1)*nsatprm_t) = satics_t(i-iICparm+1+1,isat)
   enddo
   do i = iScale, iScale+nScale_t-1
      apr_prm(i+(isat-1)*nsatprm_t) = apr_ScaleBias(i-iScale+1,isat)
   enddo
   do i = iBias, iBias+nBias_t-1
      apr_prm(i+(isat-1)*nsatprm_t) = apr_ScaleBias(i-iBias+3+1,isat)
   enddo
   ! PT140821: added code for a priori values for once- and twice-per-rev along-track acceleration
   do i = ionepr, ionepr+n_onepr_t-1
      apr_prm(i+(isat-1)*nsatprm_t) = satics_t(i-ionepr+iICparm+ngprm_t+1,isat)
   enddo
   do i = itwopr, itwopr+ntwopr_t-1
      apr_prm(i+(isat-1)*nsatprm_t) = satics_t(i-itwopr+iICparm+ngprm_t+3,isat)
   enddo
   do i = iantoff, iantoff+nantoff_t-1
      apr_prm(i+(isat-1)*nsatprm_t) = GPSantoff(i-iantoff+1+4,isat)
   enddo
enddo

! store the original pointer of the input record number (= the "END OF HEADER" record)

!*********************** READ GTORB FILES ***********************

! Read in positions, velocities and partials as well as quaternions for each satellite every epoch
ioerr = 0
do isat = 1, nsat_t
      call status_update('STATUS','GRACEFIT','input_readGTORBh5_Amat',' ','reading orbit coords and partials',0)
      call orb_read(orb_files(isat+local_offset), time =epoch(isat,:), coord = rvec(:,isat,:),  &
        sciframe_quat = sciframe_quat(:,isat,:))
      !  call orb_read_in_part(orb_files(isat+local_offset), (/6*(isat-1),12*(isat-1),0/), (/6,12,nepochs_t/), partials = part)
      call orb_read_in_aMat(orb_files(isat+local_offset), (/6*(isat-1),12*(isat-1)/), (/1, 1/),  (/nobs_t, 1/), (/6, 12/), nepochs_t,  part, 1)

   if (nmascons_t > 0 .and. lu_offset == 0) then
      call status_update('STATUS','GRACEFIT','input_readGTORBh5_Amat',' ','reading mascon partials',0)
      !call orb_read_in_part(orb_files(isat+local_offset), (/6*(isat-1),24,0/), (/6,nmascons_t,nepochs_t/), msc_part = part)
      call orb_read_in_aMat(orb_files(isat+local_offset), (/6*(isat-1),24/), (/1, 1/),  (/nobs_t, 1/), (/6, nmascons_t/), nepochs_t,  part, 2)


   endif
   if (nmsc_tid_constit_t > 0) then
      call status_update('WARNING','GRACEFIT','input_readGTORBh5_Amat',' ','reading tidal amplitude partials not working yet',0)
   endif

   call status_update('STATUS','GRACEFIT','input_readGTORBh5_Amat',' ','now quaternians',0)
   do iepoch = 1 , nepochs_t
      call rotation_quat2mat_3d (sciframe_quat(:,isat,iepoch),srf2trf_rotmat(:,:,iepoch,isat))
      ! Validate recorded data once both satellites have been read in

      if(isat == nsat_t) call input_isValid(iepoch, epoch(1,iepoch), epoch(nsat_t,iepoch),'input_readGTORBh5_Amat','FATAL')
   enddo


!!! Apriori transfer to be done here! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! ******** Mascon  a priori values  ************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Read in mascon a priori values if necessary
   if(nmascons_t > 0)then
     call status_update('STATUS','GRACEFIT','input_readGTORBh5_Amat',' ','now msc_apriori values',0)
      call orb_read(orb_files(isat+local_offset), msc_apr =  apr_prm(imascons:imascons+nmascons_t-1) )
      do i=1,nmascons_t
         if(apr_prm(imascons+i-1) > 900.d0)then   ! it is an ocean mascon
            apr_prm(imascons+i-1) = apr_prm(imascons+i-1) - 1000.d0
            mcon_ocean(i) = .true.
         else
            mcon_ocean(i) = .false.
         endif
      enddo
   endif

   if(est_msc_tides > 0) then

      write(message,'(a)')' Tidal mascons not implemented yet'
      call status_update('WARNING','GRACEFIT','input_readGTORBh5_Amat',' ',message,0)
   endif
end do

!Numerically differentiate the SRF to TRF rotation matrix elements
call status_update('STATUS','GRACEFIT','input_readGTORBh5_Amat',' ','differentiate SRF to TRF matrices',0)
!call status_update('WARNING','GRACEFIT','input_readGTORBh5_Amat',' ','differentiate SRF to TRF matrices disabled!!!!',0)
do isat = 1, nsat_t
   do i = 1, 3
      do j = 1, 3
         call noise_robust_deriv(srf2trf_rotmat(i,j,:,isat),srf2trf_deriv_rotmat(i,j,:,isat) &
              ,epoch_interval,nepochs_t,5)
      enddo
   enddo
enddo
end subroutine  input_readGTORBh5_Amat
