!********************************************************************************************************************************
!  File: header_lib.f90
!
!  Purpose: Set of subroutines used to read headers
!
!  Author: Thomas Greenspan
!
!  API:
!       header_readGTORB     : Reads the GTORB headers
!       header_readbinGTORB  : Reads the binary GTORB headers
!       header_readkbr1b     : read the kbr1b header information?
!
!  August 2, 2013
!
!********************************************************************************************************************************
!********************************************************************************************************************************
! header_readbinGTORB: Reads in data from the header of the binary GTORB file(s), checks to make sure files
!                   are consistent with each other and with the command file and sets values as
!                   needed
! Author: P. Tregoning
!         4 November 2013
!********************************************************************************************************************************

subroutine header_readh5GTORB(irec,satics_t, apr_ScaleBias, GPSantoff, mcon_tide1, combined_mascon_file &
     , msc_hdr_code,ocean_mascon_file,msc_ocn_hdr_code,total_ocean_prim, local_offset, prmnam,prmnam_size)

  !use gracefit_mod
  use gracefit_mod
  implicit none
   
  include 'input.h90' ! not sure we need it. 

  integer, intent(out) :: irec(2)                                 !< Number of record in GTORB
  double precision, intent(out) :: satics_t(maxrvecprm+1,nsat_t)  !< Apriori values for IC's
  double precision, intent(out) :: apr_ScaleBias(max_SB,nsat_t)   !< Apriori values for bias and scale
  double precision, intent(out) :: GPSantoff(7,nsat_t)            !< Quaternion and offsets for antenna
  character (len=*),     intent(out) :: combined_mascon_file,ocean_mascon_file  ! name of temporal and ocean mascon files used in GRACEORB
  character (len=*),      intent(out) :: msc_hdr_code,msc_ocn_hdr_code           ! mascon file codes of temporal and ocean mascon files used in GRACEORB
  integer  ,      intent(out) :: total_ocean_prim
  integer,        intent(inout) :: mcon_tide1(4000,2)              ! bit-mapped info about which tidal mascons to estimate
  integer , intent(in):: local_offset
  integer*4,intent(in):: prmnam_size
  character*30,intent(inout):: prmnam(prmnam_size)

  ! Data read in from header
  character(80)    :: gt_agency(nsat_t),gt_runby(nsat_t),gt_rundate(nsat_t) 
  character(80)    :: gt_datarecref(nsat_t),gt_models(nsat_t),gt_swver(nsat_t)
  character(8)     :: gt_filefmt(nsat_t),gt_satnam(nsat_t),ics_name(maxrvecprm+1,nsat_t)
  character(8)     :: ScaleBias_name(max_SB,nsat_t),GPSantoff_name(7,nsat_t),gt_datafmt(30,nsat_t)
  integer*4        :: gt_hdrrecs(nsat_t),gt_datarecs(nsat_t),gt_nmcs(nsat_t)
  integer*4        :: gt_stime(nsat_t),gt_etime(nsat_t)
  double precision :: gt_datarecint(nsat_t),satics_i(maxrvecprm+1,nsat_t)
  character (80):: msc_hdr_code2,msc_ocn_hdr_code2
  integer, dimension(7) :: timevec
  integer :: isat,i
  character (len=100) :: message

  do isat = 1, nsat_t
     call ReadAttribute(orb_files(isat+local_offset), "Agency", 1, strScalar= gt_agency(isat))
     call ReadAttribute(orb_files(isat+local_offset), "Created by", 1, strScalar=gt_runby(isat))
     call ReadAttribute(orb_files(isat+local_offset), "Creation date",1,strScalar=gt_rundate(isat)) 
     call ReadAttribute(orb_files(isat+local_offset), "Software version",1,strScalar=gt_swver(isat)) 
     call ReadAttribute(orb_files(isat+local_offset), "File format",1,strScalar=gt_filefmt(isat)) 
     call ReadAttribute(orb_files(isat+local_offset), "Satellite name",1,strScalar=gt_satnam(isat)) 
     call ReadAttribute(orb_files(isat+local_offset), "NRec",1,IntegerScalar=gt_datarecs(isat))
     call ReadAttribute(orb_files(isat+local_offset), "Time first obs (Sec past epoch)",7,IntegerArray = timevec)
     gt_stime(isat) = timevec(1)
     call ReadAttribute(orb_files(isat+local_offset), "Time last obs (Sec past epoch)",7,IntegerArray = timevec)
     gt_etime(isat) = timevec(1)
     call ReadAttribute(orb_files(isat+local_offset), "Sampling",1,RealScalar=gt_datarecint(isat))   
     !call ReadAttributeorb_file(isat)e%file_id, "Static_field",1 , strScalar=gt_statfieldmod(1)) 
     !call ReadAttributeorb_file(isat)e%file_id, "Tide_model",1,  strScalar=gt_oceantidemod(1)) 
     !call ReadAttributeorb_file(isat)e%file_id, "Atm_tide_model",1,  strScalar=gt_atmtidemod(1)) 
     !call ReadAttributeorb_file(isat)e%file_id, "dealiasing",1 , strScalar=gt_dealiasmod(1)) 
     !call ReadAttributeorb_file(isat)e%file_id, "dealiasing", 1, strScalar=gt_dealiasmod(1)) 
     !call ReadAttributeorb_file(isat)e%file_id, "reference framce",1,strScalar=coorspace) 
     call ReadAttribute(orb_files(isat+local_offset), "reference frame", 1, strScalar = gt_datarecref(isat))
     !call ReadAttribute(orb_files(isat), "t_in",1,  IntegerScalar=tin) 
     call ReadAttribute(orb_files(isat+local_offset), "ICs terrestrial",10 ,  RealArray=satics_t(2:11,isat)) 
     !call ReadAttribute(orb_files(isat), "ICs inertial",10, RealArray=incoor) 
     call ReadAttribute(orb_files(isat+local_offset), "Acc Scale", 3, RealArray=apr_scalebias(1:3,isat)) 
     call ReadAttribute(orb_files(isat+local_offset), "Acc Bias",3, RealArray=apr_scalebias(4:6,isat)) 
     call ReadAttribute(orb_files(isat+local_offset), "GPS mag",1, RealScalar=GPSantoff(1,isat)) 
     call ReadAttribute(orb_files(isat+local_offset), "GPS cos",3, RealArray=GPSantoff(2:4,isat)) 
     call ReadAttribute(orb_files(isat+local_offset), "GPS off",3, RealArray=GPSantoff(5:7,isat))
     call ReadAttribute(orb_files(isat+local_offset), "NMSC",1,  IntegerScalar=gt_nmcs(isat)) 
     call ReadAttribute(orb_files(isat+local_offset), "NTMSC",1,  IntegerScalar=total_ocean_prim) 
     call ReadAttribute(orb_files(isat+local_offset), "MSC fname", 1,  strScalar = combined_mascon_file)
     !print*, combined_mascon_file
     call ReadAttribute(orb_files(isat+local_offset), "MSC code", 1, strScalar = msc_hdr_code2)
     call ReadAttribute(orb_files(isat+local_offset), "Ocean fname", 1,   strScalar = ocean_mascon_file)
     call ReadAttribute(orb_files(isat+local_offset), "Ocean Code", 1,   strScalar = msc_ocn_hdr_code2)
  enddo

  msc_ocn_hdr_code=trim(msc_ocn_hdr_code2)
  msc_hdr_code=trim(msc_hdr_code2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !******************* VALIDATE AND STORE INPUT *******************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Number of epochs
  !    print *, gt_datarecs(1) , gt_datarecs(nsat_t)
  if(gt_datarecs(1) /= gt_datarecs(nsat_t)) call status_update('FATAL','GRACEFIT','header_readGTORB',' ',&
       'Inconsistency in epoch number between the GTORB files',0)
  nepochs_t = gt_datarecs(1)
  nKBobs_t = nepochs_t ! Set maximum number of Kband observations
  ! Epoch interval
  if(gt_datarecint(1) /= gt_datarecint(nsat_t)) call status_update('FATAL','GRACEFIT','header_readGTORB',' ',&
       'Inconsistency in epoch interval between the GTORB files',0)
  if(gt_datarecint(1) /= epoch_interval) then
     write(message,'(a,f6.2,a,f6.2,a)')'Inconsistency in epoch interval (',epoch_interval,') and JPL data setup (' &
          , gt_datarecint(1),')'
     call status_update('FATAL','GRACEFIT','header_readGTORB',' ', message,0)
  endif
  ! Check GTORB ref frame
  if(gt_datarecref(1) /= gt_datarecref(nsat_t)) call status_update('FATAL','GRACEFIT','header_readGTORB',' ',&
       'Inconsistency in reference frame between the GTORB files',0)
  ! Start and end times
  if(gt_stime(1) /= gt_stime(nsat_t)) call status_update('FATAL','GRACEFIT','header_readGTORB',' ',&
       'Inconsistency in start times between the GTORB files',0)
  if(gt_etime(1) /= gt_etime(nsat_t)) call status_update('FATAL','GRACEFIT','header_readGTORB',' ',&
       'Inconsistency in end times between the GTORB files',0)
  if((gt_etime(1)-gt_stime(1))/gt_datarecint(1)+1 /= nepochs_t) call status_update('FATAL','GRACEFIT','header_readGTORB',' ',&
       'Inconsistency in start and end times and number of epochs in GTORB files',0)
  starting_epoch = dble(gt_stime(1))
  spanhrs = nepochs_t*epoch_interval/3600.d0  !   Set number of hours
  nGPSobs_t = 0
  do isat = 1 ,nsat_t
     nGPSobs_t(isat) = nepochs_t/gnv_interval + 1      !   Set maximum number of GPS observations (updated in input_readGNV1B)
  enddo
  ! Mascons
  if(gt_nmcs(1) /= gt_nmcs(nsat_t)) call status_update('FATAL','GRACEFIT','header_readGTORB',' ',&
       'Inconsistency in number of mascons between the GTORB files',0)
  nmascons_t = nmascons_t*int(gt_nmcs(1)) ! NOTE: Before computation: mascons on -> nmascons_t = 1; mascons off -> nmascons_t = 0
  ! PT140616: now we can calculate the pointer for the tidal mascons
  if(est_msc_tides == 1)imsctide = imascons + nmascons_t

  ! Total number of parameters and parameter partials stored in rvec
  nrvecprm_t = 6 !@# Hardcoded until GTORB header indicates how many there are (not parameter because this may become dependent on GTORB file)
  if(nrvecprm_t < 6) call status_update('FATAL','GRACEFIT','header_readGTORB',' ',&
       'GTORB files must have at least positions and velocities for every coordinate',0)
  if(nrvecprm_t < ngobs_t) call status_update('FATAL','GRACEFIT','header_readGTORB',' ',&
       'GTORB files must have at least as many observations as g-observations considered',0)
  if(nrvecprm_t < nCoM_t) call status_update('FATAL','GRACEFIT','header_readGTORB',' ',&
       'GTORB files must have at least as many observations as antenna offset conditions considered',0)
  ! PT140616: added the mascon tidal amplitudes to nrvec_t (2 amplitudes per constituent per active mascon x number of ICs)
  nrvec_t = nrvecprm_t + norbprm_t*nrvecprm_t + nmascons_t*nrvecprm_t + nmsc_tid_constit_t*2*nrvecprm_t
  ! Update nparam_t to include mascons
  nparam_t = nparam_t + nmascons_t
  !   set mascon names (to be used in output)
  ! RM190715: changed mascon numbers to 5-digit
  do i = imascons, nparam_t
     write(prmnam(i),'(i5,a,i5.5,a)') i,'. MC',i-imascons+1,' (m)'
  enddo
  ! PT140616: update nparam_t to include mascon tidal amplitudes
  !   if(est_msc_tides == 1 .and. nmascons_t > 0 )then 
  !     nparam_t = nparam_t + nmsc_tid_constit_t*2
  !     char_mcon_tides(1) = "M2"
  !     char_mcon_tides(2) = "O1"
  !     char_mcon_tides(3) = "S2"
  !     char_mcon_tides(4) = "K1"
  !     char_mcon_tides(5) = "K2"

  ! and set mascon tidal amplitude names
  !     count_mcon_tides = imsctide
  !     prmnam(imsctide:imsctide+nmsc_tid_constit_t*2 -1) = " "
  !     do i=1, nmascons_t
  !print*,'header_lib mascon and mcon_tide1',i,mcon_tide1(i,1),LUGT(1+lu_offset),lu_offset
  !       if(mcon_tide1(i,1) > 0)then
  !         do j=1,max_mcon_tides
  !print*,'mascon,j,bitmap',i,j,kbit(mcon_tide1(i,1),j)
  !           if(bitmap(mcon_tide1(i,1),j))then
  !             write(prmnam(count_mcon_tides)  ,'(i5,a,i4.4,a1,a2,a)')count_mcon_tides,'. TMC',i,' ',char_mcon_tides(j),'sin (m)'
  !             write(prmnam(count_mcon_tides+1),'(i5,a,i4.4,a1,a2,a)')count_mcon_tides+1,'. TMC',i,' ',char_mcon_tides(j),'cos (m)'
  !             count_mcon_tides = count_mcon_tides + 2
  !           endif
  !         enddo
  !       endif
  !     enddo
  !   else
  call status_update('WARNING','GRACEFIT','readbinGTORBH5',' ','   **** No tidal mascon amplitudes to estimate',0)
  !    endif 



end subroutine header_readh5GTORB



!********************************************************************************************************************************
!********************************************************************************************************************************
! header_readGNV1B: Reads in data from the header of the GNV1B file(s) and checks to make sure
!                   files are consistent with each other, with the command file and with the
!                   GTORB file(s)
! Author: Unknown
!         (modified by Thomas Greenspan)
!********************************************************************************************************************************

subroutine header_readGNV1B()
  use gracefit_mod

  !use gracefit_mod
  implicit none
   
  include 'input.h90'

  !********************  Variable declarations ********************

  integer*4 :: isat           ! Counter variable that runs through satellites
  integer*4 :: i              ! Counter variable
  integer*4 :: ierr,jerr      ! IOSTAT errors
  integer*4 :: indx           ! Position in string
  integer*4 :: trimlen        ! Function that returns the length of a string
  character(256)   :: line    ! Line read from file
  character(256)   :: cval    ! String read from line
  double precision :: val     ! Value read from line
  character(32)    :: cmd     ! Commands read from file header lines.

  ! Data read in from header
  character(80)    :: gn_agency(nsat_t),gn_runby(nsat_t),gn_rundate(nsat_t),gn_models(nsat_t)
  character(8)     :: gn_filefmt(nsat_t),gn_satnam(nsat_t)
  integer*4        :: gn_hdrrecs(nsat_t),gn_datarecs(nsat_t)
  integer*4        :: gn_stime(nsat_t),gn_etime(nsat_t)
  !****************************************************************

  ! Initialize some variables
  gn_hdrrecs  = 0
  gn_datarecs = 0
  gn_stime    = 0.d0
  gn_etime    = 0.d0
  !****************************************************************

  !************************* PARSE HEADER *************************

  do isat = 1, nsat_t

     do while (ierr == 0)
        !   Try to read a header record
        read (LUGN(isat),'(a)',iostat=ierr) line
        if( ierr /= 0 .or. trimlen(line) <= 0 ) cycle

        !   Find the index of the first character past the : separator
        indx = 0
        do while(indx < 256)
           indx = indx + 1
           if(line(indx:indx) == ':') exit
        enddo
        indx = indx + 1  ! index one past ':'

        ! Store header record label and make sure it is all uppercase
        cmd = adjustl(line(1:indx-2))
        call lwr_to_upper(cmd)

        !   See which header line is found
        if (cmd(1:15)      == 'PRODUCER AGENCY'             ) then
           !   Get name of agency
           gn_agency(isat) = line(indx:trimlen(line))
        else if (cmd(1:20) == 'PRODUCER INSTITUTION'        ) then
           !   Get name of run_by
           gn_runby(isat) = line(indx:trimlen(line))
        else if (cmd(1:19) == 'PRODUCT CREATE TIME'         ) then
           gn_rundate(isat) = line(indx:trimlen(line))
        else if (cmd(1:16) == 'SOFTWARE VERSION'            ) then
        else if (cmd(1:11) == 'FILE FORMAT'                 ) then
           read(line(32:33),*)i
           if ( i == 1 ) then
              gn_filefmt(isat) = 'ASC'
           else if( i == 0 ) then
              gn_filefmt(isat) = 'BIN'
           else
              call status_update('WARNING','GRACEFIT','header_readGNV1B',' ','Format of GNV1B file unknown',i)
           endif
        else if (cmd(1:14) == 'SATELLITE NAME'              ) then
           gn_satnam(isat) = adjustl(line(indx:trimlen(line)))
        else if (cmd(1:24) == 'NUMBER OF HEADER RECORDS'    ) then       !   Get the number of header lines in the file
           read(line(32:38),*)gn_hdrrecs(isat)
        else if (cmd(1:22) == 'NUMBER OF DATA RECORDS'      ) then       !   Get the number of data records in the file
           read(line(32:38),*)gn_datarecs(isat)
        else if (cmd(1:11) == 'MODELS USED'                 ) then
           gn_models(isat) = line(indx:trimlen(line))
        else if (cmd(1:21) == 'TIME EPOCH (GPS TIME)'       ) then       !   Not decoded
        else if (cmd(1:30) == 'TIME FIRST OBS(SEC PAST EPOCH)') then     !   Get the start time of the integration
           read(line(32:48),*)val
           gn_stime(isat) = nint(val)
        else if (cmd(1:29) == 'TIME LAST OBS(SEC PAST EPOCH)') then      !   Get the end time of the integration
           read(line(32:48),*)val
           gn_etime(isat) = nint(val)
        else if (cmd(1:18) == 'DATA RECORD FORMAT'          ) then       !   Format not read into a variable
        else if (cmd(1:15) == 'INPUT FILE NAME'             ) then       !   Files cols not read into a variable!
        else if (cmd(1:8) == 'FILENAME'                     ) then       !   Not decoded
        else if (cmd(1:13) == 'END OF HEADER'               ) then
           exit
        endif
     enddo  ! End of header loop
  enddo  ! End of satellite loop
  !****************************************************************

  !************************ VALIDATE INPUT ************************

  ! Satellites
  if((gn_satnam(1)(7:7) /= satnam(1)).or.(gn_satnam(nsat_t)(7:7) /= satnam(nsat_t))) &
       call status_update('FATAL','GRACEFIT','header_readGNV1B',' ',&
       'Inconsistency in satellites between command file and GNV1B file(s)',0)
  ! Number of epochs
  if(gn_datarecs(1) /= gn_datarecs(nsat_t)) call status_update('FATAL','GRACEFIT','header_readGNV1B',' ',&
       'Inconsistency in record number between the GNV1B files',0)
  if(gn_datarecs(1) < nepochs_t) call status_update('FATAL','GRACEFIT','header_readGNV1B',' ',&
       'Not enough records in GNV1B file(s)',0)

  ! Start and end times
  if(gn_stime(1) /= gn_stime(nsat_t)) call status_update('FATAL','GRACEFIT','header_readGNV1B',' ',&
       'Inconsistency in start times between the GNV1B files',0)
  if(gn_etime(1) /= gn_etime(nsat_t)) call status_update('FATAL','GRACEFIT','header_readGNV1B',' ',&
       'Inconsistency in end times between the GNV1B files',0)
  if(gn_stime(1) > starting_epoch) call status_update('FATAL','GRACEFIT','header_readGNV1B',' ',&
       'Inconsistency between start time in GNV1B file(s) and in GTORB file(s)',0)
  if(gn_etime(1) < starting_epoch+(nepochs_t-1)*epoch_interval) call status_update('FATAL','GRACEFIT','header_readGNV1B',' ',&
       'Inconsistency between end time in GNV1B file(s) and in GTORB file(s)',0)
  !****************************************************************

  return
end subroutine header_readGNV1B

!********************************************************************************************************************************
!********************************************************************************************************************************
! header_readGNV1B: Reads in data from the header of the GNV1B file(s) and checks to make sure
!                   files are consistent with each other, with the command file and with the
!                   GTORB file(s)
! Author: Unknown
!         (modified by Thomas Greenspan)
!********************************************************************************************************************************

subroutine header_readKBR1B()
  use gracefit_mod

  !use gracefit_mod
  implicit none
   
  include 'input.h90'

  !********************  Variable declarations ********************

  integer*4 :: isat           ! Counter variable that runs through satellites
  integer*4 :: i              ! Counter variable
  integer*4 :: ierr,jerr      ! IOSTAT errors
  integer*4 :: indx           ! Position in string
  integer*4 :: trimlen        ! Function that returns the length of a string
  character(256)   :: line    ! Line read from file
  character(256)   :: cval    ! String read from line
  double precision :: val     ! Value read from line
  character(32)    :: cmd     ! Commands read from file header lines.

  ! Data read in from header
  character(80)    :: kb_agency,kb_runby,kb_rundate
  integer*4        :: kb_hdrrecs,kb_datarecs
  character(8)     :: kb_filefmt,kb_proclevl
  integer*4        :: kb_stime,kb_etime
  !****************************************************************

  ! Initialize some variables
  kb_hdrrecs  = 0
  kb_datarecs = 0
  kb_stime    = 0
  kb_etime    = 0
  !****************************************************************

  !************************* PARSE HEADER *************************

  do while (ierr == 0)
     !   Try to read a header record
     read (LUKB,'(a)',iostat=ierr) line
     if( ierr /= 0 .or. trimlen(line) <= 0 ) cycle

     !   Find the index of the first character past the : separator
     indx = 0
     do while(indx < 256)
        indx = indx + 1
        if(line(indx:indx) == ':') exit
     enddo
     indx = indx + 1  ! index one past ':'

     ! Store header record label and make sure it is all uppercase
     cmd = adjustl(line(1:indx-2))
     call lwr_to_upper(cmd)

     !   See which header line is found
     if (cmd(1:15)     == 'PRODUCER AGENCY'             ) then    !   Get name of agency
        kb_agency = line(indx:trimlen(line))
     else if (cmd(1:20) == 'PRODUCER INSTITUTION'        ) then  !   Get name of run_by
        kb_runby = line(indx:trimlen(line))
     else if (cmd(1:23) == 'PRODUCT CREATE END TIME'     ) then
        kb_rundate = line(indx:trimlen(line))
     else if (cmd(1:16) == 'SOFTWARE VERSION'            ) then  !    Not decoded
     else if (cmd(1:11) == 'FILE FORMAT'                 ) then
        read(line(32:33),*)i
        if (i == 1) then
           kb_filefmt = 'ASC'
        else if( i == 0 ) then
           kb_filefmt = 'BIN'
        else
           call status_update('WARNING','GRACEFIT','header_readKBR1B',' ','Format of KBR1B file unknown',i)
        endif
     else if (cmd(1:24) == 'NUMBER OF HEADER RECORDS'    ) then   !   Get the number of header lines in the file
        read(line(32:36),*)kb_hdrrecs
     else if (cmd(1:22) == 'NUMBER OF DATA RECORDS'      ) then   !   Get the number of data records in the file
        read(line(32:38),*)kb_datarecs
     else if (cmd(1:21) == 'TIME EPOCH (GPS TIME)'       ) then   !   Not decoded
     else if (cmd(1:30) == 'TIME FIRST OBS(SEC PAST EPOCH)') then !   Get the start time of the integration
        read(line(32:48),*)val
        kb_stime = int(val)
     else if (cmd(1:29) == 'TIME LAST OBS(SEC PAST EPOCH)') then  !   Get the end time of the integration
        read(line(32:48),*)val
        kb_etime = int(val)
     else if (cmd(1:15) == 'INPUT FILE NAME'             ) then   !   Files cols not read into a variable!
     else if (cmd(1:8) == 'FILENAME'                     ) then   !   Not decoded
     else if (cmd(1:13) == 'PROCESS LEVEL'               ) then
        kb_proclevl = line(indx:trimlen(line))
     else if (cmd(1:13) == 'END OF HEADER'               ) then
        exit
     endif
  enddo  ! End of header loop
  !****************************************************************

  !************************ VALIDATE INPUT ************************

  ! Number of epochs
  if(kb_datarecs < nepochs_t) call status_update('WARNING','GRACEFIT','header_readKBR1B',' ',&
       'Incomplete records in KBR1B file',0)
  ! Start and end times
  if(kb_stime > starting_epoch) call status_update('WARNING','GRACEFIT','header_readKBR1B',' ',&
       'Inconsistency between start time in KBR1B file and in GTORB file(s)',0)
  if(kb_etime < starting_epoch+(nepochs_t-1)*epoch_interval) call status_update('WARNING','GRACEFIT','header_readKBR1B',' ',&
       'Inconsistency between end time in KBR1B file and in GTORB file(s)',0)
  !****************************************************************

  return
end subroutine header_readKBR1B

!********************************************************************************************************************************
