!********************************************************************************************************************************
!  File: output_lib.f90
!
!  Purpose: Set of subroutines used to write out the solution or data that should be
!           when running the program. The subroutines should be self contained and 
!           weakly linked to rest of code.
!
!  Author: Thomas Greenspan
!          (some subroutines taken or inspired from subroutines by Simon McClusky or Paul Tregoning)
!
!  API:
!       output_openFile      : Opens a files
!       output_openAllFiles  : Opens all files relevant to GRACEFIT
!       output_writeVCV      : Write solution to .vcv file
!       output_writeFIT      : Write pre- and post-fit residuals as well as parameter information to .fit file
!       output_writeRMS      : Write general residual information to .rms file
!       output_writeSVS      : Write adjusted IC's to .svs file
!       output_writeJPL      : Write jpl kbrr residuals to .resid file
!       output_writeCOR      : Write correlations to .corr file
!       output_writeNORM     : Write normal equations and RHS to .norm file   PT131028
!       outout_codes2labels  : convert between abbreviated and full labels for parameter names PT131028
!       output_wtiteAmat     : write binary A matrix to .Amat file RM210203
!       output_writebmat     : write binary b matrix to .bmat file RM210203
!
!  July 29, 2013
!
!********************************************************************************************************************************
!********************************************************************************************************************************
! output_openFile: Opens a file given the file unit number, name and status as well as the
!                  name of the program and subroutine, the type of error if there is one and
!                  the message written in case of an error.
! Author: Thomas Greenspan
!********************************************************************************************************************************

subroutine output_openFile(unit_num,file_name,file_status,program_name,error_type,error_message,access_type,rec_len)

  implicit none

  !********************  Variable declarations ********************

  integer*4,    intent(in) :: unit_num       ! File unit number to be opened
  character(*), intent(in) :: file_name      ! Name of file to be opened
  character(*), intent(in) :: file_status    ! Status of file to be opened
  character(*), intent(in) :: program_name   ! Program name
  character(*), intent(in) :: error_type     ! Type of error if any occurs
  character(*), intent(in) :: error_message  ! Message in case of error
  ! PT131018: add info to permit the opening of binary files
  character(*), intent(in) :: access_type    ! sequential (ascii) or direct (binary)
  integer*4,    intent(in) :: rec_len        ! record length for binary file

  integer*4 :: trimlen
  integer*4 :: ioerr
  character*200 message
  !****************************************************************

  ! Open file and report stat error if there is one
  if(trim(access_type) == "sequential") then
     open(unit=unit_num,file=trim(file_name),status=file_status,iostat=ioerr)
  else if (trim(access_type) == "direct")then
     open(unit=unit_num,file=trim(file_name),status=file_status,iostat=ioerr,access="direct",recl=rec_len)
  else
     write(message,'(a,a,a)')'Unknown file access type ("',access_type,'"). Must be "sequential" or "direct".'
     call status_update('FATAL',program_name,'output_openFILE',trim(file_name),message,0)
  endif

  if (ioerr /= 0) call status_update(error_type,program_name,'output_lib',trim(file_name),error_message,ioerr)
  !****************************************************************

  return
end subroutine output_openFile

!********************************************************************************************************************************
!********************************************************************************************************************************
! output_openAllFiles: Opens all the files needed for output of the program gracefit given the 
!                      the name of the rms file and the name of the plot file
! Author: Thomas Greenspan
!********************************************************************************************************************************

subroutine output_openAllFiles(rmsnam)
  use gracefit_mod
  implicit none

   
  include 'output.h90'

  !********************  Variable declarations ********************

  character(*), intent(inout) :: rmsnam    ! Name of rms file taken from command-line

  integer*4     :: isat               ! Counter that runs through satellites
  character(64) :: fitnam             ! Name of fit file
  character(64) :: svsnam             ! Name of svs.apr file
  character(64) :: vcvnam             ! Name of vcv file
  character(64) :: residnam           ! Name of residual file
  character(64) :: corrnam            ! Name of correlations file
  character(64) :: pltkbnam           ! Name of plot kb file
  character(64) :: pltkbnacc         ! Name of plot kb file acc
  character(64) :: pltsatnam(nsat_t)  ! Name of plot sat file(s)
  integer*4 :: iwhere   ! Counter variable uset to indicate '.' in rms file in command line argument
  integer*4 :: index    ! Function that returns the index of a certain character in a string
  integer*4 :: last_nonblank    ! Function that returns the length of a string
  !****************************************************************

  !  Create file names for the fit file, svs.apr file, rms file, vcv file, resid file, plot file and correlation file
  iwhere = index(rmsnam,'.')
  !   cover case where user leaves off "."
  if (iwhere == 0 ) then
     iwhere = last_nonblank(rmsnam) + 1
     rmsnam(iwhere:iwhere) = "."
     rmsnam(iwhere+1:iwhere+3) = "rms"
  endif
  fitnam   = rmsnam(1:iwhere)//"fit"
  svsnam   = rmsnam(1:iwhere)//"svs"
  vcvnam   = rmsnam(1:iwhere)//"vcv"
  residnam = rmsnam(1:iwhere)//"resid"
  corrnam  = rmsnam(1:iwhere)//"corr"
  pltkbnam = "plt_"//rmsnam(1:iwhere)//"kb"
  pltkbnacc = "plt_"//rmsnam(1:iwhere)//"kba"
  do isat = 1, nsat_t
     pltsatnam(isat) = "plt_"//rmsnam(1:iwhere)//satnam(isat)
  enddo

  !************************** OPEN FILES **************************

! PT220125: no longer output the files we never use: .svs .rms .resid
  ! Open  ouput files.
  !   open the rms file
!  call output_openFile(LURMS,rmsnam,'unknown','GRACEFIT','FATAL','Error opening the rms file: ','sequential',0)
  !   open the SVS apr file
!  call output_openFile(LUSVS,svsnam,'unknown','GRACEFIT','FATAL','Error opening the svs.apr file: ','sequential',0)
  !   open the solution file
  call output_openFile(LUVCV,vcvnam,'unknown','GRACEFIT','FATAL','Error opening the vcv file: ','sequential',0)
  !   open the fit-file
  call output_openFile(LUFIT,fitnam,'unknown','GRACEFIT','FATAL','Error opening the fit file: ','sequential',0)
  !   open the jpl kbrr residual file
!  call output_openFile(LUJPL,residnam,'unknown','GRACEFIT','WARNING','Error opening the jpl kbrr residual file: ','sequential',0)
  !   open the correlation file
  call output_openFile(LUCOR,corrnam,'unknown','GRACEFIT','WARNING','Error opening the correlation file: ','sequential',0)

  !   open plot kb file
  call output_openFile(LUGN_KB,pltkbnam,'unknown','GRACEFIT','WARNING','Error opening the plot kb file: ','sequential',0)
  call output_openFile(LUGN_KBA,pltkbnacc,'unknown','GRACEFIT','WARNING','Error opening the plot kb file: ','sequential',0)
  !   open plot sat file(s)
  do isat = 1, nsat_t
     call output_openFile(LUGN_PLT(isat),pltsatnam(isat),'unknown','GRACEFIT','WARNING',&
          'Error opening the plot sat '//satnam(isat)//' file: ','sequential',0)
  enddo
  !****************************************************************

  return
end subroutine output_openAllFiles

!********************************************************************************************************************************
!********************************************************************************************************************************
! output_writeVCV: Writes out the vcv file given the a priori and adjustment vectors, the inverted
!                  form of the right-hand side of the LS equation, the name of the GTORB files and
!                  the total number of epochs run through
! Author: Lydie Lescarmontier
!         (Modified by Thomas Greenspan)
!********************************************************************************************************************************

subroutine output_writeVCV(gt_fnam,apr_prm,kbrr_misfit_num,adjust,normeq,prmnam,prmnam_size)
  use gracefit_mod
  implicit none

   
  include 'output.h90'

  !********************  Variable declarations ********************

  character(*), intent(in)     :: gt_fnam(maxsat)           ! Name of GTORB files
  double precision, intent(in) :: apr_prm(nparam_t)         ! A priori values for parameters
  integer*4, intent(in)        :: kbrr_misfit_num           ! Number of kbrr misfits
  double precision, intent(in) :: adjust(nparam_t)          ! Adjustments to apriori values
  double precision, intent(in) :: normeq(nparam_t,nparam_t) ! The inverted left side of the LS algorithm. P_t*W*P
  integer*4        ,intent(in)    :: prmnam_size
  character*30,     intent(inout) :: prmnam(prmnam_size)

  integer*4        :: isat            ! Counter that runs through satellites
  integer*4        :: i,j             ! Counter variables
  double precision :: middle_epoch    ! The epoch (in grace seconds) in the middle of the span of epochs looped over
  double precision :: sigma(nparam_t) ! Sigmas
  character(32)    :: dattim          ! Date the program was run
  integer*4        :: date(5)         ! Date being looked at by program
  double precision :: sec             ! Seconds part of date (being looked at)
  !****************************************************************

  ! Set values
  !   Set what the middle epoch is in grace seconds
  middle_epoch = starting_epoch+(nepochs_t-1)*epoch_interval/2.d0
  !   Set sigmas
  do i=1,nparam_t
     sigma(i) = dsqrt(normeq(i,i))
  enddo

  !************************* WRITE HEADER *************************

  ! Write out headers
  call datetime(dattim)
! RM190329: add in option to write out V2 or V3 VCV file
  write(LUVCV, '(a,a,a,a19)') 'V2 ',calling_prog,' VCV solution:  run ',dattim
!  write(LUVCV, '(a,a,a,a19)') 'V3 ',calling_prog,' VCV solution:  run ',dattim
  write(LUVCV, '(a,a,1x,a)') 'Reference GTORB files: ', (trim(gt_fnam(i)), i = 1, 2)
  write(LUVCV, '(a,f6.1,a)') 'Reference GTORB files span: ',spanhrs,' hrs'
  call gsec_to_ymdhms(middle_epoch, date, sec)
  write(LUVCV,'(a,a,a,i4,4(a1,i2.2),a1,f05.2)') 'Reference ',calling_prog,' solution epoch: ' &
       ,date(1),"-",date(2),"-",date(3)," ",date(4),"-",date(5),"-",sec
  !PT150509: don't write out number of tidal mascon amplitudes unless we estimated them
  if(est_msc_tides > 0)then
     write(LUVCV,'(a32,3i7)')'Number of estimated parameters: ',nparam_t,nmascons_t,nmsc_tid_constit_t
  else
     write(LUVCV,'(a32,3i7)')'Number of estimated parameters: ',nparam_t,nmascons_t,0
  endif
  ! Write out information on uncertainties and number of epochs that have been used for various observations
  ! Uncertainties do not have any appropriate values as of now
  !    write(LUVCV,'(a,f8.6)') 'KBRR uncertainty (um): ',
  !    write(LUVCV,'(a,f8.6)') 'GPS position uncertainty (m): ',
  write(LUVCV,'(a,i4)')   'Number of KB observations and used: ',nKBobs_t
  write(LUVCV,'(a,i4)')   'Number of missing KB observations: ',nepochs_t - nKBobs_t
  write(LUVCV,'(a,i4)')   'Number of KBRR misfits: ',kbrr_misfit_num
  write(LUVCV,'(a,i4)')   'Number of KBRR observations used: ',nKBobs_t - kbrr_misfit_num
  do isat = 1, nsat_t
     write(LUVCV,'(a,a1,a,i4)') 'Number of GPS observations used for GRACE ',satnam(isat),': ',nGPSobs_t(isat)
     write(LUVCV,'(a,a1,a,i4)') 'Number of relevant missing GPS observations for GRACE ',satnam(isat),': ',&
          nepochs_t/gnv_interval + 1 - nGPSobs_t(isat)
  enddo
  write(LUVCV,'(a,i4)')   'Total Number of epochs: ',nepochs_t
  !****************************************************************

  !************************ WRITE SOLUTION ************************

  ! NOTE: If we want to make a separate subroutine that prints out the solutions for gracefit and the gracekal forwards and backwards filters, only the solution should be taken care of by the subroutine in question. The subroutine would be called here.

  write(LUVCV,'(/,a)') ' SOLUTION A PRIORI AND VECTOR:'
  write(LUVCV,'(a)') ' PARAMETER                     A PRIORI             VECTOR            SIGMA'   
  ! Write out names of parameters, a priori, vector and sigmas
  do i = 1, nparam_t
     write(LUVCV,'(a,f17.7,f17.7,f17.7)') prmnam(i),apr_prm(i),apr_prm(i)+adjust(i),sigma(i)
  enddo
  ! Write the vcv solution            
  write(LUVCV,*)' VCV SOLUTION '
  ! PT220127: turn this off for now - we are not using it at present
  call status_update('STATUS',calling_prog,'output_writeVCV',' ','Not outputting the whole VCV. Code commented out',0)
  !do i = 1, nparam_t
  !   write(LUVCV,*) ((normeq(j, i)), j = 1, nparam_t)   
  !enddo
  !****************************************************************

  close (LUVCV)

  return
end subroutine output_writeVCV

!********************************************************************************************************************************
!********************************************************************************************************************************
! output_writeFIT: Calculate the pre- and post-fit rms values and print these as well as the
!                  parameters to a fit file. Input is the GTORB file names, the names of the
!                  observations, the a priori values for paramters and their adjustments, the
!                  matrix of partials as well as the observed minus calculated values and the
!                  inverted left-hand side of the LS algorithm
! Author: Thomas Greenspan
!         (taken in large part from a section of the subroutine write_summary by unknown author)
!
! MODS:
! PT160229: passed in and wrote out the a priori sigmas on the observations
! PT190917: write name of regularisation file to .fit file header
!********************************************************************************************************************************

subroutine output_writeFIT(gt_fnam,apr_prm,apr_wght,adjust,kbrr_misfit_num,Amat,pre_omc,normeq,post_omc &
                           ,prmnam,prmnam_size,regularisation_file)
  use gracefit_mod
  implicit none

   
  include 'output.h90'

  !********************  Variable declarations ********************

  character(*), intent(in)      :: gt_fnam(maxsat)                 ! Name of GTORB files
  double precision, intent(in)  :: apr_prm(nparam_t)               ! A priori values for parameters
  double precision, intent(in)  :: apr_wght(nobs_t)                ! A priori weights on observations
  double precision, intent(in)  :: adjust(nparam_t)                ! Adjustments to a priori values for parameters
  integer*4, intent(in)         :: kbrr_misfit_num                 ! Number of kbrr misfits
  double precision, intent(in)  :: Amat(nobs_t*nepochs_t,nparam_t) ! Partials of observables with respect to the parameters
  double precision, intent(in)  :: pre_omc(nobs_t,nepochs_t)       ! prefit Observed Minus Computed values
  double precision, intent(in)  :: normeq(nparam_t,nparam_t)       ! The inverted left-hand side of the LS algorithm. P_t*W*P
  double precision, intent(in) :: post_omc(nobs_t,nepochs_t)      ! postfit Observed Minus Computed values
  integer*4        ,intent(in)    :: prmnam_size
  character*30,     intent(inout) :: prmnam(prmnam_size)
  character*100    ,intent(in)    :: regularisation_file           ! name of regularisation file used in the inversion


  integer*4        :: iepoch, ordered(nepochs_t)             ! Counter that runs through epochs
  integer*4        :: isat                ! Counter that runs through satellites
  integer*4        :: i,j, k              ! Counter variables
  character(32)    :: dattim              ! Date the program was run
  character(256)   :: message             ! Message to be printed out concerning RR postfit rms
  integer*4        :: tmp_indx

  ! Variables for writing the pre- and postfit rms values
  double precision :: prefit_sum(ngobs_t*nsat_t+nkobs_t)   ! Sum of prefit rms values for each observation
  double precision :: postfit_sum(ngobs_t*nsat_t+nkobs_t)  ! Sum of postfit rms values for each observation
  double precision :: prefit_obs_sum(2),postfit_obs_sum(2) ! Total prefit and postfit rms values for pos and vel
  double precision :: prefit_rms,postfit_rms              ! Variable in which prefit and posfit are placed in and printed out from

  ! Variables for writing the parameters
  double precision :: sigma(nparam_t)     ! Sigmas
  double precision :: unit_conv(nparam_t) ! Factor necessary to apply wanted unit conversions
  double precision :: fract               ! Fraction of the adjustment over the uncertainty for a parameter
  double precision :: post_prm            ! A priori value for a parameter after the adjustment is applied
  !****************************************************************

  ! Set sigmas
  do i=1, nparam_t
     sigma(i) = dsqrt(normeq(i,i))
  enddo
  !************************* WRITE HEADER *************************

  call datetime(dattim)
  write(LUFIT, '(a,a,a,a19)') '# V2 ',calling_prog,' FIT solution: run ',dattim
  write(LUFIT,'(a,a19,a,a)') '# Parameters adjustments:  run ',dattim," Regularisation file: ",regularisation_file
  write(LUFIT,'(a,a,1x,a)') '# Reference GTORB files: ',(trim(gt_fnam(i)),i=1,2)
  write(LUFIT,'(a,f6.1,a)') '# Total overlapped time range: ',spanhrs,' hrs'
  ! PT170516: reformat this to be one set of uncertainties per line. Make conditional on obs having been used.
  write(LUFIT,'(a,6f10.4)')                 '# Obs sigma GPS pos/vel:',1.0/sqrt(apr_wght(1:6))
  if(nkbrr_t /= 0)write(LUFIT,'(a,f15.5,a)')"# Obs sigma        KBRR:",1.0/sqrt(apr_wght(ikbrr))*1.e6,' um/s'
  if(nkbra_t /= 0)write(LUFIT,'(a,e15.6,a)')"# Obs sigma        KBRA:",1.0/sqrt(apr_wght(ikbra))*1.e9,' nm/s^2'
  if(ncond_t /= 0)write(LUFIT,'(a,4f15.5,a)')"# Obs sigma   shad cond:",1.0/sqrt(apr_wght(icond:icond+3)),' nm/s^2'
  if(naobs_t /= 0)write(LUFIT,'(a,3e15.5,a)')"# Obs sigma    mean acc:",1.0/sqrt(apr_wght(iacobs:iacobs+2)),' nm/s^2'

  !    write(LUFIT,100)'Obs Sigmas (RR, pos/vel, shad constraints, mean acc):'  &
  !             ,1.0/sqrt(apr_wght(ikbrr))           &  ! a priori sigmas on range rate observations
  !             ,1.0/sqrt(apr_wght(1:6))             &  ! a priori sigmas on GPS position and velocity
  !             ,1.0/sqrt(apr_wght(icond:icond+3))   &  ! a priori sigmas on Y/Z accelerometer values during shadow
  !             ,1.0/sqrt(apr_wght(iacobs:iacobs+2))    ! a priori sigmas on Y/Z accelerometer mean values being the same
  !100 format(a,14e12.3)
  !****************************************************************

  call status_update('STATUS',calling_prog,'output_writeFIT',' ', 'Calculating residuals',0)

  !******************** WRITE PREFIT STATISTICS *******************

  prefit_sum = 0.d0

  ! Calculate intermediate values
  ! PT140901: use only the epochs that were stacked in the normal equations
  do iepoch = start_neq,end_neq
     do isat = 1, nsat_t
        do i = 1, ngobs_t
           !       Calculate rms values
           prefit_sum(i+(isat-1)*ngobs_t) = prefit_sum(i+(isat-1)*ngobs_t) + pre_omc(i+igobs-1+(isat-1)*ngobs_t,iepoch)**2
           if (i <= 3)then       !   Combined GRACE A/B position prefit rms
              prefit_obs_sum(1) =  prefit_obs_sum(1) + pre_omc(i+igobs-1+(isat-1)*ngobs_t,iepoch)**2
           else if (i <= 6)then  !   Combined GRACE A/B velocity prefit rms
              prefit_obs_sum(2) =  prefit_obs_sum(2) + pre_omc(i+igobs-1+(isat-1)*ngobs_t,iepoch)**2
           endif
        enddo  ! End of observations loop
     enddo  ! End of satellite loop
     if(nkbr_t /= 0)  prefit_sum(ikbr)  = prefit_sum(ikbr)  + pre_omc(ikbr,iepoch)**2   ! Combined GRACE prefit R  rms
     if(nkbrr_t /= 0) prefit_sum(ikbrr) = prefit_sum(ikbrr) + pre_omc(ikbrr,iepoch)**2  ! Combined GRACE prefit RR rms
     if(nkbra_t /= 0) prefit_sum(ikbra) = prefit_sum(ikbra) + pre_omc(ikbra,iepoch)**2  ! Combined GRACE prefit RA rms
  enddo  ! End of epoch loop

  ! Write out prefit rms values
  write(LUFIT,'(a)') '# PreFit Statistics: '

  do isat = 1, nsat_t
     do i = igobs, igobs+ngobs_t-1
        prefit_rms = dsqrt(prefit_sum(i-igobs+1+(isat-1)*ngobs_t)/dble(nGPSobs_t(isat)))
        if (i <= 3)then
           ! write(LUFIT,'(a,a,1x,a3,1x,f11.7)') '# Prefit RMS (m)     GRACE:  ',satnam(isat),obsnam(i)(14:15),prefit_rms
           write(LUFIT,'(a,a,a,a,a,1x,f11.6)') '# Prefit RMS GRACE ', satnam(isat),' Pos ', obsnam(i)(14:15), ' (m)    : ',prefit_rms
        else if (i <= 6)then
           !write(LUFIT,'(a,a,1x,a3,1x,f11.7)') '# Prefit RMS (m/s)   GRACE:  ',satnam(isat),obsnam(i)(14:15),prefit_rms
            write(LUFIT,'(a,a,a,a,a,1x,f11.6)') '# Prefit RMS GRACE ', satnam(isat),' Vel ', obsnam(i)(14:15), ' (m/s)  :  ',prefit_rms
        endif
     enddo  ! End of observations loop
  enddo  ! End of satellite loop

  !   Write total position, velocity and acceleration rms
  prefit_rms = dsqrt(prefit_obs_sum(1)/dble(nGPSobs_t(1)*3+nGPSobs_t(2)*3))
  if(ngobs_t > 0) write(LUFIT,'(a,1x,f7.4)') '# Prefit RMS Pos (m)      :',prefit_rms
  prefit_rms = dsqrt(prefit_obs_sum(2)/dble(nGPSobs_t(1)*3+nGPSobs_t(2)*3))
  if(ngobs_t > 3) write(LUFIT,'(a,1x,f7.4)') '# Prefit RMS Vel (mm/s)   :',prefit_rms*m_mm

  !   Write Range, Range Rate and Range Acceleration rms
  if(nkbr_t /= 0)then
     prefit_rms = dsqrt(prefit_sum(ikbr)/dble(nKBobs_t))
     write(LUFIT,'(a,1x,f10.5)')               '# Prefit RMS R   (m)     :',prefit_rms
  endif
  if(nkbrr_t /= 0)then
     prefit_rms = dsqrt(prefit_sum(ikbrr)/dble(nKBobs_t-kbrr_misfit_num))
     write(LUFIT,'(a,1x,f10.5)')               '# Prefit RMS RR  (um/s)   :',prefit_rms*m_um
  endif
  if(nkbra_t /= 0)then
     prefit_rms = dsqrt(prefit_sum(ikbra)/dble(nKBobs_t))
     write(LUFIT,'(a,1x,f10.5)')               '# Prefit RMS RA  (nm/s^2) :',prefit_rms*m_nm
  endif
  !****************************************************************

  !******************* WRITE POSTFIT STATISTICS *******************

  postfit_sum = 0.d0


  ! PT140901: for the statistics, use only the epochs that were stacked in the normal equations
  do iepoch = start_neq,end_neq
     do isat = 1, nsat_t
        do i = 1, ngobs_t
           !       Calculate rms values
           postfit_sum(i+(isat-1)*ngobs_t) = postfit_sum(i+(isat-1)*ngobs_t) + post_omc(i+igobs-1+(isat-1)*ngobs_t,iepoch)**2
           if (i <= 3)then       !   Combined GRACE A/B position postfit rms
              postfit_obs_sum(1) =  postfit_obs_sum(1) + post_omc(i+igobs-1+(isat-1)*ngobs_t,iepoch)**2
           else if (i <= 6)then  !   Combined GRACE A/B velocity postfit rms
              postfit_obs_sum(2) =  postfit_obs_sum(2) + post_omc(i+igobs-1+(isat-1)*ngobs_t,iepoch)**2
           endif
        enddo
     enddo
     !   Combined GRACE postfit RR rms
     ! PT130510: eliminate the first and last 6 points to exclude the tails of our AORC generation
     if( iepoch > 6 .and. iepoch < (nepochs_t - 6) )then
        if(nkbr_t /= 0)  postfit_sum(ikbr)  = postfit_sum(ikbr)  + post_omc(ikbr,iepoch)**2   ! Combined GRACE postfit R  rms
        if(nkbrr_t /= 0) postfit_sum(ikbrr) = postfit_sum(ikbrr) + post_omc(ikbrr,iepoch)**2  ! Combined GRACE postfit RR rms
        if(nkbra_t /= 0) postfit_sum(ikbra) = postfit_sum(ikbra) + post_omc(ikbra,iepoch)**2  ! Combined GRACE postfit RA rms
     endif
  enddo

  ! Write out postfit rms values
  write(LUFIT,'(/,a)') '# PostFit Statistics: '

  do isat = 1, nsat_t
     do i = igobs, igobs+ngobs_t-1
        postfit_rms = dsqrt(postfit_sum(i-igobs+1+(isat-1)*ngobs_t)/dble(nGPSobs_t(isat)))
        if ( i <= 3 )then
           !write(LUFIT,'(a,a,1x,a3,1x,f8.6)') '# Postfit RMS Pos (m)     GRACE:  ',satnam(isat),obsnam(i)(14:15),postfit_rms
           write(LUFIT,'(a,a,a,a,a,1x,f11.6)') '# Postfit RMS GRACE ', satnam(isat),' Pos ', obsnam(i)(14:15), ' (m)    :  ', postfit_rms
        else if (i <= 6)then
           !write(LUFIT,'(a,a,1x,a3,1x,f8.6)') '# Postfit RMS (m/s)   GRACE:  ',satnam(isat),obsnam(i)(14:15),postfit_rms
           write(LUFIT,'(a,a,a,a,a,1x,f11.6)') '# Postfit RMS GRACE ', satnam(isat),' Vel ', obsnam(i)(14:15), ' (m/s)  :  ', postfit_rms
        endif
     enddo  ! End of observations loop
  enddo  ! End of satellite loop

  !   Write total position, velocity and acceleration rms
  postfit_rms = dsqrt(postfit_obs_sum(1)/dble(nGPSobs_t(1)*3+nGPSobs_t(2)*3))
  if(ngobs_t > 0)then
     write(LUFIT,'(a,1x,f8.6)') '# Postfit RMS Pos (m)      : ',postfit_rms
     write(message,'(a,a,f15.5,a)') gt_fnam(1)(7:19),'# Postfit Pos  RMS = ',postfit_rms, ' (m)'
     call status_update('STATUS',calling_prog,'output_writefit',' ',message,0)
  endif
  postfit_rms = dsqrt(postfit_obs_sum(2)/dble(nGPSobs_t(1)*3+nGPSobs_t(2)*3))
  if(ngobs_t > 3)then
     write(LUFIT,'(a,1x,f8.6)') '# Postfit RMS Vel (mm/s)   : ',postfit_rms*m_mm
     write(message,'(a,a,f15.5,a)') gt_fnam(1)(7:19),'# Postfit Vel  RMS = ',postfit_rms*m_mm, ' (mm)'
     call status_update('STATUS',calling_prog,'output_writefit',' ',message,0)
  endif

  !   Write Range, Range Rate and Range Acceleration rms
  if(nkbr_t /= 0)then
     postfit_rms = dsqrt(postfit_sum(ikbr)/dble(nKBobs_t-12))
     write(LUFIT,'(a,1x,f11.5)')              '# Postfit RMS R   (m)     : ',postfit_rms
  endif
  if(nkbrr_t /= 0)then
     postfit_rms = dsqrt(postfit_sum(ikbrr)/dble(nKBobs_t-12-kbrr_misfit_num))
     write(LUFIT,'(a,1x,f11.5)')              '# Postfit RMS RR  (um/s)   : ',postfit_rms*m_um
     ! Write message to standard output
     write(message,'(a,a,f15.5,a)') gt_fnam(1)(7:19),' Postfit   RR RMS = ',postfit_rms*m_um,' (um/sec)'
     call status_update('STATUS',calling_prog,'output_writeFIT',' ', message,0)
  endif
  if(nkbra_t /= 0)then
     postfit_rms = dsqrt(postfit_sum(ikbra)/dble(nKBobs_t-12))
     write(LUFIT,'(a,1x,f11.5)')              '# Postfit RMS RA  (nm/s^2) : ',postfit_rms*m_nm
     ! Write message to standard output
     write(message,'(a,a,f15.5,a)') gt_fnam(1)(7:19),' Postfit   RA RMS = ',postfit_rms*m_nm,' (nm/sec^2)'
     call status_update('STATUS',calling_prog,'output_writeFIT',' ', message,0)
  endif


  !****************************************************************

  !****************** WRITE PARAMETER ADJUSTMENTS *****************

  ! Set unit conversions as necessary
  do isat = 1, nsat_t
     do i = 1, nsatprm_t
        if(i <= 3)then
           unit_conv(i+(isat-1)*nsatprm_t) = m_mm
           if(index(prmnam(i+(isat-1)*nsatprm_t),'(m)') /= 0)then
             tmp_indx = index(prmnam(i+(isat-1)*nsatprm_t),'(m)')
             prmnam(i+(isat-1)*nsatprm_t)(tmp_indx:tmp_indx+3) = "(mm)"
           endif
        else if(i <= ngprm_t)then
           unit_conv(i+(isat-1)*nsatprm_t) = m_mm
           if(index(prmnam(i+(isat-1)*nsatprm_t),'(m/s)') /= 0)then
             tmp_indx = index(prmnam(i+(isat-1)*nsatprm_t),'(m/s)')
             prmnam(i+(isat-1)*nsatprm_t)(tmp_indx:tmp_indx+5) = "(mm/s)"
           endif
        else if(i <= norbprm_t)then
           unit_conv(i+(isat-1)*nsatprm_t) = 1.d0  ! No unit conversion
        else if(i <= nsatprm_t)then
           unit_conv(i+(isat-1)*nsatprm_t) = m_mm
           if(index(prmnam(i+(isat-1)*nsatprm_t),'(m)') /= 0)then
             tmp_indx = index(prmnam(i+(isat-1)*nsatprm_t),'(m)')
             prmnam(i+(isat-1)*nsatprm_t)(tmp_indx:tmp_indx+3) = "(mm)"
           endif
        endif
     enddo  ! End if satellite-specific parameter loop
  enddo  ! End of satellite loop

  !   Set unit conversion for mascons (and mascon tidal amplitudes)
  unit_conv(imascons:nparam_t) = 1.d0  ! No unit conversion

  ! Write parameters
  write(LUFIT,'(/,a)') '    PARAMETER                        A PRIORI       ADJUST         POSTFIT        SIGMA    FRACT'
  write(LUFIT,'(1x)')

  do i = 1, nparam_t
     post_prm = apr_prm(i) + adjust(i)
     fract = adjust(i)/sigma(i)
     write(LUFIT,'(a30,f17.5,f12.5,f19.5,f12.5,f9.1)') prmnam(i),apr_prm(i)*unit_conv(i),adjust(i)*unit_conv(i),&
          post_prm*unit_conv(i),sigma(i)*unit_conv(i),fract
     if(mod(i,nsatprm_t) == 0) write(LUFIT,'(1x)') ! Add divider between sections
  enddo
  !****************************************************************

  close (LUFIT)

  return
end subroutine output_writeFIT

!********************************************************************************************************************************
!********************************************************************************************************************************
! output_writeRMS: Calculates and writes out rms values (xyz, rac, means..) to RMS file given 
!                  the names of the GTORB files, the observed minus computed array and the
!                  vector of positions and velocities
! Author: Thomas Greenspan
!         (taken in large part from a section of the subroutine write_summary by unknown author)
!********************************************************************************************************************************

subroutine output_writeRMS(gt_fnam,post_omc,rvec)
  use gracefit_mod
  implicit none

   
  include 'output.h90'

  !********************  Variable declarations ********************

  character(*), intent(in)     :: gt_fnam(maxsat)                ! Name of GTORB files
  double precision, intent(in) :: post_omc(nobs_t,nepochs_t)     ! Observed Minus Computed values to be used in LS algorithm
  double precision, intent(in) :: rvec(6,nsat_t,nepochs_t)       ! Vector of positions, velocities

  integer*4        :: iepoch              ! Counter that runs through epochs
  integer*4        :: isat                ! Counter that runs through satellites
  integer*4        :: i                   ! Counter variables
  character(32)    :: dattim              ! Date the program was run
  character(256)   :: message             ! Message to be printed out

  double precision :: dot                        ! Function that returns the dot product of a vector of length 3 or less 
  double precision :: drac(3,nsat_t,nepochs_t)           ! Radial residuals
  double precision :: sumxyz(3,nsat_t),sumrac(3,nsat_t)  ! Intermediate values for xyz and rac rms calculations
  double precision :: rmsxyz(3,nsat_t),rmsrac(3,nsat_t)  ! xyz and rac rms values for each satellite
  double precision :: rmstot(nsat_t)                     ! Total rms values for each satellite
  double precision :: mean_rmsxyz(3),mean_rmsrac(3)      ! Mean xyz and rac rms values
  double precision :: mean_rmstot                        ! Mean total rms values
  !****************************************************************

  !************************* WRITE HEADER *************************

  call datetime( dattim )
  write(LURMS,'(a,a19)') 'Summary of rms differences: GRACEFIT run ',dattim
  write(LURMS,'(a,a,1x,a)') 'Reference GTORB files: ',(trim(gt_fnam(i)),i=1,2)
  write(LURMS,'(a,f6.1,a)') 'Reference GTORB files span: ',spanhrs,' hrs'
  write(LURMS,'(/,a,/,a)') 'SAT   Total     delta-X   delta-Y   delta-Z   d-Radial  d-Along   d-Cross' &
       ,'============================================================================'
  !****************************************************************

  !*********************** CALCULATE VALUES ***********************

  sumxyz = 0.d0
  sumrac = 0.d0
  ! Get residuals and statistics for each satellite
  do isat = 1, nsat_t
     do iepoch = 1, nepochs_t
        if (mod(iepoch-1,gnv_interval) == 0) then   ! Only apply when GPS observations were used
           !       Get radial, along-track,cross-track
           call xyz2rac(rvec(1:6,isat,iepoch),post_omc(1:3,iepoch),drac(1:3,isat,iepoch) )
           !       Calculate sum
           do i = 1, 3
              sumxyz(i,isat) = sumxyz(i,isat) + post_omc(i+(isat-1)*ngobs_t,iepoch)**2
              sumrac(i,isat) = sumrac(i,isat) + drac(i,isat,iepoch)**2
           enddo
        endif  ! End of "mod(gnv_interval)" if statement
     enddo  ! End of epoch loop
     !   Calculate rms values
     do i = 1, 3
        rmsxyz(i,isat) = dsqrt(sumxyz(i,isat)/dble(nGPSobs_t(1)))
        rmsrac(i,isat) = dsqrt(sumrac(i,isat)/dble(nGPSobs_t(1)))
     enddo
     !   Total rms values
     rmstot(isat) = dsqrt(dot(rmsxyz(:,isat),rmsxyz(:,isat))/3)
  enddo  ! End of satellite loop

  !   Write the mean values at the bottom of the rms file
  do i = 1, 3
     mean_rmsxyz(i) = dsqrt(sum(rmsxyz(i,1:nsat_t)**2)/dble(nsat_t))*m_mm
     mean_rmsrac(i) = dsqrt(sum(rmsrac(i,1:nsat_t)**2)/dble(nsat_t))*m_mm
  enddo
  mean_rmstot = dsqrt(sum(rmstot(1:nsat_t)**2)/dble(nsat_t))*m_mm
  !****************************************************************

  !*********************** WRITE TO RMS FILE **********************

  do isat = 1, nsat_t
     !   Write rms values for each satellite to file
     write(LURMS,'(i3,7(1x,f9.5))') isat,rmstot(isat)*m_mm,(rmsxyz(i,isat)*m_mm,i=1,3),(rmsrac(i,isat)*m_mm,i=1,3)
  enddo  ! End of satellite loop

  write(LURMS,'(a)') '============================================================================'
  write(LURMS,'(a,7(f9.5,1x),/)') 'MEAN',mean_rmstot,mean_rmsxyz,mean_rmsrac
  !   Write the total in the STATUS file for sh_gamit grep'ing
  write(message,'(a,f11.5,a)') 'Overall fit (rms) to GNV1B orbits =',mean_rmstot,' (mm)'
  call status_update('STATUS','GRACEFIT','orbits/output_writeRMS',' ', message,0)
  !****************************************************************

  close (LURMS)

  return
end subroutine output_writeRMS

!********************************************************************************************************************************
!********************************************************************************************************************************
! output_writeSVS: Calculates and writes out IC values to svs file given the vector of a priori
!                  values and the vector of adjustments to these values
! Author: Unknown
!         (Modified by Thomas Greenspan)
!********************************************************************************************************************************

subroutine output_writeSVS(apr_prm,adjust)
  use gracefit_mod
  implicit none

   
  include 'output.h90'

  !********************  Variable declarations ********************

  double precision, intent(in)    :: apr_prm(nparam_t)  ! A priori values for parameters
  double precision, intent(in)    :: adjust(nparam_t)   ! Adjustments to a priori values for parameters

  integer*4        :: isat                         ! Counter that runs through satellites
  integer*4        :: i                            ! Counter variables
  double precision :: adj_satICs(nsatprm_t*nsat_t) ! Adjusted satellite IC's
  character(32)    :: dattim                       ! Date the program was run
  !****************************************************************

  !************************* WRITE HEADER *************************

  call datetime(dattim)
  write(LUSVS,'(a,a19)') '#  svs generated by GRACEFIT run ',dattim
  write(LUSVS,'(a)') '#'
  write(LUSVS,'(a,a)') '#  Time   X_pos Y_pos Z_pos  X_dot  Y_dot   ', &
       'Z_dot acc_scl_X acc_scl_Y acc_scl_Z acc_bs_X  acc_bs_Y  acc_bs_Z'
  write(LUSVS,'(a,a)') '# (gsec)   (m)   (m)   (m)  (m/sec)(m/sec) (m/sec) ', & 
       '                             (m/sec^2) (m/sec^2) (m/sec^2)'
  write(LUSVS,'(a)') '#'
  !****************************************************************

  !*********************** WRITE TO SVS FILE **********************

  do i = 1, nsatprm_t*nsat_t
     adj_satICs(i) = apr_prm(i) + adjust(i)
  enddo  ! End of loop over parameters

  do isat = 1, nsat_t
     write(LUSVS,*) SATNAM(isat),int(starting_epoch),(adj_satics(i+(isat-1)*nsatprm_t),i=1,norbprm_t)
  enddo  ! End of satellite loop

  call status_update('STATUS','GRACEFIT','grace/output_writeSVS',' ','Successfully wrote SVS ',0)
  !****************************************************************

  close (LUSVS)

  return
end subroutine output_writeSVS

!********************************************************************************************************************************
!********************************************************************************************************************************
! output_writeJPL: Writes Range Rate residuals to JPL resid file
! Author: Thomas Greenspan
!********************************************************************************************************************************

subroutine output_writeJPL(post_omc)
  use gracefit_mod
  implicit none

   
  include 'output.h90'

  !********************  Variable declarations ********************

  double precision, intent(in) :: post_omc(nobs_t,nepochs_t)  ! Postfit Observed Minus Computed values

  integer*4        :: iepoch  ! Counter that runs through epochs
  !****************************************************************

  !*********************** WRITE TO JPL FILE **********************
  if(nkobs_t > 0)then
     do iepoch = 1, nepochs_t
        write(LUJPL,*)iepoch,post_omc(ikbrr,iepoch)*1.d3,post_omc(ikbra,iepoch)*1.d6
     enddo
  endif
  !****************************************************************

  close (LUJPL)

  return
end subroutine output_writeJPL

!********************************************************************************************************************************
!********************************************************************************************************************************
! output_writeCOR: Calculates and writes out correlation matrix for GRACE A to Correlation file
! Author: Thomas Greenspan
!         (based on code by Paul Tregoning)
!********************************************************************************************************************************

subroutine output_writeCOR(normeq)
  use gracefit_mod
  implicit none

   
  include 'output.h90'

  !********************  Variable declarations ********************

  double precision, intent(in) :: normeq(nparam_t,nparam_t)  ! The inverted left-hand side of the LS algorithm. P_t*W*P

  integer*4        :: i,j                              ! Counter variables
  double precision :: corr(nparam_t,nparam_t,nsat_t) ! Matrix of correlations
  !****************************************************************

  call status_update('STATUS','GRACEFIT','gracefit',' ',"Writing correlation file",0)
  do i = 1, nparam_t
     do j = 1, nparam_t
        corr(i,j,1) = normeq(i,j)/( dsqrt(normeq(i,i))*dsqrt(normeq(j,j)) ) !Put in satellite loop if correlations for both sats is wanted
        !  if(i > 24) print*,'i,j,corr(i,j)',i,j, corr(i,j,1)
     enddo
  enddo

  !*********************** WRITE TO COR FILE **********************

  do i = 1, norbprm_t    
     write(LUCOR,'(i4,12e15.4)')i,(normeq(i,j),j=1,norbprm_t)
  enddo
  write(LUCOR,*)'Correlation matrix for GRACE '//satnam(1)
  do i = 1, norbprm_t
     write(LUCOR,'(i4,12e15.4)')i,(corr(i,j,1),j = 1, norbprm_t)
  enddo
  ! write the correlations between tidal amplitude parameters
  !    if(imsctide > 0 )then
  !    write(LUCOR,*)'Correlation matrix for tidal mascon amplitudes '
  !      do i=1,nmsc_tid_constit_t*2
  !        write(LUCOR,'(i4,20e15.4)')i,(corr(i,j,1),j = 1, nmsc_tid_constit_t*2)
  !      enddo
  !    endif
  ! write the whole damn thing
  write(LUCOR,*)'Correlation matrix'
  do i = 1, nparam_t
     write(LUCOR,*)i,(corr(i,j,1),j = 1, nparam_t)
  enddo

  write(LUCOR,*)'Correlation matrix for mascons'          ,nparam_t,nmascons_t,nparam_t-nmascons_t+1
  do i = nparam_t-nmascons_t+1,nparam_t
     !      print*,i,((nparam_t-nmascons_t+1+j-1),j=1,nmascons_t)
     write(LUCOR,*)i,(corr(i,nparam_t-nmascons_t+1+j-1,1),j = 1, nmascons_t)
  enddo

  !****************************************************************

  close (LUCOR)

  return
end subroutine output_writeCOR


!********************************************************************************************************************************
! output_writeNORM: Writes normal equations and RHS to a binary file
! Author: Paul Tregoning
!         28 October 2013
!********************************************************************************************************************************

subroutine output_writeNORM(rmsnam,apr_prm,normeq,AtWb, apr_wght, mcon_tides,prmnam,prmnam_size)

  ! MODS
  ! PT140709: add the number of mascon tidal amplitudes estimated to the header line (nmsc_tid_constit_t)
  ! PT170611: include mascons_mod for the number of ocean mascons
  ! PT190222: write out the observation uncertainties into the normal equation file

  use mascon_mod
  use gracefit_mod
  implicit none

        ! provides starting epoch information and duration of integration
  include 'output.h90'       ! provides unit number for binary normeq output file

  !********************  Variable declarations ********************

  character(*)    , intent(inout) :: rmsnam                  ! root name from which to form the normeq file name

  double precision, intent(in) :: apr_prm(nparam_t)          ! a priori parameter values
  double precision, intent(in) :: normeq(nparam_t,nparam_t)  ! Normal equations
  double precision, intent(in) :: AtWb(nparam_t)             ! RHS
  real(kind=8)    , intent(in) :: apr_wght(nobs_t)           ! observation uncertainties                 
  integer*4       , intent(in) :: mcon_tides(4556,2)         ! bit-mapped tidal amplitudes per mascon
  integer*4        ,intent(in)    :: prmnam_size
  character*30,     intent(inout) :: prmnam(prmnam_size)

  real(kind=8 )    :: version                                ! version number of binary file
  integer*4        :: i,j,irec                               ! Counters
  integer*4        :: rec_len                                ! record length of binary file (nparam_t*8)
  character(64)    :: normnam                                ! name of normal equations output file (binary)
  integer*4        :: iwhere          ! Counter variable uset to indicate '.' in rms file in command line argument
  integer*4        :: index           ! Function that returns the index of a certain character in a string
  integer*4        :: last_nonblank           ! Function that returns the length of a string
  integer*4        :: date(5)         ! Starting date of the orbit integration
  double precision :: sec             ! Seconds part of date 
  logical          ::  IC_solve, msc_solve, tmsc_solve  ! logicals that describe for subroutine codes2labels the type of solution that was estimated
  integer*4        :: iflag
  integer*4, allocatable  :: param_type(:)
  character*1 :: sat1, Sat2

  allocate (param_type(nparam_t))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! define satellite names based on mission
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(mission == 0)then
  sat1 = "A"
  sat2 = "B"
else if (mission == 1)then
  sat1 = "C"
  sat2 = "D"
else
  sat1 = "X"
  sat2 = "Y"
endif

!****************************************************************

! define the record length. Since we write one row of the normal equations at a time, this will be the number of parameters * 8 bytes for each number.
  rec_len = nparam_t*8

!  Create file names for the normeq file.
  iwhere = index(rmsnam,'.')
!   cover case where user leaves off "."
  if (iwhere == 0 ) then
     iwhere = last_nonblank(rmsnam) + 1
     rmsnam(iwhere:iwhere) = "."
     rmsnam(iwhere+1:iwhere+3) = "rms"
  endif
  normnam   = rmsnam(1:iwhere)//"norm"

! define the starting epoch in terms of yr/mm/dd/hr/mn/sec
  call gsec_to_ymdhms(starting_epoch, date, sec)

!********************** OPEN THE BINARY FILE ********************
  call output_openFile(LUNORM,normnam,'unknown','GRACESIM','FATAL','Error opening normal equation file'  &
       ,'direct',rec_len)

!*************** WRITE HEADER INO TO BINARY FILE ****************
!
  call status_update('STATUS','GRACEFIT','writeNORM',normnam,'Writing binary normal equations file',0)
! version number, number of parameters, and number of mascons
  version = 1.0
  irec=1
! PT140715: added the number of IC parameters to this list
! PT190222: include observation uncertainties in this list
  write(LUNORM,rec=irec)version,nparam_t,nmascons_t,nparam_t-nmascons_t-nmsc_tid_constit_t*2,nmsc_tid_constit_t*2,date,int(sec) &
       ,spanhrs !,(apr_wght(i),i=1,nobs_t)
! a priori parameter values
  irec=2
  write(LUNORM,rec=irec)(apr_prm(i),i=1,nparam_t)
! convert the parameter codes into character labels to be written out in gracefit vcv format
  IC_solve   = .false.
  msc_solve  = .false.
  tmsc_solve = .false.
  if(nsatprm_t > 0  ) IC_solve   = .true.
  if(nmascons_t > 0 ) msc_solve  = .true.
  if(nmsc_tid_constit_t > 0) tmsc_solve = .true.

  iflag = 1    ! flag to convert from ascii labels to parameter codes for output file
  ! PT210217: added the required offset of (nparam_t-nmascons_t-nmsc_tid_constit_t*2) ICs to put mascon codes after the ICs
  call output_codes2labels(iflag,1,IC_solve,msc_solve,tmsc_solve,nparam_t,nmascons_t,nmsc_tid_constit_t,nsatprm_t*2 &
       ,param_type,sat1,sat2,1,nparam_t-nmascons_t-nmsc_tid_constit_t*2,prmnam)
! Record 3: coded info related to parameter types. 
  irec=3
  write(LUNORM,rec=irec)sat1,param_type(1:nsatprm_t),sat2,param_type(nsatprm_t+1:nsatprm_t*2) &
       ,param_type(nsatprm_t*2 + 1:nparam_t)


! PT140610: add another line here, which is a series of logicals that indicate whether each mascon is an ocean ("T") or land ("F") mascon.
!  irec = 4
!  if(nmascons_t > 0 )then
!    write(LUNORM,rec=irec)(mcon_ocean(i),i=1,nmascons_t)  ! mocn_ocean is not defined in this subroutine !!!!
!  endif

! PT140709: also write out the bit-mapped mascon tidal amplitudes that were estimated. We need this if we are to stack these parameters in addnorm
  irec = 5
! PT170611: we write out the number of ocean mascons, not nmascons_t
  write(LUNORM,rec=irec)(mcon_tides(i,1),i=1,total_ocean_prim_ampl)

! For version 1.0 the normal equations commence on record 4
  ! PT140610: this is now going to be record 5, since I've added the mcon_ocean information into record 4.
  ! PT140709: now it is record 6,               since I've added the bit-mapped tidal amplitude information
  !************* WRITE normal equations to normeq FILE ***********
  do i= 1, nparam_t
     irec=irec+1
     write(LUNORM,rec=irec)(normeq(i,j),j=1,nparam_t)
! DEBUG:
!print*,(normeq(i,j),j=1,nparam_t)
  enddo
  !****************************************************************

  !******************* WRITE RHS to normeq FILE *******************
  irec=irec+1
  write(LUNORM,rec=irec)(AtWb(i),i=1,nparam_t)
!DEBUG:
!print*,"AtB",(AtWb(i),i=1,nparam_t)
  !****************************************************************
  close (LUNORM)

  return
end subroutine output_writeNORM

!********************************************************************************************************************************
! output_codes2labels: converts parameter codes (1-17) into parameter labels for output files
! Author: Paul Tregoning
!         31 October 2013
!********************************************************************************************************************************

subroutine output_codes2labels(iflag,nfile,IC_solve,msc_solve,tmsc_solve,nparam_t,nmascon_t,nmsc_tidal_t,nIC_t &
     ,param_type,sat1,sat2,iIC,msc_offset,prmnam)

  implicit none

  !********************  Variable declarations ********************

  integer*4,    intent(in) :: iflag                                ! label to code (1) or code to label (2)
  integer*4,    intent(in) :: nfile                                ! the number  in the list of the binary norm file (to identify the IC)
  logical     , intent(in) :: IC_solve,msc_solve,tmsc_solve        ! tells us what type of estimation is being done
  integer*4   , intent(in) :: nparam_t,nmascon_t,nmsc_tidal_t,nIC_t! number of parameters, mascons, tidal mascon amplitudes and ICs
  character*1 , intent(in) :: sat1,sat2                            ! satellite labels for the parameters
  integer*4   , intent(in) :: iIC                                  ! pointer to the start of the next ICs into the prmnam array
  integer*4   , intent(in) :: msc_offset                           ! offset possibly required to shift mascons from after to before the ICs (or zero entered otherwise)
  integer*4   , intent(inout) :: param_type(*)                     ! parameter codes from binary norm file
  character*30, intent(inout):: prmnam(nparam_t)                   ! character labels for each parameter

  integer*4    :: i
  character*14 :: label
  character*1  :: tmpsat

  !****************************************************************

  if(iflag == 1) then         ! label to code (for gracefit/gracesim)
     do i=1,nparam_t
        if(prmnam(i)(13:16) == "X0  ")param_type(i) = 1
        if(prmnam(i)(13:16) == "Y0  ")param_type(i) = 2
        if(prmnam(i)(13:16) == "Z0  ")param_type(i) = 3
        if(prmnam(i)(13:16) == "XV0 ")param_type(i) = 4
        if(prmnam(i)(13:16) == "YV0 ")param_type(i) = 5
        if(prmnam(i)(13:16) == "ZV0 ")param_type(i) = 6
        if(prmnam(i)(13:16) == "sclx")param_type(i) = 7
        if(prmnam(i)(13:16) == "scly")param_type(i) = 8
        if(prmnam(i)(13:16) == "sclz")param_type(i) = 9
        if(prmnam(i)(13:16) == "bsx ")param_type(i) = 10
        if(prmnam(i)(13:16) == "bsy ")param_type(i) = 11
        if(prmnam(i)(13:16) == "bsz ")param_type(i) = 12
        if(prmnam(i)(13:16) == "GanX")param_type(i) = 13
        if(prmnam(i)(13:16) == "GanY")param_type(i) = 14
        if(prmnam(i)(13:16) == "GanZ")param_type(i) = 15

        ! static or time-varying mascon (they don't exist yet but some day we may have mascons estimated more than once within an orbit integration)
        if(prmnam(i)(7:9) == " MC")param_type(i) = 16     ! constant mascon
        if(prmnam(i)(8:10) == "TMC")param_type(i) = 17    ! tidal amplitude  mascon

        ! PT140821: once-per-rev along-track accelerations
        if(prmnam(i)(13:16) == "1prS")param_type(i) = 18
        if(prmnam(i)(13:16) == "1prC")param_type(i) = 19

        ! PT140821: twice-per-rev along-track accelerations
        if(prmnam(i)(13:16) == "2prS")param_type(i) = 20
        if(prmnam(i)(13:16) == "2prC")param_type(i) = 21

        ! PT170420: rpy parameters
        if(prmnam(i)(13:16) == "rpyX")param_type(i) = 22
        if(prmnam(i)(13:16) == "rpyY")param_type(i) = 23
        if(prmnam(i)(13:16) == "rpyZ")param_type(i) = 24

     enddo

  else if (iflag == 2 )then   ! code to label (for addnorm)
     if(nIC_t > 0 .and. (IC_solve) ) then       ! we need IC labels
        tmpsat = " "
        do i=1,nparam_t
           label = " "
           if(param_type(i) == 1 )label = " X0   (mm)    "
           if(param_type(i) == 2 )label = " Y0   (mm)    "
           if(param_type(i) == 3 )label = " Z0   (mm)    "
           if(param_type(i) == 4 )label = " XV0  (mm/s)  "
           if(param_type(i) == 5 )label = " YV0  (mm/s)  "
           if(param_type(i) == 6 )label = " ZV0  (mm/s)  "
           if(param_type(i) == 7 )label = " sclx (n/a)   "
           if(param_type(i) == 8 )label = " scly (n/a)   "
           if(param_type(i) == 9 )label = " sclz (n/a)   "
           if(param_type(i) == 10)label = " bsx  (nm/s^2)"
           if(param_type(i) == 11)label = " bsy  (nm/s^2)"
           if(param_type(i) == 12)label = " bsz  (nm/s^2)"
           if(param_type(i) == 13)label = " GanX (mm)    "
           if(param_type(i) == 14)label = " GanY (mm)    "
           if(param_type(i) == 15)label = " GanZ (mm)    "
           if(param_type(i) == 18)label = " 1prS (um/s^2)"
           if(param_type(i) == 19)label = " 1prC (um/s^2)"
           if(param_type(i) == 20)label = " 2prS (um/s^2)"
           if(param_type(i) == 21)label = " 2prC (um/s^2)"
           if(param_type(i) == 22)label = " rpyX (mrad)  "
           if(param_type(i) == 23)label = " rpyY (mrad)  "
           if(param_type(i) == 24)label = " rpyZ (mrad)  "

           ! PT210217: this code assumes that we have a 12-parameter satellite model. It will fail if this is not the case!
           if(label(1:2) /= " ")then
             ! set the satellite name. Assume that it is ALWAYS sat1 first and then sat2 in the list of ICs
             if( (tmpsat == " " .or. tmpsat == sat2) .and. label == " X0   (mm)    ")then
               tmpsat = sat1
             else if (tmpsat == sat1 .and. label == " X0   (mm)    ")then
               tmpsat = sat2
             endif
             write(prmnam(iIC+i),'(i3,i3,1x,a3,1x,a1,a1,1x,a14)')nfile,param_type(i),"SAT",tmpsat,":",label


           endif
        enddo
     endif

     if (nmascon_t > 0 ) then                                  ! we need mascon labels
        do i=1,nparam_t
           if(param_type(i) == 16)then
              write(prmnam(i+msc_offset),'(i6,a,i5.5,a)')i,". MC",i+msc_offset," (m)     "
           else if (param_type(nIC_t+i) == 17)then
              ! PT140711: do nothing - assign these labels back in addnorm where we know something about which mascon and which tidal element
              !          write(prmnam(i),'(i6,a24)')i," tidal ampl mascon (m)  "
           endif
        enddo
     endif

  else
     call status_update('FATAL','GRACEFIT','codes2labels',' ','Wrong direction flag (must be 1 or 2). Check calling arguments',0)
  endif

  return
end subroutine output_codes2labels


!********************************************************************************************************************************
! output_writeAmat: Writes A matrix to a binary file
! Author: Rebecca McGirr
!         03 February 2021
!********************************************************************************************************************************

subroutine output_writeAmat(rmsnam,Amat)

  use mascon_mod
  use gracefit_mod
  implicit none

        ! provides starting epoch information and duration of integration
  include 'output.h90'       ! provides unit number for binary normeq output file

  !********************  Variable declarations ********************
  character(*)    , intent(inout) :: rmsnam                           ! root name from which to form the Amat file name
  double precision, intent(in)    :: Amat(nepochs_t*nobs_t,nparam_t)  ! partial equations

  real(kind=8 )    :: version                                ! version number of binary file
  integer*4        :: i,j,irec                               ! Counters
  integer*4        :: rec_len                                ! record length of binary file (nparam_t*8)
  character(64)    :: amatnam                                ! name of normal equations output file (binary)
  integer*4        :: iwhere          ! Counter variable uset to indicate '.' in rms file in command line argument
  integer*4        :: index           ! Function that returns the index of a certain character in a string
  integer*4        :: last_nonblank           ! Function that returns the length of a string
  integer*4        :: date(5)         ! Starting date of the orbit integration
  double precision :: sec             ! Seconds part of date 

!****************************************************************

! define the record length. Since we write one row of the partial equations at a time, this will be the number of observations * epochs * 8 bytes for each number.
  rec_len = nparam_t*8

!  Create file names for the Amat file.
  iwhere = index(rmsnam,'.')
!   cover case where user leaves off "."
  if (iwhere == 0 ) then
     iwhere = last_nonblank(rmsnam) + 1
     rmsnam(iwhere:iwhere) = "."
     rmsnam(iwhere+1:iwhere+3) = "rms"
  endif
  amatnam   = rmsnam(1:iwhere)//"Amat"

! define the starting epoch in terms of yr/mm/dd/hr/mn/sec
  call gsec_to_ymdhms(starting_epoch, date, sec)

!********************** OPEN THE BINARY FILE ********************
  call output_openFile(LUAMAT,amatnam,'unknown','GRACESIM','FATAL','Error opening Amat file'  &
       ,'direct',rec_len)

!*************** WRITE HEADER INFO TO BINARY FILE ****************
!
  call status_update('STATUS','GRACEFIT','writeAmat',amatnam,'Writing binary Amat file',0)
! version number, number of observations, number of parameters
  version = 1.0
  irec=1
  write(LUAMAT,rec=irec)version,nobs_t,nepochs_t,nparam_t,nmascons_t,date,int(sec),spanhrs

!*************** WRITE A MATRIX TO BINARY FILE ****************
  do i= 1, nobs_t*nepochs_t
     irec=irec+1
     write(LUAMAT,rec=irec)(Amat(i,j),j=1,nparam_t)
! DEBUG:
!print*,(Amat(i,j),j=1,nparam_t)
  enddo

  close (LUAMAT)

  return
end subroutine output_writeAmat

!********************************************************************************************************************************
! output_writebmat: Writes b matrix to a binary file
! Author: Rebecca McGirr
!         03 February 2021
!********************************************************************************************************************************

subroutine output_writebvec(rmsnam,bvec)

  use mascon_mod
  use gracefit_mod
  implicit none

        ! provides starting epoch information and duration of integration
  include 'output.h90'       ! provides unit number for binary normeq output file

  !********************  Variable declarations ********************
  character(*)    , intent(inout) :: rmsnam                           ! root name from which to form the Amat file name
  double precision, intent(in)    :: bvec(nepochs_t*nobs_t)  ! partial equations

  real(kind=8 )    :: version                                ! version number of binary file
  integer*4        :: i,irec                                 ! Counters
  integer*4        :: rec_len                                ! record length of binary file (nparam_t*8)
  character(64)    :: bvecnam                                ! name of normal equations output file (binary)
  integer*4        :: iwhere          ! Counter variable uset to indicate '.' in rms file in command line argument
  integer*4        :: index           ! Function that returns the index of a certain character in a string
  integer*4        :: last_nonblank           ! Function that returns the length of a string
  integer*4        :: date(5)         ! Starting date of the orbit integration
  double precision :: sec             ! Seconds part of date 

!****************************************************************

! define the record length. Since we write one row of the partial equations at a time, this will be the number of  4 * 8 bytes for each number.
  rec_len = 4*8

!  Create file names for the Amat file.
  iwhere = index(rmsnam,'.')
!   cover case where user leaves off "."
  if (iwhere == 0 ) then
     iwhere = last_nonblank(rmsnam) + 1
     rmsnam(iwhere:iwhere) = "."
     rmsnam(iwhere+1:iwhere+3) = "rms"
  endif
  bvecnam   = rmsnam(1:iwhere)//"bvec"

! define the starting epoch in terms of yr/mm/dd/hr/mn/sec
  call gsec_to_ymdhms(starting_epoch, date, sec)

!********************** OPEN THE BINARY FILE ********************
  call output_openFile(LUBVEC,bvecnam,'unknown','GRACESIM','FATAL','Error opening bvec file'  &
       ,'direct',rec_len)

!*************** WRITE HEADER INFO TO BINARY FILE ****************
!
  call status_update('STATUS','GRACEFIT','writebvec',bvecnam,'Writing binary bvec file',0)
! version number, number of observations, number of parameters
  version = 1.0
  irec=1
  write(LUBVEC,rec=irec)version,nobs_t,nepochs_t
! split head into to two lines so rec_len can be shorter
  irec=2
  write(LUBVEC,rec=irec)date,int(sec),spanhrs
!*************** WRITE A MATRIX TO BINARY FILE ****************
  do i= 1, nobs_t*nepochs_t
     irec=irec+1
     write(LUBVEC,rec=irec)bvec(i)
! DEBUG:
!print*,(bvec(i,j),j=1,nparam_t)
  enddo

  close (LUBVEC)

end subroutine output_writebvec



!********************************************************************************************************************************
