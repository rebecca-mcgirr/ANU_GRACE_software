!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! various subroutines to read L1B headers and data ....
! 
! level1B_open    :  opens a L1B file
! sca_read_hdr    :  reads the header of a SCA1B file (GRACE or GRACE FO)
! sca_read_data   :  reads the quaternions and camera flag of SCA1B data (GRACE or GRACE FO)
! thr_read_hdr    :  reads the header of a SCA1B file (GRACE or GRACE FO)
! thr_read_data   :  reads the thrust timing, duration and thruster data from THR1B file (GRACE or GRACE FO)
! gnv_read_hdr    :  reads the header of a GNV1B file (GRACE or GRACE FO)
! gnv_read_data   :  reads the position/velocity data from the GNV1B file (GRACE or GRACE FO)
! kbr_read_hdr    :  reads the header of a KBR1B file (GRACE or GRACE FO)
! kbr_read_data   :  reads the KBR, KBRR, KBRA data from KBR1B file (GRACE or GRACE FO)
! lri_read_hdr    :  reads the header of a LRI1B file (GRACE FO)
! lri_read_data   :  reads the LRR, LRRR, LRRA data from LRR1B file (GRACE FO)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine level1B_open(luin,calling_prog,level1B_file)

! open the accelerometer file
  
  implicit none

! passed variables
  character(*),intent(in)   :: level1B_file     ! name of accelerometer file
  character(*),intent(in)   :: calling_prog ! name of calling program
  integer*4   ,intent(in)   :: luin         ! unit number of file

! local variables
  integer*4     :: ioerr

! P. Tregoning
! 31 May 2018

! open the file
  open(luin,file=level1B_file,status='old',iostat=ioerr)

! output information regarding opening the file
  if(ioerr == 0)then
    call status_update('STATUS',calling_prog,'level1B_open',level1B_file,'Have opened Level-1B file',0)
  else
    call status_update('FATAL',calling_prog,'level1B_open',level1B_file,'Error opening Level-1B file',ioerr)
  endif

  return
  end subroutine level1B_open
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine gnv_read_data(luin,calling_prog,gnv_file,mission,maxgprm,nsat_t,nepochs_t,isat,starting_epoch,epoch_interval,n_gnv &
                           ,gnv_step,gnv_epochs,max_gnv,gnv_obs )

! subroutine to read positions/velocities from a GNV1B file. The data are actually the same format for GRACE and GRACE FO.
!
! P. Tregoning
! 28 June 2018

  implicit none

! passed variables
  integer*4   ,intent(in)    :: luin                                 ! unit number of file
  character(*),intent(in)    :: calling_prog                         ! name of calling program
  character(*),intent(in)    :: gnv_file                             ! name of star camera file
  integer*4   ,intent(inout) :: mission                              ! 0: GRACE, 1: GRACE FO, 2: GRACE II  
  integer*4   , intent(in)   :: maxgprm,nsat_t,nepochs_t             ! dimensions for the pos/vel array
  integer*4   ,intent(in)    :: isat                                 ! number of actual satellite data to be read and stored
  real(kind=8),intent(in)    :: starting_epoch                       ! graceseconds of the starting epoch of the range of data we want
  real(kind=8),intent(in)    :: epoch_interval                       ! interval of the data to be used and analysed in GRACEFIT/SIM
  integer*4   ,intent(in)    :: n_gnv                                ! number of gnv observations in the GNV1B file
  real(kind=8),intent(in)    :: gnv_step                             ! #seconds between each GNV1B observation
  integer*4,   intent(in)    :: max_gnv                              ! dimension of GNV1B epoch array
  integer*4   ,intent(out)   :: gnv_epochs(max_gnv)                  ! array of epochs for which there are GNV1B obs
  real(kind=8),intent(out)   :: gnv_obs(maxgprm,nsat_t,nepochs_t)    ! epochs, quaternions, camera flag

! local variables
  integer*4     :: gracesec,iepoch,i,ioerr,counter
  character*1   :: sat,frame
  character*250 :: line,message
  real(kind=8)  :: pos(3),pos_sig(3),vel(3),vel_sig(3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the position/velocity data. The header should already have been read, so we just need to read in the data.
! 
! gracesec sat frame pos(3) pos_sig(3) vel(3) vel_sig(3) 8-bits_of_flags e.g:
! 271868390 D E -1346567.5776538 86679.84337818051 -6706961.24194364 0 0 0 7394.27443317575 -1048.34110525824 -1494.87440176657 0 0 0  00000000
  ioerr = 0
  counter = 0
  do i=1,n_gnv
    read(luin,*,iostat=ioerr,end=1000)gracesec,sat,frame,pos,pos_sig,vel,vel_sig
! determine where to fit the epoch into the gvec array (which gos from 1:nepochs_t)
    iepoch = (gracesec - nint(starting_epoch))/gnv_step + 0
!PT200624: change logic so that we store only 5 second data for both GRACE and GRACE-FO
    if(iepoch >= 0 .and. counter < nepochs_t  )then
      if( (gnv_step < 5.d0 .and. mod(dble(iepoch),epoch_interval) == 0.d0) .or. gnv_step >= 5.d0)then
        counter = counter + 1
        gnv_obs(1:3,isat,counter) = pos
        gnv_obs(4:6,isat,counter) = vel
        gnv_epochs(counter) = gracesec
      endif
    endif 
  enddo

1000 continue

  return
  end subroutine gnv_read_data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine gnv_read_hdr(luin,calling_prog,gnv_file,mission,n_gnv,gnv_step,gnv_span )

! subroutine to read through the header information of an SCA1B file
!
! P. Tregoning
! 19 June 2018

  implicit none

! passed variables
  integer*4   ,intent(in)    :: luin         ! unit number of file
  character(*),intent(in)    :: calling_prog ! name of calling program
  character(*),intent(in)    :: gnv_file     ! name of GNV1B file
  integer*4,  intent(inout)  :: mission      ! 0: GRACE, 1: GRACE FO, 2: GRACE II
  integer*4,  intent(out)    :: n_gnv        ! number of star camera records in file
  real(kind=8), intent(out)  :: gnv_step     ! interval between SCA1B observations (5s for GRACE, 1s for GRACE FO)

! local variables
  real(kind=8)  :: seconds
  integer*4     :: gnv_start,gnv_end,gnv_span
  integer*4     :: ioerr
  character*250 :: line,message

! check whether the mission variable has been set
  if(mission < 0) then

  ! not set. Determine the mission (GRACE, GRACE FO, GRACE II) from the first line of the header
    read(luin,'(a)')line
    if(line(1:7) == "PRODUCE")then
      mission = 0    ! it is a GRACE star camera file
      call status_update('STATUS',calling_prog,'gnv_read_hdr',' ','GRACE format found ',0)
    elseif(line(1:7) == "header:")then
      mission = 1    ! it is a GRACE FO star camera file
      call status_update('STATUS',calling_prog,'gnv_read_hdr',' ','GRACE FO format found',0)
    endif
    backspace(luin)
  endif

! now, read the header information.
  line = ""
  ioerr = 0
  if(mission == 0)then
    do while (line(1:13) /= "END OF HEADER" .and. ioerr == 0)
      read(luin,'(a)',iostat=ioerr)line
      if(line(1:31) == "NUMBER OF DATA RECORDS        :")read(line(32:41),*)n_gnv
      if(line(1:14) == "TIME FIRST OBS")then
         read(line(32:42),*) seconds        ! HM190608 read seconds as real in case it is <10^8  
         gnv_start=int(seconds)
         read(luin,'(a)',iostat=ioerr) line
         read(line(32:42),*) seconds
         gnv_end=int(seconds)
         gnv_span = gnv_end-gnv_start+1
      endif
      gnv_step = 5  ! 5 sec data for GRACE
    enddo

  else if (mission == 1)then
    do while (line(1:20) /= "# End of YAML header" .and. ioerr == 0)
      read(luin,'(a)',iostat=ioerr)line
      if(line(1:13) == "  dimensions:")read(luin,'(17x,i7)')n_gnv
      if(line(1:17) == "    summary: 1-Hz")then
        call status_update('STATUS',calling_prog,'gnv_read_hdr',gnv_file,line(1:50),0)
        gnv_step = 1  ! 1Hz data for GRACE FO (at least, the simulated data is 1 s)
      endif

! PT180730: get the start and end times of the data to know the data span
      if(line(1:26) == "    start_time_epoch_secs:")then
         read(line(27:37),*)gnv_start
         read(luin,'(25x,i11)')gnv_end
         gnv_span = gnv_end-gnv_start+1
      endif

    enddo    

  endif

  write(message,'(i8,a,i8,a)')n_gnv," pos/vel records found in GNV1B file. Data span of ",gnv_span,' seconds'
  call status_update('STATUS',calling_prog,'gnv_read_hdr',gnv_file,message,0)

  return
  end subroutine gnv_read_hdr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sca_read_hdr(luin,calling_prog,sca_file,mission,n_sca,sca_step )

! subroutine to read through the header information of an SCA1B file
!
! P. Tregoning
! 19 June 2018

  implicit none

! passed variables
  integer*4   ,intent(in)    :: luin         ! unit number of file
  character(*),intent(in)    :: calling_prog ! name of calling program
  character(*),intent(in)    :: sca_file     ! name of star camera file
  integer*4,  intent(inout)  :: mission      ! 0: GRACE, 1: GRACE FO, 2: GRACE II
  integer*4,  intent(out)    :: n_sca        ! number of star camera records in file
  real(kind=8), intent(out)  :: sca_step     ! interval between SCA1B observations (5s for GRACE, 1s for GRACE FO)

! local variables
  integer*4     :: ioerr
  character*250 :: line,message
  integer*4     :: i,sca_epoch

! check whether the mission variable has been set
  if(mission < 0) then

  ! not set. Determine the mission (GRACE, GRACE FO, GRACE II) from the first line of the header
    read(luin,'(a)')line
    if(line(1:7) == "PRODUCE")then
      mission = 0    ! it is a GRACE star camera file
      call status_update('STATUS',calling_prog,'sca_read_hdr',' ','GRACE format found ',0)
    elseif(line(1:7) == "header:")then
      mission = 1    ! it is a GRACE FO star camera file
      call status_update('STATUS',calling_prog,'sca_read_hdr',' ','GRACE FO format found',0)
    endif
    backspace(luin)
  endif

! now, read the header information.
  line = ""
  ioerr = 0
  if(mission == 0)then
    do while (line(1:13) /= "END OF HEADER" .and. ioerr == 0)
      read(luin,'(a)',iostat=ioerr)line
      if(line(1:31) == "NUMBER OF DATA RECORDS        :")read(line(32:41),*)n_sca
      sca_step = 5  ! 5 sec data for GRACE
    enddo

! SA190812 Read the number of data from the file. not the header... I used the same structure of what Paul did for the RL04 data.
! I am sure there is a much better way to do it... 
    do i=1,n_sca
       read(luin,*,end=1000)sca_epoch
    enddo
    rewind(luin)
    line = " "
    do while (line(1:13) /= "END OF HEADER" .and. ioerr == 0)
       read(luin,'(a)',iostat=ioerr)line
    enddo


  else if (mission == 1)then
    do while (line(1:20) /= "# End of YAML header" .and. ioerr == 0)
      read(luin,'(a)',iostat=ioerr)line
      if(line(1:13) == "  dimensions:")then
        read(luin,'(17x,i7)')n_sca
      endif
      if(line(1:36) == "    summary: 1-Hz processed SCA data")then
        call status_update('STATUS',calling_prog,'sca_read_hdr',sca_file,"   1-Hz processed SCA1B data",0)
        sca_step = 1  ! 1Hz data for GRACE FO (at least, the simulated data is 1 s)
      endif
    enddo    

! PT190525: one cannot rely on the header telling us how many records are in the file. 2018-12-30 file
!           says 87000 records but there are only 86400 - and the data starts at 00UT but the header
!           says 23:55 of the day before. So, we are going to have to read all the records beyond the
!           header and confirm that the value of n_sca is correct!
    do i=1,n_sca
      read(luin,*,end=1000)sca_epoch  
    enddo
    rewind(luin)
    line = " "
    do while (line(1:20) /= "# End of YAML header" .and. ioerr == 0)
      read(luin,'(a)',iostat=ioerr)line
    enddo

  endif

  write(message,'(i8,a)')n_sca," star camera records found in SCA1B file"
  call status_update('STATUS',calling_prog,'sca_read_hdr',sca_file,message,0)

  return

1000 continue
! if we get to here it is because the yaml header said that there were more obs than there actually are!!!
  write(message,'(a,i8,a,i8)')"yaml header is wrong. There are not", n_sca &
                          ," star camera records found in SCA1B file. There are only",i-1
  call status_update('STATUS',calling_prog,'sca_read_hdr',sca_file,message,0)
  n_sca = i - 1
  rewind(luin)
  line = " "
  ! SA This line should work, don't know why
  !do while ((line(1:20) /= "# End of YAML header" .or. line(1:13) /= "END OF HEADER" )  .and. ioerr == 0)
  do while (line(1:20) /= "# End of YAML header"   .and. ioerr == 0)
    read(luin,'(a)',iostat=ioerr)line
    if (line(1:13) == "END OF HEADER" ) ioerr =5
  enddo

  end subroutine sca_read_hdr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sca_read_data(luin,calling_prog,sca_file,mission,n_sca,sca_obs )

! subroutine to read quaternions in a SCA1B file. The data are actually the same format for GRACE and GRACE FO.
!
! P. Tregoning
! 19 June 2018

  implicit none

! passed variables
  integer*4   ,intent(in)    :: luin                ! unit number of file
  character(*),intent(in)    :: calling_prog        ! name of calling program
  character(*),intent(in)    :: sca_file            ! name of star camera file
  integer*4   ,intent(inout) :: mission             ! 0: GRACE, 1: GRACE FO, 2: GRACE II  
  integer*4   ,intent(in)    :: n_sca               ! dimensioning of rows in sca_obs array
  real(kind=8),intent(out)   :: sca_obs(n_sca,6)    ! epochs, quaternions, camera flag

! local variables
  integer*4     :: ioerr,i,j,gracesec
  character*1   :: sat
  character*250 :: line,message
  real(kind=8)  :: tmp_obs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the star cambera data. The header should already have been read, 
! so we just need to read in the data.
! 
! gracesec camera_flag 4xquaternion_values
  ioerr = 0
  gracesec = 0
  do i=1,n_sca
    read(luin,*,iostat=ioerr,end=1000)sca_obs(i,1),sat,sca_obs(i,6),(sca_obs(i,j),j=2,5)
! PT180619: put a fatal catch if there is a gap in the quaternion data
    if (i > 1 ) then
      if ((sca_obs(i,1) - sca_obs(i-1,1)) > 5.d0)then
        write(message,'(a,i10,a,i10,a,i5,a)')"Missing quaternion data from",int(sca_obs(i-1,1))," to",int(sca_obs(i,1)) &
                                       ,". ",int(sca_obs(i,1)-sca_obs(i-1,1))/5,' epochs missing'
        call status_update('WARNING',calling_prog,'sca_read_data',sca_file,message,0)
      endif
    endif
  enddo
  return

1000 continue
  write(message,'(a,i8,a,i8)')"Error: reached end of SCA file after ",i," records instead of ",n_sca
  call status_update('FATAL',calling_prog,'sca_read_data',sca_file,message,0)

  return
  end subroutine sca_read_data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine thr_read_hdr(luin,calling_prog,thr_file,mission,n_thrusts )

! subroutine to read through the header information of an THR1B file
!
! P. Tregoning
! 31 May 2018

  implicit none

! passed variables
  integer*4   ,intent(in)    :: luin         ! unit number of file
  character(*),intent(in)    :: calling_prog ! name of calling program
  character(*),intent(in)    :: thr_file     ! name of accelerometer file
  integer*4,  intent(inout)  :: mission      ! 0: GRACE, 1: GRACE FO, 2: GRACE II
  integer*4,  intent(out)    :: n_thrusts    ! number of thrusts found in THR1B file

! local variables
  integer*4     :: ioerr
  character*250 :: line,message

! check whether the mission variable has been set
  if(mission < 0) then

  ! not set. Determine the mission (GRACE, GRACE FO, GRACE II) from the first line of the header
    read(luin,'(a)')line
    if(line(1:7) == "PRODUCE")then
      mission = 0    ! it is a GRACE thrust file
      call status_update('STATUS',calling_prog,'acc_read_hdr',' ','GRACE format found ',0)
    elseif(line(1:7) == "header:")then
      mission = 1    ! it is a GRACE FO thrust file. We don't know yet what these will look like!
      call status_update('STATUS',calling_prog,'acc_read_hdr',' ','GRACE FO format found',0)
    endif
    backspace(luin)
  endif

! now, read the header information.
  line = ""
  ioerr = 0
  if(mission == 0)then
    do while (line(1:13) /= "END OF HEADER" .and. ioerr == 0)
      read(luin,'(a)',iostat=ioerr)line
      if(line(1:31) == "NUMBER OF DATA RECORDS        :")read(line(32:41),*)n_thrusts
    enddo

  else if (mission == 1)then
    do while (line(1:20) /= "# End of YAML header" .and. ioerr == 0)
      read(luin,'(a)',iostat=ioerr)line
      if(line(1:13) == "  dimensions:")read(luin,'(17x,i7)')n_thrusts
    enddo    

  endif

  write(message,'(a,i8,a)')' read',n_thrusts," thrusts obs from THR1B file"
  call status_update('STATUS',calling_prog,'thr_read_hdr',thr_file,message,0)

  return
  end subroutine thr_read_hdr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine thr_read_data(luin,calling_prog,thr_file,mission,n_thrusts,max_thrusts,thr_obs )

! subroutine to read all the timings of thrusts in a THR1B file. 
!
! P. Tregoning
! 31 May 2018

  implicit none

! passed variables
  integer*4   ,intent(in)    :: luin                  ! unit number of file
  character(*),intent(in)    :: calling_prog          ! name of calling program
  character(*),intent(in)    :: thr_file              ! name of accelerometer file
  integer*4   ,intent(inout) :: mission               ! 0: GRACE, 1: GRACE FO, 2: GRACE II  
  integer*4   ,intent(in)    :: n_thrusts             ! dimensioning of rows in thr_obs array
  integer*4   ,intent(in)    :: max_thrusts           ! the max number of thrusts of the two satellites (used to dimension the array)
  integer*4   ,intent(out)   :: thr_obs(max_thrusts,8)  ! thrust info (gracesec, fractional sec, ms of 6 different thrusters)

! local variables
  integer*4     :: ioerr,i,j,gracesec
  character*1   :: sat
  character*250 :: line,message
  real(kind=8)  :: tmp_obs
  integer*4     :: int_array(14)                      ! dummy array to read all the integer values in the thrust line
  character*2   :: char_array(2)                      ! dummy array to read the two single character variables in the thrust line
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the thruster data. The header should already have been read, 
! so we just need to read in the data
  ioerr = 0
  gracesec = 0
  do i=1,n_thrusts
    read(luin,*,iostat=ioerr,end=1000)thr_obs(i,1),thr_obs(i,2),char_array,(int_array(j),j=1,14),(thr_obs(i,j),j=3,8)
  enddo

1000 continue

  return
  end subroutine thr_read_data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kbr_read_hdr(luin,calling_prog,kbr_file,mission,n_kbr )

! subroutine to read through the header information of an THR1B file
!
! P. Tregoning
! 27 July 2018

  implicit none

! passed variables
  integer*4   ,intent(in)    :: luin         ! unit number of file
  character(*),intent(in)    :: calling_prog ! name of calling program
  character(*),intent(in)    :: kbr_file     ! name of accelerometer file
  integer*4,  intent(inout)  :: mission      ! 0: GRACE, 1: GRACE FO, 2: GRACE II
  integer*4,  intent(out)    :: n_kbr        ! number of kbr1b obs found in KBR1B file

! local variables
  integer*4     :: ioerr
  character*250 :: line,message

! check whether the mission variable has been set
  if(mission < 0) then

  ! not set. Determine the mission (GRACE, GRACE FO, GRACE II) from the first line of the header
    read(luin,'(a)')line
    if(line(1:7) == "PRODUCE")then
      mission = 0    ! it is a GRACE thrust file
      call status_update('STATUS',calling_prog,'kbr_read_hdr',' ','GRACE format found ',0)
    elseif(line(1:7) == "header:")then
      mission = 1    ! it is a GRACE FO thrust file. We don't know yet what these will look like!
      call status_update('STATUS',calling_prog,'kbr_read_hdr',' ','GRACE FO format found',0)
    endif
    backspace(luin)
  endif

! now, read the header information.
  line = ""
  ioerr = 0
  if(mission == 0)then
    do while (line(1:13) /= "END OF HEADER" .and. ioerr == 0)
      read(luin,'(a)',iostat=ioerr)line
      if(line(1:31) == "NUMBER OF DATA RECORDS        :")read(line(32:41),*)n_kbr
    enddo

  else if (mission == 1)then
    do while (line(1:20) /= "# End of YAML header" .and. ioerr == 0)
      read(luin,'(a)',iostat=ioerr)line
      if(line(1:13) == "  dimensions:")read(luin,'(17x,i7)')n_kbr
    enddo    

  endif

  write(message,'(i8,a)')n_kbr," kbr observations found in KBR1B file"
  call status_update('STATUS',calling_prog,'kbr_read_hdr',kbr_file,message,0)

  return
  end subroutine kbr_read_hdr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kbr_read_data(luin,calling_prog,kbr_file,mission,nepochs_t,starting_epoch,epoch_interval,kbrange,kblt,kbant &
                           ,num_epochs)

! subroutine to read through the header information of an KBR1B file
!
! P. Tregoning
! 19 November 2018

  implicit none

! passed variables
  integer*4   ,intent(in)    :: luin                      ! unit number of file
  character(*),intent(in)    :: calling_prog              ! name of calling program
  character(*),intent(in)    :: kbr_file                  ! name of accelerometer file
  integer*4,  intent(inout)  :: mission                   ! 0: GRACE, 1: GRACE FO, 2: GRACE II
  real(kind=8) , intent(in ) :: starting_epoch            ! number of kbr1b obs found in KBR1B file
  real(kind=8) , intent(in ) :: epoch_interval            ! as it says

  ! data variables
  integer*4    ,intent(in)   :: nepochs_t                 ! number of epochs
  real(kind=8) , intent(out) :: kbrange(3,nepochs_t)      ! Range value (biased), rate, and acceleration
  real(kind=8) , intent(out) :: kblt(3,nepochs_t)         ! Light time corrections for range, rate and acceleration
  real(kind=8) , intent(out) :: kbant(3,nepochs_t)        ! Antenna phase center corrections for range, rate and acceleration
  integer*4    , intent(out) :: num_epochs                ! number of epochs of KBR data found

! PT180916: variables to enable reading of the error flags in the KBR1B files
  integer*4 :: dummy_values(4)
  integer*4 :: error_flags

! local variables
  integer*4     :: ioerr,i,iepoch
  character*250 :: line,message
  real(kind=8)  :: dummy_var,epoch
  integer*4     :: kb_iepoch,missing_interval  ! variables related to checking for missing data epochs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT181119: declare the arrays to have a zero value
  epoch = 0
  kbrange = 0.d0
  kblt = 0.d0
  kbant = 0.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Skip the lines that are before the starting epoch
  do while(epoch < starting_epoch)
     read(luin,*,iostat=ioerr,end=101) epoch,(kbrange(i,1),i=1,3),dummy_var,&
          (kblt(i,1),i=1,3),(kbant(i,1),i=1,3)
  enddo
101 continue !  assumed to have at least one line
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT181119: find the first measurement in the file that corresponds to the
!           starting epoch to be used in gracefit/gracesim
  kb_iepoch = (epoch - starting_epoch) / epoch_interval + 1
  num_epochs = 1
  ! Index correction if necessary
  if(kb_iepoch /= 1)then  ! Put values in correct place if this isn't already the case
     if(kb_iepoch <= nepochs_t)then
        kbrange(:,kb_iepoch) = kbrange(:,1)
        kblt(:,kb_iepoch) = kblt(:,1)
        kbant(:,kb_iepoch) = kbant(:,1)
     else
        kb_iepoch = nepochs_t
        num_epochs = 0
     endif
     missing_interval = int((kb_iepoch-1)*epoch_interval)
     kbrange(:,1) = 0.d0
     kblt(:,1) = 0.d0
     kbant(:,1) = 0.d0
     write(message,'(a,i6,a,f12.1,a,f12.1,a,a,i6,a,i6,a)')'Missing data from KBR1B file:',missing_interval,&
          ' seconds (missing epochs from ', starting_epoch,' to ', starting_epoch+missing_interval-epoch_interval,' included)'&    
            ,' (epochs ',iepoch,' to ',kb_iepoch-1,')'
     call status_update('WARNING','GRACEFIT','kbr_read_data',' ',message,0)
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read in the rest of the files
  do iepoch = kb_iepoch+1, nepochs_t  ! start at kb_iepoch+1 because the first line was already read in

     if(iepoch-1 /= kb_iepoch) cycle

! PT180916: also read in the error bits and use them to eliminate bad kbr observations
     read(luin,*,iostat=ioerr,end=102) epoch,(kbrange(i,iepoch),i=1,3),dummy_var,&
          (kblt(i,iepoch),i=1,3),(kbant(i,iepoch),i=1,3),dummy_values,error_flags
     ! Assert there was no error reading file
     if(ioerr /= 0) call status_update('FATAL',calling_prog,'kbr_read_data',' ','Error reading KBR1B input file',ioerr)

       kb_iepoch = ( epoch - starting_epoch ) / epoch_interval + 1
       ! index correction if necessary
       if(kb_iepoch /= iepoch ) then !.or. error_flags /= 0)then  ! Put values in correct place if this isn't already the case
         if(kb_iepoch <= nepochs_t)then
           kbrange(:,kb_iepoch) = kbrange(:,iepoch)
           kblt(:,kb_iepoch) = kblt(:,iepoch)
           kbant(:,kb_iepoch) = kbant(:,iepoch)
         else
           kb_iepoch = nepochs_t
         endif
         missing_interval = int((kb_iepoch-iepoch)*epoch_interval)
         kbrange(:,iepoch) = 0.d0
         kblt(:,iepoch) = 0.d0
         kbant(:,iepoch) = 0.d0
         write(message,'(a,i6,a,a,f12.1,a,f12.1,a,a,i6,a,i6,a)')'Missing data from KBR1B file:',missing_interval,&
             ' seconds',' (missing epochs from ', starting_epoch+(iepoch-1)*epoch_interval,' to ',&
             starting_epoch+(iepoch-1)*epoch_interval+missing_interval-epoch_interval,' included)' &    ! PT211020: fixed bug in end time here
            ,' (epochs ',iepoch,' to ',kb_iepoch-1,')'
         call status_update('WARNING',calling_prog,'kbr_read_data',' ',message,0)
       endif

       num_epochs = num_epochs + 1
  enddo ! end of epoch loop

102 continue

! Close file, it should no longer be necessary
  close(luin)



  return
  end subroutine kbr_read_data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


