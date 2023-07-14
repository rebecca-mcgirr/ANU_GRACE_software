  program extend_L1B

! program to extend the temporal span of a L1B file to include a specified amount of data before/after the 
! data for any particular day. This involves:
!
! 1. Changing the header information
!    o  number of epochs of data
!    o  start and end time records
! 2. Adding data records before the data for the specified day
! 3. Adding data records after  the data for the specified day
!
! The program will create a new file - in L1B format - in directory "extended_L1B". The file will have the
! same name as the information inputted to the program from the command line
!
! P. Tregoning
! 9 August 2018
!
! PT180816: fixed the problem of getting the correct yYYY MM DD for preceding and subsequent days
! PT190807: modify code to work for ACT1B files as well as ACC1B. 
! PT221025: modify it again to work for ACH1B as well

  implicit none

! command line arguments
  character*3 :: L1B_type              ! ACC or KBR or LRI
  integer*4   :: year, month, day        ! requested date for which observation records will be extended
  integer*4   :: t_extra               ! amount of time (in seconds) to extend the observations before/after
  character*1 :: sat                   ! satellite name
  integer*4   :: mission               ! GRACE: 0, GRACE FO: 1
  character*10:: arg
  character*3 :: data_type             ! GRACE: "asc", GRACE FO: "txt"

! L1B file names
  character*26 :: L1Bfile(3)           ! names of L1B files
  character*39 :: L1B_extended         ! name of extended L1B file, including subdirectory name ("extended_L1B/")

! unit numbers 
  integer*4     :: luacc(4),luout

! Accelerometer variables
  integer*4     :: nvals_acc(3)        ! the number of observations in each file of accelerometer data
  integer*4     :: acc_span(3,2)       ! start/stop time (in gracesec?) of each file
  integer*4     :: iday                ! loop counter
  integer*4     :: max_obs             ! maximum number of observations in any of the three days of data
  real(kind=8), allocatable :: acc_obs(:,:,:)        ! the actual accelerometer observations for each day
  real(kind=8), allocatable :: acc_obs_extended(:,:) ! the extended set of observations

! time variables and counters
  integer*4     :: new_start,new_end   ! new start/stop time of extended file
  integer*4     :: iobs
  integer*4     :: n_obs               ! the number of observations in the extended set
  integer*4     :: ncols               ! number of columns to be read from the L1B file
  integer*4     :: date(5)             ! passed to/from gamit routines to vonvert dates
  real*8        :: jd                  ! julian day. We use this as the means to add/subtract a day
  real*8        :: sec                 ! seconds of day (not used but required by the subroutines)

! variables for making output file, header transfer etc
  logical       :: end_of_header
  character*100 :: line

! local variables
  character*15  :: calling_prog        ! pass to subroutines if necessary
  character*150 :: message
  integer*4     :: RL_num              ! release number based on mission
 

! define some values
  calling_prog = "extend_L1B"

!------------------------------------------------------------ 
! decode runstring
  call getarg(1,L1Bfile(2))
  if(L1Bfile(2)(1:1) == "")then
    call status_update('STATUS','UTIL','extend_L1B',' ',"extend_L1B L1Bfile t_extra",0)
    call status_update('FATAL','UTIL','extend_L1B',' ' &
           ,"e.g. extend_L1B ACC1B_2017-05-04_A_02.asc 7200 (where last value is extension time in secods)",0)
  endif

! extension time
  call getarg(2,arg)
  read(arg,*)t_extra
!------------------------------------------------------------ 


!------------------------------------------------------------ 
! extract information from input file name
! L1B type
  read(L1Bfile(2)(1:3),*)L1B_type
! PT190807: added type "ACT1B", the transplanted accelerometer data
  if (L1B_type /= "ACC" .and. L1B_type /= "ACT"  .and. L1B_type /= "KBR" .and. L1B_type /= "LRI" .and. L1B_type /= "ACH")then
    write(message,'(a,a,a)')"L1B type '",L1B_type,"' not coded. Please add it yourself !"
    call status_update('FATAL','UTIL','extend_L1B',' ',message,0)
  endif
 
! epoch
  read(L1Bfile(2)(7:10),*)year
  read(L1Bfile(2)(12:13),*)month  
  read(L1Bfile(2)(15:16),*)day

! satellite
  read(L1Bfile(2)(18:18),'(a1)')sat

! data release number
  read(L1Bfile(2)(20:21),*)RL_num

! mission
  read(l1Bfile(2)(23:25),'(a3)')data_type
  if(data_type == "asc")then
    mission = 0
  else if (data_type == "txt")then
    mission = 1
  else
    write(message,'(a)')"Uncoded mission: "
    call status_update('FATAL','UTIL','extend_L1B',' ',message,0)
  endif

!------------------------------------------------------------ 


!------------------------------------------------------------ 
! construct the L1B filenames for the previous and subsequent days (#&Q%(# Need to do this smarter !!!
  date = 0
  sec = 0.d0
  date(1) = year
  date(2) = month
  date(3) = day
  call ymdhms_to_jd(date,sec,jd)
  jd = jd-1.d0
  call jd_to_ymdhms(jd,date,sec)
  year = date(1)
  month = date(2)
  day = date(3)
  write(L1Bfile(1),100)L1B_type,"1B_",year,"-",month,"-",day,"_",sat,"_",RL_num,".",data_type       
  ! now the subsequent day (add 2 to the "previous day")
  jd = jd+2.d0
  call jd_to_ymdhms(jd,date,sec)
  year = date(1)
  month = date(2)
  day = date(3)
  write(L1Bfile(3),100)L1B_type,"1B_",year,"-",month,"-",day,"_",sat,"_",RL_num,".",data_type       
100 format(a3,a3,i4.4,a1,i2.2,a1,i2.2,a1,a1,a1,i2.2,a1,a3)

  write(message,'(a,a,a,i6,a)')"Will extend file ",L1Bfile(2)," by ",t_extra," seconds at either end"
  call status_update('STATUS','UTIL','extend_L1B',' ',message,0)
!------------------------------------------------------------ 


!-------------------------------------------------------------------------------------------------------------- 
! now, open and read the data for each of the three days
  do iday = 1,3
    luacc(iday) = 10 + (iday-1)
    ! open file
    call level1B_open(luacc(iday),calling_prog,L1Bfile(iday))
    ! read the header information
    call acc_read_hdr(luacc(iday),calling_prog,L1Bfile(iday),mission,acc_span(iday,:))
    nvals_acc(iday) = acc_span(iday,2) - acc_span(iday,1) + 1
  enddo

  ! allocate the array, using the maximum obs of the three days as the dimension
  max_obs = nvals_acc(1)
  if (nvals_acc(2) > max_obs) max_obs = nvals_acc(2)
  if (nvals_acc(3) > max_obs) max_obs = nvals_acc(3)
  allocate(acc_obs(3,max_obs,11))

  ! now read the actual data
  do iday = 1,3
    ! read the ACC1B linear accelerations
    call acc_read_data(luacc(iday),calling_prog,L1Bfile(iday),mission,max_obs,10,acc_span(iday,:),acc_obs(iday,:,:))
  enddo

! open the new extended file
  write(L1B_extended,'(a,a)')"extended_L1B/",L1Bfile(2)
  open(luacc(4),file=L1B_extended,status='unknown')
!-------------------------------------------------------------------------------------------------------------- 


!-------------------------------------------------------------------------------------------------------------- 
! work out the required start/stop times for the extended file
  new_start = acc_obs(2,1,1) - t_extra
  new_end   = acc_obs(2,nvals_acc(2),1) + t_extra
  write(message,'(a,i12,a,i12)')"Will try to create extended file from",new_start," to",new_end
  call status_update('STATUS','UTIL','extend_L1B',L1B_extended,message,0)
!-------------------------------------------------------------------------------------------------------------- 


!-------------------------------------------------------------------------------------------------------------- 
! loop through all the observations and save off the ones in the required time period
  if(L1B_type == "ACC" .or. L1B_type == "ACT" .or. L1B_type == "ACH")then
    ncols = 10
  else
    print*,L1B_type," not yet coded - please do it yourself !"
    stop
  endif
  allocate(acc_obs_extended(sum(nvals_acc(:)),ncols+1))  ! this will be the size of all the obs for all three days. Too big, but at least big enough!

  n_obs = 0
  do iday=1,3
    do iobs = 1,nvals_acc(iday)
      if(acc_obs(iday,iobs,1) >= new_start .and. acc_obs(iday,iobs,1) <= new_end)then
        n_obs = n_obs + 1
        acc_obs_extended(n_obs,:) = acc_obs(iday,iobs,:)
      endif
    enddo
  enddo
  write(message,'(a,i8,a)')"There are now",n_obs," L1B observation epochs in extended file"
  call status_update('STATUS','UTIL','extend_L1B',L1B_extended,message,0)
!-------------------------------------------------------------------------------------------------------------- 

!-------------------------------------------------------------------------------------------------------------- 
! transfer the header of the centre day to the output file, changing the values of #obs, start/stop time etc
  rewind(luacc(2))
  end_of_header = .false.
  do while (.not. end_of_header)
    read(luacc(2),'(a)')line

! GRACE
    if(mission == 0)then
      ! update the gracesecs of first epoch
      if(line(1:31) == "TIME FIRST OBS(SEC PAST EPOCH):")write(line(32:48),'(f17.6)')dble(new_start)  ! NOTE: this doesn't update the date on the line
      ! update the gracesecs of first epoch
      if(line(1:31) == "TIME LAST OBS(SEC PAST EPOCH) :")write(line(32:48),'(f17.6)')dble(new_end)    ! NOTE: this doesn't update the date on the line
      ! update number of observations
      if(line(1:31) == "NUMBER OF DATA RECORDS        :")write(line(32:38),'(i7)')n_obs

      ! end of header
      if(line(1:) == "END OF HEADER")end_of_header = .true.

! GRACE FO
    else if (mission == 1)then
      ! update the gracesecs of first epoch
      if(line(1:26) == "    start_time_epoch_secs:")write(line(27:36),'(i10)')new_start  
      ! update the gracesecs of first epoch
      if(line(1:25) == "    stop_time_epoch_secs:") write(line(26:35),'(i10)')new_end    
      ! update number of observations
      if(line(1:16) == "    num_records:")write(line(17:26),'(i10)')n_obs

      ! end of header
      if(line(1:20) == "# End of YAML header")end_of_header = .true.
    endif

    ! write to output file the line of the header
    write(luacc(4),'(a)')line
  enddo
!-------------------------------------------------------------------------------------------------------------- 


!-------------------------------------------------------------------------------------------------------------- 
! now write all the data records to the output file
  do iobs=1,n_obs

! ACC1B data
    if(L1B_type == "ACC" .or. L1B_type == "ACT" .or. L1B_type == "ACH")then
      write(luacc(4),101)nint(acc_obs_extended(iobs,1)),sat,acc_obs_extended(iobs,2:ncols),nint(acc_obs_extended(iobs,ncols+1))
101   format(i9,1x,a1,9e23.15,1x,i8.8)

! KBR1B data
    else if (L1B_type == "KBR")then
      stop "KBR1B not coded yet - please add it yourself!"
    else if (L1B_type == "KBR")then
      stop "KBR1B not coded yet - please add it yourself!"
    endif
  enddo
!-------------------------------------------------------------------------------------------------------------- 
  end



!

