  program square_thrust_ACC1B

! This program will remove (as best we can) the presence of filtered orientation thrust events in the ACC1B 
! observations and replace them with a square thrust instead. The process involves knowing scale factors of
! the admittance of angular thrusts into linear accelerations, and we use scale factors derived at ANU. 
!
! The program works as follows:
!
! INPUT INFORMATION
! 1. A Level-1B thruster file provides information on timing and duration of thrust events
! 2. The total nominal force is scaled by scale factors (determined externally of this program) 
!    to replicate the right amount of force that leaks into linear accelerations
! 3. A synthetic, 200 Hz time series is generated, with the beginning of the thrust event at the time
!    specified in the THR1B input file. The thrust is represented as a square pulse.
! 4. A 3rd-order Butterworth filter is applied, output time series is then 10 Hz
! 5. A CRN filter (140.7 sec) filter is then applied, output time series is 1 Hz and yields the synthetic
!    ACC1B data for that particular thrust event.
! 6. This is then subtracted from the ACC1B observed accelerations to effectively remove the filtered thrust signals
! 7. A time series of square pulse thrusts is generated, based on the time of the thrust and the scaling factors to
!    provide the admittance. 
! 8. The square pulse thrusts are then added to the ACC1B observations, thus re-inserting the leaked thrusts
!
! This process effectively replaces the 140 second long, filtered thrust information with a 5sec square pulse. It will
! still not be right but it will be a lot closer to what the satellite actually underwent than the filtered version will be.

! Given that we have a 5sec time step in our orbit integrator, we choose the height of our square pulse so that the
! area under the curve matches the (scaled) amplitude times the duration.  
!
!
! Paul Tregoning (Paul.Tregoning@anu.edu.au)
! 4 April 2012
!
! MODS:
!
!  PT120405: the first version of this program didn't work. We must generate first a time series of the Butterworth filtered 10 Hz data, 
!                     then run the CRN filter across it. This is necessary because many of the thrusts overlap and will appear within the window 
!                     length of the CRN filter. Therefore, it is not correct to sum the CRN-filtered thrusts afterwards - the CRN filter needs to
!                     be run on 1047 samples already containing all information of thrusts during that period. 
!
!                     I need to restructure the program to 
!                        1) generate for the whole day a time series of Butterworth filtered thrusts, 10 Hz sampling
!                        2) run the CRN filter over the whole day, with a 1047 sample sliding window stepping forward 1 sec at a time
!                              - within each 1047 samples, see whether there is a non-zero element. If so, run the filter. If not, output is zero.
!
! PT130814: the original program, synthetic_ACC1B, has been beefed up to also do the following steps:
!           1. read in the actual ACC1B observations
!           2. subtract the modelled filtered thrusts from the actual obs
!           3. model and add back in the square pulse thrusts.

  implicit none

  integer*4 maxepoch             ! maximum number of epochs
  parameter (maxepoch = 86400)  ! number of epochs in a day at 1 Hz sampling

  character THRfile*100          !   input Level1B thrust file
  character in_ACC1Bfile*100     !   input Level1B accelerometer file 
  character out_ACC1Bfile*100    !   output accelerometer file (in ACC1B format) with filtered thrusts replaced by square pulses
  character arg*100,line*100,char1*1  ! bits and pieces
  character message*150          ! string for status_update
  character*8 acc_code(maxepoch) ! bit flags at the end of the ACC1B observation records
  character debug*5              ! debug flag from the command line (6th argument)
  character*6 thrust_type(6)     ! character names for each thrust type
  data thrust_type /    "  +yaw", "  -yaw", "+pitch", "-pitch", " +roll", " -roll" /

  integer*4 flg                  ! flag to set the output option (has value 1 to 4)
  integer*4 hr,mn,isec           ! as the name suggests
  integer*4 start_secs           ! time of first thrust in the THR1B file
  integer*4 nthrust              ! number of thrusts in the THR1B file
  integer*4 axis                 ! the accelerometer axis being modelled
  integer*4 count_thrusts(6)     ! number of each type of thrust
  integer*4 av_thrust(6)         ! average duration of each thrust type
  integer*4 ioerr,line_num,ithrust,i,j ! bits and pieces

  double precision :: acc_data(maxepoch,10)      ! input data from the ACC1B file
  double precision accscale(6,2,3) ! six thrust pairs, 2 satellites, 3 axes
  double precision synth_10Hz(maxepoch*10,6,3)   ! 10Hz synthetic square thrust observations
  double precision :: synth_acc(maxepoch,6,3)    ! for each epoch, 6 thrust pairs, 3 axes
  double precision :: square_acc(maxepoch,6,3)   ! 1Hz scaled synthetic square thrust observations           
  double precision epoch

! satellite identifier. GRACE A is 1, GRACE B is 2
  integer isat

! variables for reading the THR1B file
  integer tmpint(12),tmpint2(2), tmpdur(12), tmpsec, tmpmicrosec

! CRN filter code variables
  logical do_crn_filt,skipFn
  integer counter,iepoch
  double precision tmp_crn(1407),filtered_val

! initialise some variables
  synth_acc = 0.d0
  accscale = 0.d0
  nthrust = 0
  epoch = 0.d0

! decode command line arguments
  call getarg(1,THRfile)
  if(THRfile(1:1).eq.' ')then
    print*,"Runstring: square_thrust_ACC1B THR1B_2010-09-15_A_02.asc A/B flg ACC1B_2010-09-15_A_02.asc ACC1B_2010-09-15_A_Sy.asc"
    print*,'where flg setting will output  1: filtered, modelled thrusts'
    print*,'                               2: square pulse modelled thrusts'
    print*,'                               3: ACC1B obs minus filtered modelled thrusts'
    print*,'                               4: ACC1B obs - filtered + square pulse thrusts'
    print*,'                               5: original ACC1B obs only'
    stop 
  endif
  call getarg(2,char1)
  if (char1.eq."A")isat = 1
  if (char1.eq."B")isat = 2
  call getarg(3,arg)
  read(arg,*)flg
  call getarg(4,in_ACC1Bfile)
  call getarg(5,out_ACC1Bfile)
  call getarg(6,debug)

  call status_update('STATUS','UTIL','square_thrust_ACC1B',THRfile,"Generating synthetic filtered thrusts from file ",0)
  call status_update('STATUS','UTIL','square_thrust_ACC1B',out_ACC1Bfile,"Output will be written to file ",0)

! open the output ACC1B file
  open(20,file=out_ACC1Bfile,status='unknown')

! Open the THR1B file
  open(10,file=THRfile,status='old',iostat=ioerr)
  if(ioerr.ne.0)call status_update('FATAL','UTIL','square_thrust_ACC1B',THRfile,"Error opening input file.",0)

! read the header of the THR1B file. Do this quickly rather than rigorously for now, but extract out the time of first obs
  line = " "
  do while (line(1:13).ne."END OF HEADER")
    read(10,'(a)',iostat=ioerr)line
    if(ioerr.ne.0)then
      call status_update('FATAL','UTIL','square_thrust_ACC1B',THRfile,"Did not find the END OF HEADER line. Error reading file: ",0)
      stop
    endif
! extract out the time of first obs. Use this to get the number of seconds at UT0000 for this day. This becomes our time counter for storing all the obs.
! TIME FIRST OBS(SEC PAST EPOCH): 373291353.654915 (2011-10-31 00:02:33.65)       
    if(line(1:14).eq."TIME FIRST OBS")then
       read(line(32:41),'(i10)')start_secs
       read(line(62:69),'(i2,1x,i2,1x,i2)')hr,mn,isec
       start_secs = start_secs - hr*3600 - mn*60 - isec
       write(message,'(a11,a,i11)')line(51:61),' 00UT (in grace seconds):',start_secs
       call status_update('STATUS','UTIL','square_thrust_ACC1B',' ',message,0)
    endif
  enddo

! open the ACC1B file
!  call input_openFile(11,in_ACC1Bfile,'old','UTIL','square_thrust_ACC1B','FATAL','Error opening ACC1B file')
  open(11,file=in_ACC1Bfile,status='old')

! read in all the ACC1B observations. skip down to after the header information.
  line = " "
  do while (line(1:13) /= "END OF HEADER")
    read(11,'(a)')line
  enddo
! now read in 86400 lines of information
  do iepoch = 1,86400
    read(11,*)acc_data(iepoch,1),char1,(acc_data(iepoch,j),j=2,10),acc_code(iepoch)
  enddo
  call status_update('STATUS','UTIL','square_thrust_ACC1B',in_ACC1Bfile,"Read acceleration obs from file ",0)

! define the scale factors to convert the nominal thrust force into accelerations. There
! are 6 pairs of thrusters (+yaw, +pitch, -yaw, -pitch, +roll, -roll) and each has a different
! scaling factor. The ones used here were determined by S. McClusky and P. Tregoning from
! an assessment of a few different months of actual GRACE B ACC1B observations. We will need
! to do the same for GRACE A, so set it up here accordingly.
!
! the array "scale" is as follows:
! column 1: rotation type (+Y, +P, -Y, -P, +R, -R)
! column 2: satellite (1: GRACE A; 2: GRACE B)
! column 3: axis in science reference frame (1: X, 2: Y, 3: Z)
!!!!!!!!!!!!!!!!
!!!! GRACE A !!!  Use average (or thereabouts) of March2008, Dec 20101 and Oct2011
!!!!!!!!!!!!!!!!
  accscale(1,1,1) = -0.0555d3   ! leakage of +yaw into x-axis of GRACE A
  accscale(3,1,1) = -0.0696d3   ! leakage of -yaw into x-axis of GRACE A
  accscale(1,1,2) =  0.0d3    ! no discernable leakage of +yaw into y-axis of GRACE A
  accscale(3,1,2) = -0.0d3    ! no discernable leakage of -yaw into y-axis of GRACE A
  accscale(1,1,3) =  0.0327d3   ! leakage of +yaw into z-axis of GRACE A
  accscale(3,1,3) =  0.0244d3   ! leakage of -yaw into z-axis of GRACE A

  accscale(2,1,1) = 0.0d0     ! no discernable leakage of +pitch into x-axis of GRACE A
  accscale(4,1,1) = 0.0d0     ! no discernable leakage of -pitch into x-axis of GRACE A
  accscale(2,1,2) = 0.035d3     ! leakage of +pitch into y-axis of GRACE A
  accscale(4,1,2) = 0.035d3     ! leakage of -pitch into y-axis of GRACE A
  accscale(2,1,3) = 0.0d0     ! no discernable leakage of +pitch into z-axis of GRACE A
  accscale(4,1,3) = 0.009d3     ! leakage of -pitch into z-axis of GRACE A

  accscale(5,1,1) =  -0.1373d3  ! leakage of +roll into x-axis of GRACE A
  accscale(6,1,1) =  -0.1227d3  ! leakage of -roll into x-axis of GRACE A
  accscale(5,1,2) =  0.0d3    ! no discernable leakage of +roll into y-axis of GRACE A
  accscale(6,1,2) =  0.0d3    ! no discernable leakage of -roll into y-axis of GRACE A
  accscale(5,1,3) =  0.0285d3   ! leakage of +roll into z-axis of GRACE A
  accscale(6,1,3) =  0.0361d3   ! leakage of -roll into z-axis of GRACE A

!!!!!!!!!!!!!!!!
!!!! GRACE B !!!  Use average (or thereabouts) of of March2008, Dec2010 and Oct2011
!!!!!!!!!!!!!!!!
  accscale(1,2,1) = 0.0d3   ! no discernable leakage of +yaw into x-axis of GRACE B
  accscale(3,2,1) = 0.0d3   ! no discernable leakage of -yaw into x-axis of GRACE B
  accscale(1,2,2) =  0.0152d3   ! leakage of +yaw into y-axis of GRACE B
  accscale(3,2,2) = -0.0108d3   ! leakage of -yaw into y-axis of GRACE B
  accscale(1,2,3) =  0.0220d3   ! leakage of +yaw into z-axis of GRACE B
  accscale(3,2,3) =  0.0246d3   ! leakage of -yaw into z-axis of GRACE B

  accscale(2,2,1) = 0.0d0   ! no discernable leakage of +pitch into x-axis of GRACE B
  accscale(4,2,1) = 0.0d0   ! no discernable leakage of -pitch into x-axis of GRACE B
  accscale(2,2,2) = 0.0d0   ! no discernable leakage of +pitch into y-axis of GRACE B
  accscale(4,2,2) = 0.0d0   ! no discernable leakage of -pitch into y-axis of GRACE B
  accscale(2,2,3) =  0.015d3   ! leakage of +pitch into z-axis of GRACE B
  accscale(4,2,3) = -0.023d3   ! leakage of -pitch into z-axis of GRACE B

  accscale(5,2,1) =  0.1009d3   ! leakage of +roll into x-axis of GRACE B
  accscale(6,2,1) = -0.2132d3   ! leakage of -roll into x-axis of GRACE B
  accscale(5,2,2) =  0.0120d3   ! leakage of +roll into y-axis of GRACE B
  accscale(6,2,2) =  0.0320d3   ! leakage of -roll into y-axis of GRACE B
  accscale(5,2,3) =  0.0192d3   ! leakage of +roll into z-axis of GRACE B
  accscale(6,2,3) =  0.0390d3   ! leakage of -roll into z-axis of GRACE B

! ok, so we are at the first line of the thrust entries in the file. Read the file one line
! at a time, compute the leakage of the thrust into linear accelerations on each accelerometer
! component and then store the values against the relevant epoch.
  ioerr = 0
  line_num = 1
  do while (ioerr.eq.0)
    read(10,*,iostat=ioerr,end=1000)tmpsec,tmpmicrosec,char1,char1,(tmpint(i),i=1,12),(tmpint2(i),i=1,2),(tmpdur(i),i=1,12)
    nthrust = nthrust + 1
    if(debug == "debug")write(*,'(a,i5,a,i10,a1,i6.6,6i6)')'Thrust',nthrust,' at time : ',tmpsec,".",tmpmicrosec,(tmpdur(i),i=1,6)
!100 format(i9,i7,1x,a1,1x,a1,2x,12i7,
    do ithrust = 1,6
      if(tmpdur(ithrust).gt.0)then  ! a thrust occurred on this pair
        count_thrusts(ithrust) = count_thrusts(ithrust) + 1
        av_thrust(ithrust) = av_thrust(ithrust) + tmpdur(ithrust)

        do axis = 1,3
          if(accscale(ithrust,isat,axis).ne.0.d0)then
! Generate 10Hz data for this event
             epoch = dble(tmpsec) + dble(tmpmicrosec)/1.d6

! PT130815: for unknown reasons, the observed thrusts are time-shifted by -1 second for thrusts of duration > 1 sec. Account
!           for that here, by subtracting one second when generating the square pulse for filtering. We will insert the square
!           pulse thrust at the right epoch, though.
             if(tmpdur(ithrust).gt.1000.d0)then
                write(message,'(a,i10,a,a6,a,i4)')'Epoch ',tmpsec,' Thrust ',thrust_type(ithrust),' duration ',tmpdur(ithrust)
                call status_update('STATUS','UTIL','square_thrust_ACC1B',' ',message,0)
                call generate_synth_10Hz(start_secs, ithrust, axis, maxepoch, epoch-1.d0 &
                   , tmpdur(ithrust) , accscale(ithrust,isat,axis), synth_10Hz )
            else
               call generate_synth_10Hz(start_secs, ithrust, axis, maxepoch, epoch &
                   , tmpdur(ithrust) , accscale(ithrust,isat,axis), synth_10Hz )
            endif

! Generate a scaled, square pulse that will be reinserted later on.
             iepoch = nint(epoch)-start_secs
             square_acc( iepoch,ithrust,axis ) = 20.d0*tmpdur(ithrust)*1.d-3 /accscale(ithrust,isat,axis)
          endif     ! end of check on whether scale factor is non-zero
        enddo     ! end of loop over the three axes
      endif      ! end of check on length of thrust duration
    enddo       ! end of loop over the 6 thruster pairs

  enddo          

1000 call status_update('STATUS','UTIL','square_thrust_ACC1B',' ',"10 Hz observations for all thrusts have been computed",0)


! We now have in synth_10Hz array a time series of butterworth-filtered thrusts for each thrust pair for each axis from 0000UT to 2400UT. Now need to filter each pair/axis with the CRN filter, stepping along 1 second at a time with a sliding window of 1047 samples width.
!
! The filter requires 73 seconds either side of the first point, so we start here at second 74 and finish at 86326

! set this to false if you don't have the file Fn.dat. It takes time to run, so unset it once you've generated the file Fn.dat
  skipFn = .true. 
  do ithrust = 1,6
    write(message,'(a,a6,a,i4,a,f6.1,a)')'CRN filter for thrusts: ',thrust_type(ithrust),' (total ',count_thrusts(ithrust) &
                     ,' thrusts). Av duration:',av_thrust(ithrust)/dble(count_thrusts(ithrust)),' milliseconds.'
    call status_update('STATUS','UTIL','square_thrust_ACC1B',' ',message,0)

    do axis = 1,3
      do iepoch = 70,86329    ! loop from the first second to the end
! generate the 1047 length 10Hz samples in a vector
        do_crn_filt = .false.
        counter = iepoch*10 + 4 - 703 - 1
        do i=1,1407
          counter = counter + 1
!  print*,'ithrust, axis, counter,iepoch',ithrust,axis,counter,iepoch,synth_10Hz(counter,ithrust,axis)
          tmp_CRN(i) = synth_10Hz(counter,ithrust,axis)
!          print*,i,counter,ithrust,axis,synth_10Hz(counter,ithrust,axis)
! check whether there are any non-zero values in this 1047 sample window
          if(dabs(tmp_CRN(i)).gt.1.0d-9)do_crn_filt = .true.
        enddo
! invoke the filter if there were non-zero values
        if(do_crn_filt)then

          call CRN_filter(skipFn,1407,tmp_CRN,filtered_val)
!          print*,"CRN filter for thrust",ithrust," axis",axis," seconds of day",iepoch,filtered_val
! PT120412: the epoch is 2 seconds early for some reason. Set it back 2 secs
          synth_acc(iepoch+2,ithrust,axis) = filtered_val
        else
!          print*,'no non-zero values centred on ',iepoch,' thrust',ithrust,' axis',axis
!          synth_acc(iepoch+2,ithrust,axis) = 0.d0
          synth_acc(iepoch+3,ithrust,axis) = 0.d0
        endif
      enddo
    enddo
  enddo
          

  call write_ACC1B(20, flg, isat, out_ACC1Bfile, start_secs, synth_acc, square_acc, acc_data, acc_code, maxepoch)

  print*,"Have a nice space flight ... !"
  end
!*****************************************************************************************

!**************  Subroutine write_ACC1B ***********************

  subroutine write_ACC1B(iout, flg, isat, outfile, start_secs, synth_acc, square_acc, acc_data, acc_code, maxepoch)


! subroutine to write out the synthetic acceleration values caused by modelled thrust. Output format is to match the ACC1B files.
!
! IN:
!       iout                                   :  unit number of output file
!       flg                                    :  flag to set what data will be output
!       isat                                   :  1 (GRACE A) or 2 (GRACE B)
!       start_secs                             :  seconds since Jan2000 of the 0000UT epoch of the day in question. This will be the start of the file
!       synth_acc (maxepoch,6,3)  :  synthetic observations. 1 row per epoch. Columns are the thruster pairs (+Y, +P, -Y, -P, +R, -R), 3rd dimension are axes (X,Y,Z)
!       square_acc (maxepoch,6,3)
!       acc_data   (maxepoch,10)
!
! OUT:  writes formatted, ascii output to the output file
!
! Paul Tregoning (Paul.Tregoning@anu.edu.au)
! 5 April 2012

  implicit none

  integer*4 flg                                                    ! output option flag (value 1 to 5)
  integer iout,iepoch,ithrust,axis,maxepoch,start_secs,isat,j
  character csat*1,  outfile*80, message*100
  double precision, intent(in) :: synth_acc(maxepoch,6,3)          ! filtered thrust model
  double precision, intent(in) :: square_acc(maxepoch,6,3)         ! square thrust model
  double precision, intent(in) :: acc_data(maxepoch,10)            ! original ACC1B data
  character*8, intent(in)      :: acc_code(maxepoch)

  double precision             :: accel_obs(3)                     ! total acceleration to output
  double precision             :: filt_obs(3)                      ! total acceleration to output
  double precision             :: square_obs(3)                    ! total acceleration to output
  double precision             :: zeroes(6)

  zeroes = 0.d0

! why did I ever store GRACE A/B as an integer .... ?
  if (isat.eq.1)csat = "A"
  if (isat.eq.2)csat = "B"

  if (flg == 1 ) then
    write(message,'(a)')"Writing out filtered thrust model"
  else if (flg == 2 ) then
    write(message,'(a)')"Writing out square thrust model"
  else if (flg == 3 ) then
    write(message,'(a)')"Writing out ACC1B minus filtered thrust model"
  else if (flg == 4) then
    write(message,'(a)')"Writing out ACC1B minus filtered thrust + square thrust model"
  else if (flg == 5) then
    write(message,'(a)')"Writing out original ACC1B observations"
  endif
  call status_update('STATUS','UTIL','square_thrust_ACC1B',outfile,message,0)

! write the file header
  call write_ACC1B_header(iout,csat,outfile )

! now, for 86400 epochs (being 1/second for an entire day), sum up the total accelerations per axis, then write them out along with the values per pair per axis
  do iepoch = 1,86400
     accel_obs = 0.d0
     filt_obs = 0.d0
     square_obs = 0.d0
     do axis = 1,3
        do ithrust = 1,6
          filt_obs(axis)   =   filt_obs(axis) +  synth_acc(iepoch,ithrust,axis)
          square_obs(axis) = square_obs(axis) + square_acc(iepoch,ithrust,axis)
        enddo
        if(flg == 1) then         ! output just the filtered thrust model values
          accel_obs(axis) = filt_obs(axis)
        else if (flg == 2) then   ! output just the square thrust model values
          accel_obs(axis) = square_obs(axis)
        else if (flg == 3) then   ! ACC1B minus filtered
          accel_obs(axis) = acc_data(iepoch,axis+1)*1.d6 - filt_obs(axis)
        else if (flg == 4) then   ! ACC1B minus filtered
          accel_obs(axis) = acc_data(iepoch,axis+1)*1.d6 - filt_obs(axis) + square_obs(axis)
        else if (flg == 5) then   ! ACC1B minus filtered
          accel_obs(axis) = acc_data(iepoch,axis+1)*1.d6
        else
          call status_update('FATAL','UTIL','square_thrust_ACC1B',' ',"Unknown flg value. Must be integer 1 to 4",0)
        endif

     enddo

! now just write it out
! PT120405: the computations done generated accelerations in micrometres/sec^2. Write out in m/s^2
!    write(*,1000)start_secs+iepoch-1, csat, accel_obs,zeroes,"00000000"
    write(iout,1000)start_secs+iepoch-1, csat, accel_obs*1.d-6,(acc_data(iepoch,j),j=5,10),acc_code(iepoch)

1000 format(i9,1x,a1,1x,9e23.14,1x,a8)
  enddo

  return
  end

  

!!!!!!!!&&&&&&&&&&&&&&&&&!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine generate_synth_10Hz(start_secs, ithrust, axis, maxepoch, epoch, duration, tscale, synth_10Hz )

! subroutine to set up the synthetic, square thrust pulse at 200 Hz, run
! a 3rd-order Butterworth filter over it and downsample to 10 Hz. The synthetic observations that
! are non-zero are then stored in the array synth_10Hz at the appropriate epoch.
!
! Paul Tregoning (Paul.Tregoning@anu.edu.au)
! 4 April 2012.
!
! IN:
! start_secs : seconds since Jan2000 for 0000UT of the day being processed
!   ithrust  :    thruster pair number (from 1-6, they are +Y,+P,-Y,-P,+R,-R)
!   axis     :    science frame axis number (1: X, 2: Y, 3: Z)
!   maxepoch :    dimensions output array
!   epoch    :    epoch in seconds (floating point value) of start of thrust
!   duration :    duration of thrust event (integer value in milliseconds)
!   tscale   :    scale factor for this thruster pair/axis
!   synth_10Hz :   output array of synthetic accelerometer observations at 10Hz

  implicit none

! passed variables
  integer maxepoch, axis, ithrust, duration,start_secs
  double precision  tscale, synth_10Hz(maxepoch*10,6,3),epoch

! local variables
  double precision t_start,t_thrust,fs
  integer nsecs,thrust_counter,i,counter
  parameter (nsecs = 5)
  double precision ts_200Hz(2*nsecs*200)
  double precision ts_butt(2*nsecs*200),ts_10Hz((2*nsecs*200)/10)

! set the 200 Hz array to zero
  ts_200Hz = 0.d0
  ts_10Hz = 0.d0
  ts_butt = 0.d0

! first, set up a 200 Hz time series with a square thrust pulse starting at time "epoch" with a width
! of "duration". For the Butterworth filter, 5 seconds before and after the thrust is more than enough. the time series must start at least 150 sec before the start of the
  fs = 1.d0/200.d0         ! 200 Hz sampling
  t_thrust = int(epoch)

! 200 Hz sampling gives us a resolution of 5 milliseconds, but the thrust times and durations go down to 1 ms and 1 ns, so we'll have to round it off here.
  t_start = t_thrust - nsecs
  thrust_counter = int(float(nsecs)/fs) + int(int((epoch-t_thrust)*1.d3)/5.d0)
!  print*,"t_start, epoch, thrust_counter and duration (ms):",t_start,epoch,thrust_counter,duration,int(t_thrust-start_secs)

! insert the square pulse, starting at time "epoch" which is in location "thrust_counter"
!  print*,"inserting square pulse in 200 Hz elements",thrust_counter,thrust_counter+int(duration/fs)
  do i = thrust_counter,thrust_counter+int(duration*1.d-3/fs)    ! duration in seconds x 200 obs/sec 
!    print*,"inserting square pulse in array row",i
    ts_200Hz(i) = 20.d0  / tscale ! * 1.06d0    ! PT120405: the scale factors were done using "20" rather than using 10mN/mass of sat. 1.06 corrects for this.
  enddo

! ok, now filter it with the Butterworth filter.
  call butt_filter(2*nsecs*200,ts_200Hz,ts_butt,ts_10Hz)
! debug
!  do i=1,2*nsecs*200
!    print*,i,ts_butt(i)
!  enddo
!  stop
! The Butterworth filter has introduced a 140 ms delay (the Butterworth delay), which equates to 14 epochs
! in the 10 Hz time series. So, the correct entry in the array for our start epoch is now element 15 of ts_10 Hz.

! Now we have the square pulse signal once-filtered and downsampled to 10 Hz in array "ts_10Hz".
! The thrust is the 1130th entry in the 200Hz array, so it is the 58th entry in the 10Hz array

! add the values into the synth_10Hz array.
  
  do i = 1,nsecs*20 - 14
    counter = (t_start - start_secs)*10  + (i -1)    ! this is the number of seconds since the start of the day, and is the row in the synth_10Hz array
   synth_10Hz(counter,ithrust,axis) = synth_10Hz(counter,ithrust,axis) + ts_10Hz(i+14)
!    print*,i,counter,start_secs,t_thrust,t_start, ithrust, axis,synth_10Hz(counter,ithrust,axis)
  enddo

!  stop


! that's all folks.
  return

  end



! !----------------------------------------------------------------------
!   SUB-ROUTINE "butt_filter"
!
!   This sub-routine applies a 3rd order low pass Butterworth filter to
!   the input array to emulate the onboard filtering of the GRACE accelerometer data.
! the Butterworth coefficients come from Christopher Watson from his matlab script
!  butterA(1) =  1.0000   
!  butterA(2) = -2.8116    
!  butterA(3) =  2.6405   
!  butterA(4) = -0.8281
!  butterB(1) = 0.0954d-3    
!  butterB(2) = 0.2863d-3    
!  butterB(3) = 0.2863d-3    
!  butterB(4) = 0.0954d-3

!
!   Input: datain, "nepoch" element array
!   Ouptut dataout, "nepoch" element array
!
! PT080930: modified to accept any number of elements, rather than hardwired
!           to 250
! -----------------------------------------------------------------------      
      SUBROUTINE butt_filter(nepoch,datain,dataout,data10hZ)
       
      IMPLICIT NONE
      INTEGER i,j,nepoch,nepoch10Hz
      REAL*8 datain(nepoch), dataout(nepoch), data10Hz(nepoch), a(21), b(21)  
      
      DATA a/1.00000000000000,-2.8116,2.6405,-0.8281, & 
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0/
     
      DATA b/0.0000954,0.0002863,0.0002863, 0.0000954,&
              0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0/

     
 
      dataout(1)=b(1)*datain(1)
      
      DO 10 i = 2,21
      
        dataout(i)=0.0
        DO 11 j = 1,i
          dataout(i)=dataout(i)+b(j)*datain(i-j+1)
 11     CONTINUE
        DO 12 j = 1,i-1
          dataout(i)=dataout(i)-a(j+1)*dataout(i-j)
 12     CONTINUE
 
 10   CONTINUE
      

      DO 13 i = 22,nepoch
      
        dataout(i)=0.0
        DO 14 j = 1,21
          dataout(i)=dataout(i)+b(j)*datain(i-j+1)
 14     CONTINUE
        DO 15 j = 1,20
          dataout(i)=dataout(i)-a(j+1)*dataout(i-j)
 15     CONTINUE
 
 13   CONTINUE


! now need to resample the data at 10 Hz (rather than the 200 Hz that we chose to start with). So we take every 20th point and throw the rest away ....
      nepoch10Hz = 0
      do i = 1,nepoch,20
        nepoch10Hz = nepoch10Hz + 1
        data10Hz(nepoch10Hz) = dataout(i)
!        print*,i,data10Hz(nepoch10Hz)
      enddo
      

      return
      END
      
! -----------------------------------------------------------------------
!   END OF SUB-ROUTINE "butt_filter"
! -----------------------------------------------------------------------



      SUBROUTINE CRN_filter(skipFn,nepoch,datain,filtered_val)

! CRN filter to emulate the step between GRACE Level-1A and GRACE Level-1B accelerometer data
!
! Filter coefficients adapted from matlab code provided by Ulrich Meyer (AIUB)
!
! P. Tregoning
! 1 April 2012

  implicit none

  logical skipFn
  integer nepoch, junk
  double precision datain(nepoch),filtered_val

! variables for the filter
  integer Nf, Nh, NB,i,j,k,m, n
  double precision f0, fs, Tf, B, Nc, Hk(1407), sumHkn, Fnorm, sumHki, pi, Fn(1407)

  pi = 4.d0*datan(1.d0)

  fs=10.d0    ! samples per s
  Tf=140.7d0   ! filterlength in s
  Nf=fs*Tf   ! samples in filterlength: 1407
  Nh=(Nf-1)/2    ! summation limit: 703

  B=0.035d0    ! target low-pass bandwith
  NB=nint(B*Tf)   ! frequency bins in passband: 5

  Nc=7.d0    ! self convolution number

  if (.not.skipFn) then
    Hk = 0
    do k = -Nh,Nh
      do m = -NB,NB
        if(m.ne.k)then
          Hk(k+Nh+1)=Hk(k+Nh+1)+(dsin(pi*(k-m)/Nc)/dsin(pi*(k-m)/Nf))**Nc
        else
          Hk(k+Nh+1)=Hk(k+Nh+1)+(Nf/Nc)**Nc
        endif
      enddo
    enddo
    print*,'computed Hk'

    f0=0.00037;   ! basic frequency from J2 (ca. 2/5400s)


    do n = -Nh, Nh
    print*,'n-loop',n
      sumHkn=0.d0
      do k = -Nh, Nh
        sumHkn=sumHkn+Hk(k+Nh+1)*dcos(2.0*pi*k*n/Nf) 
      enddo
      Fnorm=0.d0  
      do i = -Nh, Nh
        sumHki=0.d0
        do k = -Nh, Nh
          sumHki=sumHki+Hk(k+Nh+1)*dcos(2.0*pi*k*i/Nf) 
        enddo
        Fnorm=Fnorm+dcos(2.d0*pi*f0*i/fs)*sumHki
      enddo
      Fn(n+Nh+1)=sumHkn/Fnorm
    enddo
    print*,'computed Fn'

    open(30,file='Fn.dat',status='unknown')
    do i=-Nh,Nh
      write(30,*)i+Nh+1,Fn(i+Nh+1)
    enddo  
    close(30)
  else  ! Fn exists already
    open(30,file='Fn.dat',status='unknown')
    do i=-Nh,Nh
      read(30,*)junk,Fn(i+Nh+1)
!      print*,'In CRN_filter  i,Fn(i+Nh+1): ',i,Fn(i+Nh+1)
    enddo  
    close(30)
  endif

! ok, so now I have to write the application of the CRN filter, given Fn and the series to be filtered ..... we only loop through n=-Nh,Nh because we only compute this for the centre of the 1047 values (the main program loops through all the seconds and sends the appropriate 1047 values down to this subroutine)

! debug
!  do i=1,1407
!     print*,"in CRN_filter tmpvals: ",i,datain(i)
!  enddo

    i = 704
    filtered_val = 0.d0
    do n = -Nh, Nh
         filtered_val = filtered_val + Fn(n+Nh+1)*datain(i-n)
!         print*,'in CRN_filter  n,Fn(n+Nh+1), datain(i-n):',n,Fn(n+Nh+1), i-n,datain(i-n),filtered_val
    enddo

  return
  end

!!!!!**************** write_header ******************
!
  subroutine write_ACC1B_header(iout, csat,outfile )

! subroutine to write the ACC1B header
! P. Tregoning
! 11 April 2012

  implicit none

  integer iout
  character csat*1,c_date*8,c_time*10,outfile*80

! a sample header (from 2011-10-31 GRACE A file)
!PRODUCER AGENCY               : NASA                                            
!PRODUCER INSTITUTION          : JPL                                             
!FILE TYPE ipACC1BF            : 8                                               
!FILE FORMAT 0=BINARY 1=ASCII  : 1                                               
!NUMBER OF HEADER RECORDS      : 23                                              
!SOFTWARE VERSION              : $Id: ACC_compress.c 1.75 10/24/06 22:27:10 gl $ 
!SOFTWARE LINK TIME            : @(#) 2007-01-03 19:38:30 glk  j2                
!REFERENCE DOCUMENTATION       : GRACE Level 1 Software Handbook                 
!SATELLITE NAME                : GRACE A                                         
!SENSOR NAME                   : ACC  GRACE_ICU_cal.txt 1.5 03/10/02             
!TIME EPOCH (GPS TIME)         : 2000-01-01 12:00:00                             
!TIME FIRST OBS(SEC PAST EPOCH): 373291200.000000 (2011-10-31 00:00:00.00)       
!TIME LAST OBS(SEC PAST EPOCH) : 373377599.000000 (2011-10-31 23:59:59.00)       
!NUMBER OF DATA RECORDS        : 86400                                           
!PRODUCT CREATE START TIME(UTC): 2011-11-13 00:32:41 by l0tol1                   
!PRODUCT CREATE END TIME(UTC)  : 2011-11-13 00:33:26 by l0tol1                   
!FILESIZE (BYTES)              : 19316156                                        
!FILENAME                      : ACC1B_2011-10-31_A_01.asc                       
!PROCESS LEVEL (1A OR 1B)      : 1B                                              
!INPUT FILE NAME               : ACC1A<-ACC1A_2011-10-31_A_01.dat                
!INPUT FILE TIME TAG (UTC)     : ACC1A<-2011-11-01 01:41:47 by l0tol1            
!INPUT FILE NAME               : CLK1B<-CLK1B_2011-10-31_A_01.dat                
!INPUT FILE TIME TAG (UTC)     : CLK1B<-2011-11-13 00:27:36 by l0tol1            
!END OF HEADER                                                                   

  write(iout,'(a)')"PRODUCER AGENCY               : ANU"
  write(iout,'(a)')"PRODUCER INSTITUTION          : RSES"
  write(iout,'(a)')"FILE TYPE ipACC1BF            : 8"
  write(iout,'(a)')"FILE FORMAT 0=BINARY 1=ASCII  : 1"
  write(iout,'(a)')"NUMBER OF HEADER RECORDS      : 19"
  write(iout,'(a)')"SOFTWARE VERSION              : square_thrust_ACC1B"
  write(iout,'(a)')"SOFTWARE LINK TIME            : 2013-08-14" 
  write(iout,'(a)')"REFERENCE DOCUMENTATION       : McClusky et al (JGR paper)"
  write(iout,'(a,a1)')"SATELLITE NAME                : GRACE ",csat
  write(iout,'(a)')"SENSOR NAME                   : ACC  GRACE_ICU_cal.txt 1.5 03/10/02"
  write(iout,'(a)')"TIME EPOCH (GPS TIME)         : 2000-01-01 12:00:00"
  write(iout,'(a)')"TIME FIRST OBS(SEC PAST EPOCH):"
  write(iout,'(a)')"TIME LAST OBS(SEC PAST EPOCH) :"
  write(iout,'(a)')"NUMBER OF DATA RECORDS        : 86400"

! get the system date/time
     call date_and_time(date=c_date,  &  ! character(len=8) ccyymmdd
                      time=c_time )                ! character(len=10) hhmmss.sss
!                      zone=c_zone, &               ! character(len=10) +/-hhmm (time zone)
!                      values=ivalues)                ! integer ivalues(8) all of the above 

  write(iout,'(a,a4,2(a1,a2),4x,2(a2,a1),a6)')"PRODUCT CREATE START TIME(UTC): " &
                                                               ,c_date(1:4),"-",c_date(5:6),"-",c_date(7:8)  &
                                                               ,c_time(1:2),":",c_time(3:4),":",c_time(5:10)
  write(iout,'(a,a4,2(a1,a2),4x,2(a2,a1),a6)')"PRODUCT CREATE END TIME(UTC)  : " &
                                                               ,c_date(1:4),"-",c_date(5:6),"-",c_date(7:8)  &
                                                               ,c_time(1:2),":",c_time(3:4),":",c_time(5:10)
  write(iout,'(a,a)')"FILENAME                      : ",outfile
  write(iout,'(a)')"PROCESS LEVEL (1A OR 1B)      : 1B"
  write(iout,'(a)')"END OF HEADER"

  return
  end

