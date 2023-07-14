  program addnorm2vcv

! program to extract a set of ICs and mascons from an addnorm.vcv file and create a single vcv file that can
! be run back through GRACEORB. This will permit multi-day mascon estimates to be fed back a priori through the
! orbit integrator.
!
! P. Tregoning
! 30 May 2017
!
! PT180314: provide three options for output solution type: a priori, adjustment or updated value
! JP191121: correct format to read input vcv file generated using addnorm

  implicit none

! command line arguments
  character*100  :: vcv_in                 ! input vcv file (in addnorm.vcv format)
  character*100  :: vcv_out                ! output vcv file name
  character*5    :: pos,vel,bias,scl,msc   ! flags for apriori or estimated value for each parameter type
  integer*4      :: year, month, day       ! date of ICs required

! local variables
  character*256  :: message,arg
  integer*4      :: luin,luout,lutmp       ! input/output unit numbers
  integer*4      :: ioerr,i
  character*100  :: line
  logical        :: found

! addnorm file variables
  integer*4      :: n_orbits               ! number of orbits stacked in the input addnorm vcv  file
  integer*4      :: tmporbit,tmpyear,tmpmonth,tmpday
  integer*4      :: iorbit                 ! pointer to the particular orbit that matches the requested IC date
  integer*4      :: iparam                 ! IC parameter number
  character*7    :: csat                   ! IC satellite string
  character*24   :: cparam                 ! IC parameter character string
  real(kind=8)   :: apriori,vector,sigma   ! IC parameter values
  integer*4      :: n_ICs,n_msc            ! number of ICs and mascons to be output  
  character*100  :: junk(6)                ! junk variable to read in things we don't want and then ignore               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! decode command line
! input file
  call getarg(1,vcv_in)
  if(vcv_in(1:1) == "")then
    write(message,'(a,a)')"Runstring: addnorm2vcv addnorm.vcv output.vcv 2012 07 03 pos_A vel_A scl_A bs__A msc_E " &
      ,"(A:apriori, C:correction, E:estimated)"
    call status_update('FATAL','UTIL','addnorm2vcv',' ',message,0)
  endif
! output file
  call getarg(2,vcv_out)
! date
  call getarg(3,arg)
  read(arg,*)year
  call getarg(4,arg)
  read(arg,*)month
  call getarg(5,arg)
  read(arg,*)day
! parameter flags
  call getarg(6,pos)
  call getarg(7,vel)
  call getarg(8,scl)
  call getarg(9,bias)
  call getarg(10,msc)
print*,'flags to use:',pos,' ',vel,' ',scl,' ',bias,' ',msc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! open the files
  luin = 10
  luout = 11
  !lutmp = 12
  open(luin,file=vcv_in,status='old',iostat=ioerr)
  if(ioerr /= 0)then
    call status_update('FATAL','UTIL','addnorm2vcv',vcv_in,"Error opening input addnorm file",0)
  endif  
  open(luout,file=vcv_out,status='unknown',iostat=ioerr)
  open(lutmp,file="vcv.tmp",status='unknown')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the header of the addnorm vcv file
  read(luin,'(a)')line
! PT180913: check whether it is an old-style VCV or one that stipulates which program created it on the first line
  if(line(1:2) == "V2")then
    ! PT180913: check that it is a VCV file
    ! RM190402: changed line(14:16) to line(15:17)
    if(line(15:17) /= "VCV")then
      call status_update('FATAL','UTIL','addnorm2vcv',vcv_in,"Input file is not a VCV file format!",0)
    else
      write(message,'(a,a)')"Input file created by ",line(4:12)
      call status_update('STATUS','UTIL','addnorm2vcv',vcv_in,message,0)
    endif
  else if (line(1:8) == "Solution")then
      call status_update('STATUS','UTIL','addnorm2vcv',vcv_in,"Old-style VCV file with incomplete header information",0)
      rewind(luin)
  endif



! second line tells us how many orbits are stacked in the addnorm.vcv file
! JP: replace '(69x,i5)' by '(74x,i7)' to read n orbits
  read(luin,'(74x,i7)')n_orbits
  write(message,'(a,i5,a)')"Found ",n_orbits," orbits in the input file"
  call status_update('STATUS','UTIL','addnorm2vcv',vcv_in,message,0)

! PT180913: read the 2 lines containing regularisation information
! RM190402: should only read 1 line, was skipping over first orbit
  read(luin,'(a)')line
  !read(luin,'(a)')line

! now read the IC times and check for a match with requested date. 
  found = .false.
  do i=1,n_orbits
    read(luin,'(a5,i4,a9,3i5)')junk(1),tmporbit,junk(2),tmpyear,tmpmonth,tmpday
    if(.not. found .and. (year == tmpyear .and. month == tmpmonth .and. day == tmpday) )then
      found = .true.
      iorbit = tmporbit
      write(message,'(a,i3,a,3i5,a)')"Orbit ",iorbit," is a match for requested IC date (",year,month,day,")"
      call status_update('STATUS','UTIL','addnorm2vcv',' ',message,0)
    endif
  enddo   
! read down to the first IC line.
  do while (line(2:10) /= "PARAMETER")
    read(luin,'(a)')line
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now, read through the ICs and extract out the ICs for the particular orbit that matches the requested IC date
  tmporbit = 0
  do while (tmporbit /= iorbit)
    read(luin,*)tmporbit
  enddo
  backspace(luin)

! now we are at the start of the GRACE A ICs for the requested orbit epoch. Transfer them to the output file
  call status_update('STATUS','UTIL','addnorm2vcv',vcv_out,"Writing out ICs",0)
  n_ICs = 0
  do while (tmporbit == iorbit)
    read(luin,'(a)')line
    if(line(8:10) == "SAT")then    ! it is still an IC line
      read(line(1:3),*)tmporbit
      if(tmporbit == iorbit)then
        n_ICs = n_ICs + 1
        ! strip off the orbit number and write out either the a priori or the estimated value
! PT180226: changed line(4:81) to line(5:81)
        read(line(5:81),'(i3,a7,1x,a16,3f17.7)')iparam,csat,cparam,apriori,vector,sigma
        ! Position IC
        if(cparam(2:3) == "0 ")then
          ! we need to change the units from mm to m. It is actually metres in the addnorm vcv file, but is wrongly labelled
          cparam(6:9) = "(m) "
          if(pos == "pos_A")write(lutmp,'(i3,a2,a7,a16,2x,f17.7,f19.9,f15.7)')iparam,". ",csat,cparam,apriori,apriori,sigma
         if(pos == "pos_C")write(lutmp,'(i3,a2,a7,a16,2x,f17.7,f19.9,f15.7)')iparam,". ",csat,cparam,apriori,vector-apriori,sigma
          if(pos == "pos_E")write(lutmp,'(i3,a2,a7,a16,2x,f17.7,f19.9,f15.7)')iparam,". ",csat,cparam,apriori,vector,sigma
        endif
        ! Velocity IC
        if(cparam(2:3) == "V0")then
          if(vel == "vel_A")write(lutmp,'(i3,a2,a7,a16,2x,f17.7,f19.9,f15.7)')iparam,". ",csat,cparam,apriori,apriori,sigma
         if(vel == "vel_C")write(lutmp,'(i3,a2,a7,a16,2x,f17.7,f19.9,f15.7)')iparam,". ",csat,cparam,apriori,vector-apriori,sigma
          if(vel == "vel_E")write(lutmp,'(i3,a2,a7,a16,2x,f17.7,f19.9,f15.7)')iparam,". ",csat,cparam,apriori,vector,sigma
        endif
        ! Scale IC
        if(cparam(1:3) == "scl")then
          if(scl == "scl_A")write(lutmp,'(i3,a2,a7,a16,2x,f17.7,f19.9,f15.7)')iparam,". ",csat,cparam,apriori,apriori,sigma
         if(scl == "scl_C")write(lutmp,'(i3,a2,a7,a16,2x,f17.7,f19.9,f15.7)')iparam,". ",csat,cparam,apriori,vector-apriori,sigma
          if(scl == "scl_E")write(lutmp,'(i3,a2,a7,a16,2x,f17.7,f19.9,f15.7)')iparam,". ",csat,cparam,apriori,vector,sigma
        endif
        ! Bias IC
        if(cparam(1:2) == "bs")then
          if(bias == "bs__N")print*,'no bias in vcv'
          if(bias == "bs__A")write(lutmp,'(i3,a2,a7,a16,2x,f17.7,f19.9,f15.7)')iparam,". ",csat,cparam,apriori,apriori,sigma
          if(bias == "bs__C")write(lutmp,'(i3,a2,a7,a16,2x,f17.7,f19.9,f15.7)')iparam,". ",csat,cparam,apriori,vector-apriori,sigma
          if(bias == "bs__E")write(lutmp,'(i3,a2,a7,a16,2x,f17.7,f19.9,f15.7)')iparam,". ",csat,cparam,apriori,vector,sigma
        endif
      endif
    else
      tmporbit =100000
    endif
  enddo

! now skip over the rest of the ICs 
  do while (line(7:9) /= " MC")
    read(luin,'(a)')line
  enddo
  backspace(luin)

! now, transfer the mascon values to the output file
  call status_update('STATUS','UTIL','addnorm2vcv',vcv_out,"Writing out mascons",0)
  n_msc = 0
  do while (line(3:14) /= "VCV SOLUTION")
    read(luin,'(a)')line
    if(line(3:14) /= "VCV SOLUTION")then
      read(line(1:81),'(i5,1x,a12,12x,3f17.7)')iparam,cparam,apriori,vector,sigma
      ! Mascon
      if(cparam(2:3) == "MC")then
        n_msc = n_msc + 1
        if(msc == "msc_A")write(lutmp,'(i5,a1,a12,12x,f17.7,f19.9,f15.7)')iparam,".",cparam,apriori,apriori,sigma
       if(msc == "msc_C")write(lutmp,'(i5,a1,a12,12x,f17.7,f19.9,f15.7)')iparam,".",cparam,apriori,vector-apriori,sigma
       if(msc == "msc_N")write(lutmp,'(i5,a1,a12,12x,f17.7,f19.9,f15.7)')iparam,".",cparam,apriori,apriori-vector,sigma
        if(msc == "msc_E")write(lutmp,'(i5,a1,a12,12x,f17.7,f19.9,f15.7)')iparam,".",cparam,apriori,vector,sigma
      endif
    endif
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now, write the header for the file, then copy the temporary info over to the output file
! PT220914: we should write "GRACEFIT" here instead of "ADDNORM" because it is not in addnorm format. This also pushes
!           the "VCV" one column to the right, which is needed for it to be read properly
  write(luout,'(a,3i5)')"V2 GRACEFIT   VCV solution extracted solution: ",year,month,day
  write(luout,'(a)')"Reference GTORB files: (not relevant)"
  write(luout,'(a)')"Reference GTORB files span: (not relevant)"
  write(luout,'(a)')"Reference GRACEFIT solution epoch: (not relevant)"
  write(luout,'(a,3i7)')"Number of estimated parameters: ",n_ICs+n_msc,n_msc,0
  write(luout,'(a)')"Number of KB observations and used: (not relevant)"
  write(luout,'(a)')"Number of missing KB observations:  (not relevant)"
  write(luout,'(a)')"Number of KBRR misfits:   (not relevant)"
  write(luout,'(a)')"Number of KBRR observations used:  (not relevant)"
  write(luout,'(a)')"Number of GPS observations used for GRACE A:  (not relevant)"
  write(luout,'(a)')"Number of relevant missing GPS observations for GRACE A:  (not relevant)"
  write(luout,'(a)')"Number of GPS observations used for GRACE B:  (not relevant)"
  write(luout,'(a)')"Number of relevant missing GPS observations for GRACE B:  (not relevant)"
  write(luout,'(a)')"Total Number of epochs:  (not relevant)"
  write(luout,'(a)')" "
  write(luout,'(a)')" SOLUTION A PRIORI AND VECTOR: "
  write(luout,'(a)')" PARAMETER                     A PRIORI             VECTOR            SIGMA"

! copy over the solution
  ioerr = 0
  rewind(lutmp)
  do while(ioerr == 0)
    read(lutmp,'(a)',iostat=ioerr, end=1000)line
    write(luout,'(a)')line
  enddo

1000 continue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  close(lutmp,status='delete')
  end












