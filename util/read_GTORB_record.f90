  program read_GTORB_record

! program to read a particular record in a binary GTORB file. Input variables are 
! a) name of GTORB file
! b) the record number to be read
! c) type of read ("C", "I", "R" for character, integer or real*8)
! d) for "I" or "R", the number of values to read
!

  use gtorb_mod   ! PT170609: GTORB record length information variables defined in here

  implicit none

  integer*4 irec,nvals
  character arg*100,gtorbfile*100,read_type*1,line*200
  integer*4, allocatable :: ivals(:)
  real*8,    allocatable :: rvals(:)
  real*4,    allocatable :: svals(:)
! variables that match actual GTORB records
  integer*4  timeout
  real*8     ICs(6)
  real*8     quat(4)

! first and last column for extracting partials
  integer*4 :: first_val, last_val
! number of mascons in the GTORB file
  integer*4 nmascons  

! local variables
  integer*4     :: ioerr
  character*200 :: message

! PT170609: make this bigger than likely necessary at this stage
  nmascons = 10000

! input GTORB file to be read
  call getarg(1,gtorbfile)

! help if no arguments on command line
  if (gtorbfile(1:1) == " ")then
    print*,"Runstring: read_GTORB_record <GTORB file> <record to read> <flag>"
    print*,"flag = C(haracter), I(nteger), R(eal*8), S(ingle precision R*4), G(= GRACE record)" &
             ," P(= particular column after the GRACE seconds, ICs and quats) [add"
    stop
  endif

! record to be interrogated
  call getarg(2,arg)
  read(arg,*)irec

! type of read of record
  call getarg(3,read_type)
  if(read_type /= "C" .and. read_type /= "T")then
!   number of values to read
    call getarg(4,arg)
    if(arg(1:1) == " ")then
      call status_update('FATAL','UTIL','read_GTORB_record',' ' &
                  ,"Must indicate how many records to read as 4th command line argument",0)
    endif
    read(arg,*)nvals
! PT140811: if we are extracting partials, allow to indicate the end column as well
    if(read_type == "P") then
      call getarg(5,arg)
      if(arg(1:1) /= " ")then
        first_val = nvals
        read(arg,*)last_val
        nvals = last_val
        if(last_val < first_val)then
          print*,'last column value (',last_val,') must be greater than or equal to first column value (',first_val,')'
          stop
        endif
      else
        print*,'Must enter the end column as well as start column when using option "P"'
        stop
      endif
    endif

    allocate(ivals(nvals))
    allocate(rvals(nvals))
    allocate(svals(nvals))
  endif

  if(read_type == "T")then
    allocate(rvals(nmascons))
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT170609: determine the correct record length. It is returned in file_recl
  call GTORB_record_length(10,'reda_GTORB_record',gtorbfile,file_recl)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now open the file again, using the correct record length for this particular file
  open(10,file=gtorbfile,access='direct',form='unformatted',recl=file_recl)

  if(read_type == "C")then
! read the record as a character line
    read(10,rec=irec)line(1:200)
    write(*,'(a,i6,a,a)')"Record ",irec,' ',line
  elseif(read_type == "I")then
! read the record as a series of integers
    read(10,rec=irec)ivals(1:nvals)
    print*,'Record ',irec,ivals
  elseif(read_type == "R")then
! read the record as a series of double precision floating point
    read(10,rec=irec)rvals(1:nvals)
    print*,'Record ',irec,rvals
  elseif(read_type == "S")then
! read the record as a series of double precision floating point
    read(10,rec=irec)svals(1:nvals)
    print*,'Record ',irec,svals
  elseif(read_type == "G")then    ! read a line containing GRACE seconds, ICs, quaternions, partials
    read(10,rec=irec)timeout,ICs,quat,rvals(1:nvals)
    print*,'Record ',irec,timeout,ICs,quat,rvals(1:nvals)
  elseif(read_type == "P")then    ! read and output the value of a particular column(s) (i.e. to check a particular partial(s))
    read(10,rec=irec)timeout,ICs,quat,rvals(1:nvals)
    print*,'Record ',irec,' column ',first_val,rvals(first_val:last_val)
! PT160622: add code to read the lines containing information on the tidal amplitudes
  else if (read_type == "T")then  ! read and output 4556 values, being the a priori tidal amplitudes for one constituent for either cosine or sine (depending on what line was chosen)
    read(10,rec=irec)rvals(1:nmascons)
    print*,'Record ',irec,' tidal amplitudes ',rvals(1:nmascons)
  endif

  end






