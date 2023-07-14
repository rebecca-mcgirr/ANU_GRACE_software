  program scale_vcv

! program to read in a VCV file and output a file where the EWH values are scaled by a user-defined amount
!
! P. Tregoning
! 20 March 2020

  implicit none

  character*100 :: line,arg
  character*100 :: infile,outfile
  real(kind=8)  :: EWH
  real(kind=8)  :: scale_factor
  integer*4     :: i,ioerr
  integer, parameter :: luin=10,luout=20

! decode command line arguments
  call getarg(1,infile)
  if(infile(1:1) == "")then
    !call status_update('FATAL','UTIL','scale_vcv',' ',"Runstring: scale_vcv input.vcv output.vcv scale_factor",0)
    print*,"Runstring: scale_vcv input.vcv output.vcv scale_factor"
    stop
  endif
  call getarg(2,outfile)
  call getarg(3,arg)
  read(arg,*)scale_factor

  print*,'file ',infile,' will be scaled to ',scale_factor*100,' percent of the original mascon values'

! open the files
  open(luin,file=infile)
  open(luout,file=outfile)


! read through the file. Simply copy non-mascon lines as they are. Scale the mascon lines by the input scale factor
  ioerr = 0
  do while (ioerr == 0)
    read(luin,'(a)',iostat=ioerr,end=1000)line
    if(ioerr == 0)then
      if(line(7:9) == " MC")then
        ! it is a mascon line. Need to convert the mascon value to floating point, scale it then write it back
        read(line(48:64),*)EWH
        write(line(48:64),'(f17.7)')EWH*scale_factor
      endif

      write(luout,'(a)')line
    endif
  enddo

1000 close(luin)
     close(luout)

  print*,'End of program'
  end

