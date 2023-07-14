  program reformat_TN13

! program to take the TN13 file format (with data in GRCOF2 lines) and output as:
! start_yyyymmdd end_yyyymmdd     C_1,0           C_1,1             S_1,1
! 20210901.0000 20211001.0000 -8.098178278e-10 -2.330012540e-10 +3.081041730e-10
!
! P. Tregoning
! 16 February 2022

  implicit none
  
  character*100 :: TN13_file       ! e.g. TN-13_GEOC_JPL_RL06.txt
  character*101 :: line
  integer*4     :: ioerr
  character*100 :: epochs,C10,C11,S11
  

! get file name from command line
  call getarg(1,TN13_file)
  
! open the TN13 file
  open(10,file=trim(TN13_file),status='old',iostat=ioerr)
  if(ioerr /= 0)then
    call status_update('FATAL','UTIL','reformat_TN13',TN13_file,"Error opening input TN13 file",0)
  endif
  
! read through the file until we find a GRCOF2 line
  line = " "
  do while (line(1:5) /= "GRCOF")
    read(10,'(a)')line
  enddo
  backspace(10)
  
  
! now, loop through to the end of the file
  ioerr = 0
  do while (ioerr == 0)
    ! read the C_1,0 line
    read(10,'(a)',iostat=ioerr,end=1000)line
    if(ioerr ==0)then
      C10 = line(17:32)
      
      ! read the deg_1, order_1 line
      read(10,'(a)')line
      C11 = line(17:32)
      S11 = line(34:49)
      
      ! read the start/stop date string
      epochs = line(75:101)
      
      ! output them all in the required order on the same line
      write(*,'(a,1x,a,1x,a,1x,a)')trim(epochs),trim(C10),trim(C11),trim(S11)
    endif
    
  enddo
  
1000 continue
  call status_update('STATUS','UTIl','reformat_TN13',' ',"End of program",0)
  
  end
  
