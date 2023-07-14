  subroutine GTORB_record_length(luin,calling_prog,gtorbfile,tmp_file_recl)

! subroutine to determine the record length of a GTORB file based on the information contained within the file
! 
! a) Prior to 09/06/2017, there was no information on this in the files. From 27/08/2015 to 08/06/2017 the record
!    length was 318788.
!
! b) From 09/06/2017 to -->  the record length has been the first I*4 variable in the file

! P. Tregoning
! 9 June 2017
!
! PT180907: add fatal stop if error in opening the GTORB file

  use gtorb_mod

  implicit none

! argument variables
  integer*4    , intent(in) :: luin          ! unit number of GTORB file
  character*(*), intent(in) :: calling_prog  ! name of calling program
  character*80 , intent(in) :: gtorbfile     ! name of GTORB file
  integer*4     ,intent(out):: tmp_file_recl ! value decided upon for record length of GTORB file

! local variables 
  integer*4                 :: ioerr
  character*200             :: message


! PT170609: open the GTORB file generically to then read the correct record length. It is an I*4 number in the first record
  open(luin,file=gtorbfile,access='direct',status='old',form='unformatted',recl=1000,iostat=ioerr)
  if(ioerr /= 0)then
    call status_update('FATAL',calling_prog,'gtorb_record_length',gtorbfile,"Error opening file. Does it exist?",ioerr)
  endif

! read the first number in the file
  read(luin,rec=1,iostat=ioerr)tmp_file_recl
  close(luin)

! check that a valid number was read
! PT190306: increase this check number from 1 million to 4 million to allow over 40000 mascons
  if(ioerr == 0 .and. abs(tmp_file_recl) < 4000000)then
    write(message,'(a,i8)')"Record length of binary GTORB file: ",tmp_file_recl
    call status_update('STATUS','UTIL','read_gtorb_file',gtorbfile,message,0)
  else
! assume the old default record length
    tmp_file_recl = GTORB_recl_old
    write(message,'(a,i8)')"Old-style GTORB file. Record length set to: ",tmp_file_recl
! PT190306: fix info in error message call
    call status_update('STATUS','LIB','gtorb_record_length',gtorbfile,message,0)
  endif


  return
  end subroutine GTORB_record_length
