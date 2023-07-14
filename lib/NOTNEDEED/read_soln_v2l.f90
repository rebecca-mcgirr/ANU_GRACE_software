 subroutine read_soln_v2 (gracenam, lu, apriori, soln, prm_input, VCV_obs, VCV_flag,imascon,imsc_tide)
!  subroutine read_soln_v2 (lu, N)

! Program read_soln reading the solution from graceful to gracekal and gracekal to gracekal
! LL 14 January 2013 modified the 29th of January
! APP130304: Modified to use write_soln include file and to pass data as input parameters
! PT130528:  add the GPS antenna offset parameter values to the read of the vcv file
! PT131107:  read and pass out the ascii descriptors of the estimated parameters

 implicit none

! include '../includes/gracefit.h90'
 include '../includes/grace.param'
 include '../includes/write_soln.h90'

!-- PARAMETERS --

 character*30 , intent(out) :: prm_input(maxparm)   ! the descriptors of the estimated parameters included in the input file
 integer*4    , intent(out) :: imascon              ! pointer to the parameter number of the first mascon parameter
 integer*4    , intent(out) :: imsc_tide            ! pointer to the first tidal mascon parameter
 character*8 gracenam
 character*200 message 
 character*10 message2
 character*60 message3
 integer*4 :: lu, i, N, j, VCV_flag, indx, len !LL nparam !length
 DOUBLE PRECISION :: apriori(maxparm), soln(maxparm), VCV_obs(maxparm, maxparm)
 integer*4 trimlen

! skip 3 lines
 do i=1, 3
   read(lu, '(a)') message
 enddo

! Now read the number of parameters
! PT150618: replace the (32x,i15) with (32x,i9)
  read(lu, '(32x,i9)') nparam
  write(message,'(a,i6)')'number of parameters' , nparam
  call report_stat('STATUS',gracenam,'lib/read_soln_v2',' ',message,0)

print*, 'TOTOTOTOTO'

  N = nparam

!! SOLUTION A PRIORI !!
! PT131107: change this to read in the parameter character descriptor as well as the numerical values
 do while ( message2(2:10) .ne. "PARAMETER" )
   read(lu, '(a)', end = 1000) message2
   if ( message2(2:10) .eq. "PARAMETER" ) then
     imascon = 0
     imsc_tide = 0
     do i = 1, N
! PT140610: identify and save the value of the first mascon parameter and first mascon tide parameter
       read(lu, '(a30, f17.7, f17.7, f17.7)') prm_input(i), apriori(i), soln(i), VCV_obs(i, i) 
       if(imascon   == 0 .and. prm_input(i)(8:9)  == "MC" )imascon   = i
       if(imsc_tide == 0 .and. prm_input(i)(8:10) == "TMC")imsc_tide = i
     enddo
   endif
 enddo

! APP130321: Only read in VCV matrix if VCV_flag is non-zero
 if ( VCV_flag .ne. 0 ) then
   do while ( message2(3:5) .ne. "VCV" )
     read(lu, '(a)', end = 1001) message2
     if ( message2(3:5) .eq. "VCV" ) then
       do i = 1, N
         read(lu, *) (VCV_obs(i, j), j = 1, N)
       enddo
     endif
   enddo
 endif

 return

! PT130322: what do these two lines do?
 1000 write(message,'(a)')'Error reading VCV file. Line "SOLUTION A PRIORI" not found'
      call report_stat('FATAL',gracenam,'read_soln_v2',' ',message,0)
 1001 write(message,'(a)')'Error reading VCV file. Line "SOLUTION VCV" not found'
      call report_stat('FATAL',gracenam,'read_soln_v2',' ',message,0)

 return
 end subroutine read_soln_v2
