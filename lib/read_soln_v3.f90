 subroutine read_soln_v3 (gracenam, lu, VCV_flag &
                          , sTid, sScl, sBias, sEmp, sMasc, n_emp, nScl, nBias)
!  subroutine read_soln_v2 (lu, N)

! Subroutine to read a GRACEFIT/KAL/ solution VCV file and return the a priori, solution and VCV.
!
! NOTES: This version (version 3) doesn't use any include files, thus making it more portable
  use soln_mod

 implicit none

  character*8 , intent(in)  :: gracenam                   ! name of calling program
  integer*4,    intent(in)  :: lu                         ! unit number of VCV file to be read
!  integer*4,    intent(in ) :: maxparam                   ! dimensioning variable of the arrays
  logical,      intent(in ) :: VCV_flag                   ! flag to indicate whether to read the full VCV or just the solution

!  integer*4,    intent(out) :: nparam                     ! number of parameters in this particular solution
!  real(kind=8), intent(out) :: apriori(maxparam)          !   apriori values of all parameters in this solution
!  real(kind=8), intent(out) :: soln(maxparam)             ! estimated values of all parameters in this solution
!  character*30, intent(out) :: prm_input(maxparam)        ! the descriptors of the estimated parameters included in the input file
!  real(kind=8), intent(out) :: VCV_obs(maxparam,maxparam) ! VCV of the solution
  integer*4,    intent(out) :: sScl(2)                    ! pointers to the start of the accelerometer scale parameters in the input solution
  integer*4,    intent(out) :: sBias(2)                   ! pointers to the start of the accelerometer bias parameters in the input solution
  integer*4,    intent(out) :: n_emp                      ! number of empirical "per-rev" parameters in the solution in the input solution
  integer*4,    intent(out) :: nScl                       ! number of scale parameters per satellite in the solution in the input solution
  integer*4,    intent(out) :: nBias                      ! number of bias parameters  per satellite in the solution in the input solution
  integer*4,    intent(out) :: sEmp(2)                    ! pointers to the start of the accelerometer bias parameters in the input solution
  integer*4,    intent(out) :: sMasc                      ! pointers to the parameter number of the first mascon parameter in the input solution
  integer*4,    intent(out) :: sTid                       ! pointer to the first tidal amplitude parameter in the input solution

! temporary parameter numbers
  integer*4 :: tmp_nmasc, tmp_ntid

  character*200 message 
  character*10 message2
  integer*4 ::  i, j, indx
  integer*4 :: trimlen
  logical   :: debug

!print*,'in read_soln_v3',gracenam,lu,vcv_flag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            declare some variables               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  debug = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         start reading the input file            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! skip 3 lines (the first line has been read already, so need skip only lines 2-4)
 do i=1, 3
   read(lu, '(a)') message
if(debug)print*,message
 enddo

! Now read the total number of parameters
  read(lu, '(32x,3i7)') nparam,tmp_nmasc,tmp_ntid
  write(message,'(a,i6,a,i6,a,i6)')'Total parameters: ',nparam,' Number of mascons: ',tmp_nmasc &
                                   ,' Number of tidal amplitudes: ',tmp_ntid*2   ! PT160429: there are two tidal amplitudes per tide per tidal mascon 
  call status_update('STATUS',gracenam,'lib/read_soln_v3',' ',message,0)

!! PT190312: now that we have the number of parameters, allocate the arrays
   maxparm = nparam
   allocate(apriori(maxparm))
   allocate(soln(maxparm))
   allocate(VCV_obs(maxparm,maxparm))
!   allocate(prmnam(maxparm))
   prmnam = " "

!! SOLUTION A PRIORI !!
! PT131107: change this to read in the parameter character descriptor as well as the numerical values
!debug=.false.
 do while ( message2(2:10) .ne. "PARAMETER" )
   read(lu, '(a)', end = 1000) message2
   !print *,message2
   if ( message2(2:10) .eq. "PARAMETER" ) then
if(debug)print*,'found PARAMETER line. nparam=',nparam,maxparm
     sMasc = 0
     sTid = 0
     sScl = 0
     sBias = 0
     n_emp = 0

if(debug)then
!  do i=1,nparam
!    print*,'debug:',i,prmnam(i),'xx'
!  enddo
!  stop 'stopped in read_soln_v3 because debug is on'
endif
     do i = 1, nparam
       ! PT170713: changed the format here to have more significant figures on the solution vector values
!if(debug)print*,i,' about to read it',maxparm
       read(lu, '(a30, f17.7, f19.9, f15.7)') prmnam(i), apriori(i), soln(i), VCV_obs(i, i) 
       !print *,i, prmnam(i), apriori(i), soln(i), VCV_obs(i, i)
if(debug)print*,"from .vcv file:",i,apriori(i), soln(i), VCV_obs(i, i)
! PT140610: identify and save the value of the first occurence of several parameter types
!      mascons
       if(sMasc   == 0 .and. prmnam(i)(8:9)  == "MC" )then
        sMasc   = i
       endif
!      tidal amplitudes
       if(sTid == 0 .and. prmnam(i)(8:10) == "TMC")sTid = i
!      accelerometer scales
       if(sScl(1) == 0 .and. prmnam(i)(6:10) == "SAT A" .and. prmnam(i)(13:15) == "scl")then
         sScl(1) = i
         nScl = 3
       endif
       if(sScl(2) == 0 .and. prmnam(i)(6:10) == "SAT B" .and. prmnam(i)(13:15) == "scl")then
         sScl(2) = i
         nScl = 3
       endif
!      accelerometer biases
       if(sBias(1) == 0 .and. prmnam(i)(6:10) == "SAT A" .and. prmnam(i)(13:14) == "bs")then
         sBias(1) = i
         nBias = 3
       endif
       if(sBias(2) == 0 .and. prmnam(i)(6:10) == "SAT B" .and. prmnam(i)(13:14) == "bs")then
         sBias(2) = i
         nBias = 3
       endif
!      empirical parameters
       if(prmnam(i)(14:15) == "pr")then
         n_emp = n_emp + 1
         if(sEmp(1) == 0 .and. prmnam(i)(6:10) == "SAT A" )sEmp(1) = i
         if(sEmp(2) == 0 .and. prmnam(i)(6:10) == "SAT B" )sEmp(2) = i
       endif
!if(debug)print*,i, ' end of do loop'
     enddo
! convert total number of empirical force parameters into the number per satellite
     n_emp = n_emp/2
   endif
 enddo


! APP130321: Only read in VCV matrix if VCV_flag is non-zero
 if ( VCV_flag ) then
  call status_update('STATUS',gracenam,'lib/read_soln_v3',' ',"Reading VCV matrix",0)
   do while ( message2(3:5) .ne. "VCV" )
     read(lu, '(a)', end = 1001) message2
     if ( message2(3:5) .eq. "VCV" ) then
       do i = 1, nparam
! PT141124: speed this up by reading only the lower triangular part
           read(lu, *) (VCV_obs(i, j), j = 1, i)
! PT170308: fix bug here by putting this in a loop from 1 to i
           do j=1,i
             VCV_obs(j,i) = VCV_obs(i,j)
           enddo
       enddo
     endif
   enddo
 endif

! DEBUG
! print*,'read_soln_v3 maxparam =',maxparam
!  do i=25,45
!    print*,i,vcv_obs(i,25:45)
!  enddo
! call status_update('STATUS',gracenam,'lib/read_soln_v3',' ',"Finished reading in solution information for this file",0)
!stop 'stopped after reading vcv'

 return

! error messages if we've not properly read the input VCV file
 1000 write(message,'(a)')'Error reading VCV file. Line "SOLUTION A PRIORI" not found'
      call status_update('FATAL',gracenam,'read_soln_v3',' ',message,0)
 1001 write(message,'(a)')'Error reading VCV file. Line "SOLUTION VCV" not found'
      call status_update('WARNING',gracenam,'read_soln_v3',' ',message,0)

 return
 end subroutine read_soln_v3
