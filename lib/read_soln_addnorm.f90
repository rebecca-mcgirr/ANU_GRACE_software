 subroutine read_soln_addnorm (gracenam,soln_type,soln_number,lu,maxparam,nparam,apriori,soln,prm_input,VCV_obs,VCV_flag &
                          , sTid, sScl, sBias, sEmp, sMasc, n_emp, nScl, nBias)

! Subroutine to read a addnorm solution FIT/VCV file and return the a priori, solution and VCV.
!
! PT190301: subroutine modified from read_soln_v3 to accommodate the minor differences in the addnorm output files

 implicit none

  character*8 , intent(in)  :: gracenam                   ! name of calling program
  character*5 , intent(in)  :: soln_type                  ! type of addnorm solution file ("FIT" or "VCV")
  integer*4   , intent(in)  :: soln_number                ! the particular set of ICs that are to be returned from the addnorm soln file
  integer*4,    intent(in)  :: lu                         ! unit number of VCV file to be read
  integer*4,    intent(in ) :: maxparam                   ! dimensioning variable of the arrays
  logical,      intent(in ) :: VCV_flag                   ! flag to indicate whether to read the full VCV or just the solution
  integer*4,    intent(out) :: nparam                     ! number of parameters in this particular solution
  real(kind=8), intent(out) :: apriori(maxparam)          !   apriori values of all parameters in this solution
  real(kind=8), intent(out) :: soln(maxparam)             ! estimated values of all parameters in this solution
  character*30, intent(out) :: prm_input(maxparam)        ! the descriptors of the estimated parameters included in the input file
  real(kind=8), intent(out) :: VCV_obs(maxparam,maxparam) ! VCV of the solution
  integer*4,    intent(out) :: sScl(2)                    ! pointers to the start of the accelerometer scale parameters in the input solution
  integer*4,    intent(out) :: sBias(2)                   ! pointers to the start of the accelerometer bias parameters in the input solution
  integer*4,    intent(out) :: n_emp                      ! number of empirical "per-rev" parameters in the solution in the input solution
  integer*4,    intent(out) :: nScl                       ! number of scale parameters per satellite in the solution in the input solution
  integer*4,    intent(out) :: nBias                      ! number of bias parameters  per satellite in the solution in the input solution
  integer*4,    intent(out) :: sEmp(2)                    ! pointers to the start of the accelerometer bias parameters in the input solution
  integer*4,    intent(out) :: sMasc                      ! pointers to the parameter number of the first mascon parameter in the input solution
  integer*4,    intent(out) :: sTid                       ! pointer to the first tidal amplitude parameter in the input solution

! local variables
  integer*4 :: tmp_soln,nmsc,n_files

  character*200 message , message2
  integer*4 ::  i, j, indx
  integer*4 :: trimlen
  logical   :: debug

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            declare some variables               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  prm_input = " "

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            define pointers for ICs              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  sMasc = 1    ! we only need this one at this stage .... diff_ternarys requires only this !
  nparam = 0   ! we will return in this variable the number of mascons

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         start reading the input file            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the first line has been read, so continue reading ...
! line 2 contains the number of files
! PT191231: Julia has has changed this to "n_files" beginning after column 74. Not sure why, but make it the same here ...
!  read(lu, '(a68,i7)')message, n_files
  read(lu, '(a74,i7,22x,i6)')message, n_files,nparam

! skip over the file information
 do i=1, n_files
   read(lu, '(a)') message
 enddo

!! SOLUTION A PRIORI. Move to the appropriate line in the file !!
  debug=.false.
  do while ( message2(2:10) .ne. "PARAMETER" )
    read(lu, '(a)') message2
    if(debug)print*,"skipping line:",message2
  enddo

! PT190301: in GRACESIM/GRACEFIT we know in the file headers how many parameters there are. In ADDNORM output files we don't have
!           that information. Therefore, we need different logic here.
!
!           find the set of ICs that relate to the solution requested
  tmp_soln = 0
  i=0
  sMasc = 0    ! initialise the counter for the pointer to the first mascon parameter
  do while (tmp_soln /= soln_number)
    read(lu,'(i3)')tmp_soln
   enddo
! now read in all the ICs for each satellite for this requested solution
  backspace(lu)
  do while (tmp_soln == soln_number)

    sMasc = sMasc + 1
    if(soln_type(1:4) == " VCV")then
      read(lu, '(i3,a27, f17.7, f19.9, f15.7)')tmp_soln, prm_input(sMasc), apriori(sMasc), soln(sMasc), VCV_obs(sMasc, sMasc)
    else if (soln_type(1:3) == "FIT")then
      read(lu, '(i3,a27, f17.7,12x, f19.9, f12.7)')tmp_soln, prm_input(sMasc), apriori(sMasc), soln(sMasc), VCV_obs(sMasc, sMasc)
    else
      write(message,'(a,a)')"Unknown solution file type: ",soln_type
      call status_update('FATAL',gracenam,'read_soln_addnorm',' ',message,0)
    endif

    if(tmp_soln /= soln_number)then
      write(message,'(a,i4)')"Extracted ICs for orbit number: ",soln_number
      call status_update('STATUS',gracenam,'read_soln_addnorm',' ',message,0)
      backspace(lu)
    endif
  enddo

! now, skip over the remaining IC lines until we reach the mascon entries
  !message = " "
  do while (message(8:9) /= "MC")
    read(lu,'(a)')message
  enddo
  backspace(lu)

! now, simply read the mascons
  nmsc = sMasc - 1 

  do while (message(8:9) == "MC")
    read(lu,'(a)',end=1000)message
    if(message(8:9) == "MC")then
      nmsc = nmsc + 1
      if(soln_type(1:4) == " VCV")then
        read(message,'(a30, f17.7, f19.9, f15.7)') prm_input(nmsc), apriori(nmsc), soln(nmsc), VCV_obs(nmsc, nmsc)
      else
        read(message,'(a30, f17.7,12x, f19.9, f15.7)') prm_input(nmsc), apriori(nmsc), soln(nmsc), VCV_obs(nmsc, nmsc)
      endif
    endif

  enddo

1000  call status_update('STATUS',gracenam,'lib/read_soln_addnorm',' ',"Not reading the VCV matrix ...",0)

  return
  end subroutine read_soln_addnorm



