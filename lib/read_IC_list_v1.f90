  subroutine read_IC_list_v1(luin,tin,extra_params)

! subroutine to read all ICs for both satellites from a file that contains all info for a single day on a single line
!
! P. Tregoning
! 9 May 2017
!
! HM190502: modified subroutine to use soln_mod and allocate apriori arrays (consistent with read_soln_v3)
! PT190917: include the number of  mascons in the dimensioning of "apriori" and "soln", if required
  use soln_mod		 ! defines the prmnam array, soln and apriori vectors 
  use sat_mod            ! provides name of satellite being integrated
 
  implicit none

! passed variables
  integer*4,intent(in)      :: luin            ! unit number of input file (opened in INPUT_read_apriori)
  integer*4,intent(in)      :: tin             ! GRACE seconds of requested ICs
  integer*4,intent(in)      :: extra_params    ! either 0 or total_prim, depending on whether we want to use a priori mascon EWH values

! local variables
  logical       :: found
  real(kind=8)  :: ICs(19,2)               ! pos/vel/scl/bs/1pr/2pr/rpy  = 3+3+3+3+2+2+3 per satellite
  real(kind=8)  :: seconds
  integer*4     :: ioerr,i
  integer*4     :: date(5)
  integer*4     :: year, month, day        ! requested epoch for which we want ICs
  integer*4     :: tmpyr,tmpmonth,tmpday
  character*100 :: line  
  character(30) :: text(38)		   ! prmnam template

  data text/&
"  1. SAT A: X0   (m)          ",&
"  2. SAT A: Y0   (m)          ",&
"  3. SAT A: Z0   (m)          ",&
"  4. SAT A: XV0  (m/s)        ",&
"  5. SAT A: YV0  (m/s)        ",&
"  6. SAT A: ZV0  (m/s)        ",&
"  7. SAT A: sclx (n/a)        ",&
"  8. SAT A: scly (n/a)        ",&
"  9. SAT A: sclz (n/a)        ",&
" 10. SAT A: bsx  (um/s^2)     ",&
" 11. SAT A: bsy  (um/s^2)     ",&
" 12. SAT A: bsz  (um/s^2)     ",&
" 13. SAT A: 1prS (um/s^2)     ",&
" 14. SAT A: 1prC (um/s^2)     ",&
" 15. SAT A: 2prS (um/s^2)     ",&
" 16. SAT A: 2prC (um/s^2)     ",&
" 17. SAT A: rpyx (mrad)       ",&
" 18. SAT A: rpyy (mrad)       ",&
" 19. SAT A: rpyz (mrad)       ",&
"  1. SAT B: X0   (m)          ",&
"  2. SAT B: Y0   (m)          ",&
"  3. SAT B: Z0   (m)          ",&
"  4. SAT B: XV0  (m/s)        ",&
"  5. SAT B: YV0  (m/s)        ",&
"  6. SAT B: ZV0  (m/s)        ",&
"  7. SAT B: sclx (n/a)        ",&
"  8. SAT B: scly (n/a)        ",&
"  9. SAT B: sclz (n/a)        ",&
" 10. SAT B: bsx  (um/s^2)     ",&
" 11. SAT B: bsy  (um/s^2)     ",&
" 12. SAT B: bsz  (um/s^2)     ",&
" 13. SAT B: 1prS (um/s^2)     ",&
" 14. SAT B: 1prC (um/s^2)     ",&
" 15. SAT B: 2prS (um/s^2)     ",&
" 16. SAT B: 2prC (um/s^2)     ",&
" 17. SAT B: rpyx (mrad)       ",&
" 18. SAT B: rpyy (mrad)       ",&
" 19. SAT B: rpyz (mrad)       "/

! PT/HMcQ190917: change prmnam satellite names from A/B to C/D if required
  if(sat == "C" .or. sat == "D")then
    do i = 1,19
      text(i)(10:10) = "C"
      text(i+19)(10:10) = "D"
    enddo
  endif
     

! convert the start time (in GRACE seconds) to yyyy mm dd
  call gsec_to_ymdhms(tin,date,seconds)
  year=date(1)
  month=date(2)
  day=date(3)

! setup IC arrays
  maxparm=38
  allocate(apriori(maxparm+extra_params))
  allocate(soln(maxparm+extra_params))
  allocate(VCV_obs(maxparm+extra_params,maxparm+extra_params))
! now already allocated in soln_mod
!  allocate(prmnam(maxparm))
  prmnam = " "

! read the IC file until we find a line that matches the requested date
  ioerr = 0
  found = .false.

  do while (ioerr == 0 .and. .not. found)
    read(luin,*,iostat=ioerr,end=1000)tmpyr,tmpmonth,tmpday,(ICs(:,i),i=1,2)
    if(tmpyr == year .and. tmpmonth == month .and. tmpday == day)then
      found = .true.
      print*,'Found ICs for ',year,month,day," They are :",ICs
      close(luin)
      apriori(1:19) = ICs(1:19,1)
      apriori(20:38) = ICs(1:19,2)
      soln = apriori
    
! define the character names for the parameters
 
      do i=1,19
        prmnam(i)=text(i)
        prmnam(i+19)=text(i+19)
      enddo
      return
    endif
  end do

1000 continue
  if(ioerr /= 0)then
    call status_update('STATUS',"LIB",'read_IC_list_v1',' ' &
            ,"Error reading the IC file. Perhaps the requested date does not exist in the file.",ioerr)
  endif
  return
  
  end subroutine read_IC_list_v1
