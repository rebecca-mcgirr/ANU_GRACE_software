  program diff_ternarys

! program to calculate and output the difference in mascon estimates on the ternary mascons between two files. The input files 
! need to be VCV file format (so that I can read them with read_soln_v3) but the primary mascon geometry doesn't need to be the
! same in the two files. 
!
! PT190301: modified to enable input files to be any of a GRACEFIT/GRACESIM .fit/.vcv or an ADDNORM .fit/.vcv
!
! P. Tregoning
! 5 December 2018
!
! MODS:
! PT200302: added the uncertainty of the primary mascon of each solution to each output ternary line
! PT200915: add the capability to use the a priori values of each field rather than the estimates
! PT220621: add a scale factor option, so that rates can be determined from the difference

  use mascon_mod
  use soln_mod

  implicit none

  character*200 :: message,message2
  character*150 :: soln_file1,soln_file2,mascon_file1,mascon_file2,output_file
  real(kind=8)  :: delta_time
  
! variables for determining whether to use a priori or estimated mascon values
  character*2 :: which_soln(2)

! unit numbers
  integer*4 :: lumsc1,lumsc2,lusoln1,lusoln2,luout

! variables for reading solution VCV info
  integer*4 :: sTid, sScl, sBias, sEmp, sMasc, n_emp, nScl, nBias
  integer*4,parameter :: total_params=52000
  character*30, allocatable :: prm_input(:)      ! no idea how big this needs to be
  real(kind=8),allocatable :: VCV_obs_local(:,:)
  real(kind=8),allocatable :: soln_sigma(:)

! variables for extracting ternary info
  integer*4 :: tern_num
  integer*4,parameter :: maxterns=1485118
  real(kind=8),allocatable :: soln1_ternarys(:,:),soln2_ternarys(:,:)

! decimate ternaries by this value
  integer*4 :: n_decimate
  
! counters
  integer*4 :: iprim,itern,i

! other variables
  real(kind=8)  :: pi
  character*100 :: line,arg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! allocate variables
  allocate(prm_input(total_params))
  allocate(VCV_obs_local(total_params,total_params))
  allocate(soln1_ternarys(1485118,5))  ! PT200302: added a 5th variable, being the sigma of the primary mascon estimate
  allocate(soln2_ternarys(1485118,5))
!  allocate(apriori(total_params))
!  allocate(soln(total_params))
!  allocate(vcv_obs(total_params,total_params))

! define values
  pi = 4.d0*datan(1.d0)

! decode runstring
  call getarg(1,mascon_file1)
  if(mascon_file1(1:1) == " ")then
    print*,"Runstring: diff_ternarys mascon_file1 soln_file1 mascon_file2 soln_file2 output_file [A/E A/E]"
    stop
  endif
  call getarg(2,soln_file1)
  call getarg(3,mascon_file2)
  call getarg(4,soln_file2)
  call getarg(5,output_file)

! PT200915: check whether we want apriori or estimated. These are optional command line arguments.
  which_soln = "E"         ! set the default to be "estimated"
  call getarg(6,arg)
  if(arg(1:1) == "A" .or. arg(1:1) == "Y")then
    print*,"File 1 mascon values to use: ", arg(1:1)
    which_soln(1) = arg(1:1)
  endif
  call getarg(7,arg)
  if(arg(1:1) == "A" .or. arg(1:1) == "Y")then
    print*,"File 2 mascon values to use: ", arg(1:1)
    which_soln(2) = arg(1:1)
  endif

! PT220608: 8th command line argument is the value by which to decimate the output ternaries
  arg = " "
  call getarg(8,arg)
  if(arg(1:1) == " ")then
    n_decimate = 1
  else
    read(arg,*)n_decimate
  endif

! PT220621: add a scale factor (or delta_time) to permit differences to be converted to rates
  arg = " "
  call getarg(9,arg)
  if(arg(1:1) == " ")then
    delta_time = 1.d0
  else
    read(arg,*)delta_time
  endif
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! open the files
  lumsc1 = 10
  lumsc2 = 11
  lusoln1 = 12
  lusoln2 = 13
  luout  = 14

  open(luout,file=output_file)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the first mascon solution file
  open(lusoln1,file=soln_file1)
  read(lusoln1,'(a26)')message   ! need to read the first line to determine the type of file

! VCV file
! V2 GRACEFIT   VCV solution:
  if(message(4:8) == "GRACE" .and. message(15:17) == "VCV")then
    call read_soln_v3('UTIL    ', lusoln1, .false., sTid, sScl, sBias, sEmp, sMasc, n_emp, nScl, nBias) 

    ! PT200302: save off the uncertainty of the mascons
    allocate(soln_sigma(total_params))
    do iprim=1,total_prim
      soln_sigma(iprim+sMasc-1) = dsqrt(VCV_obs_local(iprim+sMasc-1,iprim+sMasc-1))
    enddo

! FIT file
! V2 GRACEFIT   FIT solution
  else if(message(6:10) == "GRACE" .and. message(17:19) == "FIT")then  ! PT191228 changed to 17:19
    call status_update('FATAL','UTIL','diff_ternarys',' ',"Can't read GRACEFIT/GRACESIM.FIT files. Please write code yourself :-)",0)
  
! ADDNORM VCV file
! V2 ADDNORM   VCV solution
  else if(message(4:10) == "ADDNORM")then 
    ! PT190909: we need to allocate apriori, soln 
    allocate(apriori(total_params))
    allocate(soln(total_params))
    allocate(soln_sigma(total_params))
    ! PT190301: for now, this returns all params but not the details of the pointers to the different ICs. Mascons are fine though.
    call read_soln_addnorm('UTIL    ',message(14:18), 1, lusoln1, total_params, nparam, apriori, soln, prm_input &
                          , VCV_obs_local, .false. , sTid, sScl, sBias, sEmp, sMasc, n_emp, nScl, nBias) 

    ! PT200302: save off the uncertainty of the mascons. "nparam" has returned from read_soln_addnorm the number of mascons.
    do iprim=1,nparam
      soln_sigma(iprim+sMasc-1) = dsqrt(VCV_obs_local(iprim+sMasc-1,iprim+sMasc-1))
    enddo

! a priori VCV file generated by util/apriori_mascons 
! PT220316: or by h5_to_VCV.py
  else if (message(4:15) == "APRIORI  VCV" .or. message(1:26) == "VCV file written from hdf5")then
    read(lusoln1,*)nparam
    allocate(apriori(nparam))
    allocate(soln(nparam))
    ! and read in the a priori and the solution for the mascons.
    do i=1,4
      read(lusoln1,'(a)')line
    enddo
    do iprim=1,nparam
      read(lusoln1,'(31x,2f17.7)')apriori(iprim),soln(iprim)
    enddo
    sMasc = 1

    ! PT200302: there are no uncertainties for the mascons in the case of just an apriori file. Set to 0.1 mm
    allocate(soln_sigma(nparam))
    do iprim=1,nparam
      soln_sigma(iprim+sMasc-1) = 0.1d-3
    enddo

  else
    print*,'got to here message(1:26) ',message(1:26),' and code now stops!'
    stop
  endif
  close(lusoln1)
! read the first mascon file
  open(lumsc1,file=mascon_file1)
  call read_msc_hdr(lumsc1,mascon_file1,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                                ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
  call allocate_mascon_arrays
  call read_mascon_file(lumsc1,mascon_file1)  
  close(lumsc1)

! now, assign to each ternary mascon its primary mascon value. To do this we loop through all primarys and assign to each
! ternary within the primary
  do iprim = 1,total_prim
    do itern = 1, nint(mcon_prim(iprim,8))
      tern_num = mcon_tern(iprim,itern,7)
      soln1_ternarys(tern_num,1) = mcon_tern(iprim,itern,1)
      soln1_ternarys(tern_num,2) = mcon_tern(iprim,itern,2)
      ! PT200915: use either the a priori or the estimated solution, depending on what has been asked for
      if(which_soln(1) == "E")then
        soln1_ternarys(tern_num,3) = soln(iprim+sMasc-1)
      else if (which_soln(1) == "A")then
        soln1_ternarys(tern_num,3) = apriori(iprim+sMasc-1)
      endif
      soln1_ternarys(tern_num,4) = dble(iprim)
      soln1_ternarys(tern_num,5) = soln_sigma(iprim+sMasc-1)
    enddo
  enddo

! RM190403: deallocate arrays allocated within read_soln_v3
  if(message(4:8) == "GRACE" .and. message(15:17) == "VCV" )then
    deallocate(apriori)
    deallocate(soln)
    deallocate(soln_sigma)
    deallocate(vcv_obs)
!    deallocate(prmnam)
  endif
! PT190909: deallocate the arrays used in reading the addnorm files
  if(message(4:10) == "ADDNORM" .or. message(4:10) == "APRIORI")then
    deallocate(apriori)
    deallocate(soln)
    deallocate(soln_sigma)
  endif

! PT220316: deallocate if the input is from VCV files created from hdf5
  if(message(1:26) == "VCV file written from hdf5")then
    deallocate(apriori)
    deallocate(soln)
    deallocate(soln_sigma)
  endif
  
! PT190301: don't do this unless the second file is not the same file as the first
  if(mascon_file1 /= mascon_file2)call deallocate_mascon_arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the second mascon solution file
!  soln = 0.d0
  open(lusoln2,file=soln_file2)
  read(lusoln2,'(a)')message   ! need to read the first line to determine the type of file
! VCV file
! V2 GRACEFIT   VCV solution:
  if(message(4:8) == "GRACE" .and. message(15:17) == "VCV")then
    call read_soln_v3('UTIL    ', lusoln2, .false., sTid, sScl, sBias, sEmp, sMasc, n_emp, nScl, nBias)
    print*,'2nd file first 10 mascon solution values:',soln(sMasc:sMasc+10) 

    ! PT200302: save off the uncertainty of the mascons
    allocate(soln_sigma(total_params))
    do iprim=1,total_prim
      soln_sigma(iprim+sMasc-1) = dsqrt(VCV_obs_local(iprim+sMasc-1,iprim+sMasc-1))
    enddo

! FIT file
! V2 GRACEFIT   FIT solution
  else if(message(4:8) == "GRACE" .and. message(15:17) == "FIT")then
    call status_update('FATAL','UTIL','diff_ternarys',' ',"Can't read GRACEFIT/GRACESIM.FIT files. Please write code yourself :-)",0)
  
! ADDNORM VCV file
! V2 ADDNORM   VCV solution
  else if(message(4:10) == "ADDNORM")then 
    ! PT190909: we need to allocate apriori, soln 
    allocate(apriori(total_params))
    allocate(soln(total_params))
    allocate(soln_sigma(total_params))
    ! PT190301: for now, this returns all params but not the details of the pointers to the different ICs. Mascons are fine though.
    call read_soln_addnorm('UTIL    ',message(14:18), 1, lusoln2, total_params, nparam, apriori, soln, prm_input &
                          , VCV_obs_local, .false. , sTid, sScl, sBias, sEmp, sMasc, n_emp, nScl, nBias) 

    ! PT200302: save off the uncertainty of the mascons. "nparam" has returned from read_soln_addnorm the number of mascons
    do iprim=1,nparam
      soln_sigma(iprim+sMasc-1) = dsqrt(VCV_obs_local(iprim+sMasc-1,iprim+sMasc-1))
    enddo

! a priori VCV file generated by util/apriori_mascons
! PT220316: or by h5_to_VCV.py
  else if (message(4:15) == "APRIORI  VCV" .or. message(1:26) == "VCV file written from hdf5")then
    read(lusoln2,*)nparam
    allocate(apriori(nparam))
    allocate(soln(nparam))
    ! and read in the a priori and the solution for the mascons.
    do i=1,4
      read(lusoln2,'(a)')line
    enddo
    do iprim=1,nparam
      read(lusoln2,'(31x,2f17.7)')apriori(iprim),soln(iprim)
    enddo
    sMasc = 1

    ! PT200302: save off the uncertainty of the mascons
    allocate(soln_sigma(total_params))
    do iprim=1,nparam
      soln_sigma(iprim) = 0.1d-3
    enddo

  else
    write(message2,'(a,a,a)')'unknown vcv file type: ',message(4:10),' Code reads GRACEFIT, GRACESIM, ADDNORM files'
    call status_update('FATAL','UTIL','diff_ternarys',soln_file2,message,0)
  endif
  close(lusoln2)

! read the first mascon file - but only if it is not the same mascon file as mascon_file1
  if(mascon_file2 /= mascon_file1)then
    open(lumsc2,file=mascon_file2)
    call read_msc_hdr(lumsc2,mascon_file2,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                                ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
    call allocate_mascon_arrays
    call read_mascon_file(lumsc2,mascon_file2)  
    close(lumsc2)
  else
    call status_update('STATUS','UTIL','diff_ternarys',' ',"Mascon file is the same, so do not need to re-read it",0)
  endif


! now, assign to each ternary mascon its primary mascon value. To do this we loop through all primarys and assign to each
! ternary within the primary
  do iprim = 1,total_prim
    do itern = 1, nint(mcon_prim(iprim,8))
      tern_num = mcon_tern(iprim,itern,7)
      soln2_ternarys(tern_num,1) = mcon_tern(iprim,itern,1)
      soln2_ternarys(tern_num,2) = mcon_tern(iprim,itern,2)
      ! PT200915: use either the a priori or the estimated solution, depending on what has been asked for
      if(which_soln(2) == "E")then
        soln2_ternarys(tern_num,3) = soln(iprim+sMasc-1)
      else if(which_soln(2) == "A")then
        soln2_ternarys(tern_num,3) = apriori(iprim+sMasc-1)
      endif
      soln2_ternarys(tern_num,4) = dble(iprim)
      soln2_ternarys(tern_num,5) = soln_sigma(iprim+sMasc-1)
    enddo
  enddo
  call deallocate_mascon_arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Output the difference
  write(message,'(a,i5)')'Writing out solution differences on each ternary. Decimated by ',n_decimate
  call status_update('STATUS','UTIL','diff_ternarys',output_file,message,0)
  do itern=1,maxterns,n_decimate
    write(luout,'(2f10.5,5f15.2)')soln2_ternarys(itern,1:2)*180.d0/pi &
             ,(soln1_ternarys(itern,3)-soln2_ternarys(itern,3))*1.e3 / delta_time &  ! converts to rates if dt /= 1.d0
             ,soln1_ternarys(itern,3)*1.e3,soln2_ternarys(itern,3)*1.e3  &
             ,soln1_ternarys(itern,5)*1.e3,soln2_ternarys(itern,5)*1.e3

  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end




