 subroutine read_soln_v1 (lu, VCV_flag )
!  subroutine read_soln_v1 (lu,N)  

! Program read_soln, reading the solution from gracefit to gracekal

! L. Lescarmontier
! 14 January 2013 from P. Tregoning 4 January 2013
! APP130304: Modified to use write_soln include file and to accept and return apriori
!            soln & VCV_obs as input parameters
! APP130321: Removed user-passed value N defining the number of parameters
!            Instead we are using nparam which is read from the input file.
! PT190312 : fixed up problems of multiple dimensioning of the number of parameters. Now all done using mod_subs/soln_mod.f90

   use soln_mod

  implicit none

! --PARAMETERS--
  
!  include '../includes/grace.param'
!  include '../includes/write_soln.h90'

!  integer*8 length
  
  integer*8 lu, i, j, nlines, supp_line, count, nvals, N, VCV_flag !nparam

! PT190312: these are now declared in mod_subs/soln_mod.f90
!  double precision apriori(maxparm), soln(maxparm), VCV_obs(maxparm,maxparm)
  character message*200, line*400
  
  count = 0

! skip three lines                                                                                             
  do i=1,3
    read(lu,'(a)')message
    print*, 'message' ,message
  enddo

! now read the number of parameters                                                                            
  read(lu,'(32x,i15)')nparam
  call status_update('STATUS','GRACEKAL','gracekal/read_soln',' ','Reading number of parameters',0)
  print*, 'number of parameters' ,nparam
! work out how many lines this will be, given that there are 25 values written per line                        
  nlines = int(nparam/25)
  if ( mod(nparam, 25) .ne. 0 ) then
    supp_line = 1
  else
    supp_line = 0
  endif


!! PT190312: now that we have the number of parameters, allocate the arrays
   maxparm = nparam
   allocate(apriori(maxparm))
   allocate(soln(maxparm))
   allocate(VCV_obs(maxparm,maxparm))


!!!! SOLUTION A PRIORI                  
! read through until we find the "SOLUTION A PRIORI" line                                                      
   message = " "
   do while ( message(1:8) .ne. "SOLUTION" )
     read(lu, '(a)', end = 1000) message
     if ( message(1:8) .eq. "SOLUTION") then
       if ( message(12:17) .eq. "PRIORI") then
         call status_update('STATUS','GRACEKAL','gracekal/read_soln',' ','Reading Solution A Priori',0)
         count = 0
         do i = 1, nlines
           read(lu, *) (apriori(count + j), j = 1, 25)
           count = count + 25
         enddo

! read any values hanging on the last line                                                                     
        if ( supp_line .eq. 1 ) then
          nvals = mod(nparam, 25)
          read(lu,*)(apriori(count+j),j=1,nvals)
        endif
      else
        call status_update('WARNING','GRACEKAL','gracekal/read_soln',' ','No SOLUTION A PRIORI line, set apriori to zero',0)
        apriori = 0.d0
        backspace(lu)
      endif
    endif

  end do

!!!  SOLUTION VECTOR                                                                                           
! now read the solution vector (in the gracefit vcv file this is apriori + adjustment  )                       
  call status_update('STATUS','GRACEKAL','gracekal/read_soln',' ','Reading Solution Vector',0)
  read(lu, '(a)') message
  count = 0
  do i = 1, nlines
    read(lu, *) (soln(count + j), j = 1, 25)
!    print*, i, count, (soln(count + j), j = 1, 25)            
    count = count + 25
  enddo
! read any values hanging on the last line                                                                     
  if ( supp_line .eq. 1 ) then
    nvals = mod(nparam, 25)
    read(lu, *) (soln(count + j), j = 1, nvals)
  endif

!!! VCV part !!!
! read the variance covariance matrix
! APP130321: Insert VCV flag option to skip reading of VCV matrix
  if (VCV_flag .ne. 0 )  then
    call status_update('STATUS','GRACEKAL','gracekal/read_soln',' ','Reading Observation VCV',0)
    read(lu, '(a)') message
    do i = 1, nparam
!      print*, 'reading VCV line' ,i, 'of' ,nparam
      read(lu, *) (VCV_obs(i, j), j = 1, nparam)
    enddo
  endif

  return

1000 write(message,'(a)')'Error reading gracefit VCV file. Line "SOLUTION A PRIORI" not found'
  call status_update('FATAL','GRACEKAL','gracekal/read_soln',' ',message,0)

  return
  end subroutine read_soln_v1
