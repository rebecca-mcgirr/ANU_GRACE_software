 program model_acc_erin

! program to compute a model (deg 2 polynomial) to flatten out the non-linear shapes that we see in the
! uncalibrated accelerometer data, especially in the cross-track
!
! P. Tregoning
! 15 May 2017
!
! PT170724: output bogus 999 values if the observations are missing at any particular epoch. GRACEORB will then
!           have to deal with it.

  implicit none

  character*100 :: accfiles(2),accfiles_out(2),thrfiles(2)
  integer*4     :: yr,month,day, PLT_ORIG
  integer*4     :: luacc(2),luthr(2),luacc_out(2)
  integer*4     :: poly_degree        ! maximum degree for polynomial to fit to the data
  integer*4     :: bit_map            ! 1: X, 2: Y, 4: Z (i.e. 7 = all three components)
  logical       :: bitmap               ! function to decipher the bit-mapped value

! Level1B  data
  real(kind=8)  :: acc_data(90000,4,2)   ! ACC1B define it for 90000 seconds (more than necessary) for 4 components (xyz+flag), 2 satellites
  real(kind=8)  :: thr_data(18000,2)     ! THR1B define it for 18000 epochs (more than necessary) for epoch of thrusts, 2 satellites
  integer*4     :: gracesec
  integer*4     :: start_gracesec(2)

! temp variables for reading accelerometer values
  real(kind=8)  :: tmpsec,tmpacc(3)
  integer*4     :: dt    , q                ! seconds between time of a thrust and the graceseconds of the start of an orbit

! least squares variables
  real(kind=8),allocatable :: A(:,:),B(:,:),At(:,:),W(:,:),AtW(:,:),AtWA(:,:), VCV(:,:),AtWB(:,:),soln(:,:)  
  integer*4     :: nepochs,nobs,nparam,iter
  real(kind=8)  :: t                     ! time (in decimal days) since the start of the accelerometer file
  real(kind=8)  :: apriori(3,3)          ! apriori values for offset, rate and acceleration terms of LS polynomial model

! counters and character variables
  integer*4     :: iepoch,ithrust,ioerr,isat,i,axis
  character*100 :: line,arg
  character     :: csat*2
  character*256 :: message

! parse command line arguments
  call getarg(1,arg)
  if(arg(1:1) == "")then
    call status_update('FATAL','UTIL','model_acc',' ',"Runstring: model_acc 2010 09 07 poly_degree bit-mapped-value",0)
  else
    read(arg,*)yr
    call getarg(2,arg)
    read(arg,*)month
    call getarg(3,arg)
    read(arg,*)day
  endif
  call getarg(4,arg)
  read(arg,*)poly_degree
  call getarg(5,arg)
  read(arg,*)bit_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! define file names
! accelerometer files
  write(accfiles(1),'(a,i4,2(a1,i2.2),a)')"ACC1B_",yr,'-',month,'-',day,'_A_02.asc'
  write(accfiles(2),'(a,i4,2(a1,i2.2),a)')"ACC1B_",yr,'-',month,'-',day,'_B_02.asc'
  write(thrfiles(1),'(a,i4,2(a1,i2.2),a)')"THR1B_",yr,'-',month,'-',day,'_A_02.asc'
  write(thrfiles(2),'(a,i4,2(a1,i2.2),a)')"THR1B_",yr,'-',month,'-',day,'_B_02.asc'

  write(accfiles_out(1),'(a,i4,2(a,i2.2),a)')"ACC_A_",yr,'-',month,'-',day,'.asc'
  write(accfiles_out(2),'(a,i4,2(a,i2.2),a)')"ACC_B_",yr,'-',month,'-',day,'.asc'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! open the files
! ACC1B
  luacc(1) = 10
  luacc(2) = 11
  open(luacc(1),file=accfiles(1),status='old')
  open(luacc(2),file=accfiles(2),status='old')

! THR1B
  luthr(1) = 12
  luthr(2) = 13
  open(luthr(1),file=thrfiles(1),status='old')
  open(luthr(2),file=thrfiles(2),status='old')

! output files
  luacc_out(1) = 14
  luacc_out(2) = 15
  open(luacc_out(1),file=accfiles_out(1),status='unknown')
  open(luacc_out(2),file=accfiles_out(2),status='unknown')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read all the ACC1B data
  line = " "
  do isat = 1,2
    do while (line(1:13) /= "END OF HEADER")
      read(luacc(isat),'(a)')line
    enddo

! get the starting epoch in grace seconds
    read(luacc(isat),*,iostat=ioerr,end=123)start_gracesec(isat)
123 if(ioerr /= 0)then
      write(message,'(a)')"It appears that there are no accelerometer data!"
      call status_update('FATAL','UTIL','model_acc',accfiles(isat),message,0)
    else
      backspace(luacc(isat))
    endif
open(PLT_ORIG, FILE='acceleration.txt')

! read the accelerometer data and store in the appropriate epoch in the matrix (there may be holes if data are missing)
!    call status_update('STATUS','UTIL','model_acc',accfiles(isat),"Reading accelerometer data",0)
    ioerr = 0
    nepochs = 0
    acc_data(:,4,isat) = 999.d0              ! set the flag a priori that ALL data are "bad"
    do while (ioerr == 0)
      read(luacc(isat),'(a)',end=1000,iostat=ioerr)line
      if(ioerr == 0)then
        nepochs = nepochs+1
        read(line,*)tmpsec,csat,tmpacc(1:3)
write(PLT_ORIG,*) tmpacc(2)
        iepoch = (tmpsec - start_gracesec(isat) )+1
        acc_data(iepoch,1:3,isat) = tmpacc
        acc_data(iepoch,4,isat) = 0.d0       ! set the flag that the accelerometer data for this epoch are "good"
      endif
    enddo
1000 continue
  enddo
 close(PLT_ORIG)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read all the THR1B data
  line = " "
  do isat = 1,2
!    call status_update('STATUS','UTIL','model_acc',thrfiles(isat),"Reading thrust data",0)
    do while (line(1:13) /= "END OF HEADER")
      read(luthr(isat),'(a)')line
    enddo

! read the accelerometer data
    ioerr = 0
    ithrust = 0
    do while (ioerr == 0)
      read(luthr(isat),'(a)',end=1001)line
      if(ioerr == 0)then
        ithrust = ithrust + 1
        read(line,*)thr_data(ithrust,isat)
! set to "bad" the epochs +/- 60 seconds of the thrust event
        dt = thr_data(ithrust,isat) - start_gracesec(isat)
        if(dt < 15)then
          acc_data(1:dt,4,isat) = 99.d0
        else if (86400 - dt < 60)then
          acc_data(dt:86400,4,isat) = 99.d0
        else
          acc_data(dt-15:dt+15,4,isat) = 99.d0
        endif

      endif
    enddo
1001 continue
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !print*,'nepochs and ithrust',iepoch,ithrust
print*, 'ERINS VERSION!'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the least squares. Solve for offset, rate and quadratic terms
  nparam = 3
  allocate(A(nepochs,nparam))
  allocate(At(nparam,nepochs))
!We don't use it  allocate(W(nepochs,nepochs))
  allocate(AtW(nparam,nepochs))
  allocate(AtWA(nparam,nparam))
  allocate(VCV(nparam,nparam))
  allocate(B(nepochs,1))
  allocate(AtWB(nparam,1))
  allocate(soln(nparam,1))
q = 1
! loop over satellites
  do isat = 1 ,2
   call status_update('STATUS','UTIL','model_acc',' ',"Least squares inversion: fill out A matrix",0)
    do axis = 1,3
      apriori(:,axis) = 0.d0
      nobs = 0
      do iter=1,4
        nobs = 0
        A = 0.d0
        B = 0.d0
        do iepoch = 1, 86400
!if(mod(iepoch,100) == 0)print*,'iepoch:',iepoch
          
if (q == 1) then     ! it is a valid observation that is unaffected by thrusts
            nobs = nobs+1
            t = dble(iepoch)/86400.d0
            A(nobs,1) = t**poly_degree
            A(nobs,2) = t**(poly_degree-1)
            A(nobs,3) = 1.d0
         
            B(nobs,1) = acc_data(nobs,axis,isat) - ( apriori(1,axis)*t**poly_degree &
		+ apriori(2,axis)*t**(poly_degree - 1) + apriori(3,axis)*t**(poly_degree - 2) )       ! just do the cross-track for now ....
          endif
        enddo

! least squares inversion
 !     call status_update('STATUS','UTIL','model_acc',' ',"Least squares inversion: transpose A matrix",0)
        call transp(A,At,nepochs,nparam)
        call matmult(At,A,AtWA,nparam,nepochs,nparam)
  !    call status_update('STATUS','UTIL','model_acc',' ',"Least squares inversion: invert AtWA",0)
        call invert(AtWA,VCV,nparam,nparam)
 
   !  call status_update('STATUS','UTIL','model_acc',' ',"Least squares inversion: AtWB",0)
        call matmult(At,B,AtWB,nparam,nepochs,1)
        call matmult(VCV,AtWB,soln,nparam,nparam,1)

     ! print*,'isat, axis, xhat: ',isat,axis,soln(:,1)
        apriori(:,axis) = apriori(:,axis) + soln(:,1)
      print*,'soln(:,1)',soln(:,1)
      enddo
    enddo
!print*,'isat,apriori(:,2)',isat,apriori(:,2)
! now, output the observed - modelled values
! PT170522: use the input bit-mapped value to determine which axes to update. Do this by setting the estimated
!           value (stored in "apriori") to zero if we don't want to update
    if( .not. bitmap(bit_map,1))then
      apriori(:,1) = 0.d0
    else
      write(message,'(a,i3)')" Remove model from accelerometer obs for X-axis for sat: ",isat
      call status_update('STATUS','UTIL','model_acc',' ',message,0)
    endif
    if(  .not. bitmap(bit_map,2))then
      apriori(:,2) = 0.d0
    else
      write(message,'(a,i3)')" Remove model from accelerometer obs for Y-axis for sat: ",isat
      call status_update('STATUS','UTIL','model_acc',' ',message,0)
    endif
    if(  .not. bitmap(bit_map,3))then
      apriori(:,3) = 0.d0
    else
      write(message,'(a,i3)')" Remove model from accelerometer obs for Z-axis for sat: ",isat
      call status_update('STATUS','UTIL','model_acc',' ',message,0)
    endif

    do iepoch = 1,86400
      if (acc_data(iepoch,4,isat) < 999.d0) then     ! it is a valid observation 
        t = dble(iepoch)/86400.d0
! output, correcting for only acceleration and rate (leave the offset bias in there .....)
        write(luacc_out(isat),*)start_gracesec(isat)+iepoch-1 &
           ,(acc_data(iepoch,axis,isat) -  (apriori(1,axis)*t**poly_degree +apriori(2,axis)*t**(poly_degree-1)) &
		,axis=1,3) &
           ,acc_data(iepoch,1:3,isat)

      else
! PT170724: write out the missing epochs as just huge observations
        write(luacc_out(isat),*)start_gracesec(isat)+iepoch-1, 999.0,999.0,999.0
	print*, luacc_out
      endif
    enddo

  enddo



  end
  
