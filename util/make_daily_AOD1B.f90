  program make_daily_AOD1B

! program to do the following:
! 1. read in the daily spherical harmonic fields for the ocean and atmosphere,
! 2. read in the S1 and S2 atmospheric tide coefficients for each epoch (the sum is read in, not the individual tides)
! 3. remove the S2 atmospheric tide component from every coefficient for the epochs of the AOD1B (00 06 12 18 24)
! 4. add the non-tidal atmosphere to the ocean
! 5. voila, we have our dealiasing product for the day that doesn't include any atmospheric tides
!
! P. Tregoning
! 15 October 2016
!
! MODS
! PT190526: for RL06 the tides have been removed, so we just need to output in the format that we need. We can
!           perhaps simply subtract tide models with zero coefficients (ie do nothing).
! PT190526: On the other hand, RL06 is 3-hourly, so need to re-dimension the arrays here

  implicit none

  character  :: input_aod1b_1*100,input_aod1b_2*100,output_aod1b*100,atm_tide_coeffs*100
  character  :: message*250,line*100

! Epochss of AOD1B coefficients (in GRACE seconds)
  integer*4, allocatable :: GRACE_epochs(:)
  integer*4              :: date(5)
  real(kind=8)           :: jd,sec
  real(kind=8)           :: grace_start_mjd = 51544.5d0

! atmospheric tide coefficients
  integer*4  :: nepochs_atm_tide
  real(kind=8),allocatable :: Catm_tide(:,:,:),Satm_tide(:,:,:)

! AOD1B atmosphere coefficients
  integer*4  :: nepochs_atm
  real(kind=8),allocatable :: Catm(:,:,:),Satm(:,:,:)

! AOD1B ocean coefficients
  integer*4  :: nepochs_ocn
  real(kind=8),allocatable :: Cocn(:,:,:),Socn(:,:,:)

! AOD1B atm + ocean coefficients
  real(kind=8),allocatable :: Ctot(:,:,:),Stot(:,:,:)

! unit numbers
  integer*4 :: luaod_in_1,luaod_in_2,luaod_out,luatm_tide

! other stuff
  integer*4 :: ioerr,i,j,iepoch,nepoch
  integer*4 :: maxdeg_atm,maxdeg_ocn,maxdeg_atm_tide
  integer*4 :: tmpdeg,tmpord
  real*8    :: tmpC,tmpS
  character*4 :: RL_num
  integer*4 :: tmp_nepochs_atm,tmp_nepochs_ocn


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get command line arguments
  call getarg(1,input_aod1b_1)
  if(input_aod1b_1(1:1) == " ")then
    write(message,'(a,a)')"Runstring:  make_daily_AOD1B AOD1B_2010-09-09_X_05.asc AOD1B_2010-09-10_X_05.asc aod_2010-09-09.asc " &
                         ," atm_tide_all.hs"
    call status_update('FATAL','UTIL','make_daily_AOD1B',' ',message,0)
  endif

  call getarg(2,input_aod1b_2)
  call getarg(3,output_aod1b)
  call getarg(4,atm_tide_coeffs)
  call getarg(5,RL_num)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! open the files
  luaod_in_1 = 10
  luaod_in_2 = 11
  luaod_out = 12
  luatm_tide = 13
  open(luaod_in_1,file=input_aod1b_1,status='old')
  open(luaod_in_2,file=input_aod1b_2,status='old')
  open(luaod_out,file=output_aod1b,status='unknown')
! PT190526: allow "none" for the tide model, for RL06 where tides have been removed already
  if(RL_num /= "RL06")then
    open(luatm_tide,file=atm_tide_coeffs,status='old')
    ! read the header information to find out the maximum degree
    ! of the dealiasing and atm tide coefficients
    read(luatm_tide,*)maxdeg_atm_tide
  else
    maxdeg_atm_tide = 0
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read the header information to find out the maximum degree
!! of the dealiasing and atm tide coefficients
!  read(luatm_tide,*)maxdeg_atm_tide
!  line=" "

  do while (line(1:14) /= "MAXIMUM DEGREE")
    read(luaod_in_1,'(a)')line
  enddo
  read(line(32:35),*)maxdeg_atm
  write(luaod_out,*)maxdeg_atm
  maxdeg_ocn = maxdeg_atm
  write(message,'(a,i4,a,i4)')"Atm tide model degree:",maxdeg_atm_tide," AOD1B degree:",maxdeg_atm
  call status_update('STATUS','UTIL','make_daily_AOD1B',' ',message,0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dimension the arrays
! PT190526: 6-hourly for RL < 06, 3-hourly for RL06
  if(RL_num /= "RL06")then
    nepochs_atm = 5
    nepochs_atm_tide = 145
  else if (RL_num == "RL06")then
    nepochs_atm = 9
    nepochs_atm_tide = 0
  endif
  allocate( Catm(0:maxdeg_atm,0:maxdeg_atm,nepochs_atm) )
  allocate( Satm(0:maxdeg_atm,0:maxdeg_atm,nepochs_atm) )

  allocate( Cocn(0:maxdeg_atm,0:maxdeg_atm,nepochs_atm) )
  allocate( Socn(0:maxdeg_atm,0:maxdeg_atm,nepochs_atm) )

  allocate( Ctot(0:maxdeg_atm,0:maxdeg_atm,nepochs_atm) )
  allocate( Stot(0:maxdeg_atm,0:maxdeg_atm,nepochs_atm) )

  allocate( Catm_tide(0:maxdeg_atm_tide,0:maxdeg_atm_tide,0:nepochs_atm_tide-1) )
  allocate( Satm_tide(0:maxdeg_atm_tide,0:maxdeg_atm_tide,0:nepochs_atm_tide-1) )

  allocate(GRACE_epochs(nepochs_atm))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read in the atmospheric tide coefficients
  if(RL_num /= "RL06")then
    call status_update('STATUS','UTIL','make_daily_AOD1B',' ',"Reading in atm tide coefficients",0)
    do iepoch = 0,nepochs_atm_tide-1
      do i = 1,(maxdeg_atm_tide+1)*(maxdeg_atm_tide+2)/2
        read(luatm_tide,*)tmpdeg,tmpord,tmpC,tmpS
        Catm_tide(tmpdeg,tmpord,iepoch) = tmpC
        Satm_tide(tmpdeg,tmpord,iepoch) = tmpS
      enddo
    enddo
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read in the AOD1B data

! skip the rest of the header
  do while (line(1:13) /= "END OF HEADER")
    read(luaod_in_1,'(a)')line
  enddo

  ioerr = 0
  tmp_nepochs_atm = 0
  tmp_nepochs_ocn = 0
  do while (ioerr == 0)
    read(luaod_in_1,'(a)',iostat=ioerr,end=1000)line
    if(line(58:68) == "OF TYPE atm")then   ! this is a header line of a set of atm coefficieints
      tmp_nepochs_atm = tmp_nepochs_atm + 1
      call status_update('STATUS','UTIL','make_daily_AOD1B',' ',line,0)
! convert YMDHHMMSS to GRACE seconds
      read(line(38:56),'(i4,4(1x,i2),1x,f2.0)')date,sec
      call ymdhms_to_jd(date, sec, jd)
      grace_epochs(tmp_nepochs_atm) = ( (jd - 2400000.5d0) - grace_start_mjd)*86400.d0

      do i = 0, maxdeg_atm
        do j=0,i
          read(luaod_in_1,*)tmpdeg,tmpord,Catm(i,j,tmp_nepochs_atm),Satm(i,j,tmp_nepochs_atm)
        enddo
      enddo
    else if (line(58:68) == "OF TYPE ocn")then   ! this is a header line of a set of ocean coefficients
      tmp_nepochs_ocn = tmp_nepochs_ocn + 1
      call status_update('STATUS','UTIL','make_daily_AOD1B',' ',line,0)
      do i = 0, maxdeg_ocn
        do j=0,i
          read(luaod_in_1,*)tmpdeg,tmpord,Cocn(tmpdeg,tmpord,tmp_nepochs_ocn),Socn(tmpdeg,tmpord,tmp_nepochs_ocn)
        enddo
      enddo
    else  if (line(58:65) == "OF TYPE ")then
!      call status_update('STATUS','UTIL','make_daily_AOD1B',' ',line,0)
      do i=0,maxdeg_atm
        do j=0,i
          read(luaod_in_1,*)tmpdeg,tmpord,tmpC,tmpS
        enddo
      enddo
    endif

  enddo
1000 continue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the second AOD1B file to get the 24th hr epoch for the day
! that we want
! skip the header
  line = " "
  do while (line(1:13) /= "END OF HEADER")
    read(luaod_in_2,'(a)')line
  enddo

! extract out the atm and the ocn fields
  do while (tmp_nepochs_atm < nepochs_atm .or. tmp_nepochs_ocn < nepochs_ocn)
    read(luaod_in_2,'(a)',iostat=ioerr,end=1000)line
    if(line(58:68) == "OF TYPE atm")then   ! this is a header line of a set of atm coefficieints
      tmp_nepochs_atm =  tmp_nepochs_atm + 1
      call status_update('STATUS','UTIL','make_daily_AOD1B',' ',line,0)
! convert YMDHHMMSS to GRACE seconds
      read(line(38:56),'(i4,4(1x,i2),1x,f2.0)')date,sec
      call ymdhms_to_jd(date, sec, jd)
      grace_epochs(tmp_nepochs_atm) = ( (jd - 2400000.5d0) - grace_start_mjd)*86400.d0

      do i = 0, maxdeg_atm
        do j=0,i
          read(luaod_in_2,*)tmpdeg,tmpord,Catm(i,j,tmp_nepochs_atm),Satm(i,j,tmp_nepochs_atm)
        enddo
      enddo
    else if (line(58:68) == "OF TYPE ocn")then   ! this is a header line of a set of ocean coefficients
      tmp_nepochs_ocn = tmp_nepochs_ocn + 1
      call status_update('STATUS','UTIL','make_daily_AOD1B',' ',line,0)
      do i = 0, maxdeg_ocn
        do j=0,i
          read(luaod_in_2,*)tmpdeg,tmpord,Cocn(tmpdeg,tmpord,tmp_nepochs_ocn),Socn(tmpdeg,tmpord,tmp_nepochs_ocn)
        enddo
      enddo
    else  if (line(58:65) == "OF TYPE ")then
!      call status_update('STATUS','UTIL','make_daily_AOD1B',' ',line,0)
      do i=0,maxdeg_atm
        do j=0,i
          read(luaod_in_2,*)tmpdeg,tmpord,tmpC,tmpS
        enddo
      enddo
    endif
  enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! OK, so we now have all the data read in. We now need to 
! subtract the tide coefficients from each of the AOD1B
! atm coefficients at each epoch
!
! PT190526: don't do this for RL06
  if(RL_num /= "RL06")then
    do iepoch = 1,nepochs_atm
      Catm(:,:,iepoch) = Catm(:,:,iepoch) - Catm_tide(0:maxdeg_atm,0:maxdeg_atm,(iepoch-1)*36)
      Satm(:,:,iepoch) = Satm(:,:,iepoch) - Satm_tide(0:maxdeg_atm,0:maxdeg_atm,(iepoch-1)*36)
    enddo
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now we add the non-tidal atm coefficients to the ocean
! coefficients to generate and write out our dealiasing
! non-tidal product
  Ctot = Catm + Cocn
  Stot = Satm + Socn

  call status_update('STATUS','UTIL','make_daily_AOD1B',output_aod1b,"Writing non-tidal AOD1B file",0)
  do iepoch = 1,nepochs_atm
! write out the GRACE seconds for this epoch
    write(luaod_out,*)int(GRACE_epochs(iepoch))

! now the spherical harmonic coefficients
    do i=0,maxdeg_atm
      do j=0,i
        write(luaod_out,*)i,j,Ctot(i,j,iepoch),Stot(i,j,iepoch),"  0.000 0.000"
      enddo
    enddo
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  end

 
