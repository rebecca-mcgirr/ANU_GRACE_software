  program norm_diag

! program to extract the diagonal elements for the mascons out of the entire binary normal equations 
!
! A quick and dirty program for a specific, quick problem. Not meant to have a long lifetime and not written with structure and longevity in mind.

! P. Tregoning
! 29 August 2018

  implicit none

  integer*4 :: k,j,nparams,lufile,nmascons,nICs,nmsc_tidal,rec_len,year,month,day,hr,minute,irec
  real(kind=8),allocatable :: tmpvals(:),apr_tmp(:)
  real(kind=8) :: sec,duration,version,junk(24)
  character*150 :: normeq_file

! local variables
  real(kind=8),allocatable :: msc_diag(:),AtWb(:)

  call getarg(1,normeq_file)

  lufile = 10

  open(lufile,file=normeq_file,status='old',form='unformatted',access='direct',recl=24)
  read(lufile,rec=1)version,nparams,nmascons,nICs,nmsc_tidal
  close(lufile)
  rec_len = nparams*8
  open(lufile,file=normeq_file,status='old',form='unformatted',access='direct',recl=rec_len)


  allocate(tmpvals(nparams))
  allocate(apr_tmp(nparams))
  allocate(msc_diag(nmascons))
  allocate(AtWb(nmascons))


!  line 1 contains version, # parameters, date and duration of orbit integration
  read(lufile,rec=1)version,nparams,nmascons,nICs,nmsc_tidal,year,month,day,hr,minute,sec &
       ,duration
  print *, "line 1:",version,nparams,nmascons,nICs,nmsc_tidal,year,month,day,hr,minute,sec &
       ,duration

! line 2 of version 1.0 files are the a priori values of the parameters
  read(lufile,rec=2)(apr_tmp(j),j=1,nparams)

! now, read the normal equations. Just the lines of the mascons
  do k=1,nparams
     irec = 5+k
     read(lufile,rec=irec)(tmpvals(j),j=1,nparams)
     if(k > 24 .and. k < nmascons+24)then
      msc_diag(k-24) = tmpvals(24+k)
      ! write(*,*)k,tmpvals(24+k)
     endif
  enddo

! now read the AtWb
  read(lufile,rec=irec+1)(junk(j),j=1,nICs),(AtWb(j),j=1,nmascons)


! now, output various values to see which one will tell me something about the magnitudes of the mascon adjustments
  do k=1,nmascons
    write(*,*)k,AtWb(k)/msc_diag(k),msc_diag(k),AtWb(k)
  enddo

  end



