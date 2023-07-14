  program norm_diag_ascii

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

  open(lufile,file=normeq_file,status='old')
  read(lufile,*)version,nparams,nmascons,nICs,nmsc_tidal

  allocate(tmpvals(nparams))
  allocate(apr_tmp(nparams))
  allocate(msc_diag(nmascons))
  allocate(AtWb(nmascons))


! now, read the normal equations. Just the lines of the mascons
print*,'reading the normal equations'
  do k=1,nparams
     read(lufile,*)(tmpvals(j),j=1,nparams)
     if( k <= nmascons)then
      msc_diag(k) = tmpvals(k)
      ! write(*,*)k,tmpvals(24+k)
     endif
  enddo

! now read the AtWb
print*,'reading AtWb'
  do k=1,nmascons
    read(lufile,*)junk(1),AtWb(k)
  enddo

! now, output various values to see which one will tell me something about the magnitudes of the mascon adjustments
  do k=1,nmascons
    write(*,*)k,dabs(AtWb(k)/msc_diag(k)),msc_diag(k),AtWb(k)
  enddo

  end



