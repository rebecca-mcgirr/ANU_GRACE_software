  program sum_components

  implicit none

  character*100 :: infile, arg
  integer*4 :: n_excl,counter,ioerr,icomp,i,j
  real*8 :: values(13,200000),sum_total(200000),sum_use(200000)
  integer*4 :: excl(13),epoch(200000)
  logical use_comp(13)

  use_comp = .true.

! command arguments
  call getarg(1,infile)
  call getarg(2,arg)
  read(arg,*)n_excl

  do i=1,n_excl
    call getarg(2+i,arg)
    read(arg,*)excl(i)
    use_comp(excl(i)) = .false.
  enddo

!do i=1,13
!  print*,"Component states: ",use_comp(i),i
!enddo

! read the file
  counter = 1
  open(10,file=infile,status='old')
  ioerr = 0
  do while (ioerr == 0)
    read(10,*,iostat=ioerr,end=1000)epoch(counter),(values(j,counter),j=1,13)
!print*,counter,epoch(counter),(values(j,counter),j=1,13)
    if(ioerr == 0)counter = counter+1
  enddo

1000 continue

! ok, add them all up and add up only the ones that we want
  sum_total = 0.d0
  sum_use = 0.d0
  do i=1,counter-1
    do icomp=1,13
      sum_total(i) = sum_total(i) + values(icomp,i)
!print*,i,icomp,values(icomp,i),sum_total(i)
      if(use_comp(icomp)) then
!print*,'i,icomp. Use component',i,icomp
        sum_use(i) = sum_use(i) + values(icomp,i)
      endif
    enddo
  
    write(*,*)epoch(i),sum_total(i),sum_use(i)
  enddo



  end

