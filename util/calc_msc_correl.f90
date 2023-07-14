  program calc_msc_correl

! program to read a vcv file and calculate the correlations between all the mascons. For now, the program
! is set up to read an addnorm vcv file but this could be expanded to reading a gracefit version later on.
!
! P. Tregoning
! 27 September 2017

  implicit none
  
  character   :: infile*100,outfile*100
  integer*4   :: n_orbs
  integer*4   :: n_msc,n_IC
  real*8, allocatable :: VCV(:,:)
  real*8, allocatable :: msc_correl(:,:)
  real*8, allocatable :: ICs(:)

! counters
  integer*4   :: i,j

! unit numbers
  integer*4, parameter   :: luin=10,luout=11

! other stuff
  character   :: line*100

! decode command line
  call getarg(1,infile)
  call getarg(2,outfile)
  open(luin,file=infile,status='old')
  open(luout,file=outfile,status='unknown')

! read and ignore lines until we get to the PARAMETER line
  line = " "
  do while (line(1:10)/=" PARAMETER")
    read(luin,'(a)')line
  enddo

! now read through to the first mascon parameter. Count them so we know how many IC parameters there are
  line = " "
  n_IC = 0
  do while (line(6:9)/=". MC")
    read(luin,'(a)')line
    n_IC = n_IC + 1
  enddo
  n_IC = n_IC - 1
  print*,'There are ',n_IC,' IC parameters'

! read and ignore lines until we get to the VCV line
  line = " "
  do while (line(1:14)/="  VCV SOLUTION")
    read(luin,'(a)')line
  enddo


! backspace one line to find out how many mascons there are
  backspace(luin)
  backspace(luin)
  read(luin,'(a)')line
  read(line(10:14),*)n_msc

  print*,'There are ',n_msc,' mascons in the vcv file'

! skip the VCV line
  read(luin,'(a)')line

! ok, so now we need to read in the VCV matrix
  allocate(ICs(n_IC))
  allocate(VCV(n_msc,n_msc))
  allocate(msc_correl(n_msc,n_msc))
  msc_correl = 0.d0

  call status_update('STATUS','UTIL','calc_msc_correl',' ',"Reading VCV matrix",0)
! PT170928: don't do this! Tthe VCV in addnorm.vcv is written out with parameters in the order: mascons, tidal_amplitudes, ICs
!  do i=1,n_IC
!    read(luin,'(a)')line
!  enddo
  do i=1,n_msc
    read(luin,*)(VCV(i,j),j=1,i)  !,ICs(:)
  enddo
  do i=1,n_msc
    do j=1,i
      VCV(j,i) = VCV(i,j)
    enddo
  enddo
  call status_update('STATUS','UTIL','calc_msc_correl',' ',"Have read in the lower triangular VCV matrix",0)


! compute the correlations between the all mascons
  do i=1,n_msc
    if(mod(i,500) == 0)then
      write(line,'(a,i7)')'Correlations for mascon: ',i
      call status_update('STATUS','UTIL','calc_msc_correl',' ',line,0)
    endif

    do j=1,n_msc
      msc_correl(i,j) = VCV(i,j) / (dsqrt(VCV(i,i)) * dsqrt(VCV(j,j)) )
!print*,i,j,msc_correl(i,j)
!      msc_correl(j,i) = msc_correl(i,j)
    enddo
    write(luout,'(10421f15.11)')msc_correl(i,:)
  enddo


  end 
