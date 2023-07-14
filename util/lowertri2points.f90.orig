  program lowertri2points

! program to read a lower-triangular correlation matrix and output it as a series of points that can then be plotted

  implicit none

  integer*4 :: npoints,junk
  real*8,allocatable :: values(:,:)
  character*100 :: infile,line
  integer*4     :: irow,icol,luin


! get the input file from the command line
  call getarg(1,infile)

! get the total number of parameters
  call getarg(2,line)
  read(line,*)npoints

! allocate arrays
  allocate(values(npoints,npoints))

! open the file
  luin = 10
  open(luin,file=infile,status='old')

! skip the first line
  read(luin,'(a)')line

! read in the lower-triangular matrix
  do irow = 1,npoints
    read(luin,*)junk,(values(irow,icol),icol=1,irow)
  enddo
!  print*,'have read in the lower-triangular matrix'

! output it one point per line
  do irow=1,npoints
    do icol=1,irow
      write(*,'(i7,i7,f12.4)')irow,icol,values(irow,icol)
    enddo
  enddo

  end

