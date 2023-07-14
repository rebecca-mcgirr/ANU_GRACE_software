  subroutine read_precnut

! subroutine to open and read the new file that contains nutation angles, dpsi, deps, and the precession/nutation matrix elements
! every 5 seconds for the required day of integration. Values are returned via the common variables in rotation_mod.f90
!
! P. Tregoning
! 24 February 2014

  use rotation_mod

  implicit none

  integer*4 maxepoch,i,j,k,junk
  character line*100
  integer*4 date(5)            ! date from the header of the input file (nutabl.info)
  real(kind=8) :: pi
  
  pi = 4.d0*datan(1.d0)

! open the file
  open(39,file='nutabl.info',status='old')

! read the date from the header line. We may add a check later that this date is appropriate for the requested integration span.
  read(39,'(73x,5i4)')date
  write(line,'(a,5i5)')" Read nutation/precession information for date",date
  call status_update('STATUS','GRACEORB','read_recnut',' ',line,0)

! now, read 17280 lines, being every 5 seconds for an entire day. Each line contains:
! seconds  dpsi  deps  rotation_matrix (1,1  1,2  1,3  2,1  2,2  2,3  3,1  3,2  3,3 )
  i = 0
  do while (i <= 86400/5)
    read(39,*)nut_epoch(i),dpsi(i),deps(i), ((rot_nutprec(j,k,i),k=1,3),j=1,3)
    i = i + 1
  enddo

  return
  end
   
