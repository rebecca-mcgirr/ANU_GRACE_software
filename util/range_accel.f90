  program range_accel

! program to calculate the range acceleration (rather than range rate) from a set of range rate residuals (either pre- or post-fit)
!
! P. Tregoning
! 20 October 2016
!
! PT170726: don't output the first two points at the start/end nor just before/after any gap in the data.

  implicit none

  character :: infile*100,outfile*100 , resid_type*7 ,message*200,line*100

! number of observations
  integer*4 :: maxobs
  real(kind=8),allocatable :: range_rate(:),range_acc(:)

! observations from input file
  integer*4    :: epoch
  real(kind=8) :: postfit_RR,prefit_RR
  real(kind=8),allocatable :: crds(:,:),rpyA(:,:),rpyB(:,:)

! unit numbers
  integer*4,parameter :: luin=10,luout=11

! local variables
  integer*4 :: iobs



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! decode command line arguments
  call getarg(1,infile)
  if(infile(1:1) == " ")then
    write(message,'(a)')"Runstring: range_accel kb_file output_file prefit/postfit"
    call status_update('FATAL','UTIL','range_accel',' ',message,0)
  endif

  call getarg(2,outfile)
  call getarg(3,resid_type)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! open the files and get maxobs from input file
  open(luin,file=infile,status='old')
  open(luout,file=outfile,status='unknown')

  read(luin,'(a)')line
  write(luout,'(a)')line
  read(line,'(18x,i8)')maxobs
print*,'maxobs=',maxobs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! allocate the arrays
  allocate(range_rate(maxobs))
  allocate(range_acc(maxobs))
  allocate(crds(3,maxobs))
  allocate(rpyA(3,maxobs))
  allocate(rpyB(3,maxobs))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the kbrr file to get the residuals
! skip the header lines
  do iobs=1,5
    read(luin,'(a)')line
    write(luout,'(a)')line
  enddo

! read the obs
  do iobs = 1,maxobs
    read(luin,*)epoch,postfit_RR,crds(:,iobs),rpyA(:,iobs),rpyB(:,iobs),prefit_RR
    if(resid_type == "prefit ")then
      range_rate(iobs)=prefit_RR
    elseif(resid_type == "postfit")then
      range_rate(iobs)=postfit_RR
    else
      write(message,'(a,a,a)')'residual type ',resid_type,' unknown. Must be "prefit " or "postfit"'
      call status_update('FATAL','UTIL','range_accel',' ',message,0)
    endif
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate the derivative (which is the range acceleration)
  call noise_robust_deriv(range_rate*1.e3,range_acc,5.d0,maxobs,5)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write out the epoch, range rate and range acceleration
  call status_update('STATUS','UTIL','range_accel',outfile,"Writing out range accelerations",0)
! PT170726: don't output the first two and last two points. We should also remove two points before/after any gap in the data ....
  do iobs = 3,maxobs-3
    write(luout,'(i7,f16.10,2f12.5,f16.5,6f9.4,f21.10)')iobs,range_acc(iobs)*1.d0,crds(:,iobs),rpyA(:,iobs),rpyB(:,iobs) &
                    ,range_acc(iobs)*1.d0
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  end

