  program xyz_to_llr

! convert GRACE/GRACE-FO position coords to spherical lat/lon/radius given a GNV1B file
!
! R. McGirr
! 10 December 2020

  implicit none

  character*100               :: arg
  integer*4                   :: year, month, day    
  integer*4                   :: mission
  character*150               :: gnv_file
  real(kind=8),allocatable    :: llr(:,:)

! Data read in from GNV1B files
  real(kind=8), allocatable   :: gvec(:,:,:)    ! Positions and velocities for each satellite at every epoch
  real(kind=8), allocatable   :: gvec_epochs(:) ! Epochs of positions and velocities for each satellite at every epoch
  integer*4                   :: num_gps_used   ! number of GPS observations that will be used
  real(kind=8)                :: gnv_step       ! epoch step size for each GNV1B file
  integer*4                   :: gnv_span       ! span of epochs of each GNV1B file ? 
  integer*4                   :: epoch_interval
  real(kind=8)                :: starting_epoch

  real(kind=8)                :: pi
  integer*4                   :: nepochs_t
  integer*4                   :: iepoch
  integer*4                   :: ioerr
  integer*4                   :: tmp_int
  integer*4,parameter         :: lugnv=10, lullr=11
  character*20                :: calling_prog 
  character*150               :: out_file
  integer*4                   :: date(5)
  

  calling_prog = "xyz_to_llr"
  date = 0
  pi = 4.d0*datan(1.d0)

! decode command line
  call getarg(1,arg)
  if(arg(1:1) == "")then
    print*,"Runstring: xyz_to_llr year month day mission"
    stop
  endif
  read(arg,*)year
  call getarg(2,arg)
  read(arg,*)month
  call getarg(3,arg)
  read(arg,*)day
  ! mission
  call getarg(4,arg)
  read(arg,*)mission

  ! calculate (in GRACE seconds) the starting epoch at 00UT for this day
  date(1) = year
  date(2) = month
  date(3) = day
  call ymdhms_to_gracesec(date,0.d0,tmp_int)
  starting_epoch = dble(tmp_int)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! form up the names of the GNV1B files
  if (mission == 0)then
    write(gnv_file,'("GNV1B_",i4,"-",i2.2,"-",i2.2,"_A_02.asc")')year,month,day
  else
    write(gnv_file,'("GNV1B_",i4,"-",i2.2,"-",i2.2,"_C_04.txt")')year,month,day
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! open the GNV1B files
  open(lugnv,file=gnv_file,status='old',iostat=ioerr)
  if(ioerr /= 0)then
    call status_update('FATAL','UTIL','xyz_to_llr',gnv_file,"Error opening GNV1B file",0)
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the GNV1B files. Use the standard subroutines for reading L1B files
  call gnv_read_hdr(lugnv,calling_prog,gnv_file,mission,num_gps_used,gnv_step,gnv_span)

! allocate the array of epochs for gvec (unused in this subroutine but required when calling gnv_read_data)
  nepochs_t = int(gnv_span / gnv_step)

! Data read in from GNV1B files
  allocate(gvec(6,1,nepochs_t))
  allocate(gvec_epochs(nepochs_t))

  epoch_interval = 5.0
  call gnv_read_data(lugnv,calling_prog,gnv_file,mission,6,1,nepochs_t,1,starting_epoch &
                     ,gnv_step,num_gps_used,1.d0,gvec_epochs,nepochs_t,gvec)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! allocate array for llr converison
  allocate(llr(nepochs_t,3))
  call status_update('STATUS','UTIL','xyz_to_llr','',"Calculating longitude, latitude and radius from positions",0)

  do iepoch=1,nepochs_t
  ! Do the conversion to calculate the radius
    llr(iepoch,3) = dsqrt(gvec(1,1,iepoch)**2 + gvec(2,1,iepoch)**2 + gvec(3,1,iepoch)**2)

  ! latitude
    llr(iepoch,1) = dasin(gvec(3,1,iepoch)/llr(iepoch,3))*180.d0/pi

  ! longitude
    llr(iepoch,2) = datan(gvec(2,1,iepoch)/gvec(1,1,iepoch))*180.d0/pi

  ! get the right quadrant
    if(gvec(1,1,iepoch) > 0.d0 .and. gvec(2,1,iepoch) < 0.d0)then
      llr(iepoch,2) = 360.d0+llr(iepoch,2)
    else if (gvec(1,1,iepoch) < 0.d0 .and. gvec(2,1,iepoch) < 0.d0)then
      llr(iepoch,2) = 180.d0 + llr(iepoch,2)
    else if (gvec(1,1,iepoch) < 0.d0 .and. gvec(2,1,iepoch) > 0.d0)then
      llr(iepoch,2) = 180.d0 + llr(iepoch,2)
    endif
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write out epoch lat lon rad
  if (mission == 0)then
    write(out_file,'("LLR_",i4,"-",i2.2,"-",i2.2,"_A.asc")')year,month,day
  else
    write(out_file,'("LLR_",i4,"-",i2.2,"-",i2.2,"_C.asc")')year,month,day
  endif
  open(lullr,file=out_file,status='unknown')
  write(lullr,'(a)')"# latitude longitude radius"

  do iepoch=1,nepochs_t
    if (llr(iepoch,3) /= 0.d0) write(lullr,'(f12.6,f12.6,f16.6)')llr(iepoch,1),llr(iepoch,2),llr(iepoch,3)
  enddo
  call status_update('STATUS','UTIL','xyz_to_llr',out_file,'Latitude longitude radius written to file',0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end
