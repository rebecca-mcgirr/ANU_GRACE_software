  program transplant
  
! try to work out the time separation and, therefore, the roll/pitch/yaw corrections required to be made to the transplanted accelerometer
! obs for the GRACE-FO mission
!
! P. Tregoning
! 8 March 2022

  use gracefit_mod    ! brings in the "mission" variable
  
  implicit none

! command line arguments
  integer*4     :: date(5)
  character     :: arg*100,message*250

! Data read in from GNV1B files
  double precision, dimension(2) :: gnv_step     ! epoch step size for each GNV1B file
  double precision, allocatable  :: epochs(:)     ! matrix of angles between LOS vectors
  double precision, allocatable  :: gvec(:,:,:)     ! Positions and velocities for each satellite at every epoch
  double precision, allocatable  :: gvec_epoch(:,:) ! Epochs of positions and velocities for each satellite at every epoch

  integer*4:: n_gnv            ! span of epochs of each GNV1B file ?
  integer*4:: num_epochs_used  ! span of epochs of each GNV1B file ?

  integer*4, dimension(2) :: num_gps_used        ! number of GPS observations that will be used (after decimating the 5sec data)
  integer*4, dimension(2) :: gnv_span            ! span of epochs of each GNV1B file ?


! ******************************************************************
! decode runstring
  call getarg(1,arg)

  if ( arg(1:1) == " ") then
    write(message,'(a)')"Runstring: transplant 2019 10 02 mission "
    call status_update('FATAL','UTIL','transplant',' ',message,0)
  else
!   Do Nothing
  end if

! epoch: date(1) = year, date(2) = month, date(3) = day
  date = 0

  read(arg,*) date(1)

  call getarg(2,arg)
  read(arg,*) date(2)

  call getarg(3,arg)
  read(arg,*) date(3)

! mission
  call getarg(4,arg)
  read(arg,*) mission

  write(6,'(a,i5,2i3,/)') "Input date: ", date(1), date(2), date(3)
! ******************************************************************


! ******************************************************************
! form up the names of the GNV1B files
  if ( mission == 0 ) then
    write(gnv(1),'("data/GNV1B_",i4,"-",i2.2,"-",i2.2,"_A_02.asc")') date(1), date(2), date(3)
    write(gnv(2),'("data/GNV1B_",i4,"-",i2.2,"-",i2.2,"_B_02.asc")') date(1), date(2), date(3)
  else
    write(gnv(1),'("data/GNI1B_",i4,"-",i2.2,"-",i2.2,"_C_04.txt")') date(1), date(2), date(3)
    write(gnv(2),'("data/GNI1B_",i4,"-",i2.2,"-",i2.2,"_D_04.txt")') date(1), date(2), date(3)
  end if

  write(6,'(3(a,/))') "Input files: ", gnv(1), gnv(2)
! ******************************************************************

  end
  
