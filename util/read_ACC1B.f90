  subroutine read_ACC1B(starting_epoch,nepochs_t,csat,LUACC,accobs)

! reads the ACC1B data for one satellite (only)
!
! P. Tregoning
! 19 November 2015

    implicit none


!********************  Variable declarations ********************
    integer*4      :: nepochs_t                                      ! number of epochs (probably 86400/5)
    character*1    :: csat                                           ! satellite A or B
    integer*4      :: starting_epoch                                 ! first epoch (in GRACE seconds?)

    real(kind=8),     intent(inout) :: accobs(nepochs_t,5) ! non-gravitational accelerations (XYZ)

    integer*4      :: iepoch, kepoch ! Counters for epochs
    integer*4      :: i              ! Counter variable
    integer*4      :: dummy_ints(14) ! Dummy variable used to store unwanted data from THR1B
    character(2)   :: dummy_char     ! Dummy variable used to store unwanted data from THR1B
    character(150) :: message        ! Message printed when done reading file
    integer*4      :: ioerr          ! Standard error variable

! PT151118: I've had to add some declarations here so that I don't have
!           to use the gracefit.h90 declarations
    integer*4               :: LUACC
    real(kind=8),parameter  :: epoch_interval=5.0
    real(kind=8)            :: junk(4)
    integer*4               :: j

!****************************************************************

 
! Skip the header in file(s)
    call input_skipHeader(LUACC,'GRACE '//csat//' ACC1B')

!*********************** READ ACC1B FILES ***********************

! TG130627: set all values to 0 (end of array will remain that way)
    accobs(:,1:3) = 0
    ioerr = 0

    do iepoch = 1,nepochs_t  ! Upper bound. Should never reach this number (or even close)
      read(LUACC,*,iostat=ioerr,end=10)accobs(iepoch,1),dummy_char(1:1),(accobs(iepoch,i),i=2,4)
! PT151119: we want to store only every 5th observation so that it matches with the GNV1B sampling. Therefore, throw away the next 4 observations
      do j=1,4
        read(LUACC,*,iostat=ioerr,end=10)junk(1),dummy_char,junk(2:4)
      enddo

! Assert there was no error reading file
      if(ioerr /= 0) call status_update('FATAL','UTIL','read_ACC1B',' ','Error reading ACC1B input file',ioerr)

    enddo ! end of epoch loop

10    continue

    write(message,'(a,i6,a,a)')"read ",iepoch-1,' acceleration obs from ACC1B_',csat
    call status_update('STATUS','UTIL','read_ACC1B',' ',message,0)

!****************************************************************


! Close files, they should no longer be necessary
    close(LUACC)

    return
  end subroutine read_ACC1B

