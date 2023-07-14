  subroutine read_THR1B(starting_epoch,nepochs_t,csat,LUTH,thr,thrusters_off)

! reads the THR1B data for one satellite (only) and fills out an array to indicate whether epochs
! are affected by orientation thrusts or not.
!
! P. Tregoning
! 19 November 2015

    implicit none


!********************  Variable declarations ********************
    integer*4      :: nepochs_t                                      ! number of epochs (probably 86400/5)
    character*1    :: csat                                           ! satellite A or B
    integer*4      :: starting_epoch                                 ! first epoch (in GRACE seconds?)

    double precision, intent(out) :: thr(nepochs_t,16)        ! Thruster data for each satellite at every epoch
    real(kind=8),     intent(out) :: thrusters_off(nepochs_t) ! Data on whether thrusters are having any effect on accelerations

    integer*4      :: iepoch, kepoch ! Counters for epochs
    integer*4      :: i              ! Counter variable
    integer*4      :: dummy_ints(14) ! Dummy variable used to store unwanted data from THR1B
    character(2)   :: dummy_char     ! Dummy variable used to store unwanted data from THR1B
    character(150) :: message        ! Message printed when done reading file
    integer*4      :: ioerr          ! Standard error variable

! PT151118: I've had to add some declarations here so that I don't have
!           to use the gracefit.h90 declarations
    integer*4               :: LUTH
    integer*4,parameter     :: BUTTERWORTH_OFFSET=14
    real(kind=8),parameter  :: epoch_interval=5.0
    integer*4               :: end_loop
!****************************************************************

 
! Skip the header in file(s)
    call input_skipHeader(LUTH,'GRACE '//csat//' THR1B')

!*********************** READ THR1B FILES ***********************

! TG130627: set all values to 0 (end of array will remain that way)
    thr = 0
    ioerr = 0

    do iepoch = 1,nepochs_t  ! Upper bound. Should never reach this number (or even close)
      read(LUTH,*,iostat=ioerr,end=10)(thr(iepoch,i),i=1,2),dummy_char,dummy_char,(dummy_ints(i),i=1,14),&
                                              (thr(iepoch,i),i=3,16)

! Assert there was no error reading file
      if(ioerr /= 0) call status_update('FATAL','UTIL','read_THR1B',' ','Error reading THR1B input file',ioerr)

!  Only read in relevant thrusts (first nepochs_t epochs + the 70 second residue from butterworth filters)
      if(starting_epoch+(nepochs_t-1)*epoch_interval+70 < thr(iepoch,1)) exit

    enddo ! end of epoch loop

10    continue

    write(message,'(a,i6,a,a)')"read ",iepoch-1,' thrusts from THR1B_',csat
    call status_update('STATUS','UTIL','read_THR1B',' ',message,0)

!****************************************************************

    thrusters_off = 0
! Set thrusters_off according to data in the THR1B files
    kepoch = 1

    do iepoch = 1,nepochs_t*int(epoch_interval)
      if((thr(iepoch,1) == 0)) exit

! find the appropriate row to start making changes in thrusters_off
      do while ((starting_epoch+(kepoch-1)*epoch_interval < thr(iepoch,1)-70).and.(kepoch < nepochs_t*int(epoch_interval)))
        kepoch = kepoch + 1
      enddo

      if (kepoch == nepochs_t*int(epoch_interval)) exit
! There are 140s/5=28 epochs affected by each butterworth filter on each thrust (so the epochs k to k+BUTTERWORTH_OFFSET*2-1)
! PT151120: define the maximum of this loop to be either the width of the thrust or the end of the day (i.e. 17280)
      if(kepoch + BUTTERWORTH_OFFSET*2 - 1 > 17280)then
        end_loop = 17280
      else
        end_loop = kepoch + BUTTERWORTH_OFFSET*2 - 1
      endif
      do kepoch = kepoch,end_loop
        thrusters_off(kepoch) = 1
        if (kepoch == nepochs_t*int(epoch_interval)) exit
      enddo
! Set kepoch back to what it was (note that do loop goes one beyond what it computes)
      kepoch = kepoch - BUTTERWORTH_OFFSET*2
    enddo  ! end of epoch loop
!****************************************************************

! Close files, they should no longer be necessary
    close(LUTH)

    return
  end subroutine read_THR1B

