  program atm_tide_model

! program to take the atmospheric tide constituent files associated with the AOD1B product and generate an atmospheric tide model
! file, with values every 10 minutes, for digestion into graceorb.
!
! There are 12 tidal consituents that have been removed from the 3-hourly AOD1B fields. Restoring these 12 constituents, therefore, will create
! a complete atmospheric mass variation model. The tide models are stored in the 12 separate files as degree 180 spherical harmonic models of 
! the cosine and sine components of the tides.
!
! The 12 tidal constituents were sourced from ftp://isdcftp.gfz-potsdam.de/grace/Level-1B/GFZ/AOD/RL06/TIDES/ as advised by Henryk Dobeslaw on
! 5 February 2020.
!
! P. Tregoning
! 6 February 2020

  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the number of tidal constituents
  integer*4,parameter      :: ntides=12

! the size of the spherical harmonic models
  integer*4                :: maxdeg,maxord

! arrays to store the tidal models
  real(kind=8),allocatable :: tides(:,:,:,:,:)      !  (tide number, deg, order, cos/sine)
  real(kind=8),allocatable :: sum_tides(:,:,:,:)    ! the sum of all tides every 10 minutes (epoch, deg, order, cos/sine)

! variables related to the tidal files
  character*100            :: tide_dir
  character*100            :: tide_files(ntides)
  character*2              :: tide_codes(ntides)
  data tide_codes/'P1','S1','K1','N2','M2','L2','T2','S2','R2','T3','S3','R3'/
  real(kind=8)             :: tide_freq(ntides)                ! tidal frequencies in degrees/hour
  data tide_freq/ 14.9589314d0,15.0d0      ,15.0410686d0 &
                 ,28.4397295d0,28.9841042d0,29.5284789d0 &
                 ,29.9589333d0,30.0d0      ,30.0410667d0 &
                 ,44.9589300d0,45.0d0      ,45.0410700d0 /
   
  character*2,parameter    :: RL="06"
  integer*4                :: interval                         ! required output interval (in seconds) of computed tides

! counters, loop indices
  real(kind=8)             :: t                                ! time (in seconds) since start of day
  integer*4                :: iepoch,itide,icomp,ideg,iord
  character*200            :: message
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! allocate the tide variables
  maxdeg = 180
  maxord = maxdeg
  allocate(tides(ntides,0:maxdeg,0:maxord,2,2))
  allocate(sum_tides(0:144,0:maxdeg,0:maxord,2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write some information to the screen to start the program off
  write(*,100)"atm_tide_model: will read the tide files provided with the RL06 AOD1B model and" &
             ," output a file of spherical harmonic coefficients to represent the sum of all " &
             ," tidal consituents" &
             ,"The tidal constituents are: K1, L2, M2, N2, P1, R2, R3, S1, S2, S3, T2, T3"
100 format(//,3a,//,a,//)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set up the names of the tide files
  tide_dir = "/Gdata/GRAV/GRACE/data/podacc_L1B/GFZ/AOD1B/RL06/TIDES/"
  do itide = 1, ntides
    tide_files(itide) = trim(tide_dir)//"AOD1B_ATM_"//tide_codes(itide)//"_"//RL//".asc"
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read in each tide constituent cosine/sine spherical harmonic model
  do itide = 1,ntides
    call read_atm_tide_constituent(maxdeg,tide_files(itide),tides(itide,:,:,:,:))
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! evaluate each constituent cosine/sine every 10 minutes for 24 hours, starting at time t=0
  interval = 600                 ! 10 minute output interval for tide evaluations
  sum_tides = 0.d0
print*,'calculate tides on ',86400/interval + 1,' epochs'
  do iepoch = 0,86400/interval
    t = dble(iepoch*interval)
    write(message,'(a,f10.1)')" Summing constituents for epoch: ",t
    if(mod(t,6000.d0) == 0.d0)call status_update('STATUS','UTIL','atm_tide_model',' ',message,0)
    do itide=1,ntides
      call calc_atm_tide_constituent(tide_codes(itide),tide_freq(itide),t,maxdeg,tides(itide,:,:,:,:),sum_tides(iepoch,:,:,:))
!print*,iepoch,itide,sum_tides(iepoch,0,0,1)
    enddo
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output the file for use in graceorb. Output will be a spherical harmonic model for the cosine
! and then sine component of the sum of the tides, every 10 minutes starting at 00UT
  call output_sum_tides(86400/interval + 1,maxdeg,sum_tides)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call status_update('STATUS','UTIL','atm_tide_model',' ',"End of ATM_TIDE_MODEL",0)
  end

