  subroutine read_acc_sca_thr(tin,efic)

! A complete rewrite of the subroutine accred.f90 to enable reading the ascii files of the L1B data, for any of the 
! space gravity missions.
!
! P. Tregoning
! 19 June 2018
!
! Emma-Kate Potter, 18 April, 2011 
! This subroutine reads in the accelerometer data spanning the 
! time range of the integration

! MODIFIED: APP 121004
! To allow user-defined names for the star camera and accelerometer files
! 
! PT180619: a substantial rewrite was undertaken to generalise this routine
! PT200827: add thermal corrections of low_freq cross-track signals onto along-track and radial

  use accred_mod
  use inmod_mod
  use inmod_mod     ! provides model options for removing non-linear effects in accelerometer obs
  use gm_mod        ! gm mod of the Earth
  use sat_mod       ! provides the 1-char name of the satellite. Use this to apply thermal corrections to along-track

  implicit none
  real(kind=8), intent(in)  :: tin          ! start of integration (in gracesecs)
  real(kind=8), intent(in)  :: efic(3)      ! earth-fixed position

  integer :: ACCnum, SCAnum 
  integer*4 :: i,j, ioerr  
  real(kind=8) :: temp1, temp2, temp3, temp4, temp5 
  character*8  :: calling_program

! PT180619: variables added to interact with the generic reads of L1B ascii data
  integer*4    :: n_acc, n_sca          ! number of accelerometer and star camera observations
  integer*4    :: n_sca_full            ! number of star camera records including epochs of missing data
  integer*4    :: acc_span(2)           ! start/end epoch (in gracesec) of the star camera data
  integer*4    :: mission               ! 0: GRACE, 1: GRACE FO, 2: GRACE II
  real(kind=8),allocatable :: acc_obs_tmp(:,:)  ! temporary storage of ACC data
  real(kind=8),allocatable :: sca_obs_tmp(:,:)  ! temporary storage of SCA data, concatenated over gaps
  integer*4    :: start_gap,end_gap     ! start and end epochs around a gap in the data
  integer*4    :: npoints               ! number of points to use either side of a gap
  real(kind=8) :: params(10)            ! parameters (that may be) needed to linearise the accelerometer obs
     
! PT180705: variables for reading the thrust file (to know when to exclude ACC1B obs)
  !character(80) :: thr_file                   ! THR1B files to be used
  integer*4     :: n_thrusts                  ! number of thrusts in the file
  integer*4     :: THRnum                     ! unit number for THR1B file
  integer*4,allocatable :: thr_obs(:,:,:)     ! array for the thrust observations
  integer*4     :: thrust_width               ! number of points to remove either side of a thrust event
  integer*4     :: start_flag,end_flag        ! temporary values of the epoch for the start/end of the thrust-affected data
  integer*4     :: ithrust                    ! loop counter for thrusts

! PT180809: variables related to the EMD decomposition
  integer*4     :: start_emd,end_emd          ! start/end components for the reconstruction of the signal
  logical       :: debug
  integer*4     :: ndecomp                    ! maximum number of components for the EMD decomposition
  integer*4     :: icomp                      ! loop counter for accelerometer axes
  integer*4     :: iepoch                     ! loop counter for epochs
  real(kind=8),allocatable :: acc_decomp(:,:) ! temporary array to hold the components of the EMD decomposition

! RM190628: variables related to the FFT and cos filter
  integer*4                :: nfft
  real(kind=4)             :: rad             ! radius from IC positions
  real(kind=4)             :: T               ! period of orbit in seconds
  real(kind=4)             :: flo             ! frequencies below flo are removed
  real(kind=4)             :: fhi             ! frequencies above fhi are kept
  real(kind=4),allocatable :: pspec(:)        ! power spectrum
  real(kind=4),allocatable :: fspec(:)        ! frequencies
  real(kind=8),allocatable :: acc_fft_out(:)  ! output from acc_fft
  real(kind=8),allocatable :: acc_diff(:)     ! diff acc_obs and acc_fft

! local counters and stuff
  integer*4    :: irow
  integer*4    :: acc_model_start,acc_model_end  ! start/end epochs, relative to start of acc obs, of the requested acc obs to be modelled

! variables to hardwire the addition of a saw-tooth model for day 2016-07-26
  real*8 :: acc_extra(86400,2)

! PT200819: variables to scale low_freq cross-track signals to along-track and radial thermal corrections
  real(kind=8) :: thermal_scl_factors(3,2)   
  real(kind=8),allocatable :: acc_CT(:)        ! backup copy of the original cross-track observations  

  ACCnum=61
  SCAnum=60
  THRnum=62
  calling_program = "GRACEORB"
! PT180621: set the mission variable to "-1" so that the program will determine the type of data from the data format itself
  mission = -1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT200819: set the cross-track to along-track scaling factors
  ! GRACE A
  thermal_scl_factors(1,1) = -0.042d0   ! along-track thermal scale factor for GRACE A
  
  ! GRACE B
  ! PT220429: the transfer value for GRACE B depends on whether it is transplant data or not. Decide this based on the epoch
  if(tin < 531144000.d0)then   ! before end of October 2016 it is GRACE-B non-transplant data
    thermal_scl_factors(1,2) = -0.017d0   ! along-track thermal scale factor for GRACE B
  else if (tin >= 531144000.d0 .and. tin < 546955200.d0)then   ! Nov'16 to 2 May 2017 it is transplant data for GRACE B so use the GRACE A transfer function
    thermal_scl_factors(1,2) = -0.042d0   ! along-track thermal scale factor for GRACE B
  else if (tin >= 547041600.d0 .and. tin < 549460800.d0)then   ! May'17 is NOT transplant                                                        
    thermal_scl_factors(1,2) = -0.017d0   ! along-track thermal scale factor for GRACE B
  else if (tin >= 549460800.d0 .and. tin < 552139200.d0)then   ! June'17 is again transplant                                                        
    thermal_scl_factors(1,2) = -0.042d0   ! along-track thermal scale factor for GRACE B
  else    ! for GRACE-FO don't use anything
    thermal_scl_factors(1,2) = 0.d0       ! along-track thermal scale factor for GRACE B
  endif                                                            
                                                              
  thermal_scl_factors(3,1) =  0.d0      ! radial thermal scale factor for GRACE A
  thermal_scl_factors(3,2) =  0.d0      ! radial thermal scale factor for GRACE B
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! open L1B files for accelerometer (ACC) and star camera (SCA) and thruster (THR)
  call level1B_open(ACCnum,calling_program,ACC_file)
  call level1B_open(SCAnum,calling_program,SCA_file)
  call level1B_open(THRnum,calling_program,THR_file)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Read thrust data 
! read the THR1B header
  call thr_read_hdr(THRnum,calling_program,thr_file,mission,n_thrusts)

! now we can dimension the arrays to hold the thrust observations.
  allocate(thr_obs(1,n_thrusts,8))              ! gracesec of times of thrusts (we don't care which thruster fired, just that something fired) 

! read the THR1B linear thrust data
  call thr_read_data(THRnum,calling_program,thr_file,mission,n_thrusts,n_thrusts,thr_obs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Read accelerometer data 
  call acc_read_hdr(ACCnum,calling_program,ACC_file,mission,acc_span )
! the number of observations can be deduced from the span of the data
  n_acc = acc_span(2) - acc_span(1) + 1
  allocate(acc_obs(1,n_acc,5))
  call acc_read_data(ACCnum,calling_program,ACC_file,mission,n_acc,3,acc_span,acc_obs(1,:,:) )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! flag accelerometer observations within 120 seconds of a thrust to eliminate obs that
! may be affected by the CRN filtering of the square thrusts
  acc_obs(1,:,5) = 1.0
  do ithrust = 1,n_thrusts
    ! remove points "N" seconds before the thrust
    thrust_width = 60
    start_flag = (thr_obs(1,ithrust,1)-acc_span(1))-thrust_width
    if( (thr_obs(1,ithrust,1)-thrust_width) <= acc_span(1) )start_flag = 1    ! check that we don't start before the data

    ! remove points 60 seconds after the thrust
    end_flag = (thr_obs(1,ithrust,1)-acc_span(1)) + thrust_width
    if( (thr_obs(1,ithrust,1)+thrust_width) > acc_span(2) ) end_flag = (acc_span(2)-acc_span(1)+1)    ! check that we don't end after the data
    ! set the flags to false
    acc_obs(1,start_flag:end_flag,5) = 0.0
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! PT180705: remove the relevant model(s) from the different components to get rid of any non-linear features
! PT180702: now apply the relevant model to remove non-linear behaviour in the accelerometer obs
      ! the user input values (acc_fit_start(i), acc_fit_end(i)) control how much of the accelerometer data is
      ! to be used when fitting the model. Trouble is, these values are in epochs of day, whereas the accelerometer
      ! data is in GRACE seconds - and may start before the start of the day. So we need to convert the user input
      ! epochs into a counter relative to the start of the accelerometer obs, not the start of the day.

! PT200827: make a copy of the original cross-track obs. We need this for the thermal model corrections to AT and Rad
  allocate(acc_CT(n_acc))
  acc_CT = acc_obs(1,:,3)

  do icomp=1,3
    acc_model_start = (nint(tin) - acc_span(1))/5 + acc_fit_start(icomp)
    acc_model_end   = (nint(tin) - acc_span(1))/5 + acc_fit_end(icomp)
! ----------------------------------------------
    if( acc_fit_model(icomp)(1:9) == "quadratic")then
      call acc_fit_quadratic(.false.,calling_program,n_acc,acc_obs(1,:,icomp+1),acc_obs(1,:,5) &
             ,acc_model_start,acc_model_end,params(1:3))
      write(message,'(a,3e20.8)')"  Model parameter values (accel, rate, offset) (nm/s^2):",params(1:3)*1.e9
      call status_update('STATUS',calling_program,'read_acc_sca_thr',' ',message,0)
! ----------------------------------------------


! ----------------------------------------------
    else if( acc_fit_model(icomp)(1:11) == "exponential")then
      call acc_fit_exponential(.true.,calling_program,n_acc,acc_obs(1,:,icomp+1),acc_obs(1,:,5) &
             ,acc_model_start,acc_model_end,params(1:4))
      write(message,'(a,5e20.8)')"  Exponential model parameters  (nm/s^2):",params(1:4)*1.e9
      call status_update('STATUS',calling_program,'read_acc_sca_thr',' ',message,0)
! ----------------------------------------------


! ----------------------------------------------
! PT180809: use the EMD to deompose and remove the long-wavelength components (that are the temperature-related effects ?)
    elseif(acc_fit_model(icomp)(1:4) == "EMD1")then
      ndecomp    = 13
      debug = .false.
      allocate(acc_decomp(ndecomp,n_acc))
      call acc_EMD1(debug,calling_program,icomp,n_acc,ndecomp,nint(acc_obs(1,n_acc,icomp+1)-acc_obs(1,1,icomp+1)) &
                    ,acc_obs(1,:,icomp+1),acc_decomp,end_emd)

      ! reconstruct with only the components that don't contain the non-linear signal
      do iepoch = 1,n_acc
         acc_obs(1,iepoch,icomp+1) = sum(acc_decomp(1:end_emd-1,iepoch))
      enddo
      write(message,'(a,i1,a,i3,a,i3,a)')"  Axis ",icomp,": removed long-wavelength components",end_emd," to",ndecomp &
                                          ," of EMD filter."
      call status_update('STATUS','GRACEORB','read_acc_sca_thr',' ',message,0)
      deallocate(acc_decomp)
! ----------------------------------------------

! ----------------------------------------------
! RM180914: use the EMD2 to deompose and remove the long-wavelength components (that are the temperature-related effects ?)
    elseif(acc_fit_model(icomp)(1:4) == "EMD2")then
      ndecomp    = 13
      debug = .false.
      allocate(acc_decomp(ndecomp,n_acc))
      call acc_EMD2(debug,calling_program,icomp,n_acc,ndecomp,nint(acc_obs(1,n_acc,icomp+1)-acc_obs(1,1,icomp+1)) &
                    ,acc_obs(1,:,icomp+1),acc_decomp,end_emd)

      ! reconstruct with only the components that don't contain the non-linear signal
      do iepoch = 1,n_acc
         acc_obs(1,iepoch,icomp+1) = sum(acc_decomp(1:end_emd-1,iepoch))
      enddo
      write(message,'(a,i1,a,i3,a,i3,a)')"  Axis ",icomp,": removed long-wavelength components",end_emd," to",ndecomp &
                                          ," of EMD2 filter."
      call status_update('STATUS','GRACEORB','read_acc_sca_thr',' ',message,0)
      deallocate(acc_decomp)
! ----------------------------------------------

! ----------------------------------------------
! RM181203: use the EMD3 to deompose and remove the long-wavelength components (that are the temperature-related effects ?)
    elseif(acc_fit_model(icomp)(1:4) == "EMD3")then
      ndecomp    = 13
      debug = .false.
      allocate(acc_decomp(ndecomp,n_acc))
      call acc_EMD3(debug,calling_program,icomp,n_acc,ndecomp,nint(acc_obs(1,n_acc,icomp+1)-acc_obs(1,1,icomp+1)) &
                    ,acc_obs(1,:,icomp+1),acc_decomp,end_emd)

      ! reconstruct with only the components that don't contain the non-linear signal
      do iepoch = 1,n_acc
         acc_obs(1,iepoch,icomp+1) = sum(acc_decomp(1:end_emd-1,iepoch))
      enddo
      write(message,'(a,i1,a,i3,a,i3,a)')"  Axis ",icomp,": removed long-wavelength components",end_emd," to",ndecomp &
                                          ," of EMD3 filter."
      call status_update('STATUS','GRACEORB','read_acc_sca_thr',' ',message,0)
      deallocate(acc_decomp)
! ----------------------------------------------

! ----------------------------------------------
! PT180809: use the EMD to deompose and remove the long-wavelength components (that are the temperature-related effects ?)
    elseif(acc_fit_model(icomp)(1:5) == "ACC_Z" .and. icomp == 3)then
      ndecomp    = 13
      debug = .false.

! PT180912: the GRACE A accelerometer data in July 2016 needs a saw-tooth model added to it so that the obs capture the total
!           maximum SRP signal. Here I will code a hard-wired version of a model (linear increase + exponential decay) to see
!           whether it makes any difference
      open(342,file='model_acc_2016_07_26_A_Z.dat',status='old')
      do i=1,86400
        read(342,*)acc_extra(i,:)
      enddo
      ! now add the model to the obs, matching up on GRACESEC ....
      do iepoch=1,86400
        ! just add it on, starting at epoch 86401
         acc_obs(1,86400+iepoch,icomp+1) = acc_obs(1,86400+iepoch,icomp+1) - acc_extra(iepoch,2)
      enddo

      allocate(acc_decomp(ndecomp,n_acc))
      call acc_EMD1(debug,calling_program,icomp,n_acc,ndecomp,nint(acc_obs(1,n_acc,icomp+1)-acc_obs(1,1,icomp+1)) &
                    ,acc_obs(1,:,icomp+1),acc_decomp,end_emd)

      ! reconstruct with only the components that don't contain the non-linear signal
      do iepoch = 1,n_acc
         acc_obs(1,iepoch,icomp+1) = sum(acc_decomp(1:end_emd-1,iepoch))
      enddo
      write(message,'(a,i1,a,i3,a,i3,a)')"  Axis ",icomp,": removed long-wavelength components",end_emd," to",ndecomp &
                                          ," of EMD filter."
      call status_update('STATUS','GRACEORB','read_acc_sca_thr',' ',message,0)
      deallocate(acc_decomp)
! ----------------------------------------------


! ----------------------------------------------
    elseif(acc_fit_model(icomp)(1:3) == "fft")then
      debug = .false.
      nfft = 2**20
      ! calculate the period of the orbit
      rad = (efic(1)**2 + efic(2)**2 + efic(3)**2)**0.5
      T = ((4*pi**2*rad**3)/(gm(1)))**0.5 
      ! set low and high frequency cutoffs
! PT/RMcG: set these to specific frequencies
!      flo = 1.0/(T*1.5) ! freq of 2/3 per rev signal (i.e 1/~135 cycles/min)
!      fhi = 1.0/(T*1.25) ! freq of once per rev signal (i.e. 1/~90 cycles/min)
! PT/RMcG190927: these frequencies were erroneously set to d-9 instead of d-3. Changed on 2019-09-27
      !flo = 0.09d-3
      !fhi = 0.11d-3
! RMcG200512: changed values to 0.45d-4 and 0.55d-4
      flo = 0.45d-4
      fhi = 0.55d-4
      allocate(acc_fft_out(n_acc))

      call acc_cos_fft(debug,calling_program,icomp,n_acc,nfft,acc_obs(1,:,icomp+1),acc_fft_out,flo,fhi) 
      acc_obs(1,:,icomp+1) = acc_fft_out

      write(message,'(a,i1,a,e12.5,a)')"  Axis ",icomp,": filtered low frequency components below ",fhi," cycles/sec"
      call status_update('STATUS','GRACEORB',calling_program,' ',message,0)
      deallocate(acc_fft_out)
! ----------------------------------------------


! ----------------------------------------------
    elseif(acc_fit_model(icomp)(1:5) == "therm")then
      debug = .false.
      nfft = 2**20
      ! calculate the period of the orbit
      rad = (efic(1)**2 + efic(2)**2 + efic(3)**2)**0.5
      T = ((4*pi**2*rad**3)/(gm(1)))**0.5 
      ! set low and high frequency cutoffs
! PT/RMcG: set these to specific frequencies
!      flo = 1.0/(T*1.5) ! freq of 2/3 per rev signal (i.e 1/~135 cycles/min)
!      fhi = 1.0/(T*1.25) ! freq of once per rev signal (i.e. 1/~90 cycles/min)
! PT/RMcG190927: these frequencies were erroneously set to d-9 instead of d-3. Changed on 2019-09-27
      !flo = 0.09d-3
      !fhi = 0.11d-3
      flo = 0.45d-4
      fhi = 0.55d-4
      allocate(acc_fft_out(n_acc))
      allocate(acc_diff(n_acc))

! run the fft on the cross-track component
      call acc_cos_fft(debug,calling_program,icomp,n_acc,nfft,acc_CT,acc_fft_out,flo,fhi) 
      acc_diff = acc_CT-acc_fft_out

! RM/PT200819: apply the relevant scaling from cross-track to other component, depending on which satellite
      if(sat == "A")then
        acc_obs(1,:,icomp+1) = acc_obs(1,:,icomp+1)-(acc_diff*thermal_scl_factors(icomp,1)) 
      elseif(sat == "B")then
        acc_obs(1,:,icomp+1) = acc_obs(1,:,icomp+1)-(acc_diff*thermal_scl_factors(icomp,2))
      endif

      write(message,'(a,i1,a)')"  Axis ",icomp,": applied thermal-based correction from cross-track low_freq thermal signal"
      call status_update('STATUS','GRACEORB',calling_program,' ',message,0)
      deallocate(acc_fft_out)
      deallocate(acc_diff)
! ----------------------------------------------


! ----------------------------------------------
    elseif(acc_fit_model(icomp)(1:4) == "none")then
      ! do nothing
      write(message,'(a,i2)')"  No model removed from accelerometer obs for axis",icomp
      call status_update('STATUS',calling_program,'read_acc_sca_thr',' ',message,0)
! ----------------------------------------------

! ----------------------------------------------
    else
      write(message,'(a,a,a)')"Model `",acc_fit_model(icomp),"' not coded. Please add it yourself"
      call status_update('FATAL',calling_program,"read_acc_sca_thr",' ',message,0)
    endif
! ----------------------------------------------
  enddo

  ! PT180621: for now, just transfer the values to the previous arrays, so that the code should run
  allocate(ACCtime(n_acc))
  allocate(accr(1:3,1:n_acc))

  ACCtime = acc_obs(1,:,1)
  do iepoch=1,n_acc
    do icomp=1,3
      accr(icomp,iepoch) = acc_obs(1,iepoch,icomp+1)
    enddo
    !print*,"DEBUG:",iepoch,accr(:,iepoch)
  enddo
  ACCn = n_acc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Read star camera data 
  call sca_read_hdr(SCAnum,calling_program,SCA_file,mission,n_sca,sca_step )
  allocate(sca_obs_tmp(n_sca,6))
! sca_obs contains:
! gracesec, 4xquaternions  star_camera_flag
  call sca_read_data(SCAnum,calling_program,SCA_file,mission,n_sca, sca_obs_tmp)
  sca_step = sca_obs_tmp(2,1)-sca_obs_tmp(1,1)
  do iepoch=3,n_sca
     sca_step = min(sca_step, sca_obs_tmp(iepoch,1)-sca_obs_tmp(iepoch-1,1))
  enddo
! DEBUG
! jump across the infill and just use as read
!  allocate(sca_obs(n_sca,6))
!  sca_obs = sca_obs_tmp
!  goto 9876


! update the number of star camera observations to now include missing data
  n_sca_full = nint( (sca_obs_tmp(n_sca,1)-sca_obs_tmp(1,1))/sca_step) + 1

! PT180620: now, fill out the SCA data, leaving zero entries where there are data gaps (to be filled in by interpolation later)
  allocate(sca_obs(n_sca_full,6))
  sca_obs = 0.d0
  sca_obs(1,:) = sca_obs_tmp(1,:)
  do i=2,n_sca
    irow = nint( (sca_obs_tmp(i,1)-sca_obs_tmp(1,1))/sca_step) + 1
    sca_obs( irow,:) = sca_obs_tmp(i,:)
  enddo

  write(message,'(a,i10)')" Number of star camera observations (including infilled gaps): ",n_sca_full
  call status_update('STATUS','GRACEORB','read_acc_sca',' ',message,0)

! now, fill in any missing values using a 2nd-order quadratic interpolation
  n_sca = n_sca_full
  do i=1,n_sca-1
    if(nint(sca_obs(i,1)) == 0)then      ! it is a missing epoch
      if(i == 1)then
        sca_obs(i,:) = sca_obs(i+1,:)    ! just make it all the same as epoch 2
      else
        ! using 800 points either side of the gap, fit a quadratic to the data and use modelled values to infill missing obs
        ! first, find the end of the data gap
        start_gap = i
        end_gap = i
        do while (end_gap < n_sca .and. nint(sca_obs(end_gap,1)) == 0)
          end_gap = end_gap + 1
        enddo
        ! now adjust by 800 epochs before/after
        npoints = 100
        if (start_gap > npoints)then
          start_gap = start_gap - npoints
        else
          start_gap = 1
        endif
        if(int((n_sca - end_gap)/sca_step) > npoints)then
          end_gap = end_gap + npoints
        else
          end_gap = n_sca
        endif
        ! fit the quadratic model to the span of data
        call infill_TS(.false.,calling_program,"quadratic ","SCA",end_gap-start_gap+1,5,sca_obs(start_gap:end_gap,:),sca_step)
      endif
    endif
  enddo

! PT180621: for now, store the new data in the names of the old arrays, just to get the code working
9876 continue
  allocate(SCAtime(1:n_sca))
  allocate(varr(1:4,1:n_sca))

  SCAtime = sca_obs(:,1)
  do i=1,n_sca
    do j=1,4
      varr(j,i) = sca_obs(i,j+1)
    enddo
  enddo

  STARn = n_sca
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    return
    end subroutine read_acc_sca_thr
