  subroutine mask_outliers(calling_prog,use_outlier_mask,nepochs_t,date,use_obs)

! subroutine to read file outlier.mask and identify any epochs of KBR/LRI observations that should be eliminated from 
! the gracefit/gracesim inversion.
!
! P. Tregoning
! 15 October 2020
!
! MODS
! PT220518: increase the dimensioning of use_obs to 8 observations to include position and velocity
! PT220518: always mask out range acceleration if range rate is masked out
! PT220624: mask three extra range acceleration points either side of range rate data outage

  implicit none

! passed variables
  character(*),  intent(in)    :: calling_prog           ! name of calling program
  logical     ,  intent(in)    :: use_outlier_mask       ! logical to either use or not use the outlier mask file
  integer*4,     intent(in)    :: nepochs_t              ! number of epochs in the gracefit/gracesim solution
  integer*4,     intent(in)    :: date(3)                ! year, month, day of observations being processed
  logical  ,     intent(inout) :: use_obs(nepochs_t,8)   ! flag to use/not use KBR/LRI obs of range/range rate/range accel

! local variables
  integer*4,parameter   :: lumask=600                  ! unit number of file outlier.mask
  integer*4             :: ioerr
  character*3           :: instrument                  ! "KBR" or "LRI"
  character*2           :: obs_type                    ! "R_", "RR", "RA"
  integer*4             :: start_outlier,end_outlier   ! start/end of outlier segment (inclusive)
  character*250         :: line,message
  integer*4             :: obs_num                     ! 1: kbr; 2: kbrr; 3: kbra; 4: LR; 5: LRR; 6: LRA  (first 3: kbr, last 3: LRI)
  integer*4             :: tmpdate(3)                  ! date read temporarily from the file outliers.mask
  logical               :: found_outlier               ! set to true if at least one outlier is found on the day
  real(kind=8)          :: tmp_mask                    ! 0: don't use outlier.mask, 1: try to read outliers from the file

  found_outlier = .false.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(use_outlier_mask)then
    ! open the file containing outliers to mask
    open(lumask,file='outlier.mask',status='old',iostat=ioerr)
    if(ioerr /= 0)then
      call status_update('WARNING',calling_prog,'mask_outliers','outlier.mask',"Outlier mask file not found",0)
      return
    else
       call status_update('STATUS',calling_prog,'mask_outliers','outlier.mask',"Opened outlier mask file",0)
   endif
  else
    call status_update('STATUS',calling_prog,'mask_outliers',' ',"User did not request that outliers be masked",0)
    return
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the file to find any entries for the requested date
  ! header line
  read(lumask,'(a)')line

  ioerr = 0
  do while (ioerr == 0)
    read(lumask,*,iostat=ioerr,end=1000)tmpdate(1:3),instrument, obs_type, start_outlier,end_outlier
    if(ioerr == 0 .and. tmpdate(1) == date(1) .and. tmpdate(2) == date(2) .and. tmpdate(3) == date(3))then
      write(message,'(a,a3,a,a2,a,i6,a,i6,i8,2(a1,i2.2))')"Remove outlier of instrument ",instrument," observation type " &
                         , obs_type ," from epoch ",start_outlier," to",end_outlier, date(1),"-",date(2),"-",date(3)
      call status_update('STATUS',calling_prog,'mask_outliers','outlier.mask',message,0)

      ! set logical to true
      found_outlier = .true.

      if(instrument == "KBR")then
        if(obs_type == "R_")obs_num=1    ! KBR range
        if(obs_type == "RR")obs_num=2    ! KBR range-rate
        if(obs_type == "RA")obs_num=3    ! KBR range-accel
      else if (instrument == "LRI")then
        if(obs_type == "R_")obs_num=4    ! KBR range
        if(obs_type == "RR")obs_num=5    ! KBR range-rate
        if(obs_type == "RA")obs_num=6    ! KBR range-accel
! PT220518: add GPS observation types
      else if (instrument == "GPS")then
        if(obs_type == "PO")obs_num=7
        if(obs_type == "VE")obs_num=8
      else
        write(message,'(a)')"Unknown instrument type ",instrument
        call status_update('FATAL',calling_prog,'mask_outliers','outlier.mask',message,0)
      endif

      ! now mask out the observations
      if(start_outlier > nepochs_t)start_outlier = nepochs_t
      if(end_outlier > nepochs_t)end_outlier = nepochs_t

      use_obs(start_outlier:end_outlier,obs_num) = .false.
! PT220518: also mask out the range acceleration if the range rate was masked out
      if(obs_num == 2 .or. obs_num == 5)then

        ! PT220624: also mask out an additional three RA obs each side of a RR data outage - the numiercal differentiator needs it thus
        if(start_outlier > 3)then
          start_outlier = start_outlier - 5
        else
          start_outlier = 1
        endif
        if(end_outlier < nepochs_t - 5)then
          end_outlier = end_outlier - 5
        else
          end_outlier = nepochs_t
        endif
        
        use_obs(start_outlier:end_outlier,obs_num+1) = .false.
        write(message,'(a,a3,a,a2,a,i6,a,i6,i8,2(a1,i2.2))')"Also   outlier of instrument ",instrument," observation type " &
                         , "RA" ," from epoch ",start_outlier," to",end_outlier, date(1),"-",date(2),"-",date(3)
        call status_update('STATUS',calling_prog,'mask_outliers','outlier.mask',message,0)
      endif
      
  
    endif

  enddo

1000 continue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! print a message if no outliers were found on this day
  if(.not. found_outlier)then
    write(message,'(a,i8,2(a1,i2.2))')"No outliers found in file for date: ", date(1),"-",date(2),"-",date(3)
    call status_update('STATUS',calling_prog,'mask_outliers','outlier.mask',message,0)
  endif

  return
  end







