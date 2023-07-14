  program apriori_model_timeseries

! program to read in the low-pass filter a priori model file and output a time series at a single lat/lon location
!
! P. Tregoning
! 25 March 2022

  use mascon_mod          ! defines all the mascon arrays

  implicit none
  
! command line arguments
  character*150 :: mascon_file        ! file containing mascon geometry
  character*150 :: msc_model_file     ! file containing the quadratic+annual model parameters for the corresponding mascons
  real(kind=8)  :: user_lat, user_lon ! coords of requested point on which to generate the time series
  integer*4     :: model              ! type of model required (bit-mapped values)
  
! mascon file and user mascon geometry variables
  integer*4                :: tern_number
  real(kind=8),allocatable :: mcon_prim_usr(:,:)
  integer*4                :: total_prim_usr
  logical                  :: same_msc_geom
  real(kind=8)             :: lat_spacing
  integer*4                :: model_msc_number
  
  ! variables to read in mascon metric information
  real(kind=8),allocatable :: msc_crds(:,:)          ! mascon longitude/latitude
  real(kind=8),allocatable :: msc_annual(:,:)        ! cosine/sine amplitudes of annual signal for each mascon
  real(kind=8),allocatable :: msc_wrms(:,:)          ! the wrms values (**list them)
  real(kind=8),allocatable :: msc_monthly(:,:)       ! monthly low-frequency signal values from mascon time series
  real(kind=8),allocatable :: epochs(:)

  integer*4                :: nmonths                ! number of monthly values of low-freq signal
  integer*4                :: nmsc                   ! number of mascons in the WRMS file
  integer*4                :: nrms                   ! number of RMS values in file (between annual ampl and low-freq values)
  integer*4                :: ntot                   ! total number of columns
  character*100            :: line

! variables to interpolate the low-frequency msc metrics to the solution epoch
  real(kind=8)             :: dt,dt2
  integer*4                :: i1,i2,imonth,iRMS
  real(kind=8)             :: lowfrq

! external functions
  logical                  :: bitmap

! local variables
  real(kind=8)  :: dec_yr           ! decimal year of required epoch
  integer*4     :: ioerr
  integer*4     :: imsc
  real(kind=8)  :: pi,w
  character*110 :: arg,message
  integer*4     :: trimlen
  character*15  :: calling_prog
  real(kind=8)  :: msc_apr_value
  real(kind=8)  :: soln_epoch
  
! unit numbers
  integer*4, parameter :: lumsc=10, lumsc_model=11,lumsc_usr=12,luout=13
  
  pi = 4.d0*datan(1.d0)
  calling_prog = "apriori_model_timeseries"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the command line
  call getarg(1,mascon_file)
  if(mascon_file(1:1) == " ")then
    print*,"apriori_model_timeseries mascons_stage5_V006_200km msc_model.all lat lon bit-mapped_model"
    print*," where bit-mapped models are: "
    print*,"       1: low-frequency component"
    print*,"       2: annual amplitude"
    print*,"       4: WRMS of total time series"
    print*,"       8: WRMS of land time series, 0.1m ocean"
    print*,"      16: WRMS of high-frequecy component"
    print*,"      32: annual amplitude if > 0.2 m"
    print*," "
    stop
  endif
  
  call getarg(2,msc_model_file)  ! mascon file of modelled values
  open(lumsc_model,file=msc_model_file,status='old',iostat=ioerr)
  if(ioerr /= 0)then
    call status_update('FATAL','UTIL',calling_prog,msc_model_file &
                     ,"Error opening mascon apriori model file. Does it exist?",ioerr)
  else
    call status_update('STATUS','UTIL',calling_prog,msc_model_file,"Opened mascon apriori model file. ",ioerr)
  endif

  call getarg(3,arg)
  read(arg,*)user_lat
  call getarg(4,arg)
  read(arg,*)user_lon
  write(message,'(a,2f6.2,a)')"Will compute time series of apriori values at location (",user_lat,user_lon,")"
  call status_update('STATUS','UTIL',calling_prog,msc_model_file,message,0)
    ! bit-mapped model required
  call getarg(5,arg)
  read(arg,*)model

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the mascon file that relates to the mascon EWH model
  call read_msc_hdr(lumsc,mascon_file,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                                ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
  call allocate_mascon_arrays
! PT180220: using trimlen here doesn't work - the subroutine is expecting a C*150, so pass the whole thing through
!  call read_mascon_file(lumsc_in,mascon_file(1:trimlen(mascon_file)))
  call read_mascon_file(lumsc,mascon_file)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! open and read the apriori model file
! read the file. First line is a comment, giving column headers. Second line contains number of mascons.
  open(lumsc_model,file=msc_model_file,status='old',iostat=ioerr)
  if(ioerr /= 0)call status_update('FATAL','UTIL',calling_prog,msc_model_file,"Error opening file",0)

  read(lumsc_model,'(a)')line
  read(lumsc_model,*)nmsc,nmonths,ntot

! the "ntot" value in the file includes lat/lon/density/cosine_ampl/sine_ampl
  nrms = ntot - nmonths - 5
  allocate(msc_crds(nmsc,3))          ! lat, lon, density
  allocate(msc_annual(nmsc,2))        ! cosine/sine amplitudes
  allocate(msc_WRMS(nmsc,nrms))
  allocate(msc_monthly(nmsc,nmonths))
  allocate(epochs(nmonths))

! read the epochs contained in the file. They are in decimal years
  read(lumsc_model,'(a)')line
  read(lumsc_model,*)(epochs(imonth),imonth=1,nmonths)
print*,'Epochs in apriori model file: ',epochs

! ok, now loop through mascons and read in all the information
  read(lumsc_model,'(a)')line
  print*,'reading values for each mascon'
  do imsc=1,nmsc
    read(lumsc_model,*)msc_crds(imsc,:),msc_annual(imsc,:),(msc_wrms(imsc,iRMS),iRMS=1,nrms) &
                ,(msc_monthly(imsc,imonth),imonth=1,nmonths)

  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! determine in which mascon the requested lat/lon resides
  ! set the ternary latitudinal spacing
  lat_spacing = 10.d0/60.d0     ! default is 10' ternary mascons in latitude

  ! establish the array that contains how many ternary latitude bands and how many ternary mascons per band
  call tern_lat_bands_ell(lat_spacing)
  call calc_which_ternary(.false.,user_lat,user_lon,lat_spacing,tern_number)
  imsc = mcon_tern_ptr(tern_number,1)
  print*,'mascon at ',user_lat,user_lon,' is in primary mascon ',imsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now, starting at 2003 (dt = 0.d0), calculate the a priori model each month
  dt = 1.d0 !8.d0/12.d0
  do while (2002.d0 + dt < 2023.d0)
    soln_epoch = 2002.d0+dt
    ! find the two monthly values that straddle the required epoch
    imonth = 1
    do while ( imonth <= nmonths)
      if(epochs(imonth) < 2002.d0 + dt)then
        imonth = imonth+1
      else
        exit ! this will stop it looping
      endif
    enddo

    ! low-frequency component of model.
    if(bitmap(model,1) )then  
!print*,'in the low-freq part of the code',model,bitmap(model,1)
      ! interpolate between i2 and i1
      if(imonth == 1)then
        write(message,'(a,2(f15.4,a))')"Requested epoch (",soln_epoch,") before first epoch in file (",epochs(1) &
                                        ,"). Using model of the first epoch."
        !call status_update('STATUS',calling_prog,'read_msc_model',msc_model_file,message,0)
        lowfrq = msc_monthly(imsc,1)
      else if (imonth > nmonths)then
        write(message,'(a,2(f15.4,a))')"Requested epoch (",soln_epoch,") after last epoch in file (",epochs(nmonths) &
                                        ,"). Using model of the last epoch."
        !call status_update('STATUS',calling_prog,'read_msc_model',msc_model_file,message,0)
        lowfrq = msc_monthly(imsc,nmonths)
      else
        i2 = imonth
        i1 = imonth - 1
        dt2 = epochs(i2)-epochs(i1)
        write(message,'(3(a,f15.4))')"Requested epoch (",soln_epoch,") between epochs",epochs(i1)," and",epochs(i2)
        !call status_update('STATUS',calling_prog,'read_msc_model',msc_model_file,message,0)
        lowfrq = msc_monthly(imsc,i1) + (soln_epoch-epochs(i1)) * (msc_monthly(imsc,i2)-msc_monthly(imsc,i1))/dt2
      endif
      msc_apr_value =  lowfrq

    endif

    ! annual variation. Used differently for a priori value vs regularisation constraint. 
    if(bitmap(model,2) )then  
        msc_apr_value = msc_apr_value + msc_annual(imsc,1)*dcos(2.d0*pi*soln_epoch) &
                                            + msc_annual(imsc,2)*dsin(2.d0*pi*soln_epoch)
    endif

    ! output the computed value
    write(*,'(f10.4,f10.4,2f10.4,a)')soln_epoch,msc_apr_value,user_lat,user_lon," apriori model value (mm EWH)"
    
    ! increment dt by one month
    dt = dt + 1.d0/12.d0
      
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



 
  
  
  end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
