  program norm2netcdf
  
! an attempt to convert a binary .norm file into a netcdf. Hopefully, we can write it out as a netCDF4 compressed !
!
! P. Tregoning
! 23 August 2022

  use netcdf
  use soln_mod   ! PT221122: declares the arrays to be read from the .norm file
  
  implicit none

! command line arguments
  character*100 :: norm_in,netcdf_out

! norm file header variables
  integer*4     :: nparams               ! total number of paramters
  integer*4     :: nICs                  ! number of IC parameters in a binary file
  integer*4     :: nmascons              ! number of mascon parameters in a binary file
  integer*4     :: nmsc_tidal            ! number of tidal mascon parameters in a binary file
  real(kind=8 ) :: version               ! version number of binary normal equation file
  integer*4     :: rec_len               ! required record length for the binary norm file
  
  ! first line
  integer*4,allocatable     :: year(:), month(:), day(:), hr(:), minute(:), sec(:)

! arrays for data from the input file
  real(kind=8), allocatable :: tmpvals(:)
  character*1               :: sats(2)

! unit numbers
  integer*4, parameter :: luin = 10, luout = 11

! variables to read daily solution dates from a snorm file
  integer*4,parameter    :: max_epochs=400

! variables to reorder the normal equations etc to put mascon parameters at the top
  real(kind=8),allocatable  :: apriori_tmp(:),AtWA_tmp(:,:),AtWb_tmp(:)
  integer*4,    allocatable :: param_type_tmp(:)

! other stuff
  character*200 :: message
  integer*4     :: i,j,k,irec,iMsc,nsolns,isoln,ndates


! decode the runstring
  call getarg(1,norm_in)
  if(norm_in(1:1) == "")then
    stop 'Runstring: norm2cdf input.norm output_norm.nc'
  endif
  call getarg(2,netcdf_out)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

! open the input binary norm file
  open(luin,file=norm_in,access='direct',recl=24)

! read the first line of the input norm file
  read(luin,rec=1)version,nparams,nmascons,nICs,nmsc_tidal

  if(version == 1.0)then
    ! it is a .norm file from gracefit/sim. There is only a single date inside it. Allocate arrays
    nsolns = 1
  else
    nsolns = max_epochs
  endif
  allocate(year(nsolns))
  allocate(month(nsolns))
  allocate(day(nsolns))
  allocate(hr(nsolns))
  allocate(minute(nsolns))
  allocate(sec(nsolns))
  allocate(duration(nsolns))
  allocate(file_names(nsolns))

!!!!!! allocate arrays to fill from the input binary file
  allocate(tmpvals(nparams))
  allocate(apriori(nparams))
  allocate(param_type(nparams))
  allocate(AtWA(nparams,nparams))
  allocate(AtWb(nparams))
  AtWA = 0.d0
  AtWb = 0.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! close and reopen with the correct record length
  close(luin)
  rec_len = nparams*8
  open(luin,file=norm_in,access='direct',recl=rec_len)


! read the header information again, this time all of it
  nsolns = 1    ! default value for a binary .norm file. If it is a snorm then the solution dates are read from record 4.
  read(luin,rec=1)version,nparams,nmascons,nICs,nmsc_tidal,year(nsolns),month(nsolns),day(nsolns),hr(nsolns) &
                 ,minute(nsolns),sec(nsolns),duration(nsolns)
  write(message,'(a,f3.1,4(a,i6),a,6i5,a,f7.3 )') &
            "File version: ",version," Num params: ",nparams," Num mascons: ",nmascons," Num ICs: "&
            ,nICs, " Num tidal mascons: ",nmsc_tidal," Date: ",year(nsolns),month(nsolns),day(nsolns) &
          ,hr(nsolns),minute(nsolns),sec(nsolns)," Duration: ",duration(nsolns)
  call status_update('STATUS','UTIL','norm2cdf',' ',message,0)

  ! line 2 of version 1.0 and version 2.0 files are the a priori values of the parameters
  read(luin,rec=2)(apriori(j),j=1,nparams)

  ! line 3 of version 1.0 files contains (coded) descriptors for the parameters
  irec=3  
  if(version == 1.0)then
    if(nmascons == 0) then                ! the parameters will be only ICs
       read(luin,rec=irec)sats(1),(param_type(j),j=1,nICs/2),sats(2),(param_type(j),j=nICs/2+1,nICs)
    else if (nICs == 0) then      ! the parameters will be only mascons and possibly tidal amplitudes (ie no satellite indicator first).
       read(luin,rec=irec)(param_type(j),j=1,nmascons+nmsc_tidal)
    else
       read(luin,rec=irec)sats(1),(param_type(j),j=1,nICs/2),sats(2),(param_type(j),j=nICs/2+1,nICs) &
              ,(param_type(nICs+j),j=1,nmascons+nmsc_tidal)
    endif
  else 
    read(luin,rec=irec)param_type,sats(:)
  endif

!! Line 4
  if(version == 2.0)then 
    ! it is a snorm version 2 file. It can have multiple input daily files that made the solution within it
    call status_update('STATUS','UTIL','norm2netcdf',' ','Reading dates of stacked normal equations',0)
    irec = 4
    read(luin,rec=irec)ndates,(year(isoln),month(isoln),day(isoln),hr(isoln),minute(isoln),sec(isoln) &
                          ,duration(isoln), file_names(isoln),isoln=1,ndates) 
    nsolns = ndates
  else if (version == 1.0)then
    file_names(1) = trim(norm_in)   
  endif

! record 5 is the tidal mascons to be estimated. Ignore this - we no longer estimate tidal mascons.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

! DEBUG
do i=1,ndates
  print*,i,year(i),month(i),day(i),hr(i),minute(i),sec(i),duration(i)," ",trim(file_names(i))
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! Read in the AtWA matrix
  irec = 6
  call status_update('STATUS','UTIL','norm2netcdf',norm_in,'Reading input AtWA',0)
  do k=1,nparams
    read(luin,rec=irec)AtWA(k,:)
    irec = irec + 1
  enddo

! Read in the AtWb matrix
  call status_update('STATUS','UTIL','norm2netcdf',norm_in,'Reading input AtWb',0)
  read(luin,rec=irec)AtWb(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! reconfigure the order of parameters for .norm files to match the .snorm files. That is, we
! want it as mascons first, then ICs.
  iMsc = nICs+1
  if(version == 1.0)then
    call status_update('STATUS','UTIL','norm2netcdf',norm_in,'Reordering the apriori parameter information',0)
    allocate(apriori_tmp(nparams))
    apriori_tmp(1:nmascons) = apriori(iMsc:iMsc+nmascons-1)
    apriori_tmp(nmascons+1:nparams) = apriori(1:nICs)
    apriori = apriori_tmp
    deallocate(apriori_tmp)
    
    call status_update('STATUS','UTIL','norm2netcdf',norm_in,'Reordering the parameter type information',0)
    allocate(param_type_tmp(nparams))
    param_type_tmp(1:nmascons) = param_type(iMsc:iMsc+nmascons-1)
    param_type_tmp(nmascons+1:nparams) = param_type(1:nICs)
    param_type = param_type_tmp
    deallocate(param_type_tmp)
    
    call status_update('STATUS','UTIL','norm2netcdf',norm_in,'Reordering AtWb',0)
    allocate(AtWb_tmp(nparams))
    AtWb_tmp(1:nmascons) = AtWb(iMsc:iMsc+nmascons-1)
    AtWb_tmp(nmascons+1:nparams) = AtWb(1:nICs)
    AtWb = AtWb_tmp
    deallocate(AtWb_tmp)

    call status_update('STATUS','UTIL','norm2netcdf',norm_in,'Reordering AtWA',0)
    allocate(AtWA_tmp(nparams,nparams))
    !AtWA_tmp = AtWA
    ! the diagonal block of the mascon parameters
    AtWA_tmp(1:nmascons,1:nmascons) = AtWA(iMsc:iMsc+nmascons-1,iMsc:iMsc+nmascons-1)
    ! the diagonal block of the ICs
    AtWA_tmp(nmascons+1:nparams,nmascons+1:nparams) = AtWA(1:nICs,1:nICs)
    ! the rows of off-diagonal information. Store it in the upper-triangular part of the matrix
    do j=1,nmascons
        AtWA_tmp(j,nmascons+1:nmascons+nICs) = AtWA(1:nICs,nICs+j)
    enddo
    AtWA = AtWA_tmp
    deallocate(AtWA_tmp)
  endif    
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! debug:
!do i=1,10
!  print*,i,AtWA(i,1:10)
!enddo  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! detail of the ymdhms/duration for the output netcdf file
  allocate(epoch(nsolns,6))
  do i=1,nsolns
    epoch(i,1)=year(i)
    epoch(i,2)=month(i)
    epoch(i,3)=day(i)
    epoch(i,4)=hr(i)
    epoch(i,5)=minute(i)
    epoch(i,6)=sec(i)
    !file_names(i) = trim(norm_in)
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now write out a netcdf file
  version = 3.0
  call write_norm_netcdf('norm2netcdf',netcdf_out,version,nparams,nmascons,nICs,nmsc_tidal,nsolns, AtwA,AtWb, param_type &
                                   ,apriori,file_names,epoch,duration,sats)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  end


