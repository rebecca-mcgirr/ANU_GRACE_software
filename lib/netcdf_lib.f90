!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutines related to netcdf files
!
! read_norm_netcdf    :  reads the information of a .norm file but from its netcdf format
! write_norm_netcdf   :  writes the info of a binary .norm file into a netcdf file
! check_netcdf_status :  error checking for netcdf operations within fortran
! reorder_normeq      :  flip the ordering of normal equations etc so mascons are at the top, then ICs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_norm_netcdf(calling_prog,what_to_read,infile,version,nparams,nmascons,nICs,nmsc_tidal,nsolns,satellites)

! read normal equations and RHS from a netcdf file. There are multiple versions of the file:
!
! Version 1.0 : the entire AtWA is written out, order in the matrix is mascons then ICs
! Version 2.0 : the "mascons", "ICs" and "msc/ICs" blocks are written out and need to be formed into the AtWA matrix
!
! in all cases, the "apriori", "AtWb" and "param_type" vectors are mascons then ICs

  use netcdf
!  use norm_netcdf_mod
  use soln_mod  ! PT221122: variables that were in norm_netcdf_mod are now declared in soln_mod
  
  implicit none

  character(*)  , intent(in) :: calling_prog    ! name of calling program
  character(*) , intent(in)  :: what_to_read    ! "dimensions","file_names" or "all"
  integer*4     :: nmascons,nICs,nmsc_tidal     ! numbers of each of these parameters in file "infile"
  character(*)  :: infile
  character*150 :: message
  integer*4     :: nstatus,i,j,k
  integer*4     :: nparams,nsolns      ! number of parameters and daily solutions in stacked netcdf file
  character(*)  :: satellites(2)       ! names of the two space gravity satellites

  character*50  :: vname
  integer*4     :: dim1
  real(kind=8)  :: version
  integer*4     :: nver

! PT221125: arrays to read different blocks of the AtWA matrix (for V2.0)
  real(kind=8),allocatable :: AtWA_mascons(:,:)
  real(kind=8),allocatable :: AtWA_ICs(:,:)
  real(kind=8),allocatable :: AtWA_mscIC(:,:)
  
! netcdf variables to read out the file names character array into the netcdf
!  character*150, allocatable :: file_names(:)
  integer*4 :: namelen_dimid,n_strings_dimid,file_name_varid,n_strings,namelen,nepoch_var,satslen

! stuff
  integer(kind=4) :: ncid, xtype, ndims, varid,dimids(2)
  integer*4,parameter :: maxvar=11                          ! maximum number of variables to loop through to read the file
  integer*4       :: ivar                                   ! loop counter to loop through to read the file
  logical :: debug


  write(message,'(a,a,a)')'Reading "',trim(what_to_read),'" from input netcdf normal equations file'
  call status_update('STATUS',calling_prog,'read_norm_netcdf',trim(infile),message,0)

! open the file
  nstatus = nf90_open(infile, nf90_nowrite, ncid)
  call check_netcdf_status (nstatus,"open file")

  
debug = .false.
  
  if(what_to_read(1:3) == "dim" .or. what_to_read(1:3) == "all")then
  
    ! get the dimensions of the variables
    nstatus = nf90_inquire_dimension(ncid,1,vname,nparams)
    call check_netcdf_status (nstatus,"nparams")

    nstatus = nf90_inquire_dimension(ncid,2,vname,nmascons)
    call check_netcdf_status (nstatus,"nmascons")

    nstatus = nf90_inquire_dimension(ncid,3,vname,nmsc_tidal)
    call check_netcdf_status (nstatus,"nmsc_tidal")

    nstatus = nf90_inquire_dimension(ncid,4,vname,nICs)
    call check_netcdf_status (nstatus,"nICs")
  
    nstatus = nf90_inquire_dimension(ncid,5,vname,nsolns)
    call check_netcdf_status (nstatus,"nsolns")

    nstatus = nf90_inquire_dimension(ncid,6,vname,nepoch_var)
    call check_netcdf_status (nstatus,"nepoch_var")

    nstatus = nf90_inquire_dimension(ncid,8,vname,namelen)
    call check_netcdf_status (nstatus,"namelen")

    nstatus = nf90_inquire_dimension(ncid,9,vname,n_strings)
    call check_netcdf_status (nstatus,"n_strings")

    nstatus = nf90_inquire_dimension(ncid,10,vname,satslen)
    call check_netcdf_status (nstatus,"satslen")
  endif
  
  if(what_to_read(1:3) == "all")then
   
    ! allocate arrays
    allocate(AtWA(nparams,nparams))
    allocate(AtWb(nparams))
    allocate(epoch(nsolns,6))
    allocate(apriori(nparams))
    allocate(param_type(nparams))
    allocate(duration(nsolns))

    ! now read in the data
    ! AtWA
    ! PT221125: Version 1.0 is the whole AtWA, mascons then ICs
    do ivar=1,maxvar
      nstatus = nf90_inquire_variable(ncid,ivar,vname,xtype,ndims,dimids)
if(debug)print*,ivar," vname = ",trim(vname)
      ! the full AtWA matrix (from version 1.0 netcdf file formats)
      if(trim(vname) == "AtWA")then
        nstatus = nf90_inq_varid(ncid,vname,varid)
        nstatus = nf90_get_var(ncid,varid,AtWA)
        call check_netcdf_status (nstatus,"Reading AtWA")
      endif
      
      if(trim(vname) == "AtWA_mascons")then
      ! mascons block from V2.0 file
        nstatus = nf90_inq_varid(ncid,vname,varid)
        nstatus = nf90_get_var(ncid,varid,AtWA(1:nmascons,1:nmascons))
        call check_netcdf_status (nstatus,"Reading AtWA_mascons")
      endif
      
      if(trim(vname) == "AtWA_ICs")then
      ! IC block from V2.0 file
        nstatus = nf90_inq_varid(ncid,vname,varid)
        nstatus = nf90_get_var(ncid,varid,AtWA(nmascons+1:nmascons+nICs,nmascons+1:nmascons+nICs))
        call check_netcdf_status (nstatus,"Reading AtWA_IC")
      endif
      if(trim(vname) == "AtWA_mscIC")then
      ! msc/IC block from V2.0 file
        nstatus = nf90_inq_varid(ncid,vname,varid)
        nstatus = nf90_get_var(ncid,varid,AtWA(1:nmascons,nmascons+1:nmascons+nICs))
        call check_netcdf_status (nstatus,"Reading AtWA_mscIC")

if(debug)print*,'last row of AtWA_mscIC',nmascons,nICs,AtWA(nmascons,nmascons+1:nmascons+nICs)
      endif

      if(trim(vname) == "AtWb")then
        nstatus = nf90_inq_varid(ncid,vname,varid)
        nstatus = nf90_get_var(ncid,varid,AtWb)
        call check_netcdf_status (nstatus,"Reading AtWb")
      endif

      if(trim(vname) == "apriori")then
        nstatus = nf90_inq_varid(ncid,vname,varid)
        nstatus = nf90_get_var(ncid,varid,apriori)
        call check_netcdf_status (nstatus,"Reading apriori")
        if(debug)print*,'ICs from apriori',apriori(nmascons+1:nmascons+nICs)
      endif

      if(trim(vname) == "param_type")then
        nstatus = nf90_inq_varid(ncid,vname,varid)
        nstatus = nf90_get_var(ncid,varid,param_type)
        call check_netcdf_status (nstatus,"Reading parameter types")
      endif

      if(trim(vname) == "epoch")then
        nstatus = nf90_inq_varid(ncid,vname,varid)
        nstatus = nf90_get_var(ncid,varid,epoch)
        call check_netcdf_status (nstatus,"Reading epochs")
      endif

      if(trim(vname) == "duration")then
        nstatus = nf90_inq_varid(ncid,vname,varid)
        nstatus = nf90_get_var(ncid,varid,duration)
        call check_netcdf_status (nstatus,"Reading daily solution durations")
      endif

      if(trim(vname) == "version")then
        nstatus = nf90_inq_varid(ncid,vname,varid)
        nstatus = nf90_get_var(ncid,varid,version)
        call check_netcdf_status (nstatus,"Reading version of the file format")
      endif

      if(trim(vname) == "sats")then
        nstatus = nf90_inq_varid(ncid,vname,varid)
        nstatus = nf90_get_var(ncid,varid,satellites)
        call check_netcdf_status (nstatus,"Reading satellite names")
      endif
    enddo

  endif
  
  if(what_to_read(1:3) == "fil" .or. what_to_read(1:3) == "all")then
    allocate(file_names(nsolns))
    do ivar=1,maxvar
      nstatus = nf90_inquire_variable(ncid,ivar,vname,xtype,ndims,dimids)

      if(trim(vname) == "file_names")then 
        ! daily file names
        call check_netcdf_status (nstatus,"inquire variable file names")
        nstatus = nf90_inq_varid(ncid,vname,varid)
        call check_netcdf_status (nstatus,"inquire varid file names")
        nstatus = nf90_get_var(ncid,varid,file_names)
        call check_netcdf_status (nstatus,"Reading file names")
        if(debug)print*,"file names = : ",file_names
      endif
    enddo
    
  endif

! close the file
  nstatus = nf90_close(ncid)
  call check_netcdf_status(nstatus, 'close netcdf')
  

  return
  end subroutine read_norm_netcdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_norm_netcdf(calling_prog,netcdf_out,version,nparams,nmascons,nICs,nmsc_tidal,nsolns, AtwA,AtWb &
                                   , param_type,apriori,file_names,epoch,duration,sats)
  
! write out the information of a binary .norm file inot netcdf format
!
! P. Tregoning
! 23 August 2022

  use netcdf

  implicit none

! passed arguments
  character(*) , intent(in)  :: calling_prog            ! name of calling program
  character(*) , intent(in)  :: netcdf_out              ! name of output netcdf file
  real(kind=8) , intent(in)  :: version                 ! version of the output file. netCDF4 versions start at 3.0
  integer*4    , intent(in)  :: nparams                 ! number of parameters
  integer*4    , intent(in)  :: nICs                    ! number of IC parameters
  integer*4    , intent(in)  :: nmascons                ! number of mascon parameters
  integer*4    , intent(in)  :: nmsc_tidal              ! number of tidal mascon parameters (hardwired to zero - we don't estimate them anymore)
  integer*4    , intent(in)  :: nsolns                  ! number of daily solutions contained within input file
  real(kind=8) , intent(in)  :: AtWA(nparams,nparams)   ! normal equations
  real(kind=8) , intent(in)  :: AtWb(nparams)           ! RHS
  integer*4    , intent(in)  :: param_type(nparams)     ! identifier of type of parameter
  real(kind=8) , intent(in)  :: apriori(nparams)        ! a priori values of the parameters
  character(*) , intent(in)  :: file_names(nsolns)      ! file names of solutions contained within input file
  integer*4    , intent(in)  :: epoch(nsolns,6 )        ! epochs of the solutions (ymdhms)
  real(kind=8) , intent(in)  :: duration(nsolns)        ! decimal hours of the duration of each solution
  character*1  , intent(in)  :: sats(2)                 ! 1-char names of satellites


! netcdf variables
  integer*4 :: ncid            ! id for netcdf file (aka "unit number")
  integer*4 :: nstatus         ! integer i/o value from a nnetcdf action
  integer*4 :: dimid_nparams,dimid_nICs,dimid_nmascons,dimid_nmsc_tidal,dimid_version   ! dimension id for number of parameters and number of IC parameters

  integer*4 :: dimid_nsolns    ! dimension id for number of daily solutions in file
  integer*4 :: dimid_nepochvar ! dimension id for number of ymdhms (=6)
  integer*4, allocatable :: param_nums(:)
  integer*4 :: varid_param_nums,varid_nsolns,varid_AtWA,varid_AtWb,varid_param_type,varid_apriori
  integer*4 :: varid_epoch,varid_duration,varid_version

! netcdf variables to write out the file names character array into the netcdf
  integer*4 :: namelen_dimid,n_strings_dimid,file_name_varid
  
! variables to write out the satellite names
  integer*4 :: satslen_dimid,n_strings2_dimid,sats_varid

! other variables
  integer*4 :: i
  logical   :: debug
  
  debug = .false.

! open the netcdf file
  nstatus = nf90_create(trim(netcdf_out),NF90_NETCDF4, ncid)
  call check_netcdf_status(nstatus,"open netcdf file")
  nstatus = nf90_enddef(ncid)

! dimensions of the array variables
  nstatus = nf90_def_dim(ncid, 'nparams', nparams, dimid_nparams)
  nstatus = nf90_def_dim(ncid, 'nmascons', nmascons, dimid_nmascons)
  nstatus = nf90_def_dim(ncid, 'nmsc_tidal', nmsc_tidal, dimid_nmsc_tidal)
  nstatus = nf90_def_dim(ncid, 'nICs', nICs, dimid_nICs)
  nstatus = nf90_def_dim(ncid, 'nsolns', nsolns, dimid_nsolns)
  nstatus = nf90_def_dim(ncid, 'nepoch_var', 6, dimid_nepochvar)
  nstatus = nf90_def_dim(ncid, 'version',1,dimid_version)
  
! write the header info of the variables
!  nstatus = nf90_def_var(ncid, 'AtWA', NF90_FLOAT, [dimid_nparams,dimid_nparams], varid_AtWA)
  nstatus = nf90_def_var(ncid, 'AtWA', NF90_DOUBLE, [dimid_nparams,dimid_nparams], varid_AtWA)
  call check_netcdf_status(nstatus,"define AtWA")
!  nstatus = nf90_def_var(ncid, 'AtWb', NF90_FLOAT, [dimid_nparams], varid_AtWb)
  nstatus = nf90_def_var(ncid, 'AtWb', NF90_DOUBLE, [dimid_nparams], varid_AtWb)
  call check_netcdf_status(nstatus,"define AtWb")
  nstatus = nf90_def_var(ncid, 'apriori', NF90_DOUBLE, [dimid_nparams], varid_apriori)
  call check_netcdf_status(nstatus,"define apriori")
  nstatus = nf90_def_var(ncid, 'param_type', NF90_INT2, [dimid_nparams], varid_param_type)
  call check_netcdf_status(nstatus,"define param_type")
  nstatus = nf90_def_var(ncid, 'epoch', NF90_INT2, [dimid_nsolns,dimid_nepochvar], varid_epoch)
  call check_netcdf_status(nstatus,"define epoch")
  nstatus = nf90_def_var(ncid, 'duration', NF90_FLOAT, [dimid_nsolns], varid_duration)
  call check_netcdf_status(nstatus,"define duration")
  nstatus = nf90_def_var(ncid, 'version', NF90_FLOAT, [dimid_version], varid_version)
  call check_netcdf_status(nstatus,"define version")
if(debug)print*,'written header info'


! PT221123: only chunk the data if mascons have been estimated
  if(nparams > 1000)then
    ! prepare "chunking" (a process of compressing the data output in the netcdf file)
    nstatus = nf90_def_var_chunking(ncid, varid_AtWA, NF90_CHUNKED, [10, 101])
    call check_netcdf_status(nstatus,"chunked")

    ! write out AtWA
    nstatus = nf90_def_var_deflate(ncid, varid_AtWA,          &
                                  shuffle = 1,                &
                                  deflate = 1,                &
                                  deflate_level = 5  )
    call check_netcdf_status(nstatus,"deflated AtWA")
    if(debug)print*,'chunked the AtWA'
  endif

! write the attributes
  nstatus = nf90_put_att(ncid, NF90_GLOBAL, 'note', 'netcdf version of daily .norm file')
  call check_netcdf_status(nstatus,"note attribute")
  nstatus = nf90_put_att(ncid, varid_AtWA, 'AtWA', 'unregularised normal equations')
  call check_netcdf_status(nstatus,"AtWA attribute")
  nstatus = nf90_put_att(ncid, varid_AtWb, 'AtWb', 'RHS')
  call check_netcdf_status(nstatus,"AtWb attribute")
  nstatus = nf90_put_att(ncid, varid_apriori, 'apriori', 'a priori values of the parameters')
  call check_netcdf_status(nstatus,"apriori attribute")
  nstatus = nf90_put_att(ncid, varid_param_type, 'param_type', 'integer code for type of parameter')
  call check_netcdf_status(nstatus,"param_type attribute")
  nstatus = nf90_put_att(ncid, varid_epoch, 'epoch', 'ymdhms of each solution in file')
  call check_netcdf_status(nstatus,"epoch attribute")
  nstatus = nf90_put_att(ncid, varid_duration, 'duration', 'decimal hours of duration of each solution in file')
  call check_netcdf_status(nstatus,"duration attribute")
  nstatus = nf90_put_att(ncid, varid_version, 'version', 'File format version')
  call check_netcdf_status(nstatus,"version attribute")
if(debug)print*,'written attributes'

! write it out
  call status_update('STATUS',calling_prog,'write_norm_netcdf',netcdf_out,'writing netcdf file',0)
  nstatus = nf90_put_var(ncid, varid_version, version)
  call check_netcdf_status(nstatus,"writing file format version")
  if(debug)print*,'written version'

  nstatus = nf90_put_var(ncid, varid_AtWA, AtWA)
  call check_netcdf_status(nstatus,"writing AtWA")
  if(debug)print*,'written AtWA'

  nstatus = nf90_put_var(ncid, varid_AtWb, AtWb)
  call check_netcdf_status(nstatus,"writing AtWb")
  if(debug)print*,'written AtWb'

  nstatus = nf90_put_var(ncid, varid_apriori, apriori)
  call check_netcdf_status(nstatus,"writing apriori")
  if(debug)print*,'written apriori'

  nstatus = nf90_put_var(ncid, varid_param_type, param_type)
  call check_netcdf_status(nstatus,"writing param_type")
  if(debug)print*,'written param_type'

  nstatus = nf90_put_var(ncid, varid_epoch, epoch)
  if(debug)print*,'epoch = ',epoch
  call check_netcdf_status(nstatus,"writing epoch info")
  if(debug)print*,'written epoch'

  nstatus = nf90_put_var(ncid, varid_duration, duration)
  call check_netcdf_status(nstatus,"writing soln durations")
  if(debug)print*,'written durations'




! try to write the character array "file_names" to the netcdf file. This code is based on information from
! http://computer-programming-forum.com/49-fortran/af9bfa2653fff83f.htm
  ! define the maximum length of the string
  nstatus = nf90_def_dim (ncid, 'namelen',len(file_names(1)),namelen_dimid)
  ! "same" as for regular dimension of an array
  nstatus = nf90_def_dim (ncid, 'n_strings',size(file_names),n_strings_dimid)
  ! now define char array variable with string length as the first dimension. That is, write it out as though it is a 2-dimensional character array
  nstatus = nf90_def_var(ncid,'file_names',NF90_CHAR &
                        ,dimIDs=(/ namelen_dimid,n_strings_dimid /) &
                        , VarID = file_name_varid )
  ! now write the character array
  nstatus = nf90_put_var(ncid, file_name_varid, file_names)
  call check_netcdf_status(nstatus,"writing soln file names")
  if(debug)print*,'written file_names'
  
  
! write out the satellite names as 1-char character codes (A, B, C, D)
  ! define the maximum length of the string (here it is char*1)
  nstatus = nf90_def_dim (ncid, 'satslen',len(sats(1)),satslen_dimid)
  call check_netcdf_status(nstatus,"writing satellite names 1")
  ! "same" as for regular dimensions
  nstatus = nf90_def_dim (ncid, 'n_strings2',size(sats),n_strings2_dimid)
  call check_netcdf_status(nstatus,"writing satellite names 2")
  ! now define char array variable with string length as the first dimension
  nstatus = nf90_def_var(ncid,'sats',NF90_CHAR &
                        ,dimIDs=(/ satslen_dimid,n_strings2_dimid /) &
                        , VarID = sats_varid )
  call check_netcdf_status(nstatus,"writing satellite names 3")
  ! now write the character array
  nstatus = nf90_put_var(ncid, sats_varid, sats)
  call check_netcdf_status(nstatus,"writing satellite names 4")
  if(debug)print*,'written satellite names'

! close the file
  nstatus = nf90_close(ncid)
  call check_netcdf_status(nstatus, 'close netcdf')

  return
  end subroutine write_norm_netcdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_norm_netcdf_v2(calling_prog,netcdf_out,version,nparams,nmascons,nICs,nmsc_tidal,nsolns &
                                   , AtwA_mascons, AtWA_ICs, AtWA_mscIC, AtWb &
                                   , param_type,apriori,file_names,epoch,duration,sats,flip,nrows,ncols)
  
! write out the information of a binary .norm file inot netcdf format
!
! P. Tregoning
! 23 August 2022
!
! PT221123: changed routine to write out three "blocks" of the AtWA normal equations, as separate bits
!           1. the mascons square
!           2. the ICs square
!           3. the column of off-diagonal mascons/IC entries

  use netcdf

  implicit none

! passed arguments
  character(*) , intent(in)  :: calling_prog            ! name of calling program
  character(*) , intent(in)  :: netcdf_out              ! name of output netcdf file
  real(kind=8) , intent(in)  :: version                 ! version of the output file. netCDF4 versions start at 3.0
  integer*4    , intent(in)  :: nparams                 ! number of parameters
  integer*4    , intent(in)  :: nICs                    ! number of IC parameters
  integer*4    , intent(in)  :: nmascons                ! number of mascon parameters
  integer*4    , intent(in)  :: nmsc_tidal              ! number of tidal mascon parameters (hardwired to zero - we don't estimate them anymore)
  integer*4    , intent(in)  :: nsolns                  ! number of daily solutions contained within input file
  real(kind=8) , intent(in)  :: AtWA_mascons(nmascons,nmascons) ! mascons normal equations
  real(kind=8) , intent(in)  :: AtWA_ICs(nICs,nICs)             ! ICs normal equations
  integer*4    , intent(in)  :: nrows,ncols                     ! dimensioning of the input off-diagonal block between ICs and mascons
  real(kind=8) , intent(in)  :: AtWA_mscIC(nrows,ncols)         ! msc/IC normal equations (but it could be IC/msc normal equations!)
  real(kind=8) , intent(in)  :: AtWb(nparams)           ! RHS
  integer*4    , intent(in)  :: param_type(nparams)     ! identifier of type of parameter
  real(kind=8) , intent(in)  :: apriori(nparams)        ! a priori values of the parameters
  character(*) , intent(in)  :: file_names(nsolns)      ! file names of solutions contained within input file
  integer*4    , intent(in)  :: epoch(nsolns,6 )        ! epochs of the solutions (ymdhms)
  real(kind=8) , intent(in)  :: duration(nsolns)        ! decimal hours of the duration of each solution
  character*1  , intent(in)  :: sats(2)                 ! 1-char names of satellites
  logical      , intent(in)  :: flip                    ! does the msc/IC block need to be transposed ?


! netcdf variables
  integer*4 :: ncid            ! id for netcdf file (aka "unit number")
  integer*4 :: nstatus         ! integer i/o value from a nnetcdf action
  integer*4 :: dimid_nparams,dimid_nICs,dimid_nmascons,dimid_nmsc_tidal,dimid_version   ! dimension id for number of parameters and number of IC parameters

  integer*4 :: dimid_nsolns    ! dimension id for number of daily solutions in file
  integer*4 :: dimid_nepochvar ! dimension id for number of ymdhms (=6)
  integer*4, allocatable :: param_nums(:)
  integer*4 :: varid_param_nums,varid_nsolns,varid_AtWb,varid_param_type,varid_apriori
  integer*4 :: varid_epoch,varid_duration,varid_version
  integer*4 :: varid_AtWA_mascons,varid_AtWA_ICs,varid_AtWA_mscIC

! netcdf variables to write out the file names character array into the netcdf
  integer*4 :: namelen_dimid,n_strings_dimid,file_name_varid
  
! variables to write out the satellite names
  integer*4 :: satslen_dimid,n_strings2_dimid,sats_varid

! other variables
  integer*4 :: i
  logical   :: debug
  real(kind=8),allocatable :: tmp_mscIC(:,:)
  
  debug = .false.

  if(debug)print*,'inside write_norm_netcdf_v2'
  
! PT221123: check whether we need to transpose the block of msc/IC information
  if (flip)then
    allocate(tmp_mscIC(ncols,nrows))
    tmp_mscIC = transpose(AtWA_mscIC)
  endif
  
  


! open the netcdf file
  if(debug)print*,'open the netcdf file to be created'
  nstatus = nf90_create(trim(netcdf_out),NF90_NETCDF4, ncid)
  call check_netcdf_status(nstatus,"open netcdf file")
  nstatus = nf90_enddef(ncid)

! dimensions of the array variables
  nstatus = nf90_def_dim(ncid, 'nparams', nparams, dimid_nparams)
  nstatus = nf90_def_dim(ncid, 'nmascons', nmascons, dimid_nmascons)
  nstatus = nf90_def_dim(ncid, 'nmsc_tidal', nmsc_tidal, dimid_nmsc_tidal)
  nstatus = nf90_def_dim(ncid, 'nICs', nICs, dimid_nICs)
  nstatus = nf90_def_dim(ncid, 'nsolns', nsolns, dimid_nsolns)
  nstatus = nf90_def_dim(ncid, 'nepoch_var', 6, dimid_nepochvar)
  nstatus = nf90_def_dim(ncid, 'version',1,dimid_version)
  
! write the header info of the variables
  ! AtWA components
  nstatus = nf90_def_var(ncid, 'AtWA_mascons', NF90_DOUBLE, [dimid_nmascons,dimid_nmascons], varid_AtWA_mascons)
  call check_netcdf_status(nstatus,"define AtWA_mascons")
  nstatus = nf90_def_var(ncid, 'AtWA_ICs', NF90_DOUBLE, [dimid_nICs,dimid_nICs], varid_AtWA_ICs)
  call check_netcdf_status(nstatus,"define AtWA_ICs")
  nstatus = nf90_def_var(ncid, 'AtWA_mscIC', NF90_DOUBLE, [dimid_nmascons,dimid_nICs], varid_AtWA_mscIC)
  call check_netcdf_status(nstatus,"define AtWA_mscIC")

  ! now the others
  nstatus = nf90_def_var(ncid, 'AtWb', NF90_DOUBLE, [dimid_nparams], varid_AtWb)
  call check_netcdf_status(nstatus,"define AtWb")
  nstatus = nf90_def_var(ncid, 'apriori', NF90_DOUBLE, [dimid_nparams], varid_apriori)
  call check_netcdf_status(nstatus,"define apriori")
  nstatus = nf90_def_var(ncid, 'param_type', NF90_INT2, [dimid_nparams], varid_param_type)
  call check_netcdf_status(nstatus,"define param_type")
  nstatus = nf90_def_var(ncid, 'epoch', NF90_INT2, [dimid_nsolns,dimid_nepochvar], varid_epoch)
  call check_netcdf_status(nstatus,"define epoch")
  nstatus = nf90_def_var(ncid, 'duration', NF90_FLOAT, [dimid_nsolns], varid_duration)
  call check_netcdf_status(nstatus,"define duration")
  nstatus = nf90_def_var(ncid, 'version', NF90_FLOAT, [dimid_version], varid_version)
  call check_netcdf_status(nstatus,"define version")
  if(debug)print*,'written header info'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! prepare "chunking" (a process of compressing the data output in the netcdf file)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nstatus = nf90_def_var_chunking(ncid, varid_AtWA_mascons, NF90_CHUNKED, [10, 101])
  call check_netcdf_status(nstatus,"chunked AtWA_mascons")
  nstatus = nf90_def_var_deflate(ncid, varid_AtWA_mascons,          &
                                shuffle = 1,                &
                                deflate = 1,                &
                                deflate_level = 5  )
  call check_netcdf_status(nstatus,"deflated AtWA_mascons")
  if(debug)print*,'chunked the AtWA_mascons'

  ! PT221125: these other sections are fairly small compared to the mascons, so don't chunk it!
  !nstatus = nf90_def_var_chunking(ncid, varid_AtWA_ICs, NF90_CHUNKED, [10, 101])
  !call check_netcdf_status(nstatus,"chunked AtWA_ICs")
  !nstatus = nf90_def_var_deflate(ncid, varid_AtWA_ICs,          &
  !                              shuffle = 1,                &
  !                              deflate = 1,                &
  !                              deflate_level = 5  )
  !call check_netcdf_status(nstatus,"deflated AtWA_ICs")
  !if(debug)print*,'chunked the AtWA_ICs'

  !nstatus = nf90_def_var_chunking(ncid, varid_AtWA_mscIC, NF90_CHUNKED, [10, 101])
  !call check_netcdf_status(nstatus,"chunked AtWA_mscIC")
  !nstatus = nf90_def_var_deflate(ncid, varid_AtWA_mscIC,          &
  !                              shuffle = 1,                &
  !                              deflate = 1,                &
  !                              deflate_level = 5  )
  !call check_netcdf_status(nstatus,"deflated AtWA_mscIC")
  !if(debug)print*,'chunked the AtWA_mscIC'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! write the attributes
  nstatus = nf90_put_att(ncid, NF90_GLOBAL, 'note', 'netcdf version of normal equations information')
  call check_netcdf_status(nstatus,"note attribute")
  ! the blocks of the AtWA matrix
  nstatus = nf90_put_att(ncid, varid_AtWA_mascons, 'AtWA_mascons', 'unregularised mascons normal equations')
  call check_netcdf_status(nstatus,"AtWA_mascons attribute")
  nstatus = nf90_put_att(ncid, varid_AtWA_ICs, 'AtWA_ICs', 'unregularised ICs normal equations')
  call check_netcdf_status(nstatus,"AtWA_ICs attribute")
  nstatus = nf90_put_att(ncid, varid_AtWA_mscIC, 'AtWA_mscIC', 'unregularised msc/IC normal equations')
  call check_netcdf_status(nstatus,"AtWA_mscIC attribute")

  nstatus = nf90_put_att(ncid, varid_AtWb, 'AtWb', 'RHS')
  call check_netcdf_status(nstatus,"AtWb attribute")
  nstatus = nf90_put_att(ncid, varid_apriori, 'apriori', 'a priori values of the parameters')
  call check_netcdf_status(nstatus,"apriori attribute")
  nstatus = nf90_put_att(ncid, varid_param_type, 'param_type', 'integer code for type of parameter')
  call check_netcdf_status(nstatus,"param_type attribute")
  nstatus = nf90_put_att(ncid, varid_epoch, 'epoch', 'ymdhms of each solution in file')
  call check_netcdf_status(nstatus,"epoch attribute")
  nstatus = nf90_put_att(ncid, varid_duration, 'duration', 'decimal hours of duration of each solution in file')
  call check_netcdf_status(nstatus,"duration attribute")
  nstatus = nf90_put_att(ncid, varid_version, 'version', 'File format version')
  call check_netcdf_status(nstatus,"version attribute")
  if(debug)print*,'written attributes'

! write it out
  call status_update('STATUS',calling_prog,'write_norm_netcdf',netcdf_out,'writing netcdf file',0)
  nstatus = nf90_put_var(ncid, varid_version, version)
  call check_netcdf_status(nstatus,"writing file format version")
  if(debug)print*,'written version'

  ! the blocks of the AtWA matrix
  nstatus = nf90_put_var(ncid, varid_AtWA_mascons, AtWA_mascons)
  call check_netcdf_status(nstatus,"writing AtWA_mascons")
  nstatus = nf90_put_var(ncid, varid_AtWA_ICs, AtWA_ICs)
  call check_netcdf_status(nstatus,"writing AtWA_ICs")
  if(flip)then
    nstatus = nf90_put_var(ncid, varid_AtWA_mscIC, tmp_mscIC)
    call check_netcdf_status(nstatus,"writing AtWA_mscIC transposed")
  else
    nstatus = nf90_put_var(ncid, varid_AtWA_mscIC, AtWA_mscIC)
    call check_netcdf_status(nstatus,"writing AtWA_mscIC")
  endif
  if(debug)print*,'written AtWA blocks'

  nstatus = nf90_put_var(ncid, varid_AtWb, AtWb)
  call check_netcdf_status(nstatus,"writing AtWb")
  if(debug)print*,'written AtWb'

  nstatus = nf90_put_var(ncid, varid_apriori, apriori)
  call check_netcdf_status(nstatus,"writing apriori")
  if(debug)print*,'written apriori'

  nstatus = nf90_put_var(ncid, varid_param_type, param_type)
  call check_netcdf_status(nstatus,"writing param_type")
  if(debug)print*,'written param_type'

  nstatus = nf90_put_var(ncid, varid_epoch, epoch)
  if(debug)print*,'epoch = ',epoch
  call check_netcdf_status(nstatus,"writing epoch info")
  if(debug)print*,'written epoch'

  nstatus = nf90_put_var(ncid, varid_duration, duration)
  call check_netcdf_status(nstatus,"writing soln durations")
  if(debug)print*,'written durations'




! try to write the character array "file_names" to the netcdf file. This code is based on information from
! http://computer-programming-forum.com/49-fortran/af9bfa2653fff83f.htm
  ! define the maximum length of the string
  nstatus = nf90_def_dim (ncid, 'namelen',len(file_names(1)),namelen_dimid)
  ! "same" as for regular dimension of an array
  nstatus = nf90_def_dim (ncid, 'n_strings',size(file_names),n_strings_dimid)
  ! now define char array variable with string length as the first dimension. That is, write it out as though it is a 2-dimensional character array
  nstatus = nf90_def_var(ncid,'file_names',NF90_CHAR &
                        ,dimIDs=(/ namelen_dimid,n_strings_dimid /) &
                        , VarID = file_name_varid )
  ! now write the character array
  nstatus = nf90_put_var(ncid, file_name_varid, file_names)
  call check_netcdf_status(nstatus,"writing soln file names")
  if(debug)print*,'written file_names'
  
  
! write out the satellite names as 1-char character codes (A, B, C, D)
  ! define the maximum length of the string (here it is char*1)
  nstatus = nf90_def_dim (ncid, 'satslen',len(sats(1)),satslen_dimid)
  call check_netcdf_status(nstatus,"writing satellite names 1")
  ! "same" as for regular dimensions
  nstatus = nf90_def_dim (ncid, 'n_strings2',size(sats),n_strings2_dimid)
  call check_netcdf_status(nstatus,"writing satellite names 2")
  ! now define char array variable with string length as the first dimension
  nstatus = nf90_def_var(ncid,'sats',NF90_CHAR &
                        ,dimIDs=(/ satslen_dimid,n_strings2_dimid /) &
                        , VarID = sats_varid )
  call check_netcdf_status(nstatus,"writing satellite names 3")
  ! now write the character array
  nstatus = nf90_put_var(ncid, sats_varid, sats)
  call check_netcdf_status(nstatus,"writing satellite names 4")
  if(debug)print*,'written satellite names'






  


! close the file
  nstatus = nf90_close(ncid)
  call check_netcdf_status(nstatus, 'close netcdf')

  return
  end subroutine write_norm_netcdf_v2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine check_netcdf_status(nstatus, operation)

  use netcdf 

  implicit none

  integer, intent(in) :: nstatus
  character(len=*), intent(in) :: operation

  if (nstatus == NF90_NOERR)then
    return
  else
    print *, "Error encountered during ", operation
    print *, nf90_strerror(nstatus)
    STOP 
  endif

  end subroutine check_netcdf_status
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine reorder_normeq(calling_prog,nparams,nmascons,AtWb,apriori,param_type &
                            ,AtWb_flip,apriori_flip,param_type_flip)

! subroutine to change the order of parameters in the AtWB and parameter lists
! so that mascons come before ICs
!
! P. Tregoning
! 23 November 2022

  implicit none
  
  character(*)  :: calling_prog
  integer*4     :: nparams                       ! total number of parameters
  integer*4     :: nmascons                      ! number of mascon parameters
  real(kind=8)  :: AtWb(nparams)                 ! original AtWb matrix, with ICs then mascons
  real(kind=8)  :: apriori(nparams)              ! original apriori matrix, with ICs then mascons
  integer*4     :: param_type(nparams)           ! original ordering of param_types, with ICs then mascons
  real(kind=8)  :: AtWb_flip(nparams)               ! reordered AtWb matrix, with mascons then ICs
  real(kind=8)  :: apriori_flip(nparams)            ! reordered apriori matrix, with mascons then ICs
  integer*4     :: param_type_flip(nparams)         ! reordered param_types, with mascons then ICs
  
! local variables
  real(kind=8),allocatable :: tmp_matrix(:)
  integer*4,allocatable    :: tmp_matrix_int(:)
  integer*4                :: j,iMsc,nICs
  
! mascons start at row (nparams - nmascons + 1) in original file
  nICs = nparams - nmascons
  iMsc = nICs + 1

  call status_update('STATUS',calling_prog,'reorder_normeq',' ','Reordering the apriori parameter information of normal equations',0)
  allocate(tmp_matrix(nparams))
  tmp_matrix(1:nmascons) = apriori(iMsc:iMsc+nmascons-1)
  tmp_matrix(nmascons+1:nparams) = apriori(1:nICs)
  apriori_flip = tmp_matrix
  deallocate(tmp_matrix)


  call status_update('STATUS',calling_prog,'reorder_normeq',' ','Reordering the parameter type information',0)
  allocate(tmp_matrix_int(nparams))
  tmp_matrix_int(1:nmascons) = param_type(iMsc:iMsc+nmascons-1)
  tmp_matrix_int(nmascons+1:nparams) = param_type(1:nICs)
  param_type_flip = tmp_matrix_int
  deallocate(tmp_matrix_int)
    
  call status_update('STATUS',calling_prog,'reorder_normeq',' ','Reordering AtWb',0)
  allocate(tmp_matrix(nparams))
  tmp_matrix(1:nmascons) = AtWb(iMsc:iMsc+nmascons-1)
  tmp_matrix(nmascons+1:nparams) = AtWb(1:nICs)
  AtWb_flip = tmp_matrix
  deallocate(tmp_matrix)
  
  return
  end subroutine reorder_normeq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




