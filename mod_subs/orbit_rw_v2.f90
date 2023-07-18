!------------------------------------------------------------------------------
! ANU/GRACE: GRACE Simulation Software
!------------------------------------------------------------------------------
!
! MODULE: orbits_mod
!
!> @author
!> Sebastien Allgeyer, ANU
!
! DESCRIPTION
!> This module deals with the new version of the orbits files (HDF5 format)
!> It is planed that the module will read, write
!
! REVISION HISTORY:
! 2017-08-15 - Initial Version - SA
!
!------------------------------------------------------------------------------

module orbits_rw
    use HDF5  !This is the module for HDF5 writing. 

    implicit none

    type orbits_data
        private
        character (len=80) :: fname
        integer (HID_T) :: file_id !File identifier
        integer(HID_T)  :: coord_id ! id for lat lons
        integer(HID_T)  :: quat_id ! id for quaternions
        integer (HID_T) :: part_id ! id for partials
        integer (HID_T) :: msc_id ! id for mascons partials
        integer (HID_T) :: tmsc_id ! id for tidal mascons partials
        integer (HID_T) :: time_id
        integer (HSIZE_T) :: idx
    end type orbits_data

    

contains
    subroutine orb_create(fname, self)
        character (len=*), intent(in) :: fname
        type(orbits_data), intent(inout):: self
        !integer :: iError
        integer :: date(3), time(3)
        self%fname=trim(fname)
        call idate(date)
        call itime(time)

        call OpenH5file(self%fname, .true., self%file_id)

        self%idx = 0
    end subroutine orb_create


    subroutine orb_open(fname, self)
      character(len=*), intent(in) :: fname
      type(orbits_data), intent(inout) :: self
      self%fname = trim(fname)
      call OpenH5file(self%fname,.false., self%file_id)
    end subroutine orb_open


    subroutine orb_close(self)
        type(orbits_data), intent(inout) :: self
        integer :: ierr
         CALL h5dclose_f(self%coord_id, ierr)
         CALL h5dclose_f(self%quat_id, ierr)
         CALL h5dclose_f(self%part_id, ierr)
         CALL h5dclose_f(self%msc_id, ierr)
         CALL h5dclose_f(self%tmsc_id, ierr)
        CALL h5fclose_f(self%file_id, ierr)
        CALL h5close_f(ierr)
    end subroutine orb_close

  subroutine orb_read(self , time,  coord , sciframe_quat,partials, msc_part, tmsc_part, msc_apr , tmsc_apr)
    type(orbits_data), intent(in) :: self
    double precision, optional, dimension(:,:) :: coord
    double precision, optional, dimension(:,:) :: sciframe_quat
    double precision, optional, dimension(:,:) :: partials
    double precision, optional, dimension(:,:) :: msc_part
    double precision, optional, dimension(:,:) :: tmsc_part
    double precision, optional, dimension(:) :: msc_apr
    double precision, optional, dimension(:,:,:) :: tmsc_apr
    double precision, optional, dimension(:) :: time 
    INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
    integer(hsize_t), dimension(1) :: data1d_dims
    integer(hsize_t), dimension(3) :: data3d_dims
    INTEGER(HID_T) :: dset_id       ! Dataset identifier
    integer :: error


    if (present(coord)) then 
    data_dims(1) = size(coord,1)
    data_dims(2) = size(coord,2)
    CALL h5dopen_f(self%file_id, "coords", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, coord, data_dims, error)
    endif
    if(present(time)) then
    data1d_dims(1) = size(time,1)
    CALL h5dopen_f(self%file_id, "time", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, time, data1d_dims, error)
    endif 

    if(present(sciframe_quat)) then
    data_dims(1) = size(sciframe_quat,1)
    data_dims(2) = size(sciframe_quat,2)
    CALL h5dopen_f(self%file_id, "quat", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, sciframe_quat, data_dims, error)
    endif 

    if (present(partials)) then
    data_dims(1) = size(partials,1)
    data_dims(2) = size(partials,2)
    CALL h5dopen_f(self%file_id, "partials", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, partials, data_dims, error)
    endif

        if (present(msc_part)) then
    data_dims(1) = size(msc_part,1)
    data_dims(2) = size(msc_part,2)
    CALL h5dopen_f(self%file_id, "msc_part", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, msc_part, data_dims, error)
    endif


    if (present(tmsc_part)) then
    data_dims(1) = size(tmsc_part,1)
    data_dims(2) = size(tmsc_part,2)
    CALL h5dopen_f(self%file_id, "tmsc_part", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tmsc_part, data_dims, error)
    endif

    if (present(msc_apr)) then
    data1d_dims(1) = size(msc_apr,1)

    CALL h5dopen_f(self%file_id, "msc_apr", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, msc_apr, data1d_dims, error)
    
    endif 



    if (present(tmsc_apr)) then
    data3d_dims(1) = size(tmsc_apr,1)
    data3d_dims(2) = size(tmsc_apr,2)
    data3d_dims(3) = size(tmsc_apr,3)

    CALL h5dopen_f(self%file_id, "msc_tide_apr", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tmsc_apr, data3d_dims, error)
    
    endif 


  end subroutine orb_read


  subroutine orb_read_in_part(self , start, count, msc_part , partials)
    type(orbits_data), intent(in) :: self
    double precision, optional, dimension(:,:,:) :: partials
    double precision, optional, dimension(:,:,:) :: msc_part
    integer, intent(in) :: start(3)
    integer, intent(in) :: count(3)
    integer(hsize_t) :: start_out(3)
    integer(hsize_t) :: count_out(3)
    integer(hsize_t) :: dim_size(3), start_in(3)
    integer(hsize_t) :: data3d_dims(3)
    integer(hsize_t) :: maxdims(3)
    integer(hsize_t) :: dataspace, memspace
    INTEGER(HID_T) :: dset_id       ! Dataset identifier
    integer :: error
    integer :: memrank
    ! data_dims(1) = size(msc_part,1)
    ! data_dims(2) = size(msc_part,2)

    if (present(msc_part)) then
        CALL h5dopen_f(self%file_id, "msc_part", dset_id, error)
        data3d_dims(1) = size(msc_part,1)
        data3d_dims(2) = size(msc_part,2)
        data3d_dims(3) = size(msc_part,3)
    endif 

    if (present(partials)) then
        CALL h5dopen_f(self%file_id, "partials", dset_id, error)
        data3d_dims(1) = size(partials,1)
        data3d_dims(2) = size(partials,2)
        data3d_dims(3) = size(partials,3)
    endif

    print*, "output mat", data3d_dims
    ! Get the dataspace
    CALL h5dget_space_f(dset_id, dataspace, error)
    ! get the dimension of the dataspace
    CALL H5Sget_simple_extent_dims_f(dataspace, dim_size, maxdims, error)
    
    if (present(partials)) then 
        dim_size(2)=12
    endif

    start_in = (/0,0,0 /) ! Initialise to the begining of the dataspace to read
    CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start_in, dim_size,  error)

    print *, "s, c", start_in, dim_size
    memrank = 3
    CALL h5screate_simple_f(memrank, data3d_dims, memspace, error)
    start_out(1)=start(1)
    start_out(2)=start(2)
    start_out(3)=start(3)
    count_out(1)=count(1)
    count_out(2)=count(2)
    count_out(3)=count(3)
!    print *, "s,c", start_out, count_out
    CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, start_out, count_out, error)

    if (present(partials)) then
        CALL H5dread_f(dset_id, H5T_NATIVE_DOUBLE, partials, data3d_dims, error, memspace, dataspace)
    endif 
    if (present(msc_part)) then
        CALL H5dread_f(dset_id, H5T_NATIVE_DOUBLE, msc_part, data3d_dims, error, memspace, dataspace)
    endif 

  end subroutine orb_read_in_part

    subroutine writeslice(self,  tin, efpos, quat, total_prim, num_mcon_tide,write_msc_partials)
       type (orbits_data), intent(inout) :: self
       integer, intent(in) :: tin
       double precision, intent(in), dimension(:) :: efpos
       double precision, intent(in), dimension(0:3) :: quat
       integer, intent(in) :: total_prim
       integer, intent(in) :: num_mcon_tide
       character*1, intent(in) :: write_msc_partials     ! N or O: don't write them; Y: write them

       INTEGER(HSIZE_T), DIMENSION(1:2) :: count = (/4,3/)  ! Size of hyperslab
       INTEGER(HSIZE_T), DIMENSION(1:2) :: start = (/2,1/) ! Hyperslab offset
       
       INTEGER(HSIZE_T), DIMENSION(1) :: count1= (/1/)  ! Size of hyperslab
       INTEGER(HSIZE_T), DIMENSION(1) :: start1= (/2/) ! Hyperslab offset
       INTEGER(HID_T) :: dataspace     ! Dataspace identifier 
       INTEGER(HID_T) :: memspace      ! memspace identifie
       integer :: ierror


       call h5dget_space_f(self%time_id, dataspace, ierror)
       start1 (1) =  self%idx
       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start1, count1, ierror) 
       call h5screate_simple_f(1,count1,memspace,ierror)
       call h5dwrite_f(self%time_id, H5T_NATIVE_INTEGER, tin , count1 , ierror, memspace, dataspace )
       call h5sclose_f(dataspace,ierror)
       call h5sclose_f(memspace,ierror)

       ! where to write coords
       call h5dget_space_f(self%coord_id, dataspace, ierror)
       count(1:2) = (/6,1/)
       start(1) =0
       start (2) =  self%idx
       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, ierror) 
       call h5screate_simple_f(2,count,memspace,ierror)
       call h5dwrite_f(self%coord_id, H5T_NATIVE_DOUBLE, efpos(1:6) , count , ierror, memspace, dataspace )

       call h5sclose_f(dataspace,ierror)
       call h5sclose_f(memspace,ierror)

       call h5dget_space_f(self%quat_id, dataspace, ierror)
       count(1:2) = (/4,1/)
       start(1) =0
       start (2) =  self%idx
       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, ierror) 
       call h5screate_simple_f(2,count,memspace,ierror)
       call h5dwrite_f(self%quat_id, H5T_NATIVE_DOUBLE, quat , count , ierror, memspace, dataspace )
       call h5sclose_f(dataspace,ierror)
       call h5sclose_f(memspace,ierror)

!partials
       if (size(efpos) > 6 ) then 
       call h5dget_space_f(self%part_id, dataspace, ierror)
       count(1:2) = (/114,1/)
       start(1) =0
       start (2) =  self%idx
       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, ierror) 
       call h5screate_simple_f(2,count,memspace,ierror)
       call h5dwrite_f(self%part_id, H5T_NATIVE_DOUBLE, efpos(7:7+114) , count , ierror, memspace, dataspace )
       call h5sclose_f(dataspace,ierror)
       call h5sclose_f(memspace,ierror)
       endif

!msc
       if ( (write_msc_partials == "Y") .and. ( total_prim /= 0 ) .and. (size(efpos)>6 )) then

       call h5dget_space_f(self%msc_id, dataspace, ierror)
       count(1:2) = (/total_prim*6,1/)
       start(1) =0
       start (2) =  self%idx
       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, ierror) 
       call h5screate_simple_f(2,count,memspace,ierror)
       call h5dwrite_f(self%msc_id, H5T_NATIVE_DOUBLE, efpos(7+114:7+114+total_prim*6) , count , ierror, memspace, dataspace )
       call h5sclose_f(dataspace,ierror)
       call h5sclose_f(memspace,ierror)
       endif

!tmsc

     if (num_mcon_tide /= 0 ) then

        call h5dget_space_f(self%tmsc_id, dataspace, ierror)
       count(1:2) = (/num_mcon_tide*6*2,1/)
       start(1) =0
       start (2) =  self%idx
       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, ierror)
       call h5screate_simple_f(2,count,memspace,ierror)
       call h5dwrite_f(self%tmsc_id, H5T_NATIVE_DOUBLE, efpos(7+114+total_prim*6:7+114+total_prim*6+num_mcon_tide*12) , &
             count , ierror, memspace, dataspace )
        call h5sclose_f(dataspace,ierror)
       call h5sclose_f(memspace,ierror)
    endif



       self%idx = self%idx + 1
   end subroutine writeslice


   subroutine h5_write_apr(orb_file, msc_tide, msc)
      implicit none

      type(orbits_data), intent(in) :: orb_file
      double precision, optional, dimension(:,:,:), intent(in) :: msc_tide
      double precision, optional, dimension(:), intent(in) :: msc


      INTEGER :: hdferr
      INTEGER(HID_T) :: file, space, dset ! handles
      INTEGER(HSIZE_T), DIMENSION(1:1)           :: dims2 ! size buffer
      INTEGER(HSIZE_T), DIMENSION(1:3)           :: dims3 ! size buffer

      if(present(msc_tide)) then
        dims3 (3) = size(msc_tide,3)
        dims3 (2) = size(msc_tide,2)
        dims3 (1) = size(msc_tide,1)
        CALL h5screate_simple_f(3, dims3, space, hdferr)
        CALL h5dcreate_f(orb_file%file_id, "msc_tide_apr", H5T_NATIVE_DOUBLE, space,  dset, hdferr)
        CALL h5dwrite_f(dset, H5T_NATIVE_DOUBLE, msc_tide, dims3, hdferr)
        CALL h5dclose_f(dset, hdferr)
      endif

      if(present(msc)) then
        dims2 (1) = size(msc,1)
        CALL h5screate_simple_f(1, dims2, space, hdferr)
        CALL h5dcreate_f(orb_file%file_id, "msc_apr", H5T_NATIVE_DOUBLE, space, dset, hdferr)
        CALL h5dwrite_f(dset, H5T_NATIVE_DOUBLE, msc, dims2, hdferr)
        CALL h5dclose_f(dset, hdferr)
      endif

   end subroutine h5_write_apr

    subroutine hdwrt( orb_file ,tin,date,sectag, tinend, dateend, sectagend, tim, step, mjd, efic, incoor, time_array, frame, &
                     GPSant_mag, GPSant_cos,gpsantoff)

! writes header file for GRACE orbit integrator output
!
! MODS
! PT130902: write out as an unformatted, binary file instead of ascii
! PT170609: add inclusion of file gtorb_mod (contains record length of GTORB file) and write recl on first line of GTORB file
    use sat_mod
    use inmod_mod
    use bsscl_mod
    use accel_mod
    use mascon_mod
  !  use gtorb_mod    ! PT170609: this declares the record length information for the GTORB file
  !  use orbits_rw

    implicit none
    type(orbits_data), intent(inout) :: orb_file
    integer :: headlines
    integer*4 :: datarecs
    integer*4 mjd
    integer :: time_array(8)
    integer*4 :: step
!    integer*4 :: num_mcon
    integer*4 :: date(5), dateend(5)
    real(kind=8) :: tim, tin, tinend
    real*8 :: sectag, sectagend
!    real*8 :: bs(3)
    real(kind=8), dimension(6) :: incoor
    real(kind=8), dimension(10) :: efic      ! earth-fixed ICs (pos, vel, twice-per-rev along-track acc
    character*11 :: coorspace, partspace
    !character (len = *):: inputfile                ! PT170127: increased from C*20
    character*10 :: fileform
    character*10 :: agency, uname
    character (len=*), intent(in) :: frame
    real(kind=8) :: GPSant_mag
    real(kind=8), dimension(3) :: GPSant_cos, gpsantoff

    INTEGER(HSIZE_T), DIMENSION(1:2) :: dimsm = (/4,3/) ! dimension of the array
    integer :: rank, ierror
    integer (HID_T) :: dataspace, data_id
    integer , dimension(7) :: date_vec

    softversion = "alpha"      ! software version
    fileform = 'hdf5_v0.1'           ! output format of file
    headlines = 24             ! number of header lines
    datarecs = step+1          ! number of data records
    coorspace = 'terrestrial'  ! coordinates are in 'terrestrial' or 'inertial' space
    partspace = 'terrestrial'  ! IC partials are in 'terrestrial' or 'inertial' space
! PT140616: change to the actual name (rather than hardwired 'GRACE.input')
 !   inputfile = input_file  ! input files with models and components etc



    call getenv('USER',uname)  ! user name running programme
    call getenv('INSTITUTE',agency) ! agency where code is run

    call WriteAttribute(orb_file, "Agency", 1, strScalar= agency)
    call WriteAttribute(orb_file,"Created by", 1, strScalar=uname)

    write (message, 145) time_array(1),time_array(2),time_array(3),time_array(5),time_array(6), time_array(7)
145 format (i4,'-',i2.2,'-',i2.2,' ',i2.2,':',i2.2,':',i2.2 )
    call WriteAttribute(orb_file, "Creation date",1,strScalar=message) 
    call WriteAttribute(orb_file, "Software version",1,strScalar=softversion) 
    call WriteAttribute(orb_file, "File format",1,strScalar=fileform) 

    call WriteAttribute(orb_file, "Satellite name",1,strScalar=sat) 
    call WriteAttribute(orb_file, "NRec",1,IntegerScalar=datarecs)
    !call WriteAttribute(orb_file%file_id, "NRec",1,RealScalar=tim)

    date_vec(1)=0
    date_vec(2)=2000
    date_vec(3)=1
    date_vec(4)=1
    date_vec(5)=12
    date_vec(6)=0
    date_vec(7)=0
    call WriteAttribute(orb_file, "Time epoch (GPS TIME)",7,IntegerArray=date_vec)
    date_vec(1)=nint(tin)
    date_vec(2)=date(1)
    date_vec(3)=date(2)
    date_vec(4)=date(3)
    date_vec(5)=date(4)
    date_vec(6)=date(5)
    date_vec(7)=int(sectag)
    call WriteAttribute(orb_file, "Time first obs (Sec past epoch)",7,IntegerArray=date_vec)
    date_vec(1)=nint(tinend)
    date_vec(2)=dateend(1)
    date_vec(3)=dateend(2)
    date_vec(4)=dateend(3)
    date_vec(5)=dateend(4)
    date_vec(6)=dateend(5)
    date_vec(7)=int(sectagend)
    call WriteAttribute(orb_file, "Time last obs (Sec past epoch)",7,IntegerArray=date_vec)
    call WriteAttribute(orb_file, "Sampling", 1, RealScalar = tim)

    if (gt_statfield(1).eq."Y") then
        call WriteAttribute(orb_file, "Static_field",1 , strScalar=gt_statfieldmod(1)) 
    else
        call WriteAttribute(orb_file, "Static_field", 1, strScalar="None") 
    endif

    if (gt_oceantide(1).eq."Y") then
        call WriteAttribute(orb_file, "Tide_model",1,  strScalar=gt_oceantidemod(1)) 
    else
        call WriteAttribute(orb_file, "Tide_model", 1,  strScalar="None") 
    endif

    if (gt_atmtide(1).eq."Y") then
        call WriteAttribute(orb_file, "Atm_tide_model",1,  strScalar=gt_atmtidemod(1)) 
    else
        call WriteAttribute(orb_file, "Atm_tide_model",1,  strScalar="None") 
    endif

    if (gt_dealias(1).eq."Y") then
        call WriteAttribute(orb_file, "dealiasing",1 , strScalar=gt_dealiasmod(1)) 
    else
        call WriteAttribute(orb_file, "dealiasing", 1 , strScalar="None") 
    endif

    call WriteAttribute(orb_file, "dealiasing", 1, strScalar=gt_dealiasmod(1)) 
    call WriteAttribute(orb_file, "reference frame",1, strScalar=coorspace) 
    call WriteAttribute(orb_file, "t_in",1,  IntegerScalar=nint(tin)) 
    call WriteAttribute(orb_file, "ICs terrestrial",10 ,  RealArray=efic) 
    call WriteAttribute(orb_file, "ICs inertial",10,  RealArray=incoor) 
    call WriteAttribute(orb_file, "Acc Scale", 3,  RealArray=scl) 
    call WriteAttribute(orb_file, "Acc Bias",3, RealArray=bs*1.0d6) 

    call WriteAttribute(orb_file, "GPS mag",1,  RealScalar=GPSant_mag) 
    call WriteAttribute(orb_file, "GPS cos",3,  RealArray=GPSant_cos) 
    call WriteAttribute(orb_file, "GPS off",3,  RealArray=gpsantoff) 


    call WriteAttribute(orb_file, "NMSC",1, IntegerScalar=total_prim) 
    call WriteAttribute(orb_file, "NTMSC",1,IntegerScalar=num_mcon_tide) 
    call WriteAttribute(orb_file, "MSC fname", 1, strScalar=combined_mascon_file)
    call WriteAttribute(orb_file, "MSC code", 1,  strScalar = msc_hdr_code)
    call WriteAttribute(orb_file, "Ocean fname", 1,  strScalar = ocean_mascon_file)
    call WriteAttribute(orb_file, "Ocean Code", 1,  strScalar = msc_ocn_hdr_code)
    call WriteAttribute(orb_file, "Input File", 1, strScalar = input_file )


    rank = 1 
    dimsm(1) = total_ocean_prim
    CALL h5screate_simple_f(rank, dimsm, dataspace, ierror)
    CALL h5dcreate_f(orb_file%file_id, "bitmap_ocean", H5T_NATIVE_INTEGER, dataspace, data_id, ierror)
    CALL h5dwrite_f(data_id, H5T_NATIVE_INTEGER, mcon_ocean_prim(:,2), dimsm, ierror)
                           CALL h5dclose_f(data_id, ierror)

    dimsm(1) = datarecs
    CALL h5screate_simple_f(rank, dimsm, dataspace, ierror)
    CALL h5dcreate_f(orb_file%file_id, "time", H5T_NATIVE_INTEGER, dataspace, orb_file%time_id, ierror)

    rank = 2
    dimsm (2) = datarecs
    dimsm (1) = 6
    CALL h5screate_simple_f(rank, dimsm, dataspace, ierror)
    CALL h5dcreate_f(orb_file%file_id, "coords", H5T_NATIVE_DOUBLE, dataspace, orb_file%coord_id, ierror)

    dimsm (2) = datarecs
    dimsm (1) = 4
    CALL h5screate_simple_f(rank, dimsm, dataspace, ierror)
    CALL h5dcreate_f(orb_file%file_id, "quat", H5T_NATIVE_DOUBLE, dataspace, orb_file%quat_id, ierror)

    dimsm (2) = datarecs
    dimsm (1) = 114
    CALL h5screate_simple_f(rank, dimsm, dataspace, ierror)
    CALL h5dcreate_f(orb_file%file_id, "partials", H5T_NATIVE_DOUBLE, dataspace, orb_file%part_id, ierror)

    dimsm(2) = datarecs
    dimsm(1) = 6*total_prim
    CALL h5screate_simple_f(rank, dimsm, dataspace, ierror)
    CALL h5dcreate_f(orb_file%file_id, "msc_part", H5T_NATIVE_DOUBLE, dataspace, orb_file%msc_id, ierror)   

    dimsm(2) = datarecs
    dimsm(1) = 6*num_mcon_tide
    CALL h5screate_simple_f(rank, dimsm, dataspace, ierror)
    CALL h5dcreate_f(orb_file%file_id, "tmsc_part", H5T_NATIVE_DOUBLE, dataspace, orb_file%tmsc_id, ierror) 
    return 
    end

    subroutine hdread( orb_file , Nrec, statics_t, apr_scalebias, GPSantoff, combined_mascon_file, msc_hdr_code, total_prim, &
     ocean_mascon_file, msc_ocn_hdr_code, total_ocean_prim,  display)


    implicit none
    type(orbits_data), intent(inout) :: orb_file
    integer, intent(out) :: Nrec ! Number of data
    double precision, intent(out), dimension(10) :: statics_t
    double precision, intent(out), dimension(6) :: apr_scalebias
    double precision, intent(out), dimension(7) :: GPSantoff
    character(len=*), intent(out) :: combined_mascon_file
    character(len=*), intent(out) :: msc_hdr_code
    character(len=*), intent(out) :: ocean_mascon_file
    character(len=*), intent(out) :: msc_ocn_hdr_code
    integer, intent(out) :: total_prim
    integer, intent(out) :: total_ocean_prim
    logical, optional, intent(in) :: display

   ! integer*4 :: i
   ! integer :: headlines
   ! integer*4 mjd
    !integer :: time_array(8)
   ! integer*4 :: step
!    integer*4 :: num_mcon
  !  integer*4 :: date(5), dateend(5)
    !real(kind=8) :: tim, tinend
    integer :: tin
   ! real*8 :: sectag, sectagend
!    real*8 :: bs(3)
    real(kind=8), dimension(6) :: incoor
   ! real(kind=8), dimension(10) :: efic      
    !character*11 :: coorspace, partspace
    character*10 :: fileform
    character*10 :: agency, uname
    !character (len=128) :: frame
    character (len=128) :: softversion
    character (len=1) :: sat
   ! INTEGER(HSIZE_T), DIMENSION(1:2) :: dimsm = (/4,3/) ! dimension of the array
    !integer :: rank, ierror

    character(len=256) :: create_date, GPS_epoch, Firstobs, LastObs


    call ReadAttribute(orb_file, "Agency", 1, strScalar= agency)
    call ReadAttribute(orb_file,"Created by", 1, strScalar=uname)
    call ReadAttribute(orb_file, "Creation date",1,strScalar=create_date) 
    call ReadAttribute(orb_file, "Software version",1,strScalar=softversion) 
    call ReadAttribute(orb_file, "File format",1,strScalar=fileform) 
    call ReadAttribute(orb_file, "Satellite name",1,strScalar=sat) 
    call ReadAttribute(orb_file, "NRec",1,IntegerScalar=Nrec)
    call ReadAttribute(orb_file, "Time epoch (GPS TIME)",1,strScalar=GPS_epoch) 
    call ReadAttribute(orb_file, "Time first obs (Sec past epoch)",1,strScalar=Firstobs) 
    call ReadAttribute(orb_file, "Time last obs (Sec past epoch)",1,strScalar=LastObs) 
    !call ReadAttribute(orb_file, "Static_field",1 , strScalar=gt_statfieldmod(1)) 
    !call ReadAttribute(orb_file, "Tide_model",1,  strScalar=gt_oceantidemod(1)) 
    !call ReadAttribute(orb_file, "Atm_tide_model",1,  strScalar=gt_atmtidemod(1)) 
    !call ReadAttribute(orb_file, "dealiasing",1 , strScalar=gt_dealiasmod(1)) 
    !call ReadAttribute(orb_file, "dealiasing", 1, strScalar=gt_dealiasmod(1)) 
    !call ReadAttribute(orb_file, "reference framce",1,strScalar=coorspace) 
    call ReadAttribute(orb_file, "t_in",1,  IntegerScalar=tin) 
    call ReadAttribute(orb_file, "ICs terrestrial",10 ,  RealArray=statics_t) 
    call ReadAttribute(orb_file, "ICs inertial",10, RealArray=incoor) 
    call ReadAttribute(orb_file, "Acc Scale", 3, RealArray=apr_scalebias(1:3)) 
    call ReadAttribute(orb_file, "Acc Bias",3, RealArray=apr_scalebias(4:6)) 
    call ReadAttribute(orb_file, "GPS mag",1, RealScalar=GPSantoff(1)) 
    call ReadAttribute(orb_file, "GPS cos",3, RealArray=GPSantoff(2:4)) 
    call ReadAttribute(orb_file, "GPS off",3, RealArray=GPSantoff(5:7)) 
    call ReadAttribute(orb_file, "NMSC",1,  IntegerScalar=total_prim) 
    call ReadAttribute(orb_file, "NTMSC",1,  IntegerScalar=total_ocean_prim) 
    call ReadAttribute(orb_file, "MSC fname", 1,  strScalar=combined_mascon_file)
    call ReadAttribute(orb_file, "MSC code", 1, strScalar = msc_hdr_code)
    call ReadAttribute(orb_file, "Ocean fname", 1,   strScalar = ocean_mascon_file)
    call ReadAttribute(orb_file, "Ocean Code", 1,   strScalar = msc_ocn_hdr_code)

    if(present(display)) then
      if (display) then
        print *, "***********************  File INFO ***********************************"
        print *, "Agency", agency
        print *, "Created by", uname, " on ", create_date
        print *, "Software version", softversion, " File version", fileform
        print *, "Satellite name", sat
      endif
    endif


    end subroutine hdread





    subroutine OpenH5file(FileName, create, id)
        Implicit none
        !Input variables
        character(len=*), intent(in) :: FileName
        logical, intent(in) :: create
        ! Output variables
        Integer (HID_T), intent(out) :: id
        !Local varaibles
        integer iError
        INTEGER(HID_T) :: dset_id       ! Dataset identifier 
        INTEGER(HID_T) :: dataspace     ! Dataspace identifie
        !Initialise HDF5
        call h5open_f(iError)
        if(create)   then
            call h5fcreate_f(trim(FileName),H5F_ACC_TRUNC_F, id, iError)
            call h5dcreate_f(id, "IC", H5T_NATIVE_INTEGER, dataspace, dset_id, iError)
            call h5dcreate_f(id, "MASCONS", H5T_NATIVE_INTEGER, dataspace, dset_id, iError)
        else
            call h5fopen_f(trim(FileName), H5F_ACC_RDONLY_F, id,iError)
        end if
    end subroutine OpenH5file


    subroutine WriteAttribute(orb_file, Att_Name, nVal, DataSetName, RealScalar, IntegerScalar, strScalar, RealArray, IntegerArray )
        implicit none

        !input variables
        type(orbits_data), intent(inout) :: orb_file
        Integer(HID_T)  :: File_id
        Character (len=*), intent(in) :: Att_Name
        Integer, intent(in) :: nVal
        Character (len=*), optional, intent(in) :: DataSetName
        double precision, optional, intent(in)  :: RealScalar
        Integer, optional, intent(in) ::  IntegerScalar
        Character (len=*), optional, intent(in) ::  strScalar
        double precision, optional, intent(in) :: RealArray(nVal)
        Integer, optional, intent(in) :: IntegerArray(nval)
        !output
        ! NAN

        !local
        integer :: Rank, ierror
        Integer(HID_T) :: DataSpace, Attr_ID, Loc_ID, aType_ID
        integer(HSIZE_T), dimension(1) :: Dimsf
        Integer(SIZE_T)               :: attr_len

        file_id = orb_file%file_id
        !if DataSetName
        if(present(DataSetName)) then
            call h5dopen_f(File_id, TRIM(DataSetName), Loc_id, iError)
        else 
            Loc_id = File_id
        end if

        Rank = 1
        Dimsf(:) = 0
        Dimsf(1) = nVal
        call h5screate_simple_f(Rank, Dimsf, DataSpace, iError)

        if (present(IntegerArray)) then
            call H5ACREATE_F(Loc_id, trim(Att_Name), H5T_NATIVE_INTEGER, Dataspace, Attr_ID, iError)
            call H5AWRITE_F(Attr_ID,H5T_NATIVE_INTEGER, IntegerArray, Dimsf, iError)
        endif
        if (present(RealArray)) then
            call H5ACREATE_F(Loc_id, trim(Att_Name), H5T_NATIVE_DOUBLE, Dataspace, Attr_ID, iError)
            call H5AWRITE_F(Attr_ID,H5T_NATIVE_DOUBLE, RealArray, Dimsf, iError)
        endif
        if (present(IntegerScalar)) then
            call H5ACREATE_F(Loc_id, trim(Att_Name), H5T_NATIVE_INTEGER, Dataspace, Attr_ID, iError)
            call H5AWRITE_F(Attr_ID,H5T_NATIVE_INTEGER, IntegerScalar, Dimsf, iError)
        endif
        if (present(RealScalar)) then
            call H5ACREATE_F(Loc_id, trim(Att_Name), H5T_NATIVE_DOUBLE, Dataspace, Attr_ID, iError)
            call H5AWRITE_F(Attr_ID,H5T_NATIVE_DOUBLE, RealScalar, Dimsf, iError)
        endif
        if(present(strScalar)) then
            call H5TCOPY_F(H5T_NATIVE_CHARACTER, aType_ID, iError)
            attr_len = len(trim(strScalar))
            call H5TSET_SIZE_F(aType_ID, attr_len, iError)
            call H5ACREATE_F(Loc_id, trim(Att_Name), aType_Id, DataSpace, Attr_ID, iError)
            call H5AWRITE_F(Attr_ID, aType_Id, trim(strScalar), Dimsf, iError)
        endif

        call h5sclose_f(Dataspace, iError)
        call h5aclose_f(Attr_ID, iError)
        if (present(DataSetName)) then
            call h5dclose_f(Loc_id, iError)
        endif
    end subroutine WriteAttribute


    subroutine ReadAttribute(orb_file, Att_Name, nVal, DataSetName, RealScalar, IntegerScalar, strScalar, RealArray, IntegerArray )
        implicit none
        type(orbits_data), intent(inout) :: orb_file

        !input variables
        Integer(HID_T)    :: File_id
        Character (len=*), intent(in) :: Att_Name
        Integer, intent(in) :: nVal
        Character (len=*), optional, intent(in) :: DataSetName

        !output
        ! NAN
        double precision, optional, intent(out)  :: RealScalar
        Integer, optional, intent(out) ::  IntegerScalar
        Character (len=*), optional, intent(out) ::  strScalar
        double precision, optional, intent(out) :: RealArray(nVal)
        Integer, optional, intent(out) :: IntegerArray(nval)
        character (len=200) :: tmp

        !local
        Integer(HID_T) :: Attr_ID, Loc_ID, Type_ID
        integer(HSIZE_T), dimension(1) :: Dimsf
        integer :: iError
        File_id = orb_file%file_id

        !print *, "try to get attribute", Att_Name

        Dimsf(1)=nVal
        !if DataSetName
        if(present(DataSetName)) then
            call h5dopen_f(File_id, TRIM(DataSetName), Loc_id, iError)
        else 
            Loc_id = File_id
        end if
        !open the attribute
        CALL H5AOPEN_F(Loc_ID, TRIM(Att_Name), Attr_ID, iError) 

        ! read the attribute data.
        if(present(RealArray))then
          RealArray=0.
          CALL H5AREAD_F(Attr_ID, H5T_NATIVE_DOUBLE, RealArray, Dimsf, iError)
        endif 
        IF(PRESENT(RealScalar))THEN
        RealScalar=0.
        CALL H5AREAD_F(Attr_ID, H5T_NATIVE_DOUBLE, RealScalar, Dimsf, iError)
        END IF
        IF(PRESENT(IntegerArray))THEN
        IntegerArray=0
        CALL H5AREAD_F(Attr_ID, H5T_NATIVE_INTEGER, IntegerArray, Dimsf, iError)
        END IF
        IF(PRESENT(IntegerScalar))THEN
        IntegerScalar=0
        CALL H5AREAD_F(Attr_ID, H5T_NATIVE_INTEGER, IntegerScalar, Dimsf, iError)
        END IF
        IF(PRESENT(StrScalar))THEN
          tmp=''
          CALL H5AGET_TYPE_F(Attr_ID, Type_ID, iError)  ! Get HDF5 data type for character string
          CALL H5AREAD_F(Attr_ID, Type_ID, tmp, Dimsf, iError)
          CALL H5TCLOSE_F(Type_ID, iError)
          StrScalar = trim(tmp)
        END IF

      IF(Loc_ID.NE.File_id)THEN
        ! Close the dataset and property list.
        CALL H5DCLOSE_F(Loc_ID, iError)
      END IF

      end subroutine ReadAttribute  


end module orbits_rw
