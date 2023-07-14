subroutine output_3d(fname, dim1, dim2, dim3, data)
  use netcdf
  ! Inputs
  integer, intent(in) :: dim1, dim2, dim3
  double precision, intent(in), dimension(dim1,dim2,dim3) :: data
  character(len=*), intent (in) :: fname
  !NetCDF vars
   integer :: ncid, dimid1, dimid2, dimid3
   integer, dimension(3) :: sdimid
   integer :: varid
   call check( nf90_create(trim(fname), NF90_CLOBBER, ncid) )
   call check( nf90_def_dim(ncid, "dim1", dim1, dimid1) )
   call check( nf90_def_dim(ncid, "dim2", dim2, dimid2) )
   call check( nf90_def_dim(ncid, "dim3", dim3, dimid3) )
   sdimid = (/dimid1, dimid2, dimid3/)
   call check( nf90_def_var(ncid, "data", NF90_Double, sdimid, varid ))
   call check( nf90_enddef(ncid) )

   call check( nf90_put_var(ncid, varid, data))
   call check( nf90_close(ncid) )
   print *, " *** Success in writing file *** ", fname
end subroutine output_3d

subroutine output_2d(fname, dim1, dim2,  array)
  use netcdf
  ! Inputs
  integer, intent(in) :: dim1, dim2
  double precision, intent(in), dimension(dim1,dim2) :: array
  character(len=*), intent (in) :: fname
  !NetCDF vars
   integer :: ncid, dimid1, dimid2, dimid3
   integer, dimension(2) :: sdimid
   integer :: varid
   call check( nf90_create(trim(fname), OR(NF90_NETCDF4, NF90_CLOBBER) , ncid) )
   call check( nf90_def_dim(ncid, "dim1", dim1, dimid1) )
   call check( nf90_def_dim(ncid, "dim2", dim2, dimid2) )
   sdimid = (/dimid1, dimid2/)
   call check( nf90_def_var(ncid, "data", NF90_Double, sdimid, varid, deflate_level=4 ))
   call check( nf90_enddef(ncid) )
   call check( nf90_put_var(ncid, varid, array))
   call check (nf90_close(ncid))
   print *, " *** Success in writing file *** ", fname
end subroutine output_2d


subroutine check(status)
  use netcdf
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop 2
    end if
end subroutine check
