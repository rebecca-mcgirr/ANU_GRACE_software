  module soln_mod
      include '../includes/grace.param'
      include '../includes/write_soln.h90'

! PT190312: declare the maxparm variable (moved from includes/grace.param to this mod file). maxparm will be given
!           a value dynamically - by reading a VCV file - rather than hardwiring it to any particular value.
      integer*4                 :: maxparm

! PT190927: declare a variable of the total number of a priori mascons
      integer*4                 :: total_msc
      integer*4                 :: n_ICs

! PT1902312: make these arrays allocatable
      real(kind=8), allocatable :: soln(:), apriori(:), soln_IC(:),apriori_IC(:),apriori_msc(:,:)  
      real(kind=8), allocatable :: VCV_obs(:,:)

! PT190312: declare the parameter name array
! PT221122: make this allocatable and declare the size in the main program, whichever it is
   integer*4,parameter :: prmnam_size = (6+3+3+3+3)*2 + 46000 + 1000*5*2
!   character*30    :: prmnam(prmnam_size)
   character*30,allocatable :: prmnam(:)

! PT221122: transfer from norm_netcdf_mod.f90 the following arrays required for reading/writing normal equations in netcdf format
  character*100, allocatable :: file_names(:)           ! names of daily files in a .norm or *norm.nc file
  real(kind=8),  allocatable :: AtWA(:,:), AtWb(:)      ! normal equation arrays
  real(kind=8),  allocatable :: duration(:)             ! duration of daily solns in .norm or *norm.nc file
  
  integer*4,     allocatable :: epoch(:,:)              ! (nsolns,6) ymdhms of daily solutions
  integer*4,     allocatable :: param_type(:)           ! array to permit conversion to character parameter descriptions


  end module soln_mod

