!------------------------------------------------------------------------------
! ANU/GRACE, Softwares for GRACE data
!------------------------------------------------------------------------------
!
! MODULE: addnorm_mod
!
!> @author
!> Paul Tregoning, ANU
!
! DESCRIPTION: 
!> Module containing variable declaration need across addnorm.
!
! REVISION HISTORY:
! 2021-12-06 - Initial Version
! PT220608: added output_correl and use_correl_in_reg

module addnorm_mod

  logical     :: IC_solve          ! estimate ICs?
  logical     :: msc_solve         ! estimate mascons?
  logical     :: tmsc_solve        ! estimate tidal mascons?
  logical     :: output_snorm      ! output stacked normal equations file?
  logical     :: output_SVD        ! output files to perform a SVD of the solution?
  logical     :: output_ascii      ! output ascii versions of AtWA and AtWb ?
  logical     :: output_correl     ! output correlation matrix ?
  logical     :: tikhonov          ! use/don't use tikhonov regularisation
  logical     :: use_correl_in_reg ! use unregularised mascon correlations when building regularisation matrix?
  logical     :: use_adaptive_reg  ! use adaptive regularisation when building regularisation matrix?
  logical     :: norm_netcdf       ! input normal equation files are netcdf format?

  real(kind=8):: msc_scl_factor    ! scale factor with which to multiply the regularisation matrix
  real(kind=8):: cons_mass_sigma   ! uncertainty value (in kg) for the conservation of mass conditional equation
  real(kind=8):: mass_offset       ! value to which conservation of mass should be done (usually zero but some simulations are otherwise)
  real(kind=8):: msc_param_const   ! constraint to be applied to all mascons, in addition to the spatially variable regularisation
  real(kind=8):: lambda            ! scaling factor for tikhonov regularisation
  real(kind=8):: adapt_sigma       ! scaling factor for tikhonov adaptive regularisation 

  integer*4   :: nfiles            ! number of normal equation files read from the file

  character*80  :: msc_constraint_file ! mascon constraint file (read from command file)
  character*100,allocatable :: normeq_files(:)   ! list of normal equation files to be used by addnorm

end module addnorm_mod
