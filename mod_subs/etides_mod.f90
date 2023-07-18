  module etides_mod

! declaration of variables related to the solid body tides
!
! P. Tregoning
! 29 October 2018

! IERS2010 Love number variables
  real(kind=8) :: k2mr(3)                                ! degree 2, anelastic real component
  real(kind=8) :: k2mi(3)                                ! degree 2, anelastic imaginary component
  real(kind=8) :: k2mp(3)                                ! degree 2, anelastic real component for deg 4 effect
  real(kind=8) :: k3m(4)                                 ! degree 3, elastic real component

! IERS 2010 variables for the number of tidal constituents
  integer*4    :: num_deg20                              ! number of tidal components for degree 2,0
  integer*4    :: num_deg21                              ! number of tidal components for degree 2,1
  integer*4    :: num_deg22                              ! number of tidal components for degree 2,2

! IERS 2010 coefficients for C20, deg 21 and deg 22
  character*7 ,allocatable          :: ztid_doodson(:)   ! array for Doodson numbers for C20
  character*7 ,allocatable          :: dtid_doodson(:)   ! array for Doodson numbers for degree 2,1
  character*7 ,allocatable          :: sdtid_doodson(:)  ! array for Doodson numbers for degree 2,2

  real(kind=8),allocatable          :: ztidip(:)         ! C20 in-phase amplitudes
  real(kind=8),allocatable          :: ztidop(:)         ! C20 out-of-phase amplitudes
  real(kind=8),allocatable          :: dtidip(:)         ! diurnal in-phase amplitudes
  real(kind=8),allocatable          :: dtidop(:)         ! diurnal out-of-phase amplitudes
  real(kind=8),allocatable          :: sdtidip(:)        ! semi-diurnal in-phase amplitudes
  real(kind=8),allocatable          :: sdtidop(:)        ! semi-diurnal out-of-phase amplitudes

! variables for the Legendre polynomials for the tidal potential
! PT210122: increased from 3 and 5 to 4 and 9, to make compatible with iers2010 solid body tide code
  real(kind=8) :: Pn_mon(4),Pnm_mon(9),Pn_sun(4),Pnm_sun(9) 



  end module etides_mod


