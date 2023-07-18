    module spherhar_mod

!  Legendre polynomials
!    real (kind=8), save, allocatable :: Plm(:)
!    real(kind=8), save, allocatable :: PDlm(:)
    real (kind=8) :: Plm(1:20301)
    real(kind=8) :: PDlm(1:20301)
    integer*4 :: ndeg
    integer*4 :: jd_old
    real(kind=8) :: t_old
    real(kind=8) :: coefC_old(0:256, 0:256)
    real(kind=8) :: coefS_old(0:256, 0:256)

    integer max_statgrav_degree

! PT170721: logical to indicate whether the atm_tide_sph routine has been called before
    logical :: atm_tide_first_call

    common / maxstatdeg/max_statgrav_degree
    

      end module
