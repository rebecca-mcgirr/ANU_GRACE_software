! Module to store raw accelerometer and star camera data
  module accred_mod
     integer*4 :: ACCn 
     integer*4 :: STARn 
     real(kind=8) :: sca_step  ! time interval between SCA1B observations (5s for GRACE, 1s for GRACE FO)
     real(kind=8), dimension(:), allocatable :: ACCtime, SCAtime
     real(kind=8), dimension(:,:), allocatable :: accr, varr 
     real(kind=8), dimension(:,:), allocatable :: sca_obs
     real(kind=8), dimension(:,:,:), allocatable :: acc_obs

! PT190212: define here accel offset values (used to reduce the uncalibrated obs down to the nm/s^2 level)
    real(kind=8) :: acc_offset(4,3)                   ! each satellite (GRACE A/B, C,D) has an offset for each orthogonal direction.

! PT191104: add another star camera array declaration to this mod file. This is needed to determine the pitch
!           corrections for GRACE FO (we need simultaneously the information from both star cameras)
     real(kind=8), dimension(:,:), allocatable :: sca_obs_2

  end module accred_mod

