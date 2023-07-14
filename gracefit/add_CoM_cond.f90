!********************************************************************************************************************************

subroutine add_CoM_cond(part, pre_omc, rvec, gvec, gps_antoff_trf, srf2trf_rotmat, srf2trf_deriv_rotmat)

  ! subroutine to impose conditional equations on the system:
  !
  ! The difference between the GPS derived position and velocity and the theoretical position and
  ! velocity for each satellite is assumed to be due to the CoM offset. Here we use the values stored
  ! in the gn_part array to set the corresponding partials for use in the normal equations.
  !
  ! Each satellite has three possible components of CoM offset (x, y, and z)
  ! These values are stored in the part() matrix in columns 13, 14, and 15 for satellite A, and 28, 29, and 30 for satellite B
  !
  ! The conditional equation observations are stored immediately after the GPS and K-band observations
  ! There are twelve such observations, position and velocity in three dimensions for each of two satellites
  !
  ! A. Purcell (Modified by Thomas Greenspan)
  ! 25 June 2013
  use gracefit_mod
  implicit none



  double precision, intent(out) :: part(nobs_t,nparam_t)             ! Partials of observables with respect to the parameters
  double precision, intent(out) :: pre_omc(nobs_t)                   ! Observed Minus Computed values to be used in LS algorithm
  double precision, intent(in)  :: srf2trf_rotmat(3,3,nsat_t)        ! SRF to TRF rotation matrix for GRACE A and B
  double precision, intent(in)  :: srf2trf_deriv_rotmat(3,3,nsat_t)  ! Differentiated SRF to TRF rotation matrix for GRACE A and B
  double precision, intent(in)  :: gvec(maxgprm,nsat_t)              ! observed GPS positions and velocities
  double precision, intent(in)  :: gps_antoff_trf(6,nsat_t)          ! CoM offset of GPS antenna in the terrestrial reference frame
  double precision, intent(in)  :: rvec(nCoM_t,nsat_t)               ! Vector of positions, velocities and partials

  integer*4 :: isat      ! Counter variable for satellites
  integer*4 :: i,j       ! Counter variables
  !=====================================================================================
  !     The partials are just the SRF - TFR rotation matrix and its derivative
  !=====================================================================================
  do isat = 1, nsat_t
     if(gvec(4,isat) == 0.d0) cycle ! Make sure given epoch existed in GNV1B file
     do i = 1, 3                ! Loop over x, y, and z components of CoM offset for current satellite
        do j= 1, nCoM_t - 3
           part(j+iCoM-1+(isat-1)*nCoM_t,   iantoff-1+i+nsatprm_t*(isat-1)) = srf2trf_rotmat(i,j,isat)
           part(j+iCoM-1+3+(isat-1)*nCoM_t, iantoff-1+i+nsatprm_t*(isat-1)) = srf2trf_deriv_rotmat(i,j,isat)
        enddo
     enddo
     ! Calculate the pre_omc vector for the conditional equation
     do i = 1, nCoM_t
        pre_omc(i+iCoM-1+(isat-1)*nCoM_t) = gvec(i,isat) - rvec(i,isat) + gps_antoff_trf(i,isat)
     enddo
  enddo

  return
end subroutine add_CoM_cond!

!********************************************************************************************************************************
