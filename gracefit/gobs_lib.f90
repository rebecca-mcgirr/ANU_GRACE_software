!********************************************************************************************************************************
!  File: gobs_lib.f90
!
!  Purpose: Set of subroutines concerning gobs calculations.
!
!  Author: Thomas Greenspan
!
!  API:
!       gobs_computeOMC       : Computes the omc (observed minus calculated)
!       gobs_computePartials  : Computes the partials
!
!   July 23, 2013
!
!********************************************************************************************************************************
!********************************************************************************************************************************
! gobs_computeOMC: Computes the prefit omc (observed minus calculated) given the observed pos, vel 
!                  and accel, the calculated pos, vel and accel and the GPS antenna offset in the 
!                  terrestrial reference frame. 
! Author: Thomas Greenspan
!********************************************************************************************************************************

subroutine gobs_computeOMC(pre_omc,gvec,rvec,gps_antoff_trf)

  use gracefit_mod
  implicit none

   

  !********************  Variable declarations ********************

  double precision, intent(out) :: pre_omc(nobs_t)                ! Observed Minus Computed values to be used in LS algorithm
  double precision, intent(in)  :: gvec(maxgprm,nsat_t)           ! Observed pos, vel and accel for each satellite
  double precision, intent(in)  :: rvec(ngobs_t,nsat_t)           ! Calculated pos, vel and accel for each satellite
  double precision, intent(in)  :: gps_antoff_trf(ngobs_t,nsat_t) ! GPS antenna offset for each obs in terrestrial reference frame

  integer*4 :: isat     ! Counter for satellites
  integer*4 :: i        ! Counter variable
  !****************************************************************

  !************************* COMPUTE OMC **************************
  do isat = 1, nsat_t
     if(gvec(4,isat) == 0.d0) cycle ! Make sure given epoch existed in GNV1B file
     do i = igobs, igobs+ngobs_t-1
        pre_omc(i+(isat-1)*ngobs_t) = gvec(i,isat) - rvec(i,isat) + gps_antoff_trf(i,isat)
!                print *, gvec(i,isat), " -",  rvec(i,isat) , '=', pre_omc(i+(isat-1)*ngobs_t), '...' ,gps_antoff_trf(i,isat)
     enddo
  enddo
  !****************************************************************

  return
end subroutine gobs_computeOMC

!********************************************************************************************************************************
!********************************************************************************************************************************
! gobs_computePartials: Computes the partial derivatives of the gobs given the part matrix
!                       in which the values are stored and a vector containing pos, vel,
!                       acc and partials.
! Author: Thomas Greenspan
!********************************************************************************************************************************

subroutine gobs_computePartials(part,rvec,GPS_ind)
  use gracefit_mod

  implicit none

   

  !********************  Variable declarations ********************

  double precision, intent(out) :: part(nobs_t,nparam_t) ! Partials of observables with respect to the parameters
  double precision, intent(in)  :: rvec(nrvec_t,nsat_t)  ! Positions, velocities and partials for each satellite
  double precision, intent(in)  :: GPS_ind(nsat_t)       ! Indicator as to whether there were observations read in (0.0 if not)

  integer*4 :: isat   ! Counter for satellites
  integer*4 :: i,j    ! Counter variable
  !****************************************************************

  !*********************** COMPUTE PARTIALS ***********************
  do isat = 1, nsat_t
     if(GPS_ind(isat) == 0.d0) cycle  ! Make sure given epoch existed in GNV1B file
     do i = 1, ngobs_t
        ! Satellite specific parameters
        do j = 1, norbprm_t
           part(i+igobs-1+(isat-1)*ngobs_t,j+nsatprm_t*(isat-1)) = rvec(nrvecprm_t+(j-1)*nrvecprm_t+i,isat)
        enddo
        ! Mascon parameters
        do j = imascons, imascons+nmascons_t-1
           part(i+igobs-1+(isat-1)*ngobs_t,j) = rvec((norbprm_t+1)*nrvecprm_t+(j-imascons)*nrvecprm_t+i,isat)
        enddo
        ! PT140618: mascon tidal amplitude partials
        if(imsctide > 0 )then
           do j = imsctide, imsctide+nmsc_tid_constit_t*2-1
              ! SA/PT160719: add nmascons_t to the offset in the rvec row
              part(i+igobs-1+(isat-1)*ngobs_t,j) = rvec((norbprm_t+nmascons_t+1)*nrvecprm_t+(j-imsctide)*nrvecprm_t+i,isat)
           enddo
        endif
     enddo  ! End of observations loop
  enddo  ! End of satellite loop
  !****************************************************************


  return
end subroutine gobs_computePartials

!********************************************************************************************************************************
