!********************************************************************************************************************************
!  File: kb_lib.f90
!
!  Purpose: Set of subroutines concerning kband calculations.
!
!  Author: Thomas Greenspan
!
!  API:
!       kb_computeTheoretical      : Computes the theoretical range, range rate and range acceleration
!       kb_computeOMC              : Computes the omc
!       kb_computePartials         : Computes the partials of R, RR and RA with respect to observations
!       kb_validate                : Validates that the values of omc are within prefit tolerance
!       kb_kbrrOutliers            : Detects outliers in kbrr postfit residuals
!       kb_kbrrDecompose           : runs an empirical decomposition on the kbrr residuals to separate noise, signal and orbit error
!       kb_computeSim              : compute range rate from pos/vel from "truth" GTORB files. Also move truthvec to gvec
!       kb_sponge                  : fit a 1/rev with linearly varying amplitude to the prefit omc kbrr   PT131121
!       kb_computeKBRA             : numerically differentiate the kbrr theoretical and partials for kbra PT170428
!       kb_differentiate_kbr       : numerically differentiate kbr to generate our own kbrr and kbra      PT170831
!       kb_filtered
!       kb_filtered_v2             : handles KBRR data with no obs during shadow (but doesn't handle short gaps!) PT190502
!       kb_cos_fft                 : filters high frequency noise from kbrr/ra residuals RM190722
!       kb_computePartials_vec     : Seb's new routine to compute the kbrr partials
!       kb_computePartials_Amat    : Routine to compute the kbrr partials and store in the 2D Amat matrix
!
!   July 23, 2013
!
!********************************************************************************************************************************
!********************************************************************************************************************************
!  File: kb_decomp_v2.f90
!
!  Purpose: Makes the model for filter_prefit and sends the synthetic data back
!
!  Author:  Erin Gray
!
! Contact info:   E. Gray esgray@princeton.edu
!                 P. Tregoning paul.tregoning@anu.edu.au
!
! Input files:  plt_unfilter.kb with missing epochs (or none as well) into filter_prefit
!
!
! Output files:  sends full list of data with synthetic points back to filter_prefit
!
!********************************************************************************************************************************

SUBROUTINE kb_decomp_v2(epoch, prefit, missing, missingepochs, totalrows, goodprefit, mean)
  USE gracefit_mod

  IMPLICIT NONE




  !********************  Variable declarations ********************
  INTEGER*4           :: ioerr                               ! Error number
  INTEGER*4           :: LUGN_KB                             ! File name
  CHARACTER(64)       :: header                              ! Header Part
  INTEGER*4           :: totalrows, missingepochs            ! Rows in the File
  INTEGER*8           :: i, k, j                             ! Counter variables
  INTEGER*4           :: epoch(totalrows)                    ! Stores epoch numbers
  DOUBLE PRECISION    :: prefit(totalrows)                   ! Stores prefit residuals
  INTEGER*4    :: missing(500)              ! List that stores missing epochs
  INTEGER*4 :: tester(500)
  INTEGER*4           :: goodepochs(17270)
  DOUBLE PRECISION           :: goodprefit(17270)
  DOUBLE PRECISION :: mean
  !***************** Initialize Values ************************
  j = 1
  DO  k = 1, totalrows
     DO j = 1, missingepochs
        IF (epoch(k) == missing(j)) THEN
           tester(j) = epoch(k)
        ENDIF
     ENDDO
  ENDDO


  DO k = 1, totalrows
     IF (prefit(k) /= 0) THEN
        goodepochs(j) = epoch(k)
        goodprefit(j) = prefit(k)
        j = j +1
     ENDIF
  ENDDO

  RETURN
END SUBROUTINE kb_decomp_v2

!********************************************************************************************************************************
! kb_computeTheoretical: Computes the theoretical range, range rate and range acceleration given a vector
!                        with positions, velocities and accelerations for each satellite
!                        NOTE: it will only calculate things as far as is indicated by nrvecprm_t (the length
!                              of the vector read in)
! Author: Unknown
!         (modified by Thomas Greenspan)
!********************************************************************************************************************************

SUBROUTINE kb_computeTheoretical(pos_vel_acc,kb_theor)
  USE gracefit_mod

  IMPLICIT NONE



  !********************  Variable declarations ********************

  DOUBLE PRECISION, INTENT(out) :: kb_theor(3)                     ! Theoretical Range, Range Rate and Range Acceleration
  DOUBLE PRECISION, INTENT(in)  :: pos_vel_acc(nrvecprm_t,nsat_t)  ! Positions, velocities and accelerations for each satellite

  INTEGER*4 :: i                           ! Counter variable
  DOUBLE PRECISION :: drv,drm,dra,dvv      ! Intermediate steps in calculations
  REAL(kind=8)     :: dxdax,dterm2         ! some more intermediate steps in calcs
  DOUBLE PRECISION :: r12(3),v12(3),a12(3) ! Differences in position, velocity and acceleration in x, y and z coordinates
  !****************************************************************
  ! Set initial values
  drv = 0.d0
  dvv = 0.d0
  drm = 0.d0
  dra = 0.d0
  dxdax = 0.d0
  dterm2= 0.d0

  !****************************************************************
  ! NOTE:
  !  Equations for theoretical range rate and
  !  range acceleration found by taking equation
  !  for range and differentiating it once and
  !  twice respectively
  !
  !  We get the following equations:
  !  (For notation purposes, x represents any coord)
  !  R  = [SUM{Dx^2}]^1/2
  !  RR = [SUM{DxDvx}]/R
  ! PT170427: this range acceleration equation is wrong, I think.
  !!  RA = [SUM{DxDax + Dvx^2}]/R^2 + [RR/R]^2
  ! try
  !  RA = (vx)^2/R - dx*RR*(vx_a+vx_b)/R^2 + sum(dxdax)/R

  ! this came from some paper, but doesn't work ...
  !  RA = 1/R [ SUM[DxDax + Dvx^2] - RR^2 ]
  !****************************************************************


  ! Set intermediate values
  DO i = 1, 3
     IF (nrvecprm_t > 0) r12(i) = pos_vel_acc(i,nsat_t)   - pos_vel_acc(i,1)
     IF (nrvecprm_t > 3) v12(i) = pos_vel_acc(i+3,nsat_t) - pos_vel_acc(i+3,1)
     IF (nrvecprm_t > 6) a12(i) = pos_vel_acc(i+6,nsat_t) - pos_vel_acc(i+6,1)
     drv = drv + v12(i) * r12(i)
     dvv = dvv + v12(i) * v12(i)
     drm = drm + r12(i) * r12(i)
     dra = dra + r12(i) * a12(i) + v12(i) * v12(i)
     dxdax = dxdax + r12(i)*a12(i)
  ENDDO

  !*********************** CALCULATE VALUES ***********************

  ! Range
  kb_theor(iRANGE) = dsqrt(drm)
  ! Range Rate
  kb_theor(iRRATE) = drv / dsqrt(drm)
  ! Range Acceleration
  ! PT170427: this is not at all the equation for RA above.....
  !    kb_theor(iRACC)  = dra / dsqrt(drm) - drv**2 / dsqrt(drm)**3
  kb_theor(iRACC) = dra/kb_theor(iRANGE)**2 + (kb_theor(iRRATE)/kb_theor(iRANGE))**2
  !

  ! try this:
  !    do i=1,3
  !      dterm2 = dterm2 + r12(i)*kb_theor(iRRATE)/kb_theor(iRANGE)**2 * (pos_vel_acc(i+3,1)+pos_vel_acc(i+3,nsat_t))
  !    enddo
  !    kb_theor(iRACC) = dvv**2/kb_theor(iRANGE) - dterm2 + dxdax/kb_theor(iRANGE)

  ! this is from some paper and doesn't work
  !    kb_theor(iRACC) = 1.d0/kb_theor(iRANGE) * (dra - kb_theor(iRRATE)**2)
  !****************************************************************

  RETURN
END SUBROUTINE kb_computeTheoretical

!********************************************************************************************************************************
!********************************************************************************************************************************
! kb_computeOMC: Computes the omc (observed minus computed) given the pre_omc matrix in which the values
!                are to be stored, given a vector of position, velocity and acceleration, and the
!                observation matrices, kbrange, kblt and kbant (taken from KBR1B file)
!                NOTE: Currently using AOC instead of the kbant correction !@#
! Author: Thomas Greenspan
!********************************************************************************************************************************

SUBROUTINE kb_computeOMC(iepoch,pre_omc,pos_vel_acc,kbrange,kblt,kbant,AOC)
  USE gracefit_mod

  IMPLICIT NONE



  !********************  Variable declarations ********************

  INTEGER*4       , INTENT(in)  :: iepoch                          ! epoch number (for debugging purposes)
  DOUBLE PRECISION, INTENT(out) :: pre_omc(nobs_t)                 ! Observed Minus Computed values to be used in LS algorithm
  DOUBLE PRECISION, INTENT(in)  :: pos_vel_acc(nrvecprm_t,nsat_t)  ! Positions, velocities and accelerations for each satellite
  DOUBLE PRECISION, INTENT(in)  :: kbrange(3)                      ! Range value (biased), rate, and acceleration
  DOUBLE PRECISION, INTENT(in)  :: kblt(3)                         ! Light time corrections for range, rate and acceleration
  DOUBLE PRECISION, INTENT(in)  :: kbant(3)                        ! Antenna phase center corrections for range, rate and acceleration
  DOUBLE PRECISION, INTENT(in)  :: AOC(9,3)

  INTEGER*4 :: i                  ! Counter variable
  DOUBLE PRECISION :: kb_theor(3) ! Theoretical Range, Range Rate and Range Acceleration
  !****************************************************************

  CALL kb_computeTheoretical(pos_vel_acc,kb_theor)
  !************************* COMPUTE OMC **************************
  !   Range
  IF(nkbr_t /= 0)THEN
     IF (kb_theor(iRANGE) == 0.0) CALL status_update('FATAL','GRACEFIT','kb_computeOMC',' ',&
          'theoretical range not set',0)
     pre_omc(ikbr) = (kbrange(iRANGE)+kblt(iRANGE)+AOC(1,iRANGE)+AOC(nsat_t,iRANGE)) - kb_theor(iRANGE) !@# is this right????
  ENDIF
  !   Range rate
  IF(nkbrr_t /= 0)THEN
     IF (kb_theor(iRRATE) == 0.0) CALL status_update('FATAL','GRACEFIT','kb_computeOMC',' ',&
          'theoretical range rate not set',0)
     ! DEBUG
     ! PT140523: use the Level1B kbant correction instead of our one
! DEBUG: PT190423: test out our own computation for 2017-05-04
     !      pre_omc(ikbrr) = (kbrange(iRRATE)+kblt(iRRATE)+AOC(1,iRRATE)+AOC(nsat_t,iRRATE)) - kb_theor(iRRATE)
     pre_omc(ikbrr) = (kbrange(iRRATE)+kblt(iRRATE)+kbant(iRRATE)) - kb_theor(iRRATE)

  ENDIF
  !   Range acceleration
  IF(nkbra_t /= 0)THEN
     IF (kb_theor(iRACC) == 0.0) CALL status_update('FATAL','GRACEFIT','kb_computeOMC',' ',&
          'theoretical range acceleration not set',0)
     ! PT170427: use the Level1B kbant correction instead of our one
     !      pre_omc(ikbra) = (kbrange(iRACC)+kblt(iRACC)+AOC(1,iRACC)+AOC(nsat_t,iRACC)) - kb_theor(iRACC)
     ! PT170428: in fact, since we don't know the analytical expression for range acceleration, just store the corrected obs for now in omc
     pre_omc(ikbra) = (kbrange(iRACC)+kblt(iRACC)+kbant(iRACC)) ! - kb_theor(iRACC)
  ENDIF
  !****************************************************************

  RETURN
END SUBROUTINE kb_computeOMC

!********************************************************************************************************************************
!********************************************************************************************************************************
! kb_computePartials: Compute the partial derivatives of R, RR and RA wrt each parameter to
!                     be estimated given the part matrix in which to store the values and a
!                     vector of position, velocities, accelerations and partials.
! Author: Unknown
!         (modified by Thomas Greenspan)
!********************************************************************************************************************************

SUBROUTINE kb_computePartials(iepoch,part,rvec,kbrr_part,kbra_part)
  USE gracefit_mod
  IMPLICIT NONE



  !********************  Variable declarations ********************

  integer         , intent(in)  :: iepoch                 ! epoch number
  double precision, intent(out) :: part(nobs_t,nparam_t)  ! Partials of observables with respect to the parameters
  double precision, intent(in)  :: rvec(nrvec_t,nsat_t)   ! Positions, velocities and partials for each satellite

  integer          :: isat                       ! Counter for satellites
  integer          :: i,j                        ! Counter variables
  double precision :: r12(3),v12(3),a12(3)       ! Differences in position, velocity and acceleration in x, y and z coordinates
  double precision :: kb_theor(3)                ! Theoretical Range, Range Rate and Range Acceleration
  double precision :: temp_partials(3,6,nsat_t)  ! Partials of range, range rate and range acceleration wrt to obs
  double precision :: kbrr_part(6,2)             ! temporary range rate         partials passed back for debug
  double precision :: kbra_part(6,2)             ! temporary range acceleration partials passed back for debug

  ! DEBUG
  ! PT140519: store then output every dRR/dP * dP/dP_o partial. I can then multiply them by delta(P_o) to get a delta(RR)
  double precision :: all_kb_partials(nsatprm_t*6,2)   ! 6 per pos/vel x 6 params, 18 per bias x 3 biases, 18 per scale x 3 scales, 6 x twice-per-rev x 2 components  = 84 in total per satellite
  double precision :: all_kb_tmsc_partials(nmsc_tid_constit_t*2)  ! two per tidal constituent per mascon

  !****************************************************************

  ! Compute intermediate values
  DO i = 1, 3
     IF (nrvecprm_t > 0) r12(i) = rvec(i,nsat_t)   - rvec(i,1)
     IF (nrvecprm_t > 3) v12(i) = rvec(i+3,nsat_t) - rvec(i+3,1)
     IF (nrvecprm_t > 6) a12(i) = rvec(i+6,nsat_t) - rvec(i+6,1)
  ENDDO

  CALL kb_computeTheoretical(rvec(1:nrvecprm_t,:),kb_theor)

  !****************************************************************
  ! NOTE:
  !  paritials are found using chain rule:
  !  dR/dP  = SUM{dR/dx.dx/dP}
  !  dRR/dP = SUM{dRR/dx.dx/dP}
  !  dRA/dP = SUM{dRA/dx.dx/dP}
  !  where:
  !    x is any observation (pos, vel, acc for each coord)
  !    dx/dP = SUM{dx/dp} where p is any parameter
  !  The {dx/dP} are found in the GTORB file and stored in rvec
  !
  !  Our theoretical equations are:
  !  (For notation purposes, x represents any coord)
  !  R  = [SUM{Dx^2}]^1/2
  !  RR = [SUM{DxDvx}]/R
  !  RA = [SUM{DxDax + Dvx^2}]/R^2 + [RR/R]^2
  !
  !  The partials are calculated from these equations
  !****************************************************************

  ! Compute the Satellite IC parameter partials of R, RR, RA wrt pos, vel and acc
  DO isat = 1, nsat_t
     DO i = 1, 3  ! Loop over coordinates
        !   Position partials
        IF (nrvecprm_t > 0)THEN
           !         Partials of R
           temp_partials(iRANGE,i,isat) = (-1.d0)**(isat) * r12(i) / kb_theor(iRANGE)
           !         Partials of RR
           temp_partials(iRRATE,i,isat) = (-1.d0)**(isat) * v12(i)/kb_theor(iRANGE)  &
                - (temp_partials(iRANGE,i,isat) * kb_theor(iRRATE) ) / kb_theor(iRANGE)
           !         Partials of RA
           temp_partials(iRACC,i,isat) = ( (-1.d0)**(isat) * a12(i) - 2*kb_theor(iRRATE) * temp_partials(iRRATE,i,isat) &
                - 2*temp_partials(iRACC,i,isat) * kb_theor(iRACC)*kb_theor(iRANGE) ) / kb_theor(iRANGE)**2
        ENDIF
        !   Velocity partials
        IF (nrvecprm_t > 3)THEN
           !         Partials of R
           temp_partials(iRANGE,i+3,isat) = 0.d0
           !         Partials of RR
           temp_partials(iRRATE,i+3,isat) = (-1.d0)**(isat) * r12(i) / kb_theor(iRANGE) !! / 1.d3  !! PT140509: temporary change of units
           !         Partials of RA
           temp_partials(iRACC,i+3,isat) = ( (-1.d0)**(isat) * 2*v12(i) - temp_partials(iRRATE,i+3,isat) &
                * kb_theor(iRRATE) ) / kb_theor(iRANGE)**2
        ENDIF
        !   Acceleration partials
        IF (nrvecprm_t > 6)THEN
           !         Partials of R
           temp_partials(iRANGE,i+6,isat) = 0.d0
           !         Partials of RR
           temp_partials(iRRATE,i+6,isat) = 0.d0
           !         Partials of RA
           temp_partials(iRACC,i+6,isat) = (-1.d0)**(isat) * r12(i) / kb_theor(iRANGE)**2 !@# or / drm if problem not fixed and important
        ENDIF
     ENDDO  ! End of coord loop
  ENDDO  ! End of satellite loop
  !****************************************************************

  ! DEBUG
!!! PT140511: store off the dRR/dX, dRR/dY etc. Pass back to gracefit and use there in debug.
  DO isat = 1,2
     DO i=1,6
        kbrr_part(i,isat) = temp_partials(iRRATE,i,isat)
        kbra_part(i,isat) = temp_partials(iRACC,i,isat)
     ENDDO
  ENDDO

  ! Calculate partials of obs and place into part() matrix
  ! print*,'norbprm_t  = ',norbprm_t      !    this = 12  for estimating pos/vel/scale/bias, 16 if also estimating once- and twice-per-rev
  ! print*,'nrvecprm_t = ',nrvecprm_t     !    this =  6
  ! print*,'nsatprm_t  = ',nsatprm_t      !    this = 12  for estimating pos/vel/scale/bias, 16 if also estimating once- and twice-per-rev
  !  stop
  DO isat = 1, nsat_t
     DO i = 1, norbprm_t
        DO j = 1, nrvecprm_t
           IF(nkbr_t /= 0)  part(ikbr,i+nsatprm_t*(isat-1)) &
                = part(ikbr,i+nsatprm_t*(isat-1)) + temp_partials(iRANGE,j,isat) * rvec(i*nrvecprm_t+j,isat)
           IF(nkbrr_t /= 0) part(ikbrr,i+nsatprm_t*(isat-1)) &
                = part(ikbrr,i+nsatprm_t*(isat-1)) + temp_partials(iRRATE,j,isat) * rvec(i*nrvecprm_t+j,isat)
           IF(nkbra_t /= 0) part(ikbra,i+nsatprm_t*(isat-1)) &
                = part(ikbra,i+nsatprm_t*(isat-1)) + temp_partials(iRACC,j,isat) * rvec(i*nrvecprm_t+j,isat)

           ! DEBUG
           !  PT140520: store the multiplied partials   dRR/dP * dP/dP_o
           all_kb_partials((i-1)*nrvecprm_t+j ,isat) = temp_partials(iRRATE,j,isat) * rvec(i*nrvecprm_t+j,isat)
        ENDDO
     ENDDO
     ! Calculate partials of mascons and place into part() matrix
     DO j = 1, nrvecprm_t
        DO i = imascons, imascons+nmascons_t-1
           IF(nkbr_t /= 0)  part(ikbr,i)  = part(ikbr,i)  &
                + temp_partials(iRANGE,j,isat) * rvec((norbprm_t+1)*nrvecprm_t+(i-imascons)*nrvecprm_t+j,isat)
           IF(nkbrr_t /= 0) part(ikbrr,i) = part(ikbrr,i) &
                + temp_partials(iRRATE,j,isat) * rvec((norbprm_t+1)*nrvecprm_t+(i-imascons)*nrvecprm_t+j,isat)
           IF(nkbra_t /= 0) part(ikbra,i) = part(ikbra,i) &
                + temp_partials(iRACC,j,isat)  * rvec((norbprm_t+1)*nrvecprm_t+(i-imascons)*nrvecprm_t+j,isat)
        ENDDO
     ENDDO
     ! Calculate partials of mascon tidal amplitudes and place into part() matrix
     IF(imsctide > 0)THEN          ! it means that we want to estimate mascon tidal amplitudes
        DO j = 1, nrvecprm_t
           DO i = imsctide, imsctide+nmsc_tid_constit_t*2 - 1
              IF(nkbr_t /= 0)  part(ikbr,i)  = part(ikbr,i)  &
                   + temp_partials(iRANGE,j,isat) * rvec((norbprm_t+nmascons_t+1)*nrvecprm_t+(i-imsctide)*nrvecprm_t+j,isat)
              IF(nkbrr_t /= 0) part(ikbrr,i) = part(ikbrr,i) &
                   + temp_partials(iRRATE,j,isat) * rvec((norbprm_t+nmascons_t+1)*nrvecprm_t+(i-imsctide)*nrvecprm_t+j,isat)
              IF(nkbra_t /= 0) part(ikbra,i) = part(ikbra,i) &
                   + temp_partials(iRACC,j,isat)  * rvec((norbprm_t+nmascons_t+1)*nrvecprm_t+(i-imsctide)*nrvecprm_t+j,isat)
              ! DEBUG
              ! PT160718: store the multiplied partials for the dKBRR/d_tidal
              all_kb_tmsc_partials(i-imsctide+1) = temp_partials(iRRATE,j,isat)  &
                   * rvec((norbprm_t+nmascons_t+1)*nrvecprm_t+(i-imsctide)*nrvecprm_t+j,isat)
           ENDDO
        ENDDO
     ENDIF

  ENDDO
  ! DEBUG: do not remove !!!!
  ! PT140520: now write out all the partials that end up in the range rate obs row in the part() matrix.
  !      write(*,'(a,i8,2(72e22.13))')"dRR/dP",iepoch,(all_kb_partials(:,i),i=1,2)
  ! PT160718: write out the partials of kbrr wrt tidal amplitudes
  !       write(*,*)"dRR/dTMSC",iepoch,(all_kb_tmsc_partials(j),j=1,nmsc_tid_constit_t*2 - 1)


  !****************************************************************

  RETURN
END SUBROUTINE kb_computePartials

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE kb_computeKBRA(use_obs,use_kb_NR,ikbrr,ikbra,nepochs_t,nobs_t,nparam_t,nrvec_t,nsat_t,nrvecprm_t,pre_omc,Amat &
     ,rvec,kbrr_part,kbra_part,n_kbr_gaps, kbrange,edge1,edge2)

! subroutine to create the obs and partials for range acceleration:
! 1. regenerate the theoretical range rate
! 2. numerically differentiate it to generate the range acceleration theoretical
! 3. numerically differentiate all the kbrr partials to generate the kbra partials
!
! P. Tregoning
! 28 April 2017
!  use gracefit_mod
!
! PT180823: difference the ranges, then compute the double derivative to get the range acceleration prefit residuals
! PT201014: make use of the use_obs logical array to know whether to use observations or not.
!           replace "part" with "ATransp"
! PT220428: include the sig_process_mod 
! PT220524: increase the dimension of use_obs to include pos and vel

  use sig_process   ! provides info on number of points to set to unusable after numerical differentiation

  IMPLICIT NONE

! passed variables
  logical     ,intent(inout) :: use_obs(nepochs_t,8)              ! flag to use KBR/LRI range/range rate/range accel + POS/VEL
  real(kind=8),intent(in)    :: use_kb_NR                         ! option for how to apply the NR filter to generate prefir kbra residuals
  INTEGER*4,   INTENT(in)    :: ikbrr,ikbra                       ! pointers to range rate and range acceleration observations in pre_omc
  INTEGER*4,   INTENT(in   ) :: nepochs_t                         ! total number of epochs
  INTEGER*4,   INTENT(in   ) :: nobs_t                            ! total number of obs
  INTEGER*4,   INTENT(in   ) :: nparam_t                          ! total number of parameters
  INTEGER*4,   INTENT(in   ) :: nrvec_t                           ! total number of pos/vel
  INTEGER*4,   INTENT(in   ) :: nsat_t                            ! total number of satellites
  INTEGER*4,   INTENT(in)    :: nrvecprm_t                        ! what is this?
  INTEGER*4,   INTENT(in   ) :: n_kbr_gaps                        ! number of gaps in the kbr observations
  REAL(kind=8),INTENT(inout) :: kbrange(3,nepochs_t)              ! Level-1B range, range rate, range acceleration obs
  INTEGER*4,   INTENT(in   ) :: edge1(n_kbr_gaps),edge2(n_kbr_gaps) ! start/stop epochs of the gaps in kbr obs

  REAL(kind=8),INTENT(inout) :: pre_omc(nobs_t,nepochs_t)         ! observed minus computed array
  REAL(kind=8),INTENT(inout) :: Amat(nobs_t*nepochs_t,nparam_t)   ! matrix of partials
! PT201007: change first dimension of rvec from "nrvec_t" to "6" to match Seb's changes in gracefit
  REAL(kind=8),INTENT(in)    :: rvec(6,nsat_t,nepochs_t)    ! pos/vel information of satellites
  REAL(kind=8),INTENT(in)    :: kbrr_part(6,2,nepochs_t)          ! range rate partials
  REAL(kind=8),INTENT(out)   :: kbra_part(6,2,nepochs_t)          ! range acceleration partials

! local variables
  REAL(kind=8),ALLOCATABLE   :: kb_theor(:,:)                     ! theoretical range, range rate, range acceleration
  INTEGER*4                  :: i,j,iepoch
  character*100              :: message


! allocate arrays
  ALLOCATE(kb_theor(3,nepochs_t))


! first, we need to regenerate the theoretical range rate
  DO iepoch=1,nepochs_t
     CALL kb_computeTheoretical(rvec(1:nrvecprm_t,:,iepoch),kb_theor(:,iepoch))
  ENDDO

! now, numerically differentiate the range rate to get the range acceleration
  CALL noise_robust_deriv(kb_theor(2,:),kb_theor(3,:),5.d0,nepochs_t,5)

! now do the same for each of the range rate partials to get range acceleration partials
! all parameters (pos,vel,bias,scale, 1/rev, 2/rev, rpy, mascons, tidal amplitudes)
  DO i=1,nparam_t
   
     ! CALL noise_robust_deriv(part(ikbrr,i,:),part(ikbra,i,:),5.d0,nepochs_t,5)
     CALL noise_robust_deriv(Amat(ikbrr:nobs_t*nepochs_t:nobs_t,i),Amat(ikbra:nobs_t*nepochs_t:nobs_t,i),5.d0,nepochs_t,5)

  ENDDO
  ! set the first 5 and last 5 epochs to zero (edge effects of the nimerical differentiation
  ! PT220828: use the variable n_edge_unusable instead of hardwiring to "5"
  use_obs(1:n_edge_unusable,3) = .false.
  use_obs(nepochs_t-n_edge_unusable+1:nepochs_t, 3) = .false.

! compute the kbra OMC.
! The pre_omc for the kbra contains only the "obs" at this stage, which is the kbra + kblt + kbant for range acceleration. So we
! need to subtract from this our numerical derivative of the theoretical kbrr

  pre_omc(ikbra,:) = pre_omc(ikbra,:) - kb_theor(3,:)
! CALL status_update('STATUS','GRACEFIT','kb_computeKBRA',' ','Using Level-1B range acceleration',0)
! PT170428: alternatively, just numerically differentiate the kbrr pre_omc
  if(nint(use_kb_NR) == 1)then
    call status_update('STATUS','GRACEFIT','kb_computeKBRA',' ' &
                     ,'Numerically differentiate prefit range rate residuals to get range acceleration',0)
    call noise_robust_deriv(pre_omc(ikbrr,:),pre_omc(ikbra,:),5.d0,nepochs_t,5)

  else if (nint(use_kb_NR) == 2)then
! PT180823: difference the ranges, then numerically differentiate twice
    call status_update('STATUS','GRACEFIT','kb_computeKBRA',' ' &
                            ,'Numerically differentiate (twice) range residuals to get range acceleration residuals',0)
    call noise_robust_deriv(kbrange(1,:)-kb_theor(1,:),pre_omc(ikbrr,:),5.d0,nepochs_t,5)
    call noise_robust_deriv(pre_omc(ikbrr,:),pre_omc(ikbra,:),5.d0,nepochs_t,5)
  else if (nint(use_kb_NR) == 3)then
! PT/RMcG190920: just use the difference of the observation (which is L1B at this stage) from the theoretical (which is numerically differentiated range rate)
!                Computation is done above the if statement, so we don't need to do anything here.
    write(message,'(a)')"pre_omc for KBRA is L1B-d_kbrr/dt "
    call status_update('STATUS','GRACEFIT','kb_computeKBRA',' ',message,0)
    
  else
    write(message,'(a,i3,a)')"option for kb_NR ",nint(use_kb_NR)," not coded."
    call status_update('FATAL','GRACEFIT','kb_computeKBRA',' ',message,0)
  endif

! again, remove the first 5 and last 5 points
! PT220428: use the variable defined in sig_process_mod
  use_obs(1:n_edge_unusable,3) = .false.
  use_obs(nepochs_t-n_edge_unusable:nepochs_t,3) = .false.
  pre_omc(ikbra,1:n_edge_unusable) = 0.d0
  pre_omc(ikbra,nepochs_t-n_edge_unusable:nepochs_t) = 0.d0
  DO j = 1,n_kbr_gaps
     IF(edge1(j) > n_edge_unusable .AND. edge2(j) < nepochs_t-n_edge_unusable)THEN
       use_obs(edge1(j)-n_edge_unusable+1:edge2(j)+n_edge_unusable-1,2:3) = .false.
       pre_omc(ikbrr,edge1(j)-n_edge_unusable+1:edge2(j)+n_edge_unusable-1) = 0.d0
       pre_omc(ikbra,edge1(j)-n_edge_unusable+1:edge2(j)+n_edge_unusable-1) = 0.d0
     ENDIF
  ENDDO
  RETURN
END SUBROUTINE kb_computeKBRA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE kb_differentiate_kbr(replace,kbrange,kblt,kbant)

! subroutine to replace the kbrr and kbra obs with those derived from numerically differentiating the KBR1B range
! 1. fix cycle slips in the range measurements
! 2. numerically differentiate it to generate the range rate observation (this is accurate to 0.05 um/s RMS, with some spikes to 0.4 um/s)
! 3. numerically differentiate the newly formed range rate to create a new range acceleration (accurate to ~2 nm/s^2 RMS, with spikes to 40 nm/s^2)
!
! P. Tregoning
! 31 August 2017
!
! MODS:
! PT201014: set observations to unusable in the use_obs array if necessary
  
  USE gracefit_mod

  IMPLICIT NONE



  ! passed variables
  LOGICAL         , INTENT(in   )    :: replace                   ! flag to replace Level-1B obs or not
  double precision, INTENT(inout)    :: kbrange(3,nepochs_t)      ! Range value (biased), rate, and acceleration
  double precision, INTENT(inout)    :: kblt(3,nepochs_t)         ! Light time corrections for range, rate and acceleration
  double precision, INTENT(inout)    :: kbant(3,nepochs_t)        ! Antenna phase center corrections for range, rate and acceleration

  ! number of points in differentiator
  integer:: num_points

  ! DEBUG
  double precision, allocatable      :: debug_vector(:)
  integer  :: i

  ALLOCATE(debug_vector(nepochs_t))
  IF(replace) THEN
!!!!!!!!!!!!!!!!!!
     !!  light time  !!
!!!!!!!!!!!!!!!!!!
     ! KBR1B light time corrections for range rate and range acceleration are the same (to nm/s and pm/s^2) as our differentiated values, so don't do it.

!!!!!!!!!!!!!!!!!
     !! Ant Offset  !!
!!!!!!!!!!!!!!!!!
     num_points = 7
     ! KBR1B antenna offset corrections for range rate and range acceleration need replacing
     debug_vector = 0.d0
     CALL noise_robust_deriv(kbant(1,:),debug_vector,5.d0,nepochs_t,num_points)  ! use a 5-point numerical derivative
     kbant(2,4:nepochs_t-3) = debug_vector(4:nepochs_t-3)
     debug_vector = 0.d0
     CALL noise_robust_deriv(kbant(2,:),debug_vector,5.d0,nepochs_t,num_points)  ! use a 5-point numerical derivative
     kbant(3,4:nepochs_t-3) = debug_vector(4:nepochs_t-3)  ! SA correct copy into acceleration kbant(3 ... not velocity kbant(2
!!!!!!!!!!!!!!!!!
     !!  R A N G E  !!
!!!!!!!!!!!!!!!!!
     ! first, identify and repair any cycle slips in the range measurements
     CALL kb_cycle_slips(kbrange,nepochs_t)

     !! next, numerically differentiate
     debug_vector = 0.d0
     CALL noise_robust_deriv(kbrange(1,:),debug_vector,5.d0,nepochs_t,num_points)  ! use a 5-point numerical derivative
     kbrange(2,4:nepochs_t-3) = debug_vector(4:nepochs_t-3)
     debug_vector = 0.d0
     CALL noise_robust_deriv(kbrange(2,:),debug_vector,5.d0,nepochs_t,num_points)  ! use a 5-point numerical derivative
     kbrange(3,4:nepochs_t-3) = debug_vector(4:nepochs_t-3)

     ! PT170901: check the kbr obs and eliminate points (ie set the range obs to zero) for any range rate/range accel. obs within num_points/2.0 of
     !           a data gap. Such points will not have accurate time differentiated values.
     !!!!SACOmment CALL kb_data_gaps(kbrange,num_points,nepochs_t)

     CALL status_update('STATUS','GRACEFIT','kb_differentiate_kbr',' '," Have replaced Level-1B KBRR and KBRA",0)
  ELSE
     CALL status_update('STATUS','GRACEFIT','kb_differentiate_kbr',' '," Did NOT replace Level-1B KBRR and KBRA",0)
  ENDIF

  DEALLOCATE(debug_vector)
  RETURN

END SUBROUTINE kb_differentiate_kbr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!********************************************************************************************************************************
!********************************************************************************************************************************
! kb_cycle_slips: Checks Leve-1B KBR values and fixes obvious step functions in range. Uses the Level-1B range rate to extrapolate
! Author: Paul Tregoning
! Date  : 4 September 2017
!********************************************************************************************************************************
SUBROUTINE kb_cycle_slips(kbrange,nepochs_t)

  IMPLICIT NONE

  ! passed variables
  INTEGER*4,    INTENT(in )   :: nepochs_t              ! number of epochs
  REAL(kind=8), INTENT(inout) :: kbrange(3,nepochs_t)   ! KBR1B observations

  ! local variables
  INTEGER*4    :: iobs,last_obs_epoch
  REAL(kind=8) :: kbr_predicted, d_range, dt
  CHARACTER    :: message*250

  ! loop through observations
  last_obs_epoch = 1
  DO iobs = 2,nepochs_t
     IF(kbrange(1,iobs) /= 0.d0)THEN    ! there is a K-band observation
        dt = (iobs - last_obs_epoch)*5.d0
        kbr_predicted = kbrange(1,last_obs_epoch) + kbrange(2,last_obs_epoch)*dt + 0.5d0*kbrange(3,last_obs_epoch)*dt**2 ! s + ut + 1/2 at^2
        d_range = kbrange(1,iobs)-kbr_predicted
        IF(dabs(d_range) .GT. 0.1d0)THEN   ! discrepancy of more than 100 mm in observed and predicted range
           WRITE(message,'(f15.3,a,i10,2f15.3)')d_range," step found in range at epoch ",iobs,kbrange(1,iobs),kbr_predicted
           CALL status_update('STATUS','GRACEFIT','kb_cycle_slips',' ',message,0)
           kbrange(1,iobs:nepochs_t) = kbrange(1,iobs:nepochs_t) - d_range

        ENDIF
        !      if(iobs > 3500)print*,'i,kbrange(1,iobs)',iobs,kbrange(1,iobs)

        ! set this epoch as the most recent one to have had observations
        last_obs_epoch = iobs
     ENDIF

  ENDDO

END SUBROUTINE kb_cycle_slips
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!********************************************************************************************************************************
!********************************************************************************************************************************
! kb_data_gaps: Finds gaps in the k-band range observations and sets to zero the range obs that are within n_points/2
!               of a data gap. Such points will not have accurate time differentiated range rate and range accelerations.
! Author: Paul Tregoning
! Date  : 1 September 2017
!********************************************************************************************************************************
SUBROUTINE kb_data_gaps(kbrange,n_points,nepochs_t)

  IMPLICIT NONE

  ! passed variables
  INTEGER*4,    INTENT(in )   :: nepochs_t              ! number of epochs
  INTEGER*4,    INTENT(in )   :: n_points               ! number of points used in the time differentiation
  REAL(kind=8), INTENT(inout) :: kbrange(3,nepochs_t)   ! KBR1B observations

  ! local variables
  INTEGER*4   :: iobs,igap
  INTEGER*4   :: start_pt,end_pt,n_pts_to_remove,n_gaps
  INTEGER*4   :: gap_epochs(1000,2)                           ! array of start/end epochs of data gaps
  LOGICAL     :: in_gap

  ! loop through observations
  in_gap = .FALSE.
  DO iobs = 1,nepochs_t
     IF(kbrange(1,iobs) == 0.d0)THEN    ! no K-band observations for this epoch
        IF(iobs > 1 ) THEN
           IF(kbrange(1,iobs-1) /= 0.d0) THEN   ! it is the start of a gap
              n_gaps = n_gaps + 1
              gap_epochs(n_gaps,1) = iobs        ! store the epoch of the start of the gap
           ENDIF
        ELSE
           PRINT*,'DEBUG: data gap starts at epoch ',iobs
           n_gaps = 1
           gap_epochs(n_gaps,1) = 1
        ENDIF
        in_gap = .TRUE.
     ELSE                               ! have k-band observations for this epoch
        ! check whether it is the first epoch after a gap
        IF(in_gap)THEN
           gap_epochs(n_gaps,2) = iobs - 1
           in_gap = .FALSE.
        ENDIF
     ENDIF
  ENDDO

  ! now, remove data points either side of gaps and start/end of orbits because the time derivatives won't be good there ...
  n_pts_to_remove = INT(DBLE(n_points)/2.d0) + 1
!  PRINT*,'number of points to remove:',n_pts_to_remove
  DO igap = 1,n_gaps
     start_pt = gap_epochs(igap,1) - n_pts_to_remove
     end_pt   = gap_epochs(igap,2) + n_pts_to_remove
     !   check that neither goes before or after the first/last point
     IF(start_pt < 1)start_pt = 1
     IF(end_pt > nepochs_t)end_pt = nepochs_t

     !   set the range obs to zero (used as a flag that there are no observations at all)
     kbrange(1,start_pt:end_pt) = 0.d0
  ENDDO
!  print *, '** n_pts_to_remove', n_pts_to_remove, 'nepochs_t', nepochs_t
  ! finally, remove the points at the start and end of the orbit
  !FIXME: BUG1: Followed if the extrem values are at 0 the EMD isn't working.
  kbrange(1,1:n_pts_to_remove) = 0.d0
  kbrange(1,nepochs_t-n_pts_to_remove+1:nepochs_t) = 0.d0

  !   stop 'stopped in kb_data_gaps'

END SUBROUTINE kb_data_gaps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!********************************************************************************************************************************
!********************************************************************************************************************************
! kb_validate: Checks to make sure that the omc values calculated are within the prefit tolerance
!              read in from the command file
! Author: Thomas Greenspan
!********************************************************************************************************************************

SUBROUTINE kb_validateRR(part_ikbrr,pre_omc_ikbrr,kbrr_misfit_num,kbrr_prefit_tol,iepoch)
  USE gracefit_mod

  IMPLICIT NONE



  !********************  Variable declarations ********************
  DOUBLE PRECISION, INTENT(inout) :: part_ikbrr(nparam_t) ! Partials for range rate
  DOUBLE PRECISION, INTENT(inout) :: pre_omc_ikbrr        ! Observed Minus Computed values for range rate
  INTEGER*4, INTENT(inout)        :: kbrr_misfit_num      ! Number of kbrr misfits up to point when subroutine is called
  DOUBLE PRECISION, INTENT(in)    :: kbrr_prefit_tol      ! KBRR prefit tolerance
  INTEGER*4, INTENT(in)           :: iepoch               ! Place in epoch loop in which subroutine is called

  CHARACTER(128) :: message    ! Message in case of bad observations
  !****************************************************************

  ! Check to see if the prefit kbrr residule exceeds the assigned tolerance
  IF( dabs(pre_omc_ikbrr) > kbrr_prefit_tol ) THEN
     WRITE(message,'(a,i5,a,f12.1,a,f6.1,a)')'Bad KBRR observation at iepoch',iepoch,' (epoch ', &
          starting_epoch+(iepoch-1)*epoch_interval,'): ',pre_omc_ikbrr*m_um,'(um) misfit'
     CALL status_update('WARNING','GRACEFIT','kb_validateRR',' ',message,0)
     pre_omc_ikbrr = 0.d0
     part_ikbrr = 0.d0
     kbrr_misfit_num = kbrr_misfit_num + 1
  ENDIF
  !****************************************************************

  RETURN
END SUBROUTINE kb_validateRR

!********************************************************************************************************************************
!********************************************************************************************************************************
! kb_kbrrOutliers: Prints out the kbrr prefit outliers using residuals smoothed by a 20 point running mean
!                  and indicates whether the ikbrr postfit residual at a certain epoch is an outlier
! Author: Paul Tregoning
!********************************************************************************************************************************

SUBROUTINE kb_kbrrOutliers(post_omc,kbrr_prefit_tol,kbrr_outliers)
  USE gracefit_mod

  IMPLICIT NONE



  !********************  Variable declarations ********************
  DOUBLE PRECISION, INTENT(in) :: post_omc(nobs_t,nepochs_t) ! Postfit Observed Minus Computed values
  DOUBLE PRECISION, INTENT(in) :: kbrr_prefit_tol            ! KBRR prefit tolerance
  LOGICAL, INTENT(out)         :: kbrr_outliers(nepochs_t)   ! Indicates whether the kbrr postfit values at a given epoch are outliers

  INTEGER*4        :: iepoch                      ! Counter variable that runs through epochs
  DOUBLE PRECISION :: kb_omc_smoothed(nepochs_t)  ! kbrr residuals smoothed by a 20 point running mean
  CHARACTER(128)   :: message                     ! Message in case of bad observations
  !****************************************************************

  ! Initialize array
  kbrr_outliers = .FALSE.

  ! Outlier detection on the kbrr residuals. For this, we generate a 20-point smoothing and then look for outliers > 3 um/s away from the smoothed value.
  CALL running_mean(nepochs_t,post_omc(ikbrr,:),20,kb_omc_smoothed)

  ! ok, so we step through all the epochs and, if the difference between kbrr residual and smoothed
  ! is > kbrr_prefit_tol
  DO iepoch = 1, nepochs_t
     IF(dabs(post_omc(ikbrr,iepoch) - kb_omc_smoothed(iepoch)) > kbrr_prefit_tol)THEN
        WRITE(message,'(a,i6,f12.1,f10.3,a)')'KBRR postfit residual outlier: ',iepoch,starting_epoch+(iepoch-1)*epoch_interval, &
             (kb_omc_smoothed(iepoch)-post_omc(ikbrr,iepoch))*m_um,' um/s'
        CALL status_update('STATUS','GRACEFIT','gracefit/kb_kbrrOutliers',' ',message,0)
        kbrr_outliers(iepoch) = .TRUE.
     ENDIF
  ENDDO
  !****************************************************************

  RETURN
END SUBROUTINE kb_kbrrOutliers

!******************************************************************************************************
! kb_KBRRdecompose: empirical decomposition of kbrr residuals to separate noise, signal and orbital error
!
! Author: Paul Tregoning  12 August 2013
!
! MODS
! PT40529: compress the pre_omc vector to remove data gaps before filtering, then re-expand afterwards
!
!******************************************************************************************************
SUBROUTINE kb_kbrrDecompose(ikbrr,kbRANGE_obs,nobs_t,nepochs_t,obs_omc,filt_omc)

  IMPLICIT NONE

  ! ********************  Variable declarations ********************
  INTEGER*4       , INTENT(in)  :: nobs_t,nepochs_t            ! number of rows,columns of the passed variables
  INTEGER*4       , INTENT(in)  :: ikbrr                      ! row number of the kbrr observation
  DOUBLE PRECISION, INTENT(in)  :: kbRANGE_obs(nepochs_t)      ! actual k-band range observations (zero if no observation)
  DOUBLE PRECISION, INTENT(in)  :: obs_omc(nepochs_t)          ! Prefit Observed Minus Computed values (range, range rate or range acceleration)
  DOUBLE PRECISION, INTENT(out) :: filt_omc(nepochs_t)         ! reconstructed using the empirical decomposition (Montillet et al., 2013)

  DOUBLE PRECISION :: prefit_omc(0:nepochs_t)                 ! KBRR omc (reformatted to a vector)
  INTEGER*4        :: ncomp , PLT1, PLT2                                   ! number of components for the decomposition
  PARAMETER (ncomp = 13)
  DOUBLE PRECISION :: decomp(ncomp,nepochs_t)                      ! decomposed values
  INTEGER*4    ::  iepoch                                      ! Counter variable that runs through epochs
  INTEGER*4    ::  i,j                                         ! do loop counters
  DOUBLE PRECISION :: mean                                     ! mean value of the input time series

  ! variables to enable the filter to handle data gaps
  INTEGER*4               :: numobs                            ! counter for the number of non-zero kbrr observations
  INTEGER*4, ALLOCATABLE  :: epoch_pntr(:)                     ! array to indicate which epoch is stored where when the missing-obs are removed

  ALLOCATE(epoch_pntr(nepochs_t))

  ! ****************************************************************
  ! This routine will decompose the input kbrr residuals into 8 components. We think they are as follows:
  ! 1. noise
  ! 2. mostly noise (but < 0.2 um/s). Perhaps a tiny bit of signal at the poles.
  ! 3. signal (< 0.3 um/s)
  ! 4. signal
  ! 5. signal
  ! 6. signal - in particular, the apparent "reverberation" of large-scale mass changes into the kbrr residuals
  ! 7. orbit error
  ! 8 and above. orbit error

  CALL status_update('STATUS','GRACEFIT','gracefit','kb_kbrrDecompose ','Decompose kbrr prefit residuals',0)

  filt_omc = 0.d0
  prefit_omc = 0.d0

  ! PT140529: the EMD filter doesn't work if there are data gaps. To overcome this, compress the data (ie simply remove the
  !           epochs without observations from the vector), filter and then expand it out. We just need to keep a vector of
  !           pointers so that we can re-establish the correct epochs
  !
  ! PT170427: removed the pointer to iKBRR, so that this routine can work for range, range rate or range acceleration
  numobs = 0
  epoch_pntr = 0
  DO i=1,nepochs_t
     IF (kbRANGE_obs(i) /= 0.d0)THEN
        numobs = numobs + 1
        prefit_omc(numobs-1) = obs_omc(i)
        epoch_pntr(i) = numobs-1
     ENDIF
  ENDDO


  ! reformat the input omc to be a column vector
  mean = SUM(prefit_omc(0:numobs))/DBLE(numobs)

  DO i=1,numobs
     prefit_omc(i-1) = prefit_omc(i-1) - mean
  ENDDO

  ! call Jean-Phillipe Montillet's subroutine to perform the decomposition
  CALL emd(prefit_omc,numobs,ncomp,decomp)

  ! now, add together the components that we want to include in the new "filtered range rate omc". Hardwire to be 2-6 for now
  DO i=1,nepochs_t
     ! Pt170906: include component 2 if using numerically differentiated range rate and range acceleration
     DO j=2,ncomp
        IF(epoch_pntr(i) /= 0)THEN
           filt_omc(i) = filt_omc(i) + decomp(j,epoch_pntr(i))
        ENDIF
     ENDDO
     IF(epoch_pntr(i) /= 0)filt_omc(i) = filt_omc(i) + mean
  ENDDO
!  DO 34 i = 1, nepochs_t
!      print*, filt_omc(i), prefit_omc(i), i
! 34   CONTINUE


     RETURN
   END SUBROUTINE kb_kbrrDecompose


   !********************************************************************************************************************************

   !********************************************************************************************************************************
   ! kb_computeSim: compute range rate from pos/vec of "truth" GTORB files for program gracesim
   !
   ! Author: Paul Tregoning  15 October 2013
   !******************************************************************************************************
   SUBROUTINE kb_computeSim(truthvec,gvec,kbrange)
     USE gracefit_mod

     IMPLICIT NONE



     ! ********************  Variable declarations ********************

! PT201007: change truthvec from (nparam_t,nsat_t,nepochs_t) to (6,nsat_t,nepochs_t)
     DOUBLE PRECISION, INTENT(in)    :: truthvec(6,nsat_t,nepochs_t)          ! positions and velocities (and partials) of truth GTORB
     DOUBLE PRECISION, INTENT(inout) :: gvec(maxgprm,nsat_t,nepochs_t)        ! array to store the GNV1B positions/velocities
     DOUBLE PRECISION, INTENT(inout) :: kbrange(3,nepochs_t)                  ! array to store the KBR1B range, range rate

     INTEGER*4 iepoch, i
     ! ****************************************************************

     ! loop through the epochs
     DO iepoch = 1,nepochs_t

        ! compute the range rate
        CALL kb_computeTheoretical(truthvec(1:nrvecprm_t,:,iepoch),kbrange(:,iepoch) )

        ! transfer from truthvec to gvec
        DO i=1,nsat_t
           gvec(1:6, i, iepoch) = truthvec(1:6, i, iepoch)
        ENDDO

     ENDDO

     RETURN

   END SUBROUTINE kb_computeSim


   !********************************************************************************************************************************


   !********************************************************************************************************************************
   ! kb_sponge: derive a model of 1/rev with quadratically varying amplitude as a fit to the kbrr prefit omc
   !
   ! Author: Paul Tregoning  21 November 2013 !PT170703 Runs data w/ gaps and makes accurate (Erin Gray)
   !******************************************************************************************************
   SUBROUTINE kb_sponge(nvals,values,model,model_params)

     IMPLICIT NONE

     INTEGER*4    , intent(in)  :: nvals                       ! number of input values
     REAL(kind=8 ), intent(in)  :: values(nvals)               ! input time series
     REAL(kind=8 ), intent(out) :: model(nvals)                ! modelled values
! PT190807: model equation is (a*dt+b+c*dt^2)cos(w*dt+phi)+offset+rate*dt
     real(kind=8 ), intent(out) :: model_params(7)             ! model parameter values (params are a,b,c,phi,offset,rate). omega passed back but fixed.
     INTEGER*4     :: niter,iter
     INTEGER*4  :: i,j,nparam

     ! variables for the equation to model the unwanted signal in the kbrr prefit residuals
     REAL(kind=8)  :: a,b,c,w,phi,dt,pi,offset, rate


     ! variables for LS
     REAL(kind=8),ALLOCATABLE   :: amat(:,:), omc(:,:), AtA(:,:), Atb(:,:), At(:,:), VCV(:,:),soln(:,:)
     CHARACTER message*250

     pi = 4.d0*datan(1.d0)

!!! for the formulation of a phase and offset/rate of amplitude, period, nparam = 4. Adding offset and rate makes it 6.
     ! PT140203: try adding a quadratic term to the amplitude of the 1/rev signal. This makes nparam = 7
     nparam = 7
     niter = 70

     ! allocate the variables
     ALLOCATE(amat(nvals,nparam))
     ALLOCATE(omc(nvals,1))
     ALLOCATE(AtA(nparam,nparam))
     ALLOCATE(Atb(nparam,1))
     ALLOCATE(At(nparam,nvals))
     ALLOCATE(VCV(nparam,nparam))
     ALLOCATE(soln(nparam,1))

     !  call status_update('STATUS','GRACEFIT','gracefit/sim','kb_kbsponge ','kband sponge is now a trend + quadratically varying 1/rev amplitude',0)
     CALL status_update('STATUS','GRACEFIT','kb_sponge',' ','kband sponge is now a trend + quadratically varying 1/rev amplitude',0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! assign a priori values to parameters
     a = 1.d-9
     b = 0.d0
     c = 1.d-9
     phi = 0.d0
     offset = 0.d0
     rate = 0.d0
     ! assign the period
     w = 2.d0*pi/1120.d0     ! 1860 seconds in 93 minutes which is the time of one revolution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! we  need to detrend the time series so that we can estimate the 1/rev model, so add an offset and slope to the model. The values
     ! are not passed back .... we want to leave the rate in the kbrr prefit because it comes from the accelerometer biases !

     DO iter = 1,niter
        ! form up the A matrix. columns are for parameters a, b, phi, w
        ! Y(t) = offset + rate*(t-t0) + (a(t-t0)+b) * cos(w(t-t0) + phi)    where w is the 1/rev period
        ! PT140203: try adding a quadratic term to the amplitude of the 1/rev:
        ! Y(t) = offset + rate*(t-t0) + (a(t-t0)+b + c(t-t0)^2) * cos(w(t-t0) + phi)    where w is the 1/rev period

        !if(iter>60)WRITE(*,'(a,i3,7e18.9)')'iter,a,b,c,w,phi,offset,rate',iter,a,b,c,w,phi,offset,rate
        DO i=1,nvals
           !if(values(i) /= 0) THEN
           dt = DBLE(i-1)
           amat(i,1) = dt * dcos(w*dt + phi)                   ! amplitude of rate of change of 1/rev amplitude
           amat(i,2) = dcos(w*dt + phi)                        ! constant part of amplitude of 1/rev signal
           amat(i,3) = -dsin(w*dt + phi) * (a*dt+b+c*dt**2)    ! phase partial
           amat(i,4) = -dsin(w*dt + phi) * dt*(a*dt+b+c*dt**2) ! w partial
           amat(i,5) = 1.d0                                    ! offset partial
           amat(i,6) = dt                                      ! rate partial
           amat(i,7) = dt**2 * dcos(w*dt + phi)                ! amplitude of quadratic of change of 1/rev amplitude

           omc(i,1) = values(i) - ( (a*dt+b+c*dt**2) * dcos(w*dt + phi) + offset + rate*dt)
           !endif
        ENDDO

        ! transpose A
        CALL transp(Amat,At,nvals,nparam)
        ! AtA
        CALL matmult(At,Amat,AtA,nparam,nvals,nparam)

        ! PT131206: add a really tight constraint on the rate of the amplitude - to make it just a simple 1/rev model
        !    AtA(1,1) = AtA(1,1) + 1.d12

        ! invert it
        CALL invert(AtA,VCV,nparam)

        ! Atb
        CALL matmult(At,omc,Atb,nparam,nvals,1)

        ! AtA^-1Atb
        CALL matmult(VCV,Atb,soln,nparam,nparam,1)

!!!!!!!!!!!!!!!!!!! update a priori parameter values !!!!!!!!!!!!!
        a =  a + soln(1,1)
        b = b + soln(2,1)
        phi = phi + soln(3,1)
        w = w + soln(4,1)
        offset = offset + soln(5,1)
        rate = rate + soln(6,1)
        c = c + soln(7,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ENDDO

!!!!!!!!!!!!!!!!!!! output the final solution for the parameters !!!
     WRITE(message,'(a,7e18.9)')"1/rev model A,b,c,omega,phi,offset,rate:",a,b,c,w,phi,offset,rate
     CALL status_update('STATUS','GRACEFIT','kb_sponge',' ',message,0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!! calculate model values !!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! PT131206: only return the 1/rev part of this - leave the offset and rate out of it.
     DO i=1,nvals
        dt = DBLE(i-1)
        model(i) = (a*dt+b+c*dt**2) * dcos(w*dt + phi)  
     ENDDO

! PT190807: return the model parameters in a single array
     model_params(1) = a
     model_params(2) = b
     model_params(3) = c
     model_params(4) = w
     model_params(5) = phi
     model_params(6) = offset
     model_params(7) = rate
   

   END SUBROUTINE kb_sponge

  subroutine kb_filtered_v2(calling_prog,obs_type,nepochs_t,start_neq,end_neq,start_reconstruct,prefit,filtered)

! subroutine designed to apply the EMD filtering to kbrr observations where no data are available when the satellites are in eclipse.
! This occurs in the final year of the GRACE mission and affects the EMD filter badly. Therefore, instead of running the EMD across
! the whole of the data span (including infilling gaps beforehand) I am going to try running it on each segment of kbr data (i.e. each time
! when the satellites come out of eclipse to when they go back into eclipse)
! 
! P. Tregoning
! 1 May 2019

  implicit none

!********************  Variable declarations ********************
  character*10,intent(in) :: calling_prog                           ! name of calling program
  character*4, intent(in) :: obs_type                               ! 'kbrr' or 'kbra'
  integer*4,intent(in):: nepochs_t                                  ! total number of epochs
  INTEGER*4,INTENT(in):: start_neq,end_neq                          ! start and end epochs for which we are going to stack normal equations
  INTEGER*4,INTENT(in):: start_reconstruct                          ! reconstruct EMD recomposition from this component
  real(kind=8),intent(inout) :: prefit(nepochs_t)                   ! Prefit OMC (range, range rate, or range acceleration)
  real(kind=8),intent(inout) :: filtered(nepochs_t)                 ! reconstructed prefit OMC using the empirical decomposition (Montillet et al., 2013)

! local variables
  integer*4      :: segment_start,segment_end
  integer*4      :: iepoch,i,j,nvals
  integer*4      :: gap_epochs, max_gap
  character*256  :: message
  real(kind=8)   :: mean
  integer*4      :: ncomp                                           ! number of components for the decomposition
  parameter (ncomp = 13)                                            ! Given component Value
  real(kind=8),allocatable :: decomp(:,:)                           ! decomposed values
!*****************************************************************


!*****************************************************************
! allocate and assign values
  max_gap = 10     ! a gap of 60 epochs (5 minutes) will cause the creation of a new segment in the KBR data
!*****************************************************************


! set the segment start/end to zero
  segment_start = 0
  segment_end = 0
  gap_epochs = 0

! loop through all the required epochs
  do iepoch = start_neq,end_neq
    if(prefit(iepoch) /= 0.d0)then   ! we have an observation
      if(segment_start == 0)then     ! this is the beginning of a new segment of observations
        segment_start = iepoch
        segment_end = iepoch

      else          
        segment_end = iepoch
        ! PT190531: if it is the last epoch then we need to filter the segment
        if(iepoch == end_neq)then
          write(message,'(a,2(a,i7))')obs_type,': Need segment from ',segment_start,' to ',segment_end
          call status_update('STATUS',calling_prog,'kb_filtered_v2',' ',message,0)
          mean = SUM(prefit(segment_start:segment_end))/DBLE(segment_end-segment_start+1)
          do i= segment_start,segment_end
            prefit(i) = prefit(i) - mean
          enddo
          nvals = segment_end-segment_start+1
          allocate(decomp(ncomp,nvals))
          decomp = 0.d0
          call emd(prefit(segment_start:segment_end),nvals,ncomp,decomp)                             ! CALL THE FILTER
          ! PT170728: put the mean back into prefit
          do i = segment_start,segment_end
            prefit(i) = prefit(i) + mean
            ! PT190501: reconstruct the filtered signal
            if(obs_type == "kbrr" .or. obs_type == "kbra")then
              do j=start_reconstruct,ncomp
                filtered(i) = filtered(i) + decomp(j,i-segment_start+1)         ! MOVE FILTERED DATA OVER
              enddo
              filtered(i) = filtered(i) + mean
            endif
          enddo

          deallocate(decomp)
          ! reset so that it can search for the next segment
          segment_start = 0
          segment_end = 0
          gap_epochs = 0
        endif


      endif

    else    ! we don't have an observation
      if(segment_start /= 0)then     ! this is the start of a gap, or just a missing observation
        gap_epochs = gap_epochs + 1

      ! is the gap big enough to mean that we need a new segment
        if(gap_epochs > max_gap .or. iepoch == end_neq)then   ! yes, we need a new segment
          ! PT190502: for kbra, remove the first and last 5 points
          if(obs_type == "kbra")then
          segment_start = segment_start + 5
          segment_end = segment_end - 5
          endif

          write(message,'(a,2(a,i7))')obs_type,': Need segment from ',segment_start,' to ',segment_end
          call status_update('STATUS',calling_prog,'kb_filtered_v2',' ',message,0)


      ! send just this section of data to EMD to be filtered
        mean = SUM(prefit(segment_start:segment_end))/DBLE(segment_end-segment_start+1)
! debug
!if(obs_type == "kbra")then
!  print*,obs_type," segment: ",segment_start,segment_end," mean",mean
!endif

        do i= segment_start,segment_end
          prefit(i) = prefit(i) - mean
        enddo

!        call status_update('STATUS',calling_prog,'kb_filtered_v2',' ','Filtering Data',0)
        nvals = segment_end-segment_start+1
        allocate(decomp(ncomp,nvals))

        decomp = 0.d0
        call emd(prefit(segment_start:segment_end),nvals,ncomp,decomp)                             ! CALL THE FILTER
        ! PT170728: put the mean back into prefit
        do i = segment_start,segment_end
          prefit(i) = prefit(i) + mean
          ! PT190501: reconstruct the filtered signal
          if(obs_type == "kbrr" .or. obs_type == "kbra")then
            do j=start_reconstruct,ncomp
              filtered(i) = filtered(i) + decomp(j,i-segment_start+1)         ! MOVE FILTERED DATA OVER
            enddo
            filtered(i) = filtered(i) + mean
! DEBUG
!if(obs_type == "kbra")print*,obs_type,i,prefit(i),filtered(i),mean,decomp(1:4,i-segment_start+1) 
          else
          endif

        enddo


        deallocate(decomp)
          ! reset so that it can search for the next segment
          segment_start = 0
          segment_end = 0
          gap_epochs = 0
        !stop ' stopped 1'


        endif

      endif

    endif

  enddo



  end subroutine kb_filtered_v2


   !********************************************************************************************************************************
   !
   !  kb_filtered: Model missing epochs (if there are any), synthesize missing data,filter it, remove synthesized data
   !
   !  Author:  Erin Gray
   !
   ! Contact info:   E. Gray esgray@princeton.edu
   !                 P. Tregoning paul.tregoning@anu.edu.au
   !
   ! Output files:  Correctly Filtered Data (filtered)
   !
   !********************************************************************************************************************************

   SUBROUTINE kb_filtered(obs_type,nepochs_t,start_neq,end_neq,start_reconstruct,prefit,filtered)

     IMPLICIT NONE


     !********************  Variable declarations ********************
     character*4, intent(in) :: obs_type                               ! 'kbrr' or 'kbra'
     integer*4,intent(in):: nepochs_t                                  ! total number of epochs
     INTEGER*4,INTENT(in):: start_neq,end_neq                          ! start and end epochs for which we are going to stack normal equations
     INTEGER*4,INTENT(in):: start_reconstruct                          ! reconstruct EMD recomposition from this component
     real(kind=8),intent(inout) :: prefit(nepochs_t)                   ! Prefit OMC (range, range rate, or range acceleration)
     real(kind=8),intent(inout) :: filtered(nepochs_t)                 ! reconstructed prefit OMC using the empirical decomposition (Montillet et al., 2013)


! local variables
     INTEGER*4             :: i, k, j                                  ! Counter variables
     INTEGER*4             :: beginning, ending, nvals                 ! Modelling Data values
     INTEGER*4,allocatable :: epoch(:)                                 ! Stores epoch numbers
     real(kind=8)          :: mean                                     ! mean value of prefit residuals (need to remove the mean to use EMD)
     INTEGER*4             :: missingepochs                            ! total columns in file that are missing
     INTEGER*4,allocatable :: missing(:), ordered(:)                   ! Stores Missing Epochs, Missing Order of Epochs
     INTEGER*4             :: ncomp                                      ! number of components for the decomposition
     PARAMETER (ncomp = 13)                                            ! Given component Value
     real(kind=8),allocatable :: decomp(:,:)                                ! decomposed values
     INTEGER*4             :: ioerr                                      ! Error number
     INTEGER*4             :: firstm, lastm, check1, check2, PLT_FIXY, ngaps              ! For Multiple Gapsw
     integer*4             :: edge
     real(kind=8),allocatable :: original(:)
     character*150   ::  message

! allocate variables
     allocate(epoch(nepochs_t))
     allocate(original(nepochs_t))
     allocate(missing(nepochs_t))
     allocate(ordered(nepochs_t))
     allocate(decomp(ncomp,nepochs_t))
    
     edge = 4
     !****************************** INITIALIZE VALUES **************************************
     missingepochs = 0
     ordered = 0
     filtered = 0.d0
     check1 = -1
     check2 = -2
     ngaps = 0
     missing = 0
     ordered = 0
     decomp = 0
     filtered = 0
     epoch = 0
     !****************************** CHECK FOR MISSING VALUES *******************************
     DO i = start_neq, end_neq
        epoch(i) = i
!        original(i)=prefit(i)
     ENDDO
     original = prefit


     IF (start_neq == 1) THEN
        j = 1                                                                  !CORNER CASES FOR FIRST MISSING EPOCHS

        DO k = 1, 9
           ordered(k) = epoch(k)
        ENDDO

        DO k = 10, end_neq-10
           IF ((prefit(k)) == 0) THEN     ! there are no kbr obs at this epoch

              IF(prefit(k-1) /= 0) THEN    ! epoch k is the first epoch of a gap in the kbr data
                 missing(j) = -1
                 j = j + 1
              ENDIF

              missing(j) = epoch(k)
              j = j + 1

              IF(prefit(k+1) /= 0) THEN    ! epoch K is the last epoch of a data gap in kbr obs
                 missing(j) = -2
                 j = j + 1
              ENDIF

           ELSE
              ordered(k) = epoch(k)                                              !MAKE INVERSE LIST OF MISSING

           ENDIF
        ENDDO

        DO k = end_neq -9, end_neq                                      !CORNER CASES FOR LAST MISSING EPOCHS
           ordered(k) = epoch(k)
        ENDDO

     ELSE

        j = 1

        DO k = start_neq, end_neq - 10
           IF (prefit(k)== 0) THEN
              IF((prefit(k-1)) /= 0) THEN
                 missing(j) = -1
                 j = j + 1
              ENDIF                                           !ADD MISSING EPOCHS TO A LIST FOR STORAGE

              missing(j) = epoch(k)
              j = j + 1

              IF(prefit(k+1) /= 0) THEN       ! epoch K is the last epoch of a data gap in kbr obs
                 missing(j) = -2
                 j = j + 1
              ENDIF

           ELSE
              ordered(k) = epoch(k)                                              !MAKE INVERSE LIST OF MISSING

           ENDIF
        ENDDO

        DO k = end_neq -9, end_neq                                      !CORNER CASES FOR LAST MISSING EPOCHS
           ordered(k) = epoch(k)
        ENDDO

     ENDIF

     i = 0
     CALL status_update('STATUS','GRACEFIT','kb_filtered',' ','STORED',ioerr)

     DO j = 1, end_neq
        IF (missing(j) > 0) THEN
           i = i + 1
        ENDIF
     ENDDO


     missingepochs = i                                                  ! TOTAL NUMBER OF MISSING EPOCHS
     !****************************** STATUS UPDATE ******************************************
     IF (missingepochs == 0) THEN
        CALL status_update('STATUS','GRACEFIT','kb_filtered',' ','No Missing Epochs',ioerr)

     ELSE
        write(message,'(i8,a)')missingepochs,' missing epochs found'
        CALL status_update('STATUS','GRACEFIT','kb_filtered',' ',message,ioerr)
        !****************************** MODELLING EPOCHS ******************************************
        firstm = 0
        lastm = 0
        DO j = 1, end_neq
           IF (missing(j) == -1) THEN
              firstm = missing(j + 1)
           ENDIF
           IF (missing(j) == -2) THEN
              lastm = missing(j-1)
! PT180905: I think if either lastm or firstm are zero then there are no obs to filter ....
              IF (lastm > 0 .and. firstm > 0)then
                 if(prefit(lastm) == 0 .AND. prefit(firstm) == 0) THEN
                   CALL kb_model(firstm, lastm, end_neq, ordered, prefit, ngaps)
                 endif
              ENDIF
           ENDIF
        ENDDO
     ENDIF

     prefit (1:edge) = prefit(edge+1)
     prefit (end_neq-edge:end_neq) = prefit(end_neq-edge)
     edge = 0
     !****************************** FILTERING ALL DATA *************************************
     mean = SUM(prefit(start_neq+edge:end_neq-edge))/DBLE(end_neq-start_neq-edge-edge)         ! FIND THE MEAN OF THE DATA
     DO i= start_neq+edge, end_neq-edge
        prefit(i) = prefit(i) - mean
     ENDDO

     CALL status_update('STATUS','GRACEFIT','kb_filtered',' ','Filtering Data',ioerr)
     nvals = end_neq - start_neq + 1 -edge -edge
     CALL emd(prefit (start_neq+edge: end_neq-edge),nvals,ncomp,decomp)                             ! CALL THE FILTER
     ! PT170728: put the mean back into prefit
     DO i= start_neq+edge, end_neq-edge
        prefit(i) = prefit(i) + mean
     ENDDO

!******************************
! DEBUG: do NOT remove!
if(obs_type == 'kbrr')then
   open(111,file='emd_kbrr.decomp',status='unknown')
else if (obs_type == 'kbra')then
   open(111,file='emd_kbra.decomp',status='unknown')
endif
do i=start_neq,end_neq
  write(111,*)i,decomp(:,ordered(i))
enddo
 close(111)
!******************************

     edge = 4
     !****************************** REMOVING SYNTHESIZED DATA ******************************
     ! PT180803: add a comment regarding how the data were reconstructed
     write(message,'(a,i3,a,i3)')"Reconstructing filtered data from component",start_reconstruct,' to',ncomp
     CALL status_update('STATUS','GRACEFIT','kb_filtered',' ',message,ioerr)
     DO i=start_neq+edge,end_neq-edge
        ! PT170911: use the starting component from the gracefit command file
        if(obs_type == "kbrr")then
          DO j=start_reconstruct,ncomp
            IF(ordered(i) /= 0) THEN                                    ! REMOVE THE SYNTHESIZED DATA
              filtered(i) = filtered(i) + decomp(j,ordered(i))          ! MOVE FILTERED DATA OVER
            ENDIF
          ENDDO
        else
! PT180930: try cropping out the longer-wavelength components of kbra, leaving just the gravity field signals and not the 1/rev
!          DO j=start_reconstruct,4
          DO j=start_reconstruct,ncomp
            IF(ordered(i) /= 0) THEN                                    ! REMOVE THE SYNTHESIZED DATA
              filtered(i) = filtered(i) + decomp(j,ordered(i))          ! MOVE FILTERED DATA OVER
            ENDIF
          ENDDO
        endif
        IF(ordered(i) /= 0) THEN                                      ! USES ORDERED LIST TO ENSURE
           filtered(i) = filtered(i) + mean                          ! NO SYNTHESIZED DATA IS PLACED
        ENDIF                                                         ! INTO FINAL ANSWER
     ENDDO
! PT180817: need to change the logic here, so that "filtered" contains all original data as well as the filtered data between start_new and end_neq
     filtered(1:start_neq+10) = original(1:start_neq+10)
     filtered(end_neq-10:nepochs_t) = original(end_neq-10:nepochs_t)

     WRITE(*,'(a,i5)')'FILTER/missing_epoch: Total Number of Gaps Modelled:       ', ngaps

!! PT181003: if only using coponents 3 and 4 we need to remove the mean of the filtered data (a 1/rev in range rate differentiates to a constant in RA)
!     mean = SUM(filtered(start_neq+10:end_neq-10))/DBLE(end_neq-start_neq-10-10) 
!     filtered(start_neq+10:end_neq-10) = filtered(start_neq+10:end_neq-10) - mean
       
     RETURN

   END SUBROUTINE kb_filtered
   !******************************************************************************************************************

   !********************************************************************************************************************************
   !
   !  kb_model: Model missing epochs, synthesize missing data
   !
   !  Author:  Erin Gray
   !
   ! Contact info:   E. Gray esgray@princeton.edu
   !                 P. Tregoning paul.tregoning@anu.edu.au
   !
   ! Output files:  prefit is input and output, with synthesized data in it
   !
   !********************************************************************************************************************************

   SUBROUTINE kb_model(firstm, lastm, totalrows, ordered, prefit, ngaps)

     IMPLICIT NONE



     !********************  Variable declarations ********************
     INTEGER*4           :: ioerr                                     ! Error number
     INTEGER*8           :: i, k, j                                    ! Counter variables
     INTEGER*4           :: firstm, lastm, totalrows                   ! Row the prefit is in, total columns in file that are missing
     INTEGER*4           :: beginning, ending, nvals                   ! Modelling Data values
     DOUBLE PRECISION    :: prefit(totalrows)                          ! Stores prefit residuals, Prefit OMC (range, range rate, or range acceleration)
     INTEGER*4           :: ordered(totalrows)  , PLT_FIXY , ngaps                      ! Stores Missing Epochs, Missing Order of Epochs
     DOUBLE PRECISION    :: premodel(totalrows)                        ! Prefit Data sent to the model
     DOUBLE PRECISION    :: model(totalrows)                           ! Returned Model Data
     real(kind=8)        :: kb_sponge_model(7)    ! all the parameters for a linear plus sinusoid model


     !****************************** MODEL DATA AROUND MISSING EPOCHS ***********************
     premodel = 0
     CALL status_update('STATUS','FILTER_PREFIT','missing_epochs_status',' ','Modelling Epochs',ioerr)

     beginning = firstm - (42*(lastm-firstm))                    !TAKING DATA AROUND MISSING DATA
     ending = lastm + (42*(lastm-firstm))

     IF ((beginning < 1) .AND. ((ending-beginning) > totalrows)) THEN
        beginning = firstm - 3000
        ending = lastm + 3000
     ENDIF

     IF (beginning < 1) THEN                                        ! IF THE DATA IS TOO CLOSE TO BEG/END OF FILE
        ending = ending - beginning                                  ! THE SAME NUMBER OF EPOCHS WILL STILL BE USED
        beginning = 6                                                ! JUST MOVED TO THE OTHER SIDE OF THE MISSING DATA
     ENDIF                                                          !BEG = 6 FOR CORNER CASES

     IF (ending > totalrows) THEN
        beginning = beginning + (ending - totalrows)
        ending = totalrows - 5                                       ! Total Rows -5 accounts for edge effect
     ENDIF

     ending = ending-1
     nvals = ending - beginning                                     ! TOTAL VALUES TO BE ESTIMATED BASED ON THE RATIO


     j = firstm- 10                                                  ! ACCOUNTS FOR EDGE EFFECT
     k = lastm + 10                                                 ! ACCOUNTS FOR EDGE EFFECT
     IF (j < 1) j = 1                                               ! ACCOUNTS FOR EDGE EFFECT
     IF (k > totalrows) k = totalrows                               ! ACCOUNTS FOR EDGE EFFECT
     DO i = j, k                                                 ! ACCOUNTS FOR EDGE EFFECT
        ordered(i) = 0
     ENDDO

     WRITE(*,'(a,i5)')'FILTER/missing_epoch: Missing Epochs begin at:         ', firstm
     WRITE(*,'(a,i5)')'FILTER/missing_epoch: Missing Epochs end at:           ', lastm
     WRITE(*,'(a,i5)')'FILTER/missing_epoch: Total Missing Epochs :         ', lastm-firstm
     WRITE(*,'(a,i5)')'FILTER/missing_epoch: Modelling Epochs beginning at:  ', beginning
     WRITE(*,'(a,i5)')'FILTER/missing_epoch: Modelling Epochs ending at:      ', ending
     WRITE(*,'(a,i5)')'FILTER/missing_epoch: Total Epochs Modelled:          ', nvals


     i = 1
     DO k = beginning, ending                                    ! MOVE DATA TO BE MODELLED INTO NEW STACK
        premodel(i) = prefit(k)                                      ! PREMODEL CONTAINS ONLY DATA TO BE MODELLED
        i = i +1
     ENDDO

     !****************************** SENDING KBRR TO PRE_FILTER SUBROUTINE ******************
     model = 0
     CALL kb_sponge(nvals,premodel,model,kb_sponge_model)                           ! CALL THE MODELLING PROGRAM
     CALL status_update('STATUS','FILTER_PREFIT','missing_epochs_status',' ','Epochs Successfully Modelled',ioerr)
     CALL status_update('STATUS','FILTER_PREFIT','missing_epochs_status',' ','Data Shuffling',ioerr)



     i = 1
     DO k = beginning, ending
        IF ((prefit(k))== 0) THEN                ! ADD MODELLED DATA TO THE ACTUAL DATA
           prefit(k) = model(i)
        ENDIF
        i = i + 1
     ENDDO

     ngaps = ngaps + 1   ! counter for number of modelled gaps (not total gaps)

     RETURN

   END SUBROUTINE kb_model

   !******************************************************************************************************************
  SUBROUTINE kb_cos_fft(calling_prog,use_obs,obs_type,nepochs_t,nfft,start_neq,end_neq,freqs,prefit,filtered)

! written by Bec McGirr some time in 2019 (not documented when)
!
! MODS
! PT201014: pass in use_obs array to indicate whether to use particular observations
! PT220429: include sig_process for the declaration of n_edge_unusable

  use sig_process
  
  IMPLICIT NONE


!********************  Variable declarations ********************
  character(*),intent(in)     :: calling_prog            ! name of calling program
  logical     ,intent(inout)  :: use_obs(nepochs_t)      ! flag whether to use KBR/LRI range/range rate/range accel  obs
  character*4, intent(in)     :: obs_type                ! 'kbrr' or 'kbra'
  integer*4,   intent(in)     :: nepochs_t               ! total number of epochs
  integer*4,   intent(in)     :: nfft                    ! see note in subroutine header
  integer*4,   intent(in)     :: start_neq,end_neq       ! start and end epochs
  real(kind=8),intent(in)     :: freqs(2)                ! cosine lowpass filter band
  real(kind=8),intent(inout)  :: prefit(nepochs_t)       ! Prefit OMC (range, range rate, or range acceleration)
  real(kind=8),intent(inout)  :: filtered(nepochs_t)     ! reconstructed prefit OMC using cos lowpass filter

! local variables
  integer*4                   :: i
  integer*4                   :: inc                     ! sample spacing
  character*15                :: spec_outfile
  character*17                :: filt_outfile
  character*10                :: ftype                   ! set to 'low'
  character*10                :: wtype                   ! set to 'uniform'
  real(kind=4),allocatable    :: pspec_raw(:)   ! contains unfiltered power spectrum
  real(kind=4),allocatable    :: pspec_filt(:)  ! contains filtered power spectrum
  real(kind=4),allocatable    :: fspec(:)       ! contains frequencies
  real(kind=4),allocatable    :: filt(:)        ! contains filter
  logical                     :: debug                   ! debug flag
  character*100               :: message
  real(kind=8)                :: kb_model(7)             ! all the parameters for a linear plus sinusoid model
  real(kind=8)                :: model(18000)             ! temporary vector to store (but not use) modelled kbr values from kb_sponge
  real(kind=8)                :: mean

! PT190807: define a new vector for the prefit, which adds 1 hour to each end of the prefit.
  real(kind=8),allocatable    :: tmp_prefit(:)
  real(kind=8),allocatable    :: tmp_filtered(:)
  real(kind=4),allocatable    :: window(:)
  integer*4                   :: n_extend            ! number of epochs of extrapolated data to add to each end
  real(kind=8) :: dt, omega, pi
  integer*4    :: iepoch,next_epoch

! variables for infilling data gaps
  logical, allocatable :: infilled(:)
  logical              :: in_gap
  integer*4            :: start_gap,end_gap,ep_start,ep_end
  real(kind=8)         :: offset
  integer*4            :: end_offset
  integer*4            :: ngaps

  pi = 4.d0*datan(1.d0)
! set the period of a sinusoidal model to be roughly 1/rev
  omega = 2.d0*pi/1100.d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                       !
!                Infill data gaps                       !
!                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT190808: run through the data and infill any data gaps.
  allocate(infilled(nepochs_t))
  infilled = .false.

  in_gap = .false.
  start_gap = 0
  end_gap = 0
  ngaps = 0

  iepoch = start_neq
  do while (iepoch < end_neq)
    if(.not. use_obs(iepoch))then  ! we have a data gap
      start_gap = iepoch
      in_gap = .true.
      next_epoch = iepoch+1
      
      ! check forward to find the next non-gap epoch
      do while (in_gap .and. next_epoch <= end_neq)
        !print *, " .... II ", in_gap, next_epoch, use_obs(next_epoch-1: next_epoch)
        if(use_obs(next_epoch) .or. next_epoch == end_neq)then  ! end of gap
          end_gap = next_epoch-1
          in_gap = .false.
          iepoch = next_epoch
        endif
        if (next_epoch == end_neq) then
                end_gap = end_neq
                in_gap = .false.
                iepoch = next_epoch
                ngaps = ngaps + 1
        endif
        next_epoch = next_epoch+1
      enddo

      if(start_gap /= 0 .and. end_gap /= 0  )then    ! we have a gap from start_gap to end_gap
        write(message,'(a,i6,a,i6)')"Interpolating to fill data gap from epoch ",start_gap," to ",end_gap
        call status_update('STATUS',calling_prog,'kb_cos_fft',' ',message,0)

        ! PT220706: treat differently the case of a gap from the start through to the first epoch (ie processing only 1200-2400 data)
        
        ! PT220624: treat differently the case of a gap through to nepochs_t 
        if(end_gap < nepochs_t)then
          ! fit a model to 700 epochs of data, 350 epochs either side of the gap
          ep_start = start_gap - 350
          if(ep_start < 1)ep_start = 1
          ep_end = end_gap + 350
          if(ep_end > nepochs_t)ep_end = nepochs_t
        else
          ep_start = start_gap - 1000
          if(ep_start < 1)ep_start = 1
          ep_end = start_gap - 1
        endif
        call fit_sinusoid_rate(calling_prog,ep_end-ep_start,omega,prefit(ep_start:ep_end),kb_model)
        ! now infill the prefit vector
        ! LD 220218: Removed the -1 after end_gap
        do i=start_gap,end_gap
          if(i == start_gap .and. start_gap > 1)then
          !if(i == start_gap .and. start_gap > 1)then
            dt = start_gap - ep_start - 1
            offset = prefit(start_gap-1) - (kb_model(1) +kb_model(2)*dt + kb_model(3)*dcos(omega*dt) + kb_model(4)*dsin(omega*dt))
            !offset = prefit(start_gap) - (kb_model(1) +kb_model(2)*dt + kb_model(3)*dcos(omega*dt) + kb_model(4)*dsin(omega*dt))
          endif
          dt = dble(i) - ep_start
          prefit(i) = kb_model(1) +kb_model(2)*dt + kb_model(3)*dcos(omega*dt) + kb_model(4)*dsin(omega*dt) + offset

          ! record that these epochs have been infilled          
          infilled(i) = .true.
        enddo
        !LD 22/02/2022: Add correction to force end of infilling to fit to end of gap
        do i=start_gap,end_gap
          if(start_gap < end_neq .and. end_gap < end_neq)then
            prefit(i) = prefit(i) + (((prefit(end_gap+1) - prefit(end_gap))/(end_gap-start_gap)) * (i-start_gap))
          endif
        enddo
      endif

    else
      in_gap = .false.
    endif
    
    ! PT190712: the above logic fails to capture that the value at epoch nepochs_t should not be used for kbra data. Set flag here
    if(iepoch == end_neq .and. obs_type == 'kbra')infilled(iepoch) = .true.
    iepoch = iepoch + 1
    start_gap=0
    end_gap=0
  enddo

! PT220706: there is a need to infill from epoch 1 to epoch start_neq if start_neq /= 1.
  if (start_neq /= 1)then
    write(message,'(a,i7)')'infill from epoch 1 to ',start_neq+1
    call status_update('STATUS',calling_prog,'kb_cos_fft',' ',message,0)

    ! Fit a model from start_neq+2 using 2000 epochs, then extrapolate it back. Epochs start_neq and start_neq+1 are corrupted for kbra.
    call fit_sinusoid_rate(calling_prog,2000,omega,prefit(start_neq+2:start_neq+2001),kb_model)
    dt = 0.d0
    offset = prefit(start_neq+2) - (kb_model(1) +kb_model(2)*dt + kb_model(3)*dcos(omega*dt) + kb_model(4)*dsin(omega*dt))
    do i=1,start_neq+1
      dt = dble(i) - dble(start_neq+2)
      prefit(i) = kb_model(1) +kb_model(2)*dt + kb_model(3)*dcos(omega*dt) + kb_model(4)*dsin(omega*dt) + offset
    enddo
  endif 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! PT190807: set up the temporary prefit vector to be ~2 revolutions longer at each end. Then infill it using a sine model extrapolation.
  n_extend = 2000

! allocate arrays
  allocate(tmp_prefit(-n_extend:nepochs_t+n_extend))
  allocate(tmp_filtered(-n_extend:nepochs_t+n_extend))
  allocate(window(-n_extend:nepochs_t+n_extend))
  allocate(fspec(nfft/2 + 1))
  allocate(filt(nfft/2 + 1))
  allocate(pspec_raw(nfft/2 + 1))
  allocate(pspec_filt(nfft/2 + 1))

  tmp_prefit = 0.d0
  tmp_filtered = 0.d0
  tmp_prefit(1:nepochs_t) = prefit(1:nepochs_t)      ! transfer from prefit to tmp_prefit

! fit a model to the first part of the data. Tests in python showed that using 700 epochs works ok
! PT190814: increased to 2000 epochs, to avoid unrealistically high rates being estimated in the sinusoid model
  call fit_sinusoid_rate(calling_prog,2000,omega,tmp_prefit(1:2000),kb_model)

! extrapolate the model to make the extended prefit data at the start of the data. Infill first 4 epochs, which have been set to zero after kb_NR
! PT220428: number of epochs to infill at the start/end now defined by n_edge_unusable in sig_process_mod
  do i=-n_extend,n_edge_unusable
    if(i == -n_extend)then
      if(obs_type == "kbrr")then
         dt = 1.d0
      else if (obs_type == "kbra")then
         dt = n_edge_unusable + 1     ! PT220429: made this generic rather than hardwired to "6.d0"
      endif
      offset = prefit(6) - (kb_model(1) +kb_model(2)*dt + kb_model(3)*dcos(omega*dt) + kb_model(4)*dsin(omega*dt))
      !offset = 0.d0
    endif
    dt = dble(i)
    tmp_prefit(i) = kb_model(1) +kb_model(2)*dt + kb_model(3)*dcos(omega*dt) + kb_model(4)*dsin(omega*dt) + offset
  enddo

! fit a model to the last part of the data. Tests in python showed that using 700 epochs works ok
! PT190814: increased to 2000 epochs, to avoid unrealistically high rates being estimated in the sinusoid model
! PT220624: we need to use variable "end_neq" to define the end of the observations, NOT nepochs_t, then extend out beyond it
  call fit_sinusoid_rate(calling_prog,2000,omega,tmp_prefit(end_neq-1999:end_neq),kb_model)

! extrapolate the model to make the extended prefit data at the start of the data
! PT190814: we need a different epoch for the offset computation, because kbrr data goes to nepochs_t whereas kbra stops 5 epochs beforehand
  if(obs_type == "kbrr")then
    end_offset = 0.d0
  else
    end_offset = n_edge_unusable-1
  endif
  do i=end_neq-end_offset,nepochs_t+n_extend
    if(i == end_neq-end_offset)then
      dt = 2000 - end_offset
      offset = prefit(end_neq-end_offset) - (kb_model(1) +kb_model(2)*dt + kb_model(3)*dcos(omega*dt) + kb_model(4)*dsin(omega*dt))
    endif
    dt = dble(i) - end_neq + 2000
    tmp_prefit(i) = kb_model(1) +kb_model(2)*dt + kb_model(3)*dcos(omega*dt) + kb_model(4)*dsin(omega*dt) + offset
  enddo

  write(message,'(2a)')"Using KB_COS_FFT on obs type ",obs_type
  call status_update("STATUS",calling_prog,"KB_COS_FFT","",message,0)

! RM190814: calculate mean 
  mean = sum(tmp_prefit(-n_extend:nepochs_t+n_extend),nepochs_t+2*n_extend)/dble(nepochs_t+2*n_extend)

! RM190816: apply window before sending to compute_fft
  wtype = 'hann'
  window = 0.0
  call window_function(calling_prog,nepochs_t+2*n_extend,wtype,window)

! call compute_fft use low pass cosine filter
  ftype = 'low'
  debug = .false.
  call compute_fft(debug,calling_prog,nepochs_t+2*n_extend,nfft,(tmp_prefit(-n_extend:nepochs_t+n_extend)-mean)*window(-n_extend:nepochs_t+n_extend),tmp_filtered(-n_extend:nepochs_t+n_extend),ftype,real(freqs(1),4)*5,real(freqs(2),4)*5,fspec,pspec_raw,pspec_filt,filt)

! PT190807: transfer the required part of the filtered values, without the extrapolated end parts. Don't transfer first/last 5 epochs as they 
!           were interpolated (infilled because the numerical differentiation clobbers them earlier)
!  filtered(6:nepochs_t-20) = tmp_filtered(6:nepochs_t-20)
  filtered(1:nepochs_t) = (tmp_filtered(1:nepochs_t) / window(1:nepochs_t)) + mean

  do i=-n_extend,nepochs_t+n_extend
    !if(obs_type == "kbra")print*,i,window(i),(tmp_prefit(i)-mean),(tmp_prefit(i)-mean)*window(i),tmp_filtered(i),(tmp_filtered(i)/window(i))+mean," ",obs_type
    !if(obs_type == "kbra")print*,i,tmp_prefit(i),tmp_filtered(i)," kbra_tmp_prefit"
  enddo

! PT190808: reset back to zero any infilled observations from the prefit
  do iepoch=start_neq,end_neq
    if(infilled(iepoch))filtered(iepoch) = 0.d0
  enddo

! write out frequencies, power spectrum and filter
  call status_update("STATUS",calling_prog,"KB_COS_FFT","","Writing out spectrum",0)
  write(spec_outfile,'(2a)')"kb_cos_fft.",obs_type
  open(1000,file=spec_outfile,status='unknown')
  write(1000,*)freqs(1),freqs(2)
  do i = 1,nfft/2 + 1
    write(1000,*)fspec(i)/5,pspec_raw(i),pspec_filt(i),filt(i)
  enddo
  close(1000)

! write out prefit and filtered
  call status_update("STATUS",calling_prog,"KB_COS_FFT","","Writing out prefit and filtered",0)
  write(filt_outfile,'(2a)')"kb_filtered.",obs_type
  open(1001,file=filt_outfile,status='unknown')
  write(1001,*)freqs(1),freqs(2)
  do i = 1,nepochs_t
    write(1001,*)i,tmp_prefit(i),tmp_filtered(i)/window(i)
  enddo
  close(1001)

! deallocate arrays
  deallocate(infilled)
  deallocate(tmp_prefit)
  deallocate(tmp_filtered)
  deallocate(window)
  deallocate(fspec)
  deallocate(filt)
  deallocate(pspec_raw)
  deallocate(pspec_filt)

  END SUBROUTINE kb_cos_fft



SUBROUTINE kb_computePartials_vec(rvec, part) !,kbrr_part,kbra_part)
   USE gracefit_mod
   IMPLICIT NONE
 
 
 
   !********************  Variable declarations ********************
 
   ! integer         , intent(in)  :: iepoch                 ! epoch number
   double precision, intent(inout) :: part(nobs_t,nparam_t,nepochs_t)  ! Partials of observables with respect to the parameters
   double precision, intent(in)  :: rvec(6, nsat_t,nepochs_t)    ! Positions, velocities and partials for each satellite
 
   integer          :: isat                       ! Counter for satellites
   integer          :: i,j                        ! Counter variables
   double precision, allocatable :: delta(:,:)       ! Differences in position, velocity and acceleration in x, y and z coordinates
   double precision, allocatable :: kb_theor(:,:)                ! Theoretical Range, Range Rate and Range Acceleration
   double precision :: temp_partials(3, 12, nepochs_t)  ! Partials of range, range rate and range acceleration wrt to obs
   !****************************************************************
 
   allocate(delta(6, nepochs_t))
   allocate(kb_theor(3, nepochs_t))

   delta = rvec(:, 2, :) - rvec(:, 1, :)

   kb_theor(1,:) = dsqrt(delta(1,:)*delta(1,:) + delta(2,:)*delta(2,:) + delta(3,:)*delta(3,:))
   kb_theor(2,:) = (delta(1,:)*delta(4,:) + delta(2,:) * delta(5,:) + delta(3,:) * delta(6,: )) / kb_theor(1,:)

   do isat=1, 2
      do i = 1, 3 ! axis
         ! range partials
         temp_partials(1, 6*(isat-1)+i, : ) = delta(i,:) / kb_theor(1,:) * (-1) ** (isat)
         temp_partials(1, 6*(isat-1)+i+3, : ) = delta(i,:) *0.0

         !range rate
         temp_partials(2, 6*(isat-1)+i, : )   = delta(i+3,:) / kb_theor(1,:) * (-1) ** (isat) - temp_partials(1, 6*(isat-1)+i, :) * kb_theor(2,:) /kb_theor(1,:) 
         temp_partials(2, 6*(isat-1)+i+3, : ) = delta(i,:) / kb_theor(1,:) * (-1) ** (isat)  !same as  temp_partials(1, 6*(isat-1)+ipos, : )
      enddo
   enddo


   IF(nkbr_t /= 0) then
      do i = 1, 12
         do j = 1, nparam_t
            part( ikbr, j,: ) =  part( ikbr, j,: )  + temp_partials(1, i,:) * part(i,j,:)
         ENDDO
      enddo
   endif

   IF(nkbrr_t /= 0) then
      do i = 1, 12
         do j = 1, nparam_t
            part( ikbrr, j,: ) =  part( ikbrr, j, : )  + temp_partials(2, i,:) * part(i,j,:)
         ENDDO
      enddo
   endif

   RETURN
END SUBROUTINE kb_computePartials_vec



  SUBROUTINE kb_computePartials_Amat(rvec, Amat)

! subroutine to compute the partial derivatives of the range rate with respect to all the parameters
! Written originally by Sebastien Allgeyer
! Commented by Paul Tregoning
!
! 13 October 2020
!
! MODS:
! PT201013: change the "part" array to be a two-dimensional array called ATransp. It is the transpose of the A-matrix, for ease of computations in 
!           stacking the normal equations.
!           part(nobs_t,nparams_t,nepochs_t)  changes to   Atransp(nparams_t,nobs_t*nepochs_t)


  USE gracefit_mod
  IMPLICIT NONE
 
 
 
!********************  Variable declarations ********************
 
  double precision, intent(inout) :: Amat(nobs_t*nepochs_t,nparam_t)  ! Partials of observables with respect to the parameters
  double precision, intent(in)  :: rvec(6, nsat_t,nepochs_t)          ! Positions, velocities and partials for each satellite
 
  integer*4          :: isat                         ! Counter for satellites
  integer*4          :: i,j                          ! Counter variables
  double precision, allocatable :: delta(:,:)        ! Differences in position, velocity and acceleration in x, y and z coordinates
  double precision, allocatable :: kb_theor(:,:)     ! Theoretical Range, Range Rate and Range Acceleration
  double precision :: range_partials(12, nepochs_t)  ! Partials of range wrt to obs of satellite positions and velocities
  double precision :: rrate_partials(12, nepochs_t)  ! Partials of range rate wrt to obs of satellite positions and velocities
!****************************************************************
  integer :: iepoch 
  allocate(delta(6, nepochs_t))
  allocate(kb_theor(3, nepochs_t))

! the difference in the coordinates
  delta = rvec(:, 2, :) - rvec(:, 1, :)

! kb_theor(1,:) is the range between the satellites
  kb_theor(1,:) = dsqrt(delta(1,:)*delta(1,:) + delta(2,:)*delta(2,:) + delta(3,:)*delta(3,:))

! kb_theor(2,:) is the range rate between the satellites
  kb_theor(2,:) = (delta(1,:)*delta(4,:) + delta(2,:) * delta(5,:) + delta(3,:) * delta(6,: )) / kb_theor(1,:)

  do isat=1, 2
    do i = 1, 3 ! axis
      ! range partials
      range_partials(6*(isat-1)+i, : ) = delta(i,:) / kb_theor(1,:) * (-1) ** (isat)
      range_partials(6*(isat-1)+i+3, : ) = delta(i,:) *0.0

      !range rate
      rrate_partials(6*(isat-1)+i, : )   = delta(i+3,:) / kb_theor(1,:) * (-1) ** (isat) - range_partials(6*(isat-1)+i, :) * kb_theor(2,:) /kb_theor(1,:) 
      rrate_partials(6*(isat-1)+i+3, : ) = delta(i,:) / kb_theor(1,:) * (-1) ** (isat)  !same as  temp_partials(1, 6*(isat-1)+ipos, : )
    enddo
  enddo


! PT201012: to fill out the ATransp matrix as a 2-dimensional matrix we need an epoch loop here. Use OMP_PARALLEL ?
  do iepoch = 1, nepochs_t
    IF(nkbr_t /= 0) then
      do i = 1, 12
        do j = 1, nparam_t
!          part( ikbr, j,: ) =  part( ikbr, j,: )  + temp_partials(1, i,:) * part(i,j,:)
          Amat(nobs_t*(iepoch-1)+ikbr,j)  =  Amat(nobs_t*(iepoch-1)+ikbr,j)  + range_partials(i,iepoch) * Amat(nobs_t*(iepoch-1)+ikbr,i)
        ENDDO
      enddo
    endif

    IF(nkbrr_t /= 0) then
      do i = 1, 12
        do j = 1, nparam_t
          !part( ikbrr, j,: ) =  part( ikbrr, j, : )  + temp_partials(2, i,:) * part(i,j,:)
          !Amat(nobs_t*(iepoch-1)+ikbrr,j)  = Amat(nobs_t*(iepoch-1)+ikbrr,j) + rrate_partials(i,iepoch) * Amat(nobs_t*(iepoch-1)+ikbrr,i)
          Amat(nobs_t*(iepoch-1)+ikbrr,j)  = Amat(nobs_t*(iepoch-1)+ikbrr,j) + rrate_partials(i,iepoch) *Amat(nobs_t*(iepoch-1)+i,j)
        ENDDO
      enddo
    endif

  enddo  ! end of epoch loop for filling out ATransp

  RETURN
  END SUBROUTINE kb_computePartials_Amat


