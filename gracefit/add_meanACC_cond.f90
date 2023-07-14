!********************************************************************************************************************************

subroutine add_meanACC_cond(apply_shadow_cond,acc,apr_wght,apr_ScaleBias,part, pre_omc)

  ! subroutine to impose conditional equations into the system:
  !
  ! We can expect that the calibrated accelerometer measurements are effectively the same in the same
  ! direction for GRACE A and GRACE B (after a sign swap for along-track and cross-track), at least to
  ! within some nominal value (our tests suggest < 10 nm/s^2 or thereabouts). So, we're going to impose
  ! this as a requirement on the calibrated acceleration means, applying the nominal value as an uncertainty
  ! on the "observation". Thrusts will be either irrelevant or removed with a boxcar filter.
  !
  ! PT130524: the thrusts make a difference of < 1 nm/s^2 using a sliding window of even up to 5 hrs, so they 
  ! are irrelevant (test done on day 2010 07 23)
  !
  ! PT140919: if we pass in apply_shadow_cond then we can make the mean values during shadow (no thrusts) the same
  !
  ! PT170516: the expectation that the radial are the same is flawed, because of the slight change from nominal of
  !           each satellite in pitch. This means that the lead satellite will always be slightly -ive (with an
  !           atmospheric drag component pushing the satellite away from the Earth) and the trailing satellite 
  !           will be slightly positive (with the drag pushing it towards the Earth).
  !
  ! P. Tregoning
  ! 24 May 2013
  use gracefit_mod
  implicit none

   

  !********************  Variable declarations ********************

  logical         , intent(in)  :: apply_shadow_cond(nepochs_t,nsat_t)  ! logical as to whether sat is in shadow AND not affected by thrusts
  double precision, intent(in)  :: acc(nepochs_t,4,nsat_t)      ! Accelerations read in from ACC1B files
  double precision, intent(in)  :: apr_wght(nobs_t)             ! Weights to be applied to data
  double precision, intent(in)  :: apr_ScaleBias(max_SB,maxsat) ! A priori values for parameters
  double precision, intent(out) :: part(nobs_t,nparam_t)        ! Partials of observables with respect to the parameters
  double precision, intent(out) :: pre_omc(nobs_t)              ! Observed Minus Computed values to be used in LS algorithm

  integer*4        :: isat                ! Counter that runs through satellites
  integer*4        :: i,j                 ! Counter variable
  integer*4        :: icount(3,nsat_t)    ! another counter variable
  double precision :: acc_sum(3, nsat_t)  ! Sum of accelerations over every epoch (for each coordinate and each satellite)
  double precision :: mean(3,nsat_t)      ! Mean of accelerations over every epoch (for each coordinate and each satellite)
  !****************************************************************
  ! NOTE:   The condition is as follows  
  !      The assumption is that the mean accelerations of the
  !      two satellites should be more or less the same. The 
  !      conditional equation is:
  !      Scl_Xa*ACC_Xa + Bias_Xa-(Scl_Xb*ACC_Xb + Bias_Xb) = 0
  !
  ! PT160404: In fact, there is pitch up (leading) and down (trailing) of the
  !           satellites so that they are parallel to the line-of-sight. Therefore,
  !           for the SRF_Z direction the mean of each satellite will be different 
  !           in the Z direction. We therefore cannot apply this condition that they
  !           should be the same. 
  !
  !           We can, however, apply the condition that their mean value will be zero
  !
  !      [Scl_Xa*ACC_Xa + Bias_Xa + (Scl_Xb*ACC_Xb + Bias_Xb)] / 2.0 = 0
  !
  ! PT170516: can we, rather, apply the condition that the absolute value of their difference
  !           will be, say 150 nm/s^2 ? The exact magnitude will depend on the strength of the
  !           atmospheric drag, but let's adopt 150 nm/s^2 for now ...
  !
  !      abs( Scl_Xa*ACC_Xa + Bias_Xa - (Scl_Xb*ACC_Xb + Bias_Xb) ] = 150

  !****************************************************************

  ! DEBUG
  acc_sum   = 0.d0
  do isat = 1, nsat_t
     do i = 1, 3
        icount(i,isat) = 0
        do j = 1,nepochs_t
           if(apply_shadow_cond(j,isat))then
              icount(i,isat) = icount(i,isat)+1
              acc_sum(i,isat) = acc_sum(i,isat) + acc(j,i+1,isat)
           endif
        enddo
        mean(i,isat) = acc_sum(i,isat)/dble(icount(i,isat))
     enddo
  enddo

  ! Flip the sign of the mean for GRACE B for the X and Y axes to make them compatible. Radial coord does not change sign
  !    mean(1,2) = -mean(1,2)
  !    mean(2,2) = -mean(2,2)

  !****************************************************************

  ! Set part
  do isat = 1, nsat_t
     ! PT140922: fix bug for the sign of the GRACE B partial for sclx/scly and bsx/bsy. Because the sign of the axes is opposite to GRACE A for these two axes,
     !           the partials for GRACE A and GRACE B should both be +1. The sign of sclz and bsz for GRACE B should be negative.
     !
     ! PT160404: if we include the atmospheric drag components in the z-direction (which come about because of the pitching to align the
     !           satellites to the LOS) then the mean values of Z accelerations for A and B will not be the same, because the drag adds 
     !           positively to the trailing satellite and negatively to the leading satellite. Thus, the mean of the two satellites' Z acc
     !           will be free of drag (and nearly zero). We can write the equation
     !           -10 = (Sclz_A * z_a + bsz_A) + (Sclz_B * z_B + bsz_B) / 2.0
     ! PT160405: actually, we should leave this subroutine alone and do it elsewhere. No code has been changed below.

     do i = 1, 2
        if(nScale_t > 0) part(iacobs+i-1,iScale-1+i+nsatprm_t*(isat-1)) = mean(i,isat)   ! d[equ]/dScl_Xa/b  = +ACC_X a/b (both same sign because of sign change for one satellite flying backwards)
        if(nBias_t > 0)  part(iacobs+i-1,iBias-1+i+nsatprm_t*(isat-1)) = 1.d0            ! d[equ]/dBias_XY a/b = +1.0
     enddo
     i=3
     if(nScale_t > 0)part(iacobs+i-1,iScale-1+i+nsatprm_t*(isat-1)) = (-1)**(isat+1)*mean(i,isat)  ! d[equ]/dScl_Xa/b  = +-ACC_Xa/b
     if(nBias_t > 0) part(iacobs+i-1,iBias-1+i+nsatprm_t*(isat-1))  = (-1)**(isat+1)               ! d[equ]/dBias_Z  a = +1.0, b = -1.0
  enddo  ! End of satellite loop

  ! Set pre_omc
  do i = 1, 3
     if(i < 3)then ! Change the sign of the biases for GRACE A for X and Y, since we've swapped the sign of the mean
        pre_omc(iacobs+i-1) = 0.d0 - ( apr_ScaleBias(i,1)*mean(i,1)+apr_ScaleBias(i+3,1) &
             + (apr_ScaleBias(i,2)*mean(i,2)+apr_ScaleBias(i+3,2)) ) ! Observed is 0.0
     else
        pre_omc(iacobs+i-1) = 0.d0 - ( apr_ScaleBias(i,1)*mean(i,1)+apr_ScaleBias(i+3,1) &
             - (apr_ScaleBias(i,2)*mean(i,2)+apr_ScaleBias(i+3,2)) ) ! Observed is 0.0
     endif
  enddo
  !****************************************************************

  return
end subroutine add_meanACC_cond

!********************************************************************************************************************************
