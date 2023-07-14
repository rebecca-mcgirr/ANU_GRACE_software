!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! debug_lib:  a library of debug subroutines that I might need sometimes
!
! P. Tregoning
! 11 May 2014

subroutine kbrr_LS(iepoch,X0_adj,pre_omc,part,kbrr_part)

  ! X0_adj    : adjustments to pos/vel at the first epoch (ie the "parameter" adjustments that we are seeking eventually)
  ! pre_omc   :  prefit OMC for pos/vel at this particular epoch
  ! kbrr_part : dRR/dX_t, dRR/dY_t .... dRR/dZdot_t for this epoch

  implicit none

  integer*4   , intent(in)  :: iepoch
  real(kind=8), intent(in)  :: X0_adj(12)
  real(kind=8), intent(in)  :: pre_omc(13)
  real(kind=8), intent(in)  :: part(13,24)
  real(kind=8), intent(in)  :: kbrr_part(6,2)

  integer*4 :: i,isat,j
  real(kind=8) :: delta_kbrr,delta_IC(6,2),delta_kbrr_bias,IC_adj(24),tmp_pre_omc



  ! set the IC adjust to be the prefit residuals at epoch 1 and zero for biases
  IC_adj = 0.d0
  IC_adj(1:6) = X0_adj(1:6)
  IC_adj(13:18) = X0_adj(7:12)  

  ! first, compute the "error" in range rate implied by the prefit residuals in position/velocity
  delta_kbrr = 0.d0
  ! first, compute the "error" in range rate implied by the prefit residuals in position/velocity
  delta_kbrr = 0.d0
  delta_kbrr_bias = 0.d0
  do isat = 1,2
     do i=1,6
        delta_kbrr = delta_kbrr + pre_omc((isat-1)*6+i)*kbrr_part(i,isat)
     enddo
     ! PT140612: do it again using the IC adjustments and the dP_t/dP_0 partials, also including the scale and bias ...
     do i=1,6  
        tmp_pre_omc = 0.d0
        do j=1,12
           tmp_pre_omc = tmp_pre_omc + IC_adj((isat-1)*12+j)*part((isat-1)*6+i,(isat-1)*12+j)
           delta_kbrr_bias = delta_kbrr_bias + kbrr_part(i,isat)*part((isat-1)*6+i,(isat-1)*12+j) * IC_adj((isat-1)*12+j)
        enddo
        ! also compute the position and velocity adjustment implied by the d/dP partial and "error" of the pos/vel at the first epoch
        do j=1,6
           if(i <= 3)then
              delta_IC(i,isat) = delta_IC(i,isat) + X0_adj((isat-1)*6+j)*part((isat-1)*6+i,(isat-1)*12+j)
           else
              delta_IC(i,isat) = delta_IC(i,isat) + X0_adj((isat-1)*6+j)*part((isat-1)*6+i,(isat-1)*12+j)
           endif
        enddo
     enddo
  enddo


!!! PT140613: make a LS solution to solve for the bias (and scale?) adjustments that would have made this epoch compatible with the GNV1B and the KBR1B obs
  !    if(iepoch > 1) call est_bias(iepoch,kbrr_part,part,IC_adj,pre_omc)
  !    if(iepoch > 1) call est_pos(iepoch,kbrr_part,part,IC_adj,pre_omc)

  return
end subroutine kbrr_LS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine est_bias(iepoch,kbrr_part,part,IC_adj,pre_omc)

  ! subroutine to estimate the biases - and possibly scales - of the accelerometers at a single epoch. The following assumptions are made
  ! in order to do this:
  !
  ! 1. We know accurately the IC adjustments for pos/vel, being the prefit residuals for the first epoch. These numbers seem to be consistent
  !    with the kbrr prefit residual, so let's use them.
  ! 2. All partial derivatives relating ICs to pos/vel at any epoch are correct.
  ! 3. All partials relating IC bias adjustments to pos/vel at any epoch are correct.
  ! 4. The prefit pos/vel/kbrr at each epoch then become our observations, with parameters of biases (and maybe scales)

  ! delta_X = dX/dXo delta_Xo + ... + dX/dZdot_o delta_Zdot_o + dX/dBx delta_bx + ....
  ! and for range rate
  ! delta_RR = dRR/dX x dX/dX_o x delta_Xo + ..... + dRR/dX x dX/dBx x delta_Bx + .... + dRR/dZdot x dZdot/dBx x delta_bx
  !
  ! We order the observations as XYZ velXYZ for GRACE A, XYZ velXYZ for GRACE B, KBRR. We order the parameters as Bx, By, Bz

  implicit none

  real(kind=8) :: kbrr_part(6,2)           !  partials relating instantaneous pos/vel to range rate
  real(kind=8) :: part(13,24)              !  partials relating ICs to instantaneous pos/vel
  real(kind=8) :: IC_adj(24)               !  adjustments to the ICs, being the prefit residuals at the first epoch. Don't trust the bias/scale IC adjustments here !!!!
  real(kind=8) :: pre_omc(13)              !  prefit residuals of pos/vel/RR at this epoch
  integer*4    :: iepoch

  integer*4    :: isat,i,j,k

  ! a priori values for parameters
  real(kind=8) :: bias_apr(12)             ! initially use only 6 but make it 12 in case we estimate scale later on
  real(kind=8) :: computed

  ! LS variables
  real(kind=8) :: A(13,6)
  real(kind=8) :: B(13,1)
  real(kind=8) :: At(6,13),AtW(6,13),AtWA(6,6),VCV(6,6),AtWB(6,1),soln(6,1),W(13,13)
  integer*4    :: iter
  real(kind=8) :: wgt_gps_pos,wgt_gps_vel,wgt_kbrr

  !  call printmat(part,13,24,'part ')

  ! set a priori values for biases (and scales) to be zero
  bias_apr = 0.d0
  A = 0.d0

  ! for a test, set all the IC_adj to zero, so that only the partials related to bias are included in the LS.
  !    IC_adj = 0.d0

  do iter = 1,5
!!!!!!!!!!! A  Matrix !!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! first, the pos/vel partials for each satellite
     do i=1,6                    ! number of observations to each GRACE satellite
        do j=1,3                  ! number of bias parameters for each satellite
           A(i,j) = part(i,9+j)          ! biases for GRACE A
           A(i+6,j+3) = part(i+6,21+j)   ! biases for GRACE B
        enddo
     enddo

     ! now the range rate
     do i=1,6
        do j=1,3
           A(13,j)   = kbrr_part(i,1) * part(i,9+j)    ! kbrr partial for bias GRACE A
           A(13,j+3) = kbrr_part(i,2) * part(i+6,21+j) ! kbrr partial for bias GRACE B
        enddo
     enddo

!!!!!!!!!!!!! Weight matrix !!!!!!!!!!!!!!!!!!!!!
     wgt_gps_pos = 1.d0
     wgt_gps_vel = 1.d0
     wgt_kbrr    = 1.d0
     wgt_gps_pos = 3.d-2
     wgt_gps_vel = 1.d-5
     wgt_kbrr    = 1.5d-6
     !    wgt_gps_pos = 1.d-2
     !    wgt_gps_vel = 1.d-5
     !    wgt_kbrr    = 1.d-6

     do i=1,3
        W(i,i) = 1.d0/wgt_gps_pos**2
        W(i+6,i+6) = 1.d0/wgt_gps_pos**2
        W(i+3,i+3) = 1.d0/wgt_gps_vel**2
        W(i+9,i+9) = 1.d0/wgt_gps_vel**2
     enddo
     W(13,13) = 1.d0/wgt_kbrr**2


!!!!!!!!!!!!!  OMC  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     B = 0.d0
     ! GRACE A
     do i=1,6
        computed = 0.d0
        !     pos/vel component
        do j=1,6
           computed =  computed + part(i,j)*IC_adj(j)
        enddo
        !     bias component
        do j=1,3
           computed = computed + part(i,9+j)*bias_apr(j)
        enddo
        B(i,1) = pre_omc(i) - computed
     enddo

     ! GRACE B
     do i=1,6
        computed = 0.d0
        !     pos/vel component
        do j=1,6
           computed =  computed + part(i+6,j+12)*IC_adj(j+12)
        enddo
        !     bias component
        do j=1,3
           computed = computed + part(i+6,21+j)*bias_apr(j+3)
        enddo
        B(i+6,1) = pre_omc(i+6) - computed
     enddo

     ! KBRR
     B(13,1) = pre_omc(13) 
     do i=1,6
        do j=1,6
           ! pos/vel component GRACE A
           B(13,1) = B(13,1) - kbrr_part(i,1)*part(i,j)*IC_adj(j) 
        enddo
     enddo

     ! pos/vel component GRACE B
     do i=7,12
        do j=7,12
           B(13,1) = B(13,1) - kbrr_part(i-6,2)*part(i,j+6)*IC_adj(j) 
        enddo
     enddo

     ! bias component GRACE A
     do k=1,3
        do i=1,6
           B(13,1) = B(13,1) - kbrr_part(i,1)*part(i,k+9)*bias_apr(k) 
        enddo
     enddo

     ! bias component GRACE B
     do k=1,3
        do i=1,6
           B(13,1) = B(13,1) - kbrr_part(i,2)*part(i+6,k+21)*bias_apr(k+3)
        enddo
     enddo
!!!!!!!!!!!!!  END of OMC  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!




     print*,"KBRR partials:",iepoch,A(13,:)

     ! now do the LS solution
     call transp(A, At,13,6)
     call matmult(At,W,AtW,6,13,13)
     call matmult(AtW,A,AtWA,6,13,6)
     ! DEBUG
     if(iter > 0) then
        call printmat (A,13,6,"Amat ")
        call printmat (B,13,1,"Bmat ")
        !      call printmat(At,6,13,'At ')
        !      call printmat(AtW,6,13,'AtW ')
        call printmat (AtWA,6,6,"AtWA ")
        !    call printmat(soln,6,1,"soln ")
     endif
     call invert(AtWA,VCV,6,6)
     call printmat(VCV,6,6,"VCV ")
     call matmult(AtW,B,AtWB,6,13,1)
     call printmat (AtWB,6,1,"AtWB ")
     call matmult(VCV,AtWB,soln,6,6,1)

     ! DEBUG
     if(iter > 10) then
        call printmat (A,13,6,"Amat ")
        call printmat (B,13,1,"Bmat ")
        !      call printmat(At,6,13,'At ')
        !      call printmat(AtW,6,13,'AtW ')
        call printmat (AtWA,6,6,"AtWA ")
        !    call printmat(soln,6,1,"soln ")
     endif

     do i=1,6
        if(iter > 0)print*,'Epoch ',iepoch,'iter ',iter,' Bias ',i,bias_apr(i),soln(i,1),VCV(i,i)
        bias_apr(i) = bias_apr(i) + soln(i,1)
     enddo

  enddo

  return
end subroutine est_bias
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine est_pos(iepoch,kbrr_part,part,IC_adj,pre_omc)

  ! subroutine to estimate the position adjustments at a single epoch. It is based on est_bias, but the parameters are
  ! positions of GRACE A and B rather than biases

  implicit none

  real(kind=8) :: kbrr_part(6,2)           !  partials relating instantaneous pos/vel to range rate
  real(kind=8) :: part(13,24)              !  partials relating ICs to instantaneous pos/vel
  real(kind=8) :: IC_adj(24)               !  adjustments to the ICs, being the prefit residuals at the first epoch. Don't trust the bias/scale IC adjustments here !!!!
  real(kind=8) :: pre_omc(13)              !  prefit residuals of pos/vel/RR at this epoch
  integer*4    :: iepoch

  integer*4    :: isat,i,j,k

  ! a priori values for parameters
  real(kind=8) :: bias_apr(12)             ! initially use only 6 but make it 12 in case we estimate scale later on
  real(kind=8) :: computed

  ! LS variables
  real(kind=8) :: A(13,6)
  real(kind=8) :: B(13,1)
  real(kind=8) :: At(6,13),AtW(6,13),AtWA(6,6),VCV(6,6),AtWB(6,1),soln(6,1),W(13,13)
  integer*4    :: iter
  real(kind=8) :: wgt_gps_pos,wgt_gps_vel,wgt_kbrr

  !  call printmat(part,13,24,'part ')

  ! set a priori values for biases (and scales) to be zero
  bias_apr = 0.d0
  A = 0.d0

  ! for a test, set all the IC_adj to zero, so that only the partials related to bias are included in the LS.
  IC_adj = 0.d0

  do iter = 1,5
!!!!!!!!!!! A  Matrix !!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! first, the pos/vel partials for each satellite
     do i=1,6                    ! number of observations to each GRACE satellite
        do j=1,3                  ! number of bias parameters for each satellite
           A(i,j) = part(i,9+j)          ! biases for GRACE A
           A(i+6,j+3) = part(i+6,21+j)   ! biases for GRACE B
        enddo
     enddo

     ! now the range rate
     do i=1,6
        do j=1,3
           A(13,j)   = kbrr_part(i,1) * part(i,9+j)    ! kbrr partial for bias GRACE A
           A(13,j+3) = kbrr_part(i,2) * part(i+6,21+j) ! kbrr partial for bias GRACE B
        enddo
     enddo

!!!!!!!!!!!!! Weight matrix !!!!!!!!!!!!!!!!!!!!!
     wgt_gps_pos = 1.d0
     wgt_gps_vel = 1.d0
     wgt_kbrr    = 1.d0
     wgt_gps_pos = 3.d-2
     wgt_gps_vel = 1.d-5
     !    wgt_kbrr    = 1.5d-7
     !    wgt_gps_pos = 1.d-2
     !    wgt_gps_vel = 1.d-5
     !    wgt_kbrr    = 1.d-6

     do i=1,3
        W(i,i) = 1.d0/wgt_gps_pos**2
        W(i+6,i+6) = 1.d0/wgt_gps_pos**2
        W(i+3,i+3) = 1.d0/wgt_gps_vel**2
        W(i+9,i+9) = 1.d0/wgt_gps_vel**2
     enddo
     W(13,13) = 1.d0/wgt_kbrr**2


!!!!!!!!!!!!!  OMC  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     B = 0.d0
     ! GRACE A
     do i=1,6
        computed = 0.d0
        !     pos/vel component
        do j=1,6
           computed =  computed + part(i,j)*IC_adj(j)
        enddo
        !     bias component
        do j=1,3
           computed = computed + part(i,9+j)*bias_apr(j)
        enddo
        B(i,1) = pre_omc(i) - computed
     enddo

     ! GRACE B
     do i=1,6
        computed = 0.d0
        !     pos/vel component
        do j=1,6
           computed =  computed + part(i+6,j+12)*IC_adj(j+12)
        enddo
        !     bias component
        do j=1,3
           computed = computed + part(i+6,21+j)*bias_apr(j+3)
        enddo
        B(i+6,1) = pre_omc(i+6) - computed
     enddo

     ! KBRR
     B(13,1) = pre_omc(13) 
     do i=1,6
        do j=1,6
           ! pos/vel component GRACE A
           B(13,1) = B(13,1) - kbrr_part(i,1)*part(i,j)*IC_adj(j) 
        enddo
     enddo

     ! pos/vel component GRACE B
     do i=7,12
        do j=7,12
           B(13,1) = B(13,1) - kbrr_part(i-6,2)*part(i,j+6)*IC_adj(j) 
        enddo
     enddo

     ! bias component GRACE A
     do k=1,3
        do i=1,6
           B(13,1) = B(13,1) - kbrr_part(i,1)*part(i,k+9)*bias_apr(k) 
        enddo
     enddo

     ! bias component GRACE B
     do k=1,3
        do i=1,6
           B(13,1) = B(13,1) - kbrr_part(i,2)*part(i+6,k+21)*bias_apr(k+3)
        enddo
     enddo
!!!!!!!!!!!!!  END of OMC  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     ! DEBUG
     !    call printmat (A,13,6,"Amat ")
     if(iter > 0) call printmat (B,13,1,"Bmat ")

     print*,"KBRR partials:",iepoch,A(13,:)

     ! now do the LS solution
     call transp(A, At,13,6)
     call matmult(At,W,AtW,6,13,13)
     call matmult(AtW,A,AtWA,6,13,6)
     call invert(AtWA,VCV,6,6)
     call matmult(AtW,B,AtWB,6,13,1)
     call matmult(VCV,AtWB,soln,6,6,1)

     !    call printmat(soln,6,1,"soln ")

     do i=1,6
        if(iter > 0)print*,'Epoch ',iepoch,'iter ',iter,' Bias ',i,bias_apr(i),soln(i,1)
        bias_apr(i) = bias_apr(i) + soln(i,1)
     enddo

  enddo

  return
end subroutine est_pos















