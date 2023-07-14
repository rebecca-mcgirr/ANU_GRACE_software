subroutine rpy_AOC(rvecA,rvecB,quatA,quatB,rpy,AOC)

  ! subroutine to compute the yaw, pitch and roll angles given the two satellite positions and the quaternions.
  ! It then also computes the antenna offset correction (AOC) for each satellite
  !
  ! P. Tregoning
  ! 11 April 2013
  !
  ! IN:
  !     rvecA     : GRACE A GTORB GPS position (XYZ e-fixed) and epoch (grace seconds) for each epoch) (4,maxepc)
  !     rvecB     : GRACE B GTORB GPS position (XYZ e-fixed) and epoch (grace seconds) for each epoch) (4,maxepc)
  !     quatA     : GRACE A quaternions (SRF to ERF) for each epoch (4,maxepc)
  !     quatB     : GRACE B quaternions (SRF to ERF) for each epoch (4,maxepc)
  !
  ! OUT:
  !     rpy       : yaw, pitch, roll angles at each epoch for which we have GPS positions
  !     AOC       : antenna offset range corrections for each satellite at each epoch

  ! The AOC matrix has 9 columns. The rows are the epochs.
  ! PT140807: I think it is actually 9 rows, 3 columns (range, rr, ra) and the third dimension is epoch!
  ! col 1,2  : AOC for GRACE A,B respectively
  ! col 3-5  : dAOC/d(ang) for GRACE A, where "ang" = roll, pitch, yaw respectively
  ! col 6-8  : dAOC/d(ang) for GRACE B
  ! col 9    : epoch number (not sure if we need this ...)
  !
  ! the third dimension of AOC allows us to store AOC, AORC and AOAC (range, range rate and range acceleration) corrections.
  use gracefit_mod
  implicit none


  double precision, intent(in) :: rvecA(3,nepochs_t),rvecB(3,nepochs_t),quatA(4,nepochs_t),quatB(4,nepochs_t)

  double precision, intent(out) :: rpy(3,nsat_t,nepochs_t),AOC(9,3,nepochs_t)

  double precision quattmp_A(4),quattmp_B(4),quat_rpy(4),quat_c(4),dx(2),vec1(3),vec2(3),latr,lonr,rad
  double precision LOS_mag,LOS_vec(3),ant_mag(2),dANT(3,2),tmp_ant_vec(3),rot_ant_vec(3),rot_LOS_vec(3), quat_tmp(4)
  double precision roll_deriv_A(3,3),pitch_deriv_A(3,3),yaw_deriv_A(3,3),rotmat_A(3,3)
  double precision roll_deriv_vec_A(3),pitch_deriv_vec_A(3),yaw_deriv_vec_A(3)
  double precision roll_deriv_B(3,3),pitch_deriv_B(3,3),yaw_deriv_B(3,3),rotmat_B(3,3)
  double precision roll_deriv_vec_B(3),pitch_deriv_vec_B(3),yaw_deriv_vec_B(3)
  double precision rpy_A(3),rpy_B(3),delta_rpy_A(3),delta_rpy_B(3)
  integer i,j,ioerr,iepoch
  integer num_axes
  character*100 message
  logical debug

  real*8 amag3,LOS(3)

  debug = .false.


  ! APP13-417: define the KBR antenna offsets. These are taken from the VKB1B file

  ant_mag(1) = 1.445119117290785d0
  dANT(1, 1) = 0.9999987139532216d0 * ant_mag(1)
  dANT(2, 1) = -0.0002929197980534523d0 * ant_mag(1)
  dANT(3, 1) = 0.001576797353751629d0 * ant_mag(1)

  ant_mag(2) = 1.444390929005312d0
  dANT(1, 2) = 0.9999973040502859d0 * ant_mag(2)
  dANT(2, 2) = 0.0003988672930788539d0 * ant_mag(2)
  dANT(3, 2) = 0.002287530774148089d0 * ant_mag(2)

  ! Default offsets for roll, pitch and yaw are set to 0 for both satellites, these values should be taken from input and handed off to
  ! this subroutine 
  delta_rpy_A = 0.0d0
  delta_rpy_B = 0.0d0


  ! now we need a loop through all the epochs. Both GTORB files should have the same epochs, so we can just loop through ..... (yeah, sure!)
  if(debug)print*,'rpy_AOC: rvec1 rvec2',rvecA(1:3,1),rvecB(1:3,1)
  do iepoch = 1,  nepochs_t

     if(rvecA(1,iepoch) /= 0.d0 .and. rvecB(1,iepoch) /= 0.d0)then   ! we are checking whether either satellite has a X coordinate

        ! compute the quaternion for the LOS to SRF for GRACE A
        vec1(1:3) = rvecA(1:3,iepoch)
        vec2(1:3) = rvecB(1:3,iepoch)
        LOS(:) = vec2 - vec1
        if(debug)print*,'**pos sat1**:',iepoch, vec1(1:3)
        if(debug)print*,'e-fixed LOS and unit LOS',LOS,LOS/amag3(LOS)
        call srf_2_los(debug,vec1,vec2,quattmp_A,1)   ! 1 for SRF to LOS, -1 for LOS to SRF
        if(debug)print*,'quattmp_A SRF to EF:',quattmp_A

        ! difference it from the SCA quaternion
        quat_c = quatA(1:4,iepoch)
        call quat_conj(quat_c)
        call quat_mul(quat_c,quattmp_A,quat_rpy)
        if(debug)print*,'quat_rpy (= SCA quat - SRF-to_EF quat)',quat_rpy,iepoch

        !APP130419: inserted more complete formulation of AOC

        ! calculate unit LOS vector from satellite A to satellite B
        LOS_mag = dsqrt((vec2(1)-vec1(1))**2 + (vec2(2)-vec1(2))**2 + (vec2(3)-vec1(3))**2)
        LOS_vec = (vec2 - vec1)
        LOS_vec = LOS_vec/LOS_mag

        ! apply roll pitch yaw rotation to the antenna offset vector
        tmp_ant_vec(1:3) = dANT(1:3, 1)
        quat_c = quatA(1:4,iepoch)
        call quat_rot_vect(quat_c, tmp_ant_vec, rot_ant_vec)

        ! APP130716: convert LOS vector to SRF by applying conjugate of SRF2LOS quaternion
        quat_c = quattmp_A
        call quat_conj(quat_c)
        call quat_rot_vect(quat_c, LOS_vec, rot_LOS_vec)

        ! convert to roll, pitch, yaw
        call quat_to_rpy_new(quat_rpy,rpy(1,1,iepoch),rpy(2,1,iepoch),rpy(3,1,iepoch))
        if(debug)print*,'rpy:',rpy(:,1,iepoch),iepoch
        ! incorporate any biases
        rpy_A(1:3) = rpy(1:3,1,iepoch) 
        rpy_A = rpy_A + delta_rpy_A
        rpy(1:3,1,iepoch) = rpy_A(1:3) 

        ! APP130419: calculate individual roll, pitch and yaw rotation matrices and their derivatives
        call rpy2rotmat(rpy_A(1), rpy_A(2), rpy_A(3), rotmat_A, roll_deriv_A, pitch_deriv_A, yaw_deriv_A)

        ! Apply the differentiated rotation matrices to the entenna offset vector
        call matmult( roll_deriv_A, tmp_ant_vec,  roll_deriv_vec_A, 3, 3, 1)
        call matmult(pitch_deriv_A, tmp_ant_vec, pitch_deriv_vec_A, 3, 3, 1)
        call matmult(  yaw_deriv_A, tmp_ant_vec,   yaw_deriv_vec_A, 3, 3, 1)

        ! take dot product of rotated antenna offset with unit LOS vector (note sign change required for satellite B)
        ! take dot product of differentiated antenna offsets with unit LOS vector to calculate dAOC/dangle
        AOC(1,1,iepoch) = 0.0d0
        AOC(3,1,iepoch) = 0.0d0
        AOC(4,1,iepoch) = 0.0d0
        AOC(5,1,iepoch) = 0.0d0
        do num_axes = 1, 3
           AOC(1,1,iepoch) = AOC(1,1,iepoch) + LOS_vec(num_axes)     * rot_ant_vec(num_axes)
           AOC(3,1,iepoch) = AOC(3,1,iepoch) + rot_LOS_vec(num_axes) * roll_deriv_vec_A(num_axes)/LOS_mag
           AOC(4,1,iepoch) = AOC(4,1,iepoch) + rot_LOS_vec(num_axes) * pitch_deriv_vec_A(num_axes)/LOS_mag
           AOC(5,1,iepoch) = AOC(5,1,iepoch) + rot_LOS_vec(num_axes) * yaw_deriv_vec_A(num_axes)/LOS_mag
        enddo

        ! now again for GRACE B
        ! compute the quaternion for the LOS to SRF for GRACE A
        vec1(1:3) = rvecB(1:3,iepoch)
        vec2(1:3) = rvecA(1:3,iepoch)
        call srf_2_los(debug,vec1,vec2,quattmp_B,1)      ! 1 for SRF to LOS, -1 for LOS to SRF

        ! difference it from the SCA quaternion
        quat_c = quatB(1:4,iepoch)
        call quat_conj(quat_c)
        call quat_mul(quat_c,quattmp_B,quat_rpy)

        ! PT140827: try settng this to be larger if it is too small .....
        !        if(quat_rpy(1) < 0.99999e0 )then
        !          print*,'quat_rpy(1) too small .... setting to 0.999997e0 epoch',iepoch
        !!          quat_rpy(1) = 0.999997e0
        !        endif

        ! convert to roll, pitch, yaw
        call quat_to_rpy_new(quat_rpy,rpy(1,nsat_t,iepoch),rpy(2,nsat_t,iepoch),rpy(3,nsat_t,iepoch))

        !APP130418: correct formulation of AOC

        ! apply roll pitch yaw rotation to the antenna offset vector
        tmp_ant_vec(1:3) = dANT(1:3, 2)
        LOS_vec = -LOS_vec
        quat_c = quatB(1:4,iepoch)
        call quat_rot_vect(quat_c, tmp_ant_vec, rot_ant_vec)

        ! APP130716: convert LOS vector to SRF by applying conjugate of SRF2LOS quaternion
        quat_c = quattmp_B
        call quat_conj(quat_c)
        call quat_rot_vect(quat_c, LOS_vec, rot_LOS_vec)

        ! incorporate any biases
        rpy_B(1:3) = rpy(1:3,nsat_t,iepoch) 
        rpy_B = rpy_B + delta_rpy_B
        rpy(1:3,nsat_t,iepoch) = rpy_B(1:3) 

        ! APP130419: calculate individual roll, pitch and yaw rotation matrices and their derivatives
        call rpy2rotmat(rpy_B(1), rpy_B(2), rpy_B(3), rotmat_B, roll_deriv_B, pitch_deriv_B, yaw_deriv_B)

        ! Apply the differentiated rotation matrices to the entenna offset vector
        call matmult( roll_deriv_B, tmp_ant_vec,  roll_deriv_vec_B, 3, 3, 1)
        call matmult(pitch_deriv_B, tmp_ant_vec, pitch_deriv_vec_B, 3, 3, 1)
        call matmult(  yaw_deriv_B, tmp_ant_vec,   yaw_deriv_vec_B, 3, 3, 1)

        ! take dot product of rotated antenna offset with unit LOS vector (note sign change of LOS required for satellite B)
        ! take dot product of differentiated antenna offsets with unit LOS vector to calculate dAOC/dangle
        AOC(2,1,iepoch) = 0.0d0
        AOC(6,1,iepoch) = 0.0d0
        AOC(7,1,iepoch) = 0.0d0
        AOC(8,1,iepoch) = 0.0d0

        ! DEBUG_B
        !  print*,'iepoch,    LOS_vec',iepoch,LOS_vec
        !  print*,'iepoch,rot_ant_vec',iepoch,rot_ant_vec
        !  print*,'iepoch,tmp_ant_vec',iepoch,tmp_ant_vec
        !  print*,'quatB(1:4,iepoch) ',iepoch,quatB(1:4,iepoch)

        LOS_vec(2) = LOS_vec(2) + 0.00025
        !  print*,'iepoch,    LOS_vec',iepoch,LOS_vec
        do num_axes = 1, 3
           ! PT140807: removed division by LOS_mag from the AOC correction (this was (correctly) missing for GRACE A, wrongly present for GRACE B)
           AOC(2,1,iepoch) = AOC(2,1,iepoch) + LOS_vec(num_axes)     * rot_ant_vec(num_axes)
           !  print*,'AOC(2,1,iepoch)',iepoch,num_axes,AOC(2,1,iepoch), LOS_vec(num_axes)     * rot_ant_vec(num_axes)
           AOC(6,1,iepoch) = AOC(6,1,iepoch) + rot_LOS_vec(num_axes) * roll_deriv_vec_B(num_axes)/LOS_mag
           AOC(7,1,iepoch) = AOC(7,1,iepoch) + rot_LOS_vec(num_axes) * pitch_deriv_vec_B(num_axes)/LOS_mag
           AOC(8,1,iepoch) = AOC(8,1,iepoch) + rot_LOS_vec(num_axes) * yaw_deriv_vec_B(num_axes)/LOS_mag
        enddo

        ! store the epoch, in case we need it (not sure if we will)
        AOC(9,1,iepoch) = iepoch
        !        print*,'iepoch,AOC(2,1,iepoch)',iepoch,AOC(2,1,iepoch)
     endif   ! end of check on whether there were coordinates for one satellite

  enddo  ! end of epoch loop

  !    stop ' stopped in rpy_AOC'
  ! now numerically differentiate the AOC elements
  call compute_AORC(AOC)

  return
end subroutine rpy_AOC
 
