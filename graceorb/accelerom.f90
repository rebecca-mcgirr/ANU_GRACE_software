   subroutine accelerom(jd, tin, bias_perturb,c0x,c0y,c0z,c1x,c1y,c1z,c2x,c2y,c2z,scl,scl_perturb,rpy_perturb,acc_if)

! Emma-Kate Potter, 18 August, 2010
! This subroutine calculates the contribution of the non-gravitational accelerations
! as measured by the accelerometers on board each satellite

! The accelerometer data is given in the science reference frame
!
! MODS
! PT/APP 130710: cleaned up the sign convention for the quaternions. Removed the change of sign of the input SCA quaternion, then use
!                the sign convention of the lib/quat_lib.f90 routines.
! PT140522: update the way the accelerometer bias and scale are applied so that it looks like the way the partials are
!           done. The net effect is the same, but the equation will then match code in ICpart_calc.f90
! PT140813: pass in the bias/scale information rather than using the include files, so that we can perturb the values for the perturbed orbits

    use accel_mod
    use accred_mod
! PT140813: pass in the information instead of using this
!    use bsscl_mod
    use rotation_mod

    implicit none

    integer*4,   intent(in ) :: jd               ! julian day
    real(kind=8),intent(in ) :: tin              ! seconds of day
    real(kind=8),intent(in)  :: bias_perturb(3)  ! perturbation to the bias (for calculating partials)
    real(kind=8),intent(in)  :: rpy_perturb(3)   ! perturbation to the roll, pitch, yaw (for calculating partials)
    real(kind=8),intent(in)  :: scl_perturb(3)   ! perturbation to the scale (for calculating partials)
    real(kind=8),intent(in)  :: c0x(2), c0y(2), c0z(2), c1x(2), c1y(2), c1z(2), c2x(2), c2y(2), c2z(2)   ! bias information
    real(kind=8),intent(in)  :: scl(3)           ! accelerometer scales
    real(kind=8),intent(out) :: acc_if(3)        ! non-gravitational, calibrated accelerometer values in inertial XYZ directions   

    character*356 :: message
    integer*4 :: i, j,  ioerr1, ioerr2
    integer*4,save :: t1, t2
    integer*4 :: t1i, t1f, t1tmp, t1out
    integer*4 :: t2i, t2f, t2tmp, t2out
    integer*4 :: k
    real(kind=8) ::  ACC_time 
    real(kind=8) ::  t2f_real, t2i_real 
    real(kind=8) ::  t1f_real, t1i_real 
    real(kind=8), dimension(4) ::  vari
    real(kind=8), dimension(4) ::  varf
    real(kind=8), dimension(4) ::  varout
    real(kind=8), dimension(4) ::  vartmp
    real(kind=8), dimension(0:3) :: quat_corr,quat_tmp
    real(kind=8), dimension(3) ::  acci
    real(kind=8), dimension(3) ::  accf
!    real(kind=8), dimension(3) ::  accout
! APP130328: Change name of observed accelerometer value from accout to accobs which is defined in accel_mod

    real(kind=8), dimension(3) ::  acctmp
!  Moved acc and Qrow to accel_mod so that the values can be used in partials
    real(kind=8), dimension(3,1) :: accrot 
    real(kind=8), dimension(3) :: accmtr 
!    real(kind=8), dimension(3) :: bs
    integer :: flag1 
    integer :: flag2 
    real(kind=8) ::  w, x, y, z, tmp

! PT140818: variables for adding noise to the yaw and pitch rotation angles
    real(kind=8) :: rot_ang_errors(3) 

! A quaternion rotation is applied to the ACC data using the SCA star camera data
! There are gaps in the SCA quaternion date, which require an interpolation, and there are 
! half steps in the integration, which also require interpolation.
! Interpolation is a simple linear inerpolation (average of final and initial) 

    ACC_time = dble(jd)*86400.d0+tin-2451545*86400.d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   search for appropriate ACC time step, interpolate if necessary.
    ioerr2 = 0
    flag2 = -1

    i=0
    do while (ioerr2.eq.0 .and. i < ACCn)  ! PT140608: added this extra condition to stop it reading past the end of the 24 hr day
      i=i+1

      if (i.gt.ACCn) then ! if all elements of array have been searched
        print*, "ERROR: Time is out of ACC.dat file range"
        stop
      endif

      if (ACCtime(i).eq.ACC_time) then  ! if time step of ACC data is equal to current timestep
        t2out = ACCtime(i)
        do k=1,3
          accobs(k) = accr(k,i)
! PT130815: average the 5 x 1sec values to derive the acceleration value to use for this 5 sec step (instead of just using the first value)
!           NOTE: if there are epoch gaps in the accelerometer obs then this code won't work properly !!!!!
! PT170622: this doesn't work for MDC 5sec ACC data. Put a quick fix in here to stop the averaging in this case
          if((ACCtime(i+1)-ACCtime(i))<2)then  ! only do the averaging if there is less than 2 sec between observations
            do j=1,4
              accobs(k)=accobs(k) + accr(k,i+j)
            enddo
            accobs(k) = accobs(k)/5.d0
          endif
        enddo
        ioerr2 = 1 ! to flag the end of the while loop because accobs has been calculated
      endif

      if (ACCtime(i).gt.ACC_time) then  ! the time falls between ACC time steps
                                        ! so need to set "final" values of interpolation
        t2f = ACCtime(i)
        do k=1,3
          accf(k) = accr(k,i)
        enddo
        if (flag2.eq.-1) then ! if the data time of the first array entry is greater than the integration time (flag = -1)
          call status_update('FATAL','GRACEORB','accelerom',' ',"First entry of ACC file is greater than t",0)
        else ! if (flag2.eq.0) so this is not the first array entry
          t2i = ACCtime(i-1)
          do k=1,3
            acci(k) = accr(k,i-1)
          enddo
          t2i_real = dble(t2i)
          t2f_real = dble(t2f)
          call interpol(t2i_real, t2f_real, acci, accf, ACC_time, accobs, 3)
          ioerr2 = 1  ! to flag the end of the while loop because accobs has been calculated
        endif
      endif
      flag2=0 ! to flag that the first line of the array has been passed

    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   search for appropriate SCA time step, interpolate if necessary.

    ioerr1 = 0
    flag1 = -1

    i=0
    do while (ioerr1.eq.0)
      i=i+1
      if (i.gt.STARn) then ! if all elements of array have been searched
        print*, "ERROR: Time is out of SCA.dat file range"
        stop
      endif
!print*,i,scatime(i),acc_time
      if (SCAtime(i).eq.ACC_time) then  ! if time step of ACC data is equal to current timestep
        t1out = SCAtime(i)
        do k=1,4
          varout(k)=varr(k,i)
        enddo
        ioerr1 = 1 ! to flag the end of the while loop because accobs has been calculated
      endif

      if (SCAtime(i).gt.ACC_time) then  ! the time falls between ACC time steps
                                        ! so need to set "final" values of interpolation
        t1f = SCAtime(i)
        do k=1,4
          varf(k) = varr(k,i)
        enddo

        if (flag1.eq.-1) then ! if the data time of the first array entry is greater than the integration time (flag = -1)
          write(message,'(a,f15.7,a,f15.7,a)')"First entry of SCA file (",SCAtime(i),") is greater than start of integration (" &
                                           ,ACC_time,")" 
          call status_update('FATAL','GRACEORB','accelerom',' ',message,0)
        else ! if (flag1.eq.0) so this is not the first array entry

          t1i = SCAtime(i-1)
          do k=1,4
            vari(k) = varr(k,i-1)
          enddo

          t1i_real = dble(t1i)
          t1f_real = dble(t1f)
          call interpol(t1i_real, t1f_real, vari, varf, ACC_time, varout, 4)

          ioerr1 = 1  ! to flag the end of the while loop because varout has been calculated
        endif

      endif
      flag1=0 ! to flag that the first line of the array has been passed

    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    set SCA dummy variable to real variable
! PT130510: replace this with storage in a temporary array that can then be modified with roll/pitch/yaw corrections
    quat_tmp(0) = varout(1)
    quat_tmp(1) = varout(2)
    quat_tmp(2) = varout(3)
    quat_tmp(3) = varout(4)

! PT130510: in graceorb we converted the roll/pitch/yaw biases from the command line into a matrix, rotcorr. Convert it to a quaternion
!           and then add it on now to the quaternion from the SCA file.
! PT140818: try converting random noise in the pitch and yaw directions into a quaternion to perturb the orientation
!    rot_ang_errors(1) = 0.d0
!    rot_ang_errors(2) = (rand(0)-0.5d0)*0.1d-3   ! 0.1  milliradians in pitch
!    rot_ang_errors(3) = (rand(0)-0.5d0)*0.05d-4  ! 0.05 milliradians in yaw

! PT170404: uncomment this code - rpy perturbations are now passed in as either zero (for main orbit integration) or as 
!           perturbations (for calculating partials)
    rot_ang_errors(1:3) = rpy_perturb(1:3)
    call rpy2rotmat(rot_ang_errors(1),rot_ang_errors(2),rot_ang_errors(3) , rotcorr, rotcorr_deriv_roll, rotcorr_deriv_pitch &
                    , rotcorr_deriv_yaw)

    call rotmat2quat(rotcorr, quat_corr)
    call quatcomb(quat_corr, quat_tmp, quat)

!   set ACC dummy variable to real variable
    accrom(1,1) = accobs(1)
    accrom(2,1) = accobs(2)
    accrom(3,1) = accobs(3)
    
!   use bias and scale parameters to calculate actual acceleration
    do i = 1,3
! PT190218: original equation. This produces high correlations between bias and scale
!      accrom(i,1) = bias_perturb(i) + bs(i) + (scl(i)+scl_perturb(i))*accrom(i,1)

!DEBUG
!print*,jd,tin,i,bs(i),scl(i),bias_perturb(i),scl_perturb(i),accrom(i,1) &
!       ,(scl(i)+scl_perturb(i))*(accrom(i,1) + bs(i) - bias_perturb(i))

! PT190208: this was Tony Purcell's equation. The "+ bias" is to fix a sign problem. We'll deal with that later ....
      accrom(i,1) = (scl(i)+scl_perturb(i))*(accrom(i,1) + bs(i) + bias_perturb(i))
    enddo

    call rotation_quat2mat_3d(quat, Qrow)

!   rotate the accelerometer data from the Science Reference Frame (SRF) to the inertial frame 
    call matmult(Qrow,accrom,accrot,3,3,1)

! PT170403: I think the SCA1B quaternions actually rotate to earth-fixed space although they are supposed to rotate
!           into inertial space. To test this, we will assume that accrot is actually now in e-fixed. We then need
!           to rotate by rot_e2i to get the accelerations in inertial space  
! and so we have the calibrated, non-gravitational accelerations in the inertial XYZ directions.
    acc_if(1) = accrot(1,1)
    acc_if(2) = accrot(2,1)
    acc_if(3) = accrot(3,1)

! PT170404: it is WRONG to apply the e2i transformation. The SCA quaternions already put the accelerometer obs into the 
!           inertial reference frame
!    call matmult(rot_e2i,accrot(1:3,1),acc_if,3,3,1)
!    do i=1,3
!      print*,'e-fixed and inertial acc',tin,i,accrot(i,1)*1.d9,acc_if(i)*1.d9,(accrot(i,1)-acc_if(i))*1.d9
!    enddo
!    print*," "

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    return 
    end
