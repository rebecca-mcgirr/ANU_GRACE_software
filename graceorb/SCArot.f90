   subroutine SCArot(jd, tin, SCArot_mat, quat)

! Emma-Kate Potter, 10 May, 2011
! This subroutine extracts the star camera data for this time step and returns a 
! rotation matrix that rotates from the science reference frame to the inertial frame
! There are gaps in the SCA quaternion date, which require an interpolation. 
! Interpolation is a simple two point interpolation.
!
! MODS:
! PT/APP 130710: use a consistent sign convention for the quaternions. We no longer change the sign of that read from SCA, and use the quaternion routines in lib/quat_lib.f90

!    use accel_mod
    use accred_mod ! stores the SCA data 
!    use bsscl_mod

    implicit none

    character*300 :: message
    integer*4 :: i, j, jd, ioerr1
    integer*4 :: t1i, t1f, t1tmp, t1out
    integer*4 :: SCAnum 
    integer :: k
    real(kind=8) ::  tin, ACC_time 
    real(kind=8) ::  t1f_real, t1i_real 
    real(kind=8), dimension(4) ::  vari
    real(kind=8), dimension(4) ::  varf
    real(kind=8), dimension(4) ::  varout
    real(kind=8), dimension(4) ::  vartmp
    integer :: flag1 
    real(kind=8) ::  w, x, y, z
     real(kind=8), dimension(0:3) ::  quat
     real(kind=8), dimension(0:3) ::  quattest
    real(kind=8), dimension(3,3) :: SCArot_mat

    ACC_time = dble(jd)*86400.d0+tin-2451545*86400.d0
    SCAnum=60
!*******************************

!   search for correct SCA data timestep

    ioerr1 = 0
    flag1 = -1

    i=0
    do while (ioerr1.eq.0)
      i=i+1
      if (i.gt.STARn) then ! if all elements of array have been searched
        print*, "ERROR: Time is out of SCA.dat file range"
        stop
      endif

      if (SCAtime(i).eq.ACC_time) then  ! if time step of ACC data is equal to current timestep
        t1out = SCAtime(i)
        do k=1,4
          varout(k)=varr(k,i)
        enddo
        ioerr1 = 1 ! to flag the end of the while loop because accout has been calculated
      endif

      if (SCAtime(i).gt.ACC_time) then  ! the time falls between ACC time steps
                                        ! so need to set "final" values of interpolation
        t1f = SCAtime(i)
        do k=1,4
          varf(k) = varr(k,i)
        enddo
        if (flag1.eq.-1) then ! if the data time of the first array entry is greater than the integration time (flag = -1)
          write(message,'(a,f20.7,a,f20.7,a)')"First entry of SCA file (",SCAtime(i),") is greater than start of integration (" &
                                           ,ACC_time,")" 
          call status_update('FATAL','GRACEORB','SCArot',' ',message,0)
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

!   set SCA dummy variable to real variable
    quat(0) = varout(1)
    quat(1) = varout(2)
    quat(2) = varout(3)
    quat(3) = varout(4)

!   print*, "ACCfile", sqrt(accx*accx+accy*accy+accz*accz)
!    accrom(1,1) = accout(1)
!    accrom(2,1) = accout(2)
!    accrom(3,1) = accout(3)
!    accrom(4,1) = 0.d0
    
!    call biasscale(jd, tin, bs)

!    do i = 1,3
!      accrom(i,1) = bs(i)+scl(i)*accrom(i,1)
!    enddo

    call rotation_quat2mat_3d(quat, SCArot_mat)

!test rotmat2quat subroutine
!     call rotmat2quat(SCArot_mat, quattest)
!     do j=0,3
!       print*, quat(j), quattest(j)
!     enddo
!end test

    return 
    end
