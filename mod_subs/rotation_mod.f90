    module rotation_mod

    integer*4 :: iut1pol,iut1,ipole,inut,iscrn
    real*8 :: tdtgpst
    real*8, dimension(3,3) :: rot,rotdot, rotcorr
    real*8, dimension(3,3) :: rotcorr_deriv_roll, rotcorr_deriv_pitch, rotcorr_deriv_yaw
    real*8  :: sidtm
    character*5 :: frame,precmod     
! PT140224: add the rotation matrix for frame/precession/nutation as computed by SOFA routines. Allow 17280 epochs (=1 day @ 5 sec intervals)
    real(kind=8) :: rot_nutprec(3,3,0:17280),dpsi(0:17280),deps(0:17280),nut_epoch(0:17280)

! PT140225: add xp,yp,ut1utc,mjdates, being values read from usno.finals.data
    real(kind=8),dimension(10000) :: xp,yp,ut1utc,mjdates
    integer*4 :: iusno,nvals_eop

! PT170403: move the earth-fixed to inertial (and vice versa) rotation matrices into this mod file
    real(kind=8),dimension(3,3) :: rot_i2e,rot_e2i,rotdot_i2e,rotdot_e2i,rotacc_i2e,rotacc_e2i

    end module rotation_mod

