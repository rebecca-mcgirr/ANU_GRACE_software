   subroutine biasscale(jd, tin, c0x,c0y,c0z,c1x,c1y,c1z,c2x,c2y,c2z,scl )

! Emma-Kate Potter, 18 August, 2010
! This subroutine calculates the accelerometer bias using the input coefficients 

! original values (in comments) from "Recommendation for a-priori Bias and Scale Parameters 
! for Level 1B ACC Data (Version 2)" by Srinivas Bettadpur, Jun 17 2009
! fnew = bias + scale*facc1b
! bias = c0 + c1 (Td-T0) + c2(Td-T0)^2

! expect bias coefficient "c0" to be adjusted using GRACEFIT

    use sat_mod
! PT140813: turned this off - use values passed through instead
!    use bsscl_mod
    use inmod_mod
    use accel_mod

    implicit none

    real(kind=8),intent(in)  :: c0x(2), c0y(2), c0z(2), c1x(2), c1y(2), c1z(2), c2x(2), c2y(2), c2z(2)   ! bias information
    real(kind=8),intent(in)  :: scl(3)           ! accelerometer scales

    integer :: i
    integer*4 :: jd, time, dt
    real(kind=8) :: tin
!    real(kind=8), dimension(3) :: bs 

    time= aint(jd+tin/86400.d0-2400000.5d0)

    if (time.le.52705) then ! cut off date is March 7, 2003
      dt = time-52532
      bs(1) = c0x(1) + c1x(1)*dble(dt) + c2x(1)*dble(dt*dt)
      bs(2) = c0y(1) + c1y(1)*dble(dt) + c2y(1)*dble(dt*dt)
      bs(3) = c0z(1) + c1z(1)*dble(dt) + c2z(1)*dble(dt*dt)
    else
      dt = time-53736
      bs(1) = c0x(2) + c1x(2)*dble(dt) + c2x(2)*dble(dt*dt)
      bs(2) = c0y(2) + c1y(2)*dble(dt) + c2y(2)*dble(dt*dt)
      bs(3) = c0z(2) + c1z(2)*dble(dt) + c2z(2)*dble(dt*dt)
    endif

    do i=1,3
      bs(i) = bs(i)*1.d-6
    enddo
!    print *, c0x,c0y,c0z,c1x,c1y,c1z,c2x,c2y,c2z,scl
!print *, time, dt 
!    print *, bs
!    STOP -2
    return 
    end
