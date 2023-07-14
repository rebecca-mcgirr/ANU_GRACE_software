    subroutine efix_inert(rot, rotdot, rotacc, incoor, outcoor, inacc, acc_flag)

!  E-K Potter May 2010
!  takes the inertial and calculates the ef pos
!
! MOD PT140318: pass in the rotation and rotation rate matrices, then compute the transformation (irrespective of direction)

    implicit none

    integer*4 :: i,idir
    real(kind=8), dimension(3,3), intent(in ) :: rot            ! rotation to/from inertial-efixed. Can be for either direction.
    real(kind=8), dimension(3,3), intent(in ) :: rotdot         ! rotation rate to/from inertial-efixed. Can be for either direction.
    real(kind=8), dimension(3,3), intent(in ) :: rotacc         ! rotation acceleration to/from inertial-efixed. Can be for either direction.
    real(kind=8), dimension(3)  , intent(in)  :: inacc          ! input accelerations
    real(kind=8), dimension(6),   intent(in ) :: incoor         ! input coordinates (XYZ, Xdot Ydot Zdot)
    real(kind=8), dimension(6),   intent(out) :: outcoor        ! transformed coordinates (XYZ, Xdot Ydot Zdot)
    logical                   ,   intent(in)  :: acc_flag       ! flag whether to compute accelerations or not

    integer*4 :: jd
    real(kind=8) :: t
    character*3 :: sflag
    real*8 :: PEPt
    integer*4 :: PEPjd
    real(kind=8) :: xpole, ypole

    do i=1,3
      outcoor(i) = rot(i,1)*incoor(1)+rot(i,2)*incoor(2)+rot(i,3)*incoor(3)
    enddo
    do i=4,6
      outcoor(i)=rot(i-3,1)*incoor(4)+rot(i-3,2)*incoor(5)+rot(i-3,3)*incoor(6)+ &
               rotdot(i-3,1)*incoor(1)+rotdot(i-3,2)*incoor(2)+rotdot(i-3,3)*incoor(3)
    enddo

! PT170821: include acceleration transformations if requested
    if(acc_flag) then
! X_acc = R_acc X + 2 * R_vel V + R_pos A
      do i=1,3
        outcoor(i) = rotacc(i,1)*incoor(i)+rotacc(i,2)*incoor(2)+rotacc(i,3)*incoor(3) &
                   + 2.d0*(rotdot(i,1)*incoor(4)+rotdot(i,2)*incoor(5)+rotdot(i,3)*incoor(6)) &
                   + rot(i,1)*inacc(1)+rot(i,2)*inacc(2)+rot(i,3)*inacc(3)
      enddo
    endif

    return
    end
