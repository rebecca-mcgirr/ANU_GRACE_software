    subroutine quatcomb(quatQ, quatP, quatS)

! E-K Potter May 2011
! This subroutine combines 2 sets of quaternion values 
! The resulting quaternion represents operation on vector V in frame A first by
! quaternion Q (from frame A to B) then by quaternion P (from frame B to C)

    implicit none

    real (kind=8), dimension(0:3) :: quatQ, quatP, quatS 
    
    quatS(0) = quatQ(0)*quatP(0) - quatQ(1)*quatP(1) - quatQ(2)*quatP(2) - quatQ(3)*quatP(3)
    quatS(1) = quatQ(1)*quatP(0) + quatQ(0)*quatP(1) - quatQ(3)*quatP(2) + quatQ(2)*quatP(3)
    quatS(2) = quatQ(2)*quatP(0) + quatQ(3)*quatP(1) + quatQ(0)*quatP(2) - quatQ(1)*quatP(3)
    quatS(3) = quatQ(3)*quatP(0) - quatQ(2)*quatP(1) + quatQ(1)*quatP(2) + quatQ(0)*quatP(3)

    return
    end

