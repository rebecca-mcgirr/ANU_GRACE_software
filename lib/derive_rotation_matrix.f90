!pyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      subroutine derive_rotation_matrix(angle,iaxis,r)
!
! Written by Yehuda Bock
!
! this routine computes a rotation matrix
! input: angle  - rotation angle
!        iaxis  - rotation axis
! output : r    - rotation matrix
!
      implicit none

      integer*4 iaxis,j,k
      real*8 angle,r,sinang,cosang

      dimension r(3,3)
!
      j=MOD(iaxis,3) + 1
      k=MOD(j,3) + 1
      r(iaxis,iaxis)=1.d0
      r(iaxis,j)=0.d0
      r(j,iaxis)=0.d0
      r(iaxis,k)=0.d0
      r(k,iaxis)=0.d0
      sinang=dsin(angle)
      cosang=dcos(angle)
      r(j,j)=cosang
      r(k,k)=cosang
      r(j,k)=sinang
      r(k,j)=-sinang
!
      return
      end

