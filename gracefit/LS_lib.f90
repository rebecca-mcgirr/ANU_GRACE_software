!********************************************************************************************************************************
!  File: LS_lib.f90
!
!  Purpose: Set of subroutines dealing with calculations for the Least Squares algorithm.
!
!  Author: Thomas Greenspan
!
!  API:
!       LS_norminc           : Increments the normal equations
!       LS_normSolve         : Solves the normal equations
!       LS_norminc_acc_shad  : Increments only the accelerometer shadow conditional equations
!
!   July 29, 2013
!
!********************************************************************************************************************************
!********************************************************************************************************************************
! LS_norminc: Increments normal equations (of one epoch) given the part and pre_omc matrices
!             (of specific epoch), the weights for the data, and the left and right-hand
!             sides of the Least Squares equation, normeq and AtWb respectively
! Author: Unknown
!********************************************************************************************************************************

subroutine LS_norminc_acc_shad(nobs_t,nparam_t,part,pre_omc,apr_wght,iScale,iBias,icond,iICs,normeq,AtWb)

  implicit none


  !********************  Variable declarations ********************

  integer         , intent(in)  :: nobs_t,nparam_t           ! number of observations and parameters
  double precision, intent(in)  :: pre_omc(nobs_t)           ! Observed Minus Computed values to be used in LS algorithm
  double precision, intent(in)  :: part(nobs_t,nparam_t)     ! Partials of observables with respect to the parameters
  double precision, intent(in)  :: apr_wght(nobs_t)          ! Weights to be applied to data
  integer         , intent(in)  :: iScale                    ! index of start of GRACE A scale parameters
  integer         , intent(in)  :: iBias                     ! index of start of GRACE A bias parameters
  integer         , intent(in)  :: icond                     ! index of start of conditional acceleration shadow obs
  integer         , intent(in)  :: iICs                      ! number of orbit ICs per satellite 
  real(kind=8 ), intent(out) :: normeq(nparam_t,nparam_t) ! The left side of the LS algorithm. P_t*W*P
  real(kind=8 ), intent(out) :: AtWb(nparam_t)            ! The right side of the LS algorithm. P_t*W*OMC

  integer*4 :: i,j,k,l   ! Counter variables
  character*100 message  ! character string for status_update calls
  double precision :: tempsum
  !****************************************************************

  !*********************  INCREMENT EQUATIONS *********************
  ! Left hand side
  ! GRACE A
  do j=iscale,iscale+5
     do k=iscale,iscale+5
        do i=icond,icond+1
           normeq(j,k) = normeq(j,k)  +   part(i,j)*part(i,k)*apr_wght(i)
        enddo
     enddo
  enddo

  ! GRACE B
  do j=iscale+iICs,iscale+iICs+5
     do k=iscale+iICs,iscale+iICs+5
        do i=icond+2,icond+3
           normeq(j,k) = normeq(j,k)  +   part(i,j)*part(i,k)*apr_wght(i)
        enddo
     enddo
  enddo

  ! Increment the right-hand side
  ! GRACE A
  do j=iscale,iscale+5
     do i=icond,icond+1
        AtWb(j) = AtWb(j) + part(i,j)*pre_omc(i) * apr_wght(i)
     enddo
  enddo

  ! GRACE B
  do j=iscale+iICs,iscale+iICs+5
     do i=icond,icond+1
        AtWb(j) = AtWb(j) + part(i,j)*pre_omc(i) * apr_wght(i)
     enddo
  enddo


  !****************************************************************

  return
end subroutine LS_norminc_acc_shad


!********************************************************************************************************************************
subroutine LS_norminc_Atb(nobs_t,nparam_t,part,pre_omc,apr_wght,normeq,AtWb)

   implicit none
 
 
   !********************  Variable declarations ********************
 
   integer         , intent(in)  :: nobs_t,nparam_t           ! number of observations and parameters
   double precision, intent(in)  :: pre_omc(nobs_t)           ! Observed Minus Computed values to be used in LS algorithm
   double precision, intent(in)  :: part(nobs_t,nparam_t)     ! Partials of observables with respect to the parameters
   double precision, intent(in)  :: apr_wght(nobs_t)          ! Weights to be applied to data
   real(kind=8 ), intent(out) :: normeq(nparam_t,nparam_t) ! The left side of the LS algorithm. P_t*W*P
   real(kind=8 ), intent(out) :: AtWb(nparam_t)            ! The right side of the LS algorithm. P_t*W*OMC
 
   integer*4 :: i,j,k,l   ! Counter variables
   character*100 message  ! character string for status_update calls
   double precision :: tempsum
   !****************************************************************
 
   !*********************  INCREMENT EQUATIONS *********************
   ! SA: use a "tempsum" variable to speed up the computations
   ! ! Increment the left-hand side
   ! do j = 1, nparam_t
   !    do k = j, nparam_t
   !       tempsum =0.0
   !       do i = 1, nobs_t 
   !          tempsum = tempsum + part(i,j) * part(i,k) * apr_wght(i)
   !       enddo
   !       normeq(k,j) = tempsum 
   !       normeq(j,k) = tempsum 
   !    enddo
   ! enddo
 
   ! Increment the right-hand side
   do j=1, nparam_t
      do i=1, nobs_t
         AtWb(j) = AtWb(j) + part(i,j)*pre_omc(i) * apr_wght(i)
      enddo
   enddo
   !****************************************************************
 
   return
 end subroutine LS_norminc_Atb



!********************************************************************************************************************************
subroutine LS_norminc(nobs_t,nparam_t,part,pre_omc,apr_wght,normeq,AtWb)

  implicit none


  !********************  Variable declarations ********************

  integer         , intent(in)  :: nobs_t,nparam_t           ! number of observations and parameters
  double precision, intent(in)  :: pre_omc(nobs_t)           ! Observed Minus Computed values to be used in LS algorithm
  double precision, intent(in)  :: part(nobs_t,nparam_t)     ! Partials of observables with respect to the parameters
  double precision, intent(in)  :: apr_wght(nobs_t)          ! Weights to be applied to data
  real(kind=8 ), intent(out) :: normeq(nparam_t,nparam_t) ! The left side of the LS algorithm. P_t*W*P
  real(kind=8 ), intent(out) :: AtWb(nparam_t)            ! The right side of the LS algorithm. P_t*W*OMC

  integer*4 :: i,j,k,l   ! Counter variables
  character*100 message  ! character string for status_update calls
  double precision :: tempsum
  !****************************************************************

  !*********************  INCREMENT EQUATIONS *********************
  ! SA: use a "tempsum" variable to speed up the computations
  ! Increment the left-hand side
  !do j = 1, nparam_t
  do k = 1, nparam_t
  !   do k = j, nparam_t
     do j = 1,k! j, nparam_t
        tempsum =0.0
        do i = 1, nobs_t 
           tempsum = tempsum + part(i,j) * part(i,k) * apr_wght(i)
        enddo
        normeq(k,j) = tempsum 
        normeq(j,k) = tempsum 
     enddo
  enddo

  ! Increment the right-hand side
  do j=1, nparam_t
     do i=1, nobs_t
        AtWb(j) = AtWb(j) + part(i,j)*pre_omc(i) * apr_wght(i)
     enddo
  enddo
  !****************************************************************

  return
end subroutine LS_norminc


!********************************************************************************************************************************


subroutine Chol_normSolve(nparam_t,normeq, AtWb, adjust)
   implicit none
   !********************  Variable declarations ********************
   integer :: i 
   integer*4       , intent(in)    :: nparam_t                  ! dimensioning of the number of parameters
   real(kind=8 ), intent(inout) :: normeq(nparam_t,nparam_t) ! The left side of the LS algorithm. P_t*W*P
   double precision, intent(out)   :: adjust(nparam_t)          ! Solution to normal equations
   real(kind=8 ), intent(in)    :: AtWb(nparam_t)            ! The right side of the LS algorithm. P_t*W*OMC
   integer :: ioerr
   !****************************************************************
 
   ! call status_update('STATUS','GRACEFIT','grace/LS_normSolve',' ','Solving the normal equations',0)
 
   !****************************************************************
   adjust = AtWb

   
  call dposv ( "U", nparam_t, 1 ,normeq, nparam_t, adjust, nparam_t, ioerr )
  call dpotri ( "U", nparam_t, normeq, nparam_t, ioerr )
  ! SA Fill the lower triangular part of the VCV matrix with the Upper part
  do i=1,nparam_t
        normeq(i,1:i) = normeq(1:i,i)
  enddo
end subroutine Chol_normSolve
 
!********************************************************************************************************************************

subroutine LU_normSolve(nparam_t,normeq, AtWb, adjust)
   implicit none
   !********************  Variable declarations ********************
 
   integer*4       , intent(in)    :: nparam_t                  ! dimensioning of the number of parameters
   real(kind=8 ), intent(inout) :: normeq(nparam_t,nparam_t) ! The left side of the LS algorithm. P_t*W*P
   double precision, intent(out)   :: adjust(nparam_t)          ! Solution to normal equations
   real(kind=8 ), intent(in)    :: AtWb(nparam_t)            ! The right side of the LS algorithm. P_t*W*OMC
   integer :: ioerr
   integer , allocatable :: IPV(:)
   double precision, allocatable::  WORK(:)

   allocate( IPV(nparam_t) )
   allocate( WORK(nparam_t) )
   adjust = AtWb
   call dgesv (nparam_t, 1 ,normeq, nparam_t, IPV, adjust, nparam_t, ioerr)
   call dgetri (nparam_t, normeq, nparam_t, IPV, WORK, nparam_t , ioerr)

end subroutine LU_normSolve
 
!********************************************************************************************************************************

