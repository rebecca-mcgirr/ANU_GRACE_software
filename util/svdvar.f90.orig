	SUBROUTINE svdvar_dp(v,w,cvm,ma)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	INTEGER(I4B) :: ma
	REAL(DP), DIMENSION(ma,ma), INTENT(IN) :: v
	REAL(DP), DIMENSION(ma), INTENT(IN) :: w
	REAL(DP), DIMENSION(ma,ma), INTENT(OUT) :: cvm
	REAL(DP), DIMENSION(ma) :: wti
!	ma=assert_eq((/size(v,1),size(v,2),size(w),size(cvm,1),size(cvm,2)/),&
!		'svdvar')
	where (w /= 0.d0)
		wti=1.d0/(w*w)
	elsewhere
		wti=0.d0
	end where
	cvm=v*spread(wti,dim=1,ncopies=ma)
	cvm=matmul(cvm,transpose(v))
	END SUBROUTINE svdvar_dp
