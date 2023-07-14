    subroutine reconstruct(xx,yy, otidecoef, func)

! Emma-Kate Potter, 13 Oct, 2010
! This subroutine
! 1) calls for the calculation of a grid of data that represents ocean height at current epoch
! 2) calls another subroutine to convert the grid into spherical harmonic coefficients
! 3) scales the coefficients to be used in calculation of potential due to that ocean height

    use gm_mod
    use gauleg_mod
    implicit none
    integer*4 :: m,n 
    integer*4 :: j
    real*8, dimension(0:nmax, 0:nmax) :: cfC, cfS
    real*8, dimension(maxnml) :: otidecoef
    real*8 :: legnm(0:nmax, 0:nmax)
    real*8 :: pnm(maxnm)
    real*8 :: func, xx, yy


    j=1
    do n = 0,nmax
      do m = 0, n
        if (m.ne.0) then
          cfC(n,m) = otidecoef(j)
          cfS(n,m) = otidecoef(j+1)
!print*, n,m,cfC(n,m), cfS(n,m), otidecoef(j), otidecoef(j+1)
          j=j+2
        else
          cfC(n,m) = otidecoef(j)
          cfS(n,m) = 0.d0
!print*, n,m,cfC(n,m), cfS(n,m)
          j=j+1
        endif
      enddo
    enddo
!stop
!    do n = 0,10
!      do m = 0, n
! print*, n, m, sngl(cfC(n,m)), sngl(cfS(n,m))
!      enddo
!    enddo

    pnm = 0.d0
    call calpnm(xx, pnm)

    j = 1
    do n=0,nmax
      do m=0,n
        legnm(n,m) = pnm(j)
        j=j+1
      enddo
    enddo

    func = 0.d0
    do n=0,nmax 
!    do n=0,10 
! for m=0 term:
!        print*, legnm(n,0), cfC(n,0)
        func = func + legnm(n,0)*(cfC(n,0))
      do m=1,n
        func = func + 2.d0*legnm(n,m)*(cfC(n,m)*dcos(dble(m)*yy)+cfS(n,m)*dsin(dble(m)*yy))
      enddo
    enddo
return
end



