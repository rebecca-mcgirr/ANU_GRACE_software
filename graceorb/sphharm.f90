    subroutine sphharm(xload,xharm)
 
!c   calculate the spherical harmonic expansion of the function xload
!c   and write it to the array xharm

    use gm_mod
    use gauleg_mod

    implicit none
    real(kind=8) :: xh(maxnml)
    real(kind=8) :: xload(nlon, nlat), xcop(nlon, nlat), xharm(maxnml)
    real(kind=8) :: fac
    integer*4 :: icomp, isign, i, j, n, m, nm,nml


      fac = 1.d0/(2.d0 * dble(nlon))
!c make a copy of the original function values so that they are not lost during
!c the FFT
      do i = 1, nlat
        do j = 1, nlon
          xcop(j, i) = xload(j, i)
        enddo

!c calculate the FFT for the current longitude. This effectively calculates
!c the trigonometric terms in the spherical harmonic series
        call realft(xcop(1, i), nlond2, 1)
      enddo
      
!c  now xcop(j,i) {i=1,2,...,nlat} is the Fourier transform of xload at 
!c  latitude absc(i) in the following form:
!c
!c   F0, Fnmax, Re(F1), Im(F1), Re(F2), Im(F2), ..., Re(Fnmax-1), Im(Fnmax-1)
!c
!c  Im(Fnmax)=0 and is not stored.
!c
!c  Fm = sum_{k = 0}^{2 * nmax - 1} fk exp(2 pi i k m/(2 * nmax)
!c
!c  where fk is the input function.

!c  Now use the quadrature formula for the latitude integration
!c
!c  int pnm(x) F(x,phi) cos or sin m phi
!c
!c  p00, p10, p11 cos phi, p11 sin phi, p20, p21 cos phi, p21 sin phi,
!c  p22 cos 2phi, p22 sin 2phi, ...

!c   factor 2pi/2nmax * 1/4/pi from the Fourier transform and normalisation
!c   on the sphere

!c initialise function values
      do 400 i = 1, maxnml
400     xh(i) = 0.d0


!c  loop over the latitudes - first the southern half
      do 500 i = 1, nlatd2
        icomp = nlatd2 + 1 - i


!c  loop over spherical harmonic degree
        do n = 0, nmax

!  print*,'southern hem i,n',i,n

!c isign represents the change in sign between the northern and southern abscissa
          isign = 1 - mod(n, 2) * 2

          nml = n * n + 1
          nm = n * (n + 1)/2 + 1
          xh(nml) = xh(nml) + isign * fac * wei(i) * pnmx(nm, icomp) *xcop(1, i)
          do m = 1, n
            isign = 1 - mod(n + m, 2) * 2
            nm = nm + 1
            nml = nml + 1
            if ( nm .lt. maxnm ) then
              xh(nml) = xh(nml) + isign * fac * wei(i) * pnmx(nm, icomp)* xcop(2 * m + 1, i)
              nml = nml + 1
              xh(nml) = xh(nml) + isign * fac * wei(i) * pnmx(nm, icomp)* xcop(2 * m + 2, i)
            else

!c  the cosine Fourier component for nmax is stored at xcop(2,i)
!c  the sine Fourier component for nmax is zero

              xh(nml) = xh(nml) + isign * fac * wei(i) * pnmx(nm, icomp)* xcop(2, i)
            endif
          enddo
        enddo
        xh(nml) = xh(nml) + isign * fac * wei(i) * pnmx(nm, icomp)* xcop(2, i)
500     continue

!c  loop over the latitudes - northern hemisphere
      do 600 i = nlatd2 + 1, nlat
        icomp = i - nlatd2


!c  loop over spherical harmonic degree
        do n=0,nmax
!  print*,'northern hem i,n',i,n


          nml = n * n + 1
          nm = n * (n + 1)/2 + 1
          xh(nml) = xh(nml) + fac * wei(i) * pnmx(nm, icomp) *xcop(1, i)
          do m = 1, n
            nm = nm + 1
            nml = nml + 1
            if ( nm .lt. maxnm ) then
              xh(nml) = xh(nml) + fac * wei(i) * pnmx(nm, icomp) *xcop(2 * m + 1, i)
              nml = nml + 1
              xh(nml) = xh(nml) + fac * wei(i) * pnmx(nm, icomp) *xcop(2 * m + 2, i)
            else

!c  the cosine Fourier component for nmax is stored at xcop(2,i)
!c  the sine Fourier component for nmax is zero

              xh(nml) = xh(nml) + fac * wei(i) * pnmx(nm, icomp) *xcop(2, i)
            endif
          enddo
        enddo
600     continue
  
!c copy the contents of the dummy variable xh into the output varianle xharm
      do i = 1, maxnml
        xharm(i) = xh(i)
!  print*,'end of sphharm xh(i)',i,xh(i)
      enddo
      return
      end
