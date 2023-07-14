  module coeff_mod
      integer maxsize
      parameter (maxsize = 360)
      real(kind=8), dimension(0:maxsize, 0:maxsize) :: meancoefC, meancoefS, varC, varS
      real(kind=8), dimension(0:maxsize, 0:maxsize) :: coefC, coefS
      real(kind=8), dimension((maxsize+1)*(maxsize+2)/2) :: dCdealias, dSdealias
      real(kind=8), dimension(0:maxsize, 0:maxsize) :: dCotide, dSotide
      real(kind=8), dimension(0:maxsize, 0:maxsize) :: dCatide, dSatide
      real(kind=8), dimension(0:maxsize, 0:maxsize) :: dCocepol, dSocepol
      real(kind=8), dimension(0:360, 0:360) :: AR, BR, AI, BI 
      real(kind=8) :: dCsolpol21, dSsolpol21
      integer :: flagrun 
  end module coeff_mod

