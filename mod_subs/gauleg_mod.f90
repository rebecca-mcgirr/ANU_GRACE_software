  module gauleg_mod
      integer(kind=4) :: maxnml, maxnm
      integer(kind=4) ::nmax, nlon, nlat, nlond2, nlatd2 
      parameter(nmax=256,nlon=2*nmax,nlat=nmax)
      parameter(maxnm=(nmax+1)*(nmax+2)/2)
      parameter(nlond2=nlon/2,nlatd2=nlat/2)
      parameter(maxnml=(nmax+1)*(nmax+1))
      real(kind=8) :: absc(nlat), wei(nlat)
      real(kind=8) :: pnmx(maxnm,nlatd2) 
  end module gauleg_mod

