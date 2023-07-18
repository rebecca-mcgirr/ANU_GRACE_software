    module gm_mod       
!   E-K Potter 20/05/10
!   This module declares variables for the gravitational constants 
!   of the earth, moon and sun

    real (kind=8), dimension(3) :: gm
    real (kind=8), dimension(11) :: gmpln
    real (kind=8) :: Gconst 
    real (kind=8) :: Ae
    real (kind=8) :: pi
    real (kind=8) :: Me 
    real (kind=8) :: rhow
    real (kind=8) :: Om  
    real (kind=8) :: rho_av=5515.d0      ! average density of the Earth (needed to convert from EWH to Stokes' coefficients)

! WGS84 ellipsoid parameters
    real(kind=8), parameter :: earthrad_wgs84 = 6378137.d0
    real(kind=8), parameter :: earthflat_wgs84 = 1.d0/298.257222101d0
    
    end module
