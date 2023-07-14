  module GPSant_mod
    real(kind=8), dimension(3) :: GPSant_cosL1, GPSant_cosL2
    real(kind=8)  :: GPSant_magL1, GPSant_magL2
    real(kind=8), dimension(3) :: GPSant_L1, GPSant_L2
! PT130528: offset corrections to the GPS antenna vector 
    real(kind=8), dimension(3) :: gpsantoff
  end module GPSant_mod

