  module gtorb_mod
use orbits_rw
! module to define - in only one place - the record length of a GTORB file. This value will now be written to the 
! first line of the GTORB files and used by all programs to define the record length

! Line contains:  epoch (I4), pos/vel (6xR*8), quat(4xR*8), IC partials (36xR*8), bias (18xR*8), scale (18xR*8), 1/rev and 2/rev (4*R8), GPSant (3x6xR*8)
!                 , roll/pitch/yaw (3*6*R*8), 6x num_mascxR*8, 5x2 tidal amplitudes x 1000 tidal_masconsxR*8

! PT170609: increased the GTORB record length to permit up to 7550 mascons. All numbers are considered to be R*8
!           pos/vel/bs/scl     1/rev,2/rev,GPSant,rpy       msc         tmsc
!  integer*4, parameter :: GTORB_recl = 4 + (6+4+36+18+18)*8 +      (2+2+3+3)*6*8        + (7550)*6*8 + (1000*5*2)*6*8
! RM190305: increased the GTORB record length to permit up to 48000 mascons.
! PT190403: declare the variable here but give it a numerical value in graceorb once we know how many mascons there are
!  integer*4, parameter :: GTORB_recl = 4 + (6+4+36+18+18)*8 +      (2+2+3+3)*6*8        + (48000)*6*8 + (1000*5*2)*6*8
  integer*4 :: GTORB_recl

! declare a variable to be the record length read from the GTORB file(s)
  integer*4            :: file_recl

! define the adjustable variable of the particular record of the file
  integer*4            :: nrec

! declare a variable, being the GTORB record length prior to the creation of this file. 
!    This is the record length of all GTORB files from 27/08/2015 to 08/06/2017
  integer*4, parameter :: GTORB_recl_old = 318788

! PT191126: variable to indicate that we only want to integrate and output the main orbit, no partials
  character*1  :: GTORB_partials

type(orbits_data):: orb_file

  end module gtorb_mod
