  module write_soln_mod
      include '../includes/grace.param'
      include '../includes/write_soln.h90'
! PT1902312: make these arrays allocatable
      real(kind=8), dimension(maxparm) :: soln, apriori  
      real(kind=8), dimension(maxparm, maxparm) :: VCV_obs
  end module write_soln_mod

