      program make_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! AUTHOR:  Tony Purcell
! DATE:    March 19 2013
! PURPOSE: To construct a gracefit-style input record with user-defined
!          mascon entries
! INPUT:   A series of mascon numbers and mascon heights
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

include '../includes/write_soln.h90'

real*8 :: apriori(maxparm), soln(maxparm)
real*8 :: vcv(maxparm, maxparm)

integer*4 ::
