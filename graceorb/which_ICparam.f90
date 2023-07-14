  subroutine which_ICparam(prmcode,prm_input,sat,param_offset)

! subroutine to determine which is the appropriate value in the apriori vector from the input file in 
! order to update the ICs
!
! P. Tregoning
! 7 November 2013

  implicit none

  character*4,  intent(in)   :: prmcode               ! the ascii descriptor for the parameter
  character*30, intent(in)   :: prm_input(*)          ! the list of descriptors for the ICs in the input file
  character*1,  intent(in)   :: sat                   ! A or B
  integer*4,    intent(out)  :: param_offset          ! how far through the parameters we go before finding the one of interest

  integer*4  max_ICprm       ! maximum number of ICs (currently 15 per satellite) 
  integer*4  i,j

! PT140822: increase this  from 30 to 34 so that we can include 1/rev and 2/rev terms
  max_ICprm = 34
  param_offset = 0

! loop through the available parameters and see whether we can find the one that we're looking for
  do i=1,max_ICprm
    if(prm_input(i)(10:10) == sat .and. prm_input(i)(13:16) == prmcode ) then
      param_offset = i
      return
    endif
  enddo

  return
  end

