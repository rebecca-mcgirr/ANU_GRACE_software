  subroutine which_MSCtide_param(nparam,tidecode,prm_input,imsc,param_offset)

! subroutine to determine which is the appropriate value in the apriori vector from the input file in 
! order to update the ICs
!
! P. Tregoning
! 7 November 2013

  use mascon_mod   ! needed for the mcon_ocean_prim declaration

  implicit none

  integer*4,    intent(in)   :: nparam                ! number of parameter a priori values found in input VCV file           
  character*2,  intent(in)   :: tidecode              ! the ascii descriptor for the requested tide (M2, O1, S2, K1, K2) parameter
  character*30, intent(in)   :: prm_input(*)          ! the list of descriptors for the parameters in the input file
  integer*4,    intent(in)   :: imsc                  ! the mascon number
  integer*4,    intent(out)  :: param_offset          ! how far through the parameters we go before finding the one of interest

  character*4                :: char_msc 
  integer*4  i,j

  param_offset = 0

  write(char_msc,'(i4.4)')imsc
! loop through the available parameters and see whether we can find the one that we're looking for
! PT160321: fix bug in the loop counter maximum: there are "max_msc_tides" tides per mascon
  do i=1, nparam
    if(prm_input(i)(8:10) == "TMC" .and. prm_input(i)(16:17) == tidecode  .and. prm_input(i)(11:14) == char_msc) then
      param_offset =  i
      return
    endif
  enddo

  return
  end

