  module dealias_mod
      integer max_num_epoch
      parameter ( max_num_epoch = 10 )
      real(kind=8), dimension(max_num_epoch) :: epoch_time  
!      real(kind=8), dimension(0:100, 0:100, max_num_epoch) :: CT, ST
! PT190528: RL06 AOD1B is degree 180
      integer*4,parameter  :: aod1b_degree = 180
      real(kind=8), dimension(max_num_epoch, (aod1b_degree+1)*(aod1b_degree+2)/2) :: CT, ST
!      integer dealias_counter, rnge_old
      integer num_dealias_epochs
  end module dealias_mod

