  subroutine bubble_sort_R8(nvals,index_in,vec_in,index_out,vec_out)

! subroutine to sort a bunch of floating point numbers and output the sorted values and the order of the original indices
!
! P. Tregoning
! 19 April 2019

  implicit none

  integer*4   ,intent(in)    :: nvals              ! number of values to be sorted
  real(kind=8),intent(in)    :: index_in(nvals)    ! indices of order of pre-sorted values
  real(kind=8),intent(in)    :: vec_in(nvals)      ! pre-sorted values
  real(kind=8),intent(out)   :: index_out(nvals)   ! order of indices of sorted values
  real(kind=8),intent(out)   :: vec_out(nvals)     ! sorted values

! local variables
  integer*4  :: ival,bubble,j
  real(kind=8) :: temp,temp_index


! set up the starting values
  index_out = index_in
  vec_out = vec_in
  ival = nvals

! sort the values
  do while (ival > 1)
    bubble = 0 !bubble in the greatest element out of order
    do j = 1, (ival-1)
      if (vec_out(j) > vec_out(j+1)) then
        ! adjust the indices
        temp_index = index_out(j)
        index_out(j) = index_out(j+1)
        index_out(j+1) = temp_index
        ! adjust the values
        temp = vec_out(j)
        vec_out(j) = vec_out(j+1)
        vec_out(j+1) = temp
        bubble = j
      endif 
    enddo
    ival = bubble   
  enddo   

  return

end subroutine bubble_sort_R8


! #######################################################################
  subroutine rapid_sort(less_than,specific_value,nvals,index_in,vec_in,index_out,vec_out)

! subroutine to do something like a bubble sort but only for the values passed in that are less than/more than 
! a particular value. That is, don't bother sorting values that are outside some specified range.
!
! P. Tregoning
! 10 February 2020

  implicit none

  logical     ,intent(in)    :: less_than          ! flag to indicate whether values to sort must be more or less than passed value
  real(kind=8),intent(in)    :: specific_value     ! then specified min(max) value 
  integer*4   ,intent(in)    :: nvals              ! number of values to be sorted
  real(kind=8),intent(in)    :: index_in(nvals)    ! indices of order of pre-sorted values
  real(kind=8),intent(in)    :: vec_in(nvals)      ! pre-sorted values
  real(kind=8),intent(out)   :: index_out(nvals)   ! order of indices of sorted values
  real(kind=8),intent(out)   :: vec_out(nvals)     ! sorted values

! local variables
  integer*4  :: ival,bubble,j,use_vals
  real(kind=8) :: temp,temp_index

! duplicate arrays for the temporary variables
  real(kind=8),allocatable :: vec_in_tmp(:),vec_out_tmp(:),index_in_tmp(:),index_out_tmp(:)
  allocate(vec_in_tmp(nvals))
  allocate(vec_out_tmp(nvals))
  allocate(index_in_tmp(nvals))
  allocate(index_out_tmp(nvals))


! do a first pass to separate out the values to be sorted from the values we don't care about
  use_vals = 0
  do ival =1,nvals
    if( (less_than .and. vec_in(ival) < specific_value) .or. (.not. less_than .and. vec_in(ival) > specific_value) )then
      use_vals = use_vals + 1
      index_in_tmp(use_vals) = index_in(ival)
      vec_in_tmp(use_vals) = vec_in(ival)
    endif
  enddo

! set up the starting values
  index_out_tmp = index_in_tmp
  vec_out_tmp = vec_in_tmp
  ival = use_vals

! sort the values of only the entries that meet the highest/lowest criterion
  do while (ival > 1)
    bubble = 0 !bubble in the greatest element out of order
    do j = 1, (ival-1)
      if (vec_out_tmp(j) > vec_out_tmp(j+1)) then
        ! adjust the indices
        temp_index = index_out_tmp(j)
        index_out_tmp(j) = index_out_tmp(j+1)
        index_out_tmp(j+1) = temp_index
        ! adjust the values
        temp = vec_out_tmp(j)
        vec_out_tmp(j) = vec_out_tmp(j+1)
        vec_out_tmp(j+1) = temp
        bubble = j
      endif 
    enddo
    ival = bubble   
  enddo   

! now, translate vec_out_tmp and index_out_tmp to the returned matrices
  vec_out = vec_out_tmp
  index_out = index_out_tmp

  return
  end



