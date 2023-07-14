subroutine compute_AORC(data)

  ! Tony Purcell
  ! April 12 2013
  ! Routine to take an array of data, calculate its derivative w.r.t time
  ! and store the result in the second leaf of the matrix
  use gracefit_mod
  implicit none

   

  double precision, dimension(9, 3, nepochs_t) :: data

  double precision, dimension(nepochs_t) :: vector, vectordot
  double precision :: delta

  integer :: j, k, n, count,i,jj

  ! find the length of the data records (9th column contains epoch value
  count = 1
  do while ( data(9, 1, count) .ne. 0.0d0 .and. count < nepochs_t)
     count = count + 1
  enddo

  ! PT140501: drop this down to using just three points (it was set to 5)
  delta = 5.0d0
  n = 11
  do j = 1, 8
     vector = 0.0d0
     vectordot = 0.0d0

     ! copy out j-th colum of data
     do k = 1, count
        vector(k) = data(j, 1, k) 
     enddo

     ! calculate derivative
     call noise_robust_deriv(vector, vectordot, delta, count, n)

     ! copy calculated derivative into second leaf of input array
     do k = 1, count
        data(j, 2, k) = vectordot(k)
     enddo
  enddo

  return
end subroutine compute_AORC
