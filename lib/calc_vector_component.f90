  subroutine calc_vector_component(vect1,vect2,value)

! subroutine to calculate the component of one vector in the direction of another
!
! P. Tregoning
! 11 November 2015

  implicit none

! passed variables
  real(kind=8), intent(in)  :: vect1(3),vect2(3)       ! input vectors. We want the component of vect1 in direction of vect2
  real(kind=8), intent(out) :: value                   ! output value, being the magnitude of the component of vector 1 in direction of vector 2

! local variables
  integer*4    :: i
  real(kind=8) :: dot,amag3
  real(kind=8) :: tmpvect1(3),tmpvect2(3)


! bail out if vect1 has no magnitude
  if(amag3(vect1) < 1.d-15)then
    value = 0.0
    return
  endif

! normalise the input variables
  do i=1,3
    tmpvect1(i) = vect1(i)/amag3(vect1)
    tmpvect2(i) = vect2(i)/amag3(vect2)
  enddo


! now, the component of vect1 in the direction of vect2 is the magnitude of vect1 times the dot product of the two unit vectors
  value = amag3(vect1) * dot(tmpvect1,tmpvect2)

! set really small values to zero
  if(abs(value) < 1.e-14)value = 0.0

  return
  end



