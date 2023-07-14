  subroutine infill_TS(debug,calling_prog,model,data_type,nepochs,ncols,values_in,epoch_step)

! subroutine to infill an array of quaternions for GRACE A and GRACE B so that each epoch has a quaternion value.
! Use a simple linear interpolation at this stage. If the start/end values are mission then extrapolate linearly
! using the subsequent/previous two values.
!
! P. Tregoning
! 27 March 2017
! 
! MODS:
! PT180620: this subroutine doesn't ever seem to have been used. I will take it over now and use it to infill gaps
!           in a time series using either linear or quadratic functions. Substantial mods will be made ....

  implicit none

  logical    ,  intent(in)    :: debug
  character*8,  intent(in)    :: calling_prog                ! name of calling program
  character*10, intent(in)    :: model                       ! "quadratic" or ... code it yourself !
  character*3 , intent(in)    :: data_type                   ! type of data being infilled (ACC, SCA, etc)
  integer*4,    intent(in)    :: nepochs                     ! number of epochs  of array to be infilled
  integer*4,    intent(in)    :: ncols                       ! number of columns of array to be infilled
  real(kind=8), intent(inout) :: values_in(nepochs,ncols+1)  ! array of values to be infilled. First column is epoch, rest are data.
  real(kind=8), intent(inout) :: epoch_step                  ! interval of time between epochs in input data

  real(kind=8), allocatable   :: values_tmp(:,:) ! array of infilled values. First column is epoch, rest are zero or infilled data.
  logical     , allocatable   :: missing(:)

! least squares variables
  real(kind=8),allocatable    :: A(:,:),At(:,:),B(:,:),ATA(:,:),VCV(:,:),AtB(:,:),soln(:,:)

! local variables
  integer*4    :: i,irow,icol,iter
  real*8       :: dt,comp,apr(3)
  character*256:: message

! allocate the temp values array
  allocate(values_tmp(nepochs,ncols+1))
  values_tmp = values_in
  allocate(missing(nepochs))

! allocate the LS arrays
  allocate(A(nepochs,3))
  allocate(At(3,nepochs))
  allocate(AtA(3,3))
  allocate(VCV(3,3))
  allocate(B(nepochs,1))
  allocate(AtB(3,1))
  allocate(soln(3,1))

  if(model(1:10) == "quadratic ")then
    write(message,'(a,a,a,a)')"    Infilling data type ",data_type," using model ",model
    if(debug)call status_update('STATUS',calling_prog,'infill_TS',' ',message,0)
  else
    write(message,'(a,a,a)')"Model ",model," not coded. Please code it yourself!"
    if(debug)call status_update('FATAL',calling_prog,'infill_TS',' ',message,0)
  endif

! PT180621: find the missing epochs
  do irow = 1,nepochs
    if(values_in(irow,1) == 0.d0)then
      missing(irow) = .true.
!print*,'missing epoch irow after ',values_in(irow-1,1)
    else
      missing(irow) = .false.
    endif
  enddo
      
! form up the least squares matrices
  do icol = 1,ncols
    A = 0.d0
    B = 0.d0
    At = 0.d0
    VCV =0.d0
    apr = 0.d0
  do iter = 1,3
    do irow=1,nepochs
      if(.not.missing(irow))then     ! it is an epoch that contains data. Use it in the matrices
        A(irow,1) = dble(irow)**2
        A(irow,2) = dble(irow)
        A(irow,3) = 1.d0
        comp = apr(1)*dble(irow)**2+apr(2)*dble(irow)+apr(3)
        B(irow,1) = values_in(irow,icol+1) - comp        ! if we assume our a priori parameter values are zero then the Obs-Comp is just the Obs !
      else
!        print*,'skip over record ',irow
      endif
    enddo

    ! now the LS solution
    call transp(A,At,nepochs,3)
    call matmult(At,A,AtA,3,nepochs,3)
    call invert(AtA,VCV,3)
    call matmult(At,B,AtB,3,nepochs,1)
    call matmult(VCV,AtB,soln,3,3,1)
    apr = apr + soln(:,1)
  enddo

    ! now, use the model to infill the missing values
    do irow = 1,nepochs
!print*,irow,values_in(irow,icol+1),apr(1)*dble(irow)**2+apr(2)*dble(irow)+apr(3)," all vals"
      if(missing(irow))then                                               ! it is a missing value
        if(icol == 1)then
          values_in(irow,1)  = values_in(irow-1,1) + epoch_step                                ! create the epoch (in graceseconds) of the missing epoch
          values_tmp(irow,1) = values_in(irow-1,1) + epoch_step                                ! create the epoch (in graceseconds) of the missing epoch
        endif
        values_tmp(irow,icol+1) = apr(1)*dble(irow)**2+apr(2)*dble(irow)+apr(3)  ! infill the value
      endif
!print*,irow,values_in(irow,icol+1),values_tmp(irow,icol+1)," blah"
    enddo
!stop 'stopped in infill_TS'
  enddo

! transfer values_tmp to values_in
  values_in = values_tmp

  deallocate(A)
  deallocate(At)
  deallocate(AtA)
  deallocate(VCV)
  deallocate(B)
  deallocate(AtB)
  deallocate(soln)
  deallocate(values_tmp)
  deallocate(missing)

  return
  end subroutine infill_TS
 
