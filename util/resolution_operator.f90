  program resolution_operator

! program to calculate the resolution operator and, hence, the error due to regularisation for a given set of AtWA and a regularisation file.
! We should just build this into addnorm eventually .....
!
! P. Tregoning
! 22 September 2020

  use mascon_mod    ! provides all the arrays for reading the mascon file

  implicit none

! command line variables
  character*150 :: AtWA_file                ! file of AtWA information
  character*150 :: reg_file                 ! regularisation matrix
  character*150 :: EWH_file                 ! file of EWH anomalies for each mascon
  character*150 :: output_file              ! name of output file
  real(kind=8)  :: lambda                   ! scale factor for the regularisation matrix

! matrices required here
  real(kind=8), allocatable :: AtWA(:,:)    ! normal equations, read from a file
  real(kind=8), allocatable :: P(:,:)       ! regularisation matrix
  real(kind=8), allocatable :: AtWA_P(:,:)  ! AtWA + P
  real(kind=8), allocatable :: R(:,:)       ! the resolution operator
  real(kind=8), allocatable :: VCV(:,:)     ! (AtWA + P)^-1
  real(kind=8), allocatable :: EWH(:)       ! EWH anomalies for each mascon
  real(kind=8), allocatable :: reg_error(:) ! regularisation-induced error on each mascon

! number of parameters in the AtWA file information
  integer*4    :: nparam_t,nmascon_t,nICS
  real(kind=8) :: version

! variables for lapack routines
  integer*4 :: info,lwork,lda
  real(kind=8), allocatable :: work(:)
  integer*4,    allocatable :: ipiv(:)

! local variables
  integer*4     :: imsc
  character*250 :: line,arg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! decode command line
  call getarg(1,AtWA_file)
  if(AtWA_file(1:1) == "")then
    print*,"Runstring: resolution_operator   AtWA_file   regularisation_file   EWH_file  output_file lambda"
    stop
  else
    open(10,file=AtWA_file,status='old')
  endif

  call getarg(2,reg_file)
  open(11,file=reg_file,status='old')

  call getarg(3,EWH_file)
  open(13,file=EWH_file,status='old')

  call getarg(4,output_file)
  open(12,file=output_file,status='unknown')

  call getarg(5,arg)
  read(arg,*)lambda
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! open and read the input files
  call status_update('STATUS','UTIL','resolution_operator',AtWA_file,"Reading normal equations file",0)


  ! read the number of parameters in file
  read(10,*)version,nparam_t,nmascon_t,nICs

  ! allocate the arrays
  allocate(AtWA(1:nmascon_t,1:nmascon_t))
  allocate(P(1:nmascon_t,1:nmascon_t))
  allocate(AtWA_P(1:nmascon_t,1:nmascon_t))
  allocate(VCV(1:nmascon_t,1:nmascon_t))
  allocate(R(1:nmascon_t,1:nmascon_t))
  allocate(EWH(nmascon_t))
  allocate(reg_error(nmascon_t))

  ! allocate the arrays for lapack routines
  LDA = nmascon_t
  LWORK = nmascon_t*nmascon_t
  ALLOCATE (WORK(LWORK))
  ALLOCATE (IPIV(nmascon_t))

  ! read in the AtWA matrix
  do imsc=1,nmascon_t
    read(10,*)AtWA(imsc,:)
  enddo
  close(10)

! open and read the regularisation file
  call status_update('STATUS','UTIL','resolution_operator',reg_file,"Reading regularisation file",0)

  ! skip two lines
  read(11,'(a100)')line
  print*,line(1:100)
  read(11,'(a250)')line
  print*,line

  ! now read in the regularisation matrix
  do imsc=1,nmascon_t
    read(11,*)P(imsc,:)
  enddo
  close(11)  
  ! multiply the regularisation matrix by the lambda scale factor
  P = lambda * P

  ! now read in the mascon EWH values to apply
  call status_update('STATUS','UTIL','resolution_operator',EWH_file,"Reading EWH solution file",0)

  ! now read in the regularisation matrix
  do imsc=1,nmascon_t
    read(13,*)EWH(imsc)
  enddo
  close(13)  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now, solve eq 21 of Loomis et al (2019) to derive the corresponding resolution operator, R

!  ! add the regularisation matrix
!  AtWA_P = AtWA + P
!
!  ! invert it
!  call status_update('STATUS','UTIL','resolution_operator',AtWA_file,"Inverting (AtWA + P)",0)
!  call inverse(AtWA_P,VCV,nmascon_t)
!
!  ! multiply this by AtWA
!  call status_update('STATUS','UTIL','resolution_operator',AtWA_file,"Compute VCV * AtWA",0)
!  call matmult(VCV,AtWA,R,nmascon_t,nmascon_t,nmascon_t)
!
!  ! multiply R x EWH
!  call matmult(R,EWH,reg_error,nmascon_t,nmascon_t,1)


  !***** using lapack routines ****
  VCV = AtWA + P   ! temp storage. The VCV matrix is eventually created by DGETRI below ...

  call status_update('STATUS','UTIL','resolution_operator',' ',"Calling DGETRF",0)
  CALL DGETRF( nmascon_t,nmascon_t, VCV, LDA, IPIV, INFO )
  call status_update('STATUS','UTIL','resolution_operator',' ',"Calling DGETRI",0)
  CALL DGETRI(nmascon_t, VCV, nmascon_t, IPIV, WORK, LWORK, INFO)

  ! multiply this by AtWA
  call status_update('STATUS','UTIL','resolution_operator',AtWA_file,"Compute VCV * AtWA",0)
  R = matmul(VCV,AtWA)

  ! subtract 1 from the diagonal of R
  do imsc=1,nmascon_t
    R(imsc,imsc) = R(imsc,imsc) - 1.d0
  enddo

  ! multiply modified R x EWH. That should be Eq 17 of Loomis et al (2019)
  reg_error = matmul(R,EWH)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write out computed matrix
  call status_update('STATUS','UTIL','resolution_operator',output_file,"Writing out resolution operator matrix",0)

  ! write the number of mascons
  write(12,*)nmascon_t,"   : number of mascons. Resolution operator of Loomis et al (2019)."
  
  ! write the names of the two files used to create this
  write(12,'(a,2x,a)')AtWA_file,reg_file

  ! now write out the regularisation error
  write(12,'(a)')"REGULARISATION INDUCED ERROR"
  do imsc=1,nmascon_t
    write(12,*)sum(R(imsc,:)),reg_error(imsc),EWH(imsc)
  enddo

  ! write out the R matrix
  write(12,'(a)')"RESOLUTION OPERATOR"
  do imsc=1,nmascon_t
    write(12,*)R(imsc,:)
  enddo
  close(12)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end

