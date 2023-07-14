  program mascon_threshold
  
! program to read in a vcv file and output non-zero values for only mascons that received a sufficiently large adjustment
!
! P. Tregoning
! 14 November 2022

  use mascon_mod          ! defines all the mascon arrays

  implicit none

  character :: vcv_file*150, vcv_file_out*150, arg*100
  character*150 :: mascon_file
  character*1   :: zero_ocean
  real(kind=8) :: adjust_size                 ! |adjust_size| is the threshold for a mascon to have a non-zero output value
  real(kind=8) :: incl_adjust_percent         ! percentage of adjustment to add to apriori value
  real(kind=8),allocatable :: msc(:,:)        ! mascon apriori, adjustment and estimate values
  character*81,allocatable :: mascon_lines(:)

  integer*4 :: imsc,nmascon_t
  character*130 :: line
  real(kind=8),allocatable :: msc_densities(:)
  
! unit numbers
  integer*4,parameter :: lumsc = 12
  
  ! PT221117: use the lib routines to read in the mascon file so that the densities of the mascons can be known

! decode command line
  call getarg(1,vcv_file)
  if (vcv_file(1:1) == "")then
    print*,"max_adjust vcv_file_in vcv_file_out adjustment_size %adjust_to_include mascon_file zero_ocean [y/n]"
    stop
  endif

  call getarg(2,vcv_file_out)
  
  call getarg(3,arg)
  read(arg,*)adjust_size
  call getarg(4,arg)
  read(arg,*)incl_adjust_percent
  call getarg(5,mascon_file)
  call getarg(6,zero_ocean)
  
! open the input and output files
  open(10,file=trim(vcv_file),status='old')
  open(11,file=trim(vcv_file_out),status='unknown')

! read and transfer the header part until we get to the first mascon entry
  line = " "
  do while (line(7:9) /= " MC")
    read(10,'(a)')line
    ! get the number of mascons
    if(line(1:20) == "Solution with mascon")read(line(103:108),*)nmascon_t
    if(line(7:9) /= " MC")write(11,'(a)')line
  enddo

! allocate the mascon arrays
  allocate(mascon_lines(nmascon_t))
  allocate(msc(nmascon_t,3))
  allocate(msc_densities(nmascon_t))
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! we are now at the first line of the mascon entries
  backspace(10) 
  do imsc = 1,nmascon_t
    read(10,'(a81)')mascon_lines(imsc)
    read(mascon_lines(imsc)(32:47),*)msc(imsc,1)   !   apriori mascon value
    read(mascon_lines(imsc)(48:64),*)msc(imsc,3)   ! estimated mascon value
    msc(imsc,2)  = msc(imsc,3) - msc(imsc,1)    ! adjustment to apriori mascon value
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the mascon file 
  call read_msc_hdr(lumsc,mascon_file,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                            ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
  call allocate_mascon_arrays
  ! PT180220: using trimlen here doesn't work - the subroutine is expecting a C*150, so pass the whole thing through
  !  call read_mascon_file(lumsc_in,mascon_file(1:trimlen(mascon_file)))
  call read_mascon_file(lumsc,mascon_file)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! PT221116: for now, just read a mascon density file. Ultimately, it should read the proper mascon file
!  open(12,file="msc.densities",status='old')
!  do imsc=1,nmascon_t
!    read(12,*)msc_densities(imsc)
!  enddo
  
! loop through the mascons and:
! -  set the estimated value to the apriori value if the adjustment was below the threshold input by the user
! -  change the estimated to be the apriori + adjustment*percentage input by the user

  do imsc = 1,nmascon_t
    ! small adjustments get set to zero
    if (dabs(msc(imsc,2)) < adjust_size)then
      msc(imsc,2) = 0.d0
    else if (zero_ocean == "Y" .and. mcon_prim(imsc,6) > 1010.d0)then
      msc(imsc,2) = 0.d0    
    else
      print*,"mascon",imsc,msc(imsc,:)
    endif

    ! estimate = apriori + adjustment*percentage
    msc(imsc,3) = msc(imsc,1) + msc(imsc,2) * incl_adjust_percent / 100.d0
    write(mascon_lines(imsc)(31:65),'(2f17.7)')msc(imsc,1),msc(imsc,3)
    write(11,'(a)')mascon_lines(imsc)

  enddo

! finish off the vcv file
!  write(11,'(a)')"  VCV SOLUTION"

  close(11)

  end




