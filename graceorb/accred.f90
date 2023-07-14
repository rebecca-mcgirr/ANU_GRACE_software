    subroutine accred(tin)

! Emma-Kate Potter, 18 April, 2011 
! This subroutine reads in the accelerometer data spanning the 
! time range of the integration

! MODIFIED: APP 121004
! To allow user-defined names for the star camera and accelerometer files

    use accred_mod
    use inmod_mod

    implicit none
    integer :: ACCnum, SCAnum 
    integer*4 :: i, ioerr  
    real(kind=8) :: temp1, temp2, temp3, temp4, temp5 
    real(kind=8) :: tin
    real(kind=8) :: timei, timef

    ACCnum=61
    SCAnum=60

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! open input files for accelerometer (ACC) and star camera (SCA)
    open (unit=SCAnum, file=SCA_file, status='unknown')
    open (unit=ACCnum, file=ACC_file, status='unknown')
!    open (unit=SCAnum, file='SCA.dat', status='unknown')
!    open (unit=ACCnum, file='ACC.dat', status='unknown')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Read accelerometer data in, allocate array size and save in array 
    ioerr = 0
    ACCn=0
    do while (ioerr.eq.0)  ! count how many lines in the ACC file 
      read (ACCnum, *,iostat=ioerr) temp1, temp2, temp3, temp4
      ACCn=ACCn+1
    enddo

    rewind(ACCnum)         ! rewind ACCfile

    allocate (ACCtime(1:ACCn-1))   ! allocate array to ACC file 
    allocate (accr(1:3,1:ACCn-1)) !
    
    do i=1,ACCn-1
      read (ACCnum, *) ACCtime(i), accr(1,i), accr(2,i), accr(3,i)
    enddo 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Read star camera data in, allocate array size and save in array 
    ioerr = 0
    STARn=0
    do while (ioerr.eq.0)  ! count how many lines in the ACC file
      read (SCAnum, *,iostat=ioerr) temp1, temp2, temp3, temp4, temp5
      STARn=STARn+1
    enddo

    rewind(SCAnum)         ! rewind ACCfile

    allocate (SCAtime(1:STARn-1))   ! allocate array to ACC file
    allocate (varr(1:4,1:STARn-1)) !

    do i=1,STARn-1
      read (SCAnum, *) SCAtime(i), varr(1,i), varr(2,i), varr(3,i), varr(4,i)
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    return
    end
