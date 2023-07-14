   subroutine bias_from_ACC

! subroutine to derive a priori values for the XYZ accelerometer biases directly from the accelerometer obs themselves.
! The assumptions are as follows:
!    1. the observations contain a constant part (the bias) and fluctuations with a 1/rev period. 
!    2. The fluctuations are of the order of 10-200 nm/s^2
!    3. the constant part for each bias lies somewhere between 0.5 and 30 um/s^2 (*** 1 to 15 orders of magnitude greater!)
!
!    4. for Bsx (along-track) the max calibrated ACC value would be just under zero. Range -100 to 0.
!    5. for Bsy (cross-track) the value is zero in eclipse, positive or negative out of eclipse. Range 0 to +/- 70
!    6. for Bsz (radial)      the value is just negative in eclipse, strongly positive at all other times. Range -20 to 80.
!
!  Therefore, for Bsx and Bsz, we can compute the mean value and then adjust it by 50 nm/s^2 to sit the values in the right
!  place. For Bsy we can perhaps use the shadow information to adjust the values to zero during eclipse.
!
! P. Tregoning
! 28 May 2014

    use bsscl_mod
    use accred_mod

    implicit none
    integer       :: i, i2
    real(kind=8)  :: meanbias(3)
    character*200 :: message
    integer       :: nvals(3)

! calculate the mean value of each component
    meanbias = 0.d0
    nvals = 0
    do i=1,3
! PT/RMcG190820: cut off 1000 epochs from beginning/end to eliminate fft edge effects. This will throw away 2000 epochs from data if
!                non-extended accelerometer data are used, but will leave 84400 epochs anyway, so should be ok.
      do i2 = 1000,ACCn-1000 
         if (accr(i,i2) < 1e-3) then  
             nvals(i) = nvals(i) + 1
             meanbias(i) = meanbias(i) + accr(i,i2)
         endif
      enddo 
    enddo
    
    meanbias = -1.d6 * meanbias / nvals    ! sign change so that the bias removes the mean. Biases are in um/s^2 at this stage

! transfer to the variables used for biases. Adjust Bsx and Bsz down by -50 and -30 nm/s^2. Bsy could adjust with either sign, so leave alone
    c0x = meanbias(1) - 50.d-3
    c0y = meanbias(2)
! PT140822: based on 2010-09-07 and 08, the value for bsz is too small by ~35 nm/s^2. So increase it by that much (ie add 5 rather than subtract 30
    c0z = meanbias(3) + 5.d-3

    write(message,'(a,3f15.5,a)')"Computed mean accelerometer values (Bsx Bsy Bsz):",meanbias(1)-0.05d0,meanbias(2) &
                                   ,meanbias(3)-0.03d0," um/s"
    call status_update('STATUS','GRACEORB','bias_from_ACC',' ',message,0)


    return 
    end
