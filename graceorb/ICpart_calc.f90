!CTITLE ICpart_calc

    subroutine ICpart_calc(maxnd,jd,t,x0,v0,efpos,a0) 

    use accel_mod
    use gm_mod         ! Needed for Gconst
    use inmod_mod    ! PT121221: added for inclusion/exclusion of mascon
    use coeff_mod    ! PT130802: added to get the value of C20 into this routine
    use bsscl_mod    ! PT140522: added to get bias and scale values into this routine
    use mascon_mod   ! PT160909: transfer to the new mascon arrays
    use rotation_mod ! PT170403: provides efixed-inertial matrices
    implicit none

!   Subroutine to the compure the IC partial derivative contribitions to the
!   acceleration terms.  This is d(central force acceleration)/d(xyz).
!   MOD TAH 110307: Add J2 contrribution to the partial.  Due to rotation terms
!   need to add jd and t to calling arguments to compute Inertial to Earth fixed rotation
! PT140430: added the partials of mascon accelerations (including ocean tide component) onto central body and J2 contributions
! PT150605: removed the local definition of GMe so that we use only the value defined in gmset.f90


    integer*4 , intent(in) :: maxnd           ! Dimensioning of arrays
    integer*4 , intent(in) :: jd              ! Julian day number
    real*8    , intent(in) :: t               ! Time of day (seconds)
    real*8    , intent(in) :: x0(3*maxnd)     ! Input pos'ns & IC partials (first 3 are XYZ pos'n)
    real*8    , intent(in) :: v0(3)           ! first 3 elements of v-array gives the actual velocity in inertial space
    real*8    , intent(in) :: efpos(3)        ! earth-fixed XYZ coords of satellite

    real*8    , intent(out):: a0(3*maxnd)     ! Computed partial acc'n terms (first are actual acc'ns)

 
!   Some constants should be assigned consistently later.
! PT150605: we should use the value assigned in gmset.f90 and passed in through gm_mod.f90
!    real(kind=8), parameter :: GMe = 3.986004415d14 ! m**3/s
! PT130214: make the J2 value compatible with the static gravity field valu
!    real(kind=8), parameter :: J2  =  1082.628d-6    ! J2 of Earth
! PT130802: why not compute this based on the C20 value from the input static gravity field
!    real(kind=8) , parameter :: J2  =  1082.62356d-6    ! J2 of Earth
    double precision J2
!   real(kind=8), parameter :: ae  =  6378155.0d0   ! Radius of Earth (m) Value in gm_mod

    real(kind=8), dimension(3,3) :: dadx   ! Partials of acceleration with respect to XYZ coordinates (sum of central+J2).
    real(kind=8), dimension(3,3) :: dadxc  ! Partials of central acceleration with respect to XYZ coordinates.
    real(kind=8), dimension(3,3) :: dadxm  ! Partials of J2 acceleration with respect to XYZ coordinates.
    real(kind=8), dimension(3,3) :: dadxmcon ! Partials of mascon acceleration with respect to inertial XYZ coordinates
    real*8 PF  ! Common factor in computing partials (depends on coordinates and row 3 of rotation matrix (z-component)
    real*8 GMF   ! GM factor for J2 calc.  (GM*J2*se**2)
    real*8 xr(3) ! Factors which are coordinate dependent and appear in partials.
    real*8 iz(6), ez(6)  ! Inertial z axis vector, earth fixed version.  

    real*8 rad,vmag      ! Radial distance to satellite.

    real*8 aip(3)      ! Accelerometer partial input (offsets and scales).
    real*8 aop(3)      ! Rotated accelerations elements
    real*8 posu(3)      ! unit position vector
    real*8 velu(3)      ! unit position vector
    real*8 component(3)      ! unit position vector
    real*8 poslat, poslong, vellat, vellong, Rpos, Rvel, qdec, qRA
    real*8 w2, x2, y2, z2

!  Mass con variables
    real*8 mcon_dist     ! Distance from satellite to Mass con
    real*8 mcon_efaccp(3) ! mcon earth fixed acceleration
    real*8 mcon_ifacc(3) ! mcon earth inertial acceleration
    real(kind=8) :: mcon_fact, def_fact
    integer :: i,j,k,l,m,n  ! Loop variables.
    logical bitmap

! PT130903: define some variables so that there are fewer floating point computations
    double precision rad3, rad5, rad7, rad9, GMF_5, GMF_35, PF_on_2
    integer*4 k_orig
    integer omp_get_thread_num

! PT170821: variable required in the call of efixed_inert (but not used in computations from this subroutine)
    real(kind=8) :: ifacc(3)

!  More mascon variables
    vmag = dsqrt(v0(1)**2+v0(2)**2+v0(3)**2)
    rad = dsqrt(x0(1)**2+x0(2)**2+x0(3)**2)
    rad3 = rad*rad*rad
    rad5 = rad3*rad*rad
    rad7 = rad5*rad*rad
    rad9 = rad7*rad*rad

! PT130802: derive here the value of J2 from the value of C20 that was in the input static gravity field file
    J2 = -dsqrt(5.d0)*meancoefC(2,0)

    posu(1:3) = x0(1:3)/rad
    velu(1:3) = v0(1:3)/vmag

     call cart_sph_convert(poslat, poslong, Rpos, posu)
     call cart_sph_convert(vellat, vellong, Rvel, velu)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Now add the central force partials

    do j = 1,3 
       do i = 1,3
          dadxc(i,j) = 3.d0*gm(1)*x0(i)*x0(j)/rad5
          if( i == j ) then
              dadxc(i,j) = dadxc(i,j) - gm(1)/rad3
          end if
       enddo
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Now start J2 contribution (dadxm).  We need to the rotation matrix for this
    iz = (/ 0.0d0, 0.d0, 1.d0 , 0.d0, 0.d0, 0.d0/)
!   Get Z-axis in earth fixed
    call efix_inert( rot_i2e,rotdot_i2e, rotacc_i2e, iz, ez, ifacc, .false.)

! DEBUG:
!    print*,'t,ez(1:3)',t,ez(1:3)

!   Compute common factors
    GMF = -gm(1)*J2*ae*ae 
    PF = ((3.d0*ez(1)**2-1)*x0(1)**2+(3.d0*ez(2)**2-1.d0)*x0(2)**2+(3.d0*ez(3)**2-1)*x0(3)**2  &
         + 6.d0*ez(1)*ez(2)*x0(1)*x0(2) + 6.d0*ez(2)*ez(3)*x0(2)*x0(3) + 6.d0*ez(1)*ez(3)*x0(1)*x0(3)) 

    xr(1) = (6.d0*ez(1)**2-2.d0)*x0(1)  +  6.d0*ez(1)*ez(2)*x0(2)  +  6.d0*ez(1)*ez(3)*x0(3)
    xr(2) =  6.d0*ez(2)*ez(1)*x0(1)  + (6.d0*ez(2)**2-2.d0)*x0(2)  +  6.d0*ez(2)*ez(3)*x0(3)
    xr(3) =  6.d0*ez(3)*ez(1)*x0(1)  +  6.d0*ez(3)*ez(2)*x0(2)  + (6.d0*ez(3)**2-2.d0)*x0(3)

! PT130903: define some things so that we do fewer floating point operations
    GMF_5 = 5.d0*GMF
    GMF_35 = 35.d0*GMF
    PF_on_2 = PF/2.d0

!   Now compute 3x3 matrix of contribtions 
    dadxm(1,1) = -GMF_5*x0(1)*xr(1)/rad7      + GMF_35*x0(1)**2*PF_on_2/rad9     &
                 -GMF_5*PF_on_2/rad7          + (3.d0*ez(1)**2-1.d0)*GMF/rad5
    dadxm(2,1) = -GMF_5*x0(2)*xr(1)/2.d0/rad7 + GMF_35*x0(1)*x0(2)*PF_on_2/rad9  &
                 -GMF_5*x0(1)*xr(2)/2.d0/rad7 + 3.d0*ez(1)*ez(2)*GMF/rad5
    dadxm(3,1) = -GMF_5*x0(3)*xr(1)/2.d0/rad7 + GMF_35*x0(1)*x0(3)*PF_on_2/rad9  &
                 -GMF_5*x0(1)*xr(3)/2.d0/rad7 + 3.d0*ez(1)*ez(3)*GMF/rad5

    dadxm(1,2) = -GMF_5*x0(2)*xr(1)/2.d0/rad7 + GMF_35*x0(2)*x0(1)*PF_on_2/rad9  &
                 -GMF_5*x0(1)*xr(2)/2.d0/rad7 + 3.d0*ez(1)*ez(2)*GMF/rad5
    dadxm(2,2) = -GMF_5*x0(2)*xr(2)/rad7      + GMF_35*x0(2)**2*PF_on_2/rad9     &
                 -GMF_5*PF_on_2/rad7          + (3.d0*ez(2)**2-1)*GMF/rad5
    dadxm(3,2) = -GMF_5*x0(3)*xr(2)/2.d0/rad7 + GMF_35*x0(2)*x0(3)*PF_on_2/rad9  &
                 -GMF_5*x0(2)*xr(3)/2.d0/rad7 + 3.d0*ez(2)*ez(3)*GMF/rad5

    dadxm(1,3) = -GMF_5*x0(3)*xr(1)/2.d0/rad7 + GMF_35*x0(3)*x0(1)*PF_on_2/rad9  &
                 -GMF_5*x0(1)*xr(3)/2.d0/rad7 + 3.d0*ez(1)*ez(3)*GMF/rad5
    dadxm(2,3) = -GMF_5*x0(3)*xr(2)/2.d0/rad7 + GMF_35*x0(3)*x0(2)*PF_on_2/rad9  &
                 -GMF_5*x0(2)*xr(3)/2.d0/rad7 + 3.d0*ez(2)*ez(3)*GMF/rad5
    dadxm(3,3) = -GMF_5*x0(3)*xr(3)/rad7      + GMF_35*x0(3)**2*PF_on_2/rad9     &
                 -GMF_5*PF_on_2/rad7          + (3.d0*ez(3)**2-1.d0)*GMF/rad5

! DEBUG
!    call printmat(dadxm,3,3,'dadxm ')

!   Add the J2 contribution to the Central force
! PT140514: turned off for a test
    dadx = dadxc + dadxm


! PT140430: convert the efixed mascon/ocean tide contribution to inertial
!    do i=1,3
!      do j=1,3
!        dadxmcon(i,j) = rot_e2i(i,1)*mcon_efpos_part(1,j) + rot_e2i(i,2)*mcon_efpos_part(2,j) + rot_e2i(i,3)*mcon_efpos_part(3,j)
!      enddo
!    enddo

! PT140430: add the mascon and ocean tide contribution.
!    dadx = dadx + dadxmcon

! DEBUG
!    call printmat (dadxc,3,3,'central_force_partials ')
!    call printmat (dadxm,3,3,'J2_partials ')
!    call printmat (mcon_efpos_part,3,3,'mcon_efpos_partials ')
!    call printmat (dadxmcon,3,3,'mcon_inert_partials ')
!    stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Compute contributions for each partial
! MOD TAH 110720: Increased from 3*nd to 3*(nd+1) (See loop below for bound)
    a0(4:3*(nd+1)) = 0   ! Clear vector before summing
! PT140815: changed this loop from starting at 1 to starting at 40 (to skip over pos/vel/scale/bias)
    do i = 1,nd      ! Loop over 6 partials (starts 1 up in index because first three are accelerations
                     ! (Use all partials because of change in Central and J2
                     ! forces from the accumulated motion due to the perturbations).
       do j = 1,3    ! Loop over XYZ
          do k = 1,3 ! Loop over XYZ
             a0(i*3+j) = a0(i*3+j) + dadx(j,k)*x0(i*3+k) ! Assign "acceleration" into velocity deriv slots
          end do
       enddo
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

! PT140815: KLUGE jump over the accelerometer and bias part. I should NEVER use a goto statement but I'm going to
!    goto 1234

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Now add in the contributions for accelerometers.  These are paired as scale
!   and offset for the accelerator axes directions)
! APP130328: Order of bias and scale partials swapped to agree with what
!            gracefit expects
    k = 21
    do i = 1,3   ! Loop over along, cross and radial directions
       aip = 0
!       aip(i) = accrom(i,1)  ! Scale (unit factor; maybe add multiplier later)
! APP130328: Correct error in calculation of scale partials - use raw observations
! PT140522: if we use the uncalibrated observations as the partial for scale, we end up
!           with effectively a constant number of ~1-30 um/s^2, with a temporal perturbation that 
!           is around 0.05 um/s^. The consequence is that the scale and bias partials are very
!           highly correlated. To mitigate this, I want to try first calibrating the acc obs,
!           so that the scale partial will range between 0-100 nm/s^2, while the bias stays
!           at a (calibrated) static value. So, we have the equation:
!          
!  ACC_calibrated = bias_adj + scale_apr*(ACC_obs + B_o)
! 
!  where bias_adj  = the new "bias" parameter which will be the adjustment to the a priori bias value (now zero, since already calibrated),
!        scale_adj = the scale factor for which we want to solve
!
!  therefore,  d/dScale = (ACC_obs + B_o)

       aip(i) =  accobs(i)+ bs(i)
       call matmult(Qrow,aip,aop,3,3,1)
       do j = 1,3  ! Copy values over
          k = k + 1
!         Add direct contribtion
          a0(k) = a0(k) + aop(j)
       end do
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Now do the accelerometer offsets
    do i = 1,3   ! Loop over along, cross and radial directions
       aip = 0
       aip(i) = 1  ! Offset (PT140514: in spacecraft reference frame)
       call matmult(Qrow,aip,aop,3,3,1)
!       component(i)= aop(1)*posu(1)+aop(2)*posu(2)+aop(3)*posu(3)
       do j = 1,3  ! Copy values over
          k = k + 1
!         Add the direct contribution. Make offset units be micro-m/s^2.
          a0(k) = a0(k) + aop(j)   *1.d-6 
       end do
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ADD MASS CON partials
!
!1234   k = 40   ! KULGE jump over the scale/bias partials above lands us here and we need to increment k to what it would have been

! PT140904: we need to increment k here to account for the 12 elements that relate to the 1/rev and 2/rev that are now computed in
!           the more rigorous way (not in ICpartcalc)
   k = k + 12
! PT170405: we need to increment k again to account for the 9 elements  related to the roll/pitch/yaw perturbations that are computed (using calc_perturb_acc)
   k = k + 9

if (gt_mcon(1).eq."Y") then

    k_orig = k
    do i = 1, total_prim  ! PT161005: changed counter from "num_mcon" to the new name

! Rotate efixed partials to inertial space
      do j=1,3                            
        k = k_orig + (i-1)*3 + j
        a0(k) = a0(k) + rot_e2i(j,1) * mcon_efacc_part(i,1) + &
                        rot_e2i(j,2) * mcon_efacc_part(i,2) + &
                        rot_e2i(j,3) * mcon_efacc_part(i,3)
      enddo
    enddo

endif

! END MASS CON Code          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT140613:   ADD MASS CON tidal amplitude partials
!
! PT170303: writing out the partials for the tidal amplitudes should happen irrespective of whether we use
!           non-zero apriori values or not !!!
! PT170303: but it doesn't work properly if we don't .... so I've left gt_use_apriori_msc_tide(1) set to "Y" in graceorb.f90
!           even when we haven't included apriori values !!
if (gt_use_apriori_msc_tide(1).eq."Y") then

  do i = 1, total_ocean_prim   ! PT161005: changed counter from "num_mcon" to the new name

! Rotate efixed partials to inertial space for the tidal consituents that we wanted to estimate
    if(mcon_ocean_prim(i,2) > 0)then              ! there are tidal amplitudes to estimate for this mascon
      do l = 1,max_msc_tides
! DEBUG
!   print*,'i,l,bitmap(mcon_tide1(i),j)',i,l,kbit(mcon_tide1(i),j),mcon_eftid_part(:,l,1,i)
        if(bitmap(mcon_ocean_prim(i,2),l)) then     ! we want to estimate the amplitudes for this tidal constituent PT150824: bug here on index. Was "j" should have been "l"
!print*,'tidal partial for mascon',i,' tide ',l,mcon_eftid_part(1,l,1,i)
          k_orig = k
          do j=1,3                                     
            k = k_orig + j
! sine amplitude
!print*,'tidal partial imsc,tide,component,column',i,l,j,k
            a0(k) = a0(k) + rot_e2i(j,1) * ocean_eftid_part(1,l,1,i) + &
                            rot_e2i(j,2) * ocean_eftid_part(2,l,1,i) + &
                            rot_e2i(j,3) * ocean_eftid_part(3,l,1,i)
! cosine amplitude
            a0(k+3) = a0(k+3) + rot_e2i(j,1) * ocean_eftid_part(1,l,2,i) + &
                                rot_e2i(j,2) * ocean_eftid_part(2,l,2,i) + &
                                rot_e2i(j,3) * ocean_eftid_part(3,l,2,i)

! DEBUG
!  print*,'ICpart_calc: i,l,j,k,a0(k),a0(k+3)',i,l,j,k,a0(k),a0(k+3),3*maxnd

          enddo         ! end of coordinate component loop
! PT140708: add three here to the counter to account for the cosine amplitudes that have been added into the a0 array
          k = k + 3
        endif         ! end of bitmap condition test
      enddo         ! end of tidal constituent loop
    endif         ! end of bit-map condition test
  enddo         ! end of mascon loop
endif

! END MASS CON tide amplitude Code          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!   Thats all
    return
    end
