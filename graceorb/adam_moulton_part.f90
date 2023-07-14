    subroutine adam_moulton_part(coor, onepr,twopr, jd, t, nsteps, h, quatSRF2e_init,nepochs &
                          ,rot_dates,rot_mats,rotdots,rotaccs &
! PT140319: pass the Xp, Yp values into adam_moulton_part-> gravcalc -> sphharmfield -> poltid
! PT161005: now done with usno_mod.f90
                          )

!  E-K Potter May 2010
!  This subroutine takes initial conditions for satellite position and velocity
!  and uses an adams-moulton predictor-corrector method to integrate the 
!  satellite orbit
!
!  This code adapted from adams-moulton routine found on internet (??)
!
!  The arrays below allow storage of the results for the four most recent
!  time steps.  For example y(0) contains the value of y at the beginning
!  of the current time step, y(1) contains the value of y at the beginning
!  or the previous time step, y(2) contains the value of y 2 steps back, etc.

! MODIFIED: April 3 2012 - T Purcell
! The code now uses the satellite's position to determine which mascon it
! overlies. This ternary mascon is related to a primary mascon through the
! mcon_index array. Having located the primary mascon the program then uses
! the primary (coarse) resolution for all mascons more than twenty degrees away
! the secondary resolution for intermediate sites, and the ternary (high)
! resolution for sites less than five degrees away. Rather than using the
! same list of mascons for all locations and integration steps (as assumed in
! the original version of the code). The new code constructs an updated mascon
! list each time this routine is called.

! MOD TAH 110301: Added partials for EF coordinates and velocities with
!  respect to the interial IC position and velocity.  With partials we
!  integrate 3 pos/vel and 6 sets of 3 position and velocity (PV) partials.  
!  This code will need accelerometer biases to be added to some point 
!  as well.  This is another 3 sets of 3 PV partials.
!  We extend coor to also contain partials.
!
! MOD PT130903: add some openMP lines to use parallel processing on some of the 3*nd loops
!
! MOD PT131104: write output to a binary GTORB file now (rather than the original ascii file)
!
! MOD PT140318: change the way that inertial<->efixed is calculated so that it uses IERS2010 standards
!
! MOD PT140812: change the way that the IC partials are integrated. We now integrate two orbits - the actual
!               satellite location and the position-perturbed satellite location. We then difference them to calculate
!               the effect of the perturbation at each time step. Dividing by the perturbation itself then
!               gives us the partial derivative for position.
!
! MOD PT140820: add a twice-per-rev parameterisation to the along-track acceleration (2pr cos and sine amplitudes are the new parameters)
!
! MOD PT170403: moved rot_e2i and rot_i2e into rotation_mod
!
! MOD PT170404: add perturbation of a roll offset error to compute partials for a new parameter ...
!
! MOD PT170609: define record length through mod_subs/gtorb_mod
!
! MOD SA171006:  Add HDF5
!
! MOD PT181208: replace all large static allocations with dynamic memory allocation

    use coeff_mod
    use spherhar_mod
    use gm_mod
    use rotation_mod  ! PT170403: now contains declarations for rot_e2i and rot_i2e
    use accel_mod     ! Added for mass con definitio
    use inmod_mod     ! PT121221: added to switch on/off mascon inclusion
    use bsscl_mod     ! PT140819: added to have access to the accelerometer biases
    use mascon_mod    ! PT161005: the new mascon definition module (defines total_prim to replace num_mcon)
    use gtorb_mod     ! PT170609: this declares the record length information for the GTORB file

    implicit none

! PT170609: now done in mod_subs/gtorb_mod
!    integer*4                  :: nrec     ! current line in output binary file. Need to increment before using it.
    real(kind=8), dimension(6) :: coor     ! Passed initial pos'ns & velocities.
    real(kind=8), dimension(2) :: onepr    ! amplitudes of  once-per-rev along-track acceleration empirical function
    real(kind=8), dimension(2) :: twopr    ! amplitudes of twice-per-rev along-track acceleration empirical function


!    integer*4, parameter:: maxnd = 7030 ! No. of parameters and partials to be integrated
                      ! for orbit and IC partial, nd = 7.  
                      ! With acelerometers it is 13 (3 biases and 3 scales)
                      ! With mass cons it is +number of mass cons 
! PT140722: increased maxnd from 5000 to 6300 to allow 860 (x2 components) tidal amplitudes to be estimated
! PT150827: increased to allow 1015 (x2) tidal amplitudes, which equals 5 constituents on 203 mascons.
! PT170609: increased to 7550 to account for having 7526 mascons
! PT181204: increased to 19000 to account for having 18517 mascons
! RM190305: increased to 48000 to account for having 47414 mascons
    integer*4, parameter:: maxnd = 48000 ! No. of parameters and partials to be integrated

    integer*4   , intent(in) :: nepochs               ! dimensioning of rotation and rotation rate arrays
    real(kind=8), intent(in) :: rot_mats(3,3,nepochs) ! array of rotation matrices 
    real(kind=8), intent(in) :: rotdots(3,3,nepochs)  ! arrays of rotation matrices and their rates
    real(kind=8), intent(in) :: rotaccs(3,3,nepochs)  ! arrays of rotation matrices and their rates
    real(kind=8), intent(in) :: rot_dates(nepochs)    ! vector of mjd epochs of rotation matrices and their rates

! variables related to inertial<->efixed computations
    real(kind=8) :: mjd
! PT170403: moved these declarations to rotation_mod
!    real(kind=8),dimension(3,3) :: rot_i2e,rot_e2i,rotdot_i2e,rotdot_e2i,rotacc_i2e,rotacc_e2i

    integer*4 :: jd,jds,ioerr
    integer*4 :: timeout     
    integer*4 :: nsteps, i, j,  k, m, n
!    integer :: istep
    real(kind=8) ::  xpole, ypole
    real(kind=8):: t
    real(kind=8) :: h, time, tnew, tn
    real(kind=8) :: angle
    real(kind=8), dimension(3) :: unit_SRF, unit_iner, unit_efix
! PT181208: change these to allocatable to make better use of the memory of the computer
    real(kind=8), allocatable :: anew(:), xnew(:), vnew(:), acc(:)
    real(kind=8), allocatable :: efacc(:), efx(:)
    real(kind=8), allocatable :: x0(:), x1(:), x2(:), x3(:), x4(:) 
    real(kind=8), allocatable :: efx0(:), efx1(:), efx2(:), efx3(:), efx4(:) 
    real(kind=8), allocatable :: v0(:), v1(:), v2(:), v3(:), v4(:) 
    real(kind=8), allocatable :: a0(:), a1(:), a2(:), a3(:), a4(:) 
    real(kind=8), allocatable :: xam(:), vam(:), xabt(:), vabt(:) 
    real(kind=8), allocatable :: efxabt(:) 
    real(kind=8), allocatable :: efic(:), incoor(:)
    real(kind=8), allocatable :: efpos(:) ! Contains pos'ns & velocities. 
    real(kind=8), allocatable :: x(:,:), v(:,:), a(:,:) 

    real(kind=8), dimension(0:3) :: quati2e, quatSRF2i, quatSRF2e, quattest 
    real(kind=8), dimension(3,3) :: rottest, rotSRF2e , quat_rot, srf2e_rot
    real(kind=8), dimension(0:3) :: quatSRF2e_init
    real*8 :: rlat_temp, rlong_temp, rx_temp

    real*8 xvp_ef(6)   ! IC in Earth fixed with +0.5 m change 
    real*8 xvm_ef(6)   ! IC in Earth fixed with -0.5 m change
    real*8 xvp_if(6)   ! IC in inertial with +0.5 m change
    real*8 xvm_if(6)   ! IC in inertial with -0.5 m change
    real*8 dxv(6)      ! Change in IC from Efixed to Inertial (partial)

! PT140812: declare some variables to contain the accelerations of the satellite and perturbed orbits
    real(kind=8) :: acc_sat(3)                               ! total inertial acceleration computed at satellite location
    real(kind=8) :: statgrav_perturb(3)                      ! static gravity field acceleration at the perturbed orbit location
    real(kind=8) :: statgrav_sat(3)                          ! static gravity field acceleration at the satellite location
! PT220707: we don't use the 1pr and  rpy partials anymore, so why compute them? Change from 19 to 12
    integer*4    :: norbs_perturbed=12                       ! number of perturbed orbits that we need to keep track of (X,Y,Z, XYZdot, 3 x bias, 3 x scale, 1pr(S/C), 2pr(S/C), roll/pitch/yaw offset error)
    real(kind=8) :: pos_pert_if(6),pos_pert_ef(6)            ! perturbed positions  in ef and inertial
    real(kind=8) :: vel_pert_if(6),vel_pert_ef(6)            ! perturbed velocities in ef and inertial
    real(kind=8),allocatable :: satpos_perturbed_if(:)       ! actual and 3 sets of perturbed satellite position (+ive  X,Y,Z perturbation)
    real(kind=8),allocatable :: accel_pert_if(:)             ! total accelerations at the locations of the perturbed orbit

! PT140811: set the magnitude of the parameter perturbations
    real*8, allocatable :: orbit_perturb(:)

! PT191030: variable to save the SCA-to-inert quaternion, to then calculate the SCA-to-EF quaternion
    real(kind=8) :: quat_SCA(0:3)

    character*3 :: sflag
    integer :: idir
    integer :: tempint

    real*8 :: PEPt
    integer*4 :: PEPjd

    logical output_tide

! DEBUG variables
    real(kind=8)  :: roll,pitch,yaw,d_roll,d_pitch,d_yaw
    logical :: debug
! satellite acceleration (determined roughly, as (V_i+1 - V_i)/dt
    real*8 sat_acc(3),efvel_old(3)

! PT140430: partials of satellite coords wrt mascon accelerations (computed in gravcalc, passed to and used in ICpart_calc)
!    real(kind=8) :: mcon_efpos_part(3,3), mcon_efacc_part(num_mcon,3)
! DEBUG
    real*8 :: vel_plus_ef(6),vel_minus_ef(6),vel_plus_if(6),vel_minus_if(6),dvel(6)
    real*8 :: ifacc(3)    ! temporary attempt to work out what the inertial accelerations are at a particular epoch. 

    save efvel_old

! PT181208: allocate all the really large arrays using dynamic allocation
  allocate(x(1:3*maxnd, 0:3))
  allocate(v(1:3*maxnd, 0:3))
  allocate(a(1:3*maxnd, 0:3))
  allocate(acc(3*maxnd))
  allocate(efacc(3*maxnd))
  allocate(xam(3*maxnd))
  allocate(vam(3*maxnd))
  allocate(xabt(3*maxnd))
  allocate(vabt(3*maxnd))
  allocate(efxabt(3*maxnd))
  allocate( efic(6*maxnd))
  allocate(incoor(6*maxnd))
  allocate(efpos(6*maxnd))
  allocate(anew(3*maxnd))
  allocate(xnew(3*maxnd))
  allocate(vnew(3*maxnd))
  allocate(x0(3*maxnd))
  allocate(x1(3*maxnd))
  allocate(x2(3*maxnd))
  allocate(x3(3*maxnd))
  allocate(x4(3*maxnd))
  allocate(v0(3*maxnd))
  allocate(v1(3*maxnd))
  allocate(v2(3*maxnd))
  allocate(v3(3*maxnd))
  allocate(v4(3*maxnd))
  allocate(a0(3*maxnd))
  allocate(a1(3*maxnd))
  allocate(a2(3*maxnd))
  allocate(a3(3*maxnd))
  allocate(a4(3*maxnd))
  allocate(efx0(3*maxnd))
  allocate(efx1(3*maxnd))
  allocate(efx2(3*maxnd))
  allocate(efx3(3*maxnd))
  allocate(efx4(3*maxnd))


! allocate variables
    allocate(accel_pert_if(3+norbs_perturbed*3) )
    allocate(satpos_perturbed_if(3+norbs_perturbed*3))
    allocate(orbit_perturb(norbs_perturbed+7))   ! PT220707: for now just allocate it to 19 elements
! PT161128: we need to allocate mcon_efacc_part somewhere, so I'll do it here
! PT210806: only do this if we want to include mascon partials (ie gt_mcon(1) == "Y")
    if(gt_mcon(1) == "Y")allocate(mcon_efacc_part(total_prim,3))

    if(t == 43200.d0) efvel_old = 0.d0

    call status_update('STATUS','GRACEORB','adam_moulton_part',' ',' Begin integration',0)

!   open output file for ef position
    open (unit=41, file='efpos', status='unknown') 

    flagrun = 0

    time = t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Initialize all the xnew and vnew first
    xnew(:) = 0.d0
    vnew(:) = 0.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   set initial conditions for 6 components of velocity and position 
    do j = 1,3
      xnew(j) = coor(j)
      vnew(j) = coor(j+3)
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Mapping an EFixed initial condition to Intertial.  NOTE: for 
!   position we use the rotation matrix and the derivative of the 
!   rotation matrix WRT initial conditions.  (Calculation done numerically).
    mjd = jd + t/86400.d0 - 2400000.5d0
    call inert_interp(mjd,nepochs,rot_dates,rot_mats,rotdots,rotaccs,rot_i2e,rotdot_i2e &
                      ,rot_e2i,rotdot_e2i,rotacc_e2i,rotacc_i2e)
!call \mat(rot_i2e,3,3,'adam_0_roti2e ')
!print*,'mjd=',mjd,jd,t
!stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! PT140811:      O R B I T    P E R T U R B A T I O N S
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! position(0.5m)
    orbit_perturb(1:3)  = 0.5d0     ! perturb the positions  by 0.5 m
! velocity (0.5mm)
    orbit_perturb(4:6)  = 0.05d-3    ! perturb the velocities by 5 mm/s
! scale. I want a 25 nm/s^2 effect so we need to work out an appropriate perturbation based on the size of the bias (since the scale
!        is multiplied by the bias in accelerom.f90)
    do i=1,3
! PT180927: if we want 25 nm/s^2 (as per the comment) then why not just set the perturbation to be 25.e-9 ? The bias values for 
!           cross-track are 31 um/s^2 (GRACE A) and  12 um/s^2 (GRACE B), whereas they are ~1 um/s^2 for the along-track and radial. So
!           the perturbation of scale for cross track is 12-30 times smaller than for the other two. Why? Also, we never multiply the
!           scale by the bias in accelerom.f90 .... the comment above is wrong!
!
!           set the scale perturbation to 3% - arbitrarily chosen to be an average (roughly) of what was being used for along- and radial
      orbit_perturb(6+i)  = 0.03  ! 0.025d-6/bs(i)    ! perturb the accelerometer scale so that the effect will be 25 nm/s^2
    enddo
! bias (50 nm/s^2)
    orbit_perturb(10:12)= 0.05d0    ! perturb the accelerometer biases by 50 nm/s^2
! amplitudes of once-per-rev accelerations
    orbit_perturb(13:14) = 20.d-3  ! once-per-rev amplitudes of 20 nm/s^2 (expressed here in um/s^2)
! amplitudes of twice-per-rev accelerations
    orbit_perturb(15:16) = 20.d-3  ! twice-per-rev amplitudes of 20 nm/s^2 (expressed here in um/s^2)
! roll offset error (ie star camera orientation error in the roll direction)
    orbit_perturb(17:19) = 1.d-3      ! 1 mrad orientation error around roll, pitch and yaw axes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         S T A R T - O F F     C O O R D S    
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    k = 4
    do i = 1,3   ! Loop XYZ in initial conditons
!     Compute Efixed IC
      call efix_inert(rot_i2e, rotdot_i2e, rotacc_i2e ,coor, pos_pert_ef, ifacc, .false.) 
!     Now change the i'th coordinate
      pos_pert_ef(i) = pos_pert_ef(i) + orbit_perturb(i)
!     Now convert back to inertial 
      call efix_inert(rot_e2i, rotdot_e2i, rotacc_e2i ,pos_pert_ef, pos_pert_if, ifacc, .false.) 
! and store the coordinates perturbed by the position changes
      do j = 1, 3
        xnew(k) = pos_pert_if(j)
        vnew(k) = pos_pert_if(j+3)
! increment our counter
        k = k + 1
      end do
    end do

! and now for velocities
    do i = 1,3   ! Loop XYZ velocities in initial conditons
!     Compute Efixed IC
      call efix_inert(rot_i2e, rotdot_i2e , rotacc_i2e ,coor, vel_pert_ef, ifacc, .false.) 

!     Now change the i'th coordinate
      vel_pert_ef(i+3) = vel_pert_ef(i+3) + orbit_perturb(i+3)

!     Now convert back  
      call efix_inert(rot_e2i, rotdot_e2i , rotacc_e2i ,vel_pert_ef, vel_pert_if, ifacc, .false.) 

! and store the coordinates perturbed by the velocity changes
      do j = 1, 3
        xnew(k) = vel_pert_if(j)
        vnew(k) = vel_pert_if(j+3)
! increment our counter
        k = k + 1
      end do
    end do

! PT140814: we need to fill out the rest of xnew,vnew (the coords for the 
!           perturbed orbits for bias, scale, 2-per-rev) with some starting
!           coords, i.e. those of the unperturbed orbit
    do k=21,1+norbs_perturbed*3,3  ! this covers only bias and scale at present ! PT170405: no, it covers all perturbations, surely!
      do j=1,3
        xnew(k+j) = xnew(j)
        vnew(k+j) = vnew(j)
      enddo
    enddo       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Setup output of EFixed coordinates and partials
    efpos = 0   ! Clear the Earthfixed array
    call efix_inert(rot_i2e, rotdot_i2e , rotacc_i2e ,coor, efpos, ifacc, .false.) 
    do i = 2,7  ! Only IC is set to 1 (zero from others partials)
       efpos((i-1)*6+i-1) = 1
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (1.eq.2) then ! if we require partials with respect to inertial ics
       call status_update('FATAL','GRACEORB','ADAM_MOULTON_PART',' ',"Code cannot compute partials wrt inertial ICs",0)
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! PT140904: we have 3 pos, 3 vel, 3 bias, 3 scale, 2 1/rev, 2 2/rev, 3 orientation parameters, unknown number of mascons, unknown tide amplitudes
! PT210813: don't include mascon partials if we don't want them in the GTORB files
    if(gt_mcon(1) == "Y")then
      nd = 20 + total_prim + num_mcon_tide*2  ! multiply the number of tidal constituents by 2 for a sine and a cosine amplitude
    else
      nd = 20 
    endif
!
!   END TEST Setup   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    tnew=time
    timeout = nint(dble(jd)*86400.d0+tnew-2451545*86400.d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   write to output files
!   The values on each line of the GTORB file are, in order,: Epoch, 6x pos/vel, 4x quaternions, 36x partials wrt ICs, 18x partials wrt acc
!   bias 18x partials wrt scls, 6x partials wrt mascon 
    nrec=nrec+1
!   Write out the initial IC value and partials
    if (ish5) then
        call writeslice(orb_file,  timeout, efpos, quatSRF2e_init, total_prim, num_mcon_tide, gt_mcon(1))
    else
        write(15,rec=nrec) timeout, efpos(1:6), quatSRF2e_init(0:3), efpos(7:nd*6)
    endif
! print*,timeout, efpos(1:6), quatSRF2e_init(0:3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Initialize values of previous timesteps
!   j loops over the x,y,z components of position(x) and velocity(v)
!     and partials of x and v with respect to 6-IC. (3*nd)
!   i loops over the time steps for integration
    x = 0.d0
    v = 0.d0
    a = 0.d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   do three steps of Runge Kutta
    do istep=1,3

!     shuffle past data
      do j = 1, 3*nd
        do i=3,1,-1
          x(j,i)=x(j,i-1)
          v(j,i)=v(j,i-1)
          a(j,i)=a(j,i-1) !accelern
        enddo
        x(j,0)=xnew(j)
        v(j,0)=vnew(j)
      enddo
      tn = tnew

!     Make evaluations for this time step
!      f(0)=func(xn,y(0)) + partials terms
      do j = 1, 3 * nd
        x0(j) = x(j,0)
        v0(j) = v(j,0)
      enddo

! PT140318: compute new intertial<->efixed rotation matrices for time tn

      output_tide = .false.
      mjd = jd + tn/86400.d0 - 2400000.5d0
      call inert_interp(mjd,nepochs,rot_dates,rot_mats,rotdots,rotaccs,rot_i2e,rotdot_i2e &
                      ,rot_e2i,rotdot_e2i,rotacc_e2i,rotacc_i2e)
      call gravcalc(jd, tn, x0(1:3), v0(1:3), acc_sat, efpos(1:3), output_tide)
      acc(1:3) = acc_sat(1:3)
!     Now add the central force partials to acceleration terms
      call ICpart_calc( maxnd, jd, tn, x0, v0, efpos(1:3), acc) 
! print*,'debug 1: acc for tides after ICpart_calc',acc(13720:13731)

      flagrun = 1   

      do j = 1, 3*nd 
        a(j,0) = acc(j)  
      enddo

! PT140812: to integrate the perturbed orbit we need to compute the difference in the static gravity field between perturbed and satellite location
      accel_pert_if(1:3) = a(1:3,0)
! PT140813: the variable satpos_perturbed that is passed in here contains:
!   the actual IC starting point of the satellite satpos_perturbed(1:3)
!   location when perturbed by +0.5 m in X  satpos_perturbed(4:6)
!   location when perturbed by +0.5 m in Y  satpos_perturbed(7:9)
!   location when perturbed by +0.5 m in Z  satpos_perturbed(10:12)  etc
!  call system ("date +%N")
      satpos_perturbed_if(1:3+norbs_perturbed*3) = x0(1:3+norbs_perturbed*3)
      call calc_perturbed_acc(norbs_perturbed,onepr,twopr,orbit_perturb,satpos_perturbed_if,accel_pert_if,jd,tn)

! now replace the accelerations 4:12 with those computed at the perturbed orbit locations
!  print*,'a sat_acc diff',accel_pert_if(1:3) - a(1:3,0),tn
      a(1:3+norbs_perturbed*3,0) = accel_pert_if(1:3+norbs_perturbed*3)

! DEBUG
! print*,'debug 2: acc for tides after ICpart_calc',acc(13720:13731)

!     Actual  Runge-Kutta step
!     j  loops through the 3 components of position(x) and velocity(v)

!     k1 calculation
      do j = 1, 3*nd
        x1(j) = x(j,0) 
        v1(j) = v(j,0) 
      enddo

      output_tide = .false.
      call gravcalc(jd, tn, x1(1:3),v1(1:3), acc_sat, efpos(1:3), output_tide)
      acc(1:3) = acc_sat(1:3)
      call ICpart_calc( maxnd, jd, tn, x1, v1, efpos(1:3), acc) 

      do j = 1, 3*nd
        a1(j) = acc(j)
      enddo
! PT140812: to integrate the perturbed orbit we need to compute the difference in the static gravity field between perturbed and satellite location
      accel_pert_if(1:3) = a1(1:3)
      satpos_perturbed_if(1:3+norbs_perturbed*3) = x1(1:3+norbs_perturbed*3)
      call calc_perturbed_acc(norbs_perturbed,onepr,twopr,orbit_perturb,satpos_perturbed_if,accel_pert_if,jd,tn)
!  print*,'a1 sat_acc diff',accel_pert_if(1:3) - a1(1:3),tn
      a1(1:3+norbs_perturbed*3) = accel_pert_if(1:3+norbs_perturbed*3)


!     k2 calculation
      do j = 1, 3*nd
        x2(j) = x(j,0)+0.5d0*v1(j)*h
        v2(j) = v(j,0)+0.5d0*a1(j)*h
      enddo

      output_tide = .false.
      mjd = jd + (tn+0.5d0*h)/86400.d0 - 2400000.5d0
      call inert_interp(mjd,nepochs,rot_dates,rot_mats,rotdots,rotaccs,rot_i2e,rotdot_i2e &
                      ,rot_e2i,rotdot_e2i,rotacc_e2i,rotacc_i2e)
      call gravcalc(jd, tn+0.5d0*h, x2, v2, acc, efpos(1:3),output_tide)
      call ICpart_calc( maxnd, jd, tn+0.5d0*h, x2, v2, efpos(1:3), acc) 
      do j = 1, 3*nd
        a2(j) = acc(j)
      enddo
! PT140812: to integrate the perturbed orbit we need to compute the difference in the static gravity field between perturbed and satellite location
      accel_pert_if(1:3) = a2(1:3)
      satpos_perturbed_if(1:3+norbs_perturbed*3) = x2(1:3+norbs_perturbed*3)
      call calc_perturbed_acc(norbs_perturbed,onepr,twopr,orbit_perturb,satpos_perturbed_if,accel_pert_if,jd,tn+0.5d0*h )
!  print*,'a2 sat_acc diff',accel_pert_if(1:3) - a2(1:3),tn+0.5d0*h
      a2(1:3+norbs_perturbed*3) = accel_pert_if(1:3+norbs_perturbed*3)

!     k3 calculation
      do j = 1, 3*nd
        x3(j) = x(j,0)+0.5d0*v2(j)*h
        v3(j) = v(j,0)+0.5d0*a2(j)*h
      enddo

      output_tide = .false.
      call gravcalc(jd, tn+0.5d0*h, x3, v3, acc, efpos(1:3), output_tide)
      call ICpart_calc( maxnd, jd, tn+0.5d0*h, x3, v2, efpos(1:3), acc) 

      do j = 1, 3*nd
        a3(j) = acc(j)
      enddo
! PT140812: to integrate the perturbed orbit we need to compute the difference in the static gravity field between perturbed and satellite location
      accel_pert_if(1:3) = a3(1:3)
      satpos_perturbed_if(1:3+norbs_perturbed*3) = x3(1:3+norbs_perturbed*3)
      call calc_perturbed_acc(norbs_perturbed,onepr,twopr,orbit_perturb,satpos_perturbed_if,accel_pert_if,jd,tn+0.5d0*h )
!  print*,'a3 sat_acc diff',accel_pert_if(1:3) - a3(1:3),tn+0.5d0*h
      a3(1:3+norbs_perturbed*3) = accel_pert_if(1:3+norbs_perturbed*3)

!     k4 calculation
      do j = 1, 3*nd
        x4(j) = x(j,0)+v3(j)*h
        v4(j) = v(j,0)+a3(j)*h
      enddo

      output_tide = .false.
      mjd = jd + (tn+h)/86400.d0 - 2400000.5d0
      call inert_interp(mjd,nepochs,rot_dates,rot_mats,rotdots,rotaccs,rot_i2e,rotdot_i2e &
                      ,rot_e2i,rotdot_e2i,rotacc_e2i,rotacc_i2e)
      call gravcalc(jd, tn+h, x4, v4, acc, efpos(1:3), output_tide)
      call ICpart_calc( maxnd, jd, tn+h, x4, v4, efpos(1:3), acc) 

      do j = 1, 3*nd
        a4(j) = acc(j)
      enddo
! PT140812: to integrate the perturbed orbit we need to compute the difference in the static gravity field between perturbed and satellite location
      accel_pert_if(1:3) = a4(1:3)
      satpos_perturbed_if(1:3+norbs_perturbed*3) = x4(1:3+norbs_perturbed*3)
      call calc_perturbed_acc(norbs_perturbed,onepr,twopr,orbit_perturb,satpos_perturbed_if,accel_pert_if,jd,tn+h )
!  print*,'a4 sat_acc diff',accel_pert_if(1:3) - a4(1:3),tn+h
      a4(1:3+norbs_perturbed*3) = accel_pert_if(1:3+norbs_perturbed*3)

!     end of k loops

!     Values at the end of time step(h)
      do j = 1, 3*nd
        xnew(j)=x(j,0)+h*(v1(j)+2.d0*v2(j)+2.d0*v3(j)+v4(j))/6.d0
        vnew(j)=v(j,0)+h*(a1(j)+2.d0*a2(j)+2.d0*a3(j)+a4(j))/6.d0
      enddo
      tnew=tn+h

      do k=1,nd
        do j = 1,3  ! loop over xyz
          incoor((k-1)*6  +j) = xnew((k-1)*3+j)
          incoor((k-1)*6+3+j) = vnew((k-1)*3+j)
        enddo
      enddo

! PT140813: now update our perturbed satellite positions vector
      satpos_perturbed_if(1:3+norbs_perturbed*3) = xnew(1:3+norbs_perturbed*3)

!     Now loop over Orbit and partials and convert to earthfixed.
      sflag = 'i2e'
      do k = 1, nd      
        call efix_inert(rot_i2e, rotdot_i2e , rotacc_i2e ,incoor((k-1)*6+1), efpos((k-1)*6+1), ifacc, .false.) 
      end do

      timeout = nint(dble(jd)*86400.d0+tnew-2451545*86400.d0)

! call subroutine to calculate the quaternion for the inertial to efix rotation
! matrix the matrix 'rot' is taken from the last call of rotsnp which is in the
! above call of efix_inert
!      call rotmat2quat(rot, quati2e)
      call rotation_mat2quat_3d(rot_i2e, quati2e)

! we have the SCA quaternion values quat(0:3) from the last call of accelerom
! in gravcalc combine the two quaternions (operate with quat (for SRF2i) then
! quati2e)
      call quat_mul(quati2e, quat, quatSRF2e)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  ************  write out results of Runge Kutta steps   ***********
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      nrec=nrec+1
! PT140813: for the position partials, we compute the difference in location of the actual
!           orbit and the perturbed orbit, then divide by the perturbation
! print*,'t,satpos',timeout,efpos(1:6)
! print*,'t,sclxpos',timeout,efpos(37:42)
      do k=1,norbs_perturbed
        do j=1,6
          efpos(k*6+j) = (efpos(k*6+j)-efpos(j))/orbit_perturb(k)
        enddo
      enddo
      if(ish5) then
        call writeslice(orb_file,  timeout, efpos, quatSRF2e(0:3), total_prim, num_mcon_tide, gt_mcon(1))
      else
        write(15,rec=nrec) timeout, efpos(1:6), quatSRF2e(0:3), efpos(7:nd*6)
      endif
! DEBUG:
!print*,'num_mcon_tide=',num_mcon_tide
!print*,'bottom of Runga Kutter',efpos(nd*6-num_mcon_tide*2*6+1:nd*6)
    enddo  ! end of three steps of Runga Kutter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Adam-moulton integrator uses the output of the Runge Kutta
!   Current y for Adams-Moulton
    do j = 1, 3*nd
      xam(j)=xnew(j)
      vam(j)=vnew(j)
    enddo

    do istep=4,nsteps
!SA 200130 Preprocessor flags to remove this logging information in production mode
#ifndef _NDEBUG
      if( mod(istep,1000).eq.0 ) then
        write(message,'(a,2i11)')'Time Step ',timeout, istep
        call status_update('STATUS','GRACEORB','adam_moulton_part',' ',message,0)
      endif
#endif

! Set up for next time step by shifting stored values of previous time
! steps to keep only the results of the three previous steps
      do j = 1, 3*nd 
        do i=3,1,-1
          x(j,i)=x(j,i-1)
          v(j,i)=v(j,i-1)
          a(j,i)=a(j,i-1)
        enddo 
        x(j,0)=xam(j)
        v(j,0)=vam(j)
      enddo
      output_tide = .false.
      mjd = jd + (tnew)/86400.d0 - 2400000.5d0
      call inert_interp(mjd,nepochs,rot_dates,rot_mats,rotdots,rotaccs,rot_i2e,rotdot_i2e &
                      ,rot_e2i,rotdot_e2i,rotacc_e2i,rotacc_i2e)
      call gravcalc(jd, tnew, x, v, acc, efpos(1:3), output_tide)
      call ICpart_calc( maxnd, jd, tnew, x, v, efpos(1:3), acc) 
!call printmat(rot_i2e,3,3,'adam_1_roti2e ')
!call printmat(rotdot_i2e,3,3,'adam_1_rotdoti2e ')

      do j = 1, 3*nd
        a(j,0) = acc(j)
      enddo
! PT140812: to integrate the perturbed orbit we need to compute the difference in the static gravity field between perturbed and satellite location
      accel_pert_if(1:3) = a(1:3,0)
      satpos_perturbed_if(1:3+norbs_perturbed*3) = x(1:3+norbs_perturbed*3,0)
      call calc_perturbed_acc(norbs_perturbed,onepr,twopr,orbit_perturb,satpos_perturbed_if,accel_pert_if,jd,tnew )
!  print*,'a  sat_acc diff',accel_pert_if(1:3) - a(1:3,0),tnew
      a(1:3+norbs_perturbed*3,0) = accel_pert_if(1:3+norbs_perturbed*3)

      tn = tnew
      tnew = tn+h

! Do an Adams-Bashford predictor for Adams-Moulton

      do j = 1, 3*nd
        xabt(j)=xam(j)+(55.d0*v(j,0)-59.d0*v(j,1)+37.d0*v(j,2)-9.d0*v(j,3))*h/24.d0
        vabt(j)=vam(j)+(55.d0*a(j,0)-59.d0*a(j,1)+37.d0*a(j,2)-9.d0*a(j,3))*h/24.d0
      enddo

      output_tide = .true.
      mjd = jd + (tnew)/86400.d0 - 2400000.5d0
      call inert_interp(mjd,nepochs,rot_dates,rot_mats,rotdots,rotaccs,rot_i2e,rotdot_i2e &
                      ,rot_e2i,rotdot_e2i,rotacc_e2i,rotacc_i2e)
      call gravcalc(jd, tnew, xabt, vabt, acc, efpos(1:3), output_tide)
      call ICpart_calc(maxnd, jd, tnew, xabt, vabt, efpos(1:3), acc) 

! PT191030: the call to gravcalc above means that accelerom.f90 was called. The value of "quat", therefore, gives us the coorect value
!           needed to calculate the SRF-to-EF quaternion. Save it.
      quat_SCA(0:3) = quat(0:3)

      do j = 1,3*nd
        anew(j) = acc(j)
      enddo
! PT140812: to integrate the perturbed orbit we need to compute the difference in the static gravity field between perturbed and satellite location
      accel_pert_if(1:3) = anew(1:3)
      satpos_perturbed_if(1:3+norbs_perturbed*3) = xabt(1:3+norbs_perturbed*3)
      call calc_perturbed_acc(norbs_perturbed,onepr,twopr,orbit_perturb,satpos_perturbed_if,accel_pert_if,jd,tnew )
!  print*,'anew sat_acc diff',accel_pert_if(1:3) - anew(1:3),tnew
      anew(1:3+norbs_perturbed*3) = accel_pert_if(1:3+norbs_perturbed*3)

! Now apply the Adams-Moulton correction step
      do j = 1, 3*nd
        xam(j)=xam(j)+(9.d0*vabt(j)+19.d0*v(j,0)-5.d0*v(j,1)+v(j,2))*h/24.d0
        vam(j)=vam(j)+(9.d0*anew(j)+19.d0*a(j,0)-5.d0*a(j,1)+a(j,2))*h/24.d0
      enddo

! PT170821: is this correct for the acceleration at this epoch?
      do j=1,3
!        ifacc(j) = (9.d0*anew(j)+19.d0*a(j,0)-5.d0*a(j,1)+a(j,2))/24.d0   ! the h/24 comes from above. The division by 5 is to make it m/s^2 rather than per 5 sec ... ?
         ifacc(j) = a(j,0)   ! l'avis de Seb
      enddo

      do k=1,nd
        do j = 1,3  ! loop over xyz
          incoor((k-1)*6  +j) = xam((k-1)*3+j)
          incoor((k-1)*6+3+j) = vam((k-1)*3+j)
        enddo
      enddo

! PT140813: now update our perturbed satellite positions vector
      satpos_perturbed_if(1:3+norbs_perturbed*3) = xam(1:3+norbs_perturbed*3)

! Now loop over Orbit and partials and convert to earthfixed.
       do k = 1, nd      
        call efix_inert(rot_i2e, rotdot_i2e, rotacc_i2e, incoor((k-1)*6+1), efpos((k-1)*6+1), ifacc, .false. )
      end do

! PT170818: to compute efixed accelerations, convert to efixed here
      call efix_inert(rot_i2e, rotdot_i2e, rotacc_i2e, incoor(1:6), efacc(1:3), ifacc, .true. )

! call subroutine to calculate the quaternion for the inertial to efix rotation
! matrix 
      call rotation_mat2quat_3d(rot_i2e, quati2e)

! we have the SCA quaternion values quat(0:3) from the last call of
! accelerom in gravcalc combine the two quaternions (operate with quat
! (for SRF2i) then quati2e)
! PT191030: update this to use the saved value "quat_SCA", which should be the correct value to use
      call quat_mul(quati2e, quat_SCA, quatSRF2e)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  ************  write out results of Runge Kutta steps   ***********
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      timeout = nint(dble(jd)*86400.d0+tnew-2451545*86400.d0)
      nrec=nrec+1
! PT140813: for the position partials, we compute the difference in location of the actual
!           orbit and the perturbed orbit, then divide by the perturbation
      do k=1,norbs_perturbed
        do j=1,6
          efpos(k*6+j) = (efpos(k*6+j)-efpos(j))/orbit_perturb(k)
        enddo
      enddo

      if (ish5) then
        call writeslice(orb_file,  timeout, efpos, quatSRF2e(0:3), total_prim, num_mcon_tide, gt_mcon(1))
      else
        write(15,rec=nrec) timeout, efpos(1:6), quatSRF2e(0:3), efpos(7:nd*6)
      endif
!! DEBUG
!! PT170818: print out the epoch, pos,vec and acc
!! earth-fixed (the acceleration is still wrong)
!      write(*,*)timeout, efpos(1:6),efacc(1:3)," posvelacc efixed"
!! inertial
!      write(*,*)timeout,incoor(1:6),ifacc," posvelacc inert"
!
!  quaternion info
!print*,"quat_actual_ef:",1+(timeout-252417600)/5,quatSRF2e(0:3),timeout
!print*,"quat_SCA      :",1+(timeout-252417600)/5,quat_SCA,timeout
    enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    close(11)
    return
    end subroutine adam_moulton_part



