!********************************************************************************************************************************
!  File: plot_lib.f90
!
!  Purpose: Set of subroutines
!
!  Author: Thomas Greenspan
!          (subroutines taken from parts of plt_postfit.f90 by Simon McClusky, modified by R. King)
!
!  API:
!       plot_writeKB  : Writes out kb resid to plot kb file
!       plot_writePLT : Writes out residuals and their maxima to plot sat files
!
!  July 31, 2013
!
!********************************************************************************************************************************
!********************************************************************************************************************************
! plot_writeKB: Calculates and writes out to the plot KB file information on the kband
!               postfit residuals
! Author: Thomas Greenspan
!         (taken from part of plt_postift.f90 by Simon McClusky)
!
! MODS
! PT201113: pass in postfit unfiltered kbrr residuals and output as an additional column at the end of each line
!********************************************************************************************************************************

subroutine plot_writeKB(iop,pre_omc,post_omc,post_omc_kbrr_unfilt,rvec,rpy,kbrr_prefit_tol,last_epoch)
  use gracefit_mod
  implicit none

   
  include 'output.h90'

  !********************  Variable declarations ********************

  integer*4, intent(in)        :: iop                         ! Plot resid indicator
  double precision, intent(in) :: pre_omc(nobs_t,nepochs_t)   ! Prefit Observed Minus Computed values
  double precision, intent(in) :: post_omc(nobs_t,nepochs_t)  ! Postfit Observed Minus Computed values
  double precision, intent(in) :: post_omc_kbrr_unfilt(nepochs_t)  ! Postfit OMC unfiltered KBRR residuals
  double precision, intent(in) :: rvec(3,nsat_t,nepochs_t)    ! Vector of positions
  double precision, intent(in) :: rpy(3,nsat_t,nepochs_t)     ! Roll, pitch and yaw information for both satellites at every epoch
  double precision, intent(in) :: kbrr_prefit_tol             ! KBRR prefit tolerance
  integer*4,        intent(in) :: last_epoch                  ! last epoch that was included in the analysis

  integer*4        :: iepoch                     ! Counter that runs through epochs
  integer*4        :: isat                       ! Counter that runs through satellites
  integer*4        :: i                          ! Counter variables
  logical          :: kbrr_outliers(nepochs_t)   ! Indicates whether the kbrr postfit values at a given epoch are outliers
  double precision :: orbllh(3,nsat_t,nepochs_t) ! Colatiude, longitude and radius of each satellite at every epoch
  double precision :: maxdrr                     ! Maximum value of the Range Rate omc (over every epoch)
  double precision :: rotmat(3,3)                ! Rotation matrix needed as input for XYZ_to_GEOD (not used otherwise)
  character(256)   :: message                    ! Message to be written out concerning maximum kbrr resid
  !****************************************************************

  !************************* WRITE HEADER *************************

  write(LUGN_KB,'(a,i7)') '#Number of epochs: ',nepochs_t
  write(LUGN_KB,'(a,f6.2,a)') '#(Time Interval =',epoch_interval,' sec)'
  !****************************************************************

  ! Check for kbrr outliers
  call kb_kbrrOutliers(post_omc,kbrr_prefit_tol,kbrr_outliers)

  ! If ground-track values needed, compute them now
  maxdrr = 0.d0
  if( iop == 4 )then
     do iepoch = 1, last_epoch      !nepochs_t
        if(kbrr_outliers(iepoch)) cycle  !@# Do we need the condition here as well?
        do isat = 1, nsat_t
           !       Convert X,Y,Z to Lat, Long, Height, and delta X,Y,Z, to delta N,E,U...
           call XYZ_to_GEOD(rotmat, rvec(1:3,isat,iepoch)*m_mm, orbllh(:,isat,iepoch))
           !         Convert Lat, Long, Height from radians/meters to degrees/kilometres.
           orbllh(1,isat,iepoch) = 90.d0 - orbllh(1,isat,iepoch)*rad_deg
           orbllh(2,isat,iepoch) = orbllh(2,isat,iepoch)*rad_deg
           orbllh(3,isat,iepoch) = orbllh(3,isat,iepoch)*m_km
        enddo  ! End of Satellite loop
        !     Find the maximum value for the post_omc range rate (not taking into account first and last 5 observations)
        if(iepoch.gt.5.and.iepoch.lt.(nepochs_t-5))then
           if (dabs(post_omc(ikbrr,iepoch)) > maxdrr ) maxdrr = dabs((post_omc(ikbrr,iepoch)))
        endif
     enddo  ! End of epoch loop
  endif  ! End of "ground-track values" if statement

  !********************* WRITE TO PLOT KB FILE ********************

  write(LUGN_KB,'(a,f10.6,a)') '#Maximum range rate residual:  ',maxdrr * m_mm,' (mm)'
  write(LUGN_KB,'(a)') '# '
  write(LUGN_KB,'(a,a)') '# iepoch  postfit RR (mm)    lat (deg)      lon (deg)      Az (m)     ', &
       '  Roll A   Pitch A  Yaw A     Roll B   Pitch B  Yaw B    prefit RR (mm)' 
  write(LUGN_KB,'(a)') '# '

  do iepoch = 1,last_epoch    !nepochs_t
     if(kbrr_outliers(iepoch)) cycle ! Only print to file "good" residuals
     ! PT150810: also write out the postfit position/velocity
     ! PT190509: increase the postfit pos/vel to f18.10 from f15.10
     write(LUGN_KB,'(i5,1x,f15.8,2x,2f12.5,f15.5,2(3f15.4,1x),f15.10,12f18.10,12f18.10,2f18.10)') & 
          iepoch, post_omc(ikbrr,iepoch)*m_mm, &
          (orbllh(i,1,iepoch),i=1,3),((rpy(i,isat,iepoch)*1.d3,i=1,3),isat=1,nsat_t),pre_omc(ikbrr,iepoch)*m_mm &
          , pre_omc(1:12,iepoch),post_omc(1:12,iepoch),post_omc_kbrr_unfilt(iepoch)*m_mm &
          , (post_omc(ikbrr,iepoch)-post_omc_kbrr_unfilt(iepoch))*m_mm
  enddo
  !****************************************************************

  write(message,'(a,f15.5,a)') '   Maximum range rate  residual   = ',maxdrr*m_um,' (um/sec)'
  call status_update('STATUS','GRACEFIT','plot_writeKB',' ', message,0)

  close(LUGN_KB)

  return
end subroutine plot_writeKB
!********************************************************************************************************************************
! plot_writeKBA: Calculates and writes out to the plot KB file information on the kband
!               postfit acc residuals
! Author: Sebastien Allgeyer
!********************************************************************************************************************************

subroutine plot_writeKBA(iop,pre_omc,post_omc,post_omc_kbra_unfilt,rvec,rpy,kbrr_prefit_tol,last_epoch)
  use gracefit_mod
  implicit none

   
  include 'output.h90'

  !********************  Variable declarations ********************

  integer*4, intent(in)        :: iop                         ! Plot resid indicator
  double precision, intent(in) :: pre_omc(nobs_t,nepochs_t)   ! Prefit Observed Minus Computed values
  double precision, intent(in) :: post_omc(nobs_t,nepochs_t)  ! Postfit Observed Minus Computed values
  double precision, intent(in) :: post_omc_kbra_unfilt(nepochs_t)  ! Postfit Observed Minus Computed values
  double precision, intent(in) :: rvec(3,nsat_t,nepochs_t)    ! Vector of positions
  double precision, intent(in) :: rpy(3,nsat_t,nepochs_t)     ! Roll, pitch and yaw information for both satellites at every epoch
  double precision, intent(in) :: kbrr_prefit_tol             ! KBRR prefit tolerance
  integer*4,        intent(in) :: last_epoch                  ! last epoch that was included in the analysis

  integer*4        :: iepoch                     ! Counter that runs through epochs
  integer*4        :: isat                       ! Counter that runs through satellites
  integer*4        :: i                          ! Counter variables
  logical          :: kbrr_outliers(nepochs_t)   ! Indicates whether the kbrr postfit values at a given epoch are outliers
  double precision :: orbllh(3,nsat_t,nepochs_t) ! Colatiude, longitude and radius of each satellite at every epoch
  double precision :: maxdrr                     ! Maximum value of the Range Rate omc (over every epoch)
  double precision :: rotmat(3,3)                ! Rotation matrix needed as input for XYZ_to_GEOD (not used otherwise)
  character(256)   :: message                    ! Message to be written out concerning maximum kbrr resid
  !****************************************************************

  !************************* WRITE HEADER *************************

  write(LUGN_KBA,'(a,i7)') '#Number of epochs: ',nepochs_t
  write(LUGN_KBA,'(a,f6.2,a)') '#(Time Interval =',epoch_interval,' sec)'
  !****************************************************************

  ! Check for kbrr outliers
!  call kb_kbrrOutliers(post_omc,kbrr_prefit_tol,kbrr_outliers)

  ! If ground-track values needed, compute them now
  maxdrr = 0.d0
  if( iop == 4 )then
     do iepoch = 1, last_epoch      !nepochs_t
        if(kbrr_outliers(iepoch)) cycle  !@# Do we need the condition here as well?
        do isat = 1, nsat_t
           !       Convert X,Y,Z to Lat, Long, Height, and delta X,Y,Z, to delta N,E,U...
           call XYZ_to_GEOD(rotmat, rvec(1:3,isat,iepoch)*m_mm, orbllh(:,isat,iepoch))
           !         Convert Lat, Long, Height from radians/meters to degrees/kilometres.
           orbllh(1,isat,iepoch) = 90.d0 - orbllh(1,isat,iepoch)*rad_deg
           orbllh(2,isat,iepoch) = orbllh(2,isat,iepoch)*rad_deg
           orbllh(3,isat,iepoch) = orbllh(3,isat,iepoch)*m_km
        enddo  ! End of Satellite loop
        !     Find the maximum value for the post_omc range rate (not taking into account first and last 5 observations)
        if(iepoch.gt.5.and.iepoch.lt.(nepochs_t-5))then
           if (dabs(post_omc(ikbrr,iepoch)) > maxdrr ) maxdrr = dabs((post_omc(ikbra,iepoch)))
        endif
     enddo  ! End of epoch loop
  endif  ! End of "ground-track values" if statement

  !********************* WRITE TO PLOT KB FILE ********************

  write(LUGN_KBA,'(a,f10.6,a)') '#Maximum range accel residual:  ',maxdrr * m_um*1.e3,' (nm/sec^2)'
  write(LUGN_KBA,'(a)') '# '
  write(LUGN_KBA,'(a,a)') '# iepoch  postfit RR (mm)    lat (deg)      lon (deg)      Az (m)     ', &
       '  Roll A   Pitch A  Yaw A     Roll B   Pitch B  Yaw B    prefit RR (mm)' 
  write(LUGN_KBA,'(a)') '# '

  do iepoch = 1,last_epoch    !nepochs_t
     if(kbrr_outliers(iepoch)) cycle ! Only print to file "good" residuals
     ! PT150810: also write out the postfit position/velocity
     write(LUGN_KBA,'(i5,1x,f15.8,2x,2f12.5,f15.5,2(3f15.4,1x),f15.10,12f15.10,12f15.10,2f15.10)') & 
          iepoch, post_omc(ikbra,iepoch)*m_mm*m_mm, &
          (orbllh(i,1,iepoch),i=1,3),((rpy(i,isat,iepoch)*1.d3,i=1,3),isat=1,nsat_t),pre_omc(ikbra,iepoch)*m_mm*m_mm &
          , pre_omc(1:12,iepoch),post_omc(1:12,iepoch),post_omc_kbra_unfilt(iepoch)*m_mm*m_mm &
          ,(post_omc(ikbra,iepoch)-post_omc_kbra_unfilt(iepoch))*m_mm*m_mm
  enddo
  !****************************************************************

  write(message,'(a,f15.5,a)') '  Maximum range accel residual   = ',maxdrr*m_um*1.e3,' (nm/sec^2)'
  call status_update('STATUS','GRACEFIT','plot_writeKBA',' ', message,0)

  close(LUGN_KBA)

  return
end subroutine plot_writeKBA








!********************************************************************************************************************************
!********************************************************************************************************************************
! plot_writePLT: Calculates postfit residuals and maximum postfit residuals and writes them to the 
!                plot sat files. Input is iop (resid plot indicator) the observed minus calculated
!                array and the vector of positions and velocities
! Author: Thomas Greenspan
!********************************************************************************************************************************

subroutine plot_writePLT(iop,post_omc,rvec,GPS_ind)
  use gracefit_mod
  implicit none

   
  include 'output.h90'

  !********************  Variable declarations ********************

  integer*4, intent(in)        :: iop                         ! Plot resid indicator
  double precision, intent(in) :: post_omc(nobs_t,nepochs_t)  ! Observed Minus Computed values to be used in LS algorithm
  double precision, intent(in) :: rvec(6,nsat_t,nepochs_t)    ! Vector of positions, velocities
  double precision, intent(in) :: GPS_ind(nsat_t,nepochs_t)   ! Indicator as to whether there were observations read in (0.0 if not)

  integer*4        :: iepoch                      ! Counter that runs through epochs
  integer*4        :: isat                        ! Counter that runs through satellites
  integer*4        :: i,j                         ! Counter variables
  double precision :: orbllh(3,nsat_t,nepochs_t)  ! Colatiude, longitude and radius of each satellite at every epoch
  double precision :: rotmat(3,3)                 ! Rotation matrix needed as input for XYZ_to_GEOD (not used otherwise)
  double precision :: drac(3,nsat_t,nepochs_t)    ! Radial residuals
  double precision :: maxdrac(3,nsat_t)           ! Maximum values of radial residuals
  double precision :: dxyz(3,nsat_t,nepochs_t)    ! Temporary variable containing x, y, z postion residuals of a satellite
  double precision :: dxyzvel(3,nsat_t,nepochs_t) ! Temporary variable containing x, y, z postion residuals of a satellite
  double precision :: maxdxyz(3,nsat_t)           ! Maximum value for x, y and z postion residuals for each satellite
  double precision :: dneu(3,nsat_t,nepochs_t)    ! Rotated coordinate residuals
  double precision :: maxdneu(3,nsat_t)           ! Maximum value for rotated coordinate residuals
  double precision :: az_diff(nsat_t,nepochs_t)   ! Azimuth of horizontal orbit difference components
  double precision :: len_diff(nsat_t,nepochs_t)  ! Length of horizontal orbit difference components
  !****************************************************************

  !************************* WRITE HEADER *************************

  do isat = 1, nsat_t
     if (iop == 1) then
        write(LUGN_PLT(isat),'("Delta X (mm)")')
        write(LUGN_PLT(isat),'("Delta Y (mm)")')
        write(LUGN_PLT(isat),'("Delta Z (mm)")')
     endif
     if (iop == 2) then
        write(LUGN_PLT(isat),'("D-radial (mm)")')
        write(LUGN_PLT(isat),'("D-along (mm)")')
        write(LUGN_PLT(isat),'("D-cross (mm)")')
     endif
     if (iop == 3) then
        write(LUGN_PLT(isat),'("V-X (km/s)")')
        write(LUGN_PLT(isat),'("V-Y (km/s)")')
        write(LUGN_PLT(isat),'("V-Z (km/s)")')
     endif
     if (iop == 4) then
        write(LUGN_PLT(isat),'("Latitude  (deg)")')
        write(LUGN_PLT(isat),'("Longitude (deg)")')
        write(LUGN_PLT(isat),'("Height (km)")')
        write(LUGN_PLT(isat),'("Delta N (mm)")')
        write(LUGN_PLT(isat),'("Delta E (mm)")')
        write(LUGN_PLT(isat),'("Delta U (mm)")')
        write(LUGN_PLT(isat),'("VECT AZ (deg)")')
        write(LUGN_PLT(isat),'("VECT LEN (mm)")')
     endif
     write(LUGN_PLT(isat),'("Number epochs (Time Interval=",f6.2," sec)")') epoch_interval
     write(LUGN_PLT(isat),'(a,a1,a,i4)') 'GPS observations used for GRACE ',satnam(isat),': ',nGPSobs_t(1)
     write(LUGN_PLT(isat),'(a,i4)') 'Number of epochs: ',nepochs_t
  enddo  ! End of satellite loop
  !****************************************************************

  !   If ground-track values needed, compute them now
  if( iop == 4 ) then               
     do iepoch = 1, nepochs_t
        do isat = 1, nsat_t
           !       Set dxzy to the appropriate pos x, y, z terms of post_omc (only for readability) 
           dxyz(:,isat,iepoch) = post_omc(1+(isat-1)*ngobs_t:3+(isat-1)*ngobs_t,iepoch)
           dxyzvel(:,isat,iepoch) = post_omc(4+(isat-1)*ngobs_t:6+(isat-1)*ngobs_t,iepoch)
           !       Convert X,Y,Z to Lat, Long, Height, and delta X,Y,Z, to delta N,E,U...
            call XYZ_to_GEOD(rotmat, rvec(1:3,isat,iepoch)*m_mm, orbllh(:,isat,iepoch))
           !       Convert Lat, Long, Height from radians/meters to degrees/kilometres.
           orbllh(1,isat,iepoch) = 90.d0 - orbllh(1,isat,iepoch)*rad_deg
           orbllh(2,isat,iepoch) = orbllh(2,isat,iepoch)*rad_deg
           orbllh(3,isat,iepoch) = orbllh(3,isat,iepoch)*m_km
           !       Compute Azimuth and Length of horizontal orbit difference components......
           az_diff(isat,iepoch) = datan2(dneu(2,isat,iepoch),dneu(1,isat,iepoch))
           az_diff(isat,iepoch) = az_diff(isat,iepoch)*rad_deg
           if ( az_diff(isat,iepoch) < 0.d0) az_diff(isat,iepoch) = 360.d0+az_diff(isat,iepoch)
           len_diff(isat,iepoch) = dsqrt(dneu(2,isat,iepoch)**2+dneu(1,isat,iepoch)**2)
        enddo  ! End of Satellite loop
     enddo  ! End of epoch loop
  endif  ! End of "ground-track values" if statement

  ! Find the maximum value for each component
  maxdxyz=0.d0
  maxdrac=0.d0
  maxdneu=0.d0
  do iepoch = 1, nepochs_t
     do isat= 1, nsat_t
        if (mod(iepoch-1,gnv_interval) == 0 .and. GPS_ind(isat,iepoch) /= 0.0) then   ! Only apply when GPS observations were used
           !       Get radial, along-track and cross-track (set drac)
           call xyz2rac(rvec(1:6,isat,iepoch),post_omc(1:3,iepoch),drac(1:3,isat,iepoch))
           do i = 1, 3
              if(dabs((post_omc(i+(isat-1)*ngobs_t,iepoch))) > maxdxyz(i,isat)) maxdxyz(i,isat)&
                   = dabs((post_omc(i+(isat-1)*ngobs_t,iepoch)))
              if(dabs((drac(i,isat,iepoch))) > maxdrac(i,isat)) maxdrac(i,isat) = dabs((drac(i,isat,iepoch)))
              if(iop == 4 )then
                 if (dabs((dneu(i,isat,iepoch))) > maxdneu(i,isat)) maxdneu(i,isat) = dabs((dneu(i,isat,iepoch)))
              endif
           enddo
        endif  ! End of "mod(gnv_interval)" if statement
     enddo  ! End of satellite loop
  enddo  ! End of epoch loop
  !****************************************************************

  !******************** WRITE TO PLOT SAT FILE ********************

  !   Loop over all SVs and write the GNV1B plot files
  do isat = 1, nsat_t

     !   Write time series series
     do iepoch = 1, nepochs_t
        if(mod(iepoch-1,gnv_interval) == 0 .and. GPS_ind(isat,iepoch) /= 0.0)then  ! Only apply when GPS observations were used
           if(iop == 1)then
              if(iepoch == 1)then
                 write( LUGN_PLT(isat), '(f6.2)') maxdxyz(1,isat) * m_mm
                 write( LUGN_PLT(isat), '(f6.2)') maxdxyz(2,isat) * m_mm
                 write( LUGN_PLT(isat), '(f6.2)') maxdxyz(3,isat) * m_mm
              endif
              write(LUGN_PLT(isat),'(i4,3f11.5)') iepoch,(post_omc(i+(isat-1)*6,iepoch)*m_mm,i=1,3)
           endif
           if(iop == 2)then
              if(iepoch == 1)then
                 write( LUGN_PLT(isat), '(f6.2)') maxdrac(1,isat) * m_mm
                 write( LUGN_PLT(isat), '(f6.2)') maxdrac(2,isat) * m_mm
                 write( LUGN_PLT(isat), '(f6.2)') maxdrac(3,isat) * m_mm
              endif
              write(LUGN_PLT(isat),'(i4,3f11.5)') iepoch,(drac(i,isat,iepoch)*m_mm,i=1,3)
           endif
           if(iop == 3)then
              call status_update('FATAL','GRACEFIT','grace/plot_writePLT',' ','Velocity plots not yet coded',0)
           endif
           if(iop == 4) then
              if(iepoch == 1)then
                 write( LUGN_PLT(isat), '(f6.2)') maxdneu(1,isat) * m_mm
                 write( LUGN_PLT(isat), '(f6.2)') maxdneu(2,isat) * m_mm
                 write( LUGN_PLT(isat), '(f6.2)') maxdneu(3,isat) * m_mm
              endif
              write(LUGN_PLT(isat),'(i4,3f20.7,5f14.7)') iepoch,(orbllh(i,isat,iepoch)*m_mm,i=1,3) &
                   ,(dneu(i,isat,iepoch)*m_mm,i=1,3),az_diff(isat,iepoch),len_diff(isat,iepoch)*m_mm 

           endif
        endif  ! End of "mod(gnv_interval)" if statement
     enddo  ! End of epoch loop
  enddo  ! End satellite loop
  !****************************************************************

  do isat = 1, nsat_t
     close(LUGN_PLT(isat))
  enddo

  return
end subroutine plot_writePLT

!********************************************************************************************************************************
