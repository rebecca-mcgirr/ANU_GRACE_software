   subroutine atm_tide_sph(jd, tin, rlat,rlong, rx, atm_tide_acc )  

! subroutine to calculate the acceleration caused by the atmospheric tides.
! Adapted from oceanfield.f90 that was written by Emma-Kate Potter, 13 Oct, 2010
!
! P. Tregoning
! 15 October 2016

! This subroutine 
! 1) opens a (hardwired) file that contains a spherical harmonic representation of the atmospheric tide pressures
!    (all at degree 180) every 10 minutes.
! 2) performs a linear interpolation between the two 10-min epochs to get a set of coefficients relevant for
!    the input time (jd, tin)
! 

    use gm_mod
    use gauleg_mod
    use lovenum_mod
    use coeff_mod
    use spherhar_mod
    use sat_mod         ! provides satellite name (A or B)

    implicit none

    integer*4 :: i, xxn, yyn
    integer*4 :: jd
    integer :: ideg, iord, n ,m, j ,nm

    real*8 :: tin
    real*8, dimension(nlon, nlat) :: otideh
    real*8, dimension(maxnml) :: otidecoef
    real*8 :: xx, yy, func, scale_factor
    real*8 :: klovep1 
    real*8 :: tempC, tempC2 
    double precision rlat, rlong, rx, sum_oheight

    integer, dimension(3) :: tarray
    character message*250
    logical output_groundtrack_tide

! PT160927: arrays of coefficients
    integer maxdeg_to_use
    parameter(maxdeg_to_use=180)

! PT160928: make a 3D array of all the coefficients for all epochs
    integer*4,parameter    ::  n_atm_tide_epochs=145
    real(kind=8),allocatable :: dCatm_tide_all(:,:,:)
    real(kind=8),allocatable :: dSatm_tide_all(:,:,:)

    real(kind=8) :: tmpC,tmpS
    integer*4    :: tmpdeg,tmpord,tide_1,tide_2,tide_epoch

! PT160727: counters
    integer :: counter, num_tide_epochs

! tidal acceleration computations
    real(kind=8) :: atm_tide_acc(3),maxdeg_atm_tide

    save  num_tide_epochs,tide_1,tide_2,dCatm_tide_all,dSatm_tide_all

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT160927: read the coefficients from the input (hardwired) spherical harmonic file

! SA211221: fix bug here if using 1 sec integration step
!    if(int(tin) == 43200)then
    if( .not. allocated(dCatm_tide_all)) then

! allocate the arrays
      allocate(dCatm_tide_all(0:maxdeg_to_use,0:maxdeg_to_use,n_atm_tide_epochs))
      allocate(dSatm_tide_all(0:maxdeg_to_use,0:maxdeg_to_use,n_atm_tide_epochs))

      ! open the file
      call status_update('STATUS','GRACEORB','atm_tide_sph',"atm_tide_stokes.hs",'opening atmospheric tide file ',0)
      open(125,file="atm_tide_stokes.hs",status='old')

! read all the epochs. There are 6 per hour x 24 hours = 144 epochs + 1
      read(125,*)maxdeg_atm_tide
      do tide_epoch = 1,n_atm_tide_epochs
      ! read the 00h set of coefficients
        do counter = 0,(maxdeg_to_use+1)*(maxdeg_to_use+2)/2 -1
          read(125,*)tmpdeg,tmpord,tmpC,tmpS
          dCatm_tide_all(tmpdeg,tmpord,tide_epoch) = tmpC
          dSatm_tide_all(tmpdeg,tmpord,tide_epoch) = tmpS
!if(tmpdeg < 4 .or. (tmpdeg == 180 .and.tmpord > 178))print*,"atm tide",tide_epoch,counter,tmpdeg,tmpord,tmpC,tmpS
        enddo
      enddo
      call status_update('STATUS','GRACEORB','atm_tide_sph',"atm_tide_stokes.hs",'have read all epochs in atm tide file',0)

! assign pointers to the harmonic grids
        tide_1 = 1
        tide_2 = 2
        num_tide_epochs = 1
    else                     ! just another epoch. Increment the epoch counter
        num_tide_epochs = num_tide_epochs + 1
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT160927: do we need to update our two epochs of spherical harmonic coefficients?
    if(tin > 43200 .and. mod(tin,600.) == 0)then
!      print*,'update the atm tidal spherical harmonics tin,num_tide_epochs',tin,num_tide_epochs
      tide_1 = tide_2
      tide_2 = tide_2 + 1
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT160927: now interpolate the two sets of coefficients to get values for this epoch
!           dCatide and dSatde are declared in coeff_mod.f90
    do ideg = 0,maxdeg_to_use
      do iord = 0,ideg
        dCatide(ideg,iord) = dCatm_tide_all(ideg,iord,tide_1) &
                                +(dCatm_tide_all(ideg,iord,tide_2)-dCatm_tide_all(ideg,iord,tide_1))*mod(tin,600.)/600.
        dSatide(ideg,iord) = dSatm_tide_all(ideg,iord,tide_1) &
                                +(dSatm_tide_all(ideg,iord,tide_2)-dSatm_tide_all(ideg,iord,tide_1))*mod(tin,600.)/600.
      enddo
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate the associated Legendre functions
  call legcalc_norm(sin(rlat))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT160927: limit the degree of the tide to whatever test we want
!    dCotide(181:maxsize,181:maxsize) = 0.0
!    dSotide(181:maxsize,181:maxsize) = 0.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT160929: finally, compute the tidal accelerations from the spherical harmonic coefficients for this epoch
! define the factor to convert from EWH to Stokes' coefficients (the division by 2n+1 is done in the subroutine)
  call compute_sph_acc(rlat,rlong,rx,dCatide,dSatide,maxsize,Plm,PDlm,"stokes",atm_tide_acc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT130906: calculate and store the atmospheric tide height at the satellite location for this epoch
! DEBUG
  output_groundtrack_tide = .true.
  if(output_groundtrack_tide .and. mod(tin,300.0) == 0 )then
    sum_oheight = 0.d0
    nm = 4
    do ideg = 2,maxdeg_to_use        ! limit to 200 until the graceorb code consistently goes higher
      sum_oheight = sum_oheight+Plm(nm)*dCatide(ideg,0)         ! * cos(m*lambda) = 1

      nm = nm +1

      do iord = 1, ideg    
        sum_oheight =  sum_oheight + Plm(nm)*(dCatide(ideg,iord)* &
                dcos(iord*rlong)+dSatide(ideg,iord)*dsin(iord*rlong) ) 

        nm = nm + 1
      enddo
    enddo
!    write(message,'(a,i6,a,f8.4,a,2f11.5,a)')'   Seconds of day:',int(tin-43200.0),' Tide ht (from sph har),',sum_oheight  &
!           ,' m. Coords (',rlat*180.d0/pi,rlong*180.d0/pi, ' )'
        write(message,'(a,a,a,i6,a,f8.4,a,2f11.5,a,3f10.3)')'GRACE ',sat,' Seconds of day:',int(tin-43200.)&
             ,' Atm Tide acc (from sph har),',sum_oheight  &
             ,' m. Coords (',rlat*180.d0/pi,rlong*180.d0/pi, ' )',atm_tide_acc*1.d9
!    call status_update('STATUS','GRACEORB','atm_tide_sph',' ',message,0)
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    return
    end

