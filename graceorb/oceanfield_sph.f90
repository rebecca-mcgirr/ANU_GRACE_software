   subroutine oceanfield_sph(jd, tin, rlat,rlong, rx, tide_acc )  

! subroutine essentially for debugging/testing of errors introduced by applying the ocean tide as a
! spherical harmonic model rather than through the secondary mascons. Adapted from oceanfield.f90
! that was written by Emma-Kate Potter, 13 Oct, 2010
!
! P. Tregoning
! 27 July 2016

! This subroutine 
! 1) opens a (hardwired) file that contains a spherical harmonic representation of the ocean tide heights
!    (all at degree 360) every 10 minutes.
! 2) performs a linear interpolation between the two 10-min epochs to get a set of coefficients relevant for
!    the input time (jd, tin)
! 
!
! NOTE: this subroutine is NOT meant to be robust. The input file MUST be set up by the user and be appropriate
!       for the day being integrated. It is NOT my responsibility if this goes wrong: long-term this 'S' option
!       for applying the ocean tides is NOT TO BE USED!!!! 

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
    parameter(maxdeg_to_use=200)

! PT160928: make a 3D array of all the coefficients for all epochs
    integer*4,parameter    ::  blah=145
    real(kind=8),allocatable :: dCotide_all(:,:,:)
    real(kind=8),allocatable :: dSotide_all(:,:,:)

    real(kind=8) :: tmpC,tmpS
    integer*4    :: tmpdeg,tmpord,tide_1,tide_2,tide_epoch

! PT160727: counters
    integer :: counter, num_tide_epochs

! tidal acceleration computations
    real(kind=8) :: tide_acc(3)
    real(kind=8) :: factor

    save  num_tide_epochs,tide_1,tide_2,dCotide_all,dSotide_all

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT160927: read the coefficients from the input (hardwired) spherical harmonic file

    if(int(tin) == 43200)then

! allocate the arrays
      allocate(dCotide_all(0:maxsize,0:maxsize,blah))
      allocate(dSotide_all(0:maxsize,0:maxsize,blah))

      ! open the file
      call status_update('STATUS','GRACEORB','oceanfield_sph',' ','opening file sphharm_9sep.hs',0)
        open(124,file="sphharm_9sep.hs",status='old')

! read all the epochs. There are 6 per hour x 24 hours = 144 epochs + 1
      do tide_epoch = 1,blah
      ! read the 00h set of coefficients
        do counter = 0,(maxsize+1)*(maxsize+2)/2 -1
          read(124,*)tmpdeg,tmpord,tmpC,tmpS
          dCotide_all(tmpdeg,tmpord,tide_epoch) = tmpC
          dSotide_all(tmpdeg,tmpord,tide_epoch) = tmpS
        enddo
      enddo
      call status_update('STATUS','GRACEORB','oceanfield_sph',' ','have read all epochs in sphharm_9sep.hs',0)

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
      print*,'update the tidal spherical harmonics tin,num_tide_epochs',tin,num_tide_epochs
      tide_1 = tide_2
      tide_2 = tide_2 + 1
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT160927: now interpolate the two sets of coefficients to get values for this epoch
!print*,'interpolating between tide epochs', tide_1,' and ',tide_2,tin
    do ideg = 0,maxsize
      do iord = 0,ideg
        dCotide(ideg,iord) = dCotide_all(ideg,iord,tide_1) &
                                +(dCotide_all(ideg,iord,tide_2)-dCotide_all(ideg,iord,tide_1))*mod(tin,600.)/600.
        dSotide(ideg,iord) = dSotide_all(ideg,iord,tide_1) &
                                +(dSotide_all(ideg,iord,tide_2)-dSotide_all(ideg,iord,tide_1))*mod(tin,600.)/600.
!if(ideg < 10 .and. iord < 10)then
!  print*,'int(tin-43200.),ideg,iord,mod(tin,600.),C,S',int(tin-43200.),ideg,iord,mod(tin,600.)/600. &
!                                                        ,dCotide(ideg,iord),dSotide(ideg,iord)
!
!  print*,'correction terms:',(dCotide_all(ideg,iord,tide_2)-dCotide_all(ideg,iord,tide_1))*mod(tin,600.)/600. &
!                            ,(dSotide_all(ideg,iord,tide_2)-dSotide_all(ideg,iord,tide_1))*mod(tin,600.)/600.
!  print*,"C temrs",dCotide_all(ideg,iord,tide_2),dCotide_all(ideg,iord,tide_1)
!  print*,"S temrs",dSotide_all(ideg,iord,tide_2),dSotide_all(ideg,iord,tide_1)
!endif
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
  call compute_sph_acc(rlat,rlong,rx,dCotide,dSotide,maxsize,Plm,PDlm,"EWH   ",tide_acc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT130906: calculate and store the ocean tide height at the satellite location for this epoch
! DEBUG
  output_groundtrack_tide = .true.
  if(output_groundtrack_tide .and. mod(tin,1.0) == 0 )then
    sum_oheight = 0.d0
    nm = 4
    do ideg = 2,maxdeg_to_use        ! limit to 200 until the graceorb code consistently goes higher
      sum_oheight = sum_oheight+Plm(nm)*dCotide(ideg,0)         ! * cos(m*lambda) = 1

      nm = nm +1

      do iord = 1, ideg    
        sum_oheight =  sum_oheight + Plm(nm)*(dCotide(ideg,iord)* &
                dcos(iord*rlong)+dSotide(ideg,iord)*dsin(iord*rlong) ) 

        nm = nm + 1
      enddo
    enddo
!    write(message,'(a,i6,a,f8.4,a,2f11.5,a)')'   Seconds of day:',int(tin-43200.0),' Tide ht (from sph har),',sum_oheight  &
!           ,' m. Coords (',rlat*180.d0/pi,rlong*180.d0/pi, ' )'
        write(message,'(a,a,a,i6,a,f8.4,a,2f11.5,a,3f10.3)')'   GRACE ',sat,' Seconds of day:',int(tin-43200.)&
             ,' Tide ht (from sph har),',sum_oheight  &
             ,' m. Coords (',rlat*180.d0/pi,rlong*180.d0/pi, ' )',tide_acc*1.d9
    call status_update('STATUS','GRACEORB','oceantide_sph',' ',message,0)
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   scale the coefficients of coefftide to use to calculate gravitational potential field 
!   so this converts from water height to gravitational potential
!   ****************** this should be tested further **************************
!   APP: This does not seem correct, the deformation love numbers should be scaled by a factor of 4 pi a G rho/(g * (2n + 1))
!        and they are, but shouldn't gravity be scaled by the same amount? There seems to be a factor of Ae missing.
!   APP130226: Missing factor of Ae is used when calculating the gradient function to balance the powers of (Ae/r)
!              So the original form is correct. Have introduced scale_factor term to reduce computational cost of this loop
!
    scale_factor = 4.d0*pi*Ae*Ae*rhow/Me
    do ideg = 0, nmax
      do iord = 0, ideg
        hlove(ideg) = hlove(ideg)*scale_factor*Ae/(2.d0*dble(ideg)+1.d0)
        dCotide(ideg, iord)=(1+klove(ideg))*scale_factor/((1.d0-hlove(ideg))*(2.d0*ideg+1))*dCotide(ideg,iord)
        dSotide(ideg, iord)=(1+klove(ideg))*scale_factor/((1.d0-hlove(ideg))*(2.d0*ideg+1))*dSotide(ideg,iord)
      enddo
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    return
    end

