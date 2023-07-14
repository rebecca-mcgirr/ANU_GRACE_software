module f90_kind
  ! PORTABLE KIND DEFINITION MODULES FOR FORTRAN 90
  ! Author : Laurent HILL, SIMULOG Sophia Antipolis
  ! Laurent.Hill@sophia.inria.fr

  ! WARNING : AVAILABILITY AND CHARACTERISTICS OF FLOATING POINT
  ! FOMATS ARE INHERENTLY MACHINE DEPENDENT.
  ! SOME IMPLEMENTATION MAY NOT PROVIDE SOME OF THE SPECIFIED
  ! REPRESENTATION, OR SOME REPRESENTATIONS MAY NOT HAVE THE EXPECTED 
  ! CHARACTERISTICS. (i.e. integer(INT8) is allowed to be 32 bit long !)

  ! REAL kinds : a conservative number of digits is used
  !              to get the right size across various machine float formats
  integer, parameter :: REAL32  = selected_real_kind(6)
  integer, parameter :: REAL64  = selected_real_kind(12)
  integer, parameter :: REAL128 = selected_real_kind(24)

  ! INTEGER kinds
  integer, parameter :: INT8    = selected_int_kind(2)
  integer, parameter :: INT16   = selected_int_kind(4)
  integer, parameter :: INT32   = selected_int_kind(18)
  !  integer, parameter :: INT64   = selected_int_kind(18)

  ! LOGICAL kinds
  !
  !  integer, parameter :: BOOL8  = 1   ! WARNING : works only by accident
  ! No portable way to specify 2 byte logicals :
  ! integer, parameter :: BOOL16 = 2 
  !  integer, parameter :: BOOL32 = Kind(.true.) ! WARNING : works only by accident  

  ! Aliases, to be compatible with names in NAG's module 
  ! or European meteorological organizations 

  ! REAL kinds 
  integer, parameter :: SINGLE = REAL64
  integer, parameter :: DOUBLE = REAL128
  ! ANO2.11
  integer, parameter :: S = SINGLE
  integer, parameter :: D = DOUBLE
  !  integer, parameter :: QUAD   = REAL128
  ! LOGICAL kinds
  !  integer, parameter :: BYTE    = BOOL8
  ! integer, parameter :: TWOBYTE = BOOL16  ! No portable way to specify 2 byte 
  !  integer, parameter :: WORD    = BOOL32

  ! Aliases for C type names
  !  integer, parameter :: FLOAT = SINGLE
  !  integer, parameter :: LONGLONG   = INT64
  !  integer, parameter :: LONG   = INT32
  !  integer, parameter :: INT    = INT
  !  integer, parameter :: SHORT  = INT16 

end module f90_kind

program test
! 
use f90_kind
! 
!.. Implicit Declarations .. 
implicit none


integer :: numsat ! 1 = GRACEA, 2 = GRACEB
real(SINGLE) :: xsec_grace,date1950,xtime,xsec_tmp
real(SINGLE), dimension(3) :: biais_mod
integer*4 :: iepoch,nsecs
character*100 :: arg

call getarg(1,arg)
if(arg(1:1) == " ")then
  print*,'Runstring: compute_GRACE_a_priori_accelerometer_biases_from_temperature 1/2 gracesec n_secs'
  stop
endif
read(arg,*)numsat
call getarg(2,arg)
read(arg,*)xtime
call getarg(3,arg)
read(arg,*)nsecs


if (xtime < 70000000.) then
  date1950=xtime
  xsec_grace=0.
else
  date1950=0.
  xsec_grace=xtime
end if

! PT181022: loop over 24 hours from the start of the GRACESEC entered
  do iepoch = 1, nsecs
    xsec_tmp = xsec_grace+iepoch-1
    date1950 = 0.
    call subr_compute_GRACE_a_priori_acc_biases_from_temp(numsat,xsec_tmp,date1950,biais_mod)
    write (*,'(f17.6,3e20.12,i2,"   Satellite/Sec_GRACE/Bias_XYZ")') xsec_tmp,biais_mod,numsat
  enddo


end

subroutine subr_compute_GRACE_a_priori_acc_biases_from_temp(numsat,xsec_grace,date1950,biais_mod)
!
! This routine computes the a priori acceleromater biases (biais_mod) of satellite numsat (1 = GRACEA, 2 = GRACEB) from an empirical
! law based on the date and on the temperature of the gauge closest to the center of the accelerometer.
!
! This routine accepts as input temporal parameters either the GRACE time scale in TAI seconds (xsec_grace),
! or the date in real julian days after 1/1/1950 (date1950), BUT NOT BOTH AT THE SAME TIME (one of these dates has to be 0).
! On output, both time scales are computed and returned.
!
! This routine needs the files:
! Temperature_internal_core.GRACEA.100s.dat
! Temperature_internal_core.GRACEB.100s.dat
! which are based on the internal core temperature of the accelerometers, found in the Level-1B housekeeping data files AHK1B.
! The raw temperature measurement is filtered by a low-pass filter. A time constant of approximately 27700 s for the transmission 
! of the temperature signal between the gauge closest to the center of the accelerometer and heart of the accelerometer has been
! found empirically. It is taken into account in the subroutine. 
! Due to a data gap of the temperature measurements on GRACE-A in the AHK1B files between 2013.74947778 and 2014.48790558, 
! the temperatures for this time period are replaced by a synthetic temperature computed from the accelerometer biases.
!
! La routine accepte comme paramètres temporels soit la date en jj1950 décimal (échelle TAI), 
! soit l'échelle propre GRACE en secondes (temps GPS), mais pas les deux en même temps.
! En sortie, les deux échelles de temps sont renseignées.
!
! Loi adoptée :
! y = a0+a1*ln(t-19072)+a2*(t-19072)+a3*T°(t-0.32)
! avec t en jj1950 et T° en degrés celsius. On utilise la température de t-0.32 j pour tenir compte 
! de la constante de temps de transmission de la température jusqu'au coeur de l'accéléromètre de
! l'ordre de 0.32 j, soit environ 27700 secondes. 
!
  use f90_kind
  ! 
  !.. Implicit Declarations .. 
  implicit none
  ! 
  !.. Formal Arguments .. 
  integer, intent (in) :: numsat ! 1 = GRACEA, 2 = GRACEB
  real(SINGLE), intent (inout) :: xsec_grace,date1950
  real(SINGLE), dimension(3), intent (inout) :: biais_mod
  !
  ! Parameters
  integer, parameter :: it100_deb = 0700000, it100_fin = 5700000
  !
  !.. Local Scalars ..
  integer, save :: iprem = 0
  integer :: indice_av, indice_ap, i
  real(SINGLE) :: temp,temp_av,temp_ap,f_av,f_ap,xsec_temp_decalee
  !  
  !.. Local Arrays ..
  real(SINGLE), save, dimension(it100_deb:it100_fin,2) :: temperature ! Une valeur toute les 100 secondes tsec_grace
  real(SINGLE), save, dimension(4,3,2) :: coefficients_loi_jml
  ! 
  ! ... Executable Statements ...
  !
  if (iprem == 0) then
    ! Premier appel : chargement des tableaux de température
    iprem = 1
    temperature = 0._S
    print*,'Reading file Temperature_internal_core.GRACEA.100s.dat'
    open (54321,file='Temperature_internal_core.GRACEA.100s.dat',status='old')
    do
      read (54321,*,end=100) i,temperature(i/100,1)
    end do
    100 continue
    close (54321)
    print*,'Reading file Temperature_internal_core.GRACEB.100s.dat'
    open (54321,file='Temperature_internal_core.GRACEB.100s.dat',status='old')
    do
      read (54321,*,end=200) i,temperature(i/100,2)
    end do
    200 continue
    close (54321)
    ! Premier appel : chargement du tableau des coefficients de la fonction JML
    coefficients_loi_jml(:,1,1)=(/  1.099490e-06_S,  4.222290e-08_S,  2.427620e-12_S, -7.209890e-09_S/)
    coefficients_loi_jml(:,2,1)=(/ -2.628957e-05_S, -1.080220e-06_S,  2.690460e-10_S,  1.652360e-07_S/)
    coefficients_loi_jml(:,3,1)=(/  6.962690e-07_S,  1.184700e-08_S, -8.062500e-12_S, -2.375540e-09_S/)
    coefficients_loi_jml(:,1,2)=(/  5.045210e-07_S,  2.906670e-08_S, -4.413980e-12_S, -3.595680e-09_S/)
    coefficients_loi_jml(:,2,2)=(/ -2.614180e-06_S, -1.646690e-06_S,  2.712740e-10_S,  1.265760e-07_S/)
    coefficients_loi_jml(:,3,2)=(/  1.421470e-06_S, -6.335580e-08_S,  7.552690e-12_S, -2.617490e-09_S/)
  end if

  ! Contrôle sur les dates
  if (xsec_grace .ne. 0._S .and. date1950 .ne. 0._S) then
    write(6,'(/" *** stop : Entering the subroutine with both time scales is forbidden"/)')
    stop
  end if

  ! calcul de la deuxième date
  if (xsec_grace .ne. 0._S) then
    date1950=18262._S+(xsec_grace+43219._S)/86400._S
  else
    xsec_grace=(date1950-18262._S)*86400._S-43219._S
  end if

  ! calcul de la température à l'instant voulu
  xsec_temp_decalee=xsec_grace-27700 ! constante de temps de transmission de la T° : environ 27700 secondes
  indice_av=int(xsec_temp_decalee/100._S)
  indice_ap=indice_av+1
  ! Test sur la date
  if (indice_av .lt. it100_deb .or. indice_ap .gt. it100_fin) then
    write(6,"(/'*** stop : The date is outside the interval '/)")
    stop
  end if
  temp_av=temperature(indice_av,numsat)
  temp_ap=temperature(indice_ap,numsat)
  f_av=(real(indice_ap*100,SINGLE)-xsec_temp_decalee)/100._S
  f_ap=(xsec_temp_decalee-real(indice_av*100,SINGLE))/100._S
  temp=temp_av*f_av+temp_ap*f_ap

  ! calcul des biais accélérométriques a priori
  biais_mod(:)= &
    coefficients_loi_jml(1,:,numsat) + &
    coefficients_loi_jml(2,:,numsat)*log(date1950-19072._S) + &
    coefficients_loi_jml(3,:,numsat)*(date1950-19072._S) + &
    coefficients_loi_jml(4,:,numsat)*temp

  ! A partir de l'extinction de l'accéléromètre GRAB, il faut prendre les biais de GRAA pour traiter GRAB
  ! en inversant le signe suivant X et Y, mais pas Z (radial)
  ! La date de changement est prise au milieu de l'interruption de mesures entre 20160902 (=24351) et 20161114 (=24424)
  if (numsat == 2 .and. date1950 > 24387._S) then
    temp_av=temperature(indice_av,1)
    temp_ap=temperature(indice_ap,1)
    f_av=(real(indice_ap*100,SINGLE)-xsec_temp_decalee)/100._S
    f_ap=(xsec_temp_decalee-real(indice_av*100,SINGLE))/100._S
    temp=temp_av*f_av+temp_ap*f_ap
    biais_mod(:)= &
      coefficients_loi_jml(1,:,1) + &
      coefficients_loi_jml(2,:,1)*log(date1950-19072._S) + &
      coefficients_loi_jml(3,:,1)*(date1950-19072._S) + &
      coefficients_loi_jml(4,:,1)*temp
    biais_mod(1:2) = -biais_mod(1:2)
  end if
    
end subroutine subr_compute_GRACE_a_priori_acc_biases_from_temp


