
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  read_atm_tide_constituent   :  read the RL06 AOD1B tide constituent files
!  calc_atm_tide_constituent   :  calculate the cosine/sine components of a tidal constituent at a particular epoch
!  output_sum_tides            :  output the sum of the tidal constituents in a file for reading by graceorb etc
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_atm_tide_constituent(maxdeg,tide_file,tides)

! subroutine to read in the spherical harmonic models of cosine/sine
! components of a tidal constituent
!
! P. Tregoning
! 6 February 2020

  implicit none

! passed arguments
  integer*4    :: maxdeg                       ! maximum degree of spherical harmonic models
  character(*) :: tide_file                    ! name of tide file to be opened and read
  real(kind=8) :: tides(0:maxdeg,0:maxdeg,2,2) ! cosine and sine amplitude for each spherical harmonic coefficient

! local variables
  integer*4,parameter :: lu_file=10            ! unit number to open tide constituent file
  integer*4           :: nobs                  ! total number of data records in file (cos+sine)
  integer*4           :: iobs,ideg,iord,ioerr,icomp
  character*100       :: line
  real(kind=8)        :: tmpC,tmpS

! open the file
  open(lu_file,file=trim(tide_file),status='old',iostat=ioerr)
  if(ioerr /= 0)then
    call status_update('FATAL','UTIL','read_atm_tide_constituent',tide_file,"Error opening file",0)
  else
    call status_update('STATUS','UTIL','read_atm_tide_constituent',tide_file,"Reading spherical harmonic constituents: ",0)
  endif

! skip the header information  
  line = " "
  do while (line(1:13) /= "END OF HEADER")
    read(lu_file,'(a)')line
  enddo

! read the next line to find out which "DATA SET". 01=cos, 02 = sin
  read(lu_file,'(9x,i2,1x,i7)')icomp,nobs

! read all the components of this data set
  do iobs=1,nobs
    read(lu_file,*)ideg,iord,tmpC,tmpS
    tides(ideg,iord,1,icomp) = tmpC
    tides(ideg,iord,2,icomp) = tmpS
  enddo

! now read the next data set
  read(lu_file,'(9x,i2,1x,i7)')icomp,nobs

! read all the components of this data set

  do iobs=1,nobs
    read(lu_file,*)ideg,iord,tmpC,tmpS
    tides(ideg,iord,1,icomp) = tmpC
    tides(ideg,iord,2,icomp) = tmpS
  enddo

  close(lu_file)
  return
  end subroutine read_atm_tide_constituent
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calc_atm_tide_constituent(constituent,freq,t,maxdeg,tides,sum_tides)

! subroutine to calculate the cosine/sine components of a tide at a particular epoch, t (seconds since start of day).
! Here we need to calculate the time dependence of the in-phase and out-of-phase components to produce a single C/S coefficient 
! at the required epoch.
!
! Time arguments are taken from Table 5.1 of the AOD1B product description (Dobslaw et al., 2017)
!
! P. Tregoning
! 6 February 2020

  implicit none

! passed arguments
  character*2      :: constituent
  real(kind=8)     :: freq                             ! tidal constituent frequency (in degrees/hour)
  real(kind=8)     :: t                                ! time (seconds) since start of day
  integer*4        :: maxdeg                           ! maximum degree of spherical harmonic model
  real(kind=8)     :: tides(0:maxdeg,0:maxdeg,2,2)     ! cos/sine amplides for the in-phase and out-of-phase components of tide constituent
  real(kind=8)     :: sum_tides(0:maxdeg,0:maxdeg,2)   ! sum of cos/sine amplides for the in-phase and out-of-phase components of tide constituent

! local variables
  real(kind=8)     :: time_arg
  integer*4        :: ideg,iord
  real(kind=8)     :: pi

! time variables (incomplete)
  real(kind=8) :: h,s,p,pprime,JD

  pi = 4.d0*datan(1.d0)
  
  ! calculate the mean longitudes, perigee etc for Sun and Moon
  pprime =  1.d0 ! longitude of the Sun's mean perigee 
  h = 280.460d0 + 0.9856474d0*(JD - 2451545.d0)  ! mean longitude of the Sun 

  p =  1.d0 ! longitude of the Moon's mean perigee 

! mean longitude of the Moon can be found from IERS2010 standards Eqs 5.43 (P 67)
! s = F3 - F5
!   = (L - Omega) - Omega

! here t = time in seconds since 1200UT 1/1/2000
  s =   93.27209062d0 + (1739527262.8478d0 * t - 12.7512d0 *t**2 - 0.001037d0 * t**3 + 0.00000417d0 *t**4)/3600.d0  &
      - (125.04455501d0 - (6962890.543100d0*t + 7.472200d0*t**2 + 0.00770200d0*t**3 - 0.0000593900d0*t**4)/3600.d0 ) 


 ! mean longitude of the Moon in degrees


! PT200207: write the code for one constituent at a time ...
!           solar tides, depend only on time of day
! HM200223 changed time_arg calculation from freq*24/360 to freq/3600
  
  if(constituent == "blah")then    ! will never happen - just allows me to comment out each constituent to test them all ....
    return

!  else if (constituent == "P1")then
!    time_arg = (t-h-21600.d0)               *freq/3600.d0 * pi/180.d0
  else if(constituent == "S1")then
    time_arg = (t + 43200.d0)               *freq/3600.d0 * pi/180.d0    

!  else if (constituent == "K1")then
!    time_arg = (t+h+21600.d0)               *freq/3600.d0 * pi/180.d0

!  else if (constituent == "N2")then
!    time_arg = (2.d0*t-3.d0*s+2.d0*h+p)     *freq/3600.d0 * pi/180.d0
!  else if (constituent == "M2")then   
!    time_arg = (2.d0*t-2.d0*s+2.d0*h)       *freq/3600.d0 * pi/180.d0
!  else if (constituent == "L2")then
!    time_arg = (2.d0*t-s+2.d0*h-p+43200.d0) *freq/3600.d0 * pi/180.d0


!  else if (constituent == "T2")then
!    time_arg = (2.d0*t - h + pprime)        *freq/3600.d0 * pi/180.d0
  else if (constituent == "S2")then
    time_arg = 2.d0*t                       *freq/3600.d0 * pi/180.d0
!  else if (constituent == "R2")then
!    time_arg = (2.d0*t + h - pprime + 43200.d0) *freq/3600.d0 * pi/180.d0

!  else if (constituent == "T3")then
!    time_arg = (3.d0*t-h)                   *freq/3600.d0 * pi/180.d0
  else if (constituent == "S3")then
    time_arg = 3.d0*t                       *freq/3600.d0 * pi/180.d0
!  else if (constituent == "R3")then
!    time_arg = (3.d0*t+h)                   *freq/3600.d0 * pi/180.d0

  else
!    if(t.eq.0.d0) print*,'constituent ',constituent,' is not coded'
    return
  endif

! loop through each spherical harmonic coefficient and calculate the in-phase and out-of-phase components for this epoch

  do ideg = 0,maxdeg
    do iord = 0,ideg 
      sum_tides(ideg,iord,1) = sum_tides(ideg,iord,1) + tides(ideg,iord,1,1)*dcos(time_arg) + tides(ideg,iord,1,2)*dsin(time_arg)
      sum_tides(ideg,iord,2) = sum_tides(ideg,iord,2) + tides(ideg,iord,2,1)*dcos(time_arg) + tides(ideg,iord,2,2)*dsin(time_arg)
      if(iord.eq.0.and.ideg.eq.0) write(99,*) t,time_arg,constituent,sum_tides(ideg,iord,1), tides(ideg,iord,1,1)*dcos(time_arg), tides(ideg,iord,1,2)*dsin(time_arg) 


! DEBUG
!if(ideg < 3 .and. t < 600.)then
!print*,ideg,iord,sum_tides(ideg,iord,:),t
!endif

    enddo
  enddo
!  write(99,*) sum_tides(0,0,1),sum_tides(1,1,1)

  return
  end subroutine calc_atm_tide_constituent
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine output_sum_tides(n_epochs,maxdeg,sum_tides)

! subroutine to output the spherical harmonic models for the cosine/sine components at every requested epoch
!
! P. Tregoning
! 6 February 2020

  implicit none

! passed variables
  integer*4         :: n_epochs                                ! number of epochs to output
  integer*4         :: maxdeg                                  ! maximum degree of spherical harmonic model
  real(kind=8)      :: sum_tides(0:n_epochs-1,0:maxdeg,0:maxdeg,2) ! summed tidal constituents at every requested epoch

! local variables
  integer*4           :: ioerr,iepoch,ideg,iord
  integer*4,parameter :: lu_file=10            ! unit number to open tide constituent file
  character*200       :: message
! open the output file
  open(lu_file,file="atm_AOD1B_tides.sph",status='unknown',iostat=ioerr)
  if(ioerr /= 0)then
    call status_update('FATAL','UTIL','read_atm_tide_constituent',"atm_AOD1B_tides.sph","Error opening file",0)
  else
    call status_update('STATUS','UTIL','read_atm_tide_constituent',"atm_AOD1B_tides.sph","Opened output file ",0)
  endif

print*,'n_epochs = ',n_epochs

! first line is the degree of the spherical harmonic models
  write(lu_file,*)maxdeg

! now, loop through all the epochs and output the spherical harmonic field, ordered by degree (then order)
  do iepoch = 0,n_epochs-1
    write(message,'(a,i6)')'Output epoch: ',iepoch
    if(mod(iepoch,10)==0)call status_update('STATUS','UTIL','read_atm_tide_constituent','',message,0)
    do ideg=0,maxdeg
      do iord = 0, ideg
        write(lu_file,1000)ideg,iord,sum_tides(iepoch,ideg,iord,:),"  0.01  0.01"
1000    format(2i7,2e25.17,a)
      enddo
    enddo
  enddo

  close(lu_file)
  return

  end subroutine output_sum_tides
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

