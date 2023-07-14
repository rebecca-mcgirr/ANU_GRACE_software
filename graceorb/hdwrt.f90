    subroutine hdwrt(tin,date,sectag, tinend, dateend, sectagend, tim, step, mjd, efic, incoor, time_array, frame, &
                     GPSant_mag, GPSant_cos,gpsantoff)

! writes header file for GRACE orbit integrator output

    use sat_mod
    use inmod_mod
    use bsscl_mod
    use accel_mod

    implicit none

    integer :: headlines
    integer*4 :: datarecs
    integer*4 mjd
    integer :: time_array(8)
    integer*4 :: step
!    integer*4 :: num_mcon
    integer*4 :: date(5), dateend(5)
    real(kind=8) :: tim, tin, tinend
    real*8 :: sectag, sectagend
!    real*8 :: bs(3)
    real(kind=8), dimension(6) :: incoor
    real(kind=8), dimension(6) :: efic
    character*11 :: coorspace, partspace
    character*20 :: inputfile
    character*10 :: fileform
    character*10 :: agency, uname
    character*5 frame
    real(kind=8) :: GPSant_mag
    real(kind=8), dimension(3) :: GPSant_cos, gpsantoff

    softversion = "alpha"      ! software version
    fileform = 'ASC'           ! output format of file
    headlines = 23             ! number of header lines
    datarecs = step+1          ! number of data records
    coorspace = 'terrestrial'  ! coordinates are in 'terrestrial' or 'inertial' space
    partspace = 'terrestrial'  ! IC partials are in 'terrestrial' or 'inertial' space
    inputfile = 'GRACE.input'  ! input files with models and components etc

    call getenv('USER',uname)  ! user name running programme
    call getenv('INSTITUTE',agency) ! agency where code is run


!   Print out header lines of output file:

    write (15, 130) agency
130 format (' PRODUCER AGENCY                 : ', a10)

    write (15, 140) uname
140 format (' RUN BY                          : ', a10)

    write (15, 145) time_array(1),time_array(2),time_array(3),time_array(5),time_array(6), time_array(7)
145 format (' PRODUCT CREATE TIME             : ',i4,'-',i2.2,'-',i2.2,' ',i2.2,':',i2.2,':',i2.2 )

    write (15, 150) softversion
150 format (' SOFTWARE VERSION                : ', a5)

    write (15, 160) fileform
160 format (' FILE FORMAT  (ASC or BIN)       : ', a3)

    write (15, 170) sat
170 format (' SATELLITE NAME                  : ', a1)

    write (15, 190) headlines
190 format (' NUMBER OF HEADER RECORDS        : ', i2)

    write (15, 200) datarecs
200 format (' NUMBER OF DATA RECORDS          : ', i5)

    write (15, 211)  tim
211 format (' DATA RECORD INTERVAL (SEC)      : ', f4.1)

    write (15, 212) 2000, 1, 1, 12, 0, 0
212 format (' TIME EPOCH (GPS TIME)           : ',I4,'-',I2.2,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)

    write (15, 210) nint(tin), date, nint(sectag)
210 format (' TIME FIRST OBS(SEC PAST EPOCH)  : ', i9,' (',I4,'-',I2.2,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2,')')

    write (15, 220) nint(tinend), dateend, nint(sectagend)
220 format (' TIME LAST OBS(SEC PAST EPOCH)   : ', i9,' (',I4,'-',I2.2,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2,')')

    if (gt_statfield(1).eq."Y") then
      statfieldmodname = gt_statfieldmod(1)
    else
      statfieldmodname = " "
    endif

    if (gt_oceantide(1).eq."Y") then
      oceantidemodname = gt_oceantidemod(1)
    else
      oceantidemodname = " "
    endif

    if (gt_atmtide(1).eq."Y") then
      atmtidemodname = gt_atmtidemod(1)
    else
      atmtidemodname = " "
    endif

    if (gt_dealias(1).eq."Y") then
      dealiasmodname = gt_dealiasmod(1)
    else
      dealiasmodname = " "
    endif

    write (15, 230) statfieldmodname, oceantidemodname, atmtidemodname, dealiasmodname
230 format (' MODELS USED                     : ', a5,':',a5,':',a5,':',a5)

    write (15, 231) coorspace 
231 format (' DATA RECORD REFERENCE FRAME     : ', a11)

    write (15, 240)
240 format (' INITIAL CONDITIONS FORMAT       : gsec X0 Y0 Z0 XV0 YV0 ZV0')

    write (15, 250) nint(tin), efic
250 format (' SATELLITE IC (TERRESTRIAL)      : ',i10,' ', 3f18.7, 3f18.7)

    write (15, 260) nint(tin), incoor
260 format (' SATELLITE IC (INERTIAL)         : ',i10,' ', 3f18.7, 3f18.7)

    write (15, 241)
241 format (' ACCELEROMETER BIAS SCALE FORMAT : sclx scly sclz bsx  bsy  bsz ')

!    if (mjd.le.52705) then
!      write (15, 262) scl(1), scl(2), scl(3), c0x(1), c0y(1), c0z(1), c1x(1), c1y(1), c1z(1), c2x(1), c2y(1), c2z(1) 
!    else
!      write (15, 262) scl(1), scl(2), scl(3), c0x(2), c0y(2), c0z(2), c1x(2), c1y(2), c1z(2), c2x(2), c2y(2), c2z(2)
!    endif
     write (15, 262) scl(1), scl(2), scl(3), bs(1)*1.0d6, bs(2)*1.0d6, bs(3)*1.0d6
! PT131016: increased to 15.7 to write out more significant figures on the biases
262 format (' ACCELEROMETER BIAS SCALE        : ',12d17.9)

    write (15, 263)
263 format (' GPS ANTENNA VECTOR FORMAT       : mag Xdcos Ydcos Zdcos GanX GanY GanZ ')

    write (15, 264) GPSant_mag, GPSant_cos, gpsantoff
264 format (' GPS ANTENNA VECTOR LC           : ', 1f21.13, 3f21.13, 3f10.6)


!    write (15, 265) coorspace
!265 format (' COORDINATE OUTPUT SPACE         : ',a11)

!    write (15, 266) partspace
!266 format (' PARTIAL OUTPUT SPACE            : ',a11)

    write (15, 270) num_mcon
270 format (' DATA RECORD FORMAT              : gsec X Y Z XV YV ZV 4Q &
            6dX0 6dY0 6dZ0 6dXV0 6dYV0 6dZV0 6dS1 6dS2 6dS3 6dB1 6dB2 6dB3 6dMC',i5)

    write (15, 271) inputfile
271 format (' INPUT FILE NAME                 : ',a20)

    write (15, 280) 
280 format (' END OF HEADER                 : ')

      return 
      end













