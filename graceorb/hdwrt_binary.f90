    subroutine hdwrt_binary(tin,date,sectag, tinend, dateend, sectagend, tim, step, mjd, efic, incoor, time_array, frame, &
                     GPSant_mag, GPSant_cos,gpsantoff)

! writes header file for GRACE orbit integrator output
!
! MODS
! PT130902: write out as an unformatted, binary file instead of ascii
! PT170609: add inclusion of file gtorb_mod (contains record length of GTORB file) and write recl on first line of GTORB file

    use sat_mod
    use inmod_mod
    use bsscl_mod
    use accel_mod
    use mascon_mod
    use gtorb_mod    ! PT170609: this declares the record length information for the GTORB file

    implicit none

    integer*4 :: i
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
    real(kind=8), dimension(10) :: efic      ! earth-fixed ICs (pos, vel, twice-per-rev along-track acc
    character*11 :: coorspace, partspace
    character*40 :: inputfile                ! PT170127: increased from C*20
    character*10 :: fileform
    character*10 :: agency, uname
    character*5 frame
    real(kind=8) :: GPSant_mag
    real(kind=8), dimension(3) :: GPSant_cos, gpsantoff

    softversion = "alpha"      ! software version
    fileform = 'BIN'           ! output format of file
    headlines = 24             ! number of header lines
    datarecs = step+1          ! number of data records
    coorspace = 'terrestrial'  ! coordinates are in 'terrestrial' or 'inertial' space
    partspace = 'terrestrial'  ! IC partials are in 'terrestrial' or 'inertial' space
! PT140616: change to the actual name (rather than hardwired 'GRACE.input')
    inputfile = input_file  ! input files with models and components etc

    call getenv('USER',uname)  ! user name running programme
    call getenv('INSTITUTE',agency) ! agency where code is run


! PT170609: first line now contains the record length of the file
    nrec = 1
    write(15,rec=nrec)GTORB_recl


!   Print out header lines of output file:
    nrec=nrec+1
    message = " "
    write(message,130)agency
130 format (' PRODUCER AGENCY                 : ', a10)
    write (15, rec=nrec)message

    nrec=nrec+1
    message = " "
    write (message, 140) uname
140 format (' RUN BY                          : ', a10)
    write (15, rec=nrec)message

    nrec=nrec+1
    message = " "
    write (message, 145) time_array(1),time_array(2),time_array(3),time_array(5),time_array(6), time_array(7)
145 format (' PRODUCT CREATE TIME             : ',i4,'-',i2.2,'-',i2.2,' ',i2.2,':',i2.2,':',i2.2 )
    write (15, rec=nrec)message

    nrec=nrec+1
    message = " "
   write (message, 150) softversion
150 format (' SOFTWARE VERSION                : ', a5)
    write (15, rec=nrec)message

    nrec=nrec+1
    message = " "
    write (message, 160) fileform
160 format (' FILE FORMAT  (ASC or BIN)       : ', a3)
    write (15, rec=nrec)message

    nrec=nrec+1
    message = " "
    write (message, 170) sat
170 format (' SATELLITE NAME                  : ', a1)
    write (15, rec=nrec)message

    nrec=nrec+1
    message = " "
    write (message, 190) headlines
190 format (' NUMBER OF HEADER RECORDS        : ', i2)
    write (15, rec=nrec)message

    nrec=nrec+1
    message = " "
    write (message, 200) datarecs
200 format (' NUMBER OF DATA RECORDS          : ', i5)
    write (15, rec=nrec)message

    nrec=nrec+1
    message = " "
    write (message, 211)  tim
211 format (' DATA RECORD INTERVAL (SEC)      : ', f4.1)
    write (15, rec=nrec)message

    nrec=nrec+1
    message = " "
    write (message, 212) 2000, 1, 1, 12, 0, 0
212 format (' TIME EPOCH (GPS TIME)           : ',I4,'-',I2.2,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)
    write (15, rec=nrec)message

    nrec=nrec+1
    message = " "
    write (message, 210) nint(tin), date, nint(sectag)
210 format (' TIME FIRST OBS(SEC PAST EPOCH)  : ', i9,' (',I4,'-',I2.2,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2,')')
    write (15, rec=nrec)message

    nrec=nrec+1
    message = " "
    write (message, 220) nint(tinend), dateend, nint(sectagend)
220 format (' TIME LAST OBS(SEC PAST EPOCH)   : ', i9,' (',I4,'-',I2.2,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2,')')
    write (15, rec=nrec)message

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

    nrec=nrec+1
    message = " "
! PT131121: add the name of the mascon primary file here
    write (message, 230) statfieldmodname, oceantidemodname, atmtidemodname, dealiasmodname, mascon_primary_file
230 format (' MODELS USED                     : ', a5,':',a5,':',a5,':',a5,':',a150)
    write (15, rec=nrec)message

    nrec=nrec+1
    message = " "
    write (message, 231) coorspace 
231 format (' DATA RECORD REFERENCE FRAME     : ', a11)
    write (15, rec=nrec)message

    nrec=nrec+1
    message = " "
    write (message, 240)
240 format (' INITIAL CONDITIONS FORMAT       : gsec X0 Y0 Z0 XV0 YV0 ZV0 1prS 1prC 2prS 2prC')
    write (15, rec=nrec)message

    nrec=nrec+1
    message = " "
    write (message, 250) nint(tin), efic
250 format (' SATELLITE IC (TERRESTRIAL)      : ',i10,' ', 3f18.7, 3f18.7, 4f18.7)
    write (15, rec=nrec)message

    nrec=nrec+1
    message = " "
    write (message, 260) nint(tin), incoor
260 format (' SATELLITE IC (INERTIAL)         : ',i10,' ', 3f18.7, 3f18.7)
    write (15, rec=nrec)message

    nrec=nrec+1
    message = " "
    write (message, 241)
241 format (' ACCELEROMETER BIAS SCALE FORMAT : sclx scly sclz bsx  bsy  bsz ')
    write (15, rec=nrec)message

    nrec=nrec+1
    message = " "
     write (message, 262) scl(1), scl(2), scl(3), bs(1)*1.0d6, bs(2)*1.0d6, bs(3)*1.0d6
! PT131205: increased the number of significant figures on all (y-bias was being truncated)
262 format (' ACCELEROMETER BIAS SCALE        : ',12e16.8)
    write (15, rec=nrec)message

    nrec=nrec+1
    message = " "
    write (message, 263)
263 format (' GPS ANTENNA VECTOR FORMAT       : mag Xdcos Ydcos Zdcos GanX GanY GanZ ')
    write (15, rec=nrec)message

    nrec=nrec+1
    message = " "
    write (message, 264) GPSant_mag, GPSant_cos, gpsantoff
264 format (' GPS ANTENNA VECTOR LC           : ', 1f21.13, 3f21.13, 3f10.6)
    write (15, rec=nrec)message

    nrec=nrec+1
    message = " "
    write (message, 270) total_prim, total_ocean_prim
270 format (' DATA RECORD FORMAT              : gsec X Y Z XV YV ZV 4Q &
            6dX0 6dY0 6dZ0 6dXV0 6dYV0 6dZV0 6dS1 6dS2 6dS3 6dB1 6dB2 6dB3 61prS 61prC 62prS 62prC 6dMC',2i7)
    write (15, rec=nrec)message

! PT140613: write out the total number of tidal constituents to be estimated, summed over all the primary ocean mascons
! PT170611: also include the total number of ocean primarys, since this is what is written out!
    nrec=nrec+1
    message = " "
    write(message,271)total_prim,num_mcon_tide,combined_mascon_file,msc_hdr_code,ocean_mascon_file,msc_ocn_hdr_code
271 format(' # MASCONS , TIDAL CONSTITUENTS  : ',i5,i7,2(3x,a40,1x,a8))
    write(15,rec=nrec)message

    nrec=nrec+1
    message = " "
    write (message, 272) inputfile
272 format (' INPUT FILE NAME                 : ',a40)
    write (15, rec=nrec)message

    nrec=nrec+1
    message = " "
    write (message, 280) 
280 format (' END OF HEADER                   : ')
    write (15, rec=nrec)message

    write(message,'(a,i5,a)')'  Written ',nrec,' records to binary header'
    call status_update('STATUS','GRACEORB','hdwrt_binary',' ',message,0)

! PT140613: add the bit-mapped value to indicate which tidal constituent amplitudes are being estimated at each mascon
!           (0: none, 1: M2, 2: O1, 4: S2, 8: K1, 16: K2)
  if (gt_oceantide(1).eq.'G') then  
    nrec = nrec+1
    write (15,rec=nrec)(mcon_ocean_prim(i,2),i=1,total_ocean_prim)
    write (*,*)(mcon_ocean_prim(i,2),i=1,total_ocean_prim)
    write(message,'(a,i5,a,i7,a)')'  Written bit-mapped tidal information for ',total_ocean_prim_ampl,' mascons.' &
                             ,num_mcon_tide,' (x2) tidal amplitudes to estimate'
  else
    nrec=nrec+1
    write(message,'(a)')"No bit-mapped tidal information written out"
  endif  
  call status_update('STATUS','GRACEORB','hdwrt_binary',' ',message,0)

    return 
    end













