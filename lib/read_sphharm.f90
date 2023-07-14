  subroutine read_sphharm(lu_sphharm, sphharm_file, max_deg, max_exp, low_deg, up_deg, cos_coeff, sin_coeff)

! Subroutine name: read_spharm
! Author:          Anthony Purcell
! Date:            September 24 2013 
! Purpose:         To read in  series of spherical harmonic coefficients

  implicit none

! ----------------------------- Internal Variable declaration ---------------------------------------------------------------------
  integer*4 check, count, degree, ioerr, izero, m, m1, n, n0, n1, total

  dimension check(0:max_deg, 0:max_deg)

  real*8 cos_tmp, sin_tmp
! ---------------------------------------------------------------------------------------------------------------------------------

! ----------------------------- Input Parameter declaration -----------------------------------------------------------------------
  integer*4 low_deg, lu_sphharm, max_deg, max_exp, up_deg

  real*8 cos_coeff, sin_coeff 

  dimension cos_coeff(0:max_deg,0:max_deg), sin_coeff(0:max_deg,0:max_deg)

  character sphharm_file*200, message*256
! ----------------------------------------------------------------------------------------------------------------------------------


! ----------------------------- Open & verify file of spherical harmonic coefficients ----------------------------------------------
  open(unit = lu_sphharm, file = sphharm_file, status = 'old', iostat =  ioerr)
  if (ioerr /= 0) then
    write(message,'(a,a,a)')"Error opening input file ",sphharm_file,". Does it exist?"
    call status_update('FATAL','LIB ','read_sphharm',' ',message,ioerr)
  else
    write(message,'(a,a)')"Reading file of spherical harmonic coefficients: ", sphharm_file
    call status_update('STATUS','LIB ','read_sphharm',' ',message,0)
    write(message,'(a,i4)')"Maximum Degree of expansion: ", max_exp
    call status_update('STATUS','LIB ','read_sphharm',' ',message,0)
    write(message,'(a,i4)')"Dimension of coefficient array: ", max_deg
    call status_update('STATUS','LIB ','read_sphharm',' ',message,0)
  endif
! ----------------------------------------------------------------------------------------------------------------------------------

! ----------------------------- Explicitly initialise variables --------------------------------------------------------------------
  check   = 0
  count   = 0
  ioerr   = 0
  izero   = 0
  low_deg = 10000000
  up_deg  = 0
  degree  = min(max_deg, max_exp)
! ----------------------------------------------------------------------------------------------------------------------------------

! ----------------------------- Read file - program does not assume any ordering is applied ----------------------------------------
  do while ( ioerr .eq. 0 )
    read(lu_sphharm, *, end = 110, iostat = ioerr) n1, m1, cos_tmp, sin_tmp

    if ( ( min(n1, m1) .lt. 0 ) .or. ( max(n1, m1) .gt. max_deg ) ) then
      write(message,'(a,a,a,i4,i4)') "Error reading ", sphharm_file, " invalid value for degree or order ", n1, m1
      call status_update('FATAL','LIB ','read_sphharm',' ',message,0)
    endif
    if ( m1 .gt. n1 ) then
      write(message,'(a,a,a,i3,a,i3,a)') "Error reading ", sphharm_file, " order (", m1,") larger than degree (", n1,")."
      call status_update('FATAL','LIB ','read_sphharm',' ',message,0)
    endif

    if ( n1 .le. degree ) then
      if ( check(n1, m1) .eq. 0 ) then
        up_deg = max(up_deg, n1)
        low_deg = min(low_deg, n1)
        cos_coeff(n1, m1) = cos_tmp
        sin_coeff(n1, m1) = sin_tmp
        check(n1, m1) = 1
        count = count + 1
      else
        write(message,'(a,a,a,i4,i5)') "Error reading ", sphharm_file, " multiple entries for degree & order ", n1, m1
        call status_update('FATAL','LIB ','read_sphharm',' ',message,0)
      endif
    endif
  enddo               ! End of while loop
110 close(lu_sphharm)
! ----------------------------------------------------------------------------------------------------------------------------------

! ----------------------------- Report maximum & minimum degrees read --------------------------------------------------------------
  write(message,'(a,a,a,i6)')"Maximum degree of spherical harmonic expansion: ", sphharm_file,": ",up_deg
  call status_update('STATUS','LIB ','read_sphharm',' ',message,0)
  write(message,'(a,a,a,i6)')"Minimum degree of spherical harmonic expansion: ", sphharm_file,": ",low_deg
  call status_update('STATUS','LIB ','read_sphharm',' ',message,0)
! ----------------------------------------------------------------------------------------------------------------------------------

! ----------------------------- Verify that we have read the right number of entries -----------------------------------------------
  total = (up_deg - low_deg + 1) * (up_deg + low_deg + 2)/2
  if ( count .ne. total ) then
    write(message,'(a,i6,a,i6,a)')"Total number of coefficients read: ",count,". Expected number: ",total,"."
    call status_update('FATAL','LIB ','read_sphharm',' ',message,0)
  endif
! ----------------------------------------------------------------------------------------------------------------------------------

  return
  end
