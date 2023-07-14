!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Subroutines to interact with the SPOTL routines to extract ocean height
!
! spotl_point_setup:   sets up amplitudes and phases of constituents for a 
!                        particular model at a particular location
! spotl_point_calc :   calculates the ocean height at a location, given a
!                        set of amplitudes and phases
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine spotl_point_setup(in_lat,in_lon,tidemod,cart_codes,darw_codes,ntmcon,amp,phase,debug_flag)

! subroutine to interact with the spotl routines and, eventually, generate an estimate of ocean tide height at a point at a given time
! This is based on a program by Simon McClusky that did the same but for a global grid.
!
! P. Tregoning
! 10 September 2013

  implicit none

  double precision, intent(in) :: in_lat, in_lon            ! floating point lat/lon of requested point
  character, intent(in)        :: tidemod*10                ! name of requested tide model

  integer*4, intent(out)       :: cart_codes(6,20)          ! cartwright codes as read from the input spotl codes
  character, intent(out)       :: darw_codes(20)*4          ! darwin codes
  real*4, intent(out):: amp(20)                             ! constituent amplitudes
  real*4, intent(out):: phase(20)                           ! constituent phases
  logical, intent(out)         :: debug_flag                ! whether to print to screen the phases and amplitudes 

  character                    :: tidemod_ending*10         ! end of name of requested tide model (read from otide.models)
  integer*4                    :: ioerr                     ! i/o error flag
  integer                      :: ntmcon                    ! number of tidal constituents for requested model (read from otide.models)
  character line*200
  character message*200
  logical found
  double precision             :: rlat,rlon                 ! point location in radians
  double precision             :: ocean_ht                  ! computed ocean tidal height
  integer*4                    :: i,j,trimlen
  integer*4                    :: timearr(5)                ! yr, doy, hr, min, sec

! declaration of the original spotl variables (don't allow them to be declared implicitly)
      real*8 tlat,blat,wlong,elong,xpoly,ypoly
      real*4 rlato,rlono
      integer*4 latc,longc,i1,i2,np,npoly,polynm,icte
      

! ***@@ come back and delete the variables that are not required 
      integer*4 index,numhed,llst,lgrd,nextf
      real*8 ltg(256),rad2deg,pi,version
      logical*1 ispoly,use
      character*1 modo,fingrd
      character*2 tconst(20)
      character*4 dsym
      character*8 who
      character*16 outgridfile
      character*80 otmod,mdfile,dumm
      character*80 polyf,cmd
      character*50 mdnam
      character*40 stnam
      complex camp,toloc,phasor,oamp,cz
      dimension i1(2),i2(2)

!  common block holds polygon file information
      common/polinf/xpoly(500,10),ypoly(500,10),np,npoly(10),polyf,polynm(10),use(10),ispoly
      common /modpar/mdfile,dsym,icte(6),mdnam
      common /modlim/tlat,blat,wlong,elong,latc,longc
! SCM Need this common to move to initialize next file constituent file correctly ...
      common /newmod/nextf
      phasor(toloc) = cmplx(cabs(toloc),57.2958*atan2(aimag(toloc),real(toloc)))
! PT130910: I think I need this one too ...
      real*4 rho
      common/ldens/rho

      data ispoly/.false./
      data cz/(0.,0.)/,oamp/(0.,0.)/
      data stnam/'  Interpolated ocean-tide value         '/
      
! define pi and radians-to-degree converter     
  pi = 4.d0*atan(1.d0)
  rad2deg = 180.d0/pi
      
!  Get list of tide constituents in the ocean tide model
  open(10,file="otide.models",status='old',iostat=ioerr)
  if(ioerr /= 0)then
    call report_stat('FATAL','UTIL','spotl_setup','otide.models','Cannot open file. Check link.',0)
  endif

! read through the file to find the model that we want
  found = .false.
  do while (.not. found)
    read(10,'(a)',iostat=ioerr,end=1000)line
    if (line(1:6) == tidemod(1:6))then
      found = .true.

! get the ending of the file names
      read(line(11:17),'(a)')tidemod_ending

! get number of constituents
      read(line(18:20),'(i3)')ntmcon
      do i=1,ntmcon
        tconst(i) = line(22+(i-1)*3:22+(i-1)*3+1)
      enddo
    endif
  enddo

1000 continue

! check that the model was found
  if (.not. found) then
    write(message,'(a,a,a)')"Tide model (",tidemod,") not found in file"  
    call report_stat('FATAL','UTIL','spotl_setup','otide.models',message,0)
  else
    close(10)
    if(ntmcon == 12) then
      write(message,'(i4,a,12(1x,a2),a,a)')ntmcon,' constituents (',(tconst(i),i=1,ntmcon),') found for model ',tidemod
    else
      write(message,'(i4,a,a)')ntmcon,' constituents found for model ',tidemod
    endif
    if(debug_flag)call report_stat('STATUS','UTIL','spotl_setup',' ',message,0)
  endif

   

! don't know yet what this does .....
      fingrd = 'F'
        
! Loop over constituent files  
  do i = 1,ntmcon
    write(mdfile,'(a,a2,a1,a,a1,a)')'spotl/',tconst(i),'.',tidemod(1:trimlen(tidemod)),'.',tidemod_ending(1:trimlen(tidemod_ending))
    nextf = 1     
    rlato = in_lat
    rlono = in_lon
    call ocmodl(rlato,rlono,fingrd,camp)
    if(camp.ne.cz) camp=camp/rho
    camp = conjg(camp)
    camp = phasor(camp)
    amp(i) = real(camp)
    phase(i) = aimag(camp)
! save off the cartwright and darwin codes
    do j=1,6
      cart_codes(j,i) = icte(j)
    enddo
    darw_codes(i) = dsym

    write(message,'(a,a,f12.4,a,f12.4)')darw_codes(i),'Ampl:',amp(i),'   Phase:',phase(i)
    if(debug_flag)call report_stat('STATUS','UTIL','spotl_setup',' ',message,0)
  enddo

  return

  end subroutine spotl_point_setup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!************************!!!!!!!!!!!!!!!!!!!!!!!
  subroutine spotl_point_calc(timearr,n_constit,cart_codes,darw_codes,amp,phase,ocean_ht)

! subroutine to send the required information to relevant spotl routines so that the ocean height is then calculated
!
! based on code by Simon McClusky (who pulled it out of hartid.f which is by Duncan Agnew)
!
! P. Tregoning
! 10 September 2013

  implicit none

  integer*4, intent(in)         :: timearr(5)              ! epoch for tide computation (yr, doy, hr, min, sec)
  integer*4, intent(in)         :: n_constit               ! number of tidal constituents
  real*4, intent(in )           :: amp(n_constit)          ! array of tidal amplitudes
  real*4, intent(in )           :: phase(n_constit)        ! array of tidal phases
  integer*4, intent(in)         :: cart_codes(6,20)        ! cartwright codes
  character, intent(in)         :: darw_codes(20)*4        ! darwin codes

  double precision, intent(out) :: ocean_ht                ! computed ocean height for this epoch

! local/spotl variables
  integer ncon, nt
  parameter (ncon = 20)
  integer*4   itm(5)
  integer*4 i,j,k,upperc
  logical notfound
  integer*4 idt(6,ncon)

  integer*4 irli,samp,irnt,irhi,np
  real*4, allocatable :: a(:),x(:),hc(:)
  double precision, allocatable :: scr(:),f(:),p(:),wf(:)
  double precision dr,pi      
  data  nt/342/
! we need to pass the epoch through a common (through admint to be used in tdtrph .... how convoluted is that???)
  common/date/itm


  allocate(a(nt))
  allocate(f(nt))
  allocate(p(nt))
  allocate(x(600))
  allocate(hc(2*nt))
  allocate(scr(3*nt))
  allocate(wf(nt))

!  dimension a(nt),f(nt),p(nt),x(600),hc(282),scr(423),wf(nt)


! values taken from otmgrd2arr
  data dr/.01745329252d0/,irli/1/,samp/1/,irnt/1/

  pi = 4.d0*datan(1.d0)


! transfer the epoch to the array required by spotl routines. The "it" variable is passed via commons to admint
  do i = 1, 5
    itm(i) = timearr(i)
  enddo 

!  Fill the idt (cartwright codes) array in the same order as the amplitudes and phases.  
  do i = 1,n_constit 
    do k = 1,6
      idt(k,i) = cart_codes(k,i)
    enddo
  enddo
   
!  interpolate tidal constituents to larger set of harmonics
  call admint(amp,idt,phase,a,f,p,n_constit,nt)

!  set up for first recursion, and normalize frequencies
  do i=1,nt
    p(i) = dr*p(i)
    f(i) = samp*pi*f(i)/43200.d0
    wf(i) = f(i)
  enddo
 31         irhi = min(irli+599,irnt)
            np = irhi - irli + 1
!
! set up harmonic coefficients, compute tide, and write out
  do i=1,nt
    hc(2*i-1) = a(i)*dcos(p(i))
    hc(2*i)  = -a(i)*dsin(p(i))
  enddo

  call recurs(x,np,hc,nt,wf,scr)
  ocean_ht = x(1)

  deallocate(a)
  deallocate(f)
  deallocate(p)
  deallocate(x)
  deallocate(hc)
  deallocate(scr)
  deallocate(wf)

  return

  end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!**********************************************!!!!!!!!!!!!!!!!
      subroutine spotl_point_setup_v2(in_lat,in_lon,npoints,tidemod,cart_codes,darw_codes,ntmcon,max_con,amp,phase,debug_flag)
!!!!!!!!!!!!!!!**********************************************!!!!!!!!!!!!!!!!

! subroutine to interact with the spotl routines and, eventually, generate an estimate of ocean tide height at a point at a given time
! This is based on a program by Simon McClusky that did the same but for a global grid.
!
! P. Tregoning
! 10 September 2013

  implicit none

  integer*4,        intent(in) :: npoints                   ! number of points whose amp/phase are required
  integer*4,        intent(in) :: max_con                   ! needed to dimension arrays (isn't there a better way of doing this?)
  double precision, intent(in) :: in_lat(npoints)           ! floating point arrays of lat of requested points
  double precision, intent(in) :: in_lon(npoints)           ! floating point arrays of lon of requested points
  character, intent(in)        :: tidemod*10                ! name of requested tide model

  integer*4, intent(out)  :: cart_codes(6,max_con,npoints)  ! cartwright codes as read from the input spotl codes
  character, intent(out)  :: darw_codes(max_con,npoints)*4  ! darwin codes
  real*4, intent(out)     :: amp(max_con,npoints)           ! constituent amplitudes
  real*4, intent(out)     :: phase(max_con,npoints)         ! constituent phases
  logical, intent(in)          :: debug_flag                ! whether to print to screen the phases and amplitudes 

  character                    :: tidemod_ending*10         ! end of name of requested tide model (read from otide.models)
  integer*4                    :: ioerr                     ! i/o error flag
  integer                      :: ntmcon                    ! number of tidal constituents for requested model (read from otide.models)
  character line*200
  character message*200
  logical found
  double precision             :: rlat,rlon                 ! point location in radians
  double precision             :: ocean_ht                  ! computed ocean tidal height
  integer*4                    :: i,j,k,trimlen
  integer*4                    :: timearr(5)                ! yr, doy, hr, min, sec

! declaration of the original spotl variables (don't allow them to be declared implicitly)
      real*8 tlat,blat,wlong,elong,xpoly,ypoly
      real*4 rlato,rlono
      integer*4 latc,longc,i1,i2,np,npoly,polynm,icte
      

! ***@@ come back and delete the variables that are not required 
      integer*4 index,numhed,llst,lgrd,nextf
      real*8 ltg(256),rad2deg,pi,version
      logical*1 ispoly,use
      character*1 modo,fingrd
      character*2 tconst(20)
      character*4 dsym
      character*8 who
      character*16 outgridfile
      character*80 otmod,mdfile,dumm
      character*80 polyf,cmd
      character*50 mdnam
      character*40 stnam
      complex camp,toloc,phasor,oamp,cz
      dimension i1(2),i2(2)
      
!  common block holds polygon file information
      common/polinf/xpoly(500,10),ypoly(500,10),np,npoly(10),polyf,polynm(10),use(10),ispoly
      common /modpar/mdfile,dsym,icte(6),mdnam
      common /modlim/tlat,blat,wlong,elong,latc,longc
! SCM Need this common to move to initialize next file constituent file correctly ...
      common /newmod/nextf
      phasor(toloc) = cmplx(cabs(toloc),57.2958*atan2(aimag(toloc),real(toloc)))
! PT130910: I think I need this one too ...
      real*4 rho
      common/ldens/rho

      data ispoly/.false./
      data cz/(0.,0.)/,oamp/(0.,0.)/
      data stnam/'  Interpolated ocean-tide value         '/


      
! define pi and radians-to-degree converter     
  pi = 4.d0*atan(1.d0)
  rad2deg = 180.d0/pi

!  Get list of tide constituents in the ocean tide model
  open(10,file="otide.models",status='old',iostat=ioerr)
  if(ioerr /= 0)then
    call report_stat('FATAL','UTIL','spotl_setup','otide.models','Cannot open file. Check link.',0)
  endif

! read through the file to find the model that we want
  found = .false.
  do while (.not. found)
    read(10,'(a)',iostat=ioerr,end=1000)line
    if (line(1:6) == tidemod(1:6))then
      found = .true.

! get the ending of the file names
      read(line(11:17),'(a)')tidemod_ending

! get number of constituents
      read(line(18:20),'(i3)')ntmcon
      do i=1,ntmcon
        tconst(i) = line(22+(i-1)*3:22+(i-1)*3+1)
      enddo
    endif
  enddo

1000 continue

! check that the model was found
  if (.not. found) then
    write(message,'(a,a,a)')"Tide model (",tidemod,") not found in file"  
    call report_stat('FATAL','UTIL','spotl_setup','otide.models',message,0)
  else
    close(10)
    if(ntmcon == 12) then
      write(message,'(i4,a,12(1x,a2),a,a)')ntmcon,' constituents (',(tconst(i),i=1,ntmcon),') found for model ',tidemod
    else
      write(message,'(i4,a,a)')ntmcon,' constituents found for model ',tidemod
    endif
    if(debug_flag)call report_stat('STATUS','UTIL','spotl_setup',' ',message,0)
  endif   

! don't know yet what this does .....
      fingrd = 'F'
        
! Loop over constituent files  
  do i = 1,ntmcon
    write(mdfile,'(a,a2,a1,a,a1,a)')'spotl/',tconst(i),'.',tidemod(1:trimlen(tidemod)),'.',tidemod_ending(1:trimlen(tidemod_ending))
    nextf = 1        ! this ensures that the constituent file stays in a particular state (not sure what though ...?)
    do k = 1,npoints
      rlato = in_lat(k)
      rlono = in_lon(k)
      call ocmodl(rlato,rlono,fingrd,camp)
      if(camp.ne.cz) camp=camp/rho
      camp = conjg(camp)
      camp = phasor(camp)
      amp(i,k) = real(camp)
      phase(i,k) = aimag(camp)
! save off the cartwright and darwin codes
      do j=1,6
        cart_codes(j,i,k) = icte(j)
      enddo
      darw_codes(i,k) = dsym

      write(message,'(a,a,f12.4,a,f12.4,a,2f12.4)')darw_codes(i,k),'Ampl:',amp(i,k),'   Phase:',phase(i,k),' coords',rlato,rlono
      if(debug_flag .and. k == 1 )call report_stat('STATUS','UTIL','spotl_setup',' ',message,0)
    enddo
  enddo
  if(.not. debug_flag )call report_stat('STATUS','UTIL','spotl_setup',' ',"     Computed amplitudes and phases",0)

  return

  end subroutine spotl_point_setup_v2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine spotl_get_ocean_height(npoints,first,lat,lon,timearr,tidemod,ocean_ht,debug_flag)

! subroutine to interact with the two routines that interact with spotl 
!
! P. Tregoning
! 11 September 2013

  implicit none

  integer*4,        intent(in)    :: npoints                    ! number of points in the lat/lon vectors
  double precision, intent(in)    :: lat(npoints),lon(npoints)  ! coordinates of point(s) at which to estimate tidal height
  integer*4,        intent(in)    :: timearr(5)                 ! yr, doy, hr, min,sec of epoch to estimate tide
  character,        intent(in)    :: tidemod*10                 ! requested tide model
  logical,          intent(inout) :: first                      ! flag for first time into subroutine
  logical,          intent(in)    :: debug_flag                 ! whether to print to screen the phases and amplitudes
  double precision, intent(out)   :: ocean_ht(npoints)          ! computed ocean heights

! local variables and things requred for spotl interaction
  integer*4, allocatable        :: cart_codes(:,:,:)          ! cartwright codes as read from the input spotl codes
  character, allocatable        :: darw_codes(:,:)*4          ! darwin codes
  real*4,    allocatable        :: amp(:,:)                   ! constituent amplitudes
  real*4,    allocatable        :: phase(:,:)                 ! constituent phases
  integer*4                     :: ntmcon                     ! number of tidal constituents
  integer*4                     :: max_con                    ! maximum number of constituents
  integer*4 i,j
  character message*200


  save cart_codes,darw_codes,amp,phase


! now get the amplitudes and phases for this point
  if(first) then
! allocate the arrays
  max_con = 20
    allocate(cart_codes(6,max_con,npoints))
    allocate(darw_codes(max_con,npoints))
    allocate(amp(max_con,npoints))
    allocate(phase(max_con,npoints))
    call report_stat('STATUS','UTIL','get_ocean_height',' ','Computing amplitudes/phases',0)
    call spotl_point_setup_v2(lat,lon,npoints,tidemod,cart_codes,darw_codes,ntmcon,max_con,amp,phase,debug_flag) 
  endif

! next, compute the tidal height
  write(message,'(a,5i5,i8)')'Computing ocean heights for epoch',timearr,npoints
  if(debug_flag)call report_stat('STATUS','UTIL','get_ocean_height',' ',message,0)

!$OMP PARALLEL DO  shared(timearr,ntmcon,cart_codes,darw_codes,amp,phase,ocean_ht)
  do i=1,npoints
    if(amp(1,i) > 0.0 .and. amp(2,i) > 0.0)then     ! don't compute if the point is over land
      call spotl_point_calc(timearr,ntmcon,cart_codes(:,:,i),darw_codes(:,i),amp(:,i),phase(:,i),ocean_ht(i))
    endif
  enddo
!$OMP END PARALLEL DO

! set first to false
  if(first) first = .false.

! and we're done
  return
  end subroutine spotl_get_ocean_height

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

