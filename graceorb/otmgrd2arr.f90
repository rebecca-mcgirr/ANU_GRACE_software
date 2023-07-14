   subroutine otmgrd2arr(timearr, grd_file, ocean_hts)  
!c  Compute total tide height by summing harmonics; the amp and phase of these are
!c  specified in subroutine admint
!   input:
!	it(1) = year
!	it(2) = day of year
!!       imonth = month
!!       iday = day of month
!       it(3) = hour
!       it(4) = min
!       it(5) = sec
!       grd_file = name of grd file 
!c
!c  ** calls subroutines - admint,tdfrph,recurs,spline,shells
!c  ** calls functions - mday, upperc
!c
!
! PT130903: this routine takes a total of 54" of 71" (for 2000 epochs). So let's try to speed it up !!!
      use gauleg_mod

      logical notfound,land

!      integer*4 nlat,nlon,nin,i,hindex,gindex,numhed,icount
      integer*4 nin,i,hindex,gindex,numhed,icount

      real*8 version
!      real*8, allocatable :: ocean_hts(:)
      real*8 :: ocean_hts(nlon*nlat)
!      real*8, dimension(nlat*nlon) :: ocean_hts

      character*1 ipht
      character*2 darw_code(10),tconst(20),upperc
      character*8 otmod,who
      character*16 grd_file
      character*40 dumm

      double precision f,p,dr,pi,scr

      complex camp,cz,oamp

      dimension a(141),f(141),p(141),x(600),hc(282),scr(423),wf(141)
      dimension idt(6,20),cart_code(6,10),amp(20),phase(20)

      common/date/it(5)

      data nt/141/,dr/.01745329252d0/,pi/3.1415926535898d0/,irli/1/
      data lui/5/,luo/6/,maxc/20/,nin/1/,lgrd/50/,samp/1/,irnt/1/
      data cz/(0.,0.)/,oamp/(0.,0.)/    

      integer, dimension(5) :: timearr
      integer ::  nlattmp, nlontmp
      integer, dimension(3) :: tarray

      integer ithread  !,omp_get_thread_num
!c        
!c Cartwright Tidal constituent codes
      data cart_code/ & 
              0, 2, 0, 0, 0, 0,&
              0, 1, 0,-1, 0, 0,&
              1, 1, 0, 0, 0, 0,&      
              1,-1, 0, 0, 0, 0,&
              1, 1,-2, 0, 0, 0,&
              1,-2, 0, 1, 0, 0,&   
              2, 0, 0, 0, 0, 0,&
              2, 2,-2, 0, 0, 0,&
              2,-1, 0, 1, 0, 0,&   
              2, 2, 0, 0, 0, 0 /
!c
!c Darwin Tidal constituent codes
      data darw_code/ 'Mf', 'Mm', 'K1', 'O1', 'P1', & 
                     'Q1', 'M2', 'S2', 'N2', 'K2' /
!c
!c INFO
!c     Approx Amplitudes 0.07, 0.04 m
!c     0 2 0 0 0 0 Mf
!c     0 1 0 -1 0 0 Mm
!c     Approx Amplitudes 0.37, 0.26, 0.12, 0.05 m
!c     1 1 0 0 0 0 K1
!c     1 -1 0 0 0 0 O1
!c     1 1 -2 0 0 0 P1
!c     1 -2 0 1 0 0 Q1

!c     Approx Amplitudes 0.63, 0.29, 0.12, 0.08 m
!c     2 0 0 0 0 0 M2
!c     2 2 -2 0 0 0 S2
!c     2 -1 0 1 0 0 N2
!c     2 2 0 0 0 0 K2        
!c
!c  Read the command line           
!      if(iargc().lt.7.or.iargc().gt.8) then
!	 write(luo,100)
! 100  format('Usage: ogrd2array y [d | m d] h m s grid_file',/,
!     1  '  harmonics file read from grid_file',/,
!     2  '  results written to standard output')
!	 stop
!      endif
!      call getarg(1,dumm)
!      read(dumm,102) it(1)
! 102  format(i4)
!      if(iargc().eq.6) then
!        call getarg(2,dumm)
!        read(dumm,102) it(2)
!        nb=0
!      endif
!      if(iargc().eq.7) then
!        call getarg(2,dumm)
!        read(dumm,102) imonth
!        call getarg(3,dumm)
!        read(dumm,102) iday
!        nb=1
!        it(2) = iday + mday(it(1),imonth)
!      endif
!      call getarg(nb+3,dumm)
!      read(dumm,102) it(3)
!      call getarg(nb+4,dumm)
!      read(dumm,102) it(4)
!      call getarg(nb+5,dumm)
!      read(dumm,102) it(5)
!      call getarg(nb+6,dumm)
!      read(dumm,104) grd_file
! 104  format(a)
!      print*,'it, grd_file: ',it,grd_file

    do i = 1, 5
      it(i) = timearr(i)
    enddo 
      
!c 
!c  Open the input binary ocean grid file
      open(unit=lgrd,file=grd_file,status='unknown' &
         ,access='direct',form='unformatted',recl=8,iostat=ioerr)
     
!print*, lgrd, grd_file, ioerr
 
      version = 0.0   
      who = '' 
      numhed = 0
!c      
!c  Read binay grid file header     
!       call date_and_time(VALUES=t1)
!      call itime(tarray)
!      print*, tarray
!   do j = 1, 1000
      read(lgrd,rec=1)version
      read(lgrd,rec=2)numhed    
      read(lgrd,rec=3)otmod
      read(lgrd,rec=4)who
      read(lgrd,rec=5)nlattmp,nlontmp
      read(lgrd,rec=6)nin
      do i = 1,nin
        hindex = 6 + i
        read(lgrd,rec=hindex) tconst(i)	
      enddo   
!c
!c  Fill the idt (cartwright codes) array in the same order as the
!c  tidal constituents are stored on the binary grid file.  
      do i = 1,nin 
        j = 0
        notfound = .true. 
        do while ( notfound )
          j = j + 1
          if ( upperc(tconst(i)) .eq. upperc(darw_code(j)) ) then
            do k = 1,6
              idt(k,i) = cart_code(k,j)
            enddo
            notfound = .false.
          endif
        enddo
      enddo
!c 
!c  Loop over grid file nodes summing up constituent contributions and 
!c  computing total tide height at each grid node point. 
      icount = 0
      gindex = hindex

!! !$OMP PARALLEL DO default(shared) private(icount,gindex,ilon,i,camp,land,a,f,p,irhi,np,hc,x,wf,scr,amp,phase) &
!! !$OMP&            firstprivate(nin,nt,idt,hindex,cz,dr,samp,pi,irli) 
      do ilat = 1,nlat
        do ilon = 1,nlon
          land = .true.
          irli = 1
! PT130903: we can derive icount as a function of ilat, ilon
!          icount = icount + 1
          icount = (ilat-1)*nlon + ilon
!  print*,'thread,nin,nlon,nlat,ilat,ilon,icount',ithread,nin,nlon,nlat,ilat,ilon,icount
          do i = 1,nin
            gindex = hindex + (ilat-1)*nlon*nin + (ilon-1)*nin + i
            read(lgrd,rec=gindex) camp
            if ( camp .ne. cz ) land = .false.
            amp(i) = real(camp)
            phase(i) = aimag(camp)
          enddo

!c
!c  interpolate tidal constituents to larger set of harmonics
          if ( .not. land ) then
            call admint(amp,idt,phase,a,f,p,nin,nt)
!c
!c  set up for first recursion, and normalize frequencies
            do i=1,nt
              p(i) = dr*p(i)
              f(i) = samp*pi*f(i)/43200.d0
              wf(i) = f(i)
            enddo
 31         irhi = min(irli+599,irnt)
            np = irhi - irli + 1
!c
!c set up harmonic coefficients, compute tide, and write out
            do i=1,nt
              hc(2*i-1) = a(i)*dcos(p(i))
              hc(2*i)  = -a(i)*dsin(p(i))
            enddo
! debug

            call recurs(x,np,hc,nt,wf,scr)
              ocean_hts(icount) = x(1)
          else
            ocean_hts(icount) = 0.d0
          endif

! DEBUG
!  if(ilat == 105 .and. ilon == 170)then
!    print*,timearr,-90.d0 + dble(ilat)*180.d0/dble(nlat),ilon*360.d0/512.d0,ocean_hts(icount), "-16:120 OTM"
!  endif
!  if(ilat == 131 .and. ilon == 466)then
!    print*,timearr,-90.d0 + dble(ilat)*180.d0/dble(nlat),ilon*360.d0/512.d0,ocean_hts(icount), "2:328 OTM"
!  endif
        enddo
!c  go and do the next latitude band	
      enddo
!! !$OMP END PARALLEL DO

!c  Close the binary ocean grid file      
      close(lgrd)

      return
      end

      function mday(iy,m)
!c
!c  finds day number of day before start of month m, of year iy, in
!c   Gregorian intercalation
!c
!print*, "mday"
      leap = 1 - (mod(iy,4)+3)/4
      if(mod(iy,100).eq.0.and.mod(iy,400).ne.0) leap=0
      mday = mod((367*(m-2-12*((m-14)/12)))/12+29,365) + leap*((9+m)/12)
      return
      end
     
