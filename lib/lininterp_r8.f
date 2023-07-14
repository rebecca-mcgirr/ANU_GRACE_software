      subroutine lininterp_r8( a,b,val_time,valtab,obstime,nrow,ncol
     .                    , valinterp,exttim )

c  performs a linear interpolation to calculate the value
c  at a particular time. Written specifically to interpolate
c  the tabulated 6-hourly atmospheric pressure loading 
c  displacements but could work equally well for any other
c  circumstance.
c
c  INPUT:
c    a,b         : dimensions of the matrices in the calling program
c    val_time    : array of time tags of tabulated values
c    valtab      : array of corresponding tabulated values
c    obstime     : time to which the tabulated values should be
c                  interpolated
c    nrow        : number of rows of tabulated values     
c    ncol        : number of independent variables to be
c                  interpolated 

c
c  OUTPUT:
c     valinterp  : interpolated values 
c     exttim     : if requested interpolation is outside the range
c                  of the table, this variable will have the time
c                  (neg or pos) by which the range is exceeded;
c                  if within range, extflg = 0.
c                  
c
c  P. Tregoning
c  8 January 2004

      implicit none

      integer*4 a,b,ncol,nrow,t1,t2,i

      real*8 val_time(a),valtab(a,b),exttim,valinterp(ncol)
      real*8 obstime,dval,dt,tstep,delta,eps
      character*256 message

      logical found
c      logical debug
c      data debug/.false./

                                    
c     eps is a small number tolerance to allow use of the end values
c     currently set to be 90s if times are in days
      data eps/.000001d0/
                  

c   see if the requested time is outside the range of the array   
      exttim = 0.      
      if( obstime.lt.(val_time(1)-eps) )  then
        exttim = obstime - val_time(1)  
        do i=1,ncol
         valinterp(i) = valtab(1,i)
        enddo
        return
      elseif( obstime.gt.(val_time(nrow)+eps) ) then
        exttim =  obstime - val_time(nrow) 
        do i=1,ncol
          valinterp(i) = valtab(nrow,i)
        enddo
        return
      endif

c   run through the time array to find which two times the requested
c   time lies between
      found = .false.
      t1 = 0
      do while (.not.found)
        t1 = t1+1
        if(t1+1.gt.nrow)then
          write(message,'(a,f15.7,a,2f15.7)')
     .     'Requested time',obstime,' outside tabulated values. Min/max'
     .    ,val_time(1),val_time(nrow)
          call status_update('FATAL','lib','lininterp_r8',' '
     .                  ,message,0)
        endif
c PT050415: there can be a roundoff problem if obstime is infinitessimally
c           less that the first val_time(t1) .... ! 
c RWK070111: Modified to allow interpolation within eps of last value
        if((dabs(obstime-val_time(t1)).lt.0.01d0.or.
c     .   obstime.ge.val_time(t1)).and.obstime.lt.val_time(t1+1))then
     .    obstime.ge.val_time(t1)).and.obstime.lt.(val_time(t1+1)+eps)) 
     .          then
          found = .true.
        endif
c PT061201: trap the case that the obs time is exactly the last time in the array
c RWK070111: Now allow it to be eps later than last time
c        if(dabs(obstime-val_time(t1+1)).lt.0.00001d0)then 
        if(dabs(obstime-(val_time(t1+1)+eps)).lt.0.00001d0)then 
          found = .true.
        endif
      enddo    

      t2 = t1 + 1


c  now linearly interpolate all the parameters between these two times
      do i = 1,ncol
        dt = obstime-val_time(t1)
        tstep = val_time(t2)-val_time(t1)
        dval = valtab(t2,i)-valtab(t1,i)

        delta = dval * dt / tstep
        valinterp(i) = valtab(t1,i)+delta

c DEBUG  
c      if(a .ne.9)print*,'obstime,t1,t2',obstime,t1,t2,val_time(t1)
c     .        ,val_time(t2),valtab(t1,1),valtab(t2,1)
      enddo

      return
      end

