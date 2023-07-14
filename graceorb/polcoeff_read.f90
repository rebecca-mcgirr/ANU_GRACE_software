    subroutine  polcoeff_read

! reads in the spherical harmonic coefficients for the ocean pole tide

    use coeff_mod
    use inmod_mod

    implicit none
    real*8 :: Areal, Breal, Aimg, Bimg
    integer*4 :: polcoefffile, ioerr, ierr,i, counter, ideg, iord,ndeg
    real*8 :: temp1, temp2, temp3, temp4, temp5, temp6
    character*100 :: line

!  open pole tide coefficient file - hard wired as desaicoeff.dat 
      polcoefffile = 17 
      open(unit=polcoefffile,file=pole_coeff_file,status='unknown',iostat=ioerr)
!      open(unit=polcoefffile,file='desaicoeff.dat',status='unknown',iostat=ioerr)
      if(ioerr.ne.0)then
        print*,'Error opening coefficients file: '
        stop
      else
! read and discard the header line
        read(polcoefffile,'(a)')line
      endif
      ioerr = 0
      ndeg = 0
      counter = 0

!  determine length of file and record maximum degree 
      do while (ioerr.eq.0)
          counter = counter + 1
          read(polcoefffile,*,iostat=ioerr,end=36)ideg, iord, temp1, temp2, temp3, temp4
          if (counter .eq. 1) then
            if (ideg.ne.1.or.iord.ne.0) then
              print*, "Ocean pole tide A and B coeff must start from degree 1, order 0"
              stop
            endif
          endif
          if(ideg.gt.ndeg)ndeg = ideg
      enddo
36  continue
      close(polcoefffile)
      ierr = 0

!  read in the coefficients and store
      ioerr = 0
      open(unit=polcoefffile,file=pole_coeff_file,status='unknown',iostat=ioerr)
!      open(unit=polcoefffile,file='desaicoeff.dat',status='unknown',iostat=ioerr)
      do i = 1, (ndeg+1)*(ndeg+2)/2-1
        read(polcoefffile,*, iostat=ioerr) ideg,iord,Areal, Breal, Aimg, Bimg
        if(ioerr.eq.0)then
          if(ideg.gt.ndeg)ndeg = ideg
          AR(ideg,iord) = Areal
          BR(ideg,iord) = Breal
          AI(ideg,iord) = Aimg
          BI(ideg,iord) = Bimg
        endif
      enddo
1234  continue
      write(message,'(a,i4)')"     Input ocean pole tide model is up to degree ",ndeg
      call status_update('STATUS','GRACEORB','polcoff',' ',message,0)
      close(polcoefffile)
    return
    end

    


