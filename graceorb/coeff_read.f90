    subroutine coeff_read

! Read into an array the spherical harmonic coefficients for the Earth's
! mean static gravity field

! MODIFIED: APP 121004
! To allow user defined value for the gravity coefficient file

    use coeff_mod
    use inmod_mod
    use spherhar_mod
    implicit none
    integer*4 :: lucoeff, ioerr, ierr,i, counter, ideg, iord
    real*8 :: temp1, temp2, cval, sval

!  open file with coefficients
      lucoeff = 15
! PT151014: this should be "status='old'" so that it fails if file is not there!
      open(unit=lucoeff,file=static_grav_coeff_file,status='old',iostat=ioerr)
!      open(unit=lucoeff,file='coeffs_200',status='unknown',iostat=ioerr)
      if(ioerr.ne.0)then
        write(message,'(a,a)')'Error opening coefficients file: ',static_grav_coeff_file
        call status_update('FATAL','GRACEORB','coeff_read',' ',message,0)
      endif
      ioerr = 0
      ndeg = 0
      counter = 0

!  find the end of file and record the max degree
      do while (ioerr.eq.0)
          counter = counter + 1
          read(lucoeff,*,iostat=ioerr,end=1236)ideg, iord, temp1, temp2
          if (counter .eq. 1) then
!         print*, "spherical harmonic model starting degree and order", ideg, iord
            if (ideg.ne.2.or.iord.ne.0) then
              call status_update('FATAL','GRACEORB','coeff_read',' ',"Spherical harmonic model must start from degree 2, order 0",0)
            endif
          endif
          if(ideg.gt.ndeg)ndeg = ideg
      enddo
1236  continue
      close(lucoeff)
      ierr = 0

! PT151014: add a trap for when a file of zero size has been opened
      if(ndeg < 2)then
        call status_update('FATAL','GRACEORB','coeff_read',static_grav_coeff_file &
                  ,"Static gravity field spherical harmonic model has no entries",0)
      endif
      
!     read in the coefficients and store
      ioerr = 0
      open(unit=lucoeff,file=static_grav_coeff_file,status='unknown',iostat=ioerr)
!      open(unit=lucoeff,file='coeffs_200',status='unknown',iostat=ioerr)
      do i = 1, (ndeg+1)*(ndeg+2)/2-3
        read(lucoeff,*, iostat=ioerr) ideg,iord,cval,sval
        if(ioerr.eq.0)then
          if(ideg.gt.ndeg)ndeg = ideg
          meancoefC(ideg,iord) = cval
          meancoefS(ideg,iord) = sval
        endif
      enddo
1234  continue
      write(message,'(a,a20,a,i4)')"  Input sph harm static gravity field, ",static_grav_coeff_file(1:20) &
                                 ,", is up to degree ",ndeg
      call status_update('STATUS','GRACEORB','coeff_read',' ',message,0)
      max_statgrav_degree = ndeg
      close(lucoeff)

    return 
    end
