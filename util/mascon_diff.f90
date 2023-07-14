      program mascon_diff
!     *******
!
!  calculates the difference between mascon loads in two fit/vcv format files
!
!  usage: mascon_diff <target_file> <reference_file>
!
!  input : target and reference files, both containing mascon load estimates
!               - may be in .fit file format
!                        or gracefit .vcv file format 
!                        or addnorm  .vcv file format 
!		- file formats determined by gracefit
!
!  output: difference file mdiff.fit ( =  fit-ref ) in .fit file format 
!        
!  also reports:
!          norms (rms) of mascon loads in the target and reference fit files
!          norm of differences between mascon loads in the two files
!          norm of the differences weigted by sigmas in the target fit file
!
!  H.McQueen - 180502

    implicit none

    real (kind=8), dimension(:,:), allocatable    :: sat		! target satellite ICs
    real (kind=8), dimension(:,:), allocatable    :: rsat		! reference satellite ICs
    real (kind=8), dimension(:,:), allocatable    :: val		! target mascon load
    real (kind=8), dimension(:,:), allocatable    :: rval		! reference mascon load
    character(180), dimension(:), allocatable 	  :: line	        ! output lines

    real (kind=8) sig            ! fit sigma
    real (kind=8) vmin,vmax      ! min,max of mascon load
    real (kind=8) rmin,rmax      ! min,max of reference mascon load
    real (kind=8) dmin,dmax      ! min,max of difference
    real (kind=8) wmin,wmax      ! min,max of sigma-weighted difference
    integer mic,mmasc,mtmasc     ! expected number of ICs, mascons, tidal mascons
    integer mln                  ! expected maximum no. of lines
    character*120 fit_file,ref_file,dif_file     ! filenames
    
    
! working variables
    
    real (kind=8) sum1,sum2,fnorm,rnorm
    real (kind=8) dif,wdif,sum3,sum4,sd,wsd
    integer mp
    integer nline
    integer i,j,ios,len,lun
    integer nhd,nic,nmasc,nmasc_ref
    logical set
    character*3 fext,rext,ftype,rtype
    character*16 ftag,rtag

1   format(a8,i6," mascons, ",i6," lines")
2   format(a20,f12.5,a12,f10.3,"  / ",f10.3)
3   format(a16)    ! PT180831: increased from a11 to a16
4   format(a180)
5   format(1x,f21.5,1x,f11.5,1x,f18.5,1x,f11.5,1x,f8.1)
7   format(i3)
8   format(9x,i4)
9   format(i8,$)

    write(*,'(a)') "mascon_diff: subtracting reference file"
    
! parse command line
    call getarg(1,fit_file)
    if(fit_file.eq."") then
        write(*,*) " usage: mascon_diff <target_file> <reference_file>"
        stop " no target file specified"
    else
       write(*,*) "input file: ",fit_file
        len=len_trim(fit_file)
        fext=fit_file(len-2:len)
	open(10,file=fit_file,status="old")
        read(10,3) ftag
        rewind(10)
    endif

    call getarg(2,ref_file)
    if(ref_file.eq."") then
        write(*,*) " usage: mascon_diff <target_file> <reference_file>"
	stop " no reference file specified"
    else
        write(*,*) "reference file: ",ref_file
        len=len_trim(ref_file)
        rext=ref_file(len-2:len)
	open(11,file=ref_file,status="old")
        read(11,3) rtag
        rewind(11)
     endif
!     write(*,*) fext,"|",ftag
!     write(*,*) rext,"|",rtag
     
! identify file formats
! PT180831: update to new 1st line information in .fit, .vcv and addnorm files
print*,"xx",ftag,"xx"
     ftype="nul"
     if(fext.eq."fit".and.ftag.eq."V2 GRACEFIT  FIT") ftype="fit"
     if(fext.eq."vcv".and.ftag.eq."V2 GRACEFIT  VCV") ftype="vcv"
     if(fext.eq."vcv".and.ftag.eq."V2 ADDNORM   VCV") ftype="add"
     if(fext.eq."fit".and.ftag.eq."V2 ADDNORM   FIT") ftype="fit"

print*,"xx",rtag,"xx"
     rtype="nul"
     if(rext.eq."fit".and.rtag.eq."V2 GRACEFIT  FIT") rtype="fit"
     if(rext.eq."vcv".and.rtag.eq."V2 GRACEFIT  VCV") rtype="vcv"
     if(rext.eq."vcv".and.rtag.eq."V2 ADDNORM   VCV") rtype="add"
     if(rext.eq."fit".and.rtag.eq."V2 ADDNORM   FIT") rtype="fit"

!     write(*,*) fext," ",ftag," ",ftype
     if(ftype.eq."nul") STOP "error: fit_file type unrecognised"
     if(rtype.eq."nul") STOP "error: ref_file type unrecognised" 

     write(*,*) "fit_file type: ",ftype
     write(*,*) "ref_file type: ",rtype

     
! initialize parameters

     mic=24
! mascon dimension arbitrary limit 9999
    mmasc=9999
    mtmasc=0
    mp=mic+mmasc+mtmasc
    mln=40+mp+mp/12+1
    
    i=0
    j=0
    ios=0
    sum1=0.d0
    sum2=0.d0
    sum3=0.d0
    sum4=0.d0
    set=.true.
    write(*,1) "limit ",mmasc,mln
    open(99,file="junk",status="unknown")

! allocate arrays

    allocate(sat(mic,3))
    allocate(rsat(mic,3))
    allocate(val(mmasc,3))
    allocate(rval(mmasc,3))
    allocate(line(mln))

! read ref file data

    lun=11
    call readfit(lun,rtype,line,mln,rsat,mic,rval,mmasc,nhd,nic,nmasc_ref)
    write(*,*) "ref_file n_header,n_IC,n_masc: ",nhd,nic,nmasc_ref

! read fit file data

    lun=10
    call readfit(lun,ftype,line,mln,sat,mic,val,mmasc,nhd,nic,nmasc)
    write(*,*) "fit_file n_header,n_IC,n_masc: ",nhd,nic,nmasc

    if(nmasc.ne.nmasc_ref) then
       stop "Error: mascon number mismatch - mdiff.fit not generated"
       endif
    
! output mascon difference file

    dif_file="mdiff.fit"
    open(12,file=dif_file,status="unknown")

! header lines

    do i=1,nhd
        write(12,"(a)") trim(line(i))
        enddo
     
! IC differences

     do i=1,nic
        j=nhd+i
        write(line(j)(26:99),5) 0.,0.,sat(i,2)-rsat(i,2),sat(i,3)
        write(12,"(a)") trim(line(j))
     enddo

! mascon differences and build statistics

     do i=1,nmasc
        j=nhd+nic+i
        dif=val(i,2)-rval(i,2)
        sig=val(i,3)
        wdif=dif/sig
        write(line(j)(26:99),5) 0.,0.,dif,val(i,3),wdif
        write(12,"(a)") trim(line(j))

        if(set) then
            vmin=val(i,2)
            vmax=vmin
            rmin=rval(i,2)
            rmax=rmin
            dmin=dif
            dmax=dmin
            wmin=wdif
            wmax=wmin
            set=.false.
            endif
 
        vmin=dmin1(vmin,val(i,2))
        vmax=dmax1(vmax,val(i,2))
        rmin=dmin1(rmin,rval(i,2))
        rmax=dmax1(rmax,rval(i,2))
        dmin=dmin1(dmin,dif)
        dmax=dmax1(dmax,dif)
        wmin=dmin1(wmin,wdif)
        wmax=dmax1(wmax,wdif)
                

        sum1=sum1+val(i,2)**2
        sum2=sum2+rval(i,2)**2
        sum3=sum3+dif**2
        sum4=sum4+(dif/sig)**2

        write(99,'(4f12.5)') val(i,2),rval(i,2),dif,wdif

	enddo

    close(10)
    close(11)

    nline=nhd+nic+nmasc
    write(*,1) "found  ",nmasc,nline

! write statistics
    
    fnorm=dsqrt(sum1/nmasc)
    rnorm=dsqrt(sum2/nmasc)
    sd=dsqrt(sum3/nmasc)
    wsd=dsqrt(sum4/nmasc)
    write(*,2) "fit norm:",fnorm,"range:",vmin,vmax
    write(*,2) "ref norm:",rnorm,"range:",rmin,rmax
    write(*,2) "difference norm:",sd,"range:",dmin,dmax
    write(*,2) "weighted diff norm:",wsd,"range:",wmin,wmax


    close(12)
    close(99)

!    call status_update('STATUS','UTIL','mascon_diff',"mascon__diff","End program",0)

    stop
    end


!---------------------------------------------


    subroutine readfit(lun,type,line,mln,sat,mic,val,mmasc,nhd,nic,nmasc)

    integer, intent(in)  :: lun,mln,mic,mmasc
    character(3), intent(in)  :: type

    integer, intent(out) :: nhd,nic,nmasc
    real (kind=8), intent(out), dimension(mic,3)       :: sat	! fit satellite ICs
    real (kind=8), intent(out), dimension(mmasc,3)     :: val	! fit mascon values
    character(180), intent(out), dimension(mln)        :: line	! output lines

    ! working variables

    integer i,j,k,ios
    character*180 txt
    character*25 str
    character*1 satid

4   format(a180)
7   format(2i3,5x,a1)
8   format(9x,i4)

    nhd=0
    nic=mic
    nmasc=0
 
    j=0
    ios=0

! read header
    
    do while ( ios==0 )
       j=j+1
       read(lun,4,iostat=ios) line(j)
       if(len_trim(line(j)).gt.0) then
          read(line(j),*) str
          if(str.eq."PARAMETER") go to 200
          endif
! PT180831: increase this to 1000 header lines .....
       if(j.gt.1000) stop "error: no end of header after 100 lines" 
       enddo
200    nhd=j
       

! loop over ICs and mascons
       
     read(lun,4,iostat=ios) txt
     l=j
     do while ( ios==0 )
        l=l+1

! read ICs
!   and scale position and velocity ICs if input is of type vcv
!   (pos and vel ICs are in mm and mm/s in fit files and addnorm vcv files,
!                    but in m and m/s in gracefit vcv files)

        if(txt(6:8).eq."SAT") then
! fit file or gracefit vcv file
           j=j+1
           read(txt,"(a)") line(j)(1:25)
           line(j)(26:100)=" "
           read(txt,7) i
           if(i.gt.24) then
              write(*,*) "Error: ",i," exceeded expected ICs ",24
              stop
              endif
           nic=max(nic,i)
           if(type.eq.'fit') then
              read(txt(26:100),*) sat(i,1),dum,sat(i,2),sat(i,3)
           else
              read(txt(26:100),*) sat(i,1),sat(i,2),sat(i,3)
              if(mod(i-1,12)+1.le.6) then
                 do k=1,3
                    sat(i,k)=sat(i,k)*1000
                    enddo
                 line(j)(20:25)="m"//line(j)(20:24)
                 endif
              endif

        else if(txt(8:10).eq."SAT") then
! addnorm vcv file
           if(type.ne.'add') then
              write(*,*) "Error: file type incompatibility, expecting ADDNORM vcv file"
              stop
           endif
           read(txt,7) k,i,satid
           if(k.eq.1) then
              if(satid.eq."B") i=i+12
              if(i.gt.24) then
                 write(*,*) "Error: ",i," exceeded expected ICs ",24
                 stop
                 endif
              nic=max(nic,i)
              j=j+1
              read(txt,"(a)") line(j)(1:28)
              read(txt(29:100),*) sat(i,1),sat(i,2),sat(i,3)
              write(str(1:3),"(i3)") i
              str(4:4)="."
              str(5:12)=line(j)(7:14)
              str(13:25)=line(j)(16:28)
              line(j)(1:25)=str
              line(j)(26:100)=" "              
              endif

! read mascon values
     
         else if(txt(8:9).eq."MC") then
            j=j+1   
            nmasc=nmasc+1
            read(txt,"(a)") line(j)(1:25)
            line(j)(26:100)=" "
            read(txt,8) i
            if(i.gt.mmasc) then
                write(*,*) "Error: ",i," exceeded expected mascons",mmasc
                stop
                endif
            if(type.eq."fit") then
                read(txt(26:100),*) val(i,1),dum,val(i,2),val(i,3)
            else
                read(txt(26:100),*) val(i,1),val(i,2),val(i,3)
                endif

            endif

        if(l.gt.mln) then
!           write(*,*) "line limit: exiting read at line ",l
!           write(*,*) "lines output: ",j
            return
            endif

        read(lun,4,iostat=ios) txt

        enddo

        
     return
     end subroutine readfit