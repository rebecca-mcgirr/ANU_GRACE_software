      program scan_regularisation_file
!     *******
!
!  Scans mascon regularisation file and associated mascon file
!   - counts mascons with each distinct regularisation value
!   (the mascon.reg file specifies a regularisation value for each primary mascon)
!
!  usage: scan_regfile <mascon_file>
!
!
! main arrays and variables:
! --------------------------
! pm/sm/tm - primary/secondary/ternary mascon prefix
! (for pm below, also sm...,tm...)
! pmlat - primary mascon centre latitude
! pmlon - primary mascon centre longitude
! pmr - primary mascon geocentric radius
! pma - primary mascon area
! pmh - primary mascon altitude
! pmg - primary mascon geoid height
! pmden - primary mascon water density (salt/fresh)
! pcland - percent land in primary mascon
! np,ns,nt - number of p/s/t mascons
! npex,nsex,ntex - expected number of p/s/t mascons
! sap,tap,tas - secondary associated primary etc.
! nsip,ntip,ntis - number of secondaries in a primary, etc.
! maxsip,maxtip,maxtis - maximum # of secondaries in a primary mascon, etc.
! sscol,tpcol,tscol - secondary secondary colour, etc.
! (p/s/t)region - primary region, etc.
! (p/s/t)type - primary type, etc.
! lsip,ltip,ltis - list of secondaries in primaries, etc.
!
!        
!        MHR= Max header records
!        MRC= Max char in .reg header line
!        MRV= Max regularisation values
!        rv - regularisation value array
!        irv = nint(rv)
!        nrv - regularisation value index
!        regd - diagonal element of regularisation matrix
        
    implicit none
    real (kind=8), parameter    :: pi = 3.141592653589793d0
    integer, parameter      :: MHR = 200, MRC=300, MRV=20
    
	real (kind=8), dimension(:), allocatable    :: pmlat,pmlon,pmr,pma,pmh,pmg,pmden,pcland
!	real (kind=8), dimension(:), allocatable    :: smlat,smlon,smr,sma,smh,smg,smden
!	real (kind=8), dimension(:), allocatable    :: tmlat,tmlon,tmr,tma,tmh,tmg,tmden
!	real (kind=8), dimension(:), allocatable    :: dspmax,dtpmax,dtsmax
	integer, dimension(:), allocatable          :: nsip,ntip,tidemask
!	integer, dimension(:), allocatable          :: ntis,sap,sscol
!	integer, dimension(:), allocatable          :: tap,tas,tpcol,tscol
!	integer, dimension(:), allocatable          :: isip,itip,kips,kpsc
!	integer, dimension(:), allocatable          :: itis,kist,kstc
!	integer, dimension(:), allocatable          :: hsind,hpind
!	integer, dimension(:), allocatable          :: scheck,tcheck

!    integer, dimension(:,:), allocatable 		:: lsip,ltip,ltis

	character*20, dimension(:), allocatable	   :: pregion
!	character*20, dimension(:), allocatable	   :: sregion
!	character*20, dimension(:), allocatable	   :: tregion
	character(6), dimension(:), allocatable    :: ptype
!	character(6), dimension(:), allocatable    :: stype
!	character(6), dimension(:), allocatable    :: ttype

!	character(1), dimension(:), allocatable    	:: mlevel
!	integer, dimension(:), allocatable    		:: index

	integer     MP,MS,MT,MM

	character(150), dimension(1:MHR)  :: headrec
	character*1 mtype
	character*12 hashcode
    character*150 line

    integer np,ns,nt,ip,is,it
    integer msip,mtip,mtis

! .reg parameters

    integer MP1,MS1,MT1
    integer msip1,mtip1,mtis1,nr,nreg,regt
    integer n1,n2,n3,nrv(MRV),irv(MRV)
    real(kind=8) :: rv(MRV)
    real (kind=8), dimension(:), allocatable    :: regd
    character*12 hashcode1
    character*400 rline,txt
    character*100 skip,fmt
    logical new,match
    
! working variables

    real(kind=8) :: abox(4),gbox(4)
    real(kind=8) :: lat,lon
    integer      :: unm(1000)
    character*120 infile,outfile,regfile
    integer i,j,k,nk,nu
    logical anta,grnl

    data abox/-79.d0, -71.d0, 180.d0, 310.d0/
    data gbox/ 59.d0,  84.d0, 295.d0, 350.d0/
    

1   format(i7,2x,a6,2i8,2f10.4,f11.1,f18.0,2f9.1,f6.0,f6.1,i6,2x,a20)
2   format(i7,2x,a6,i8,8x,2f10.4,f11.1,f18.0,2f9.1,f6.0,2i6,2x,a20)
3   format(i7,2x,a6,2f10.4,f11.1,f15.0,2f9.1,f6.0,4i6,4x,a20)

4    format(a150)
5    format(9x,a1)
7    format(9x,i4)
8    format(a,$)
9    format(i8,$)

! parse command line

    call getarg(1,infile)
    if(infile.eq.' ') then
		write(*,*)
        write(*,*) "usage: scanreg <mascon_file>"
		write(*,*)
		stop
        endif

    regfile=trim(infile)//".reg"
    outfile=trim(infile)//"_scan"

    write(*,*)
    write(*,*) "input mascon file ",infile
    write(*,*) "input .reg file   ",regfile



! setup

    np=0
    nr=9
    ip=0
    do j=1,nr
       nrv(j)=0
    enddo
    

! open mascon files

    open(10,file=infile,action="read",status="old")
    open(11,file=outfile,status='unknown')
    open(12,file=regfile,action="read",status="old")

! Allocate arrays

	read(10,*) hashcode,MP,MS,MT,msip,mtip,mtis
    write(*,*) 
    write(*,*) "mascon file:" 
    write(*,'(a,3i9,3i5)') hashcode,MP,MS,MT,msip,mtip,mtis
	rewind(10)

	allocate(pmlat(MP),pmlon(MP),pmr(MP),pma(MP),pmh(MP),pmg(MP),pmden(MP),pcland(MP))
	allocate(nsip(MP),ntip(MP),tidemask(MP))
	allocate(pregion(MP))
	allocate(ptype(MP))
	allocate(regd(MP))

! skip over header and read mascon file

200 read(10,4,end=920) line
    if(line(1:1).eq."#") then
!        write(11,"(a)") trim(line)
        go to 200
        endif
    backspace(10)

    write(*,8) "reading msc: "
210 read(10,4,end=220) line
    j=j+1
    read(line,5) mtype
    if(mtype=="P")then

        read(line,1) i,ptype(i),nsip(i),ntip(i),pmlat(i),pmlon(i),pmr(i),pma(i),pmh(i),pmg(i),pmden(i),&
					 pcland(i),tidemask(i),pregion(i)

        ip=ip+1
        np=max(i,np)
        endif
    if(mod(j,200000)==0) write(*,9) j
    go to 210

220 close(10)
    write(*,*)
    write(*,*)

! read mascon.reg file and load diagonal into regd

    read(12,*) hashcode1,MP1,MS1,MT1,msip1,mtip1,mtis1
    if(hashcode1.ne.hashcode) write(*,*) "hashcode mismatch: ",hashcode,hashcode1
    if(MP1.ne.MP) write(*,*) "MP mismatch: ",MP,MP1
    if(MS1.ne.MS) write(*,*) "MS mismatch: ",MS,MS1
    if(MT1.ne.MT) write(*,*) "MT mismatch: ",MT,MT1
    if(msip1.ne.msip) write(*,*) "msip mismatch: ",msip,msip1
    if(mtip1.ne.mtip) write(*,*) "mtip mismatch: ",mtip,mtip1
    if(mtis1.ne.mtis) write(*,*) "mtis mismatch: ",mtis,mtis1

! parse regularisation description
    
    read(12,'(a)') rline
    do i=3,MRC
       if(rline(i-2:i).eq."km.") then
          n1=i
          exit
       endif
    enddo
    nr=1
    do i=n1+1,MRC
       if(rline(i:i).eq."/") nr=nr+1
       if(rline(i:i).eq.":") then
          n2=i
          exit
       endif
    enddo
    do i=n2+2,MRC
       if(rline(i-1:i).eq."m.") then
          n3=i
          exit
       endif
    enddo

    write(*,*) "regularisation file:"
    write(*,*) rline(1:n1)
    write(*,*) rline(n1+1:n2)

!    write(*,8) "reading reg:"
    read(12,*) regd(1)
    do ip=2,ip
       write(skip,*) (ip-1)*26
       fmt="("//trim(adjustl(skip))//"x,f25.16)"
       read(12,trim(fmt)) regd(ip)
!       if(mod(ip,1000)==0) write(*,9) ip
    enddo

    write(*,*)
    write(*,*) 
    write(*,'(i5,a)') nr," regularisation regions (Units: cm)"
    txt=rline(n2+1:n3)
!    write(*,*) txt
    read(txt,*) (rv(i),i=1,nr)
    do j=1,nr
       rv(j)=rv(j)*100
    enddo
    write(*,'(4x,12f8.0)') (rv(j),j=1,nr)

! count regularisations

    nu=0
    do ip=1,np
       regd(ip)=100.d0*regd(ip)**(-0.5d0)
       match=.false.
       do j=1,nr
          if(abs(regd(ip)-rv(j)).lt.0.001) then
             nrv(j)=nrv(j)+1
             match=.true.
             exit
          endif
       enddo
       if(.not.match) then
          nu=nu+1
          unm(nu)=ip
          endif
    enddo

    write(*,'(3x,12i8)') (nrv(j),j=1,nr)

    write(*,*)
    write(*,*) "   reg(cm)   #primaries"
    write(*,*) "   -------   ----------"
    new=.true.
    nk=1
    irv(1)=nint(rv(1))
    write(*,'(f10.1,i11)') rv(1),nrv(1)
    nreg=nrv(1)
    do j=2,nr
       i=nint(rv(j))
       do k=1,nk
          if(i.eq.irv(k)) new=.false.
       enddo
       if(new) then
          nk=nk+1
          irv(nk)=i
          write(*,'(f10.1,i11)') rv(j),nrv(j)
          nreg=nreg+nrv(j)
       endif
       enddo

    if(nu.gt.0) then
       write(*,*)
       write(*,'(i4,a)') nu," unmatched regularizations found:"
       write(*,*) "   ip     reg(cm)   region"
       do i=1,nu
          ip=unm(i)
          write(*,'(i6,2x,f10.4,1x,a)') ip,regd(ip),pregion(ip)
       enddo
    endif       
    
    write(*,*)
    write(*,*) "mascon count: ",np
    write(*,*) "reg    count: ",nreg
    if(nreg.ne.np) write(*,'(a,i5)') " count mismatch:",np-nreg
    
    
! output rogue mascons with high regularisation level
    
    ns=0
    regt=50
    do i=1,np
       if(nint(regd(i)).eq.regt) then
          lat=pmlat(i)
          lon=pmlon(i)
          
! is it Antarctic?
          anta=(lat.gt.abox(1).and.lat.lt.abox(2).and.lon.gt.abox(3).and.lon.lt.abox(4))

! is it Greenland?
          grnl=(lat.gt.gbox(1).and.lat.lt.gbox(2).and.lon.gt.gbox(3).and.lon.lt.gbox(4))

          if(.not.(anta.or.grnl)) then
!          if(.not.anta) then
             ns=ns+1
             write(11,1)  i,ptype(i),nsip(i),ntip(i),pmlat(i),pmlon(i),pmr(i),&
                  pma(i),pmh(i),pmg(i),pmden(i),pcland(i),tidemask(i),pregion(i)
          endif
       endif
    enddo

    close(11)
    write(*,*)
    write(*,'(i5,a,i3,a)') ns," mascons with reg=",regt,"cm outside Antarctic and Greenland boxes"
!    write(*,'(i5,a,i3,a)') ns," mascons with reg=",regt,"cm outside Antarctic box"
    write(*,'(8x,a15,4f6.1)') "Antarctic box: ",(abox(i),i=1,4)
    write(*,'(8x,a15,4f6.1)') "Greenland box: ",(gbox(i),i=1,4)
    write(*,*) "  List output to: ",outfile

    stop

920 stop "mascon file error: no mascons?"
    end

!************************************************************
