      program colour_mascons
!     *******
!
!  Reads mascon file and sets mascon colours in a gmt plot file,
!   either from default colour codes or from optional fit file
!
!  input mascon file from command line argument (default: mascons_stage2)
!  output  msc_to_plot.dat
!       or msc-to_plot.txt for option tt
!
!  usage: colour_mascons -{pp|sf|sp|sh|sg|tf|tp|ts|th|tg|tt} latmin lonmin latmax lonmax
!                                  input_mascon_file [input_fit_file|N] [vcv|N] [apr|N] ]
!
!  If a fit file is supplied, the colours are fit values scaled by the multiplier indicated 
!   (currently 1000), otherwise they are index values taken from the mascon file
!    - output colour selcted by comand line option (eg -tp) as follows:
!       -p* : ip             ; primary ID#
!       -sf : pcol(sap(i))   ; fit value of secondary's associated primary
!       -sp : spcol          ; index colour of enclosing primary
!       -sh : smh            ; height of secondary
!       -sg : smg            ; geoid height of secondary
!       -tf : pcol(tap(i))   ; fit value of ternary's associated primary
!       -tp : tpcol          ; index colour of enclosing primary
!       -ts : tscol          ; index colour of enclosing secondary
!       -th : tmh            ; height of ternary
!       -tg : tmg            ; geoid height of ternary
!       -tt : tap(i)         ; creates a plot file of primary numbers on ternary points
!
! MODS
! PT181203: modified to read 5-digit mascon numbers

    implicit none
    real (kind=8), parameter    :: pi = 3.141592653589793d0

    real (kind=8), dimension(:), allocatable    :: pmlat,pmlon,pmr,pma,pmh,pmg,pmden,pcland,val
    real (kind=8), dimension(:), allocatable    :: smlat,smlon,smr,sma,smh,smg,smden
    real (kind=8), dimension(:), allocatable    :: tmlat,tmlon,tmr,tma,tmh,tmg,tmden,flag
    integer, dimension(:), allocatable          :: pmn,nsip,ntip,tidemask,pcol
    integer, dimension(:), allocatable          :: smn,ntis,sap,spcol
    integer, dimension(:), allocatable          :: tmn,tap,tas,tpcol,tscol

    character*12, dimension(:), allocatable	    :: pregion
    character*12, dimension(:), allocatable	    :: sregion
    character*12, dimension(:), allocatable	    :: tregion
    character(6), dimension(:), allocatable     :: ptype
    character(6), dimension(:), allocatable     :: stype
    character(6), dimension(:), allocatable     :: ttype


    real (kind=8) dum,rd,valmin,valmax,scale
    real (kind=8) latmin,lonmin,latmax,lonmax,mlat,mlon,dlat,dlon
    integer MP,MS,MT,maxsip,maxtip,maxtis
    integer np,ns,nt,ip,is,it
    integer nsipt,ntipt,ntist,nm
    integer nsea,nfr,ndry,nhvy,nderr,nhigh,nlow
    integer fsize,fangle,fnum,colmin,colmax
    integer i,j,ncol,col
    logical zerox

    character*1 mtype,opt1,opt2
    character*3 colopt
    character*4 fposition
    character*5 arg
    character*12 hashcode
    character*20 layout
    character*120 infile,fitfile
    character*180 line

! PT180226: allow for the reading of vcv files instead of fit files
    logical :: vcv

! PT180322: allow for the a priori mascon values to be plotted instead
    logical :: apr

! PT181203: determine the version of the input file
    integer*4 :: fit_version

1   format(i7,2x,a6,2i8,2f10.4,f11.1,f18.0,2f9.1,f6.0,f6.1,i6,2x,a20)
2   format(i7,2x,a6,i8,8x,2f10.4,f11.1,f18.0,2f9.1,f6.0,2i6,2x,a20)
3   format(i7,2x,a6,2f10.4,f11.1,f15.0,2f9.1,f6.0,4i6,12x,a6)

4    format(a180)
5    format(9x,a1)
6    format(2f9.4,i8,a7)
7    format(9x,i5)      ! PT171122: increased from 4 to 5 for 41256 mascons
8    format(2f9.4,a20,i8)

9    format(i8,$)

! parse command line

	if(iargc().lt.5) then
		write(*,*) " usage: colour_mascons -{pp|sf|sp|sh|sg|tf|tp|ts|th|tg|tt} ",&
					"latmin lonmin latmax lonmax input_mascon_file [input_fit_file|N] [vcv|N] [apr|N] "
		stop
		endif

    call getarg(1,colopt)
!    write(*,*) colopt
    if(colopt(1:1).ne."-") then
		stop "colouring option format error"
		endif
	opt1=colopt(2:2)
	opt2=colopt(3:3)

    call getarg(2,arg)
    read(arg,*) latmin
    call getarg(3,arg)
    read(arg,*) lonmin
    call getarg(4,arg)
    read(arg,*) latmax
    call getarg(5,arg)
    read(arg,*) lonmax
    if(latmin.gt.latmax) stop "!!  range error latmin > latmax  - aborting"
    if(lonmin.gt.lonmax) stop "!!  range error lonmin > lonmax  - aborting"
    if(lonmin.lt.0) then
        if(lonmax.lt.0) then
            lonmin=lonmin+360
            lonmax=lonmax+360
        else
            zerox=.true.
            endif 
        endif
    mlat=(latmin+latmax)/2
    dlat=latmax-mlat
    mlon=(lonmin+lonmax)/2
    dlon=lonmax-mlon
    write(*,'((a,2f10.3))') " plot midpoint:",mlat,mlon,"     range+/-:",dlat,dlon
!    write(*,"(4f6.1)") latmin,lonmin,latmax,lonmax
!    write(*,"(4f6.1)") mlat,dlat,mlon,dlon

! set default input values

    call getarg(6,infile)
    if(infile.eq.' ') then
        infile="mascons_stage2"
        write(*,*) "Note - mascon file not specified, using default"
        endif
    write(*,*) "  input mascon file: ",infile
    open(10,file=infile,status="old")

    call getarg(7,fitfile)
    if(fitfile.eq."N") then
        write(*,*) "No fit file specified - native mascon colours applied "
    else
        write(*,*) "input fit file: ",fitfile
        opt2="f"
        endif

! PT180226: check to see whether we want to read a vcv file (defaults is a fit file)
    call getarg(8,arg)
    if (arg(1:3) == "vcv" .or. arg(1:3) == "VCV")then
        print*,'read values from an input VCV file'
        vcv = .true.
    else
        vcv = .false.
    endif

! PT180322: check to see whether we want to read the apriori mascon values (defaults is the estimated values)
    call getarg(9,arg)
    if (arg(1:3) == "apr" .or. arg(1:3) == "APR")then
        print*,'read a priori mascon values from the input file'
        apr = .true.
    else
        apr = .false.
    endif

! initialize parameters

    rd=pi/180.d0
    np=0
    ns=0
    nt=0
!    MP=4582
!    MS=4582
!    MT=1485118
    nsipt=0
    ntipt=0
    ntist=0
    i=0
    j=0
    ip=0
    is=0
    it=0
    ncol=0
    valmin=100.
    valmax=-100.
    colmin=0
    colmax=0
    nhigh=0
    nlow=0
    fsize=8
    fangle=0
    fnum=4
    fposition=" CM "
    write(layout,'(3(i4,1x),a4)') fsize,fangle,fnum,fposition
    if(opt2=="t") then
        if(max(dlat,dlon).gt.3) write(*,*) "! range >3deg 0n a number plot - this may get messy"
        if(max(dlat,dlon).gt.10) write(*,*) "! range >10deg on a number plot - reconsider dimensions"
        endif

! Allocate arrays

    read(10,*) hashcode,MP,MS,MT,maxsip,maxtip,maxtis

    allocate(pmlat(MP),pmlon(MP),pmr(MP),pma(MP),pmh(MP),pmg(MP),pmden(MP),pcland(MP),val(MP))
    allocate(smlat(MS),smlon(MS),smr(MS),sma(MS),smh(MS),smg(MS),smden(MS))
    allocate(tmlat(MT),tmlon(MT),tmr(MT),tma(MT),tmh(MT),tmg(MT),tmden(MT))

    allocate(pmn(MP),nsip(MP),ntip(MP),tidemask(MP),pcol(MP))
    allocate(smn(MS),ntis(MS),sap(MS),spcol(MS))
    allocate(tmn(MT),tap(MT),tas(MT),tpcol(MT),tscol(MT))

    allocate(pregion(MP))
    allocate(sregion(MS))
    allocate(tregion(MT))
    allocate(ptype(MP))
    allocate(stype(MS))
    allocate(ttype(MT))
    allocate(flag(MT))

! read fit file data

    if(opt2.eq."f") then
        open(12,file=fitfile,status="old")
! PT181203: determine whether it is a V2 or a V3 file
        read(12,'(a2)')line
        if(line(1:2) == "V2")then
          print*,"V2 solution file"
          fit_version = 2
        else if (line(1:2) == "V3")then
          print*,"V3 solution file"
          fit_version = 3
        endif

110     read(12,4,end=120) line
          if(line(8:9).eq."MC") then
            read(line,7) i
!print*,'mascon: ',i
! PT180226: read either a fit or a vcv file
            if(.not.vcv)then
                if(.not. apr)read(line(70:80),*) val(i)
                if(      apr)read(line(38:47),*) val(i)
            else
                if(.not. apr)read(line(50:66),*) val(i)
                if(      apr)read(line(38:47),*) val(i)
                endif
            if (val(i).gt.1.0) nhigh=nhigh+1
            if (val(i).lt.-1.0) nlow=nlow+1
            valmin=min(val(i),valmin)
            valmax=max(val(i),valmax)
            ncol=ncol+1
!           if(ncol.le.10) write(*,"(i6,f10.4)") i,val(i)
            endif
        go to 110
120     close(12)

        write(*,*) ncol," primary mascon fits found"
        scale=1000
        write(*,"(a14,2f8.4,a8,f8.1)") "value range:",valmin,valmax,"  multiplier:",scale
        write(*,*) nlow," points below -1m,   ",nhigh,"points over +1m"
        do i=1,ncol
            pcol(i)=val(i)*scale
            enddo
        endif

! open output mascon file

    open(11,file='msc_to_plot.dat',status='unknown')
    if(opt2=="t") open(13,file='msc_to_plot.txt',status='unknown')
    write(*,*) "reading mascon file..."

! skip over header

150 read(10,4,end=920) line
    if(line(1:1).eq."#") then
!        write(11,"(a)") trim(line)
        go to 150
        endif
    backspace(10)


!!! ploting primaries !!!

if(opt1=="p") then
200 read(10,4,end=210) line
    j=j+1
    read(line,5) mtype
    if(mtype=="P")then
        read(line,*) i,ptype(i),nsip(i),ntip(i),pmlat(i),pmlon(i),pmr(i),pma(i),pmh(i),pmg(i),pmden(i),&
					pcland(i),tidemask(i),pregion(i)
!        write(11,1)  i,ptype(i),nsip(i),ntip(i),pmlat(i),pmlon(i),pmr(i),pma(i),pmh(i),pmg(i),pmden(i),&
!					pcland(i),tidemask(i),pregion(i)
        ip=ip+1
        pmn(ip)=i
        np=max(i,np)
        endif
    go to 200

210 write(*,*) "writing colour file"
    if(np.gt.ip) write(*,*) "!colour_mascons: ip<np - some primaries missing (partial mascon file?)",ip,np
220 do i=1,np
        if(abs(pmlat(i)-mlat).le.dlat.and.abs(pmlon(i)-mlon).le.dlon) then
            colmin=min(col,colmin)
            colmax=max(col,colmax)
            write(11,6) pmlat(i),pmlon(i),i
            endif
        enddo
! repeat loop if target range crosses zero meridian
    if(zerox) then
        write(*,*)
        write(*,*) "zero crossing: second pass for negative longitudes"
        zerox=.false.
        mlon=mlon+360
        go to 220
        endif



!!! ploting secondaries !!!

    else if(opt1=="s") then
300 read(10,4,end=310) line
    j=j+1
    read(line,5) mtype
    if(mtype=="S") then
        read(line,*) i,stype(i),ntis(i),smlat(i),smlon(i),smr(i),sma(i),smh(i),smg(i),smden(i),sap(i),spcol(i),sregion(i)
!        write(11,2)  i,stype(i),ntis(i),smlat(i),smlon(i),smr(i),sma(i),smh(i),smg(i),smden(i),sap(i),spcol(i),sregion(i)
        is=is+1
        smn(is)=i
        ns=max(i,ns)
        endif
    go to 300
310 write(*,*) "writing colour file"
    if(ns.gt.is) write(*,*) "!colour_mascons: is<ns - some secondaries missing (partial mascon file?)",is,ns
320 do i=1,ns
        if(abs(smlat(i)-mlat).le.dlat.and.abs(smlon(i)-mlon).le.dlon) then
            if(opt2=="f") then
                col=pcol(sap(i))
            else if(opt2=="p") then
                col=spcol(i)
            else if(opt2=="h") then
                col=smh(i)
            else if(opt2=="g") then
                col=smg(i)
            else
                col=i
                endif

            colmin=min(col,colmin)
            colmax=max(col,colmax)
            write(11,6)  smlat(i),smlon(i),col
            endif
        enddo
! repeat loop if target range crosses zero meridian
    if(zerox) then
        write(*,*)
        write(*,*) "zero crossing: second pass for negative longitudes"
        zerox=.false.
        mlon=mlon+360
        go to 320
        endif


!!! ploting ternaries !!!
    else if(opt1=="t") then
    flag = 0.d0
400 read(10,4,end=410) line
    j=j+1
    read(line,5) mtype
    if(mtype=="T") then
        read(line,*) i,ttype(i),tmlat(i),tmlon(i),tmr(i),tma(i),tmh(i),tmg(i),tmden(i),tap(i),tas(i),&
					tpcol(i),tscol(i),tregion(i)
        if(ttype(i)=="TGIA")flag(i)=1.d0
        it=it+1
        tmn(it)=i
        nt=max(i,nt)
        endif
    if(mod(j,100000)==0) write(*,9) j
    go to 400
410 nm=j
    write(*,*)
    if(nt.gt.it) write(*,*) "!colour_mascons: it<nt - some ternaries missing (partial mascon file?)",it,nt
    write(*,*) "writing colour file"
    j=1
420 do i=1,nt
        if(abs(tmlat(i)-mlat).le.dlat.and.abs(tmlon(i)-mlon).le.dlon) then
            if(opt2=="f") then
                col=pcol(tap(i))
            else if(opt2=="p") then
                col=tpcol(i)
            else if(opt2=="s") then
                col=tscol(i)
            else if(opt2=="h") then
                col=tmh(i)
            else if(opt2=="g") then
                col=tmg(i)
            else if(opt2=="t") then
                col=tmh(i)/5
                write(13,8)  tmlat(i),tmlon(i),layout,tap(i)
            else
                col=tmh(i)/5
                endif

            colmin=min(col,colmin)
            colmax=max(col,colmax)
            write(11,6)  tmlat(i),tmlon(i),col,ttype(i)   !,i   ! PT180320: added "i" as debug ......

            j=j+1
            endif

        if(mod(j,100000)==0) write(*,9) j
        enddo
    write(*,*)

! repeat loop if target range crosses zero meridian
    if(zerox) then
        write(*,*)
        write(*,*) "zero crossing: second pass for negative longitudes"
        zerox=.false.
        mlon=mlon+360
        go to 420
        endif

    else
        stop "plot option error - opt1 not recognised"
        endif

    write(*,'(a32,i6,a3,i8)') "colour range of plotted mascons:",colmin," - ",colmax

500 close(10)
    close(11)
    if(opt2=="t") close(13)

    write(*,*) "mascon count (P,S,T): ",ip,is,it
    write(*,*) "mascon range (P,S,T): ",np,ns,nt
    write(*,*) "total mascon lines found   :",nm
!    if(np.ne.ip) write(*,*) "!! primary count differs from maximum index by  ",np-ip
!    if(ns.ne.is) write(*,*) "!! secondary count differs from maximum index by",ns-is
!    if(nt.ne.it) write(*,*) "!! ternary count differs from maximum index by  ",nt-it

!    call status_update('STATUS','UTIL','rewrite_mascons',"mascon_rw","End program",0)

    stop
910 stop "fit file error"
920 stop "mascon file error"
    end

!************************************************************

subroutine read_header(lun,headrec,mhr)

!
!   reads header records and verifies hashcode
!

    implicit none

    integer lun,i,nhr,mhr
    character*12 hashcode,checkcode
    character(150), dimension(1:mhr)     :: headrec


    i=0
    rewind(lun)
    read(lun,*) hashcode
    if(hashcode(1:1).ne."#") stop "error: '#' not found on header line 1"

100 i=i+1
    if(i.gt.mhr) go to 910
    read(lun,'(a)',end=900) headrec(i)
    if(headrec(i)(1:1).eq."#") go to 100

    nhr=i-1

    backspace(lun)
    return

900 stop "error: read_header - end of file"
910 stop "error: read_header - header too long, increase MHR"

end


