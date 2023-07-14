
      program audit_mascon_file
!     *******
!
!  Analyses mascon files producing structure summary and searching for errors
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
! dspmax,dtpmax,dtsmax - maximum secondary-primary distance, etc.
! dspgmx,dtpgmx,dtsgmx - global maximum secondary-primary distance, etc.
! dspgav,dtpgav,dtsgav - global average secondary-primary distance, etc.
! lsip,ltip,ltis - list of secondaries in primaries, etc.
! isip,itp,itis -
! hpind - index file of secondaries in primaries
! hsind - index file of ternaries in secondaries
! kips(ip) - start of primary ip in hpind
! kist(is) - start of secondary is in hsind
! kpsc(ip) - count of secondaries in primary ip
! ktsc(is) - count of ternaries in secondary is
!
! Files:
        ! 10 - infile: mascon file to be audited
        ! 11 - bands_aud: P/S/T latitude bands (for regular mascons)
        ! 12 - err_aud: errors of density, %land, # of ternaries
        ! 13 - plstat_aud: primary mascon statistics
        ! 14 - index-aud: indexing errors
        ! 15 - outliers_aud: ternaries >8 deg from their pimary
        ! 16 - mismatch_aud: P/S parameters incompatible with ternaries enclosed
        ! 17 - list of small primaries
!
    implicit none

!-------

    real (kind=8), parameter      :: pi = 3.141592653589793d0
    integer, parameter      :: MHR = 500

    real (kind=8), dimension(:), allocatable    :: pmlat,pmlon,pmr,pma,pmh,pmg,pmden,pcland
    real (kind=8), dimension(:), allocatable    :: smlat,smlon,smr,sma,smh,smg,smden
    real (kind=8), dimension(:), allocatable    :: tmlat,tmlon,tmr,tma,tmh,tmg,tmden
    real (kind=8), dimension(:), allocatable    :: dspmax,dtpmax,dtsmax,dppmin
    integer, dimension(:), allocatable          :: nsip,ntip,tidemask
    integer, dimension(:), allocatable          :: ntis,sap,sscol
    integer, dimension(:), allocatable          :: tap,tas,tpcol,tscol
    integer, dimension(:), allocatable          :: isip,itip,kips,kpsc
    integer, dimension(:), allocatable          :: itis,kist,kstc
    integer, dimension(:), allocatable          :: hsind,hpind
    integer, dimension(:), allocatable          :: scheck,tcheck

    character*12, dimension(:), allocatable	   :: pregion
    character*12, dimension(:), allocatable	   :: sregion
    character*12, dimension(:), allocatable	   :: tregion
    character(6), dimension(:), allocatable    :: ptype
    character(6), dimension(:), allocatable    :: stype
    character(6), dimension(:), allocatable    :: ttype

    character(150), dimension(1:MHR)  :: headrec

!--------
    
    integer, dimension(:,:), allocatable :: lsip,ltip,ltis

    real (kind=8) dum,rad,WGS84_area
    real (kind=8) area_p_total,area_s_total,area_t_total,oceanp,afr,asea
    real (kind=8) dsp,dtp,dts,dpp,splat,cplat,sslat,cslat
    real (kind=8) dspgmx,dtpgmx,dtsgmx,dspgav,dtpgav,dtsgav
    real (kind=8) dppgmn,dppgav
    real (kind=8) xlat,xlon,dx,dxmin,s1,s2,c1,c2,dlonr,dx0
    real (kind=8) hmin,hmax,hav,gmin,gmax,gav,rmin,rmax,rav
    integer     MP,MS,MT
    character*1 mtype,choice
    character*12 hashcode

    integer, dimension(1:12)         :: pclhist
    integer np,ns,nt
    integer npex,nsex,ntex,maxsip,maxtip,maxtis
    integer nsipt,ntipt,ntist
    integer mpi,msi,mti
    integer nsea,nfr,ndry,nhvy,nderr,sindx_errs,pindx_errs,ndeft

! close and small mascon variables
    real (kind=8), dimension(:), allocatable  :: prox1d,prox2d
    integer, dimension(:,:), allocatable :: prox1,prox2
    integer, dimension(:), allocatable :: small
    real    tol1,tol2
    integer ndpp1,ndpp2,nsmall,ncrit,i1,i2,MP2

    character*180 line
    character*120 infile
    integer i,j,lun,indx,iden,ip,is,it,ic
    logical flag,sflag,tflag

! histogram
    integer hv(10),hi(10),nhi,ntipmax

! mismatch test variables

    real (kind=8) sum1,sum2,sum3,sum4,sum5,sum6
    real (kind=8) x,y,z,h,a1,a2,dif,difmin,al,as
    real (kind=8) sumx,sumy,sumz,amlat,amlon
    real (kind=8) cmlat,cmlon,cmr,cma,cmh,cmg,cmden,cpcland,tol
    integer is1,it1,k,mispt,misps,mismatch
    logical tpcol_flag,tscol_flag
    


1   format(i7,2x,a6,2i8,2f10.4,f11.1,f18.0,2f9.1,f6.0,f6.1,i6,2x,a20)
2   format(i7,2x,a6,i8,8x,2f10.4,f11.1,f18.0,2f9.1,f6.0,2i6,2x,a20)
3   format(i7,2x,a6,2f10.4,f11.1,f15.0,2f9.1,f6.0,4i6,12x,a12)
4   format(9x,a1)
5   format(6x,"Land: ",i8,6x,"Sea : ",i8,6x,"%Ocean: ",f8.2)
6   format(12i7)
7   format(i6,2f9.3,i7,f9.3,i7,2f9.3)
8   format("Max "a20," range: Global Avg:",f9.3,5x,"Global Max:",f9.3," deg")
9   format(i8,$)
10  format(a8,f9.5)
11  format(a25,3f15.3)
12  format(a22,f20.2,f10.6," % diff wrt WGS84")    

    call getarg(1,infile)
    if(infile.eq.' ') then
        infile="mascons_stage2"
        write(*,*) "Note - mascon file not specified, assuming default: ",infile
        endif
    write(*,*) "opening '",infile,"'"


! setup

    open(12,file="err_aud",status="unknown")
    open(14,file="index_aud",status="unknown")

    rad=pi/180.d0
    MP=6000
    MS=6000
    MT=1500000
    np=0
    ns=0
    nt=0
    mpi=0
    msi=0
    mti=0
    nsipt=0
    ntipt=0
    ntist=0
    ndpp1=0
    ndpp2=0
    nsmall=0
    area_p_total=0.d0
    area_s_total=0.d0
    area_t_total=0.d0
    i=0
    j=0
    sindx_errs=0
    pindx_errs=0
    do i=1,12
        pclhist(i)=0
        enddo


! Input

    open(10,file=infile,status='old')
    write(*,*) "reading..."
    read(10,*) hashcode,npex,nsex,ntex,maxsip,maxtip,maxtis

    lun=10
    call read_header(lun,headrec,MHR)
    if(hashcode(1:1).ne."#") stop "error: '#' not found on header line 1"
    if(headrec(1)(1:1).ne."#") stop "error: '#' not found on header line 2"

! check that parameters are adequate for mascon file

    if(npex.gt.MP) then
        write(*,*) "expanding max Primaries to", npex
        MP=npex
        endif
    if(nsex.gt.MS) then
        write(*,*) "expanding max Secondaries to", nsex
        MS=nsex
        endif
    if(ntex.gt.MT) then
        write(*,*) "expanding max Ternaries to", ntex
        MT=ntex
        endif

    write(*,*) "File reports MaxSIP,MaxTIP,MaxTIS: ",maxsip,maxtip,maxtis

! Allocate arrays

    allocate(pmlat(MP),pmlon(MP),pmr(MP),pma(MP),pmh(MP),pmg(MP),pmden(MP),pcland(MP))
    allocate(smlat(MS),smlon(MS),smr(MS),sma(MS),smh(MS),smg(MS),smden(MS))
    allocate(tmlat(MT),tmlon(MT),tmr(MT),tma(MT),tmh(MT),tmg(MT),tmden(MT))
    allocate(dspmax(MP),dtpmax(MP),dtsmax(MP),dppmin(MP))
    allocate(nsip(MP),ntip(MP),tidemask(MP))
    allocate(ntis(MS),sap(MS),sscol(MS))
    allocate(tap(MT),tas(MT),tpcol(MT),tscol(MT))
    allocate(isip(MP),itip(MP),kips(MP),kpsc(MP))
    allocate(itis(MS),kist(MS),kstc(MS))
    allocate(hpind(MS),scheck(MS))
    allocate(hsind(MT),tcheck(MT))

    allocate(pregion(MP))
    allocate(sregion(MS))
    allocate(tregion(MT))
    allocate(ptype(MP))
    allocate(stype(MS))
    allocate(ttype(MT))

    allocate(lsip(npex,maxsip))
    allocate(ltip(npex,maxtip))
    allocate(ltis(nsex,maxtis))

    MP2=(MP*MP+MP)/2
    allocate(prox1(MP2,2))
    allocate(prox2(MP2,2))
    allocate(prox1d(MP2))
    allocate(prox2d(MP2))
    allocate(small(MP2))

    
! read mascon file

 110  read(10,'(a180)',end=120) line
     read(line,4) mtype
    j=j+1
    if(mtype.eq."P")then
        read(line,*) i,ptype(i),nsip(i),ntip(i),pmlat(i),pmlon(i),pmr(i),pma(i),pmh(i),pmg(i),pmden(i),&
			pcland(i),tidemask(i),pregion(i)
        np=np+1
 !       write(99,*) np,j
        mpi=max(i,mpi)
        nsipt=nsipt+nsip(i)
        ntipt=ntipt+ntip(i)
        area_p_total=area_p_total+pma(i)

    else if(mtype.eq."S") then
        read(line,*) i,stype(i),ntis(i),smlat(i),smlon(i),smr(i),sma(i),smh(i),smg(i),smden(i),sap(i),sscol(i),sregion(i)
 		ns=ns+1
        msi=max(i,msi)
        ntist=ntist+ntis(i)
        area_s_total=area_s_total+sma(i)

    else if(mtype.eq."T") then
       read(line,*) i,ttype(i),tmlat(i),tmlon(i),tmr(i),tma(i),tmh(i),tmg(i),tmden(i),tap(i),tas(i),tpcol(i),tscol(i),tregion(i)
       ! repair tas connection lost during reshape
       tas(i)=tap(i)
       nt=nt+1
       mti=max(i,mti)
       area_t_total=area_t_total+tma(i)

    else
        write(*,*) line
        stop "mascon type error"
        endif

    if(mod(j,100000).eq.0) write(*,9) j
    go to 110

 120  close(10)


! write out basic counts

    write(*,*)
    write(*,*)
    write(*,*) "mascons expected (P,S,T): ",npex,nsex,ntex
    write(*,*) "mascons found    (P,S,T): ",np,ns,nt
    if(np.ne.npex) write(*,*) "! Primary count mismatch"
    if(ns.ne.nsex) write(*,*) "! Secondary count mismatch"
    if(nt.ne.ntex) write(*,*) "! Ternary count mismatch"
    write(*,*) "    total mascons: ",np+ns+nt
    write(*,*) "      check total: ",j
    write(*,*) "max index (P,S,T):    ",mpi,msi,mti
    write(*,*) "max index sum:   ",mpi+msi+mti
    write(*,*) "count of S in P: ",nsipt
    write(*,*) "count of T in P: ",ntipt
    write(*,*) "count of T in S: ",ntist
	if(mpi.gt.np.or.msi.gt.ns.or.mti.gt.nt) write(*,*) "! Warning: mascon file is incomplete"

    WGS84_area=510065621724088.44

    write(*,12) "primary total area:   ",area_p_total,100*(area_p_total-WGS84_area)/WGS84_area
    write(*,12) "secondary total area: ",area_s_total,100*(area_s_total-WGS84_area)/WGS84_area
    write(*,12) "ternary total area:   ",area_t_total,100*(area_t_total-WGS84_area)/WGS84_area
    write(*,12) "WGS84 ellipsoid area: ",WGS84_area
    if(area_t_total.lt.0.999*WGS84_area) write (*,*) "! Warning: ternary area incomplete on ellipsoid"

! build heirarchical index files

! ternaries in secondary list
    kist(1)=1
    kstc(1)=0
    do is=2,ns
        kist(is)=kist(is-1)+ntis(is-1)
        kstc(is)=0
        enddo
    do it=1,nt
        is=tas(it)
        indx=kist(is)+kstc(is)
        kstc(is)=kstc(is)+1
        hsind(indx)=it
!	if(abs(it-2).le.1) write(*,*) it,is,kist(is),kstc(is),indx,hsind(indx)
!	if(abs(it-247).le.1) write(*,*) it,is,kist(is),kstc(is),indx,hsind(indx)
        enddo
    do is=1,ns
        if(kstc(is).ne.ntis(is)) then
            write(14,*) "indexing error on secondary: ",is,kstc(is),ntis(is)
            sindx_errs=sindx_errs+1
            endif
        enddo

! secondaries in primary list
    kips(1)=1
    kpsc(1)=0
    do ip=2,np
        kips(ip)=kips(ip-1)+nsip(ip)
        kpsc(ip)=0
        enddo
    do is=1,ns
        ip=sap(is)
        if(ip.le.0) then
            write(*,*) "sap error: ",is,ip
            stop
            endif
        indx=kips(ip)+kpsc(ip)
        kpsc(ip)=kpsc(ip)+1
        hpind(indx)=is
        enddo
    do ip=1,np
        if(kpsc(ip).ne.nsip(ip)) then
            write(14,*) "indexing error on primary: ",ip,kpsc(ip),nsip(ip)
            pindx_errs=pindx_errs+1
            endif
        enddo

!  sum of tas(it) doesn't match ntis(is)
    if(sindx_errs.gt.0) write(*,*) "! ",sindx_errs," secondary indexing errors found, see index.aud"

!  sum of tap(it) doesn't match ntip(ip)
    if(pindx_errs.gt.0) write(*,*) "! ",pindx_errs," primary indexing errors found, see index.aud"

! build check arrays from the index files
	do is=1,ns
		scheck(is)=0
		enddo
	do it=1,nt
		tcheck(it)=0
		enddo
    do ip=1,np
        do i=1,nsip(ip)
            is=hpind(kips(ip)+i-1)
			scheck(is)=scheck(is)+1
            do j=1,ntis(is)
				indx=kist(is)+j-1
				it=hsind(indx)
				tcheck(it)=tcheck(it)+1
				enddo
            enddo
        enddo

! verify index file completeness
!  - log an error if any mascon appears more or less than once
	sflag=.false.
	do is=1,ns
		if(scheck(is).ne.1) then
			write(14,*) "indexing error: secondary# ",is,"  found ",scheck(is)," times"
			sflag=.true.
			endif
		enddo
	if(sflag) write(*,*) "! Warning: secondary mascon list error, see index.aud"
	tflag=.false.
	do it=1,nt
		if(tcheck(it).ne.1) then
			write(14,*) "indexing error: ternary# ",it,"  found ",tcheck(it)," times"
			tflag=.true.
			endif
		enddo
	if(tflag) write(*,*) "! Warning: ternary mascon list error, see index.aud"



! Primary mascon analysis

    nsea=0
    nfr=0
    ndry=0
    nhvy=0
    nderr=0
    afr=0.
    asea=0.
    flag=.false.
	hmin=0.d0
	hmax=0.d0
	hav=0.d0
	gmin=0.d0
	gmax=0.d0
	gav=0.d0
	rmin=1.d12
	rmax=0.d0
	rav=0.d0

    do i=1,np
        if(ptype(i).eq."      ") then
            write(*,*) "error: mascon type missing for primary #",i
            stop
            endif
        iden=nint(pmden(i))
        if(iden.eq.1000) then
            nfr=nfr+1
            afr=afr+pma(i)
        else if(iden.eq.1029) then
            nsea=nsea+1
            asea=asea+pma(i)
        else if(iden.eq.1030) then
            nhvy=nhvy+1
        else if(iden.eq.0) then
            ndry=ndry+1
        else
            write(12,*) " !density primary",i,pmden(i)
            nderr=nderr+1
            endif
        if(pcland(i).eq.0.d0) then
            indx=1
        else if(pcland(i).eq.100.d0) then
            indx=12
        else
            indx=2+pcland(i)/10
            endif
        if(indx.lt.1.or.indx.gt.12) then
            write(12,*) "pcland out of range on primary",i,pcland(i)
            flag=.true.
            endif
        pclhist(indx)=pclhist(indx)+1
		if(pmh(i).lt.hmin) hmin=pmh(i)
		if(pmh(i).gt.hmax) hmax=pmh(i)
		hav=hav+pmh(i)
		if(pmg(i).lt.gmin) gmin=pmg(i)
		if(pmg(i).gt.gmax) gmax=pmg(i)
		gav=gav+pmg(i)
		if(pmr(i).lt.rmin) rmin=pmr(i)
		if(pmr(i).gt.rmax) rmax=pmr(i)
		rav=rav+pmr(i)
        enddo

    if(flag) write(*,*) "! errors found in pcland entries"
    write(*,*)
    write(*,*) "primary mascons:"
    oceanp=100*asea/(asea+afr)
    write(*,5) nfr,nsea,oceanp
    if(ndry.gt.0) write(*,*) "   ! zero density: ",ndry
    if(nhvy.gt.0) write(*,*) "   ! high density: ",nhvy
    if(nderr.gt.0) then
        write(*,*) nderr," unidentified densities in primary mascons"
        endif
    write(*,*)
    write(*,*) "Histogram of percent land in primaries: 0, ten deciles, 100"
    write(*,6) (pclhist(i),i=1,12)
	write(*,11) "altitude min,max,mean: ",hmin,hmax,hav/np
	write(*,11) "geoid    min,max,mean: ",gmin,gmax,gav/np
	write(*,11) "radius   min,max,mean: ",rmin,rmax,rav/np

! Secondary mascon analysis

    if(ns==0) then
        write(*,*)
        write(*,*) "No secondary mascons found"
        go to 200
        endif
    nsea=0
    nfr=0
    ndry=0
    nhvy=0
    nderr=0
    afr=0.d0
    asea=0.d0
	hmin=0.d0
	hmax=0.d0
	hav=0.d0
	gmin=0.d0
	gmax=0.d0
	gav=0.d0
	rmin=1.d12
	rmax=0.d0
	rav=0.d0
    do i=1,ns
        if(stype(i).eq."      ") then
            write(*,*) "error: mascon type missing for secondary #",i
            stop
            endif

        ip=sap(i)
        isip(ip)=isip(ip)+1
        lsip(ip,isip(ip))=i

        iden=nint(smden(i))
        if(iden.eq.1000.) then
            nfr=nfr+1
            afr=afr+sma(i)
        else if(iden.eq.1029.) then
            nsea=nsea+1
            asea=asea+sma(i)
        else if(iden.eq.1030.) then
            nhvy=nhvy+1
        else if(iden.eq.0.) then
            ndry=ndry+1
        else
            write(12,*) " !density secondary",i,smden(i)
            nderr=nderr+1
            endif
		if(smh(i).lt.hmin) hmin=smh(i)
		if(smh(i).gt.hmax) hmax=smh(i)
		hav=hav+smh(i)
		if(smg(i).lt.gmin) gmin=smg(i)
		if(smg(i).gt.gmax) gmax=smg(i)
		gav=gav+smg(i)
		if(smr(i).lt.rmin) rmin=smr(i)
		if(smr(i).gt.rmax) rmax=smr(i)
		rav=rav+smr(i)

        enddo
    write(*,*)
    write(*,*) "secondary mascons:"
    oceanp=100*asea/(asea+afr)
    write(*,5) nfr,nsea,oceanp
    if(ndry.gt.0) write(*,*) "   ! zero density: ",ndry
    if(nhvy.gt.0) write(*,*) "   ! high density: ",nhvy
    if(nderr.gt.0) then
        write(*,*) nderr," unidentified densities in secondary mascons"
        endif
	write(*,11) "altitude min,max,mean: ",hmin,hmax,hav/ns
	write(*,11) "geoid    min,max,mean: ",gmin,gmax,gav/ns
	write(*,11) "radius   min,max,mean: ",rmin,rmax,rav/ns

! calculate secondary-primary inter-mascon ranges and maxima

    dspgmx=0
    dspgav=0
    do ip=1,np
        dspmax(ip)=0
        splat=sin(rad*pmlat(ip))
        cplat=cos(rad*pmlat(ip))
!        if(isip(ip).ne.nsip(ip)) write(12,*) "nsip check error ip,nsip,isip:",ip,nsip(ip),isip(ip)
        do i=1,nsip(ip)
!            is=lsip(ip,i)
            is=hpind(kips(ip)+i-1)
            dum=cplat*cos(rad*smlat(is))*cos(rad*(smlon(is)-pmlon(ip)))
            dum=dum+splat*sin(rad*smlat(is))
            dum=amin1(dum,1.)
            dsp=acos(dum)/rad
 !       write(12,*) ip,i,is,splat,cplat,rad,dsp
            dspmax(ip)=amax1(dsp,dspmax(ip))
            dspgmx=amax1(dsp,dspgmx)
            enddo
        dspgav=dspgav+dspmax(ip)
        enddo
    dspgav=dspgav/np


! Ternary mascon analysis

200    nsea=0
    nfr=0
    ndry=0
    nhvy=0
    nderr=0
	hmin=0.d0
	hmax=0.d0
	hav=0.d0
	gmin=0.d0
	gmax=0.d0
	gav=0.d0
	rmin=1.d12
	rmax=0.d0
	rav=0.d0
	ndeft=0

	do i=1,nt
		if(tmr(i).lt.6350.d3.or.tmden(i).lt.1000) ndeft=ndeft+1
		enddo
	if(ndeft.gt.0) then
		write(*,*)
		write(*,*) "! found ",ndeft," defective ternaries"
		endif

	do i=1,nt
		if(tap(i).gt.np.or.tap(i).lt.1) then
			write(*,*)
			write(*,*) "ternary associated witn nonexistent primary: it,tap,np",i,tap(i),np
			write(*,*) "  ... skipping to mascon inspection"
			write(*,*)
			go to 260
			endif
		if(tas(i).gt.ns.or.tas(i).lt.1) then
			write(*,*)
			write(*,*) "ternary associated witn nonexistent secondary: it,tas,ns",i,tas(i),ns
			write(*,*) "  ... skipping to mascon inspection"
			write(*,*)
			go to 260
			endif
		enddo

    do i=1,nt
        if(ttype(i).eq."      ") then
            write(*,*) "error: mascon type missing for ternary #",i
            stop
            endif

        ip=tap(i)
        itip(ip)=itip(ip)+1
        itip(ip)=min(itip(ip),maxtip)
        ltip(ip,itip(ip))=i

        is=tas(i)
!        if(is.eq.0) is=ip
        itis(is)=itis(is)+1
        itis(is)=min(itis(is),maxtis)
        ltis(is,itis(is))=i

        iden=nint(tmden(i))
        if(iden.eq.1000.) then
            nfr=nfr+1
        else if(iden.eq.1029.) then
            nsea=nsea+1
        else if(iden.eq.1030.) then
            nhvy=nhvy+1
        else if(iden.eq.0.) then
            ndry=ndry+1
        else
            write(12,*) " !density ternary",i,tmden(i)
            nderr=nderr+1
            endif

		if(tmh(i).lt.hmin) hmin=tmh(i)
		if(tmh(i).gt.hmax) hmax=tmh(i)
		hav=hav+tmh(i)
		if(tmg(i).lt.gmin) gmin=tmg(i)
		if(tmg(i).gt.gmax) gmax=tmg(i)
		gav=gav+tmg(i)
		if(tmr(i).lt.rmin) rmin=tmr(i)
		if(tmr(i).gt.rmax) rmax=tmr(i)
		rav=rav+tmr(i)

        enddo



    write(*,*)
    write(*,*) "ternary mascons:"

    oceanp=(nsea*100.)/nt
    write(*,5) nfr,nsea,oceanp
    if(ndry.gt.0) write(*,*) "   ! zero density: ",ndry
    if(nhvy.gt.0) write(*,*) "   ! high density: ",nhvy
    if(nderr.gt.0) then
        write(*,*) nderr," unidentified densities in ternary mascons"
        endif
	write(*,11) "altitude min,max,mean: ",hmin,hmax,hav/nt
	write(*,11) "geoid    min,max,mean: ",gmin,gmax,gav/nt
	write(*,11) "radius   min,max,mean: ",rmin,rmax,rav/nt

    write(*,*)

! calculate ternary-primary inter-mascon ranges and maxima

    dtpgmx=0
    dtpgav=0
    open(15,file="outliers_aud",status="unknown")
    write(15,*) "ternary outliers in primary mascons - range >",6
    write(15,*) "   ip     it   p-t_deg   pmlat    pmlon    tmlat    tmlon"
    
    do ip=1,np
        dtpmax(ip)=0
        splat=sin(rad*pmlat(ip))
        cplat=cos(rad*pmlat(ip))
        if(itip(ip).ne.ntip(ip)) write(12,*) "ntip check error ip,ntip,itip:",ip,ntip(ip),itip(ip)
        if(itip(ip).gt.maxtip) then
            write(*,'("! ",i8,a,i5,a)') itip(ip),"ternaries in primary",ip," -> not analysed"
            dtpmax(ip)=-1.
        else
        do i=1,ntip(ip)
            it=ltip(ip,i)
            dum=cplat*cos(rad*tmlat(it))*cos(rad*(tmlon(it)-pmlon(ip)))
            dum=dum+splat*sin(rad*tmlat(it))
            dum=amin1(dum,1.)
            dtp=acos(dum)/rad

            if(dtp.gt.6.) then
               write(15,'(i7,1x,i7,5f9.3)') ip,it,dtp,pmlat(ip),pmlon(ip),tmlat(it),tmlon(it)
             endif
            
            dtpmax(ip)=amax1(dtp,dtpmax(ip))
            dtpgmx=amax1(dtp,dtpgmx)
            enddo
            endif
        dtpgav=dtpgav+dtpmax(ip)
        enddo
     close(15)
     dtpgav=dtpgav/np
     write(*,*)

! calculate ternary-secondary inter-mascon ranges and maxima
    dtsgmx=0
    dtsgav=0
    do is=1,ns
!    do is=1,100
        dtsmax(is)=0
        sslat=sin(rad*smlat(is))
        cslat=cos(rad*smlat(is))
        if(kstc(is).le.0) then
            write(*,'(a,i6,a)') "! Warning: secondary#",is," contains no ternaries -> range analysis aborted"
            go to 250
            endif
        if(itis(is).ne.kstc(is)) write(12,*) "ntis check error is,ntis,itis,ktsc:",is,ntis(is),itis(is),kstc(is)

        if(itis(is).gt.maxtis) then
        write(*,'("! ",i8,a,i5,a)') itis(is)," ternaries in secondary#",is," -> not analysed"
            dtsmax(is)=-1.
        else

        do i=1,ntis(is)
!            it=hsind(kist(is)+i-1)
            it=ltis(is,i)
!            if(it.ne.ltis(is,i))then
!                write(12,*) "! mismatch ", is,i,ntis(is),it,ltis(is,i)
!                stop
!                endif
            dum=cslat*cos(rad*tmlat(it))*cos(rad*(tmlon(it)-smlon(is)))
            dum=dum+sslat*sin(rad*tmlat(it))
            dum=amin1(dum,1.)
            dts=acos(dum)/rad
            dtsmax(is)=amax1(dts,dtsmax(is))
            dtsgmx=amax1(dts,dtsgmx)
            enddo
            endif
        dtsgav=dtsgav+dtsmax(is)
        enddo
    dtsgav=dtsgav/ns
250 write(*,*)

! calculate primary-primary inter-mascon ranges and report close mascons
    
    dppgmn=40000  ! ~ earth circumference
    dppgav=0
    write(*,*) "Scan for close primaries:"
    tol1=0.5d0
    tol2=1.0d0
    
    !  note: ndpp1 is the cumulative number of primary pairs closer than tol1
    !        ndpp2 is the cumulative number of primary pairs closer than tol2

    do ip=1,np
        dppmin(ip)=40000 ! ~ earth circumference
        splat=sin(rad*pmlat(ip))
        cplat=cos(rad*pmlat(ip))
        do i=ip+1,np
            dum=cplat*cos(rad*pmlat(i))*cos(rad*(pmlon(i)-pmlon(ip)))
            dum=dum+splat*sin(rad*pmlat(i))
            dum=amin1(dum,1.)
            dpp=acos(dum)/rad
            if(dpp.lt.tol1) then
               ndpp1=ndpp1+1
               prox1(ndpp1,1)=ip
               prox1(ndpp1,2)=i
               prox1d(ndpp1)=dpp
               endif
            if(dpp.lt.tol2) then
               ndpp2=ndpp2+1
               prox2(ndpp2,1)=ip
               prox2(ndpp2,2)=i
               prox2d(ndpp2)=dpp
               endif
            if(ndpp2.gt.MP) stop "error: too many close mascons pairs found (more than MP)"
            dppmin(ip)=amin1(dpp,dppmin(ip))
            dppgmn=amin1(dpp,dppgmn)
            enddo
        dppgav=dppgav+dppmin(ip)
        enddo
     if(ndpp1.gt.0) then
        write(*,*) "Close primaries:"
        write(*,*) "  Separation(deg)  Mascon1      Mascon2       Lat      Lon"
        if(ndpp1.lt.50) then
           do i=1,ndpp1
              i1=prox1(i,1)
              i2=prox1(i,2)
              write(*,"(2x,f12.5,2x,2(i6,a7),2f9.2)") prox1d(i),i1,ptype(i1),i2,ptype(i1),pmlat(i1),pmlon(i1)
              enddo
        else
           write(*,*) "  -- Too many to list"
           endif
        endif
        
     write(*,"(i6,a,f5.2,a)") ndpp1," inter-primary ranges < ",tol1," deg"
     write(*,"(i6,a,f5.2,a)") ndpp2," inter-primary ranges < ",tol2," deg"
     dppgav=dppgav/np
     write(*,"(a40,f10.5,a4)") "Average closest inter-Primary distance ",dppgav," deg"
     write(*,"(a40,f10.5,a4)") "Minimum closest inter-Primary distance ",dppgmn," deg"
     write(*,*)
     
! report small primaries

    ncrit=20
    write(*,*) "Scan for small primaries:"
    write(*,"(a17,i5,a10)") "    cutoff size =",ncrit," ternaries"

    ! setup histogram
    hv=(/20,50,100,200,500,1000,2000,3000,5000,10000/)
    do i=1,10
       hi(i)=0
    enddo
    nhi=0
    ntipmax=0

    do ip=1,np
       ntipmax=max(ntip(ip),ntipmax)
       do i=1,10
          if(ntip(ip).lt.hv(i)) then
             hi(i)=hi(i)+1
             nhi=nhi+1
             go to 100
             endif
          enddo
      
 100   if(ntip(ip).lt.ncrit) then
          nsmall=nsmall+1
          small(nsmall)=ip
          endif
       enddo
       
    if(nsmall.gt.0)then
       write(*,"(a12,i6,a16)") "-> found ",nsmall," small primaries"
       open(17,file="smallpm_aud",status="unknown")
       write(17,*) "    Primary#   #ternaries  latitude    longitude  "
       write(*,*) "  - see smallpm_aud"
       do i=1,nsmall
          ip=small(i)
          write(17,"(2i12,2f12.5)") ip,ntip(ip),pmlat(ip),pmlon(ip)
          enddo
    else
       write(*,*) "   No small primaries found"
       endif
    write(*,*)
    write(*,*) "histogram of primary size ranges"
    write(*,"(10i8)") hv
    write(*,"(10i8)") hi

    write(*,*) "Histogram total: ",nhi
    write(*,*) "Largest primary: ",ntipmax
    write(*,*)

    
!-----------------------------------------------------------------------
    
! mismatch scan: check that primary mascon parameters are consistent
!                 with their contained ternaries

    write(*,*) "Primary-Ternary compatibility check"
    open(16,file="mismatch_aud",status="unknown")
    write(16,*) "P-T mismatch code - 1:lat 2:lon 4:r 10:a 20:h 40:g 100:den 200:pcland 400:region"
    write(16,*) "P-S mismatch code - 1:lat 2:lon 4:r 10:a 20:h 40:g 100:den 200:sap    400:region 800:ntip"

    mispt=0
    misps=0
    
    do ip=1,np
       is=ip
       sumx=0.d0
       sumy=0.d0
       sumz=0.d0
       sum1=0.d0
       sum2=0.d0
       sum4=0.d0
       sum5=0.d0
       sum6=0.d0

!  - calculate check parameters

       do j=1,ntis(is)
          it=ltis(is,j)
          x=dcos(tmlat(it)*rad)*dcos(tmlon(it)*rad)
          y=dcos(tmlat(it)*rad)*dsin(tmlon(it)*rad)
          z=dsin(tmlat(it)*rad)
          sumx=sumx+x*tma(it)
          sumy=sumy+y*tma(it)
          sumz=sumz+z*tma(it)
          sum1=sum1+tmlat(it)*tma(it)
          sum2=sum2+tmlon(it)*tma(it)
          sum4=sum4+tma(it)
          sum5=sum5+tmh(it)
          sum6=sum6+tmg(it)
          enddo

       x=sumx/sum4
       y=sumy/sum4
       z=sumz/sum4
       h=dsqrt(x**2+y**2)
       cmlat=datan2(z,h)/rad
       cmlon=datan2(y,x)/rad
       if(cmlon.lt.0) cmlon=cmlon+360
       
       amlat=sum1/sum4
       amlon=sum2/sum4
    
!    write(*,*) x,y,x,h
!    write(*,*) cmlat,amlat,cmlat-amlat
!    write(*,*) cmlon,amlon,cmlon-amlon

       cma=sum4
       cmh=sum5/ntis(is)
       cmg=sum6/ntis(is)

!  - get geocentric radius of contained ternary closest in latitude to secondary

       difmin=90.d0
       do j=1,ntis(is)
          it=ltis(is,j)
          dif=smlat(is)-tmlat(it)
          if(dif.lt.difmin) then
             difmin=dif
             k=it
             endif
          enddo
       if(difmin.gt.0.1d0) then
          write(*,*) ip,"closest tern lat to secondary is ",difmin," deg (too far?)"
          endif
       cmr=tmr(k)

! - sum land and sea areas of mascon
    
       al=0.d0
       as=0.d0
       do j=1,nsip(ip)
          is=lsip(ip,j)
          do k=1,ntis(is)
             it=ltis(is,k)
             if(tmden(it).le.1001.d0) then
                al=al+tma(it)
             else
                as=as+tma(it)
                 endif
             enddo
          enddo
       cpcland=100*al/(al+as)
       if(cpcland.gt.50) then
          cmden=1000.d0
       else
          cmden=1029.d0
          endif

       
 !  - test P-T parameters match
 
       mismatch=0
       if(abs(pmlat(ip)-cmlat).gt.0.333d0) mismatch=mismatch+1
       if(abs(pmlon(ip)-cmlon).gt.0.333d0) mismatch=mismatch+2
       if(abs(pmr(ip)-cmr).gt.10.d3) mismatch=mismatch+4
       if(abs((pma(ip)-cma)/cma).gt.0.1d0) mismatch=mismatch+8
       if(abs(pmh(ip)-cmh).gt.50.d0) mismatch=mismatch+16
       if(abs(pmg(ip)-cmg).gt.5.d0) mismatch=mismatch+32
       if(abs(pmden(ip)-cmden).gt.5.d0) mismatch=mismatch+64
       if(abs(pcland(ip)-cpcland).gt.5.d0) mismatch=mismatch+128
       if(pregion(ip).ne.sregion(is)) mismatch=mismatch+256

       if(mismatch.gt.0) then
          write(16,"(a13,o5)") " P-T mismatch",mismatch
          write(16,1) ip,ptype(ip),nsip(ip),ntip(ip),pmlat(ip),pmlon(ip),pmr(ip),pma(ip),pmh(ip),pmg(ip),pmden(ip),&
			pcland(ip),tidemask(ip),pregion(ip)
          write(16,1) 0,"Expect",1,ntip(ip),cmlat,cmlon,cmr,cma,cmh,cmg,cmden,&
			cpcland,tidemask(ip),sregion(ip)
          write(16,*) 
          mispt=mispt+1
          endif

 ! - check colours of ternaries in primary

       tpcol_flag=.false.
       tscol_flag=.false.
    
       is1=lsip(ip,1)
       it1=ltis(is,1)
       do j=1,nsip(ip)
          is=lsip(ip,j)
          do k=1,ntis(is)
             it=ltis(is,k)
             if(tpcol(it).ne.tpcol(it1)) tpcol_flag=.true.
             if(tscol(it).ne.tscol(it1)) tscol_flag=.true.
             enddo
          enddo

       if(tpcol_flag) write(16,*) ip, " tpcol varies in mascon"
       if(tscol_flag) write(16,*) ip, " tscol varies in mascon"



       
!  - check P/S parameters agree
  
       mismatch=0
       if(pmlat(ip).ne.smlat(is)) mismatch=mismatch+1
       if(pmlon(ip).ne.smlon(is)) mismatch=mismatch+2
!       if(pmr(ip).ne.smr(is)) mismatch=mismatch+4
       if(pma(ip).ne.sma(is)) mismatch=mismatch+8
       if(pmh(ip).ne.smh(is)) mismatch=mismatch+16
       if(pmg(ip).ne.smg(is)) mismatch=mismatch+32
       if(pmden(ip).ne.smden(is)) mismatch=mismatch+64
       if(sap(is).ne.ip) mismatch=mismatch+128
       if(pregion(ip).ne.sregion(is)) mismatch=mismatch+256
       if(ntip(ip).ne.ntis(is)) mismatch=mismatch+512
       if(nsip(ip).ne.1) mismatch=mismatch+1024
  !    if(ptype(ip).ne.stype(is)) mismatch=mismatch+2048

!      write(16,'(i7,4f7.1)') ip,(pmr(ip)-smr(is))/1000,pmlat(ip),pmr(ip)/1000-6300,smr(is)/1000-6300
       
       if(mismatch.gt.0) then
          write(16,"(i7,a13,o5)") ip," P-S mismatch",mismatch
          misps=misps+1
          endif
    
        enddo
     close(16)
     if(mispt.gt.0) write(*,*) mispt," P-T mismatches found, see mismatch_aud"
     if(misps.gt.0) write(*,*) misps," P-S mismatches found, see mismatch_aud"
     if(misps+mispt.eq.0) write(*,*) " - no mismatches found"
     write(*,*)
     
!-----------------------------------------------------------------------
   
    
! Latitude band structure analysis

    write(*,*) "Structure analysis:"
    lun=11
    open(lun,file="bands_aud",status="unknown")
    write(lun,*) "Primary bands"
    call row_stats(lun,pmlat,pmlon,np,"Primary   ")
    write(lun,*) "Secondary bands"
    call row_stats(lun,smlat,smlon,ns,"Secondary ")
    write(lun,*) "Ternary bands"
    call row_stats(lun,tmlat,tmlon,nt,"Ternary   ")
    close(lun)

    write(*,8) "secondary to primary",dspgav,dspgmx
    write(*,8) "ternary to primary ",dtpgav,dtpgmx
!    if(nxdtp.gt.0) write(*,*) "! ",nxdtp," dtp out of range"
    write(*,8) "ternary to secondary",dtsgav,dtsgmx
    !    if(nxdts.gt.0) write(*,*) "! ",nxdts," dts out of range"

    
    open(13,file="pmstat_aud",status="unknown")
    write(13,*) "   p#    lat     lon      nsip   dspmax    ntip  dtpmax   dtsmax"
    write(13,7) (i,pmlat(i),pmlon(i),nsip(i),dspmax(i),ntip(i), &
                dtpmax(i),dtsmax(i),i=1,np)
!    write(13,7) (i,pmlat(i),pmlon(i),nsip(i),dspmax(i),ntip(i),i=1,np)


    write(*,*) "audit files output:"
    write(*,*) "    bands_aud    - band structure "
    write(*,*) "    checksum     - cksum checksum of header records"
    write(*,*) "    err_aud      - error messages "
    write(*,*) "    headcheck    - trimmed header record"
    write(*,*) "    index_aud    - t-s and t-p indexing errors"
    write(*,*) "    mismatch_aud - S/T parameters incompatible with their enclosing P/S mascon"
    write(*,*) "    outliers_aud - ternaries more than 8 deg from their primary"
    write(*,*) "    pmstat_aud   - primary mascon statistics"
    close(12)
    close(13)
    close(14)

    deallocate(lsip)
    deallocate(ltip)
    deallocate(ltis)


! inspect mascons by number

260 write(*,*)
    write(*,*) "inspect mascons by number (y,n)?"
    read(*,*) choice
    if(choice.eq."y".or.choice.eq."Y") then
    
270 write(*,*) "enter P/S/T and Mascon# (x 1 to exit)"
	read(*,*,err=275,end=275) mtype,i
! convert to upper case
	ic=ichar(mtype)
	if(ic>=97.and.ic<=122) mtype=char(ic-32)
!
	if(mtype.eq."P")then
		write(*,1) i,ptype(i),nsip(i),ntip(i),pmlat(i),pmlon(i),pmr(i),pma(i),pmh(i),pmg(i),pmden(i),&
					pcland(i),tidemask(i),pregion(i)
	else if(mtype.eq."S") then
		write(*,2) i,stype(i),ntis(i),smlat(i),smlon(i),smr(i),sma(i),smh(i),smg(i),smden(i),sap(i),sscol(i),sregion(i)
	else if(mtype.eq."T") then
		write(*,3) i,ttype(i),tmlat(i),tmlon(i),tmr(i),tma(i),tmh(i),tmg(i),tmden(i),tap(i),tas(i),tpcol(i),tscol(i),tregion(i)
	else if(mtype.eq."X") then
		go to 275
	else
		write(*,*) mtype," ?"
		stop "mascon type error"
		endif
	go to 270
        endif

! inspect mascons by location

275 write(*,*) "inspect mascons by location (y/n)?"
    read(*,*) choice
    if(choice.eq."y".or.choice.eq."Y") then
       
280 write(*,*) "lat, lon (x to exit)"
       read(*,*,err=900) xlat,xlon
       if(xlat.gt.90.) stop 
       if(xlon.lt.0.) stop
!
       s1=sin(xlat*rad)
       c1=cos(xlat*rad)
       dxmin=100.d0
       ip=1
       do i=1,np
          dlonr=rad*abs(pmlon(i)-xlon)
          s2=dsin(pmlat(i)*rad)
          c2=dcos(pmlat(i)*rad)
          dx0=(s1*s2+c1*c2*dcos(dlonr))
          dx0=min(dx0,1.d0)
          dx=dacos(dx0)/rad
!	   dx=(xlat-pmlat(i))**2+(xlon-pmlon(i))**2
!          if(i.eq.9456) then
!             write(*,*) dx0,dx,xlat,pmlat(i),xlon,pmlon(i)
!             write(*,*) s1,c1,dlonr,s2,c2
!          endif
          
          if(dx.lt.dxmin) then
             dxmin=dx
             ip=i
          endif
       enddo
       i=ip
       write(*,10) "range: ",sqrt(dxmin)
       write(*,1) i,ptype(i),nsip(i),ntip(i),pmlat(i),pmlon(i),pmr(i),pma(i),pmh(i),pmg(i),pmden(i),&
				pcland(i),tidemask(i),pregion(i)
!
       dxmin=100.d0
       it=1
       do i=1,nt
!          if(tap(i).eq.ip) then
          if(abs(xlat-tmlat(i)).lt.2) then
             dlonr=rad*abs(tmlon(i)-xlon)
             s2=sin(tmlat(i)*rad)
             c2=cos(tmlat(i)*rad)
             dx0=(s1*s2+c1*c2*dcos(dlonr))
             dx0=min(dx0,1.d0)
             dx=acos(dx0)/rad
!             dx=(xlat-tmlat(i))**2+(xlon-tmlon(i))**2
             if(dx.lt.dxmin) then
                dxmin=dx
                it=i
             endif
          endif
       enddo
       i=it
       write(*,10) "range: ",sqrt(dxmin)
       write(*,3) i,ttype(i),tmlat(i),tmlon(i),tmr(i),tma(i),tmh(i),tmg(i),tmden(i),&
				tap(i),tas(i),tpcol(i),tscol(i),tregion(i)
       go to 280
endif


900 stop
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
!    do i=1,nhr
!        write(*,*) headrec(i)
!        enddo

	if(hashcode.eq."#00000") then
		write(*,*) "dummy hashcode found: ",hashcode
	else
		call hash_header(headrec,nhr,checkcode)
		if(hashcode.ne.checkcode) then
			write(*,*) "hashcode mismatch: ",checkcode," != ",hashcode
			endif
		endif

    backspace(lun)
    return

900 stop "error: read_header - end of file"
910 stop "error: read_header - header too long, increase MHR"

    end

!************************************************************

	subroutine hash_header(headrec,nhr,hashcode)

! generates a hashcode from the supplied header lines
!
!    headrec: array of header records to be hashed
!    nhr: number of header records
!    hashcode: 7 digit CRC checksum prefixed with "#"
!     - trailing blanks in header lines are ignored
!     - uses linux cksum function implementing ISO/IEC 8802-3:1989
!     - ten digit decimal hashcode compressed to base36

	implicit none

	integer i,nhr,lun
	integer (kind=8) n
	integer m,rem,base
	character*12 hashcode
	character(150), dimension(1:nhr)    :: headrec
	character*1, dimension(1:12)   ::chr

	lun=99
	open(lun,file="headcheck",status="unknown")
	write(lun,'(a)') (trim(headrec(i)),i=1,nhr)
	close(lun)

	call system("cksum headcheck|awk '{print $1}' >checksum")
	open(lun,file="checksum",status="old")
	read(lun,'(a12)') hashcode
	close(lun)

! convert hashcode to base36

	read(hashcode,*) n
	base=36
	i=0
	do while(n.gt.0)
		m=n/base
		rem=n-m*base
		n=m
		i=i+1
		if(rem.le.9) then
			chr(i)=char(rem+48)
		else
			chr(i)=char(rem+55)
			endif
		enddo
	n=i+1
	do i=1,n-1
		hashcode(i:i)=chr(n-i)
		enddo
	do i=n,12
		hashcode(i:i)=" "
		enddo


	hashcode="#"//hashcode

	return
	end

!************************************************************

subroutine row_stats(lun,lat,lon,n,level)

!
!   compiles row stats for regular mascon structures
!

implicit none

    real (kind=8), parameter      :: pi = 3.141592653589793d0
    integer, parameter      ::  MB=1081
    real (kind=8), dimension(1:n)   :: lat,lon
    real (kind=8), dimension(1:MB) :: rlat,rwid
    real (kind=8) :: scircle
    integer, dimension(1:MB)       :: lrow
    integer lun,i,j,nb,n,maxrow
    character*10 level

1   format(a20," nb>",i6," max # in band =",i6)
6   format(i6,f9.3,i6,f9.3)
    i=0

    rlat(1)=lat(1)
    j=1
    maxrow=1
    lrow(1)=0
    do i=1,n
        if (lat(i).eq.rlat(j)) then
            lrow(j)=lrow(j)+1
        else
! test for irregular mascon structure
            if(j.ge.MB) then
                if(maxrow.le.3) then
                    write(*,*) level,"lat band structure is irregular"
                else
                    write(*,1) " Too many lat bands",j,maxrow
                    endif
                return
                endif
! caculate band structure
            scircle=2*pi*6371*cos(rlat(j)*pi/180)
            rwid(j)=scircle/lrow(j)
            maxrow=max(maxrow,lrow(j))
            j=j+1
            lrow(j)=1
            rlat(j)=lat(i)
        endif
    enddo
    nb=j
    if(nb.gt.MB) write(*,*) "! nb.gt.MB in subroutine row_stats",nb,MB
! pole mascons special case
    scircle=2*pi*6371*cos(((3*rlat(1)+rlat(2))/4)*pi/180)
    rwid(1)=scircle/lrow(1)
    rwid(nb)=scircle/lrow(nb)
!
    write(*,*) nb,level,"lat bands found"
    write(lun,*) nb,level,"lat bands"
    write(lun,*) " belt    lat    nrow  width"
    write(lun,6) (j,rlat(j),lrow(j),rwid(j),j=1,nb)


    return
    end
