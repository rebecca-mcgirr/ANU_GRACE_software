      program merge_mascons
!     *******
!
!  Reads and rewrites a mascon file after merging two of its primary mascons into one

!
!  usage: merge_mascons <input mascon_file> <mascon1> <mascon2> [<newmax>]
!		- mascon2 will be merged into mascon1 and its successors renumbered
!		- if mascon1=mascon2 the mascon file is just reordered sequentially by primary mascon
!		- maxtip,maxsip increased to <newmax> if specified
!
!  output to mascons_<mascon1>m<mascon2>
!
! main arrays and variables:
! --------------------------
! ipm1 - primary mascon to grow
! ipm2 - primary mascon to merge into ipm1
! newmax - expanded array size to allow merged mascon larger than current maximum
!
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
! nc_sap,nc_tap,nc_tas - number of corrected mascon connections made

    implicit none
    real (kind=8), parameter    :: pi = 3.141592653589793d0
    integer, parameter      :: MHR = 200

	real (kind=8), dimension(:), allocatable    :: pmlat,pmlon,pmr,pma,pmh,pmg,pmden,pcland
	real (kind=8), dimension(:), allocatable    :: smlat,smlon,smr,sma,smh,smg,smden
	real (kind=8), dimension(:), allocatable    :: tmlat,tmlon,tmr,tma,tmh,tmg,tmden
	real (kind=8), dimension(:), allocatable    :: dspmax,dtpmax,dtsmax
	integer, dimension(:), allocatable          :: nsip,ntip,tidemask
	integer, dimension(:), allocatable          :: ntis,sap,sscol
	integer, dimension(:), allocatable          :: tap,tas,tpcol,tscol
	integer, dimension(:), allocatable          :: isip,itip,kips,kpsc
	integer, dimension(:), allocatable          :: itis,kist,kstc
	integer, dimension(:), allocatable          :: hsind,hpind
	integer, dimension(:), allocatable          :: scheck,tcheck

    integer, dimension(:,:), allocatable 		:: lsip,ltip,ltis

	character*20, dimension(:), allocatable	   :: pregion
	character*20, dimension(:), allocatable	   :: sregion
	character*20, dimension(:), allocatable	   :: tregion
	character(6), dimension(:), allocatable    :: ptype
	character(6), dimension(:), allocatable    :: stype
	character(6), dimension(:), allocatable    :: ttype

	integer     MP,MS,MT,MM
	integer 	ipm1,ipm2

	character(150), dimension(1:MHR)  :: headrec
	character*1 mtype
	character*12 hashcode
    character*150 line


    integer np,ns,nt,nm,nh,nn,ip,is,it
    integer maxsip,maxtip,maxtis
	integer oldmax,newmax
	integer nc_sap,nc_tap,nc_tas

! working variables

	real (kind=8) sum1,sum2,sum3,sum4,sum5,a1,a2,dif,difmin
    integer i,j,k,is1,is2,it1
    integer offset,ipo,iso
	character*120 infile,outfile
	character*14 ext
	character*6 arg1,arg2,arg3
	logical expand


1   format(i7,2x,a6,2i8,2f10.4,f11.1,f18.0,2f9.1,f6.0,f6.1,i6,2x,a20)
2   format(i7,2x,a6,i8,8x,2f10.4,f11.1,f18.0,2f9.1,f6.0,2i6,2x,a20)
3   format(i7,2x,a6,2f10.4,f11.1,f15.0,2f9.1,f6.0,4i6,4x,a20)

4   format(a150)
5   format(9x,a1)
6	format(a12,2i8,i10,3i7)
7   format(9x,i4)
8	format(a27,a15,a8,i5,a20,i5)
9   format(i8,$)

! parse command line

	if(iargc().lt.3) then
		write(*,*)
        write(*,*) "usage: merge_mascons <input mascon_file> <mascon1> <mascon2>"
        write(*,*) "	 - mascon2 will be merged into mascon1 and its successors renumbered"
		write(*,*)
		stop
        endif

	call getarg(1,infile)
	write(*,*) "Input mascon file:  ",infile

    call getarg(2,arg1)
	read(arg1,*) ipm1

	call getarg(3,arg2)
	read(arg2,*) ipm2

	expand=.false.
	if(iargc().eq.4) then
		call getarg(4,arg3)
		read(arg3,*) newmax
		expand=.true.
		endif

	ext="_"//trim(arg1)//"m"//trim(arg2)
	outfile=trim(infile)//trim(ext)
	write(*,*) "Output mascon file: ",outfile
	write(*,*) "Note: primary mascons will be sorted on output"

! open mascon files

	open(10,file=infile,status="old")
	write(*,*) "Header line 1:"

! Allocate arrays

	read(10,*) hashcode,MP,MS,MT,maxsip,maxtip,maxtis
	write(*,*) hashcode,MP,MS,MT,maxsip,maxtip,maxtis
	oldmax=maxsip
	if(expand) then
		write(*,*) "expanding maxtip,maxsip to ",newmax
		maxtip=newmax
		maxtis=newmax
		endif
	rewind(10)

	allocate(pmlat(MP),pmlon(MP),pmr(MP),pma(MP),pmh(MP),pmg(MP),pmden(MP),pcland(MP))
	allocate(smlat(MS),smlon(MS),smr(MS),sma(MS),smh(MS),smg(MS),smden(MS))
	allocate(tmlat(MT),tmlon(MT),tmr(MT),tma(MT),tmh(MT),tmg(MT),tmden(MT))

!	allocate(dspmax(MP),dtpmax(MP),dtsmax(MP))
	allocate(nsip(MP),ntip(MP),tidemask(MP))
	allocate(ntis(MS),sap(MS),sscol(MS))
	allocate(tap(MT),tas(MT),tpcol(MT),tscol(MT))
	allocate(isip(MP),itip(MP),kips(MP),kpsc(MP))
	allocate(itis(MS),kist(MS),kstc(MS))
!	allocate(hpind(MS),scheck(MS))
!	allocate(hsind(MT),tcheck(MT))

	allocate(pregion(MP))
	allocate(sregion(MS))
	allocate(tregion(MT))
	allocate(ptype(MP))
	allocate(stype(MS))
	allocate(ttype(MT))

	allocate(lsip(MP,maxsip))
!	allocate(ltip(MS,maxtip))
	allocate(ltis(MP,maxtis))

	MM=MP+MS+MT
	write(*,*) "total mascons expected:",MM

! setup

	np=0
	ns=0
	nt=0
	nm=0
	j=0
	ip=0
	is=0
	it=0
	nc_sap=0
	nc_tap=0
	nc_tas=0
	do i=1,MP
		nsip(i)=0
		ntip(i)=0
		enddo

! read header

	i=0
	read(10,4,end=920) line
	do while(line(1:1).eq."#")
		i=i+1
       	headrec(i)=trim(line)
		read(10,4,end=920) line
		enddo
	nh=i
    backspace(10)

!-----------------------------------------------------

! read mascon file

	write(*,'(a,$)') " primary read loop: "
    do i=1,MP
!        read(10,*) ip,ptype(ip),nsip(ip),ntip(ip),pmlat(ip),pmlon(ip),pmr(ip),pma(ip),pmh(ip),pmg(ip),pmden(ip),pcland(ip),tidemask(ip),pregion(ip)
        read(10,1) ip,ptype(ip),nsip(ip),ntip(ip),pmlat(ip),pmlon(ip),pmr(ip),pma(ip),pmh(ip),pmg(ip),pmden(ip),&
					 pcland(ip),tidemask(ip),pregion(ip)

		nm=nm+1
        np=max(ip,np)


		if(nsip(ip).gt.1) then
			write(*,*) "multiple secondaries in primary ",ip
			stop
			endif

    	do j=1,nsip(ip)
!       	 read(10,*) is,stype(is),ntis(is),smlat(is),smlon(is),smr(is),sma(is),smh(is),smg(is),smden(is),sap(is),sscol(is),sregion(is)
        	read(10,2) is,stype(is),ntis(is),smlat(is),smlon(is),smr(is),sma(is),smh(is),smg(is),smden(is),sap(is),&
						sscol(is),sregion(is)

			nm=nm+1
        	ns=max(is,ns)
			lsip(ip,j)=is
			if(sap(is).ne.ip) then
!				write(*,*) "secondary ",is," not assocated with its primary ",ip
				nc_sap=nc_sap+1
				endif

    		do k=1,ntis(is)
!				read(10,*) it,ttype(it),tmlat(it),tmlon(it),tmr(it),tma(it),tmh(it),tmg(it),tmden(it),tap(it),tas(it),tpcol(it),tscol(it),tregion(it)
        		read(10,3) it,ttype(it),tmlat(it),tmlon(it),tmr(it),tma(it),tmh(it),tmg(it),tmden(it),tap(it),tas(it),&
						tpcol(it),tscol(it),tregion(it)

				nm=nm+1
        		nt=max(it,nt)
				ltis(is,k)=it
!				ltip(ip,k)=it	! this assumes ony one secondary per primary

! check P/S associations of each ternary
				if(tap(it).ne.ip) then
!					write(*,*) "ternary ",it," not assocated with its primary ",ip
					nc_tap=nc_tap+1
					endif
				if(tas(it).ne.is) then
!					write(*,*) "ternary ",it," not assocated with its secondary ",is
					nc_tas=nc_tas+1
					endif

				enddo
			enddo
    	if(mod(i,1000)==0) write(*,9) i
		enddo

    close(10)
	write(*,*)
	write(*,*) "mascon range (P,S,T): ",np,ns,nt
	write(*,*) "total mascons found   :",nm

!-----------------------------------------------------

! report misconnections in sap/tap/tas arrays

	if(nc_sap+nc_tap+nc_tas.gt.0) then
		write(*,*) nc_sap," S-P misconnections found"
		write(*,*) nc_tap," T-P misconnections found"
		write(*,*) nc_tas," T-S misconnections found"
		write(*,*) " (misconnections will be repaired on output)"
		endif

! verify mascon compatibility for merger

	if(ptype(ipm1).ne.ptype(ipm2)) then
               write(*,*) "Warning: mascon types don't match"
               write(*,*) "primary ",ipm1," is ",ptype(ipm1)
               write(*,*) "primary ",ipm2," is ",ptype(ipm2)
               !stop
        endif

	if(pregion(ipm1).ne.pregion(ipm2)) then
		write(*,*) "Warning: mascon regions don't match"
		write(*,*) "primary ",ipm1," is ",pregion(ipm1)
		write(*,*) "primary ",ipm2," is ",pregion(ipm2)
		!stop
		endif

! test that merged mascon is not too large for arrays

        nn=ntip(ipm1)+ntip(ipm2)
        write(*,'(a,i5,a)') " merged mascon will contain ",nn," ternaries"
	if(nn.gt.maxtis) then
	     write(*,*) "Error: merged mascon is too large for arrays ( #ternaries > maxtis )"
 	     write(*,'(i5,a,i5,a,i5)') ntip(ipm1)," +",ntip(ipm2)," >",maxtis
		stop
		endif

!-----------------------------------------------------

! merge mascons

	if(ipm1.ne.ipm2) then
		write(*,'(a,i5,a,i5)') " merging mascon",ipm2," into",ipm1


! change colours of ternaries in retired primary (ipm2)

		is1=lsip(ipm1,1)
		it1=ltis(is1,1)
		do j=1,nsip(ipm2)
			is2=lsip(ipm2,j)
			do k=1,ntis(is2)
				it=ltis(is2,k)
				tpcol(it)=tpcol(it1)
				tscol(it)=tscol(it1)
				enddo
			enddo
! update list index for merged primary

		do i=1,ntis(is2)
			ltis(is1,ntis(is1)+i)=ltis(is2,i)
			write(99,*) is1,is2,ntis(is1),i,ltis(is2,i)
			enddo
!update merged mascon counts

!!		nsip(ipm1) is unchanged
		ntip(ipm1)=ntip(ipm1)+ntip(ipm2)
		ntis(is1)=ntip(ipm1)

! update merged secondary's parameters
! - we assume there will only ever be one secondary in each primary
!		with the same parameters as its primary

!!		stype(is1) is unchanged

		do j=1,ntis(is1)
			it=ltis(is1,j)
			sum1=sum1+tmlat(it)*tma(it)
			sum2=sum2+tmlon(it)*tma(it)
			sum3=sum3+tma(it)
			sum4=sum4+tmh(it)
			sum5=sum5+tmg(it)
			enddo

		smlat(is1)=sum1/sum3
		smlon(is1)=sum2/sum3
!		smr(is1)=
		sma(is1)=sma(is1)+sma(is2)
		smh(is1)=sum4/ntis(is1)
		smg(is1)=sum5/ntis(is1)

! secondary geocentric radius copied from ternary closest in latitude

		difmin=90.d0
		do j=1,ntis(is1)
			it=ltis(is1,j)
			dif=smlat(is1)-tmlat(it)
			if(dif.lt.difmin) then
				difmin=dif
				k=j
				endif
			enddo
		if(difmin.gt.0.1d0) then
			write(*,*) "closest tern lat to merged secondary is ",difmin," deg (too far?)"
			endif
		smr(is1)=tmr(k)

!!		smden(is1) is unchanged
!!		sap(is1) is unchanged
!!		sscol(is1) is unchanged
!!		sregion(is1) is unchanged

! update merged primary's parameters
!	- most are the same as the contained secondary

		a1=pma(ipm1)
		a2=pma(ipm2)
		pmlat(ipm1)=smlat(is1)
		pmlon(ipm1)=smlon(is1)
		pmr(ipm1)=smr(is1)
		pma(ipm1)=sma(is1)
		pmh(ipm1)=smh(is1)
		pmg(ipm1)=smg(is1)
		pcland(ipm1)=(pcland(ipm1)*a1+pcland(ipm2)*a2)/(a1+a2)

!!		pmden(ipm1) is unchanged
!!		tidemask(ipm1) is unchanged
!!		pregion(ipm1) is unchanged

		endif

!-----------------------------------------------------

! write modified mascon file
!  - P/S mascon membership is updated on output
        open(11,file=outfile,status='unknown') 
  
	if(expand) then
		maxtip=max(oldmax,ntip(ipm1))
		maxtis=maxtip
		write(*,*) "resetting maxtip,maxtis to ",maxtip
		endif
! PT180823: subtract one from the number of primary and secondary mascons written to the header (since we are killing one mascon)
	!write(headrec(1),6) hashcode,MP-1,MS-1,MT,maxsip,maxtip,maxtis
        ! PT200220: hardwire to newmax value if expand=true
        if(.not. expand)then
	  write(headrec(1),6) hashcode,MP-1,MS-1,MT,maxsip,maxtip,maxtis
        else   
	  write(headrec(1),6) hashcode,MP-1,MS-1,MT,maxsip,newmax,newmax 
        endif

	nh=nh+1
	write(headrec(nh),8) "# merge_mascons  , region: ",pregion(ipm1)," mascon ",ipm2," merged into mascon ",ipm1
	do i=1,nh
		write(11,"(a)") trim(headrec(i))
		enddo

    offset=0
	write(*,'(a,$)') " primary write loop: "
	do ip=1,MP

        if(ip.ne.ipm2) then
		ipo=ip-offset
		write(11,1)  ipo,ptype(ip),nsip(ip),ntip(ip),pmlat(ip),pmlon(ip),pmr(ip),pma(ip),pmh(ip),pmg(ip),pmden(ip),&
							pcland(ip),tidemask(ip),pregion(ip)

		do j=1,nsip(ip)
			is=lsip(ip,j)
			iso=is-offset
			write(11,2)  iso,stype(is),ntis(is),smlat(is),smlon(is),smr(is),sma(is),smh(is),smg(is),smden(is),ipo,&
							sscol(is),sregion(is)

			do k=1,ntis(is)
				it=ltis(is,k)
!				write(99,*) ip,is,k,it
        		write(11,3)  it,ttype(it),tmlat(it),tmlon(it),tmr(it),tma(it),tmh(it),tmg(it),tmden(it),ipo,iso,&
								tpcol(it),tscol(it),tregion(it)
				enddo
			enddo
        else
            offset=1
            endif
		if(mod(ip,1000)==0) write(*,9) ip
		enddo



    close(11)
    write(*,*)
    write(*,*) "new mascon count (P,S,T): ",MP-offset,MS-offset,MT

!    call status_update('STATUS','UTIL','merge_mascons',outfile,"End program",0)

    stop

920 stop "mascon file error: no mascons?"
    end

!************************************************************


