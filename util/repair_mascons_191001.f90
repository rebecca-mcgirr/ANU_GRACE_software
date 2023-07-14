
program repair_mascons
!        *************
!
!  Reads and rewrites a mascon file repairing ternary-secondary connections
!
!     Also corrects:
!        - inconsistencies between primary and mean ternary mascon parameters
!               (lat,lon,radius,area,height,geoid,pcland)
!        - inconsistent ternary colours in one primary
!        - differences between primary and secondary mascon parameters
!               (lat,lon,radius,area,height,geoid)
!
!     Header line appended and hashcode recalculated
!
!  Known errors generated when primaries are reshaped during stage 4:
!     - areas not preserved
!     - height, geoid not updated
!     - primary radius reset to 6356.7
!     - all secondaries put into primary #1  
!     - position errors near zero meridian
  
!  Output to "mascons<...>_repaired"
!     warn_log lists position corrections >1deg
!     repair_log lists all corrections
!
!  HM190725

  
! Main arrays and variables:
! --------------------------
! !
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


    use mascon_setup
  
    implicit none

    integer, dimension(:,:), allocatable 	:: lsip,ltip,ltis
    integer, dimension(:), allocatable 		:: ltin,ltex
    integer, dimension(10)   :: list
    
    real (kind=8), dimension(1:20)   :: vlon,vlat   ! polygon vertices
    integer     MP,MS,MT,MM,nv
    integer 	ipm1,ipm2,ipm3

    character*1 mtype
    character*12 hashcode
    character*150 line


    integer np,ns,nt,nm,nh,nn,ip,is,it
    integer maxsip,maxtip,maxtis
    integer oldmax,newmax
    integer nc_sap,nc_tap,nc_tas

! working variables

    real (kind=8) sum1,sum2,sum3,sum4,sum5,a1,a2,dif,difmin,area
    integer i,j,k,is1,is2,it1,nl
    integer narg
    character*120 infile,outfile,polyfile
    character*14 ext
    character*6 arg2,arg3,arg4,arg
    logical expand,inpolyc,split

    data list /3,156,1361,470,1355,16,79,1,1,1/
    data nl/8/


1   format(i7,2x,a6,2i8,2f10.4,f11.1,f18.0,2f9.1,f6.0,f6.1,i6,2x,a15)
2   format(i7,2x,a6,i8,8x,2f10.4,f11.1,f18.0,2f9.1,f6.0,2i6,2x,a15)
3   format(i7,2x,a6,2f10.4,f11.1,f15.0,2f9.1,f6.0,4i6,4x,a15)

4   format(a150)
5   format(9x,a1)
6   format(a12,2i8,i10,3i7)
7   format(9x,i4)
8   format(a27,a15,a8,i5,a20,i5)
9   format(i8,$)

! parse command line

    narg=iargc()
    if(narg.lt.1) then
       write(*,*)
       write(*,*) "usage: repair_mascons <input mascon_file>"
       write(*,*)
       stop
       endif

    call getarg(1,infile)
    write(*,*) "Input mascon file:  ",infile

    ext="_repaired"
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
       write(*,"(a40,i5,a20)") "Expanding max ternaries in primary from",maxtip," - enter new maxtip:"
       read(*,*) newmax      
!      write(*,*) "expanding maxtip,maxsip to ",newmax
       maxtip=newmax
       maxtis=newmax
       endif
!    rewind(10)

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
    allocate(ltis(MS,maxtis))
!	allocate(ltip(MS,maxtip))
!       allocate(ltin(maxtis))             ! list of ternaries inside polygon
!       allocate(ltex(maxtis))             ! list of ternaries outside polygon
    

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

! read header comments

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

! read mascon file and build connection lists

    write(*,'(a,$)') " primary read loop: "
    do i=1,MP
       read(10,1) ip,ptype(ip),nsip(ip),ntip(ip),pmlat(ip),pmlon(ip),pmr(ip),pma(ip),pmh(ip),pmg(ip),pmden(ip),&
					 pcland(ip),tidemask(ip),pregion(ip)

       nm=nm+1
       np=max(ip,np)


       if(nsip(ip).gt.1) then
          write(*,*) "multiple secondaries in primary ",ip
          stop
          endif

       do j=1,nsip(ip)
          read(10,2) is,stype(is),ntis(is),smlat(is),smlon(is),smr(is),sma(is),smh(is),smg(is),smden(is),sap(is),&
						sscol(is),sregion(is)

          nm=nm+1
          ns=max(is,ns)
          lsip(ip,j)=is
          if(sap(is).ne.ip) then
             nc_sap=nc_sap+1
             endif

          do k=1,ntis(is)
             read(10,3) it,ttype(it),tmlat(it),tmlon(it),tmr(it),tma(it),tmh(it),tmg(it),tmden(it),tap(it),tas(it),&
						tpcol(it),tscol(it),tregion(it)

             nm=nm+1
             nt=max(it,nt)
             ltis(is,k)=it
!	     ltip(ip,k)=it	! program assumes ony one secondary per primary

! check P/S associations of each ternary
             if(tap(it).ne.ip) then
                nc_tap=nc_tap+1
                endif
             if(tas(it).ne.is) then
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

!-----------------------------------------------------

! open log files

    open(12,file="repair_log",status="unknown")
    write(12,"(a)") " mascon changes from repair_mascons"
    write(12,"(a)") " ip,shift(deg)"
    open(13,file="warn_log",status="unknown")
    write(13,"(a)") " large mascon shifts (> 1deg)from repair_mascons"
    write(13,"(a)") " is,smlat,cmlat,smlon,cmlon,dlat,dlon,arcdist"

! repair mascon parameters

    write(*,'(a)') " resetting primary/secondary mascon parameters"
    do ip=1,MP  
       call repair_ps(ip,lsip,ltis,MP,MS,MT,maxsip,maxtis)
    enddo
    
!-----------------------------------------------------

! update header

    nh=nh+1
    headrec(nh)="# repair_mascons , Primary+Secondary lat,lon,radius,area,height,geoid,pcland recalculated from ternaries"
    call hash_header(headrec,nh,hashcode)
    
! write modified mascon file
!  - P/S mascon membership is updated on output

      open(11,file=outfile,status='unknown') 
  
      write(11,6) hashcode,MP,MS,MT,maxsip,maxtip,maxtis
      do i=1,nh
         write(11,"(a)") trim(headrec(i))
         enddo

      write(*,'(a,$)') " primary write loop: "
      do ip=1,MP

            write(11,1)  ip,ptype(ip),nsip(ip),ntip(ip),pmlat(ip),pmlon(ip),pmr(ip),pma(ip),pmh(ip),pmg(ip),pmden(ip),&
							pcland(ip),tidemask(ip),pregion(ip)

            do j=1,nsip(ip)
               is=lsip(ip,j)
               write(11,2)  is,stype(is),ntis(is),smlat(is),smlon(is),smr(is),sma(is),smh(is),smg(is),smden(is),ip,&
							sscol(is),sregion(is)

               do k=1,ntis(is)
                  it=ltis(is,k)
 !		write(99,*) ip,is,k,it
                  write(11,3)  it,ttype(it),tmlat(it),tmlon(it),tmr(it),tma(it),tmh(it),tmg(it),tmden(it),ip,is,&
								tpcol(it),tscol(it),tregion(it)
                  enddo
               enddo
         if(mod(ip,1000)==0) write(*,9) ip
         enddo



    close(11)
    close(12)
    close(13)
    write(*,*)
    write(*,*) "new mascon count (P,S,T ): ",MP,MS,MT
    write(*,*) "new hashcode: ",hashcode,"    number of header lines: ",nh

!    call status_update('STATUS','UTIL','merge_mascons',outfile,"End program",0)

    stop

920 stop "mascon file error: no mascons?"
    end



!************************************************************

    subroutine repair_ps(ip,lsip,ltis,MP,MS,MT,maxsip,maxtis)

! recalculate primary and secondary parameters from their contained ternaries
!      - parameters revised: lat,lon,radius,area,height,geoid,pcland

    use mascon_setup

    integer, dimension(1:MP,1:maxsip)   :: lsip
    integer, dimension(1:MS,1:maxtis)   :: ltis
    integer   MP,MS,MT,maxsip,maxtip
    integer   ip
    
    real (kind=8) sum1,sum2,sum3,sum4,sum5,sum6
    real (kind=8) x,y,z,h,a1,a2,dif,difmin,al,as
    real (kind=8) sumx,sumy,sumz,amlat,amlon,cmlat,cmlon,dist,dlat,dlon
    integer   np,ns,nt,nm,nh,nn,is,it
    integer   i,j,k,it1


!    write(*,*) ip,lsip(ip,1),ltis(ip,1),MP,MS,MT,maxsip,maxtis
    
    is=lsip(ip,1)
       sumx=0.d0
       sumy=0.d0
       sumz=0.d0
       sum1=0.d0
       sum2=0.d0
       sum4=0.d0
       sum5=0.d0
       sum6=0.d0

! set colours of ternaries to match first ternary in primary

    it1=ltis(is,1)
    do j=1,nsip(ip)
       do k=1,ntis(is)
          it=ltis(is,k)
          tpcol(it)=tpcol(it1)
          tscol(it)=tscol(it1)
          enddo
       enddo

! calculate secondary parameters
! - assume there will only ever be one secondary in each primary
! - first repair the secondary parameters, then set the primary to match

!!		stype(is) is unchanged
200    rad=pi/180
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
    
! calculate arc mean as a check on chord mean (cmlat,cmlon) of ternary positions
!  - chord mean is used for the update

    amlat=sum1/sum4
    amlon=sum2/sum4
!    write(*,*) x,y,x,h
!    write(*,'(3(a,f8.3))') "new latitude:  chord",cmlat,"   arc",amlat,"   diff",cmlat-amlat
!    write(*,'(3(a,f8.3))') "new longitude: chord",cmlon,"   arc",amlon,"   diff",cmlon-amlon

    dist=sin(smlat(is)*rad)*sin(cmlat*rad)
    dist=dist+cos(smlat(is)*rad)*cos(cmlat*rad)*cos((smlon(is)-cmlon)*rad)
    dist=acos(dist)/rad
    dlat=smlat(is)-cmlat
    dlon=smlon(is)-cmlon
    if(dist.gt.1.) then
       write(13,'(i6,7f9.3)') is,smlat(is),cmlat,smlon(is),cmlon,dlat,dlon,dist
    endif
    write(12,'(i6,f9.3)') is,dist
 
    smlat(is)=cmlat
    smlon(is)=cmlon
    sma(is)=sum4
    smh(is)=sum5/ntis(is)
    smg(is)=sum6/ntis(is)

! P/S geocentric radius copied from ternary closest in latitude

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
       write(*,*) "closest contained tern lat to revised secondary is ",difmin," deg (too far?)"
       endif
       smr(is)=tmr(k)

!!	sap(is) is unchanged
!!	sscol(is) is unchanged
!!	sregion(is) is unchanged

! sum land and sea areas of mascon
    
    al=0.d0
    as=0.d0
    do j=1,ntis(is)
       it=ltis(is,j)
       if(tmden(it).le.1001.d0) then
          al=al+tma(it)
       else
          as=as+tma(it)
       endif    
    enddo
    pcland(ip)=100*al/(al+as)

    if(smden(is).gt.1020.and.pcland(ip).gt.50) then
       write(*,"(a,i7,f7.0,a,f4.0,a)") "Ocean mascon ",is,smden(is)," is majority land ",pcland(ip),"% - density fixed"
       smden(is)=1000.d0
    endif
    if(smden(is).lt.1020.and.pcland(ip).lt.50) then
       write(*,"(a,i7,f7.0,a,f4.0,a)") "Land  mascon ",is,smden(is)," is majority ocean",pcland(ip),"% - density fixed"
       smden(is)=1029.d0
    endif
       
! update primary's parameters
!	- most are the same as the contained secondary
          
    pmlat(ip)=smlat(is)
    pmlon(ip)=smlon(is)
    pmr(ip)=smr(is)
    pma(ip)=sma(is)
    pmh(ip)=smh(is)
    pmg(ip)=smg(is)
    pmden(ip)=smden(is)

!!	tidemask(ip) is unchanged
!!	pregion(ip) is unchanged

    return
    end subroutine repair_ps


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
