program initialize_mascons
!******
!
!  generates starting set of primary and ternary mascons
!		primary and secondary mascons are identical at this stage
!		sets initial geometry of primary,secondary and ternary mascons
!		values calculated for coordinates, area and geocentric radius of mascons
!		dummy calues used for height, geoid height, density and region
!			(they will be set at stage2 by classify_mascons)
! 
!	Note this program uses co-latitude, not geographic latitude
!
implicit none
real (kind=8), parameter    :: pi = 3.141592653589793d0
real (kind=8), parameter    :: Ra = 6378.137d3, inflat = 298.257223563d0
integer, parameter      :: tlat_wid_min = 10
integer, parameter      :: plat_wid_deg = 3
integer, parameter      :: MP = 4600, MPB=61
integer, parameter      :: MS = 18000, MSB=361
integer, parameter      :: MT = 1500000, MTB=1081 ,MTIP=2000

real (kind=8), dimension(1:MP)   :: pmlat,pmlon,pmr,pma,pmh,pmden,rpn
real (kind=8), dimension(1:MS)   :: smlat,smlon,smr,sma,smh,smden
real (kind=8), dimension(1:MT)   :: tmlat,tmlon,tmr,tma,tmh,tmden
integer, dimension(1:MP)         :: pmn,nsip,ntip,tidemask
integer, dimension(1:MS)         :: smn,ntis,sap,sscol
integer, dimension(1:MT)         :: tap,tas,tpcol,tscol
character*1 mtype
character*12 hashcode
character*12 region
character*180 line
character(6), dimension(1:MP)   :: ptype
character(6), dimension(1:MS)   :: stype
character(6), dimension(1:MT)   :: ttype

integer, dimension(1:MPB)       :: npsub
integer, dimension(1:MTB)       :: npzone
integer, dimension(1:MTB)       :: ntzone,ntsub

real (kind=8) Rb,Re,Rg,area,pcland
real (kind=8) e,area_WGS84,sum,ecorr
real (kind=8) dplat,dplon,dplatkm,platb
real (kind=8) dtlat,dtlon,dtlatkm,tblat
integer ntband,ntern,npband,nprm,plat_wid_min,rho,col_range


integer, dimension(1:MTIP,1:MP)  :: itip
integer, dimension(1:MP)  :: itp
data col_range/1024/

! working variables
real (kind=8) rad,sinlat,coslat,blat,pmat,pmas,tmat,rnum,ght0
integer i,j,k,ntchk,nptm,ip,jp,it,its,itf,ip0,ipc,itc,j0,itb,ntipc
character(8)  :: date
character(10) :: time
character(5)  :: timezone

! formats

1   format(i7,2x,a6,2i8,2f10.4,f11.1,f18.0,2f9.1,f6.0,f6.1,i6,2x,a20)
2   format(i7,2x,a6,i8,8x,2f10.4,f11.1,f18.0,2f9.1,f6.0,2i6,2x,a20)
3   format(i7,2x,a6,2f10.4,f11.1,f15.0,2f9.1,f6.0,4i6,12x,a6)
4   Format(a12,2i6,i9,3i6)
5   format(a,1x,a8,1x,a10,1x,a5)
8   format(a8,$)
9   format(i8,$)


open(10,file="mascons_stage1",status="unknown")
open(11,file="junk",status="unknown")

! initialize

rad=pi/180.d0
Rb=Ra*(1-1/inflat)
Re=(Ra*Ra*Rb)**(1/3.d0)
write(*,*) Ra,Rb,Re
plat_wid_min=60*plat_wid_deg
nptm=plat_wid_min/tlat_wid_min
ipc=0
itc=0
pmat=0.d0
pmas=0.d0
tmat=0.d0
pcland=100.d0
ght0=0.d0

e=sqrt(1-(Rb/Ra)**2)
!write(*,*) " ellipticity   : ",e
area_WGS84=2*pi*(Ra**2)*(1+(1-e**2)*datanh(e)/e)
write(*,*) "WGS84 surface area ",area_WGS84
!area_WGS84=510065621724088.44

do i=1,MP
	pma(i)=0.d0
	enddo


! setup band structure

call calc_lat_bands(plat_wid_min,npband,nprm,npzone,dplat)
write(*,*) npband," primary bands",nprm," mascons"
write(*,*)
if(nprm.gt.MP) stop "error: parameter MP too small"
if(nprm.gt.MS) stop "error: parameter MS too small"
call calc_lat_bands(tlat_wid_min,ntband,ntern,ntzone,dtlat)
write(*,*) ntband," ternary bands",ntern," mascons"
write(*,*)
if(ntern.gt.MT) stop "error: parameter MT too small"

write(*,*) "Primary latitude spacing:",dplat,"deg"
write(*,*) "Ternary latitude spacing:",dtlat,"deg"
write(*,*)

! write header line
! Note: secondaries copy primaries in the stage1 mascon file
!       hashcode is just a placeholder at this stage

hashcode="#00000"
write(10,4) hashcode,nprm,nprm,ntern,0,0,0
call date_and_time(date,time,timezone)
write(10,5) "# setup by initialize_mascons ",date,time,timezone

do ip=1,nprm
    call random_number(rpn(ip))
    call random_number(rnum)
    sscol(ip)=col_range*rnum
    enddo
region="Region      "

! create ternary mascons
write(*,*) "creating ternary mascons ..."

it=1
ttype(it)="Tern "
tmlat(it)=0.d0
tmlon(it)=0.d0
tmr(it)=Rb
tma(it)=2*pi*Rb**2*(1-dcos(dtlat*rad/2))
tmh(it)=0.d0
tmden(it)=1000
tap(it)=0
tas(it)=0
tscol(it)=mod(it,1023)+1


ntsub(1)=0
do i=2,ntband-1
    tblat=(i-1)*dtlat
    dtlon=360.d0/ntzone(i)
!    write(11,*) i,ntzone(i),dtlon
    call calc_geocentric(tblat,Rg)
    area=2*pi*Rg**2*(dcos((i-1.5d0)*dtlat*rad)-dcos((i-0.5d0)*dtlat*rad))/ntzone(i)
    ntsub(i)=ntsub(i-1)+ntzone(i-1)
    do j=1,ntzone(i)
        it=it+1
        tmlat(it)=tblat
        tmlon(it)=j*dtlon
        ttype(it)="Tern "
        tmr(it)=Rg
        tma(it)=area
        tmh(it)=0.
        tmden(it)=1000
        tap(it)=0
        tas(it)=0
        tscol(it)=mod(it,1023)+1
!        if(mod(it,100000)==0) write(*,9) it

        enddo
    enddo

ntsub(ntband)=ntsub(ntband-1)+ntzone(ntband-1)

ntchk=it+1
it=ntern
ttype(it)="Tern "
tmlat(it)=180.d0
tmlon(it)=0.d0
tmr(it)=Rb
tma(it)=2*pi*Rb**2*(1-dcos(dtlat*rad/2))
tmh(it)=0.d0
tmden(it)=1000
tap(it)=0
tas(it)=0
tscol(it)=mod(it,1023)+1

! apply WGS84 ellipsoid area correction
!
do it=1,ntern
    sum=sum+tma(it)
    enddo
ecorr=area_WGS84/sum
write(*,*) " WGS84 area correction ",1/ecorr
do it=1,ntern
    tma(it)=ecorr*tma(it)
    enddo

write(*,*) "ternary mascon count:",ntern
if(ntchk.ne.ntern) write(*,*) "!!error ternary mascon check count differs:",ntchk


! assign ternaries to primaries

write(*,*)
write(*,*) "building primary mascons ..."
write(*,*) nptm,"ternary bands in a primary band"

! build north polar cap primary

ip=1
ptype(ip)="Ppolar"
pmlat(ip)=0.d0
pmlon(ip)=0.d0
pmr(ip)=Rb
pmh(ip)=0.d0
pmden(ip)=1000
tidemask(ip)=0
nsip(ip)=1
ntip(ip)=0
its=1
do j=1,nptm/2
    do k=1,ntzone(j)
        it=ntsub(j)+k
        tap(it)=ip
        tas(it)=ip
        tpcol(it)=col_range*rpn(ip)
		pma(ip)=pma(ip)+tma(it)
        pmas=pmas+tma(it)
        ntip(ip)=ntip(ip)+1
        enddo
    enddo
itf=it

! output north polar cap primary and secondary

write(10,1) ip,ptype(ip),nsip(ip),ntip(ip),90-pmlat(ip),pmlon(ip),pmr(ip),pma(ip),pmh(ip),ght0,pmden(ip),&
			pcland,tidemask(ip),region
write(10,2) ip,"S     ",ntip(ip),90-pmlat(ip),pmlon(ip),pmr(ip),pma(ip),pmh(ip),ght0,pmden(ip),&
			ip,sscol(ip),region
ipc=ipc+1
ntipc=ntip(ip)
pmat=pmat+pma(ip)

do it=its,itf
!   if(mod(it,100000)==0) write(*,9) it
    write(10,3) it,ttype(it),90-tmlat(it),tmlon(it),tmr(it),tma(it),tmh(it),ght0,tmden(it),&
				tap(it),tas(it),tpcol(it),tscol(it),region
    itc=itc+1
    tmat=tmat+tma(it)

    enddo

ip0=ip

!write(*,*)
!write(*,*) "primary count:",ipc
!write(*,*) "ternary count:",itc

! build band primaries

write(*,*)
write(*,*) "       pband   #primaries  #ternaries  #tern/#prim"
npsub(1)=0
do i=2,npband-1

    its=itf+1
    blat=(i-1)*dplat
    dplon=360.d0/npzone(i)
!    write(11,*) i,npzone(i),dplon
! PT210830: try to output the vertices of the primary rectangles
!    write(*,*)90-((i-1)*plat_wid_deg),dplon,' latitude, longitude spacing'

    npsub(i)=npsub(i-1)+npzone(i-1)

    do j=1,nptm
        itb=(i-1.5)*nptm+j
        do k=1,ntzone(itb)
            it=itf+k
! using 0.999 here avoids a roundoff error that can occur at tmlon=360
            ip=ip0+tmlon(it)/dplon+0.999
            if(ip.le.npsub(i)) write(*,*) "!! ",k,it,tmlon(it),dplon,tmlon(it)/dplon,ip
            if(ip.gt.npsub(i)+npzone(i)) write(*,*) "!! ",k,it,tmlon(it),dplon,tmlon(it)/dplon,ip
            tas(it)=ip
            tap(it)=ip
            tpcol(it)=col_range*rpn(ip)
			pma(ip)=pma(ip)+tma(it)
            pmas=pmas+tma(it)
            ntip(ip)=ntip(ip)+1
            enddo
        itf=it
        enddo

	call calc_geocentric(blat,Rg)
!    area=2*pi*Re**2*(dcos((i-1.5d0)*dplat*rad)-dcos((i-0.5d0)*dplat*rad))/npzone(i)
    do jp=1,npzone(i)
        ip=ip0+jp
        ptype(ip)="Ptess "
        pmlat(ip)=blat
        pmlon(ip)=(jp-0.5d0)*dplon
        pmr(ip)=Rg
        pmh(ip)=0.d0
        pmden(ip)=1000
        tidemask(ip)=0
        nsip(ip)=1
ntipc=ntip(ip)+ntipc
        enddo
if(itf.ne.ntipc) write(*,*) "!! ternary count mismatch ",itf,ntipc,itf-ntipc
write(*,*) i,npzone(i),itf-its+1,(itf-its+1)/npzone(i),' band, primaries, ternaries, tern_per_prim'
ip0=ip0+npzone(i)

   enddo


! build index arrays

do ip=1,nprm
    itp(ip)=0
    enddo
do it=1,itf
    ip=tap(it)
    if(ip==0) then
        write(*,*) "!! error: zero primary associated with ternary:",it
        stop
    else
        itp(ip)=itp(ip)+1
        itip(itp(ip),ip)=it
        endif
enddo
write(*,8) "tband:"

! output band primaries and secondaries

do ip=2,nprm-1
    if(itp(ip).ne.ntip(ip)) then
        write(*,*) ip,itp(ip),ntip(ip)
        stop "ntip .ne. itp"
        endif
    if(mod(ip,100)==0) write(*,9) ip

    write(10,1) ip,ptype(ip),nsip(ip),ntip(ip),90-pmlat(ip),pmlon(ip),pmr(ip),pma(ip),pmh(ip),ght0,pmden(ip),&
				pcland,tidemask(ip),region
    write(10,2) ip,"S     ",ntip(ip),90-pmlat(ip),pmlon(ip),pmr(ip),pma(ip),pmh(ip),ght0,pmden(ip),&
				ip,sscol(ip),region
    ipc=ipc+1
    pmat=pmat+pma(ip)
    do j=1,ntip(ip)
        it=itip(j,ip)
    write(10,3) it,ttype(it),90-tmlat(it),tmlon(it),tmr(it),tma(it),tmh(it),ght0,tmden(it),&
				tap(it),tas(it),tpcol(it),tscol(it),region
        itc=itc+1
        tmat=tmat+tma(it)
        enddo
    enddo

write(*,*)
!write(*,*) "primary count:",ipc
!write(*,*) "ternary count:",itc

! build south polar cap primary and secondary

ip=ip0+1
ptype(ip)="Ppolar"
pmlat(ip)=180.d0
pmlon(ip)=0.d0
pmr(ip)=Rb
pmh(ip)=0.d0
pmden(ip)=1000
tidemask(ip)=0
nsip(ip)=1
ntip(ip)=0
its=itf+1
j0=(npband-2)*nptm+nptm/2
do j=j0+1,j0+nptm/2+1
    do k=1,ntzone(j)
        it=ntsub(j)+k
        tas(it)=ip
        tap(it)=ip
        tpcol(it)=col_range*rpn(ip)
		pma(ip)=pma(ip)+tma(it)
        pmas=pmas+tma(it)
        ntip(ip)=ntip(ip)+1
        enddo
    enddo
itf=it


! output south polar cap primary

write(10,1) ip,ptype(ip),nsip(ip),ntip(ip),90-pmlat(ip),pmlon(ip),pmr(ip),pma(ip),pmh(ip),ght0,pmden(ip),&
				pcland,tidemask(ip),region
write(10,2) ip,"S     ",ntip(ip),90-pmlat(ip),pmlon(ip),pmr(ip),pma(ip),pmh(ip),ght0,pmden(ip),&
				ip,sscol(ip),region
ipc=ipc+1
pmat=pmat+pma(ip)
do it=its,itf
write(10,3) it,ttype(it),90-tmlat(it),tmlon(it),tmr(it),tma(it),tmh(it),ght0,tmden(it),&
				tap(it),tas(it),tpcol(it),tscol(it),region
    itc=itc+1
    tmat=tmat+tma(it)
    enddo


write(*,*)
write(*,*) "primary count:  ",ipc
write(*,*) "secondary count:",ipc
write(*,*) "ternary count  :",itc
write(*,*) "primary area:",pmat,pmas
write(*,*) "ternary area:",tmat
close(10)
close(11)

stop
end

! ******************************************

subroutine calc_lat_bands(lat_min,nband,nmasc,nzone,dlat)

! latitude band calculation on WGS84 ellipsoid 
! geocentric radius used to calculate small circle distance
! fixed bwid based on mean radius used for nzone calculation
! subroutine uses co-latitude, not geographic latitude

implicit none
real (kind=8), parameter         :: pi = 3.141592653589793d0
real (kind=8), parameter         :: Ra = 6378.137d3, inflat = 298.257223563d0
integer, parameter       :: MTB=1081

integer, dimension(1:MTB)       :: nzone
integer     :: lat_min

real (kind=8) dlat,dlon,bwid,scirc,blat
integer nband,nmasc

! working variables
real (kind=8) Rb,Re,Rg
real (kind=8) rad,sinlat,coslat,width,wmin,wmax,wmean,remf
integer i,j,k,nsub

1   format(i5," bands"," of ~",f6.1," km")
2   format("Mean mascon width:",f8.3,"km,  Min:",f8.3,"  Max:",f8.3)
9   format(i4,f15.7,2i8,2f12.2,f8.2)

rad=pi/180
Rb=Ra*(1-1/inflat)
Re=(Ra*Ra*Rb)**(1/3.d0)
nband=1+180*60/lat_min
dlat=lat_min/60.d0
bwid=Re*dlat*rad

write(*,1) nband,bwid/1000

nmasc=1
nzone(1)=1
nsub=1
wmax=0.d0
wmin=2*pi*Re
wmean=0.d0
do i=2,nband-1
    blat=(i-1)*dlat
    coslat=dcos(blat*rad)
    sinlat=dsin(blat*rad)
    Rg=(Ra**2*sinlat)**2+(Rb**2*coslat)**2
    Rg=sqrt(Rg/((Ra*sinlat)**2+(Rb*coslat)**2))
    scirc=2*pi*Rg*dsin(blat*rad)
    nzone(i)=nint(scirc/bwid)
    dlon=360.d0/nzone(i)
    nmasc=nmasc+nzone(i)
    width=scirc/nzone(i)
    wmax=max(width,wmax)
    wmin=min(width,wmin)
    wmean=wmean+scirc

    if(i.le.11) then
        width=scirc/nzone(i)
        remf=scirc/bwid-nzone(i)
        nsub=nsub+nzone(i)
        write(11,9) i,blat,nzone(i),nsub,scirc/1000,width/1000,remf
        endif
    enddo

nzone(nband)=1
nmasc=nmasc+1

wmean=wmean/(nmasc-2)
write(*,2) wmean/1000,wmin/1000,wmax/1000

return
end

! ******************************************

subroutine calc_geocentric(lat,Rg)

! geocentric radius at given latitude on WGS84 ellipsoid
!   input lat is co-latitude

implicit none
real (kind=8), parameter         :: pi = 3.141592653589793d0
real (kind=8), parameter         :: Ra = 6378.137d3, inflat = 298.257223563d0

real (kind=8), intent(in)  :: lat
real (kind=8), intent(out) :: Rg
real (kind=8) Rb,sinlat,coslat
real (kind=8) rad

rad=pi/180
Rb=Ra*(1-1/inflat)
coslat=dcos(lat*rad)
sinlat=dsin(lat*rad)
Rg=(Ra**2*sinlat)**2+(Rb**2*coslat)**2
Rg=sqrt(Rg/((Ra*sinlat)**2+(Rb*coslat)**2))

return
end

