      program classify_mascons
!     *******
!
!  Reads and rewrites mascon files setting mascon identities:
!     altitude, density, geoid height and region identifier
!
!   - altitudes from GEBCO14+Bedmap2
!   - density 1029 for open ocean, 1000 elsewhere
!	- geoid from
!

    implicit none
	integer, parameter      :: MP = 50000		! max primaries
	integer, parameter      :: MS = 50000		! max secondaries
	integer, parameter      :: MT = 1500000		! max ternaries
	integer, parameter      :: MPB = 181		! max primary bands
	integer, parameter      :: MSB = 181		! max secondary bands
	integer, parameter      :: MTB = 1081		! max ternary bands

    integer, parameter      :: MM = 1800000
    integer, parameter      :: MH=10

    real (kind=8), dimension(1:MP)   :: pmlat,pmlon,pmr,pma,pmh,pmg,pmden	! primary lat, lon, radius, area, geoid height & density
    real (kind=8), dimension(1:MS)   :: smlat,smlon,smr,sma,smh,smg,smden	! secondary lat, lon, radius, area, geoid height & density
    real (kind=8), dimension(1:MT)   :: tmlat,tmlon,tmr,tma,tmh,tmg,tmden	! ternary lat, lon, radius, area, geoid height & density
    integer, dimension(1:MP)         :: nsip,ntip,ntlip,tidemask		! #secondaries in primary, #ternaries in primary, tide mask
	integer, dimension(1:MS)         :: ntis,sap,sscol					! #ternaries in secondary, secondary's primary, secondary's colour
	integer, dimension(1:MT)         :: tap,tas					! ternary's primary, ternary's secondary
    integer, dimension(1:MT)         :: tpcol,tscol				! ternary's primary colour, ternary's secondary colour

	character(6),	dimension(1:MP)		:: ptype	! primary mascon type
	character(6),	dimension(1:MS)		:: stype	! secondary mascon type
	character(6),	dimension(1:MT)		:: ttype	! ternary mascon type
    character(1),   dimension(1:MM)     :: mlevel	! mascon level code (P/S/T)

    character(150), dimension(1:MH)     :: headrec	! header records
    integer npex,nsex,ntex						! expected number of primary, secondary and ternary mascons (from file header)
    integer maxsip,maxtip,maxtis				! max# of secondaries in primary, ternaries in primary, ternaries in secondary (from header)
	character(12)	:: hashcode						! hashcode identifying mascon file
	character(12)	:: region						! mascon region

! lat-lon topo grid (from GEBCO14+Bedmap2 netcdf file)

	integer, parameter :: NLATS = 21600, NLONS = 43200
	real, dimension (:),allocatable  :: lats,lons
	real, dimension (:,:),allocatable  :: ht_in

! working variables

	integer, dimension(1:MP)         :: npmh				! count of ternary points in primary for height average
	integer, dimension(1:MS)         :: nsmh				! count of ternary points in secondary for height average
	integer, dimension(1:MM)         :: index				! mascon number index

    real (kind=8) area_p_total,area_s_total,area_t_total
    real (kind=8) oceanp,pcland
    integer np,ns,nt,ip,is,it
    integer nsipt,ntipt,ntist
	integer nmascons,mhr,lun
    integer nsea,nfr,ndry,nhvy,nderr
	integer i,j
	character(1)	:: level
	character(180)	:: line
	character(8)	:: date
    character(10)	:: time
    character(5)	:: timezone

! Formats

 1   format(i7,2x,a6,2i8,2f10.4,f11.1,f18.0,2f9.1,f6.0,f6.1,i6,2x,a20)
 2   format(i7,2x,a6,i8,8x,2f10.4,f11.1,f18.0,2f9.1,f6.0,2i6,2x,a20)
 3   format(i7,2x,a6,2f10.4,f11.1,f15.0,2f9.1,f6.0,4i6,12x,a12)
 4   format(a180)
 5   format(9x,a1)
 6   format(a12,2i8,i10,3i7)
 7   format(a,1x,a8,1x,a10,1x,a5)
 9   format(i8,$)


! setup

    np=0
    ns=0
    nt=0
    npex=1607
    nsipt=0
    ntipt=0
    ntist=0
    area_p_total=0.
    area_s_total=0.
    area_t_total=0.
    i=0
    j=0
    ip=0
    is=0
    it=0

! initialize primary and secondary height and geoid 
!  - values will be set from averages of contained ternaries
	do i=1,MP
		pmh(i)=0.d0
		pmg(i)=0.d0
        ntlip(i)=0
		enddo
	do i=1,MS
		smh(i)=0.d0
		smg(i)=0.d0
		enddo

!

    open(10,file='mascons_stage1',status='old')
    open(11,file='mascons_stage2',status='unknown')
    write(*,*) "input file mascons_stage1"

! check header and storage parameters

    lun=10
    call read_header(lun,headrec,MH,hashcode,npex,nsex,ntex)
    if(npex.gt.MP) stop "error: parameter MP too small"
    if(nsex.gt.MS) stop "error: parameter MS too small"
    if(ntex.gt.MT) stop "error: parameter MS too small"

! read mascon file
	write(*,*) "reading mascon file..."
 110  read(10,4,end=120) line

    j=j+1
    read(line,5) level
    if(level=="P")then
        read(line,*) i,ptype(i),nsip(i),ntip(i),pmlat(i),pmlon(i),pmr(i),pma(i),pmh(i),pmg(i),pmden(i),&
					pcland,tidemask(i),region
        mlevel(j)=level
        index(j)=i
        ip=ip+1
		np=max(i,np)

    else if(level=="S") then
        read(line,*) i,stype(i),ntis(i),smlat(i),smlon(i),smr(i),sma(i),smh(i),smg(i),smden(i),&
					sap(i),sscol(i),region
        mlevel(j)=level
        index(j)=i
        is=is+1
        ns=max(i,ns)

    else if(level=="T") then
        read(line,*) i,ttype(i),tmlat(i),tmlon(i),tmr(i),tma(i),tmh(i),tmg(i),tmden(i),&
					tap(i),tas(i),tpcol(i),tscol(i),region
        mlevel(j)=level
        index(j)=i
        it=it+1
        nt=max(i,nt)
!        if(nt.eq.1485118) go to 120

    else
        stop "mascon type error"
        endif
    if(mod(j,100000)==0) write(*,9) j
    go to 110

 120  close(10)

    nmascons=j
    write(*,*)
    write(*,*) nmascons," mascons read"
	if(nt.le.0) then
		stop " Error: No ternary mascons - terminating"
		endif

! get geoid heights on ternary mascons

	write(*,*) "setting ternary mascon geoid heights"
	call set_Tmascon_geoid(nt,tmlat,tmlon,tmg)


! read in GEBCO/Bedmap height grid

	allocate(lats(NLATS))
	allocate(lons(NLONS))
	allocate(ht_in(NLONS,NLATS))
	call read_topo_grid(lats,lons,ht_in)

! characterise mascons


	write(*,*) "setting ternary mascons"
	call set_Tmascon_identity(nt,ttype,tmlat,tmlon,tmh,tmden,lats,lons,ht_in)
	do i=1,nt
	if(ttype(i)(1:4)=="land") ntlip(tap(i))=ntlip(tap(i))+1
		enddo

	write(*,*) "setting primary mascons"

	do j=1,nmascons
		i=index(j)
		if(mlevel(j)=="T") then
			npmh(tap(i))=npmh(tap(i))+1
			nsmh(tas(i))=nsmh(tas(i))+1
			pmh(tap(i))=pmh(tap(i))+tmh(i)
			smh(tas(i))=smh(tas(i))+tmh(i)
			pmg(tap(i))=pmg(tap(i))+tmg(i)
			smg(tas(i))=smg(tas(i))+tmg(i)
			endif

!		if(mod(j,200000)==0) write(*,9) j
		enddo
	do i=1,np
		if(npmh(i).ne.ntip(i)) then
			write(*,*) i,npmh(i),ntip(i)," ternary in primary mismatch error"
			stop
			endif
		enddo
	do i=1,ns
		if(nsmh(i).ne.ntis(i)) then
			write(*,*) i,nsmh(i),ntis(i)," ternary in primary mismatch error"
			stop
			endif
		enddo


    call set_PSmascon_identity(np,ptype,pmlat,pmlon,pmh,pmg,pmden,npmh)
    maxsip=0
    maxtip=0
    do i=1,np
        maxsip=max(nsip(i),maxsip)
        maxtip=max(ntip(i),maxtip)
        enddo

    maxtis=0
    if(ns.gt.0) then
        write(*,*) "setting secondary mascons"
        call set_PSmascon_identity(ns,stype,smlat,smlon,smh,smg,smden,nsmh)
        do i=1,ns
            maxtis=max(ntis(i),maxtis)
            enddo
        endif

	deallocate(lats)
	deallocate(lons)
	deallocate(ht_in)



! write updated mascon file

	write(*,*) "writing mascon file..."
	write(11,6) hashcode,np,ns,nt,maxsip,maxtip,maxtis
	call date_and_time(date,time,timezone)
	write(11,'(a)') trim(headrec(1))
	write(11,7) "# modified by classify_mascons ",date,time,timezone


	do j=1,nmascons
		i=index(j)
		if(mlevel(j)=="P")then
			ptype(i)="P"//ptype(i)
			call set_region(pmlat(i),pmlon(i),pmh(i),region)
			if(ptype(i).eq."Pdeep".or.ptype(i).eq."Pshelf") then
				tidemask(i)=31
			else
				tidemask(i)=0
				endif
			pcland=ntlip(i)*100./ntip(i)
			write(11,1) i,ptype(i),nsip(i),ntip(i),pmlat(i),pmlon(i),pmr(i),pma(i),pmh(i),pmg(i),pmden(i),&
						pcland,tidemask(i),region

!        call random_number(rnum)
!        sscol(i)=rnum*1600

		else if(mlevel(j)=="S") then
			stype(i)="S"//stype(i)
			call set_region(smlat(i),smlon(i),smh(i),region)
			write(11,2) i,stype(i),ntis(i),smlat(i),smlon(i),smr(i),sma(i),smh(i),smg(i),smden(i),&
						sap(i),sscol(i),region

		else if(mlevel(j)=="T") then
			ttype(i)="T"//ttype(i)
			call set_region(tmlat(i),tmlon(i),tmh(i),region)
			write(11,3) i,ttype(i),tmlat(i),tmlon(i),tmr(i),tma(i),tmh(i),tmg(i),tmden(i),&
						tap(i),tas(i),tpcol(i),tscol(i),region
			endif

		if(mod(j,100000)==0) write(*,9) j
		enddo

	close(11)

    write(*,*)
    write(*,*) "completeness check:"
    write(*,*) " mascon count (P,S,T) ",ip,is,it
    write(*,*) " mascon range (P,S,T) ",np,ns,nt
    write(*,*) " total mascon lines found   :",nmascons
    if(np.ne.ip) write(*,*) "!! primary count differs from maximum index by  ",np-ip
    if(ns.ne.is) write(*,*) "!! secondary count differs from maximum index by",ns-is
    if(nt.ne.it) write(*,*) "!! ternary count differs from maximum index by  ",nt-it

!    call status_update('STATUS','UTIL','mascon_rewrite',"mascon_rw","End program",0)

    stop
    end program

!************************************************************

subroutine read_header(lun,headrec,mhr,hashcode,npex,nsex,ntex)

!
!   reads header records and verifies hashcode
!

    implicit none

    integer lun,i,nhr,mhr
    integer npex,nsex,ntex
    character*12 hashcode,checkcode
    character(150), dimension(1:mhr)     :: headrec

    i=0
    rewind(lun)
    read(lun,*) hashcode,npex,nsex,ntex
    if(hashcode(1:1).ne."#") stop "error: '#' not found on header line 1"

100 i=i+1
    if(i.gt.mhr) go to 910
    read(lun,'(a120)',end=900) headrec(i)
    if(headrec(i)(1:1).eq."#") go to 100

    nhr=i-1
    if(nhr==0) stop "error: no '#' header lines found"

    if(hashcode=="#00000") then
        write(*,*) "dummy hashcode on input mascon file ",hashcode
    else
        call hash_header(headrec,nhr,checkcode)
        if(hashcode.ne.checkcode) then
            write(*,*) "hashcode mismatch in file header: ",checkcode," != ",hashcode
            stop "  aborting"
            endif
        endif

    backspace(lun)
    return

900 stop "error: read_header - end of file"
910 stop "error: read_header - header too long, increase MHR"

    end

!************************************************************

    subroutine hash_header(headrec,nhr,hashcode)

! generate a hashcode from the current header lines

    implicit none

    integer i,nhr
    character*12 hashcode
    character(150), dimension(1:nhr)    :: headrec

    open(99,file="headcheck",status="unknown")
    write(99,'(a120)') (trim(headrec(i)),i=1,nhr)
    close(99)

    call system("cksum headcheck|awk '{print $1}' >checksum")
    open(99,file="checksum",status="old")
    read(99,'(a12)') hashcode
    close(99)
    hashcode="#"//hashcode

    return
    end


!************************************************************

subroutine set_region(lat,lon,ht,region)


!
!   identifies terrestrial region a point resides in

!

	implicit none
    real (kind=8), parameter         :: pi = 3.141592653589793d0
	integer, parameter      :: MP=19
	real, dimension(1:MP)  :: alat,alon
	character (len=12), dimension(1:MP) :: aloc
	real (kind=8)       :: lat,lon,ht,rad,dist,dnear
	integer             :: nap,i,inear
	character (len=12)   :: region
	data nap/19/

!	data alat/90.,-90./
!	data alon/0.,0./
!	data aloc/"north","south"/

	data alat/73.,64.,33.,18.,0.,-39.,-77.,-75.,-25.,-30.,7.,5.,39.,47.,52.,22.,5.,-23.,10./
	data alon/318.,268.,260.,272.,290.,293.,287.,109.,124.,149.,133.,106.,103.,56.,6.,3.,25.,27.,190./
	data aloc/"Greenland","NAmerica","NAmerica","CAmerica","SAmerica","SAmerica","Antarctic","Antarctic",&
	"Australasia","Australasia","SEAsia","SEAsia","Asia","Asia","Europe","Africa","Africa","Africa","Oceania"/

!	open(12,file='region_anchors',status='old')
!	read(12,*) nap
!	do i=1,nap
!		read(*,*) alat(i),alon(i),aloc(i)
!		enddo
!	close(12)

    rad=pi/180.d0

	if(ht.lt.-200) then
		region="Ocean"
	else if(ht.lt.0) then
		region="Shelf"
	else
		dnear=4
		do i=1,nap
			dist=sin(lat*rad)*sin(alat(i)*rad)
			dist=dist+cos(lat*rad)*cos(alat(i)*rad)*cos((lon-alon(i))*rad)
			dist=acos(dist)
			if(dist.lt.dnear) then
				dnear=dist
				inear=i
				endif
			enddo
		region=aloc(inear)
		endif

!	write(*,*) inear,region
	return
	end subroutine set_region

!************************************************************

subroutine read_topo_grid(lats,lons,ht_in)


!
!   Reads topo surface from netCDF file
!   - data on 1 minute lat-lon grid
!

	use netcdf
	implicit none

! data file
	character (len = *), parameter :: FILE_NAME = "GEBCO14+Bedmap2_topo.nc"
	integer :: ncid

! lat-lon topo grid
!(old grid) integer, parameter :: NLATS = 10801, NLONS = 21601
	integer, parameter :: NLATS = 21600, NLONS = 43200
	real,intent(out) :: lats(NLATS), lons(NLONS)
	real,intent(out) :: ht_in(NLONS,NLATS)


	character (len = *), parameter :: LAT_NAME = "lat"
	character (len = *), parameter :: LON_NAME = "lon"
	integer :: lat_dimid, lon_dimid
	integer :: lat_varid, lon_varid

! surface height
	character (len = *), parameter :: HT_NAME = "z"
	integer :: ht_varid
!integer :: dimids(NDIMS)

! To check the units attributes
	character (len = *), parameter :: UNITS = "units"
	character (len = *), parameter :: HT_UNITS = "m"
	character (len = *), parameter :: LAT_UNITS = "degrees_north"
	character (len = *), parameter :: LON_UNITS = "degrees_east"
	integer, parameter :: MAX_ATT_LEN = 80
	integer :: att_len
	character*(MAX_ATT_LEN) :: ht_units_in
	character*(MAX_ATT_LEN) :: lat_units_in, lon_units_in

! working variables.
	integer :: ndims_in,nvars_in,ngatts_in,unlimdimid_in,ilat
	real :: tlat,tlon,tht,dlat,dlon
	character (len=20) :: location

	print *,"reading topo data from ",FILE_NAME
! Open the file and read
!	the number of netCDF variables, dimensions, and global attributes
!	and the dimension id of the unlimited dimension, if there is one.
! Expect 2 netCDF dimensions, 3 netCDF variables, 3 global attributes, 
!	and no unlimited dimension.
	call check( nf90_open(FILE_NAME, nf90_nowrite, ncid) )
	call check( nf90_inquire(ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in) )
!	write(*,*) ndims_in,nvars_in,ngatts_in,unlimdimid_in


! Read the latitude and longitude data.
	call check( nf90_inq_varid(ncid, "lat", lat_varid) )
	call check( nf90_inq_varid(ncid, "lon", lon_varid) )
	call check( nf90_get_var(ncid, lat_varid, lats) )
	call check( nf90_get_var(ncid, lon_varid, lons) )
!	write(*,*) ncid," y"," x",lat_varid,lon_varid

! Read the surface height data from the file.
	call check( nf90_inq_varid(ncid, HT_NAME, ht_varid) )
	call check( nf90_get_var(ncid, ht_varid, ht_in) )


! check variables units
	call check( nf90_get_att(ncid, lat_varid, UNITS, lat_units_in) )
	call check( nf90_inquire_attribute(ncid, lat_varid, UNITS, len = att_len) )
	if (lat_units_in(1:att_len) /= LAT_UNITS) stop "lat_units check in read_topo_grid"
	call check( nf90_get_att(ncid, lon_varid, UNITS, lon_units_in) )
	call check( nf90_inquire_attribute(ncid, lon_varid, UNITS, len = att_len) )
	if (lon_units_in(1:att_len) /= LON_UNITS) stop "lon_units check in read_topo_grid"

! Close the file, freeing up internal netCDF resources
	call check( nf90_close(ncid) )
	write(*,*) " ",FILE_NAME," topo data read complete"

	dlat=lats(2)-lats(1)
	dlon=lons(2)-lons(1)
	write(*,*) "lat range and spacing ",lats(1),lats(NLATS),dlat
	write(*,*) "lon range and spacing ",lons(1),lons(NLONS),dlon
	ilat=1+120*(90+51.477)
	write(*,*) "Greenwich altitude (46m)",ht_in(1,ilat)
	write(*,*)

	return

!**********************************

contains
	subroutine check(status)
	integer, intent ( in) :: status

	if(status /= nf90_noerr) then
		print *, trim(nf90_strerror(status))
		stop "error exit"
		end if
	end subroutine check
end subroutine read_topo_grid


!************************************************************

	subroutine set_Tmascon_identity(nm,mtype,mlat,mlon,mh,mden,lats,lons,ht_in)

!
!  sets Ternary mascon parameters representative height, type and density
!  
!	in purely land or ocean mascons the representative height is the average height
!	but in mixed mascons it is the average of the majority type (land or ocean) only

	implicit none

! formal arguments
	real (kind=8), dimension(1:nm)   :: mlat,mlon,mh,mden
	character (len=6), dimension(1:nm)  :: mtype
	integer :: nm

! lat-lon topo grid
	integer, parameter :: NLATS = 21600, NLONS = 43200
	real,intent(in) :: lats(NLATS), lons(NLONS)
	real,intent(in) :: ht_in(NLONS,NLATS)

! working variables.
    real (kind=8), parameter  :: pi = 3.141592653589793d0
	real (kind=8):: rad,tlat,tlon,tht,dlat,dlon,dmlat
	real (kind=8):: gratio,ht,sum_o,sum_l,f
	integer :: lat,lon,ilat,ilon,i,ntpm,ipole,nl,no,navp
	integer :: j,j1,j2,k,k1,k2,k3,k4,dilon,sumdi
	logical :: debug

9    format(i8,$)

    rad=pi/180.d0

! mascon grid spacing
	dmlat=abs(mlat(2)-mlat(1))

! topo grid spacing
	dlat=lats(2)-lats(1)
	dlon=lons(2)-lons(1)

! grid spacing ratio
	ntpm=dmlat/dlat
!	write(*,*) dlon,dlat,dmlat,ntpm
	open(99,file="debug",status="unknown")


! polar mascons

	do ipole=1,2

		if(ipole.eq.1) then
			i=1
			ilat=NLATS
			ilon=1
			j1=ilat-ntpm/2
			j2=ilat
		else
			i=nm
			ilat=1
			ilon=1
			j1=ilat
			j2=ilat+ntpm/2
			endif

!
! sum land and ocean areas and set mascon heights
		sum_o=0.d0
		sum_l=0.d0
		no=0
		nl=0
!		write(*,*)"A ",i,j1,j2,mlat(i),mlon(i),ilat,ilon
		do j=j1,j2
		do k=1,NLONS
			ht=ht_in(k,j)
			if(ht.le.0) then
				sum_o=sum_o+ht
				no=no+1
			else
				sum_l=sum_l+ht
				nl=nl+1
				endif
			enddo
			enddo
		if(no+nl-(j2-j1+1)*NLONS.ne.0) write(*,*) "miscount",no,nl,j1,j2,NLONS

		if(no.gt.nl) then
			tht=sum_o/no
		else
			tht=sum_l/nl
			endif


! set mascon types and densities
!		write(*,*) i,mlat(i),mlon(i),ilon,ilat,no,nl,tht
		mh(i)=tht
		if(tht.lt.-200.) then
			mtype(i)="deep"
			mden(i)=1029
		else if(tht.lt.0.) then
			mtype(i)="shelf"
			mden(i)=1029
		else
			mtype(i)="land"
			mden(i)=1000
			endif
		enddo


! loop over tesseral mascons

	sumdi=0
	do i=2,nm-1

! find nearest topo grid point to mascon centre
!   expecting mlat in (0,90) and mlon in (0,360)
		if(mod(i,100000)==0) write(*,9) i
		if(mlat(i).ge.90) then
			ilat=NLATS
		else
			ilat=int((mlat(i)+90)*120+1)
		endif
		if(mlon(i).ge.360) then
			ilon=NLONS
		else
			ilon=int(mlon(i)*120+1)
		endif

! sanity check: interpolated grid point is within half topo grid spacing of target
		if((lats(ilat)-mlat(i)).gt.dlat*0.51) then
			write(*,*) "latitude misfit:",i,mlat(i),ilat,lats(ilat),lats(ilat-1)
			stop
			endif
		if((lons(ilon)-mlon(i)).gt.dlon*0.51) then
			write(*,*) "longitude misfit:",i,mlon(i),ilon,lons(ilon),lons(ilon-1)
			stop
			endif

! dump height points used for a particular ternary
!debug=.false.
!tlat=-40.
!tlon=147.
!if(abs(mlat(i)-tlat).lt.0.1.and.abs(mlon(i)-tlon).lt.0.1) debug=.true.

!
! sum land and ocean areas and set mascon heights
		sum_o=0.d0
		sum_l=0.d0
		no=0
		nl=0
		j1=ilat-ntpm/2
		j2=ilat+ntpm/2
		gratio=1/cos(mlat(i)*rad)
		dilon=int(gratio*ntpm/2)
		f=NLONS/(2.*dilon)
		dilon=nint(f*dilon/nint(f))
		sumdi=sumdi+dilon*2
		k1=ilon-dilon
		k2=ilon
		k3=ilon+1
		k4=ilon+dilon
		if(k1.le.0) then
			k1=k1+NLONS
			k2=NLONS
			k3=1
		else if(k4.gt.NLONS) then
			k2=NLONS
			k3=1
			k4=k4-NLONS
			endif

!if(debug) write(99,*) "A ",i,j1,j2,k1,k2,k3,k4
		do j=j1,j2
		do k=k1,k2
			ht=ht_in(k,j)
			if(ht.le.0) then
				sum_o=sum_o+ht
				no=no+1
			else
				sum_l=sum_l+ht
				nl=nl+1
				endif
!if(debug) write(99,*) no,ht,sum_o
			enddo
		do k=k3,k4
			ht=ht_in(k,j)
			if(ht.le.0) then
				sum_o=sum_o+ht
				no=no+1
			else
				sum_l=sum_l+ht
				nl=nl+1
				endif
!if(debug) write(99,*) no,ht,sum_o
			enddo
			enddo

		navp=(j2-j1+1)*((k2-k1+1)+(k4-k3+1))
		if(no+nl-navp.ne.0) then
			write(*,*) "miscount",i,no,nl,j1,j2,k1,k2,k3,k4
			stop "error exit"
			endif

		if(no.gt.nl) then
			tht=sum_o/no
		else
			tht=sum_l/nl
			endif

! set mascon types and densities
!		if(i.le.30) write(*,*) i,no,nl,tht
		mh(i)=tht
		if(tht.lt.-200.) then
			mtype(i)="deep"
			mden(i)=1029
		else if(tht.lt.0.) then
			mtype(i)="shelf"
			mden(i)=1029
		else
			mtype(i)="land"
			mden(i)=1000
			endif

!if(debug) write(99,*) i,mlat(i),mlon(i),mh(i)

		enddo

	close(99)
	if(nm.gt.1.e5) write(*,*)
	write(*,*) "overlaps per line ",(sumdi/1080)-43200
	write(*,*) nm," points classified"
	return

	end subroutine set_Tmascon_identity



!************************************************************


	subroutine set_PSmascon_identity(nm,mtype,mlat,mlon,mh,mg,mden,ntim)

!
!  sets Primary/Secondary mascon parameters representative height, geoid height, type and density
!	height and geoid are averages of contained ternaries
!	type and density are defined by the height

	implicit none

! formal arguments
	real (kind=8), dimension(1:nm)   :: mlat,mlon,mh,mg,mden
	integer, dimension(1:nm) :: ntim
	character (len=6), dimension(1:nm)  :: mtype
	integer :: nm

! working variables.
	real (kind=8):: tht
	integer :: i


	do i=1,nm

		mh(i)=mh(i)/ntim(i)
		mg(i)=mg(i)/ntim(i)

! set mascon types and densities
!		if(i.le.3) write(*,*) i,mlat(i),mlon(i),mh(i),mg(i),ntim(i)

		if(mh(i).lt.-200.) then
			mtype(i)="deep"
			mden(i)=1029
		else if(mh(i).lt.0.) then
			mtype(i)="shelf"
			mden(i)=1029
		else
			mtype(i)="land"
			mden(i)=1000
			endif

		enddo

	if(nm.gt.1.e5) write(*,*)
	write(*,*) nm," points classified"
	return

	end subroutine set_PSmascon_identity

!************************************************************

	subroutine set_Tmascon_geoid(nm,mlat,mlon,mg)

!
!  sets Ternary mascon geoid height from EGM96 WGS84 geoid
!

	implicit none

! formal arguments
	real (kind=8), dimension(1:nm)   :: mlat,mlon,mg
	integer :: nm

! lat-lon geoid grid
	integer, parameter :: NGLATS = 2160, NGLONS = 4320

! working variables
	real (kind=8), parameter  :: pi = 3.141592653589793d0
	real (kind=8):: rad,tlat,tlon,ght
	integer :: i

9   format(i8,$)

	open(13,file="GEOID_points.dat",status="unknown")
	do i=1,nm
		write(13,*) mlat(i),mlon(i)
		enddo
	close(13)

	write(*,*) "calling interp_2p5m to interpolate EGU2008 2.5min geoid",&
		" onto ternary mascons"
	write(*,*) "requires local link to Und_min2.5x2.5_egm2008_isw=82_WGS84_TideFree_SE"

	call interp_2p5m

! previously interp_2p5m was run as a separate program
!	call execute_command_line ("interp_2p5m", exitstat=i)
!	write(*,*) "Exit status of interp_2p5m was ", i

	open(13,file="GEOID_heights.dat",status="old")
	do i=1,nm
		read(13,*) tlat,tlon,ght
		if((tlat-mlat(i)).le.0.1.and.(tlon-mlon(i)).le.0.1) then
			mg(i)=ght
		else
			stop "geoid point mismatch"
			endif
		enddo
	close(13)


	return

	end subroutine set_Tmascon_geoid

!************************************************************

