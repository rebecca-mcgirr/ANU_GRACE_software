      program merge_ternary
!     *******
!
!  Reads and rewrites a mascon file after moving a list of ternaries into different primary mascons

!
!  usage: merge_ternary <input mascon_file> <movement_list_file> [-x]
!	   - movements are of the form "ternary# new_primary#"
!	   - if <movement_list_file> is empty the mascon file is just reordered sequentially by primary mascon
! 	   -x : interactive increase to maxtip,maxsip
!
!  output to <input_file>_tmerge
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


!    use mascon_setup

    implicit none

!-----------------
    
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

    character*15, dimension(:), allocatable	   :: pregion
    character*15, dimension(:), allocatable	   :: sregion
    character*15, dimension(:), allocatable	   :: tregion
    character(6), dimension(:), allocatable    :: ptype
    character(6), dimension(:), allocatable    :: stype
    character(6), dimension(:), allocatable    :: ttype

    character(150), dimension(1:MHR)  :: headrec

!-----------------

    integer, dimension(:,:), allocatable 	:: lsip,ltip,ltis

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

    real (kind=8) sum1,sum2,sum3,sum4,sum5,a1,a2,dif,difmin
    integer i,j,k,is1,is2,it1,ioerr
    integer offset,ipo,iso,ipn,isn,itl,itn1,narg,ntin,ntex,nmv
    character*120 infile,mlfile,outfile
    character*14 ext
    character*6 arg2,arg3,arg4,arg
    logical expand


1   format(i7,2x,a6,2i8,2f10.4,f11.1,f18.0,2f9.1,f6.0,f6.1,i6,2x,a15)
2   format(i7,2x,a6,i8,8x,2f10.4,f11.1,f18.0,2f9.1,f6.0,2i6,2x,a15)
3   format(i7,2x,a6,2f10.4,f11.1,f15.0,2f9.1,f6.0,4i6,4x,a15)

4   format(a150)
5   format(9x,a1)
6   format(a12,2i8,i10,3i7)
7   format(9x,i4)
9   format(i8,$)

! parse command line

    narg=iargc()
    if(narg.lt.2) then
       write(*,*)
       write(*,*) "usage: merge_tern <input mascon_file> <movement_list_file> [-x]"
       write(*,*) " - movements are of the form: ternary#  new_primary# "
       write(*,*) " - if <movement_list_file> is empty the mascon file is just reordered sequentially by primary mascon"
       write(*,*) " -x : interactive increase to maxtip,maxsip"
       write(*,*)
       stop
       endif

    call getarg(1,infile)
    write(*,*) "Input mascon file:  ",infile

    call getarg(2,mlfile)
    write(*,*) "Movement list file:  ",mlfile

    expand=.false.

    if(narg.gt.2) then
       call getarg(3,arg)
       if(arg.eq."-x") then
          expand=.true.
          endif
       endif

    ext="_tmerge"
    outfile=trim(infile)//trim(ext)
    write(*,*) "Output mascon file: ",outfile
    write(*,*) "Note: primary mascons will be sorted on output"

! open mascon files

    open(10,file=infile,status="old")
    write(*,*) "Header line 1:"

! Allocate arrays

    read(10,*) hashcode,MP,MS,MT,maxsip,maxtip,maxtis
    write(*,*) hashcode,MP,MS,MT,maxsip,maxtip,maxtis
    oldmax=maxtip
    if(expand) then
! PT190712: read this off the command line rather than it being interactive
       call getarg(4,arg)
       read(arg,*)newmax
       write(*,"(a40,i5,a10, i5)") "Expanding max ternaries in primary from ",maxtip," to ",newmax
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
    allocate(ltis(MS,maxtis))
    

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
!	write(*,*) "secondary ",is," not assocated with its primary ",ip
             nc_sap=nc_sap+1
             endif

          do k=1,ntis(is)
             read(10,3) it,ttype(it),tmlat(it),tmlon(it),tmr(it),tma(it),tmh(it),tmg(it),tmden(it),tap(it),tas(it),&
						tpcol(it),tscol(it),tregion(it)

             nm=nm+1
             nt=max(it,nt)
             ltis(is,k)=it
!	     ltip(ip,k)=it	! this assumes ony one secondary per primary

! check P/S associations of each ternary
             if(tap(it).ne.ip) then
!		write(*,*) "ternary ",it," not assocated with its primary ",ip
                nc_tap=nc_tap+1
                endif
             if(tas(it).ne.is) then
!		write(*,*) "ternary ",it," not assocated with its secondary ",is
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

! main ternary shuffling loop

       
    open(12,file=mlfile,status="old")
    nmv=0
    do while (ioerr == 0)
       read(12,*,iostat=ioerr,end=200) it,ipn
       nmv=nmv+1
       ipo=tap(it)

! verify mascon compatibility for move
       write(*,*) "ipn", ipn
       write(*,*) "ipo", ipo
       if(ptype(ipn).ne.ptype(ipo)) then
          write(*,*) "! Warning: mascon types don't match",it,ttype(it),ipn,ptype(ipn)
          !stop
       endif

       if(pregion(ipn).ne.pregion(ipo)) then
          write(*,*) "! Warning: mascon regions don't match",it,ttype(it),ipn,ptype(ipn)
          !stop
       endif

! check that enlarged mascon is not too big for arrays

       nn=ntip(ipn)+1
       if(nn.gt.maxtis) then
          write(*,*) "Error: enlarged mascon is too big for arrays ( #ternaries > maxtis )"
          write(*,'(4i9)') it,ipn,": ",nn," > ",maxtis
          stop
       endif

! move ternary to new primary

       iso=lsip(ipo,1)
       isn=lsip(ipn,1)
       tap(it)=ipn
       tas(it)=isn
       write(*,*) "*move*",it,iso,isn
       ntip(ipo)=ntip(ipo)-1
       ntis(iso)=ntis(iso)-1
       ntip(ipn)=ntip(ipn)+1
       ntis(isn)=ntis(isn)+1
       maxtip=max(maxtip,ntip(ipn))
       maxtis=maxtip

! update ternary parameters to match first ternary in new primary (itn1)

       itn1=ltis(isn,1)
       tpcol(it)=tpcol(itn1)
       tscol(it)=tscol(itn1)
       tregion(it)=tregion(itn1)
       
! edit ternary lookup tables
       
       offset=0
       do itl=1,ntis(iso)
          if(ltis(iso,itl).eq.it) offset=1
          ltis(iso,itl)=ltis(iso,itl+offset)
       enddo
       ltis(isn,ntis(isn))=it      
       
    enddo

200 write(*,*) nmv,"ternaries moved"
       


!-----------------------------------------------------

! write modified mascon file
!  - P/S mascon membership is updated on output

      open(11,file=outfile,status='unknown') 
  
      if(maxtip.gt.oldmax) then
         write(*,*) "resetting maxtip,maxtis to ",maxtip
         endif

      write(headrec(1),6) hashcode,MP,MS,MT,maxsip,maxtip,maxtis
      nh=nh+1
      if(nh.gt.MHR) stop "error: exceeded maximum header lines - increase MHR"
      write(headrec(nh),'(a,i5,a)') "# merge_ternary          :",nmv," ternary mascons moved to more suitable primaries"

!      write(*,*) hashcode
!      call hash_header(headrec,nh,hashcode)
!      write(*,*) hashcode

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
    write(*,*)
    write(*,*) "new mascon count (P,S,T): ",MP,MS,MT

!    call status_update('STATUS','UTIL','merge_mascons',outfile,"End program",0)

    stop

920 stop "mascon file error: no mascons?"
    end


!--------------------------------------------------------------------------------

!!    module mascon_setup
       
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
! nc_sap,nc_tap,nc_tas - number of corrected mascon connections made

!    real (kind=8), parameter    :: pi = 3.141592653589793d0
!    integer, parameter      :: MHR = 200

!    real (kind=8), dimension(:), allocatable    :: pmlat,pmlon,pmr,pma,pmh,pmg,pmden,pcland
!    real (kind=8), dimension(:), allocatable    :: smlat,smlon,smr,sma,smh,smg,smden
!    real (kind=8), dimension(:), allocatable    :: tmlat,tmlon,tmr,tma,tmh,tmg,tmden
!    real (kind=8), dimension(:), allocatable    :: dspmax,dtpmax,dtsmax
!    integer, dimension(:), allocatable          :: nsip,ntip,tidemask
!    integer, dimension(:), allocatable          :: ntis,sap,sscol
!    integer, dimension(:), allocatable          :: tap,tas,tpcol,tscol
!    integer, dimension(:), allocatable          :: isip,itip,kips,kpsc
!    integer, dimension(:), allocatable          :: itis,kist,kstc
!    integer, dimension(:), allocatable          :: hsind,hpind
!    integer, dimension(:), allocatable          :: scheck,tcheck

!    character*15, dimension(:), allocatable	   :: pregion
!    character*15, dimension(:), allocatable	   :: sregion
!    character*15, dimension(:), allocatable	   :: tregion
!    character(6), dimension(:), allocatable    :: ptype
!    character(6), dimension(:), allocatable    :: stype
!    character(6), dimension(:), allocatable    :: ttype

!    character(150), dimension(1:MHR)  :: headrec

!!    end module mascon_setup

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
      write(lun,'(a120)') (trim(headrec(i)),i=1,nhr)
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
    end subroutine hash_header
    