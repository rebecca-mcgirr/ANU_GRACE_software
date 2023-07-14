      program cleanup_mascons
!     *******
!
!  cleans up and rewrites a mascon file after reshaping

!
!  usage: cleanup_mascons <input mascon_file>  [-list <mascon_list_file> -x]
!          - merges small primary mascons into nearest neighbour of the same type
!	   - output file is reordered sequentially by primary mascon
! 	   - maxtip,maxsip increased to <newmax> if specified
!
!  output to mascons_cleaned
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


    use mascon_setup

    implicit none

    real (kind=8), dimension(:,:), allocatable 	:: prange
    integer, dimension(:,:), allocatable 	:: lsip,ltip,ltis
    integer, dimension(:), allocatable 		:: ltin,ltex
    integer, dimension(:), allocatable 		:: mmasc

    real (kind=8), dimension(1:20)   :: vlon,vlat   ! polygon vertices
    integer     MP,MS,MT,MM,nv
    integer 	ipm1,ipm2,ipm3

    character*6 mtype
    character*12 hashcode
    character*150 line


    integer np,ns,nt,nm,nh,nn,ip,is,it
    integer maxsip,maxtip,maxtis
    integer newmax
    integer nc_sap,nc_tap,nc_tas

! working variables

    real (kind=8) sum1,sum2,sum3,sum4,sum5,a1,a2,dif,difmin
    real (kind=8) splat,cplat,dist,rmin,rad,atmp
    integer i,j,k,is1,is2,it1,it2,im
    integer offset,ipo,iso,narg,ntin,ntex
    integer nmerge,ncut,nab
    character*120 infile,outfile,listfile
    character*14 ext
    character*6 arg2,arg3,arg4,arg,val
    logical expand,inpolyc,inlist


1   format(i7,2x,a6,2i8,2f10.4,f11.1,f18.0,2f9.1,f6.0,f6.1,i6,2x,a20)
2   format(i7,2x,a6,i8,8x,2f10.4,f11.1,f18.0,2f9.1,f6.0,2i6,2x,a20)
3   format(i7,2x,a6,2f10.4,f11.1,f15.0,2f9.1,f6.0,4i6,4x,a20)

4   format(a150)
5   format(9x,a1)
6   format(a12,2i8,i10,3i7)
7   format(9x,i4)
8   format(a20,16x,i5,a)
9   format(a,$)

! parse command line

    narg=iargc()
    if(narg.lt.1) then
       write(*,*)
       write(*,*) "usage: cleanup_mascons <input mascon_file> [-list <mascon_list_file> -x -c <ncut>] "
       write(*,*) "        merges small primaries (or listed primaries) into neighbours of the same type"
       write(*,*) "          -list : override with a specific list of mascons to merge"
       write(*,*) "          -x    : expand max #ternaries in P/S mascons (value entered later)"
       write(*,*) "          -c    : set cutoff for small primaries"
       write(*,*) "                  (primaries with <=ncut ternaries will be merged)" 
       write(*,*)
       stop
       endif

    call getarg(1,infile)
    write(*,*) "Input mascon file:  ",infile

    expand=.false.
    inlist=.false.
    ncut=100
    if(narg.gt.1) then
       i=2
       do while (i.le.narg)
          call getarg(i,arg)
          if(arg.eq."-list") then
             inlist=.true.
             call getarg(i+1,listfile)
             endif
          if(arg.eq."-x") then
             expand=.true.
             endif
          if(arg.eq."-c") then
             call getarg(i+1,val)
             read(val,*) ncut
             write(*,*) "small mascon cutoff set to ",ncut
             endif
          i=i+1
          enddo
       endif
    
    outfile=trim(infile)//"_cleaned"
    write(*,*) "Output mascon file: ",outfile
    write(*,*) "Note: primary mascons will be sorted on output"

! open mascon files

    open(10,file=infile,status="old")
    write(*,*) "Header line 1:"

! Allocate arrays

    read(10,*) hashcode,MP,MS,MT,maxsip,maxtip,maxtis
    write(*,*) hashcode,MP,MS,MT,maxsip,maxtip,maxtis
    if(expand) then
       write(*,"(a40,i5,a20)") "Expanding max ternaries in primary from",maxtip," - enter new maxtip:"
       read(*,*) newmax      
!      write(*,*) "expanding maxtip,maxsip to ",newmax
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

    allocate(prange(MP,MP))            ! primary mascon ranges
    allocate(lsip(MP,maxsip))          ! list of secondaries in primaries
    allocate(ltis(MS,maxtis))          ! list of ternaries in secondaries
    allocate(ltin(maxtis))             ! list of ternaries inside polygon
    allocate(ltex(maxtis))             ! list of ternaries outside polygon
    allocate(mmasc(MP))                ! list of primaries to absorb


    MM=MP+MS+MT
    write(*,*) "total mascons expected:",MM

! setup

    rad=pi/180.d0
    np=0
    ns=0
    nt=0
    nm=0
    nab=0
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
     if(mod(i,1000)==0) write(*,9) "."
     enddo

    close(10)
    write(*,*)
    write(*,*) "mascon range (P,S,T): ",np,ns,nt
    write(*,*) "total mascons found   :",nm

!-----------------------------------------------------

! report array mismatches

    if(np.ne.MP) write(*,*) "mascon numbering mismatch,   primaries:",np,MP
    if(np.ne.MP) write(*,*) "mascon numbering mismatch, secondaries:",ns,MS
    if(np.ne.MP) write(*,*) "mascon numbering mismatch,   ternaries:",nt,MT

!-----------------------------------------------------

! report misconnections in sap/tap/tas arrays

    if(nc_sap+nc_tap+nc_tas.gt.0) then
       write(*,*) nc_sap," S-P misconnections found"
       write(*,*) nc_tap," T-P misconnections found"
       write(*,*) nc_tas," T-S misconnections found"
       write(*,*) " (misconnections will be repaired on output)"
       endif


!-----------------------------------------------------

! create primary mascon range matrix

    do ip=1,np-1
       splat=sin(pmlat(ip)*rad)
       cplat=cos(pmlat(ip)*rad)
       do j=ip+1,np
          dist=splat*sin(pmlat(j)*rad)
          dist=dist+cplat*cos(pmlat(j)*rad)*cos((pmlon(ip)-pmlon(j))*rad)
          prange(ip,j)=acos(dist)/rad
          enddo
       enddo


!-----------------------------------------------------

! read or build a list of mascons to merge

    if(inlist) then
       write(*,*) " reading list of mascons to merge"   
       open(12,file=listfile,status="old",err=900)
       do i=1,np
          read(12,*,end=200) mmasc(i)
          enddo
200    nmerge=i-1
       close(12)
    else
       nmerge=0
       write(*,'(a,i5)') " scanning for small mascons to merge, cutoff #ternaries =",ncut  
       do ip=1,np
          if(ntip(ip).lt.ncut) then
             nmerge=nmerge+1
             mmasc(nmerge)=ip
             endif
          enddo
       endif
       
     write(*,"(i5,a)") nmerge," mascons merges will be attempted"
!     write(*,"(a,i8)") " mascon ",(mmasc(i),i=1,nmerge)
     if(nmerge.eq.0) stop " - nothing to be done!"

!-----------------------------------------------------

! mascon merge loop
!  - mascon #im will be merged into target mascon #ip
     
     do i=1,nmerge
        im=mmasc(i)
        mtype=ptype(im)
  
! identify closest suitable mascon

        rmin=360.
        ip=0
        do j=1,im-1
!          if(im.eq.471)  write(99,*) j,im,prange(j,im),mtype
           if(prange(j,im).lt.rmin) then
              if(ptype(j).eq.mtype) then
                 rmin=prange(j,im)
                 ip=j
                 endif
              endif
           enddo
        do j=im+1,np
!          if(im.eq.471) write(99,*) im,j,prange(im,j),mtype
           if(prange(im,j).lt.rmin) then
              if(ptype(j).eq.mtype) then
                 rmin=prange(im,j)
                 ip=j
                 endif
              endif
           enddo

           if(ip.eq.0) then
              write(*,'(/i5,".",4(a,i6))') i," merging mascon",im," (",ntip(im)," tern)"
              write(*,*) " !! merger aborted: no matching type? ",mtype
              nab=nab+1
              go to 300
              endif

        write(*,'(/i5,".",4(a,i6),a,f7.2,a,2f9.2)') i," merging mascon",im," (",ntip(im)," tern) into",&
                    ip," (",ntip(ip)," tern) range:",rmin," degrees at",pmlat(ip),pmlon(ip)

! abort merger if primaries are too far apart
           
           if(rmin.gt.7.d0) then
              write(*,*) " !! merger aborted: primary separation too great"
              nab=nab+1
              go to 300
              endif
           
! test that allocated arrays are large enough for merged mascon

        nn=ntip(ip)+ntip(im)
        write(*,'(a,i5,2a)') " merged mascon will contain ",nn," ternaries, type: ",ptype(ip)
        if(nn.gt.maxtis) then
           write(*,*) "Error: merged mascon is too large for arrays ( #ternaries > maxtis )"
           write(*,'(i5,a,i5,a,i5)') ntip(ip)," +",ntip(im)," >",maxtis
           stop
           endif
           
            
!  load list of ternaries to be added to target mascon
        
        do j=1,nsip(im)
           is=lsip(im,j)
           do k=1,ntis(is)
              ltin(k)=ltis(is,k)
              enddo
           ntin=ntis(is)
           enddo

!  check whether the merging mascons are contiguous
           
           rmin=360.d0
           do it=1,ntin
              it1=ltin(it)
              atmp=rad*tmlat(it1)
              splat=dsin(atmp)
              cplat=dcos(atmp)
              do j=1,nsip(ip)
                 is=lsip(ip,j)
                 do k=1,ntis(is)
                    it2=ltis(is,k)
                    dist=splat*sin(tmlat(it2)*rad)
                    dist=dist+cplat*cos(tmlat(it2)*rad)*cos((tmlon(it1)-tmlon(it2))*rad)
                    dist=acos(dist)/rad
                    rmin=min(rmin,dist)
                    enddo
                 enddo
              enddo

              if(rmin.gt.(1.5d0*10.d0/60.d0)) then
                 write(*,'(a,f8.2,a)') " mascons are NON-contiguous, closest ternaries",rmin*6," grid units"  !              else
!                 write(*,'(a,f8.2,a)') " mascons are contiguous, closest ternaries",rmin*6," ternary grid units" 
                  endif
                 
! merge mascons

        call merge(ip,im,lsip,ltis,ltin,ntin,MP,MS,MT,maxsip,maxtis)
 300   enddo
      
!-----------------------------------------------------

! update max ternaries in P/S mascons
       
    newmax=0
    do ip=1,np
       newmax=max(newmax,ntip(ip))
       enddo
       
    if(expand) then
       write(*,*) "resetting maxtip,maxtis to ",newmax
    else
       if(newmax.ne.maxtip) write(*,*) " warning: maxtip,maxtis incorrect, resetting from",maxtip,maxtis,&
                                       " to ",newmax
    endif

    maxtip=newmax
    maxtis=newmax

    
!-----------------------------------------------------

! write modified mascon file
!  - P/S mascon membership is updated on output

   open(11,file=outfile,status='unknown') 
   nmerge=nmerge-nab
   write(*,*) nmerge," mergers completed"
   
! reduce the number of primary and secondary mascons written to the header by the number removed      
   write(headrec(1),6) hashcode,MP-nmerge,MS-nmerge,MT,maxsip,maxtip,maxtis
   nh=nh+1
   if(nh.gt.MHR) stop "error: exceeded maximum header lines - increase MHR"
   write(headrec(nh),8) "# cleanup_mascons,         ",nmerge," small mascons merged into neighbours "
!88 format(a27,a,a8,i5,a20,i5)
!   write(headrec(nh),88) "# merge_mascons  , region: ",pregion(1111)," mascon ",1111," merged into mascon ",2222
   do i=1,nh
      write(11,"(a)") trim(headrec(i))
      enddo

   offset=0
   write(*,'(a,$)') " primary write loop: "
   do ip=1,MP

      if(ptype(ip).ne."X") then
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
 !		write(99,*) ip,is,k,it
               write(11,3)  it,ttype(it),tmlat(it),tmlon(it),tmr(it),tma(it),tmh(it),tmg(it),tmden(it),ipo,iso,&
								tpcol(it),tscol(it),tregion(it)
               enddo
            enddo
         else
            offset=offset+1
            endif
      if(mod(ip,1000)==0) write(*,9) "."
      enddo



    close(11)
    write(*,*)
    write(*,*) "new mascon count (P,S,T): ",MP-offset,MS-offset,MT
    if(offset.ne.nmerge) then
       write(*,*) "error: mismatch in number of mascons retired: " 
       write(*,*) "  merged:  ",nmerge
       write(*,*) "  retired: ",offset
       endif
    
!    call status_update('STATUS','UTIL','merge_mascons',outfile,"End program",0)

    stop

900 stop "mascon file error: no merge_list?"
920 stop "mascon file error: no mascon file?"
    end

!************************************************************


      logical function inpolyc(x,y,xvrt,yvrt,n)
!		       *******
!                                   
!  Tests whether the point (x,y) is in the polygon defined by xv,yv
!  using the crossing algorithm (Haines, Graphics Gems 1.4)
!
!  Principle:
!	shoots a line in the +X direction from the test point and
!	counts how many polygon sides are crossed: odd=in, even=out
!  Input: x,y - test point
!	  xvrt(n),yvrt(n) - polygon vertices
!  Output: inpolyx : true if test point is in polygon
!	   nc      : number of boundary crossings (not usually needed)
!  Comments:
!	- algorithm runs slightly faster if increment of nc (crossing 
!	   count) is commented out
!	- 
!
!-------

    real (kind=8) xv0,yv0,xv1,yv1,x,y
    real (kind=8), dimension(1:n) :: xvrt,yvrt
    integer n,nc
    logical yf0,yf1,xf0


    inpolyc=.false.
    xv0=xvrt(n)
    yv0=yvrt(n)
    yf0=yv0.ge.y
    nc=0

    do j=1,n
       xv1=xvrt(j)
       yv1=yvrt(j)
       yf1=yv1.ge.y
       if(yf0.neqv.yf1) then
	  xf0=xv0.ge.x
          if(xf0.eqv.(xv1.ge.x)) then
	     if(xf0) then
                inpolyc=.not.inpolyc
	        nc=nc+1
	        endif
	     else
               sl=(xv0-xv1)/(yv0-yv1)
	       if((xv1-(yv1-y)*sl).ge.x) then
                  inpolyc=.not.inpolyc
	          nc=nc+1
	          endif
	       endif
	     endif
         yf0=yf1
 	 xv0=xv1
	 yv0=yv1
	 enddo

     return
     end

!--------------------------------------------------------------------------------

    module mascon_setup
       
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

    end module mascon_setup


!--------------------------------------------------------------------------------

    subroutine merge(ipm1,ipm2,lsip,ltis,ltern,ntern,MP,MS,MT,maxsip,maxtis)

! merge primary ipm2 into ipm1
!   - ternary mascons listed in ltern are transferred into primary mascon ipm1

    use mascon_setup

    integer, dimension(1:MP,1:maxsip)   :: lsip
    integer, dimension(1:MS,1:maxtis)   :: ltis
    integer, dimension(1:maxsip) 	:: lsec
    integer, dimension(1:maxtis) 	:: ltern
    integer   MP,MS,MT,maxsip,maxtip
    integer   ipm1,ipm2,ntern
    
    real (kind=8) sum1,sum2,sum3,sum4,sum5,sum6
    real (kind=8) x,y,z,h,a1,a2,dif,difmin,al,as
    real (kind=8) sumx,sumy,sumz,amlat,amlon,dlat,dlon
    integer   np,ns,nt,nm,nh,nn,ip,is,it
    integer   i,j,k,is1,is2,it1

!    write(*,*) ipm1,ipm2,lsip(ipm1,1),ltis(ipm1,1),ltern(1),ntern,MP,MS,MT,maxsip,maxtis


! anonymise retired primary
!  - this prevents others later being merged into it
    
    ptype(ipm2)="X"
    
! change colours of ternaries in retired primary (ipm2)

    is1=lsip(ipm1,1)
    it1=ltis(is1,1)
    do j=1,nsip(ipm2)
       is2=lsip(ipm2,j)
       do k=1,ntern
          it=ltern(k)
          tpcol(it)=tpcol(it1)
          tscol(it)=tscol(it1)
          enddo
       enddo

! update list index for merged primary

    do i=1,ntern
       ltis(is1,ntis(is1)+i)=ltern(i)
!       write(99,*) is1,ntis(is1),i,ltern(i)
       enddo

!update merged mascon counts

!!		nsip(ipm1) is unchanged
    ntip(ipm1)=ntip(ipm1)+ntern
    ntis(is1)=ntip(ipm1)
!    do i=1,ntis(is1)
!       write(98,*) ltis(is1,i)
!       enddo

! update merged secondary's parameters
! - we assume there will only ever be one secondary in each primary
!		with the same parameters as its primary

!!		stype(is1) is unchanged

    sumx=0.d0
    sumy=0.d0
    sumz=0.d0
    sum1=0.d0
    sum2=0.d0
    sum3=0.d0
    sum4=0.d0
    sum5=0.d0
    sum6=0.d0
    rad=pi/180

    do j=1,ntis(is1)
       it=ltis(is1,j)
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
    smlat(is1)=datan2(z,h)/rad
    smlon(is1)=datan2(y,x)/rad
    if(smlon(is1).lt.0) smlon(is1)=smlon(is1)+360

! amlat,amlon are arc distances (and may be in error near 360E)

    amlat=sum1/sum4
    amlon=sum2/sum4
    dlat=smlat(is1)-amlat
    dlon=smlon(is1)-amlon
    if(dlat.gt.0.2d0) write(*,'(3(a,f8.3))') "new latitude:  chord",smlat(is1),"   arc",amlat,"   diff",dlat
    if(dlat.gt.0.2d0) write(*,'(3(a,f8.3))') "new longitude: chord",smlon(is1),"   arc",amlon,"   diff",dmlon

!    sma(is1)=sma(is1)+sma(is2)
    sma(is1)=sum4
    smh(is1)=sum5/ntis(is1)
    smg(is1)=sum6/ntis(is1)
!    write(*,*) ipm1,is1,it1,ntis(is1),ntern,smh(is1),tmh(it1)
    
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
 
!!	smden(is1) is unchanged
!!	sap(is1) is unchanged
!!	sscol(is1) is unchanged
!!	sregion(is1) is unchanged

! sum land and sea areas of expanded mascon
    
    al=0.d0
    as=0.d0
    do j=1,ntis(is1)
       it=ltis(is1,j)
       if(tmden(it).le.1001.d0) then
          al=al+tma(it)
       else
          as=as+tma(it)
          endif
       enddo
!       write(*,*) al,as
       
! update merged primary's parameters
!	- most are the same as the contained secondary
          
    pmlat(ipm1)=smlat(is1)
    pmlon(ipm1)=smlon(is1)
    pmr(ipm1)=smr(is1)
    pma(ipm1)=sma(is1)
    pmh(ipm1)=smh(is1)
    pmg(ipm1)=smg(is1)
    pcland(ipm1)=100*al/(al+as)

!!	pmden(ipm1) is unchanged
!!	tidemask(ipm1) is unchanged
!!	pregion(ipm1) is unchanged

    return
    end subroutine merge
