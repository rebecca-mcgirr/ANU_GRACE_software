program repack_vcv
!     *******
!
!  Reads a multi-solution addnorm.vcv file, 
!		interactively selects a particular solution,
!		and writes it as a standard .vcv file 
!		(for input to gracefit in postfit mode)

!
!  usage: repack_vcv [input_file [output_file] ]
!
!  default output to input_file_rev.vcv
!

    implicit none
    integer, parameter      :: MP = 8000, MSP=12, MSOL=100

    real (kind=8), dimension(1:MSOL,1:2,1:MSP)	:: apr,vec,sig
	character*24,  dimension(1:MSOL,1:2,1:MSP)	:: str

    character*120 infile,outfile
	character*180 line
	character*200000 cvline
	integer i,j,k,ks,nsol,kso,jp,nhl,nmasc,ncv

! PT170713: header variable declarations
    integer*4     :: ymdhms(msol,6)         ! epoch of each set of ICs in the addnorm file
    real(kind=8)  :: duration(msol)         ! duration of the orbit (decimal hours)
    integer*4     :: nsoln                 ! number of solutions in addnorm file
    character*100 :: neq_file(msol)         ! names of normal equation files used in combined addnorm file

! PT170714: variables for how much of the VCV to write out
    real(kind=8),allocatable  :: msc_values(:)
    integer*4     :: nICs

1   format(2i3,a24,f17.7,f19.9,f15.7)
2   format(i3,".",a24,2x,f17.7,f19.9,f15.7)
3    format(a)

!5    format(9x,a1)
!9    format(i8,$)


! parse command line

    call getarg(1,infile)
    if(infile.eq.' ') then
		write(*,*)
        write(*,*) "usage: repack_vcv [input_file [output_file] ]"
		write(*,*)
		stop
        endif
        write(*,*)
        write(*,*) "  Input file ",infile

    call getarg(2,outfile)
    if(outfile.eq.' ') then
		outfile=infile(1:len_trim(infile)-4)//"_rev.vcv"
        write(*,*) "  Output file not specified, using: ",outfile
     else
         write(*,*) "  Output file ",outfile
         endif


! setup
        
	nhl=0
	nmasc=0
	ncv=0
	nsol=1

!  allocate(pmn(MP),nsip(MP),ntip(MP),tidemask(MP))

! open vcv files

	open(10,file=infile,status="old")
    open(11,file=outfile,status='unknown')
    write(*,*) "processing..."

! read and write header section unchanged

100 read(10,3,end=900) line
    if(line(2:10).ne."PARAMETER") then
		nhl=nhl+1
        write(*,"(a)") trim(line)
        if(line(1:4) == "File")then
          nsoln = nsoln + 1
          read(line(19:47),'(6i5)')(ymdhms(nsoln,i),i=1,6)    ! epoch of the solution
print*,(ymdhms(nsoln,i),i=1,6)
          read(line(61:110),*)duration(nsoln),neq_file(nsoln)  ! duration of solution and input normal equation file for this solution
        endif
        go to  100
    endif

!    write(11,"(a)") trim(line)

! read solution vectors until mascon encountered

110 read(10,3,end=910) line
    read(line,*) k
	if(k.ne.nsol) stop "k.ne.nsol"
    backspace(10)
    do i=1,2
    do j=1,MSP
        read(10,1) kso,jp,str(k,i,j),apr(k,i,j),vec(k,i,j),sig(k,i,j)
		if(kso.ne.k) then
			write(*,*) "mismatched satellite numbers: k,kso ",k,kso
			stop "error exit"
			endif
		if(jp.ne.j) then
			write(*,*) "mismatched parameter numbers: j,jp ",j,jp
			stop "error exit"
			endif
                enddo
		enddo

	read(10,3,end=910) line
	if(line(8:9).ne."MC") then
        nsol=nsol+1
        if(nsol.eq.MSOL) stop " too many solutions, increase MSOL and recompile"
		backspace(10)
		go to  110
		endif

	backspace(10)

! select desired solution

	write(*,*)
	write(*,*) "header lines      :",nhl
	write(*,*) "solutions found   :",nsol
	write(*,'(a21,$)') "  - select solution: "
	read(*,*) ks
	if(ks.lt.1.or.ks.gt.nsol) stop "out of range"

! write selected solution
! write the GRACEFIT format header for the VCV file
  write(11,'(a,a)')"V2 REPACK_VCV solution extracted from an ADDNORM VCV file",trim(infile)
  write(11,'(a,a)')"Reference GTORB files: ",neq_file(ks)
  write(11,'(a,f4.1,a)')"Reference GTORB files span:  ",duration(ks)," hrs"
  write(11,'(a,i4,"-",i2.2,"-",i2.2,i3,"-",i2.2,"-",i2.2,a)')"Reference GRACEFIT solution epoch: ",ymdhms(ks,:),".00"
  write(11,'(a)')"Number of estimated parameters:    7550   7526      0"
  write(11,'(a)')"Number of KB observations and used: ****"
  write(11,'(a)')"Number of missing KB observations:    0"
  write(11,'(a)')"Number of KBRR misfits:    0"
  write(11,'(a)')"Number of KBRR observations used:    0"
  write(11,'(a)')"Number of GPS observations used for GRACE A: ****"
  write(11,'(a)')"Number of relevant missing GPS observations for GRACE A:    1"
  write(11,'(a)')"Number of GPS observations used for GRACE B: ****"
  write(11,'(a)')"Number of relevant missing GPS observations for GRACE B:    1"
  write(11,'(a)')"Total Number of epochs: 0"
  write(11,'(a)')" "
  write(11,'(a)')" SOLUTION A PRIORI AND VECTOR:"
  write(11,'(a)')" PARAMETER                     A PRIORI             VECTOR            SIGMA"

        print*,'Writing out IC solution'
	do i=1,2
	  do j=1,MSP
! PT170714: there is one space too many in the addnorm string, so correct it here
		write(11,'(i3,".",a8,a15,3x,f17.7,f19.9,f15.7)') j+(i-1)*12,str(1,i,j)(1:8),str(1,i,j)(10:24) &
                  ,apr(ks,i,j),vec(ks,i,j),sig(ks,i,j)
	  enddo
	enddo

! read and write mascons unchanged

        print*,'Writing out mascon estimates'
120	read(10,3,end=920) line
	if(line(8:9).eq."MC") then
		nmasc=nmasc+1
		write(11,3) trim(line)
		go to 120
		endif
	backspace(10)

        nICs = 24                           ! hardwire this for now but it should be worked out from reading the ICs in the addnorm vcv file
        allocate(msc_values(nmasc+nICs))
! read and write vcv matrix unchanged

        print*,'Writing VCV matrix',nmasc,nICs
! PT170714: first, the "VCV Solution" line
        read(10,3,end=150) cvline
        write(11,'(a20)')cvline(1:20)

! now loop over the rest of the VCV
130	read(10,3,end=150) cvline
!        print*,cvline
! PT170714: read in just the mascon part
        read(cvline,*)(msc_values(i),i=1,nmasc+nICs)   ! this reads only the first IC block, not the correct IC block. Something to fix !!!!
	ncv=ncv+1
	write(11,*) msc_values
	go to 130


150	close(10)
	close(11)


	write(*,*) "mascons found     :",nmasc
	write(*,*) "variances found   :",ncv



!    call status_update('STATUS','UTIL','repack_vcv',"some.vcv","End program",0)

    stop

! error exits

900 stop "addnorm file error: no parameters?"
910 stop "addnorm file error: no mascons?"
920 stop "addnorm file error: no covariance matrix?"

    end

!************************************************************

