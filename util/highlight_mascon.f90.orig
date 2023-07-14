      program highlight_mascons
!     *******
!
!  creates a pseudo fitfile 
!
!  usage: highlight_mascons n
!  output highlight.fit
!

    implicit none
    integer, parameter      :: MP = 4582, MPB=61
    integer, parameter      :: MS = 164838, MSB=361
    integer, parameter      :: MT = 1485118, MTB=1081

	integer nh,i
	real  (kind=8) val,hival,lowval
	character*5 arg
	
1	format("       MC",i4,57x,f10.4)

    call getarg(1,arg)
    read(arg,*) nh

    open(11,file='highlight.fit',status='unknown')

	hival=0.8
	lowval=-0.5
	do i=1,MP
		if(i.eq.nh) then
			val=hival
		else
			val=lowval
			endif
    	write(11,1) i,val
		enddo
		
	close(11)
	stop
	end