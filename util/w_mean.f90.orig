       program w_mean

!  program to compute a weighted mean from single obs and their
!  means.  Data is expected in 2 column format with no blank 
!  lines in file

       implicit none
       
       integer maxobs
       parameter (maxobs = 1000000)
       real*8 values(maxobs), sigma(maxobs),mean,sdev,sumrms,junk(6)
       integer i,j,ierr,column
       logical weighted
       character infile*80,tmp*20


       call getarg(1,infile)   
       if(infile(1:1).eq." ")then
         print*,'Runstring: w_mean <infile> <column> <no_weights>'
         print*,'where: infile = input file, no_weights for unweighted'
         stop
       endif

       call getarg(2,tmp)
       read(tmp,*)column

       call getarg(3,tmp)
       if(tmp(1:9).ne."no_weight")then
         weighted = .true.
       else
         weighted = .false. 
       endif

!  HM180827: changed default sigma to 1 to suppress divide by zero errors
!             in unweighted case
       ierr = 0
       do i=1,maxobs
         values(i) = 0.d0
         sigma(i) = 1.d0
       enddo       
       sumrms = 0.d0

       open(unit=10,file=infile,status='unknown')

!  read info
       i = 1
       do while (ierr.eq.0) 
         if(weighted)then
           read(10,*,iostat=ierr)values(i),sigma(i)
         else
           if(column == 1)then
             read(10,*,iostat=ierr)values(i)
           else
             read(10,*,iostat=ierr)(junk(j),j=1,column-1),values(i)
           endif
         endif
         if(ierr.eq.0)then
!           write(*,*)i,values(i),sigma(i)

!  PT011211: calculate the wrms for the input values. This may be 
!            meaningless in some applications; however, if the input
!            values are residuals already then it will be the wrms of
!            whatever process derived them.
!
!  PT011029: wrms = sqrt [sum ( V(i)/W(i)) /n ]
           sumrms = sumrms + (values(i)/sigma(i))**2

           i = i + 1
           if(i.gt.maxobs)then
             print*,' Too many observations: ',i,'. Max is ',maxobs   
             stop
           endif
         endif
       enddo

       if(weighted)then
         print*,' computing weighted mean'
       else
         print*,' computing unweighted mean of column:',column
       endif

       j=i-1

       call meancalc(values,sigma,j,mean,sdev,weighted)
        
       if(weighted)then
         write(*,1000)mean,sdev
1000     format(//,5x,'The weighted mean is ',f15.4,' +/- ',f10.4)  
         write(*,1002)dsqrt(sumrms/(j*1.d0))
1002     format(5x,'The WRMS is ',f15.4,' mm/yr')
       else
         write(*,1001)mean,sdev/dsqrt(j*1.d0)
1001     format(//,5x,'The  mean is ',f20.15,' +/- ',f20.15) 
         write(*,'(a,f10.5)')'Sigma of single observation is ',sdev
       endif

       stop
       end



! ********************************************************

       subroutine meancalc(values,sigma,j,mean,sdev,weighted)

!  subroutine to compute the mean and std dev of a set of numbers
!
!  arguments passed in :   values(j),sigma(j)
!                      :   j    number of values
!                      :  weighted, logical variable whether to compute
!                         weighted mean and std dev.
!  arguments passed out:   mean, sdev

       implicit none
       integer j,i
       real*8 values(j),sigma(j), mean,residsq,sdev,sum,sumwt,vsqr,vwsqr
       logical weighted

       residsq = 0
       sum = 0
       sumwt = 0
       vsqr = 0

! compute the mean
       do 20 i=1,j
         if(weighted)then
           sum =  sum + values(i)/sigma(i)
           sumwt = sumwt + 1.d0/sigma(i)

         else
           sum =  sum + values(i)
         endif
20     continue
         if(weighted)then
           mean = sum/sumwt  
         else
           mean = sum/(j*1.d0)
         endif

! compute standard deviation about the mean
         if(.not.weighted)then
           mean = sum/(j*1.d0)
           do 30 i=1,j
               residsq = (mean - values(i))**2 + residsq
30         continue
         else
           mean = sum/sumwt
           do 40 i=1,j   
               residsq = ((mean - values(i))**2)/sigma(i) + residsq
40         continue
         endif

         sdev = dsqrt((residsq/(j-1)*1.d0)) 
       if(weighted)then      
!  return the standard deviation of the weighted mean
         sdev = sdev/dsqrt(sumwt)  
       endif
       return
       end

