!PROGRAM write_soln to write the output from gracefit and gracekal in the same format

! L. Lescarmontier from write_soln gracefit
! 15 January 2013
!
! MODS:
! PT130528: revamp the writing out of the solution to use the prmnam vector that contains all the estimated parameters (rather than hardwiring the names here ....)

    subroutine write_soln_v2 (gracenam, luvcv, N, apriori_vals, vector, sigma)
!     subroutine write_soln (gracenam,luvcv,N)      
         
!   Write apriori, solution vector and VCV to a file for both GRACEFIT and GRACEKAL 
     
    use soln_mod

    implicit none

!    include '../includes/grace.param'   
!    include '../includes/write_soln.h90'

! -- PARAMETERS --                 
   
    real*8 :: apriori_vals(maxparm), vector(maxparm), sigma(maxparm, maxparm)
    integer*4 :: N             ! number of estimated parameters
    integer*4 :: luvcv         ! unit number of vcv file

    integer*4 :: i,j

    real*8 :: spanhrs

    character*8 :: dattim, gracenam, gracefit, gracekal, gracebac

! PT130110: add a reference epoch to the VCV file header
    integer date(5)
    double precision sec

!LL en fonction du programme appele change l'input
! -- GRACEFIT CASE --

   if ( gracenam == "gracefit" ) then

!   Prepare and write out headers
     call datetime( dattim )
     write(luvcv, '(a,a19)') 'V2 GRACEFIT solution: GRACEFIT run ',dattim
     write(luvcv, '(a,a,1x,a)') 'Reference GTORB files: ', (trim(gt_fnam(i)), i = 1, 2)
     spanhrs = nepoch*5.d0/3600.d0
     write(luvcv, '(a,f6.1,a)') 'Reference GTORB files span: ',spanhrs,' hrs'
! PT130110: add a line, being the reference epoch of the GRACEFIT solution (ie the middle)
     call gsec_to_ymdhms( dble(gt_stime(1)+nepoch*2.d0/2.d0), date, sec )
     write(luvcv,'(a,i4,4(a1,i2.2),a1,f05.2)') 'Reference GRACEFIT solution epoch: ' &
               ,date(1),"-",date(2),"-",date(3)," ",date(4),"-",date(5),"-",sec 
     write(luvcv,*) 'Number of estimated parameters: ',N
! LL130213 to add: KBR uncertainty, GPS position uncertainty, Number of KB observations, Number of GPS
! observations, total number of epochs
     write(luvcv,'(a,f8.6)') 'KBRR uncertainty (um): ' ,postfit_rms_hk
     write(luvcv,'(a,f8.6)') 'GPS position uncertainty (m): ' ,postfit_rms_hg
     write(luvcv,'(a,i4)') 'Total Number of KBR observations: ' ,nkepochs
     write(luvcv,'(a,i4)') 'Number of missing KBR observations: ' ,fkepochs
     write(luvcv,'(a,i4)') 'Total Number of GPS observations: ' ,ngepochs
     write(luvcv,'(a,i4)') 'Number of missing GPS observations: ',fgepochs
     write(luvcv,'(a,i4)') 'Total Number of epochs: ' ,iepoch

! -- GRACEKAL CASE --
   else if (gracenam == "gracekal") then
     call datetime( dattim )
     write(luvcv,'(a,a19)') 'V2 GRACEKAL solution: GRACEKAL run', dattim
     write(luvcv,'(a,a)') 'REFERENCE GRACEFIT file:' , R_matrix_filename_initial
     write(luvcv,'(a,a)') 'Current file:' , R_matrix_filename  
     write(luvcv,'(a,i4,3(a1,i2.2),2x,a,i4,3(a1,i2.2))') 'First epoch:' , date_init(1), "-", date_init(2),"-", date_init(3),"-"&
                                              , date_init(4), 'Current epoch:' , datec(1), "-", datec(2),"-", datec(3)," ", datec(4)
     write(luvcv,*) 'Number of estimated parameters: ', N

! -- GRACEKAL BACKWARD CASE --
   else if (gracenam == "gracebac") then
     call datetime( dattim )
     write(luvcv,'(a,a19)') 'GRACEKAL solution for backward filter: GRACEKAL run',dattim
     write(luvcv,'(a,a)') 'REFERENCE GRACEFIT files:' , R_matrix_filename_initial  
     write(luvcv,'(a,a)') 'Current file:' , R_matrix_filename
     write(luvcv,'(a,i4,3(a1,i2.2),2x,a,i4,3(a1,i2.2))') 'First epoch:' ,date_init(1),"-",date_init(2),"-",date_init(3),"-"&
                                              ,date_init(4),'Current epoch:' , datec(1),"-",datec(2),"-",datec(3)," ",datec(4)
     write(luvcv,*) 'Number of estimated parameters: ',N
   endif

! -- END OF THE HEADER --

! -- WRITE THE SOLUTIONS --
   write(luvcv,'(/,a)') ' SOLUTION A PRIORI AND VECTOR:'
! PT130528: take this line out - don't distinguish anymore between orbital and mascon
!   write(luvcv,'(a)')' Orbital parameters'
   write(luvcv,'(a)') ' PARAMETER                     A PRIORI             VECTOR            SIGMA'   

! PT130528: write_soln_v2 just loops through all the parameters and uses the pre-defined names in prmnam
   do i=1,N
     write(luvcv,'(a,f17.7,f17.7,f17.7)')prmnam(i),apriori(i),vector(i),dsqrt(sigma(i,i))
   enddo

! write the vcv solution to the vcv file            
   write(luvcv,*)' VCV SOLUTION '
   do i = 1, N
     write(luvcv,*) ((sigma(j, i)), j = 1, N)   
   enddo

 return
 end Subroutine write_soln_v2
