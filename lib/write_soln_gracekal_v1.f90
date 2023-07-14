
    subroutine write_soln_gracekal_v1 (gracenam, luout, maxparm, prm_input, apriori, vector, sigma &
              ,date_init,datec, first_VCV_file)
!   Write apriori, solution vector and VCV to a file for both GRACEFIT and GRACEKAL 
!
! MODS:
! PT130528: revamp the writing out of the solution to use the prmnam vector that contains all the estimated parameters (rather than hardwiring the names here ....)
! PT150410: remove the include files ... pass in and dimension everything required here !
     
    implicit none

!    include '../includes/grace.param'   
!    include '../includes/write_soln.h90'

! -- PARAMETERS --  
    integer, intent(in)      :: luout                 ! unit number of output file
    integer, intent(in)      :: maxparm               ! number of parameters, for dimensioning of the matrices
    integer, intent(in)      :: date_init(5),datec(5) ! initial and current date for the gracekal run
    character*80, intent(in) :: first_VCV_file
    character*30, intent(in) :: prm_input(maxparm)    ! descriptors of the parameters estimated               
    real*8 :: apriori(maxparm), vector(maxparm), sigma(maxparm, maxparm)
    integer*4 :: i,j
    real*8 :: spanhrs
    character*8 :: dattim, gracenam, gracefit, gracekal, gracebac

! PT130110: add a reference epoch to the VCV file header
    integer date(5)
    double precision sec

!LL en fonction du programme appele change l'input
! -- GRACEFIT CASE --


! -- GRACEKAL CASE --
   if (gracenam == "gracekal") then
     call datetime( dattim )
     write(luout,'(a,a19)') 'V2 GRACEKAL solution: GRACEKAL run', dattim
     write(luout,'(a,a)') 'REFERENCE GRACEFIT file:' , first_VCV_file
!     write(luout,'(a,a)') 'Current file:' , R_matrix_filename  
     write(luout,'(a,i4,3(a1,i2.2),2x,a,i4,3(a1,i2.2))') 'First epoch:' , date_init(1), "-", date_init(2),"-", date_init(3),"-"&
                                              , date_init(4), 'Current epoch:' , datec(1), "-", datec(2),"-", datec(3)," ", datec(4)
     write(luout,*) 'Number of estimated parameters: ', maxparm

! -- GRACEKAL BACKWARD CASE --
   else if (gracenam == "gracebac") then
     call datetime( dattim )
     write(luout,'(a,a19)') 'GRACEKAL solution for backward filter: GRACEKAL run',dattim
     write(luout,'(a,a)') 'REFERENCE GRACEFIT files:' , first_VCV_file  
!     write(luout,'(a,a)') 'Current file:' , R_matrix_filename
     write(luout,'(a,i4,3(a1,i2.2),2x,a,i4,3(a1,i2.2))') 'First epoch:' ,date_init(1),"-",date_init(2),"-",date_init(3),"-"&
                                              ,date_init(4),'Current epoch:' , datec(1),"-",datec(2),"-",datec(3)," ",datec(4)
     write(luout,*) 'Number of estimated parameters: ',maxparm
   endif

! -- END OF THE HEADER --

! -- WRITE THE SOLUTIONS --
   write(luout,'(/,a)') ' SOLUTION A PRIORI AND VECTOR:'
! PT130528: take this line out - don't distinguish anymore between orbital and mascon
!   write(luout,'(a)')' Orbital parameters'
   write(luout,'(a)') ' PARAMETER                     A PRIORI             VECTOR            SIGMA'   

! PT130528: now we just loop through all the parameters and uses the pre-defined names in prm_input
   do i=1,maxparm
     write(luout,'(a,f17.7,f17.7,f17.7)')prm_input(i),apriori(i),vector(i),dsqrt(sigma(i,i))
   enddo

! write the vcv solution to the vcv file            
   write(luout,*)' VCV SOLUTION '
   do i = 1, maxparm
     write(luout,*) ((sigma(j, i)), j = 1, maxparm)   
   enddo

 return
 end Subroutine write_soln_gracekal_v1
