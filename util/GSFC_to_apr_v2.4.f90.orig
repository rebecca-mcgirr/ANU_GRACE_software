  program GSFCmsc_to_apr

! program to output a file of a priori mascon values based on the GSFC monthly mascon estimates.

! P. Tregoning
! 19 June 2017
!
! MODS
! PT190110: added an amplitude of white noise to be added to each mascon value
! PT190111: modified to read directly the ANU mascon file (rather than some intermediate file)
! RM191223: added option to only apply noise to mascons less than specified area

  use mascon_mod  ! provides the arrays for reading the mascon file information


  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! command line arguments
  character*150  :: msc_crd           ! file of our mascon coordinates
  integer*4      :: GSFC_epoch        ! the required epoch number of the GSFC solution required
  character*3    :: out_typ           ! fit or vcv output file format
  integer*4      :: lumsc_in          ! unit number of input mascon file

! the GSFC files
  character*100  :: GSFC_msc_crd      ! file of GSFC mascon coordinates
  character*100  :: GSFC_est          ! file of estimated GSFC mascon EWH

! various variables
  real*8, allocatable :: msc_ANU(:,:)
  real*8, allocatable :: mascons_GSFC(:,:)
  real*8, allocatable :: soln_msc(:,:)
  integer*4           :: nsoln_GSFC
  real*8, allocatable :: output_EWH(:)

  integer*4  :: imsc,jmsc,nmsc_GSFC,nmsc_ANU,j,i
  character*100 :: arg,line

! variables for adding gaussian noise
  real*8 :: noise_ampl
  real*8, allocatable :: noise(:)
  real*8 :: max_area ! RM191223: option to only add noise to mascons smaller than max_area

! variables for just getting GSFC mascon information
  real*8  :: msc_latlon(3)
  logical :: which_GSFC

! variables for ensuring conservation of mass
  real(kind=8)  :: land_area,land_vol,ocean_area,ocean_vol
  real(kind=8)  :: vol_error,ocean_vol_orig,vol_error_orig,ocean_vol_corr
  real(kind=8),allocatable :: mass_corr(:)

! variables for writing out file
  character*150 :: output_file      ! hardwired file of form "GSFC_apriori_epoch_mscfile.vcv"
  integer*4     :: luout            ! unit number of input mascon file

  real*8 :: tmplat,tmplon,lat_w,lon_w,pi

  pi = 4.d0*datan(1.d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Command line !!!!!!!!!!!
  call getarg(1,msc_crd)

  if(msc_crd(1:1) == " ")then
    print*,'Runstring: GSFC_to_apr mascons_ANU 140 fit/vcv [noise_ampl max_area]'
    print*,'Specify max_area to apply noise to mascons smaller than max_area, otherwise apply noise to all mascons'
    stop
  endif

  call getarg(2,arg)
  read(arg,*)GSFC_epoch
! PT190829: if the epoch is -99 then it is a flag to only determine in which gsfc mascon a lat/lon resides
  if(GSFC_epoch < 0)then
    which_GSFC = .true.
    print*,'will determine in which GSFC mascon lies the coordinate',msc_latlon
    call getarg(3,arg)
    read(arg,*)msc_latlon(1)
    call getarg(4,arg)
    read(arg,*)msc_latlon(2)
    msc_latlon(3) = 0.d0
  else
    which_GSFC = .false.
    call getarg(3,out_typ)
! PT190110: add an amplitude of gaussian noise to add to each mascon
    call getarg(4,arg)
    if(arg(1:1) == "")then
      noise_ampl = 0.d0
    else
      read(arg,*)noise_ampl
    endif
! RM191223: only apply noise to mascons smaller than max_area
    call getarg(5,arg)
    if(arg(1:1) == "")then
      max_area = 0.d0
    else
      read(arg,*)max_area
    endif
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(.not. which_GSFC)then
! open the input mascon coordinate file
! PT190111: do this using the standard subroutine calls
    lumsc_in = 10
    call read_msc_hdr(lumsc_in,msc_crd,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                                ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
    call allocate_mascon_arrays
    call read_mascon_file(lumsc_in,msc_crd)
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the GSFC mascon coords
  open(11,file="mascon_group.dat",status='old')
  do i=1,3
    read(11,'(a)')line
  enddo

! get the number of GSFC mascons
  read(11,'(a)')line
  read(line,'(19x,i6)')nmsc_GSFC
print*,'number of mascons:',nmsc_GSFC
! allocate arrays
  allocate(mascons_GSFC(nmsc_GSFC,4))    ! the four values per mascon are minlat,maxlat,minlon,maxlon

! skip 10 lines
  do i=1,10
    read(11,'(a)')line
  enddo

! now read in all the mascons and store as the min/max lon/lat for each one
  do imsc=1,nmsc_GSFC
    read(11,*)tmplat,tmplon,lat_w,lon_w
! make the GSFC longitudes between 0 and 360
    if (tmplon < 0)tmplon = tmplon + 360.d0

! assign the min/max extents of the mascons
    mascons_GSFC(imsc,1) = tmplat - lat_w/2.d0
    mascons_GSFC(imsc,2) = tmplat + lat_w/2.d0
    mascons_GSFC(imsc,3) = tmplon - lon_w/2.d0
    mascons_GSFC(imsc,4) = tmplon + lon_w/2.d0

 !print*,tmplat,tmplon,mascons_GSFC(imsc,1:4)

  enddo
  close(11)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get the GSFC mascon EWH values
  open(12,file="solution_group.dat",status='old')
! skip 3 lines
  do i=1,3
    read(12,'(a)')line
  enddo

! get the number of GSFC mascon solutions (columns)
  read(12,'(a)')line
  read(line,'(39x,i6)')nsoln_GSFC
print*,'number of solutions:',nsoln_GSFC

! skip 3 lines
  read(12,'(a)')line
  read(12,'(a)')line
  read(12,'(a)')line

! dimension arrays
  allocate(soln_msc(nmsc_GSFC,nsoln_GSFC))

! now read the whole thing into an array
  do imsc=1,nmsc_GSFC
    read(12,*)(soln_msc(imsc,j),j=1,nsoln_GSFC)
  enddo
  close(12)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(which_GSFC)then
! we only want to know in which GSFC mascon an input lat/lon resides. We can work this out, then exit
    imsc = 1
    call which_GSFC_mascon(jmsc,imsc,nmsc_GSFC,msc_latlon,mascons_GSFC)
    print*,'Location ',msc_latlon(1:2),' is in GSFC mascon: ',jmsc

! output all the solutions for this mascon
    do i=1,nsoln_GSFC
      write(*,'(f12.2,a,i10)')soln_msc(jmsc,i)*10,"         GSFC soln for mascon ",jmsc
    enddo
    return
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read in the ANU mascon coordinates
!  read(10,*)nmsc_ANU
! ! allocate the array
!! PT180202: increase from 2 to three columns so that we read in and store the density
!  allocate(msc_ANU(nmsc_ANU,3))
!
!! read them in
!  do imsc=1,nmsc_ANU
!    read(10,*)msc_ANU(imsc,:)
!  enddo
!  close(10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT190111: transfer the information from the mascon arrays into the msc_ANU array
  allocate(msc_ANU(total_prim,3))
  do imsc=1,total_prim
    msc_ANU(imsc,1:2) = mcon_prim(imsc,1:2)*180.d0/pi   ! decimal lat/lon
    msc_ANU(imsc,3)   = mcon_prim(imsc,6)               ! density
!print*,'msc_ANU(imsc,:)',msc_ANU(imsc,:)
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! define the gaussian noise model to be added to each mascon
  allocate(noise(total_prim))

! set up the random noise generator
! RM191223: only calculate noise if area of mascon is below some threshold if specified in command line
  call srand(total_prim)
  if(max_area == 0.d0)then
    do i=1,total_prim
      noise(i) = (rand(0)-0.5)*noise_ampl
    enddo 
  else 
    do i=1,total_prim
      if(mcon_prim(i,4) < max_area)then ! test whether mascon is smaller than some area threshold
        noise(i) = (rand(0)-0.5)*noise_ampl
      else
        noise(i) = 0.0 ! else noise is 0
      endif
    enddo 
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RM200226: write a VCV file header
  if(out_typ == "vcv")then
    luout = 13
    write(output_file,'(a,i0.3,a,a,a)')"GSFC_apriori_",GSFC_epoch,"_",trim(msc_crd),".vcv"
    open(luout,file=output_file,status='unknown')
    write(luout,'(a,i3,a)')"V2 APRIORI  VCV  created for ",GSFC_epoch," from GSFC mascon solution "
    write(luout,'(i6,a)')total_prim," mascons in file"
    write(luout,'(a,a)')"Output primary mascon geometry file: ",trim(msc_crd)
    write(luout,'(a)')"An incomplete VCV file, containing only mascon values to be read"

    write(luout,'(a)')"SOLUTION A PRIORI AND VECTOR:"
    write(luout,'(a)')"PARAMETER                     A PRIORI             VECTOR            SIGMA"
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! find the GSFC mascon in which each ANU mascon resides
  allocate(output_EWH(total_prim))
  do imsc = 1,total_prim
    call which_GSFC_mascon(jmsc,imsc,nmsc_GSFC,msc_ANU(imsc,:),mascons_GSFC)

! PT180202: there are some GSFC mascons with continental-magnitude EWH values but
!                   they fall in the ANU ocean. Set these to zero EWH value.
    if(dabs(soln_msc(jmsc,GSFC_epoch))/1.d2  > 0.05d0 .and. msc_ANU(imsc,3) > 1001.d0)then 
    !if(msc_ANU(imsc,3) > 1001.d0)then ! make all ocean mascons zero 
! it is a GSFC land mascon that lies in the ANU ocean
!       print*,'replace GSFC value of ',soln_msc(jmsc,GSFC_epoch)/1.d2,' mm with 0.000 for the ANU ocean',imsc,msc_ANU(imsc,:)
       output_EWH(imsc) = 0.d0
! RM200220: there are some continental island mascons that contain too few ternarys, set these to zero EWH
    else if(int(mcon_prim(imsc,8)) < 60)then
       output_EWH(imsc) = 0.d0
    else
      output_EWH(imsc) = soln_msc(jmsc,GSFC_epoch)/1.d2 + noise(imsc)
      !if(output_EWH(imsc) > 0.d0)then
      !   output_EWH(imsc) = output_EWH(imsc) - 0.1
      !else
      !   output_EWH(imsc) = output_EWH(imsc) + 0.1
      !endif
!print*,imsc,soln_msc(jmsc,GSFC_epoch)/1.d2,noise(imsc),msc_ANU(imsc,:)
    endif
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! need to make sure that the a priori model is mass conserving before we write it out to the vcv file ...
! calculate the sum of the volume changes across all mascons. That is, sum(a priori EWH x mascon area) should be zero
  land_area = 0.d0
  land_vol = 0.d0
  ocean_area = 0.d0
  ocean_vol = 0.d0
  do imsc = 1,total_prim
    if(mcon_prim(imsc,4) < 1010.d0)then   ! it is a land primary mascon
      land_area = land_area + mcon_prim(imsc,3)
      land_vol  = land_vol  + mcon_prim(imsc,3)*output_EWH(imsc)
    else
      ocean_area = ocean_area + mcon_prim(imsc,3)
      ocean_vol  = ocean_vol  + mcon_prim(imsc,3)*output_EWH(imsc)
    endif
  enddo

! was mass conserved?
  vol_error = land_vol + ocean_vol
  vol_error_orig = vol_error
  print*,"Original mass conservation error of ",vol_error*1.d3/(land_area+ocean_area)," mm EWH across Earth, or " &
         ,1.d3*vol_error/ocean_area," mm GSL"

! take the residual volume and assign a correction across all ocean mascons, proportional to the area of each ocean mascon
  allocate(mass_corr(total_prim))
  mass_corr = 0.d0

  do imsc = 1,total_prim
    if(mcon_prim(imsc,4) > 1010.d0)then
      mass_corr(imsc) = -1.d0*vol_error/ocean_area ! * mcon_prim(imsc,3)/ocean_area
    endif
    if(out_typ == "vcv")then
      ! write it out in GRACEFIT vcv format
      ! RM200226: write out to file
      write(luout,'(f6.0,1x,a2,i5.5,a7,9x,f17.7,f17.7,f17.7)')dble(imsc),'MC',imsc,' (Gt)   ',0.d0 &
         ,output_EWH(imsc) + mass_corr(imsc),0.d0
      !write(luout,'(f6.0,1x,a2,i5.5,a7,9x,f17.7,f17.7,f17.7)')dble(imsc),'MC',imsc,' (Gt)   ',0.d0 &
      !   ,output_EWH(imsc),0.d0

    else if (out_typ == "fit")then
      write(*,'(f6.0,1x,a2,i4.4,a7,20x,f7.5,f12.5,f19.5,f12.5,f9.1 )')dble(imsc+24),'MC',imsc,' (Gt)   ',0.d0 &
         ,output_EWH(imsc) + mass_corr(imsc),output_EWH(imsc) + mass_corr(imsc),0.01,0.1
    endif
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! check that it is now mass conserving
  ocean_vol_orig = ocean_vol
  land_area = 0.d0
  land_vol = 0.d0
  ocean_area = 0.d0
  ocean_vol = 0.d0
  ocean_vol_corr = 0.d0
  do imsc=1,total_prim
    if(mcon_prim(imsc,4) < 1010.d0)then   ! it is a land primary mascon
      land_area = land_area + mcon_prim(imsc,3)
      land_vol  = land_vol  + mcon_prim(imsc,3)*output_EWH(imsc)
    else
      ocean_area = ocean_area + mcon_prim(imsc,3)
      ocean_vol  = ocean_vol  + mcon_prim(imsc,3)*(output_EWH(imsc)+mass_corr(imsc))  ! volume in m^3
      ocean_vol_corr = ocean_vol_corr + mcon_prim(imsc,3)*mass_corr(imsc)
    endif
  enddo
! was mass conserved?
  vol_error = land_vol + ocean_vol
  print*,"Adjusted mass conservation error of ",vol_error*1.e3/(land_area+ocean_area)," mm EWH across Earth, or " &
          ,1.e3*(land_vol + ocean_vol)/ocean_area," mm GSL"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(out_typ == "vcv")then
    print*,"Written GSFC solution to apriori file: ",output_file
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end



  subroutine which_GSFC_mascon(jmsc,anu_msc,nmsc_GSFC,latlon,mascons_GSFC)

! subroutine to work out in which GSFC mascon a particular coordinate resides
!
! P. Tregoning
! 19 June 2017

  implicit none

  integer*4,  intent(in ) :: nmsc_GSFC
  integer*4,  intent(out) :: jmsc,anu_msc
  real*8   ,  intent(in ) :: mascons_GSFC(nmsc_GSFC,4)
  real*8   ,  intent(in)  :: latlon(3)

! local variables
  integer*4 :: imsc 
  logical found
  real*8  diffs(4)              
  logical debug

  debug = .false.
  found = .false.
  jmsc = 0

  imsc = 1
! loop through GSFC mascons until we find one in which the points reside
  do while (imsc <= nmsc_GSFC .and. .not. found )
!print*,'imsc=',imsc,mascons_GSFC(imsc,:)
!print*,latlon,mascons_GSFC(imsc,:)
    diffs(1) = dabs(latlon(1) - mascons_GSFC(imsc,1))
    diffs(2) = dabs(latlon(1) - mascons_GSFC(imsc,2))
    diffs(3) = dabs(latlon(2) - mascons_GSFC(imsc,3))
    diffs(4) = dabs(latlon(2) - mascons_GSFC(imsc,4))

!! DEBUG
!if(debug)then
!  if(anu_msc > 15591 .and. latlon(2) > 119.9 .and. mascons_GSFC(imsc,1) > -76. .and. mascons_GSFC(imsc,2) < -72. &
!                 .and. mascons_GSFC(imsc,3) > 118. .and. mascons_GSFC(imsc,4) < 130. ) then
!    print*,anu_msc,imsc,latlon,mascons_GSFC(imsc,:)
!  endif
!endif

    !if( dabs(latlon(1) - mascons_GSFC(imsc,1)) < 0.1d0 .and. dabs(latlon(1) - mascons_GSFC(imsc,2)) < 0.1d0 &
       !.and. dabs(latlon(2) - mascons_GSFC(imsc,3)) < 0.1d0  .and. dabs(latlon(2) - mascons_GSFC(imsc,4)) < 0.1d0 ) then

! RM190219: finds which GSFC mascon the ANU mascon coordinate is within or nearest (threshold = 0.01 deg)
! PT190314: fixed bug here. The "greater than" and the "or" need to be grouped before doing the "and"
    if( ((latlon(1) >= mascons_GSFC(imsc,1)) .or. dabs(latlon(1) - mascons_GSFC(imsc,1)) < 0.085d0) .and. &
      ((latlon(1) <= mascons_GSFC(imsc,2)) .or. dabs(latlon(1) - mascons_GSFC(imsc,2)) < 0.085d0) .and. &
      ((latlon(2) >= mascons_GSFC(imsc,3)) .or. dabs(latlon(2) - mascons_GSFC(imsc,3)) < 0.085d0) .and. &
      ((latlon(2) <= mascons_GSFC(imsc,4)) .or. dabs(latlon(2) - mascons_GSFC(imsc,4)) < 0.085d0) )then

! we found the mascon
      found = .true.
      jmsc = imsc

! PT190111: is it in the spherical cap over the South Pole?
! PT190228: make this really close to the south pole
    else if (latlon(1) < -89.9d0)then
      found = .true.
      jmsc = imsc

! check that we didn't just miss it due to roundoff error
    elseif(diffs(1) < 0.01 .or. diffs(2) < 0.01)then
      if((latlon(2) > mascons_GSFC(imsc,3))  .and. (latlon(2) <= mascons_GSFC(imsc,4)))then
        found = .true.
        jmsc = imsc
      endif
    elseif(diffs(3) < 0.01 .or. diffs(4) < 0.01)then
      if((latlon(1) > mascons_GSFC(imsc,1))  .and. (latlon(1) <= mascons_GSFC(imsc,2)))then
        found = .true.
        jmsc = imsc
      endif
 

    endif
    imsc = imsc + 1

  enddo


  if (jmsc == 0)then
    print*,'did not find a match for lat/lon:',latlon
    stop 'program halted'
  endif

  return
  end


  






