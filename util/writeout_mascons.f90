  program writeout_mascons

! program to write out a subset of continental mascons from a given input fit file solution. The mass
! change on those mascons in the selected region(s) will be written out, along with the corresponding
! change over the oceans. 
!
! The ocean signal will be the negative of the land signal and will be allocated over the ocean based
! upon a fingerprint pattern that accounts for gravitational attraction and rotational effects
!
! P. Tregoning
! 7 August 2020

  use mascon_mod

  implicit none

  character*150        :: command_file, mascon_file, soln_file        ! command line argument input files
  integer*4,parameter  :: lucmd=10,lumsc=11,lusoln=12,lufprint=13,lugia=14     ! unit numbers for input files

  integer*4                :: n_regions              ! number of regions in command file
  character*1, allocatable :: include_region(:)      ! Y/N whether to include this region or not
  character*10,allocatable :: region(:)              ! names of regions
  character*20,allocatable :: fingerprints(:)        ! names of fingerprint files for regions
  character*1              :: correct_for_gia        ! first line of input command file controls whether to apply a gia correction or not

! variables for output file
  character*150        :: output_file
  integer*4, parameter :: luout=15

! variables for fingerprint spherical harmonics
  real(kind=8),allocatable :: C_fprint(:,:),S_fprint(:,:) ! spherical harmonic coefficients of a fingerprint
  integer*4                :: maxdeg                      ! maximum degree of the fingerprint spherical harmonic file
  real(kind=8)             :: fprint_value                ! computed fingerprint value at a single location
  integer*4                :: tmpdeg,tmpord,ideg,iord
  real(kind=8)             :: tmpC,tmpS

! variables for ocean calcs + output
  real(kind=8)  :: total_area,total_volume,total_ocean_area,total_EWH,region_EWH
  integer*4     :: imsc,used_regions
  real(kind=8),allocatable :: Plm(:,:)
  character*30  :: prmnam

! GIA correction variables
  real(kind=8)  :: gia_dt
  integer*4     :: yr,month,i

  integer*4     :: iregion,nmsc,ioerr
  character*160 :: message,line
  real(kind=8)  :: msc_apr(20000),msc_est(20000)     ! mascon solution apriori and estimated EWH values
  real(kind=8)  :: msc_output(20000)
  real(kind=8)  :: msc_sigma(20000)
  real(kind=8)  :: msc_gia(20000)

  real(kind=8)  :: pi
  character*10  :: junk

  pi = 4.d0*datan(1.d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! decode command line
  call getarg(1,command_file)
  if(command_file(1:1) == "")then
    call status_update('FATAL','UTIL','writeout_mascons',' '  &
                     ,"Runstring: writeout_mascons msc.cmd mascons_stage5_V004 addnorm.fit",0)
  endif
  call getarg(2,mascon_file)
  call getarg(3,soln_file)
  call getarg(4,output_file)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! open and read the command file for this program
  open(lucmd,file=command_file,status='old')
  ! correct for GIA ?
  read(lucmd,'(a)')correct_for_gia
  ! read and store whether we want each region
  read(lucmd,*)n_regions
  allocate(include_region(n_regions))
  allocate(region(n_regions))
  allocate(fingerprints(n_regions))

  do iregion=1,n_regions
    read(lucmd,'(a1,1x,a10,1x,a20)')include_region(iregion),region(iregion),fingerprints(iregion)
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! open and read the mascon file
  open(lumsc,file=mascon_file,status='old')

  call read_msc_hdr(lumsc,mascon_file,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                                ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
  call allocate_mascon_arrays
  call read_mascon_file(lumsc,mascon_file)  
  close(lumsc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! open and read the file containing GIA rate corrections for each mascon
  if(correct_for_gia == "Y")then
    open(lugia,file="mascon_GIA_batch6.txt",status='old')
    ! read one header line
    read(lugia,'(a)')line
    call status_update('STATUS','UTIL','writeout_mascons','~/gg/grace/tables/mascon_GIA_batch6.txt',line,0)
    ! read the GIA uplift rates (in m/yr EWH)
    do imsc=1,max_prim
      read(lugia,*)junk,msc_gia(imsc)
    enddo
    close(lugia)
  else
    msc_gia = 0.d0
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the mascon EWH values from the solution file. Do this quickly and locally rather 
! than making use of existing subroutines ....
  open(lusoln,file=soln_file,status='old')
  ioerr = 0
  line = " "
  ! the fourth line of an addnorm fit file contains the date of the first .norm file. Use that as the epoch for the GIA correction
  do i=1,4
    read(lusoln,'(a)')line
  enddo
  read(line(20:28),*)yr,month
  gia_dt = dble(yr)+dble(month)/12.d0 - 2002.d0

  do while (line(7:9) /= " MC")
    read(lusoln,'(a)',iostat=ioerr)line
    if(line(7:9) /= " MC")write(luout,'(a)')line     ! transfer it to the output file
  enddo
  backspace(lusoln)

  ! we are now at the line of the first mascon
  nmsc = 1
  do while (ioerr == 0)
    read(lusoln,100,iostat=ioerr)line(1:1),msc_apr(nmsc),msc_est(nmsc),msc_sigma(nmsc)
100 format(a1,33x,f13.5,18x,f13.5,f13.5)
    if(ioerr == 0 .and. line(1:1) /= "M")then
      nmsc = nmsc + 1
    endif
  enddo
  nmsc = nmsc - 1
  write(message,'(a,i7)')"Number of mascons: ",nmsc
  call status_update('STATUS','UTIL','writeout_mascons',soln_file,message,0) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! loop through each region
  msc_output = 0.d0
  total_ocean_area = 0.d0
  used_regions = 0
  do iregion = 1,n_regions
    if(include_region(iregion) == "Y")then
      used_regions = used_regions + 1
      call status_update('STATUS','UTIL','writeout_mascons',fingerprints(iregion),'Opening fingerprint file',0)
      open(lufprint,file=fingerprints(iregion),status='old')
      read(lufprint,*)maxdeg
      allocate(C_fprint(0:maxdeg,0:maxdeg))
      allocate(S_fprint(0:maxdeg,0:maxdeg))
     
      ioerr = 0
      do while (ioerr == 0 .and. tmpdeg < maxdeg)
        read(lufprint,*,iostat=ioerr,end=1001)tmpdeg,tmpord,tmpC,tmpS
        if(ioerr == 0)then
          C_fprint(tmpdeg,tmpord) = tmpC
          S_fprint(tmpdeg,tmpord) = tmpS
        endif
      enddo
1001  close(lufprint)

      ! loop through the mascons
      total_EWH = 0.d0
      total_area = 0.d0
      do imsc = 1,max_prim

        ! we need to sum the area of the oceans
        if(used_regions == 1 .and. mcon_prim(imsc,6) > 1000.d0)then
          total_ocean_area = total_ocean_area + mcon_prim(imsc,4)
        endif

! West Antarctica
        if(region(iregion)(1:10) == "West_Antar")then
          if(mcon_prim(imsc,1) < -pi/3.d0 .and. mcon_prim(imsc,2)*180.d0/pi >= 200.d0 &
                                          .and. mcon_prim(imsc,2)*180.d0/pi <  320.d0 &
                                          .and. mcon_prim(imsc,6) < 1001.d0 )then
            ! we include this mascon in the output
            msc_output(imsc) = msc_output(imsc) + msc_est(imsc)
            total_volume = total_volume + (msc_est(imsc)-gia_dt*msc_gia(imsc))*mcon_prim(imsc,4)   ! EWH x primary mascon area
            total_area = total_area + mcon_prim(imsc,4)
          endif

! East Antarctica
        else if(region(iregion)(1:10) == "East_Antar")then
          if(mcon_prim(imsc,1) < -pi/3.d0 .and. mcon_prim(imsc,2) <= pi .and. mcon_prim(imsc,6) < 1001.d0 )then

            ! we include this mascon in the output
            msc_output(imsc) = msc_output(imsc) + msc_est(imsc)
            total_volume = total_volume + (msc_est(imsc)-gia_dt*msc_gia(imsc))*mcon_prim(imsc,4)   ! EWH x primary mascon area
            total_area = total_area + mcon_prim(imsc,4)
          endif

! Southern Amazon
        else if(region(iregion)(1:8) == "AMAZONAS")then
          if(mcon_region(imsc)(1:8) == "AMAZONAS" .and. mcon_prim(imsc,1)*180.d0/pi > -22.75d0 .and. mcon_prim(imsc,1)*180.d0/pi <= 0.d0 &
             .and. mcon_prim(imsc,6) < 1001.d0 )then

            ! we include this mascon in the output
            msc_output(imsc) = msc_output(imsc) + msc_est(imsc)
            total_volume = total_volume + (msc_est(imsc)-gia_dt*msc_gia(imsc))*mcon_prim(imsc,4)   ! EWH x primary mascon area
            total_area = total_area + mcon_prim(imsc,4)
          endif

! Northern Amazon
        else if(region(iregion)(1:8) == "Amazon_N")then
          if((mcon_region(imsc)(1:8) == "AMAZONAS" .or. mcon_region(imsc)(1:7) == "ORINOCO"  .or. mcon_region(imsc)(1:8) == "SAmerica" )  &
             .and. mcon_prim(imsc,1)*180.d0/pi > 0d0 .and. mcon_prim(imsc,1)*180.d0/pi <= 10.d0 &
             .and. mcon_prim(imsc,6) < 1001.d0 )then

            ! we include this mascon in the output
            msc_output(imsc) = msc_output(imsc) + msc_est(imsc)
            total_volume = total_volume + (msc_est(imsc)-gia_dt*msc_gia(imsc))*mcon_prim(imsc,4)   ! EWH x primary mascon area
            total_area = total_area + mcon_prim(imsc,4)
          endif

! Southern Africa 
        else if(region(iregion)(1:8) == "Africa_S")then
          if((mcon_region(imsc)(1:6) == "Africa" .or. mcon_region(imsc)(1:5) == "CONGO") &
             .and. mcon_prim(imsc,1)*180.d0/pi > -15.5d0 .and. mcon_prim(imsc,1)*180.d0/pi <= -1.5d0 &
             .and. mcon_prim(imsc,6) < 1001.d0 )then

            ! we include this mascon in the output
            msc_output(imsc) = msc_output(imsc) + msc_est(imsc)
            total_volume = total_volume + (msc_est(imsc)-gia_dt*msc_gia(imsc))*mcon_prim(imsc,4)   ! EWH x primary mascon area
            total_area = total_area + mcon_prim(imsc,4)
          endif

! Northern Africa 
        else if(region(iregion)(1:8) == "Africa_N")then
          if((mcon_region(imsc)(1:6) == "Africa" .or. mcon_region(imsc)(1:5) == "CONGO"  &
            .or. mcon_region(imsc)(1:5) == "NIGER" &
            .or. mcon_region(imsc)(1:4) == "NILE" &
            .or. mcon_region(imsc)(1:9) == "LAKE_CHAD" ) &

             .and. mcon_prim(imsc,1)*180.d0/pi > 0.d0 .and. mcon_prim(imsc,1)*180.d0/pi <= 14.5d0 &
             .and. mcon_prim(imsc,6) < 1001.d0 )then

            ! we include this mascon in the output
            msc_output(imsc) = msc_output(imsc) + msc_est(imsc)
            total_volume = total_volume + (msc_est(imsc)-gia_dt*msc_gia(imsc))*mcon_prim(imsc,4)   ! EWH x primary mascon area
            total_area = total_area + mcon_prim(imsc,4)
          endif

! Southern Asia 
        else if(region(iregion)(1:6) == "Asia_S")then
          if((mcon_region(imsc)(1:7) == "Eurasia" .or. mcon_region(imsc)(1:5) == "blaah"  &
            .or. mcon_region(imsc)(1:5) == "blah1" &
            .or. mcon_region(imsc)(1:4) == "blah" &
            .or. mcon_region(imsc)(1:9) == "blah_blah" ) &

             .and. mcon_prim(imsc,1)*180.d0/pi > 8.d0 .and. mcon_prim(imsc,1)*180.d0/pi <= 25.5d0 &
             .and. mcon_prim(imsc,2)*180.d0/pi > 67.d0 .and. mcon_prim(imsc,2)*180.d0/pi <= 109.d0 &
             .and. mcon_prim(imsc,6) < 1001.d0 )then

            ! we include this mascon in the output
            msc_output(imsc) = msc_output(imsc) + msc_est(imsc)
            total_volume = total_volume + (msc_est(imsc)-gia_dt*msc_gia(imsc))*mcon_prim(imsc,4)   ! EWH x primary mascon area
            total_area = total_area + mcon_prim(imsc,4)
          endif

! Eastern Canada 
        else if(region(iregion)(1:8) == "Canada_E")then
          if( mcon_prim(imsc,1)*180.d0/pi > 49.d0 .and. mcon_prim(imsc,1)*180.d0/pi <= 73.d0 &
             .and. mcon_prim(imsc,2)*180.d0/pi > 257.d0 .and. mcon_prim(imsc,2)*180.d0/pi <= 300.d0 &
             .and. mcon_prim(imsc,6) < 1001.d0 )then

            ! we include this mascon in the output
            msc_output(imsc) = msc_output(imsc) + msc_est(imsc)
            total_volume = total_volume + (msc_est(imsc)-gia_dt*msc_gia(imsc))*mcon_prim(imsc,4)   ! EWH x primary mascon area
            total_area = total_area + mcon_prim(imsc,4)
          endif

! Western Canada 
        else if(region(iregion)(1:8) == "Canada_W")then
          if( mcon_prim(imsc,1)*180.d0/pi > 49.d0 .and. mcon_prim(imsc,1)*180.d0/pi <= 73.d0 &
             .and. mcon_prim(imsc,2)*180.d0/pi > 219.d0 .and. mcon_prim(imsc,2)*180.d0/pi <= 257.d0 &
             .and. mcon_prim(imsc,6) < 1001.d0 )then

            ! we include this mascon in the output
            msc_output(imsc) = msc_output(imsc) + msc_est(imsc)
            total_volume = total_volume + (msc_est(imsc)-gia_dt*msc_gia(imsc))*mcon_prim(imsc,4)   ! EWH x primary mascon area
            total_area = total_area + mcon_prim(imsc,4)
          endif

! Eastern USA 
        else if(region(iregion)(1:9) == "America_E")then
          if( mcon_prim(imsc,1)*180.d0/pi > 30.d0 .and. mcon_prim(imsc,1)*180.d0/pi <= 49.d0 &
             .and. mcon_prim(imsc,2)*180.d0/pi > 260.d0 .and. mcon_prim(imsc,2)*180.d0/pi <= 300.d0 &
             .and. mcon_prim(imsc,6) < 1001.d0 )then

            ! we include this mascon in the output
            msc_output(imsc) = msc_output(imsc) + msc_est(imsc)
            total_volume = total_volume + (msc_est(imsc)-gia_dt*msc_gia(imsc))*mcon_prim(imsc,4)   ! EWH x primary mascon area
            total_area = total_area + mcon_prim(imsc,4)
          endif

! Western USA 
        else if(region(iregion)(1:9) == "America_W")then
          if( mcon_prim(imsc,1)*180.d0/pi > 30.d0 .and. mcon_prim(imsc,1)*180.d0/pi <= 49.d0 &
             .and. mcon_prim(imsc,2)*180.d0/pi > 230.d0 .and. mcon_prim(imsc,2)*180.d0/pi <= 250.d0 &
             .and. mcon_prim(imsc,6) < 1001.d0 )then

            ! we include this mascon in the output
            msc_output(imsc) = msc_output(imsc) + msc_est(imsc)
            total_volume = total_volume + (msc_est(imsc)-gia_dt*msc_gia(imsc))*mcon_prim(imsc,4)   ! EWH x primary mascon area
            total_area = total_area + mcon_prim(imsc,4)
          endif

! Alaska
        else if(region(iregion)(1:6) == "Alaska")then
          if( mcon_prim(imsc,1)*180.d0/pi > 58.d0 .and. mcon_prim(imsc,1)*180.d0/pi <= 72.d0 &
             .and. mcon_prim(imsc,2)*180.d0/pi > 192.d0 .and. mcon_prim(imsc,2)*180.d0/pi <= 219.d0 &
             .and. mcon_prim(imsc,6) < 1001.d0 )then

            ! we include this mascon in the output
            msc_output(imsc) = msc_output(imsc) + msc_est(imsc)
            total_volume = total_volume + (msc_est(imsc)-gia_dt*msc_gia(imsc))*mcon_prim(imsc,4)   ! EWH x primary mascon area
            total_area = total_area + mcon_prim(imsc,4)
          endif


! to calculate just the GIA signal everywhere ...
        else if (region(iregion)(1:7) == "giaonly")then
          msc_output(imsc) = gia_dt*msc_gia(imsc)
          total_volume = 0.d0
          total_area = total_area + mcon_prim(imsc,4)

! everything else is done using the mcon_region descriptor from the mascon_file
        else if(trim(mcon_region(imsc)) == trim(region(iregion)) )then
          ! we include this mascon in the output
          msc_output(imsc) = msc_output(imsc) + msc_est(imsc)
          total_volume = total_volume + (msc_est(imsc)-gia_dt*msc_gia(imsc))*mcon_prim(imsc,4)   ! EWH x primary mascon area
          total_area = total_area + mcon_prim(imsc,4)
        endif
      enddo

      ! now, calculate the EWH over the region
      region_EWH = -1.d0*total_volume/total_ocean_area
      write(message,'(a,f15.6,a,a,i6,i3.2)')"Eustatic sea level change of ",region_EWH*1.d3,' mm from region: ',region(iregion),yr,month
      call status_update('STATUS','UTIL','writeout_mascons',' ',message,0)

      ! now, apportion this onto each ocean mascon by multiplying it by the fingerprint value
      allocate(Plm(0:maxdeg,0:maxdeg))
      do imsc = 1,max_prim
        if(mcon_prim(imsc,6) > 1000.d0)then  ! it is an ocean mascon
          call legendre_matrix(maxdeg,dcos(pi/2.d0-mcon_prim(imsc,1)),Plm)
          fprint_value = 0.d0
! PT200810: all degrees or just from deg 2 ? Deg 0 is added lower down
          do ideg=1,maxdeg
            do iord = 0,ideg
!DEBUG
!if(ideg < 3)print*,ideg,iord,C_fprint(ideg,iord),S_fprint(ideg,iord)
              fprint_value = fprint_value + Plm(ideg,iord)*( C_fprint(ideg,iord)*dcos(dble(iord)*mcon_prim(imsc,2)) &
                                                            + S_fprint(ideg,iord)*dsin(dble(iord)*mcon_prim(imsc,2)) )
            enddo
          enddo
!          msc_output(imsc) = msc_output(imsc) + (fprint_value  +  C_fprint(0,0))  ! include degree 0 here if you want
          msc_output(imsc) = msc_output(imsc) + region_EWH*(fprint_value  +  C_fprint(0,0))  ! include degree 0 here if you want
!if(mcon_prim(imsc,6) > 1000.d0)print*,'imsc,fprint_value,region_EWH,msc_output(imsc)',imsc,fprint_value,region_EWH,msc_output(imsc)
        endif
      enddo

      ! deallocate arrays before looping for the next region
      deallocate(Plm)
      deallocate(C_fprint)
      deallocate(S_fprint)

    else
print*,'not including region: ',region(iregion)
    endif


  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! finally, write out all the mascon values
  open(luout,file=output_file,status='unknown')

  do imsc=1,max_prim
!    print*,imsc,msc_output(imsc)," ",mcon_region(imsc),' output EWH values'
    write(prmnam,'(f6.0,a3,i5.5,a4)')dble(imsc)+24.d0," MC",imsc," (m)            "
    write(luout,'(a30,f17.5,f12.5,f19.5,f12.5,f9.1)') prmnam,msc_apr(imsc),msc_est(imsc),msc_output(imsc) &
          ,msc_sigma(imsc),0.d0

!if(mcon_region(imsc)(1:6) == "Africa" .and.msc_output(imsc) /= 0.d0 )print*,imsc,msc_output(imsc)
  enddo
  write(luout,'(a)')"MASCON CONSTRAINT SIGMAS"
  close(luout)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end



