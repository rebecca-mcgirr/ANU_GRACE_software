  program extract_mascons

! program to read a vcv file and output a vcv file with non-zero mascons only for specific regions:
! 1. Antarctica
! 2. Greenland
! 3. Antarctica AND Greenland
!
! that's all I've coded for now ...
!
! P. Tregoning
! 29 May 2019
!
! MODS
! PT190906: output any mascons with |adjustments| > 0.8 m

  use mascon_mod

  implicit none

  character*250 :: message
  character*150 :: arg,infile_vcv,outfile_vcv,mascon_file,line
  character*1   :: flag
  integer*4     :: lu_msc
  parameter (lu_msc = 12)
  integer*4     :: imsc,nmsc
  real(kind=8)  :: vcv_EWH
  real(kind=8)  :: pi

  pi = 4.d0*datan(1.d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! decode the runstring
  call getarg(1,infile_vcv)
  if (infile_vcv(1:1) == " ")then
    print*,"Runstring: extract_mascons msc_iter3.vcv msc_AntGrn.vcv mascon_file flag"
    print*,"  where flag = 1 (Antarctica); 2 (Greenland); 3 (Ant+Grn) ; 4 (|adj| > 0.8m ; 5 (as for 4 but polar mascons to -15 m)&
          , ; 6 (everything but zero for oceans)"
    stop
  else
    open(10,file=infile_vcv,status='old')
  endif
  call getarg(2,outfile_vcv)
  open(11,file=outfile_vcv,status='unknown')
  call getarg(3,mascon_file)
  call getarg(4,flag)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the mascon file
  call read_msc_hdr(lu_msc,mascon_file,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                                ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)

! establish the information to identify each ternary mascon from lat/lon 
  call tern_lat_bands_ell(ternary_lat_spacing/60.d0)

! allocate the array sizes for the mascons
  call allocate_mascon_arrays

! read in the mascon information
  call read_mascon_file(lu_msc,mascon_file)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the input vcv file and transfer the header to the output file
  line = " "
  do while (line(7:9) /= " MC")
    read(10,'(a)')line
    if(line(7:9) /= " MC")write(11,'(a)')line
    if(line(1:31) == "Number of estimated parameters:")then
      read(line(40:47),'(i7)')nmsc
      write(message,'(a,i7)')"Number of mascons in VCV file:",nmsc
      call status_update('STATUS','UTIL','extract_mascons',mascon_file,message,0)
    endif
  enddo
  backspace(10)

! now read the mascon solution values
  do imsc = 1,nmsc
    read(10,'(a)')line
    ! check whether this primary mascon is in a region that we want to use a non-zero value
    if(flag == "1" .or. flag == "3")then
      if (mcon_prim(imsc,1) < -60.d0*pi/180.d0 .and. mcon_prim(imsc,6) < 1010.d0)then    ! it is a land mascon in AntarcticaA
        write(11,'(a)')line
      else if(mcon_prim(imsc,1) > 58.d0*pi/180.d0 .and. (mcon_prim(imsc,2) > 300*pi/180.d0 &
                     .and. mcon_prim(imsc,2) < 338.63d0*pi/180.d0 .and. mcon_prim(imsc,6) < 1010.d0))then   ! it is in Greenland
        write(11,'(a)')line
      else
        write(line(52:66),'(f15.9)')0.0000000
        write(11,'(a)')line
      endif
    else if (flag == "2" ) then    
      if(mcon_prim(imsc,1) > 58.d0*pi/180.d0 .and. (mcon_prim(imsc,2) > 260*pi/180.d0 &
                .and. mcon_prim(imsc,2) < 338.63d0*pi/180.d0))then   ! it is in Greenland
        write(11,'(a)')line
      else
        write(line(52:67),'(f15.9)')0.0000000
        write(11,'(a)')line
      endif

! PT190906: output any mascon with an |adjustment| > 0.8 m
    else if (flag == "4" )then
      read(line(52:66),*)vcv_EWH
      if(dabs(vcv_EWH) > 0.8d0 .or. prim_flags(imsc)(1:5) == 'PCasp')then
print*,'Assign non-zero value of',vcv_EWH,' to mascon ',prim_flags(imsc),imsc,mcon_prim(imsc,1:2)*180.d0/pi 
        write(11,'(a)')line
      else
        write(line(52:67),'(f15.9)')0.0000000
        write(11,'(a)')line
      endif

! Pt190909: put -15 m apriori mascons on polar mascons with > |0.8 m| adjustments, to make it extreme!
    else if (flag == "5")then
      read(line(52:66),*)vcv_EWH
      if (mcon_prim(imsc,1) < -60.d0*pi/180.d0 .and. mcon_prim(imsc,6) < 1010.d0 .and. &
          (vcv_EWH < -0.8d0 .or. vcv_EWH > 0.8d0) )then    ! it is a land mascon in AntarcticaA
        write(line(52:67),'(f15.9)')-15.0000000
        write(11,'(a)')line
      else if(mcon_prim(imsc,1) > 58.d0*pi/180.d0 .and. (mcon_prim(imsc,2) > 300*pi/180.d0 &
                     .and. mcon_prim(imsc,2) < 338.63d0*pi/180.d0 .and. mcon_prim(imsc,6) < 1010.d0) .and. &
                      vcv_EWH < -0.8d0)then   ! it is in Greenland
        write(line(52:67),'(f15.9)')-15.0000000
        write(11,'(a)')line
      else
        write(line(52:66),'(f15.9)')0.0000000
        write(11,'(a)')line
      endif

    else if (flag == "6")then
      read(line(52:66),*)vcv_EWH
      if (mcon_prim(imsc,6) > 1010.d0 )then
        write(line(52:66),'(f15.9)')0.0000000   ! set the ocean EWH adjustment to zero
      endif
      write(11,'(a)')line
    endif


  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call status_update('STATUS','UTIL','extract_mascons',' ','End of extract_mascons',0)
  end

