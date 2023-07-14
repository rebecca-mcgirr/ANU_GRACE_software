  program covar_vs_dist

! program to calculate (for Srinivas) the relation between covariances between mascons and their separation
!
! P Tregoning
! 25 October 2019

  use mascon_mod

  implicit none

  integer max_num_prime

  parameter (max_num_prime = 5000)

  character*150 :: mascon_file,vcv_file,line
  integer*4,parameter     :: lumsc_in=10,luvcv_in=11

  real(kind=8), allocatable :: msc_dist(:,:),VCV(:,:)
  integer*4 :: imsc,imsc2
  real(kind=8) :: mindist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print*,"covar_vs_dist ~/gg/grace/tables/mascons_stage4_V003a addnorm.vcv"
  call getarg(1,mascon_file)
  open(lumsc_in,file=mascon_file)  
  call getarg(2,vcv_file)  
  open(luvcv_in,file=vcv_file)  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call read_msc_hdr(lumsc_in,mascon_file,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                                ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
  call allocate_mascon_arrays
! PT180220: using trimlen here doesn't work - the subroutine is expecting a C*150, so pass the whole thing through
!  call read_mascon_file(lumsc_in,mascon_file(1:trimlen(mascon_file)))
  call read_mascon_file(lumsc_in,mascon_file)

! convert the mascon lat/lon/rad coords into cartesian xyz (quicker than doing it later at each epoch ...)
  call calc_mascon_xyz("N")
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate the distance between each primary mascon and each other
  allocate(msc_dist(max_prim,max_prim))
  do imsc = 1,max_prim
    do imsc2 = 1,max_prim
      msc_dist(imsc,imsc2) = dsqrt( (mcon_prim(imsc,1)-mcon_prim(imsc2,1))**2+(mcon_prim(imsc,2)-mcon_prim(imsc2,2))**2 &
                                    +(mcon_prim(imsc,3)-mcon_prim(imsc2,3))**2 )
!print*,imsc,imsc2,msc_dist(imsc,imsc2)
    enddo
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read in the VCV matrix
  allocate(VCV(max_prim,max_prim))
  ! read the file until the VCV header line is found
  line=""
  do while (line(1:14) /= "  VCV SOLUTION")
    read(luvcv_in,'(a)')line
  enddo

  ! now read in the lower-triangular VCV
  print*,'reading the VCV matrix'
  do imsc=1,max_prim
    read(luvcv_in,*)(VCV(imsc,imsc2),imsc2=1,imsc)
    VCV(1:imsc,imsc) = VCV(imsc,1:imsc)
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now, output the covariance between mascons as a function of distance. Limit the output to only mascons within 2000 km
  do imsc=1,max_prim
    do imsc2=1,max_prim
      if(msc_dist(imsc2,imsc) < 5000.d3)then
        if(imsc /= imsc2)then
            if( mcon_prim(imsc,4) <  mcon_prim(imsc2,4))then
              mindist = dsqrt(mcon_prim(imsc,4)/1.e6)
            else
              mindist = dsqrt(mcon_prim(imsc2,4)/1.e6)
            endif
            print*,msc_dist(imsc2,imsc)/1.d3,VCV(imsc2,imsc)/( dsqrt(VCV(imsc,imsc))*dsqrt(VCV(imsc2,imsc2)) ) &
                               ,mindist,dsqrt(mcon_prim(imsc,4)/1.e6),dsqrt(mcon_prim(imsc2,4)/1.e6) &
                               ,imsc,imsc2,mcon_prim(imsc,6),mcon_prim(imsc2,6) &
                               ," dist_vs_covar"

!print*,VCV(imsc2,imsc),VCV(imsc,imsc),VCV(imsc2,imsc2)
        endif
      endif
    enddo
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  end

