  program lost_ternarys

! program to do a very simple read of a mascon file and identify which, if any, ternarys are missing from the file
!
! P. Tregoning
! 7 March 2019

  implicit none

  character*100 :: mascon_file,line
  character*7   :: code,tmp_code
  integer*4     :: lumsc,nprim,nsec,ntern,sum_terns,terns_in_prim
  integer*4, allocatable :: tern_vector(:)
  integer*4     :: imsc,itern,tmp_tern,tmp_prim,tmp_sec,tern_number

  lumsc = 10

  call getarg(1,mascon_file)
  open(lumsc,file=mascon_file)

! read the number of primary and ternary mascons supposed to be in the file
  read(lumsc,*)code,nprim,nsec,ntern
  print*,mascon_file(1:50)," Line 1:  ",code,nprim,nsec,ntern

! define the ternary vector
  allocate(tern_vector(ntern))
  tern_vector = -999

! skip over the header lines
  line = "#"
  do while (line(1:1) == "#" )
    read(lumsc,'(a)')line
  enddo
  backspace(lumsc)

! now, read each primary to get prim number and the number of ternarys supposed to be in the primary
  sum_terns = 0
  do imsc = 1,nprim
    read(lumsc,'(a)')line

    if(line(8:10) /= "  P")then  ! it is NOT a primary line, but it should be!!!
      print*,'Error. Line is not a primary mascon but should be'
      print*,line
      stop

    else

      read(line(24:31),*)tmp_tern
!      print*,'Primary mascon',imsc,' should have ',tmp_tern,' ternarys'
 
      read(lumsc,'(a)')line  ! skip the secondary line
      do itern = 1,tmp_tern
        read(lumsc,'(i7,a3)')tern_number,tmp_code
 
        if(tmp_code == "  T")then  ! it is a ternary line
          tern_vector(tern_number) = imsc
          sum_terns = sum_terns + 1
        else
          print*,'problem reading ternary',itern,' in primary',imsc  !,'. There should be ',tmp_tern,'ternarys in this mascon'
          stop
        endif

      enddo

    endif
  enddo

! check whether there were enough ternarys
  if(sum_terns == ntern)then
    print*,'all ternarys present and accounted for in file: ',mascon_file
  else
    print*,'missing',ntern-sum_terns,' ternarys ...'

    do itern = 1,ntern
      if(tern_vector(itern) < 0)then
        print*,itern-1,': Primary mascon',tern_vector(itern-1)
        print*,itern,': missing'
        print*,itern+1,': Primary mascon',tern_vector(itern+1)
        print*,' '
      endif
    enddo




  endif



  end




