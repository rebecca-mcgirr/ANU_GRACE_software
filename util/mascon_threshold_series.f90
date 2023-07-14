  program mascon_threshold_series
  
! program to read in a vcv file and output non-zero values for only mascons that received a sufficiently large adjustment
! Adapted from mascon_threshold, but this version checks whether mascons have EVER exceeded
! the threshold in a given time series of vcv files
!
! P. Tregoning
! 3 February 2023
!
! PT230214: add the capability to conserve mass, distributing the apriori continental mass loss over the oceans using
!           fingerprints calculated by Tony using CALSEA

  use mascon_mod          ! defines all the mascon arrays

  implicit none

  character :: list_file*150, output_stem*50, vcv_files(1000)*150, arg*150
  character*150 :: mascon_file
  character*1   :: zero_ocean
  real(kind=8) :: adjust_size                 ! |adjust_size| is the threshold for a mascon to have a non-zero output value
  real(kind=8) :: incl_adjust_percent         ! percentage of adjustment to add to apriori value
  real(kind=8),allocatable :: msc(:,:,:)      ! mascon apriori, adjustment and estimate values
  character*81,allocatable :: mascon_lines(:)

  integer*4 :: imsc,nmascon_t,ioerr,iVCV,n_vcv,indx
  character*130 :: line,output_name
  real(kind=8),allocatable :: msc_densities(:)
  logical max_thresh,min_thresh
  real(kind=8) :: msc_max,msc_min

! variables related to reading and storing the fingerprint values
  integer*4    :: nmsc_ocean,nmsc_land
  real(kind=8),allocatable :: fingerprints(:,:)
  integer*4    :: tmp_msc
  real(kind=8) :: tmp_fing(10000)
  
! unit numbers
  integer*4,parameter :: lumsc = 12, luout = 11, luin=10, lulist = 13,lufing = 14

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! decode command line
  call getarg(1,list_file)
  if (list_file(1:1) == "")then
    print*,"max_adjust list_file_in output_stem adjustment_size %adjust_to_include mascon_file zero_ocean [y/n]"
    stop
  endif

  call getarg(2,output_stem)
  call getarg(3,arg)
  read(arg,*)adjust_size
  call getarg(4,arg)
  read(arg,*)incl_adjust_percent
  call getarg(5,mascon_file)
  call getarg(6,zero_ocean)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! PT230203: open and read the list_file to get the names of all the vcv files to be read
  open(lulist,file=trim(list_file),status='old',iostat=ioerr)
  if(ioerr == 0)then
    ! read how many vcv files there are in the list_file
    n_vcv = 0
    do while (ioerr == 0)
      read(lulist,'(a)',end=1000,iostat=ioerr)arg
      if(ioerr == 0)then
        n_vcv = n_vcv + 1
        vcv_files(n_vcv) = arg
      endif
    enddo
1000 print*,'there are ',n_vcv,' VCV file in the series to be assessed'
    close(lulist)
  else
    print*,'error opening list file'
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! open and read the fingerprint file
  open(lufing,file="mascons.fingerprint",status='old')
  ! read the primary numbers of the ocean mascons
  nmsc_ocean = 9072                         ! 9072 ocean mascons in mascon file mascons_stage5_V006_200km
  read(lufing,*)msc_ocean(1:nmsc_ocean)
  ! now, read all the ascii lines of the land mascons
  nmsc_land = 3684
  allocate(fingerprints(nmsc_land,nmsc_ocean))
  do imsc = 1,nmsc_land
    read(lufing,*)tmp_msc,tmp_fing(1:nmsc_ocean)
    fingerprints(tmp_msc,1:nmsc_ocean) = tmp_fing(1:nmsc_ocean)
  enddo
  close(lufing)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! open the first VCV file to get the number of mascons
  open(luin,file=trim(vcv_files(1)),status='old')

! read and transfer the header part until we get to the first mascon entry
  line = " "
  do while (line(7:9) /= " MC")
    read(luin,'(a)')line
    ! get the number of mascons
    if(line(1:20) == "Solution with mascon")read(line(103:108),*)nmascon_t
    if(line(8:9) /= " MC")write(luout,'(a)')line
  enddo

! allocate the mascon arrays
  allocate(mascon_lines(nmascon_t))
  allocate(msc(nmascon_t,3,n_vcv))  ! PT230203: third dimension is the number of VCV files (i.e. epochs)
  allocate(msc_densities(nmascon_t))

  close(luin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! loop through and read ALL the mascon values from each VCV file in the list_file
  do iVCV = 1,n_vcv
    ! open the first VCV file to get the number of mascons
    open(luin,file=trim(vcv_files(iVCV)),status='old')
    line = " "

    ! read and transfer the header part until we get to the first mascon entry
    line = " "
    do while (line(7:9) /= " MC")
      read(luin,'(a)')line
      ! get the number of mascons
      if(line(1:20) == "Solution with mascon")read(line(103:108),*)nmascon_t
    enddo
    
    ! we are now at the first line of the mascon entries
    backspace(luin) 
    do imsc = 1,nmascon_t
      read(luin,'(a81)')mascon_lines(imsc)
      read(mascon_lines(imsc)(32:47),*)msc(imsc,1,iVCV)        !   apriori mascon value
      read(mascon_lines(imsc)(48:64),*)msc(imsc,3,iVCV)        ! estimated mascon value
      msc(imsc,2,iVCV)  = msc(imsc,3,iVCV) - msc(imsc,1,iVCV)  ! adjustment to apriori mascon value
    enddo
    
    close(luin)
  enddo
  print*,'Have read in mascon apr/adj/est for all ',n_vcv,' VCV files'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the mascon file 
  call read_msc_hdr(lumsc,mascon_file,msc_hdr_lines,max_hdr_lines,n_hdr_lines,msc_hdr_code,max_prim,max_sec,max_tern &
                            ,max_tern_per_prim,max_tern_per_sec,max_sec_per_prim)
  call allocate_mascon_arrays
  ! PT180220: using trimlen here doesn't work - the subroutine is expecting a C*150, so pass the whole thing through
  !  call read_mascon_file(lumsc_in,mascon_file(1:trimlen(mascon_file)))
  call read_mascon_file(lumsc,mascon_file)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! loop through the mascons and:
! - determine if the mascon adjustment ever exceeded the |threshold|
! -  set the estimated value to the apriori value if the adjustment was always within the threshold input by the user
! -  otherwise, change for each epoch the estimated to be the apriori + adjustment*percentage input by the user
  print*,'Looping over ',nmascon_t,' mascons'
  do imsc = 1,nmascon_t

    max_thresh = .false.
    min_thresh = .false.

    !if(mod(imsc,100) == 0)print*,'mascon: ',imsc
    
    ! check the range of adjustments and see whether it ever exceeds the threshold
    msc_max = maxval(msc(imsc,2,:))
    msc_min = minval(msc(imsc,2,:))
    if (msc_max < adjust_size .and. dabs(msc_min) < adjust_size)then
      msc(imsc,2,:) = 0.d0

    else if (zero_ocean == "Y" .and. mcon_prim(imsc,6) > 1010.d0)then
      msc(imsc,2,:) = 0.d0    
       
    else if (msc_max > adjust_size .or. dabs(msc_min) > adjust_size)then
      ! update the adjustment for this mascon for ALL epochs
      do iVCV = 1, n_vcv
        ! PT230206: check whether this is the first epoch that exceeds the dabs(adjust) threshold
        if( msc(imsc,2,iVCV) > adjust_size .and. .not. max_thresh)then
          print*,trim(vcv_files(iVCV)),' mascon:',imsc,': first occurrence above threshold',adjust_size,msc(imsc,2,iVCV)
          max_thresh = .true.
        else if (msc(imsc,2,iVCV) < -1.d0*adjust_size .and. .not. min_thresh)then
          print*,trim(vcv_files(iVCV)),' mascon:',imsc,': first occurrence below threshold',-1.d0*adjust_size,msc(imsc,2,iVCV)
          min_thresh = .true.
        endif
        msc(imsc,2,iVCV) =  msc(imsc,2,iVCV) * incl_adjust_percent / 100.d0
      enddo
      
    endif
  enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate the new mascon values, given the adjustments
  do iVCV = 1,n_vcv

    ! loop over mascons to output adjustments
    do imsc = 1,nmascon_t
      ! add the adjustment to the apriori
      msc(imsc,3,iVCV) = msc(imsc,1,iVCV) + msc(imsc,2,iVCV) * incl_adjust_percent / 100.d0
      ! apportion the land mass change over the ocean, using the fingerprints if it is a land mascon
      if(mcon_prim(imsc,6) < 1010.d0 )then
        do iocean = 1,nmsc_ocean
          msc(iocean,3,iVCV) = msc(iocean,3,iVCV) + (-1.d0*msc(imsc,3,iVCV))*fingerprints(imsc,iocean)
        enddo
      endif
    enddo
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  OUTPUT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
    ! estimate = apriori + adjustment*percentage
  do iVCV = 1,n_vcv
    
    ! transfer the header information to the output VCV file
    open(luin ,file=trim(vcv_files(iVCV)),status='old')
    ! PT230209: assume that the input vcv files are named mscA_YYYY-MM-*. Take the YYYY and MM and replace the rest
    indx = len(trim(output_stem))
    output_name(1:12+indx+4) = "aprA_"//vcv_files(iVCV)(6:12)//trim(output_stem)//".vcv"
    print*,"writing file: ",trim(output_name)
    open(luout,file=trim(output_name),status='unknown')
    line = " "
    do while (line(7:9) /= " MC")
      read(luin,'(a)')line
      if(line(7:9) /= " MC")write(luout,'(a)')line
    enddo
    close(luin)

    ! loop over mascons to output adjustments
    do imsc = 1,nmascon_t
      ! add the adjustment to the apriori
      msc(imsc,3,iVCV) = msc(imsc,1,iVCV) + msc(imsc,2,iVCV) * incl_adjust_percent / 100.d0
    enddo

      
      ! write it to the output file
    do imsc = 1,nmascon_t
      write(mascon_lines(imsc)(31:65),'(2f17.7)')msc(imsc,1,iVCV),msc(imsc,3,iVCV)
      write(luout,'(a)')mascon_lines(imsc)
    enddo

    close(luout)
  enddo

! finish off the vcv file
!  write(11,'(a)')"  VCV SOLUTION"


  end




