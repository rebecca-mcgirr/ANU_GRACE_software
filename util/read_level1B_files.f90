  subroutine read_level1B_files(luGNV,luSCA,grace_sec,sat_efixd,quat_inert,ngaps)

! subroutine to read to the next common epoch in the GNV1B and SCA1B files. We need also to count how many epochs are skipped, then
! use that info to mask out those accobs and also to integrate the panel temperatures over the gap period.
!
! P. Tregoning
! 27 November 2015

  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! passed variables
  integer*4,   intent(in)    :: luGNV,luSCA
  real(kind=8),intent(inout) :: grace_sec          ! grace seconds of new common epoch
  real(kind=8),intent(out)   :: sat_efixd(6)       ! satellite pos/vec at new common epoch
  real(kind=8),intent(out)   :: quat_inert(4)      ! srf to inert quaternion at new common epoch
  integer*4   ,intent(out)   :: ngaps              ! number of epochs in the gap from the previous epoch

! local variables
  real(kind=8)  :: grace_sec_old,grace_sec_gnv,grace_sec_sca
  character*1   :: sat,ef_or_I
  real(kind=8)  :: junk(3)
  integer*4     :: which_camera
  real(kind=8)  :: diff_gnv,diff_sca
  integer*4     :: i,ioerr,ioerr2
  character*100 :: message
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



! store the previous epoch
  grace_sec_old = grace_sec

! read the next new epoch in the GNV1B file
  read(luGNV,*,iostat=ioerr,end=1000)grace_sec_gnv,sat,ef_or_I,sat_efixd(1:3),junk,sat_efixd(4:6)
! calculate the epoch difference (in seconds) from the previous epoch
  diff_gnv = grace_sec_gnv - grace_sec_old

! now read the next SCA1B epoch
  read(luSCA,*,iostat=ioerr2,end=1001)grace_sec_sca,sat,which_camera,quat_inert
! calculate the epoch difference (in seconds) from the previous epoch
  diff_sca = grace_sec_sca - grace_sec_old

! are the differences the same?
  if(diff_gnv == diff_sca)then
    ngaps = int(diff_gnv - diff_sca)/5
    grace_sec = grace_sec_gnv
  elseif (diff_sca > diff_gnv) then
!   we need to read more epochs of the GNV1B file to match the SCA epoch
    ngaps = int(diff_sca-diff_gnv)/5
    write(message,'(a,i3,a,i10,a,i9,a)')'Missing ',ngaps,' epochs in SCA file (' &
              ,int(grace_sec_gnv),' - ',int(grace_sec_gnv+(ngaps-1)*5),' )'
    call status_update('STATUS','UTIL','read_level1B_files',' ',message,0)
    do i=1,ngaps
      read(luGNV,*,iostat=ioerr,end=1000)grace_sec_gnv,sat,ef_or_I,sat_efixd(1:3),junk,sat_efixd(4:6)
    enddo
    grace_sec = grace_sec_gnv
  else
!   we need to read more epochs of the SCA1B file to match the GNV epoch
    ngaps = int(diff_gnv-diff_sca)/5
    write(message,'(a,i3,a,i9,a,i9,a)')'Missing ',ngaps,' epochs in GNV file (' &
               ,int(grace_sec_sca),' - ',int(grace_sec_sca+(ngaps-1)*5),' )'
    call status_update('STATUS','UTIL','read_level1B_files',' ',message,0)
    do i=1,ngaps
      read(luSCA,*,iostat=ioerr2,end=1001)grace_sec_sca,sat,which_camera,quat_inert
    enddo
    grace_sec = grace_sec_sca
  endif

  return

1000 print*,'Error: reached end of GNV1B file'
     stop

1001 print*,'Error: reached end of SCA1B file'
     stop


  end


