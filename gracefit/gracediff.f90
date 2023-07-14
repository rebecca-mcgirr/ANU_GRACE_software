!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                      !
!  PROGRAM:   GRACEDIFF   - difference two GTORB orbits and output pos/vel differences !
!                                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  program gracediff

! program to read in two binary GTORB files and compute the difference in position and velocity.
! Output will be possible in XYZ, along/cross/radial
!
! P. Tregoning
! 15 November 2013
!
! MODS:
! PT140515: change the output to include partials of bias/scale/pos/vel if requested
! PT170609: define GTORB record length information through mod_subs/gtorb_mod
! PT190328: why isn't LUGT declared anymore?

  use gracefit_mod
  use gtorb_mod

  implicit none

!  include 'gracefit.h90'
!  include 'input.h90'

  character*260   :: orbfiles(2)                      ! input GTORB files
  character*260  :: message                           ! for output using status_update()   
  integer*4      :: i,j,isat,ioerr,iepoch
  real(kind=8)   :: dposvel(6),d_rac(6)               ! delta pos/vel, delta pos/vel in along/cross/rad
  real(kind=8)   :: range_rate                        ! range rate at each epoch (works if a GRACE A and GRACE B orbit are input)
  real(kind=8)   :: conv(6)
  character*10   :: arg
  character*3    :: out_type                          ! output in along/cross/radial (ACR) or north/east/up (NEU)
  character*4    :: partial_type                      ! type of partial to output. They will be dX/dP, ... dZdot/dP where P = X(YZ), velX(YZ), Bsx(yz) or Sclx(yz)

! variables for reading the binary headers
  integer*4        :: LUGT(2)
  integer*4        :: irec,rec_length(2)
  double precision :: satics_t(maxrvecprm+1,2)                   ! A priori values for IC's
  double precision :: apr_ScaleBias(max_SB,maxsat)               ! A priori values for scale and bias
  double precision :: GPSantoff(7,maxsat)                        ! Quaternion and offsets for antenna 

! variables for reading data from the GTORB files
  double precision, allocatable :: rvec(:,:,:)               ! Vector of positions, velocities and partials for "truth" GTORB
  double precision, allocatable :: apr_prm(:)                    ! A priori values for parameters
  double precision, allocatable :: sciframe_quat(:,:,:)          ! Quaternions for each satellite for every epoch
  double precision, allocatable :: srf2trf_rotmat(:,:,:,:)       ! SRF to TRF rotation matrix for GRACE A and B
  double precision, allocatable :: srf2trf_deriv_rotmat(:,:,:,:) ! Differentiated SRF to TRF rotation matrix for GRACE A and B
  real(kind=8),   allocatable   :: epoch(:,:)
  integer*4                     :: tmp_epoch

! variables for converting to NEU, ACR, or whatever else
  real(kind=8)  :: llh(3)
  real(kind=8)  :: rotmat(3,3)
  real(kind=8)  :: d_out(6)

! PT140407: extract from the GTORB files the partials of XYZ, XYZdot wrt some parameter
  real(kind=8), allocatable :: all_partials(:,:,:)
  real(kind=8), allocatable :: Xpart(:,:,:), Ypart(:,:,:), Zpart(:,:,:), Xdotpart(:,:,:),Ydotpart(:,:,:),Zdotpart(:,:,:)
  integer*4 :: part_offset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! decode runstring
  call getarg(1,orbfiles(1))
  print *, orbfiles(1)
  if (orbfiles(1)(1:1) == '') then
    print*,"gracediff GTORB_2010-09-11_00-00-00_A_02.bin GTORB_2010-09-11_00-00-00_A_02.bin.truth nepochs [partial]"
    print*,"    where partial is a char code for the partials dX/dP .. dZdot/dp where P is X,Y,Z,Xdot,... Bsx,Bsy,.. Sclx, scly .."
    stop
  endif

  call getarg(2,orbfiles(2))
  call getarg(3,arg)
  read(arg,*)nepochs_t

! PT140515: use the 4th argument as the type of partial that we want to use. Therefore, hardwire the output type to XYZ
  out_type = "XYZ"
  partial_type = "    "
  call getarg(4,partial_type)
  if(partial_type     == "X   ") then   ! we want to output partials of pos/vel wrt Xo
    part_offset = 0
  elseif(partial_type == "Y   ") then   ! we want to output partials of pos/vel wrt Yo
    part_offset = 6
  elseif(partial_type == "Z   ") then   ! we want to output partials of pos/vel wrt Zo
    part_offset = 12
  elseif(partial_type == "Xdot") then   ! we want to output partials of pos/vel wrt Xdot_o
    part_offset = 18
  elseif(partial_type == "Ydot") then   ! we want to output partials of pos/vel wrt Ydot_o
    part_offset = 24
  elseif(partial_type == "Zdot") then   ! we want to output partials of pos/vel wrt Zdot_o
    part_offset = 30
  elseif(partial_type == "Bsx ") then   ! we want to output partials of pos/vel wrt Bsx
    part_offset = 54
  elseif(partial_type == "Bsy ") then   ! we want to output partials of pos/vel wrt Bsy
    part_offset = 60
  elseif(partial_type == "Bsz ") then   ! we want to output partials of pos/vel wrt Bsz
    part_offset = 66
  elseif(partial_type == "Sclx") then   ! we want to output partials of pos/vel wrt Sclx
    part_offset = 36
  elseif(partial_type == "Scly") then   ! we want to output partials of pos/vel wrt Scly
    part_offset = 42
  elseif(partial_type == "Sclz") then   ! we want to output partials of pos/vel wrt Sclz
    part_offset = 48
  endif

  print*,'partial requested :',partial_type,' part_offset',part_offset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dimension the problem 
! PT170609: make it bigger than we think necessary
  nmascons_t = 10000
!  nepochs_t = 17040
  nsat_t = 2
  nparam_t = 4574
  allocate(rvec(6,nsat_t,nepochs_t))
  allocate(sciframe_quat(4,nsat_t,nepochs_t))
  allocate(srf2trf_rotmat(3,3,nepochs_t,nsat_t))
  allocate(srf2trf_deriv_rotmat(3,3,nepochs_t,nsat_t))
  allocate(apr_prm(nparam_t))
  allocate(epoch(nsat_t,nepochs_t))
  allocate(all_partials(72,2,nepochs_t))
  allocate(Xpart(6,nsat_t,nepochs_t))
  allocate(Ypart(6,nsat_t,nepochs_t))
  allocate(Zpart(6,nsat_t,nepochs_t))
  allocate(Xdotpart(6,nsat_t,nepochs_t))
  allocate(Ydotpart(6,nsat_t,nepochs_t))
  allocate(Zdotpart(6,nsat_t,nepochs_t))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! conversion factors (want output in m and mm
  conv(1) = 1.d0
  conv(2) = 1.d0
  conv(3) = 1.d0
  conv(4) = 1.d3
  conv(5) = 1.d3
  conv(6) = 1.d3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT170609: determine the correct record length. It is returned in file_recl
  LUGT(1) = 13
  LUGT(2) = 14
  call GTORB_record_length(LUGT(1),'GRACEDIFF',orbfiles(1),rec_length(1))
  call input_openFile(LUGT(1),orbfiles(1),'old','GRACEDIFF','gracediff',"direct    ",rec_length(1),'FATAL',&
                          'Error opening the first GTORB file: ')
  call GTORB_record_length(LUGT(2),'GRACEDIFF',orbfiles(2),rec_length(2))
  call input_openFile(LUGT(2),orbfiles(2),'old','GRACEDIFF','gracediff',"direct    ",rec_length(2),'FATAL',&
                          'Error opening the second GTORB file: ')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




! now just read in the position/velocity for the entire file
  do isat = 1,2
! PT131121: read the 13th line to find out the models used and the type of mascon files
    if(rec_length(isat) == 318788)then
      read(LUGT(isat),rec=13)message(1:170)
    else
      read(LUGT(isat),rec=14)message(1:170)
    endif
    call status_update('STATUS','GRACEDIFF','gracediff',orbfiles(isat),message,0)

! PT140618: irec now needs to start at 25 (after adding the bit-mapped tidal mascons record )
! PT170609: no, it is either 26 (pre-20170609) or 27 (ie the value of the line BEFORE the first line of data
    irec = 26
    if(rec_length(isat) /= 318788)irec=irec+1   ! add one extra line if it is a post-20170609 GTORB file

    do iepoch = 1, nepochs_t
      irec=irec+1
      read(LUGT(isat),rec=irec,iostat=ioerr)&
         tmp_epoch,(rvec(i,isat,iepoch),i=1,6),(sciframe_quat(i,isat,iepoch),i=1,4) &
! PT140515: read in all the partials here, then extract later the ones requested
!         ,(Xpart(i,isat,iepoch),i=1,6),(Ypart(i,isat,iepoch),i=1,6),(Zpart(i,isat,iepoch),i=1,6) &
!         ,(Xdotpart(i,isat,iepoch),i=1,6),(Ydotpart(i,isat,iepoch),i=1,6),(Zdotpart(i,isat,iepoch),i=1,6)
         ,(all_partials(i,isat,iepoch),i=1,72) ! 6 x 6 partials for pos/vel, then 6 x 3 for bias, then 6 x 3 for scale
         epoch(isat,iepoch)=dble(tmp_epoch)        ! this is to account for GRACEORB writing out as integer but gracefit reading as R*8
    enddo
  enddo
  call status_update('STATUS','GRACEDIFF','gracediff',' ',"Read in the binary GTORB files for pos/vel information",0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! for each epoch, calculate and output the differences in position and velocity
  do iepoch = 1,nepochs_t
    do i=1,6
      dposvel(i) = (rvec(i,2,iepoch) - rvec(i,1,iepoch)) * conv(i)
      if(out_type == "ACR")then
! convert to along/cross/radial
        call xyz2rac(rvec(:,1,iepoch),dposvel(1:3),d_out(1:3))
        call xyz2rac(rvec(:,1,iepoch),dposvel(4:6),d_out(4:6))
      else if (out_type == "NEU")then
!       Convert delta X,Y,Z, to delta N,E,U...
        call rotate_geod( dposvel(1:3),d_out(1:3),'XYZ','NEU',rvec(1:3,1,iepoch)*m_mm, llh,rotmat)
        call rotate_geod( dposvel(4:6),d_out(4:6),'XYZ','NEU',rvec(1:3,1,iepoch)*m_mm, llh,rotmat)
      else
        d_out = dposvel
      endif
    enddo

! PT140207: we should get gracediff to compute the range rate and output it as well
    range_rate = (dposvel(1)*dposvel(4)/conv(4)+dposvel(2)*dposvel(5)/conv(5)+dposvel(3)*dposvel(6)/conv(6)) &
                / dsqrt(dposvel(1)**2+dposvel(2)**2+dposvel(3)**2)
    write(*,100)iepoch,int(epoch(1,iepoch)),dposvel,d_out,out_type,range_rate*1.d6, " dPos/Vel"
100 format(i6,i10,1x,2(3f21.10,3x),6x,2(3f21.10,3x),3x,a3,3x,f14.8, a)
    write(*,110)"Partials: d/d",partial_type,iepoch,(all_partials(part_offset+i,1,iepoch),i=1,6) &
                            ,(all_partials(part_offset+i,2,iepoch),i=1,6)          !,(Xpart(j,1,iepoch),j=1,6),(Ypart(j,1,iepoch),j=1,6),(Zpart(j,1,iepoch),j=1,6)
110 format(a,a4,i7,6e17.9,7x,6e17.9)
  enddo


  end
