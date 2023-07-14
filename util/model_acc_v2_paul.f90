 program model_acc_v2

! program to compute a model  to flatten out the non-linear shapes that we see in the
! uncalibrated accelerometer data, especially in the cross-track
!
! P. Tregoning
! 31 May 2017
!
! PT180702: modified to fit the linearising model to each component separately.
! PT180803: added EMD1 filter as an option to remove long-wavelength components
! PT180806: output original observations as well as model-removed observations
! RM180913: added EMD2 filter

  implicit none

! command line arguments
  character*150 :: acc_file         ! input ACC1B file (format can be GRACE, GRACE FO or ???? )
  character*100 :: arg
  character*20  :: model_type(3)       ! character string to define the model to fit ("quadratic", "exponential", [add more] )
  integer*4     :: bitmap           ! bit-map variable to indicate which components of accelerometer data should be modelled
  integer*4     :: start_ep(3),end_ep(3)  ! start and end epochs over which to fit the model

! thrust file (made from ACC1B file name)
  character*150 :: thr_file   

! unit numbers
  integer*4     :: luacc_in,luthr_in

! space gravity mission identifier (0: GRACE, 1: GRACE FO, 2: GRACE II)
  integer*4     :: mission

! accelerometer arrays and start/stop times
  integer*4     :: acc_span(2)                     ! start/stop time (in grace seconds) of accelerometer data
  real(kind=8),allocatable :: acc_obs(:,:,:)       ! array to store the ACC1B accelerometer observations
  real(kind=8),allocatable :: acc_obs_fixed(:,:,:) ! array to store the fixed/filtered/detrended accelerometer observations
  integer*4     :: nvals_acc
  integer*4     :: start_flag,end_flag

! thrust variables
  integer*4     :: n_thrusts                  ! number of thrusts found in THR1B file
  integer*4   ,allocatable :: thr_obs(:,:,:)  ! array to store the thrust observations
  integer*4     :: thrust_width               ! number of points to flag either side of thrust time
  logical       :: need_thrusts               ! .true. for models where we need to ignore thrust-affected observations

! variables related to setting up different estimation strategies
  integer*4     :: nparams                    ! number of parameters in the model
  real(kind=8)  :: params(10)                 ! parameter values estimated when fitting model to data

! variables related to the EMD filter
  integer*4     :: ndecomp                    ! max number of decomposition components for the EMD filter  
  integer*4     :: start_comp,end_comp        ! start/end components required for reconstructed acc obs
  real(kind=8),allocatable :: acc_decomp(:,:) ! the EMD decomposition of the signal
  real(kind=8),allocatable :: valrange(:)     ! peak-to-peak range of each EMD decomposition component

! variables for the PSD computations
  integer*4     :: start_epoch,end_epoch

! other variables
  integer*4     :: ioerr,next_row,icomp,iepoch,j,narg
  character     :: calling_prog*17
  logical       :: bitmap
  character*250 :: message

! debug
  logical debug
  integer*4 :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  luacc_in = 10
  luthr_in = 11
  calling_prog = "UTIL/MODEL_ACC_V2"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! decode command line
  call getarg(1,acc_file)
  if(acc_file(1:1) == "")then
    print*,"model_acc_v2 acccelerometer_file  model_X start end model_Y start end model_Z start end"
    print*,"e.g. model_acc_v2 ACC1B_2016-07-28_A_02.asc none quadratic start_fit end_fit quadratic start end"
    print*," "
    stop
  endif

! character string for which model to use
  narg = 2
  do i=1,3
    call getarg(narg,model_type(i))
    call getarg(narg+1,arg)
    read(arg,*)start_ep(i)
    call getarg(narg+2,arg)
    read(arg,*)end_ep(i)
    narg = narg + 3
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now, work out whether we need thrust data or not
  need_thrusts = .false.
  do icomp = 1,3
    if(model_type(icomp)(1:9) == "quadratic" .or. model_type(icomp)(1:11) == "exponential")then
      need_thrusts = .true.
    endif
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! open files
  call level1B_open(luacc_in,calling_prog,acc_file)
  if(need_thrusts)then
    thr_file = acc_file
    write(thr_file(1:3),'(a3)')"THR"
    call level1B_open(luthr_in,calling_prog,thr_file)
  else
    call status_update('STATUS',calling_prog,' ',' '," No thrust data required for requested ACC1B modelling",0)
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  A C C 1 B   D A T A
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the ACC1B header
  mission = -99
  call acc_read_hdr(luacc_in,calling_prog,acc_file,mission,acc_span)

! now we can dimension the arrays to hold the accelerometer observations.
  nvals_acc = acc_span(2) - acc_span(1) + 1
  allocate(acc_obs(1,nvals_acc,5))              ! gracesec, XYZ accelerations, 5th component is a thrust flag
  allocate(acc_obs_fixed(1,nvals_acc,5))              ! gracesec, XYZ accelerations, 5th component is a thrust flag

! set all the acc obs flags to good a priori
  acc_obs(1,:,5) = 1

! read the ACC1B linear accelerations
  call acc_read_data(luacc_in,calling_prog,acc_file,mission,nvals_acc,3,acc_span,acc_obs)
!do iepoch=1,nvals_acc
!  print*,iepoch,acc_obs(1,iepoch,2)
!enddo
!stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  T H R 1 B   D A T A
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PT180806: we don't need thrust data for the EMD ..
  if(need_thrusts)then
    ! read the THR1B header
    call thr_read_hdr(luthr_in,calling_prog,thr_file,mission,n_thrusts)

    ! now we can dimension the arrays to hold the thrust observations.
    allocate(thr_obs(1,n_thrusts,8))              ! gracesec of times of thrusts (we don't care which thruster fired, just that something fired) 

    ! read the THR1B linear thrust data
    call thr_read_data(luthr_in,calling_prog,thr_file,mission,n_thrusts,n_thrusts,thr_obs)
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! flag accelerometer observations within 120 seconds of a thrust to eliminate obs that
! may be affected by the CRN filtering of the square thrusts
  do i = 1,n_thrusts
    ! remove points "N" seconds before the thrust
    thrust_width = 60
    start_flag = (thr_obs(1,i,1)-acc_span(1))-thrust_width
    if( (thr_obs(1,i,1)-thrust_width) < acc_span(1) )start_flag = 1    ! check that we don't start before the data

    ! remove points 60 seconds after the thrust
    end_flag = (thr_obs(1,i,1)-acc_span(1)) + thrust_width
    if( (thr_obs(1,i,1)+thrust_width) > acc_span(2) ) end_flag = (acc_span(2)-acc_span(1)+1)    ! check that we don't end after the data
    ! set the flags to false
    acc_obs(1,start_flag:end_flag,5) = 0
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   R E M O V E    S O M E    S O R T    O F    M O D E L
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  here we will fit a model to one or more components of the accelerometer observations
!  and remove it from the data. This is done to (try to) create a set of observations
!  that have a quasi-repetitive once/rev pattern but have zero slope across the time series.

! Option 1: fit a quadratic and remove the rate and acceleration terms. This is what the
!           old ga/util/model_acc used to do

! PT180702: we now allow a different model to be applied to each component, with different time spans as well
  do icomp = 1,3
! *** do nothing ***
    if(model_type(icomp)(1:4) == "none")then
      write(message,'(a,i2,a)')" No model applied to component",icomp,"."
      call status_update('STATUS','UTIL',calling_prog,' ',message,0)
      acc_obs_fixed(1,:,icomp+1) = acc_obs(1,:,icomp+1)

! *** remove a quadratic ***
    else if (model_type(icomp)(1:9) == "quadratic")then
      call acc_fit_quadratic(.true.,calling_prog,nvals_acc,acc_obs(1,:,icomp+1),acc_obs_fixed(1,:,5) &
             ,start_ep(icomp),end_ep(icomp),params(1:3))
      write(message,'(a,3e20.8,a)')"  Model parameter values (accel, rate, offset) (nm/s^2):",params(1:3)*1.e9
      call status_update('STATUS','UTIL','util/model_acc_v2',' ',message,0)

! *** remove an exponential (doesn't work!!!) ***
    else if (model_type(icomp)(1:11) == "exponential")then
      call acc_fit_exponential(.false.,calling_prog,nvals_acc,acc_obs(1,:,icomp+1),acc_obs_fixed(1,:,5) &
             ,start_ep(icomp),end_ep(icomp),params(1:5))
      write(message,'(a,5e20.8)')"  Exponential model parameter (nm/s^2):",params(1:5)*1.e9
      call status_update('STATUS','UTIL','util/model_acc_v2',' ',message,0)

! *** EMD1 filter ***
    else if (model_type(icomp)(1:4) == "EMD1")then
      ! set the EMD decomposition variables
      ndecomp    = 13
      start_comp =  1
      debug = .false.
      allocate(acc_decomp(ndecomp,nvals_acc))
      allocate(valrange(ndecomp))
      call acc_EMD1(debug,calling_prog,icomp,nvals_acc,ndecomp,nint(acc_obs(1,nvals_acc,icomp+1)-acc_obs(1,1,icomp+1)) &
                 ,acc_obs(1,:,icomp+1),acc_decomp,end_comp)

!! PT180820: try removing the identified long wavelength features, then running the EMD again
!      do iepoch = 1,nvals_acc
!         acc_obs_fixed(1,iepoch,icomp+1) = sum(acc_decomp(start_comp:end_comp+1,iepoch))
!      enddo
!      call acc_EMD1(debug,calling_prog,icomp,nvals_acc,ndecomp,nint(acc_obs(1,nvals_acc,1)-acc_obs(1,1,1)) &
!                 ,acc_obs_fixed(1,:,icomp+1),acc_decomp,end_comp)

      ! reconstruct with only the components that don't contain the non-linear signal
      do iepoch = 1,nvals_acc
         acc_obs_fixed(1,iepoch,icomp+1) = sum(acc_decomp(start_comp:end_comp-1,iepoch))
      enddo
      write(message,'(a,i1,a,i3,a,i3,a)')"  Axis ",icomp,": removed long-wavelength components",end_comp," to",ndecomp &
                                          ," of EMD filter."
      call status_update('STATUS','UTIL','util/model_acc_v2',' ',message,0)
      deallocate(acc_decomp)
      deallocate(valrange)

! *** EMD2 filter ***
    else if (model_type(icomp)(1:4) == "EMD2")then
      ! set the EMD decomposition variables
      ndecomp    = 13
      start_comp =  1
      debug = .false.
      allocate(acc_decomp(ndecomp,nvals_acc))
      allocate(valrange(ndecomp))
      call acc_EMD2(debug,calling_prog,icomp,nvals_acc,ndecomp,nint(acc_obs(1,nvals_acc,icomp+1)-acc_obs(1,1,icomp+1)) &
                  ,acc_obs(1,:,icomp+1),acc_decomp,end_comp)
      ! reconstruct with only the components that don't contain the non-linear signal
      do iepoch = 1,nvals_acc
         acc_obs_fixed(1,iepoch,icomp+1) = sum(acc_decomp(start_comp:end_comp-1,iepoch))
      enddo
      write(message,'(a,i1,a,i3,a,i3,a)')"  Axis ",icomp,": removed long-wavelength components",end_comp," to",ndecomp &
                                          ," of EMD2 filter."
      call status_update('STATUS','UTIL','util/model_acc_v2',' ',message,0)
      deallocate(acc_decomp)
      deallocate(valrange)

! *** EMD3 filter ***
    else if (model_type(icomp)(1:4) == "EMD3")then
      ! set the EMD decomposition variables
      ndecomp    = 13
      start_comp =  1
      debug = .false.
      allocate(acc_decomp(ndecomp,nvals_acc))
      allocate(valrange(ndecomp))
      call acc_EMD2(debug,calling_prog,icomp,nvals_acc,ndecomp,nint(acc_obs(1,nvals_acc,icomp+1)-acc_obs(1,1,icomp+1)) &
                  ,acc_obs(1,:,icomp+1),acc_decomp,end_comp)

! PT180926: invoke a power-spectrum density test to check whether the excluded components do actually contain any signal
!           at 1/rev or higher. If so, move the end_comp value so that they are included.
      if(nvals_acc > 86400)then
        start_epoch = nint(dble(nvals_acc - 86400)/2.d0) + 1
        end_epoch = start_epoch + 86400
      else
        start_epoch = 1
        end_epoch = nvals_acc
      endif
end_comp = 0
      call acc_emd_psd(debug,calling_prog,nvals_acc,ndecomp,start_epoch,end_epoch,acc_decomp,end_comp)

      ! reconstruct with only the components that don't contain the non-linear signal
      do iepoch = 1,nvals_acc
         acc_obs_fixed(1,iepoch,icomp+1) = sum(acc_decomp(start_comp:end_comp-1,iepoch))
      enddo
      write(message,'(a,i1,a,i3,a,i3,a)')"  Axis ",icomp,": removed long-wavelength components",end_comp," to",ndecomp &
                                          ," of EMD2 filter."
      call status_update('STATUS','UTIL','util/model_acc_v2',' ',message,0)
      deallocate(acc_decomp)
      deallocate(valrange)


! *** model not coded. ***
    else
      write(message,'(a,a,a)')'Model "',model_type,'" not coded. You will have to write it yourself ...'
      call status_update('FATAL','UTIL','model_acc_v2',' ',message,0)
    endif

  enddo

! PT180702: we need this program to output something ..... how about in ACC1B format without the header?
  print*,"GRACE SECONDS        Acc_X              Acc_Y                Acc_Z     (um/s^2)      Models applied: ",model_type
  do iepoch = 1,nvals_acc
    write(*,'(i12,6e25.15)')nint(acc_obs(1,iepoch,1)),(acc_obs_fixed(1,iepoch,j),j=2,4),(acc_obs(1,iepoch,j),j=2,4)
  enddo

  end




  
      SUBROUTINE spctrm(p,m,k,ovrlap,w1,w2)
      INTEGER k,m
      REAL p(m),w1(4*m),w2(m)
      LOGICAL ovrlap
!CU    USES four1
      INTEGER j,j2,joff,joffn,kk,m4,m43,m44,mm
      REAL den,facm,facp,sumw,w,window
      window(j)=(1.-abs(((j-1)-facm)*facp))
!C     window(j)=1.
!C     window(j)=(1.-(((j-1)-facm)*facp)**2)
      mm=m+m
      m4=mm+mm
      m44=m4+4
      m43=m4+3
      den=0.
      facm=m
      facp=1./m
      sumw=0.

!print*,'in spctrm m,k:',m,k

      do 11 j=1,mm
        sumw=sumw+window(j)**2
11    continue
      do 12 j=1,m
        p(j)=0.
12    continue
      if(ovrlap)then
!print*,'reading from file 9'
        do j=1,m
          read (9,*) w2(j)
!print*,'j,w2(j)',j,w2(j)
        enddo
      endif
      do 18 kk=1,k
        do 15 joff=-1,0,1
          if (ovrlap) then
            do 13 j=1,m
              w1(joff+j+j)=w2(j)
!print*,'j,w2(j)',j,w2(j)
13          continue
!print*,'line366 reading from 9. kk,joff',kk,joff
            read (9,*) (w2(j),j=1,m)
            joffn=joff+mm
            do 14 j=1,m
              w1(joffn+j+j)=w2(j)
14          continue
          else
            read (9,*) (w1(j),j=joff+2,m4,2)
          endif
15      continue
        do 16 j=1,mm
          j2=j+j
          w=window(j)
          w1(j2)=w1(j2)*w
          w1(j2-1)=w1(j2-1)*w
16      continue
        call four1(w1,mm,1)
        p(1)=p(1)+w1(1)**2+w1(2)**2
        do 17 j=2,m
          j2=j+j
          p(j)=p(j)+w1(j2)**2+w1(j2-1)**2+w1(m44-j2)**2+w1(m43-j2)**2
17      continue
        den=den+sumw
18    continue
      den=m4*den
      do 19 j=1,m
        p(j)=p(j)/den
19    continue
      return
      END
!C  (C) Copr. 1986-92 Numerical Recipes Software .
!C------------------------------------------------------------------------------

!C------------------------------------------------------------------------------
      SUBROUTINE four1(data,nn,isign)
      INTEGER isign,nn
      REAL data(2*nn)
      INTEGER i,istep,j,m,mmax,n
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
      END
!C  (C) Copr. 1986-92 Numerical Recipes Software .
!C------------------------------------------------------------------------------

