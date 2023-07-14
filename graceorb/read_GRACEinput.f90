    subroutine read_GRACEinput(ifile)

! read in the GRACE orbit integrator input file
! contains details of which components to include in calculation
! and model names for static field, ocean and atm tides and
! non-tidal ocean and atm
! modified from Simon's input file read subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! example GRACE.input file:

!SUN_MOON          : Y
!PLANETS           : Y
!SOLID_TIDE        : Y
!POLE_TIDE         : Y
!GENERAL_REL       : Y
!ACCELEROMETER     : Y
!STATIC_FIELD      : Y
!MOD_STATIC_FIELD  : GGM02
!OCEAN_TIDE        : Y
!MOD_OCEAN_TIDE    : FES95
!ATM_TIDE          : N
!MOD_ATM_TIDE      : ATMMD
!DEALIAS           : Y
!MOD_DEALIAS       : AOD04
!MASCON            : Y
!DATA_DIR          : /Users/tony/work/gracesoft/grace_data
!MASCON_DIR        : /Users/tony/work/gracesoft/grace_data/mascon_files/3deg_sphere_1layer
!DEALIAS_DIR       : /Users/tony/work/gracesoft/grace_data/AOD1B
!APRIORI_FILE      : 
!NUM_APRIORI_PARAM : 24
!USE_APRIORI_ICS   : N
!USE_APRIORI_MCS   : N
!MSC_APRIORI_FILE  : 
!
!A1_ACC_BIAS_SCALE : -1.106   27.042   -0.5486    2.233d-4  4.46d-3   -1.139d-6   2.5d-7   1.1d-6  1.7d-7   0.9595   0.9797   0.9485
!A2_ACC_BIAS_SCALE : -1.2095  29.3370  -0.5606   -4.128d-5  6.515d-4  -2.352d-6   9.7d-9  -3.9d-7  3.8d-9   0.9595   0.9797   0.9485
!B1_ACC_BIAS_SCALE : -0.5647  7.5101   -0.8602   -7.788d-5  7.495d-3   1.399d-4   2.4d-7  -9.6d-6  2.5d-7   0.9465   0.9842   0.9303
!B2_ACC_BIAS_SCALE : -0.6049  10.6860  -0.7901   -1.982e-5  1.159d-3   4.783d-5   3.5d-9  -4.3d-7 -6.5d-9   0.9465   0.9842   0.9303
!
!GPS_ANT_A_L1      : 0.4514203544369704 -0.0008860920781892878 -0.0008860920781892878 -0.9999992148405207
!GPS_ANT_A_L2      : 0.4756503363816744 -0.0008409538886124732 -0.0008409538886124732 -0.9999992927963072
!GPS_ANT_B_L1      : 0.4517310303930869 0.001332651421967077 0.001669134837480358 -0.9999977190119395
!GPS_ANT_B_L2      : 0.475960977938318 0.001264809570329978 0.001584163481775421 -0.999997945339296
!
!INT_TIME_INTERVAL : 5
!END_OF_INPUT_FILE :
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use sat_mod
    use inmod_mod
    implicit none

!     variables for status reporting
      integer*4   ierr, jerr,  &   ! IOSTAT errors
                  indx,       & ! Position in string
                   i, j,       &   ! loop counters
                  ifile,      &   ! file number being read
                  len,rcpar,  &   ! command line arguments
                  ihdw          ! header descriptor width
!
      real*8 val
!
      character*100  cdum       ! Dummy character for read_line.
      character*30  cmd         ! Commands read from file header lines.
      character*80  fname,prog_name
! PT121211: message is defined in inred_mod.f90
!      character*256 message, &
       character*256  line,      &  ! Line read from file
                  cval        ! String read from line
!
!      integer*4 lugt(1)


      logical the_end    ! End of header found
      logical found      ! Header separator found

!Set standard gtorb header descriptor width (characters)
      ihdw = 0
      ierr = 0

!     Initialize some variables
!      gt_hdrrecs(ifile) = 0
!      gt_datarecs(ifile) = 0
!      gt_stime(ifile)    = 0.0d0
!      gt_etime(ifile)    = 0.0d0
      the_end = .false.

10    do while (ierr .eq. 0 .and. .not. the_end )

!        Try to read a header record
         read (lugt(ifile),'(a)',iostat=ierr) line
!         if( ierr.eq.0 .and. line(1:1).eq.' ' .and. trim(line).gt.0 ) then
         if( ierr.eq.0 .and.  len(trim(line)).gt.0 ) then
!           Find the index of the first character past the : separator
            indx = 1
            found = .false.
20          do while ( .not. found )
               if ( line(indx:indx+1) .eq. ":" ) then
                  found = .true.
               else
                  indx = indx + 1
               endif
            enddo
            cmd = adjustl(line(1:indx-1))
            call lwr_to_upper(cmd)

!           Move the index one space past the : separator
            ihdw = indx + 2

!           See which header line is found
            if (cmd(1:8) .eq. 'SUN_MOON'                 ) then
!           Get sun moon flag (Y/N)
               indx = ihdw
               call first_word(line, cdum, indx)
               gt_sunmoon(ifile) = trim(cdum)

            elseif (cmd(1:7) .eq. 'PLANETS'                       ) then
!            Get planet point mass flag (Y/N)
               indx = ihdw
               call first_word(line, cdum, indx)
               gt_planet(ifile) = trim(cdum)

            elseif (cmd(1:10) .eq. 'SOLID_TIDE'                       ) then
!            Get solid earth tide flag  (Y/N)
               indx = ihdw
               call first_word(line, cdum, indx)
               gt_solidtide(ifile) = trim(cdum)

            elseif (cmd(1:9) .eq. 'POLE_TIDE'                       ) then
!            Get pole tide flag  (Y/N)
               indx = ihdw
               call first_word(line, cdum, indx)
               gt_poletide(ifile) = trim(cdum)

            elseif (cmd(1:11) .eq. 'GENERAL_REL'                       ) then
!            Get general relativity flag (Y/N) 
               indx = ihdw
               call first_word(line, cdum, indx)
               gt_genrel(ifile) = trim(cdum)

            elseif (cmd(1:13) .eq. 'ACCELEROMETER'                       ) then
!            Get accelerometer flag  (Y/N)
               indx = ihdw
               call first_word(line, cdum, indx)
               gt_acc(ifile) = trim(cdum)

            elseif (cmd(1:12) .eq. 'STATIC_FIELD'                       ) then
!            Get static field flag  (Y/N)
               indx = ihdw
               call first_word(line, cdum, indx)
               gt_statfield(ifile) = trim(cdum)

            elseif (cmd(1:16) .eq. 'MOD_STATIC_FIELD'                       ) then
!            Get static field model (XXXXX)
               indx = ihdw
               call first_word(line, cdum, indx)
               gt_statfieldmod(ifile) = trim(cdum)

            elseif (cmd(1:10) .eq. 'OCEAN_TIDE'                       ) then
!            Get ocean tide flag  (Y/N) ! PT140204: also allow 'G' for ocean tides from a grid rather than from spherical harmonics
               indx = ihdw
               call first_word(line, cdum, indx)
               gt_oceantide(ifile) = trim(cdum)

            elseif (cmd(1:14) .eq. 'MOD_OCEAN_TIDE'                       ) then
!            Get ocean tide model (XXXXX)
               indx = ihdw
               call first_word(line, cdum, indx)
               gt_oceantidemod(ifile) = trim(cdum)

            elseif (cmd(1:8) .eq. 'ATM_TIDE'                       ) then
!            Get atm tide flag  (Y/N)
               indx = ihdw
               call first_word(line, cdum, indx)
               gt_atmtide(ifile) = trim(cdum)

            elseif (cmd(1:12) .eq. 'MOD_ATM_TIDE'                       ) then
!            Get atm tide model (XXXXX)
               indx = ihdw
               call first_word(line, cdum, indx)
               gt_atmtidemod(ifile) = trim(cdum)

! PT121101: distinguish between use dealias and the dealias directory
            elseif (cmd(1:9) .eq. 'DEALIAS '                       ) then
!            Get dealias flag  (Y/N)
               indx = ihdw
               call first_word(line, cdum, indx)
               gt_dealias(ifile) = trim(cdum)

            elseif (cmd(1:12) .eq. 'MOD_DEALIAS'                       ) then
!            Get dealias model (XXXXX)
               indx = ihdw
               call first_word(line, cdum, indx)
               gt_dealiasmod(ifile) = trim(cdum)

! PT121101: distinghush between use mascons and the mascon directory
            elseif (cmd(1:7) .eq. 'MASCON '                       ) then
!            Get mascon flag 
               indx = ihdw
               call first_word(line, cdum, indx)
               gt_mcon(ifile) = trim(cdum)
!               gt_mascon(ifile) = line(ihdw:trim(line))

! PT130906: get the name of the mascon directory
            elseif (cmd(1:10) .eq. 'MASCON_DIR'   ) then
              continue

! APP130227: Create option to read ICs from a Gracefit style output file
            elseif (cmd(1:12) .eq. 'APRIORI_FILE'                ) then
!            Get apriori file (in gracefit format)
               indx = ihdw
               call first_word(line, cdum, indx)
               apriori_file = trim(cdum)

! PT190917: Create option to read ICs from a Gracefit style output file
            elseif (cmd(1:16) .eq. 'MSC_APRIORI_FILE'            ) then
!            Get apriori file (in gracefit format)
               indx = ihdw
               call first_word(line, cdum, indx)
               msc_apriori_file = trim(cdum)

! APP130317: Create option to read orbit IC values from a Gracefit style output file
            elseif (cmd(1:15) .eq. 'USE_APRIORI_ICS'                ) then
!            Get apriori file (in gracefit format)
               indx = ihdw
               call first_word(line, cdum, indx)
! PT130905: make this an integer, bit-mapped to indicate which ICs to update a priori
! PT160215: increase to a 3-digit integer
               read(cdum,'(i3)')gt_use_apriori_ics(1)

! APP130321: Create option to read orbit mascon values from a Gracefit style output file
            elseif (cmd(1:15) .eq. 'USE_APRIORI_MCS'                ) then
!            Get apriori file (in gracefit format)
               indx = ihdw
               call first_word(line, cdum, indx)
               gt_use_apriori_mcs(1) = trim(cdum)

            elseif (cmd(1:17) .eq. 'NUM_APRIORI_PARAM'                       ) then
!  Get number of parameters in the apriori file
               indx = ihdw
               call first_word(line, cdum, indx)
               gt_num_apriori_par(ifile) = trim(cdum)

            elseif (cmd(1:17) .eq. 'A1_ACC_BIAS_SCALE'   ) then
              continue
            elseif (cmd(1:17) .eq. 'A2_ACC_BIAS_SCALE'   ) then
              continue
            elseif (cmd(1:17) .eq. 'B1_ACC_BIAS_SCALE'   ) then
              continue
            elseif (cmd(1:17) .eq. 'B2_ACC_BIAS_SCALE'   ) then
              continue
            elseif (cmd(1:12) .eq. 'GPS_ANT_A_L1'   ) then
              continue
            elseif (cmd(1:12) .eq. 'GPS_ANT_A_L2'   ) then
              continue
            elseif (cmd(1:12) .eq. 'GPS_ANT_B_L1'   ) then
              continue
            elseif (cmd(1:12) .eq. 'GPS_ANT_B_L2'   ) then
              continue
            elseif (cmd(1:12) .eq. 'GPS_ANT_C_L1'   ) then
              continue
            elseif (cmd(1:12) .eq. 'GPS_ANT_C_L2'   ) then
              continue
            elseif (cmd(1:12) .eq. 'GPS_ANT_D_L1'   ) then
              continue
            elseif (cmd(1:12) .eq. 'GPS_ANT_D_L2'   ) then
              continue
! PT121101: skip over the data directory entries
            elseif (cmd(1:8) .eq. 'DATA_DIR'   ) then
              continue
            elseif (cmd(1:11) .eq. 'DEALIAS_DIR'  ) then
              continue
            elseif (cmd(1:17) .eq. 'INT_TIME_INTERVAL' ) then
!  Get integration step size (sec)
               indx = ihdw
               call first_word(line, cdum, indx)
               gt_intstpsz(ifile) = trim(cdum)

            elseif (cmd(1:17) .eq. 'END_OF_INPUT_FILE' ) then
                  the_end = .true.

            else
               write(*,800) trim(line)
800            format('** WARNING ** Unknown header line: ',a)
            endif
         endif
      enddo

! PT210820: put a catch here for when "MASCONS   : Y/O" but "USE_APRIORI_MCS : N". This is a non-sensical selection.
      if( (gt_mcon(1) == "O" .or. gt_mcon(1) == "Y") .and. gt_use_apriori_mcs(1) == "N")then
        write(message,'(a,a1,a,a1,a)')"Cannot model mascon accelerations (MASCON   : ",gt_mcon(1) &
                                   ,") without a priori values (USE_APRIORI_MCS  : ",gt_use_apriori_mcs(1) &
                                   ,"). They will be set to zero."
        call status_update('STATUS','GRACEORB','read_GRACEinput',' ',message,0)
      else if (gt_mcon(1) == "N" .and. gt_use_apriori_mcs(1) == "Y")then
        write(message,'(a)')"Requested to NOT use mascons but to use apriori mascon values. The latter is ignored!"
        call status_update('STATUS','GRACEORB','read_GRACEinput',' ',message,0)
        gt_use_apriori_mcs(1) = "N"
      endif
      return 
      end subroutine read_GRACEinput













