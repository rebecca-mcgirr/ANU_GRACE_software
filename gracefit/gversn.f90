SUBROUTINE GVERSN(VERSION)

  !   Update history for GRACE utilities GRACEFIT ....

  CHARACTER(10) :: GETMAC,MACHIN
  CHARACTER(35) :: VERS
  character(45) :: libver
  CHARACTER(120) :: VERSION

  integer*4 :: last_nonblank

  !MACHIN = GETMAC(1)

  !   Get the orbit utilities version
  WRITE (VERS,5) MACHIN(1:last_nonblank(machin))
5 format ('2.00 2020/10/14 11:11:00 (',a,')')

  !   Get library version
  !CALL LVERSN(libver)

  version = ' '
  WRITE(version,10) vers,libver
10 FORMAT('ver. ',a35,' Library ver. ',a45)

!   0.10     GRACEFIT converted to F90 McClusky 120515
!            PARTL, NORMINC:  fixed bug in storing and using range-rate partials   Tregoning 121219
!            GRACEFIT, PARTL, NORMINC, READ_INPUT: added conditional equation to make calibrated accelerometer values the same Tregoning/McClusky 130524
!            NORMINC_COND : new routine to add conditional equation above  Tregoning 130524
!            GRACEFIT, READ_INPUT, KB_PART, PARTL : add estimation of GPS antenna offset values McClusky/Tregoning 130528 
!            GRACEFIT: fixed bug in computation of velocity adjustments as a result of CoM estimates Purcell/Tregoning 130617
!            GRACEFIT, SHADOW_GRACE: compute whether the satellites are in the Earth's shadow  Tregoning 130620
!            OPEN_ORBFILES: stop if ut1./pole./nutabl. not found  Tregoning 130620
!            CALC_BETA_ANGLE: compute the beta angle (elevation of sun wrt orbital plane)  Tregoning 130621
!
!   2.00     Completely rewritten software, changing names of nearly all subroutines and variables   Greenspan/Tregoning/McClusky 130801
!            GRACEFIT, CALC_BETA_ANGLE: compute mean beta angle rather than just one value   Tregoning 130814
!            LS_norminc: only compute (then mirror) upper triangle of normal equations (40% faster)  Tregoning 130814
!            output_writeVCV: re-inserted "V2" as the first two characters of the first line    Tregoning 130819
!            kb_lib/kb_sponge: added offset and rate to the 1/rev model in kband  Tregoning 140108
!            GRACESIM: added output information explaining which GTORB files are "truth" and which are "modelled"   Tregoning 140113
!            GRACEADD : new program to add adjustments to pos/vel values in a GTORB file and output new files  Tregoning 140207
!            GRACEFIT, WRITEJPL : fix bugs for case when no kbrr obs are used  Tregoning 140501
!            ADDNORM : added constraints to the scale parameters  Tregoning 140506
!            GRACEFIT : replace filtered omc for kbrr with zero value if there is no kbrr obs for the epoch Tregoning 140527
!            kb_kbrrDECOMPOSE : compress the pre_omc vector to remove data gaps before filtering, then re-expand afterwards  Tregoning 140529
!            GRACEFIT/GRACESIM : pass kRANGE obs array into kb_kbrrDECOMPOSE so that data gaps can be identified  Tregoning 140529
!            output_writeNORM  : add ocean/land mascon flags to record 4 of the binary normal equations file     Tregoning 140610
!            ADDNORM           : read ocean/land mascon flags from record 4 of the binary normal equations file  Tregoning 140610
!
!   3.00     Add ocean tide constituents as parameters for ocean mascons
!            HEADER_readBIN : added read of bit-mapped tidal mascon info from record after "END OF HEADER"  Tregoning 140613
!            GRACEFIT.h     : added variable to store bit-mapped tidal mascon info   Tregoning 140613
!            input_openGTORB : increase record length by 80000 to allow tidal amplitude partials  Tregoning 140613
!            GRACEFIT.h     : added pointer to msc tide parameters  Tregoning 140616
!            HEADER_readBIN, command_read_SETUP : added code to handle properly the estimating (or not) of tidal mascon amplitudes  Tregoning 140616
!            GRACEFIT       : output prefit obs RMS to status file   Tregoning 140618
!            command_readAprConst : added mascon tidal amplitude a priori constraint   Tregoning 140618
!            GRACEFIT, plot_writePLT, plot_writeKB, output_writeFIT : limit statistics to last used epoch (rather than the whole orbit length)  Tregoning 140618
!            GRACEDIFF : update record length and number of header records to read GTORB files with mascon tides  Tregoning 140618 
!            GRACESIM  : updated to be compatible with the recent changes made to GRACEFIT   Tregoning 140626
!            input_readGTORB : fixed bug in assigning apr_prm names when no tidal mascons estimated   Tregoning 140626
!            GRACESIM : fixed bug concerning the number of parameters when tidal mascon amplitudes are estimated  Tregoning 1406707
!            output_writeNORM : add number of mascon tidal amplitudes to header of normal equations file   Tregoning 140709
!            GRACESIM, GRACEFIT, CALC_BETA_ANGLE, shadow_satelliteStatus: added iepoch to the calling arguments   Tregoning 140801
!            RPY_AOC : fixed bug in the computation of GRACE B AOC range correction (incorrect division by LOS_mag)  Tregoning 140807
!
!   4.00     Add twice-per-rev along-track accelerations to the GTORB files
!            input_openGTORB : increase record length of GTORB file to include an additional 12 partials                         Tregoning 140821
!            command_readSetup, gracefit.h90 : add two-per-rev parameters to command file and include file. Also added
!            paremeter names to prmnam for the twice-per-rev variables                                Tregoning 140821
!            output_codes2labels : added translation of twice-per-rev labels                          Tregoning 140821
!            input_readGTORB : added transfer of partials for twice-per-revs, and increased the dimension of tmprvec accordingly Tregoning 140821
!            command_readAprConst : added code for apriori twice-per-revs                             Tregoning 140821
!            command_readSetup, GRACEFIT : read from input command file the start/end epochs to be stacked in normal equations   Tregoning 140901 
!            command_storeCommand  : fix bug when the command name is present but no value is listed after it                    Tregoning 140901
!            output_writeFIT, GRACEFIT   : write out statistics for only epochs stacked in the normal equations                  Tregoning 140901
!            GRACEFIT : initialised GPSAnt_adj to be zero                                             Tregoning 140918
!            ADD_NEAMACC_COND : fixed bug in sign of partials for GRACE B XY mean condition equations Tregoning/Purcell 140922
!            SHADOW_APPLYCONDITION : calibrate the acc obs in the scale partials by adding the bias   Tregoning 140922
!            ADDNORM  : read the mascon constraint matrix filename and mascon constraint scale factor from the input command file. Tregoning 140923
!            ADDNORM  : read (and print) a header line from the mascon constraint file  Tregoning 140923
!            OUTPUT_WRITEVCV : added number of mascons and tidal constituent amplitudes to header line of # parameters           Tregoning 141120
!            OUTPUT_WRITEVCV : write out zero tidal constituent amplitudes if mascon tides not estimated   Tregoning 150509
!            shadow_isINSHADOW : undid the change of sign of the satpos vector (it is already sat wrt earth)  Tregoning 150522
!            shadow_satelliteStatut: undid the change of sign of the sun vector wrt Earth. Also changed the order
!                   of the subtraction to calculate the sun_to_sat vector   Tregoning 150525
!            shadow_getSun   : new subroutine to return Sun coords wrt Earth in an earth-fixed frame  Tregoning 150525
!            shadow_SRP  : new routine to compute a scalar relation between sclY and scales of the other axes                    Tregoning 150525
!            shadow_applyCondition : changed the addition of the bias to the acc obs when computing the O-C                      Tregoning 150603
!            command_readSetup : increased ncond_t from 2 to 4                                        Tregoning 150609
!            GRACEFIT, shadow_ApplyCondition : pass iepoch and full arrays of accel obs and thrusters_off into subroutine        Tregoning 150609
!            shadow_find_AccObs  : new routine to find an epoch half-a-revolution away from a shadow obs but not affected by thrusts  Tregoning 150609
!            shadow_ApplyCondition : removed - once and for all - the addition of the bias to the acc obs. It doesn't work if we add the 
!                                    bias before scaling the obs, so DO NOT do it.                    Tregoning 150612
!            shadow_ApplyCondition : added more condition obs, being 1/2 revolution away from shadow periods.                    Tregoning 150612
!            GRACEFIT, command_READ_SETUP, gracefit.h90    : add mascon mass-conserving conditional equation                     Tregoning 150721
!            ADDNORM  : added mass conservation constraint plus an error message if spatial constraint file not found            Tregoning 150722
!            GRACEFIT : added code to handle header lines in the mascon spatial constraint file       Tregoning 150729
!            PLOT_LIB/plt_writeKB : output postfit pos/vel residuals to kb file                       Tregoning 150810
!            GRACESIM : transferred several changes from GRACEFIT code to GRACESIM                    Tregoning 150812
!            input_READGTORB : fixed bug in pointers for tidal mascon partials                        Tregoning 150825
!            GRACEFIT, GRACESIM, command_readCommand : added apr_tide_ampl_const variable for constraints on tidal amplitudes
!            SHADOW_applyCondition : only apply shadow scale factor constraint obs if requested in command file                  Tregoning 151015 
!            SHADOW_getSun,SHADOW_SRP, SHADO_SatelliteStatus : added starting_epoch as an argument    Tregoning 151022
!            GRACEFIT : set end_neq to be no later than the end of the orbit integration              Tregoning 160223
!            GRACEFIT, SHADOW_applyCondition_drag : change AccZ shadow condition equation             Tregoning 160405
!            SHADOW_applyCondition_drag : add new acc_Z condition that makes trailing satellite calibrated values > leading sat  Tregoning 160421
!            GRACEFIT : fixed bug in applying a priori constraints on tidal amplitude parameters      Tregoning 160428
!            command_readAprConst : fixed formatting of output print statements on tidal amplitude constraints                   Tregoning 160428
!            GRACESIM  :  change epoch loop to stop at end_neq rather than last_epoch                 Tregoning 160622
!            HEADER_readbinGTORB  : fixed bugs in generating parameter labels for tidal amplitudes    Tregoning/Swadba 160622
!            ADDNORM   :  updated to Kellie's version of addnorm.f90, which reads a command file for constraint values
!                         and also fixes a couple of bugs related to tidal amplitude parameters  Tregoning 161026
!            GRACEFIT/GRACESIM    : read more header information from the mascon regularization file  Tregoning 161128
!            ADDNORM   :  updated to read mutliple header lines in the regularisation file            Tregoning 161130
!            ADDNORM   :  fix output mascon statement so that addnorm.out aligns with gracefit/gracesim fit files                Tregoning 161130
!            HEADER_LIB:  fixed bug in reading mascon line in GTORB header                            Tregoning 170228
!            GRACEFIT  :  removed duplicate read of 1st line of regularization file                   Tregoning 170307
!            INPUT_READGTORB : transfer roll/pitch/yaw partials to rvec. Update counters for subsequent parameters               Tregoning 170410
!            OUTPUT_CODES2LABELS : added the rpy parameters                                           Tregoning 170420
!
!  5.00  allow solutions using accelerations
!            COMMAND_READSETUP     : change fatal to warning when range accelerations selected        Tregoning 170427
!            KB_COMPUTETHEORETICAL : changed equation for theoretical range acceleration              Tregoning 170427
!            KB_KBRRDECOMPOSE      : modified so that it can be called for range, rate or accel       Tregoning 170427
!            GRACEFIT              : modified to send only range, rate or accel obs to KBRRDecompose  Tregoning 170427
!            GRACEFIT              : added subroutine to compute the KBRA                             Tregoning 170428
!            OUTPUT_writeFIT       : fixed bug in dimensioning of apr_wght(nobs_t)                    Tregoning 170516
!            OUTPUT_writeFIT       : change formatting of output of obs uncertainties in fit file     Tregoning 170516
!            ADDNORM               : reordered to put input filenames at the end of output line       Tregoning 170601
!            GRACESIM              : added subroutine to compute the KBRA                             Tregoning 170428
!            GRACEDIFF             : redimensioned record length of GTORB to include more mascons     Tregoning 170609
!
!  6.00  record length of GTORB file defined within (or from) the file itself
!            INPUT_openGTORB       : decide upon GTORB record length from the file itself             Tregoning 170609
!            HEADER_readbinGTORB   : user file record length to know how to read the file             Tregoning 170609
!            HEADER_readbinGTORB   : read total_ocean_prim from GTORB header info                     Tregoning 170611
!            GRACEFIT,GRACESIM     : pass in the rec_len of the binary file and pass back from 
!                                     readbinGTORB the number of total_ocean_prim                     Tregoning 170617
!            plot_writeKB          : fixed output format for rpy                                      Tregoning 170622
!            GRACEFIT              : fix (long-standing) problem of filtering kbr when missing data   Tregoning 170703
!            GRACEFIT              : fix bugs in kbra residuals when no kbrr filtering                Tregoning 170703
!            GRACEFIT              : changed OMP handling of normeq, AtWb to increase speed by 22%    Allgeyer  170907
!            LS_norminc            : introduced tempsum variable to increase speed by 11%             Allgeyer  170907
!            KB_filtered, command_readSetup : added number of first component for KBRR/KBRA reconstruction  Tregoning 170911
!
!  7.00  added capability to read GTORB in HDF5 format
!            GRACEFIT              : detects and reads GTORB in hdf5 format                            Allgeyer 171009
!            KBLIB                 : changed handling of edge effects by EMD and NR filters            Allgeyer 180606
!            GRACEFIT,COMMAND_LIB  : added control of NR filter choices in command file                Allgeyer 180608
!            PLOTLIB               : write plt*.kba containing prefit and postfit kbra                 Allgeyer 180615
!            ADDNORM               : record regularisation file and scale in addnorm.out               McQueen  180615
!
!  8.00  read L1B ascii files for GRACE/GRACE FO/GRACE II
!            input_openAllFiles    : adapt for GRACE/GRACE FO. Added read of LRI1B for GRACE FO 
!                                    and defined the RL_NUM for each mission                          Tregoning 180625
!            INPUT.H90             : added unit number for LRI. Removed RL_NUM as a parameter         Tregoning 180625
!            GRACEFIT_MOD          : Added "mission" variable. Changed SAT_A, SAT_B to SAT_1, SAT_2   Tregoning 180625
!            COMMAND_readSETUP     : Added read of mission. Made satellites mission-dependent         Tregoning 180625
!            INPUT_openALLFILES, GRACEFIT_MODS : declare L1B file name variables in GRACEFIT_MODS     Tregoning 180625
!            INPUT_openH5GTORB, INPUT_openGTORB     : change SAT_A/SAT_B to SAT_1/SAT_2               Tregoning 180625
!            plot_writeKB, plot_writeKBA : changed format statement for output max residual           Tregoning 180626
!            output_writeFIT       : changed format statement for output postfit RMS                  Tregoning 180626
!            ADDNORM               : added Tikhonov regularisation as an option                       Tregoning 180626
!            INPUT.h90             : changed LUGT to LUGTORB                                          Tregoning 180723
!            command_readsetup     : fixed bug in passing "mission" variable                          Tregoning 180723
!   HEADER_readGTORB,header_readH5GTORB : replace LUGT with LUGTORB                                   Tregoning 180723
!   input_openGTORB, input_readGTORB : replace LUGT with LUGTORB                                      Tregoning 180723
!   GRACEFIT                       : fix logic around the use of the regularisation file              Tregoning 180723
!   shadow_leadsat                 : fixed logic and variable names                                   Tregoning 180723
!   input_readGNV1B                : fixed bug in storing number of GPS observations                  Tregoning 180723
!   output_writeFIT                : removed debug statement                                          Tregoning 180723
!   add_meanACC_cond               : removed debug statement                                          Tregoning 180723
!   command_readSetup              : removed debug statements                                         Tregoning 180723
!   GRACEFIT_MOD                   : add acc_obs_model array                                          Tregoning 180726
!   command_readSetup              : add "acc_obs_model:" for linearising the accelerometer obs       Tregoning 180726
!   input_readACC1B                : added modelling of non-linear behaviour of ACC1B data            Tregoning 180726
!   input_readTHR1B                : fixed bugs in logic of setting thrust flags                      Tregoning 180726
!   input_openAllFiles             : fixed KBR1B filename for GRACE FO files                          Tregoning 180727
!   gnv_read_hdr                   : added code to read start/stop times for GRACE data               Tregoning 180802
!   kb_filtered                    : added STATUS comment regarding reconstruction of filtered data   Tregoning 180803
!   readbingtorb                   : changed length of character print for mascon file names          Tregoning 180814
!   input_readACC1B                : add number of columns to read to arguments to read_acc_data      Tregoning 180814
!   GRACEFIT,kb_filtered           : fixed bug to now include unfiltered kb obs in postfit residuals  Tregoning 180817
!   GRACEFIT,kb_computeKBRA        : add option of double differentiation of range resids for kbra    Tregoning 180823
!   kb_computeKBRA                 : fixed bug in variable type for use_kb_nr                         Tregoning 180828
!   plot_writeKBA                  : fixed two missing header lines                                   McQueen   180828
!   ADDNORM                        : changed addnorm.out to addnorm.fit. Added new 1st header line.   Tregoning 180831
!   output_writeVCV                : Use calling_prog value plus write "VCV" to 1st header line.      Tregoning 180831
!   output_writeFIT                : added new 1st header line to make compatible with VCV file.      Tregoning 180831
!   kb_filtered                    : fixed bug when integration starts in data gap                    Tregoning 180905
!   input_readGTORB                : add calling program to gtorb_record_length arguments             Tregoning 180907
!   ADDNORM                        : increased length of regularisation file to char*80               Tregoning 180910
!   ADDNORM                        : minor mods to output statements. Remove debug.                   Tregoning 180910
!   output_writefit                : re-inserted status message for pos/vel postfit RMS               Tregoning 180925
!   input_readTHR1B                : fixed bug when thrust occurs in first 2 secs of day              Tregoning/McGirr 180927
!   command_readSetup              : remove obsolete warning about range acceleration obs             Tregoning 181003
!   makefile                       : restored GRACESIM as a program in this directory                 Tregoning 181119
!   GRACESIM                       : fix bugs in calling arguments to plot_writeKB(A)                 Tregoning 181128
!   GRACESIM                       : removed NRD on simulated KBr to remove spikes?                   HM,SA     181215
!   GRACEFIT,WRITENORM             : add observation uncertainties to .norm header information        Tregoning 190222
!   ADDNORM                        : fixed format bugs. moved "VCV" across one column in header line  Tregoning 190301
!   GRACEFIT, GRACESIM             : fix bug in calculating mean beta angle                           Tregoning/McQueen 190301
!   header_readbinGTORB            : improved error message when inconsistent mascon numbers          Tregoning 190312
!   shadow_lib                     : removed legacy comma before i/o list in write statement          McQueeen  190424
!   GRACEFIT,kb_filtered_v2        : add new capability to apply EMD filtering to sub-segments of KBRR/KBRA data Tregoning 190502
!   input_openAllFiles             : change GRACEFO RL to RL04 and accelerometer file to ACT1B        Tregoning 190526
!   input_readGNV1B                : pass variable epoch_interval to gnv_read_data                    Tregoning 190529
!   kb_filtered_v2                 : fix bug to now filter final segment of data                      Tregoning 190531
!   kb_filtered_v2                 : fix bug in filtering the last segment of kbra data. Also, reduced
!                                    the max gap width down to 10 epochs (not sure why!)              Tregoning 190531
!   GRACEFIT,header_lib            : changed mascon number output from 4 to 5-digit zero-padded       McGirr    190715
!   GRACEFIT,command_lib           : added cosine fft filtering option for kband obs                  Tregoning/McGirr 190722
!   gracefit_mod                   : added use_kb_cos_fft variable                                    Tregoning/McGirr 190722
!   gracefit,kb_lib,kb_sponge      : pass fitted model param values back from kb_sponge subroutine    Tregoning 190807
!   kb_cos_fft                     : add interpolated data in gaps, extrapolated obs at each end      Tregoning 190809
!   kb_cos_fft                     : remove offset between obs and extrapolated/interpolated obs      Tregoning/McGirr 190812
!   kb_cos_fft                     : added call to window_function and mean calc before compute_fft   McGirr    190816
!   gracefit                       : replaced LS_normsolve with Chol_normsolve                        Allgeyer/Tregoning 190820
!   ADDNORM                        : replaced LS_normsolve with Chol_normsolve                        Allgeyer/Tregoning 190820
!   ADDNORM                        : turned OFF output of the VCV and corrl matrix                    Tregoning 190821
!   input_openAllFiles             : stop fatally if L1B_files.txt doesn't exist                      Tregoning 190903
!   kb_cos_fft                     : fixed bug (wrong pointer used for "infilled" array)              Tregoning 190903
!   gracefit,command_lib,header_lib,output_lib : change declaration/dimensioning for apriori,soln,prmnam arrays to
!                                                read them from mod_subs/soln_mod                     Tregoning 190905
!   output_writeNORM               : define names of satellites properly by mission                   Tregoning 190906
!   input_openh5                   : removed debug                                                    Tregoning 190917
!   output_writeFIT,gracefit       : write name of regularisation file to .fit file header            Tregoning 190917
!   header_readh5GTORB             : removed debug                                                    Tregoning 190917
!   gracefit                       : leadsat 4th args is nepochs_t not num_gps_used                   HM,SA     190924
!   addnorm                        : Change VCV row/col variable names and loop over rows, then columns Tregoning 191029
!   LS_lib/Chol_normsolve          : return the full VCV matrix instead just the Upper part           Allgeyer  191029
!   gracefit                       : fixed bug in beta angle calculation                              McGirr    191112
!   plt_lib...                     : Modify header Ax, Ay => lat, lon of the satellite                Allgeyer  191118
!   plt_lib...                     : add "#" symbol in front of header lines                          Allgeyer  191118
!   output_lib                     : Reformat header ... starting with "#"; Uniformized R; RR; RA     Allgeyer  191118
!   output_lib                     : Correct header MC (m) not (Gt)                                   Allgeyer  191118
!   header_lib                     : Correct header MC (m) not (Gt)                                   Allgeyer  191118
!   output_writeFIT                : added missing "#" to line "Postfit RA RMS"                       Tregoning 191120
!   command_readSetup, gracefit_mod : added "gracefit_step" as optional variable in gracefit.cmd      Tregoning 191126
!   gracefit                       : use "gracefit_step" variable to stop at various places in the run Tregoning 191126
!   gracesim                       : various changes to make it compatible with gracefit              Tregoning 191230
!   gracesim,input_readGTORBh5     : fixed bugs in reading h5 truth orbits                            Tregoning 200102
!   command_storeCommand           : improve error message when insufficient values requested         Tregoning 200110
!   addnorm                        : include output file character stem on command line               Tregoning 200108
!   addnorm                        : change mass conservation to include consideration of msc areas   Tregoning 200112
!   gracefit                       : Preprocessor flags to remove some logging information in prod.   Allgeyer  200130  
!   read_msc_ts_WRMS               : modified to read cosine/sine amplitudes rather than annual amplitude Tregoning 200225
!   gracefit/input_lib             : modified subroutine to read the GTORBh5 files directly in the part 
!                                    array (lots of memory gain)                                      Allgeyer (from) 200401
!   gracefit                       : Optimisation of the stacking loop (halved the stacking time)     Allgeyer (to)   200519
!   gracesim                       : update to gracefit standards for memory saving and speed         Tregoning 201006
!   kb_computeSim                  : change dimensioning of truthvec array                            Tregoning 201007
!   input_readGTORBh5_v1           : for a TEST ONLY, turn of differentiation of rotation matrices    Tregoning 201007
!   kb_computeKBRA                 : fixed bug in dimension of rvec                                   Tregoning 201007
!   gracesim                       : updated to chol_normSolve. Comment out gobs_computePartials      Tregoning/McGirr 201009
!   kb_lib                         : added comment to top of file                                     Tregoning 201012
!
!
!  9.00  "part" matrix renamed to "Amat" and converted to a two-dimensional matrix
!   input_readKBR1B                : added new array, use_obs, and set to T if kbr range obs exist    Tregoning 201014
!   gracefit_mod                   : added new array, use_obs. Declared as allocatable.               Tregoning 201014
!   kb_computeKBRA                 : update to use Amat array and use_obs array                       Tregoning 201014
!   kb_cos_fft                     : pass in a vector of use_obs                                      Tregoning 201014
!   kb_computeKBRA                 : remove setting kbra obs to "use" by default                      Tregoning 201109
!   mask_outliers                  : new routine to read epochs to remove from file outliers.mask     Tregoning 201109
!   makefile                       : added new subroutine mask_outliers.f90                           Tregoning 201109
!   gracefit                       : removed debug, commented out/removed shadow stuff, added in the
!                                       removal of outliers identified in file outliers.mask          Tregoning 201109
!   kb_cos_fft                     : turned off debug ".... II "                                      Tregoning 201110
!   mask_outliers                  : improve status messages                                          Tregoning 201112
!   gracefit, plot_writeKB         : calculate and output unfiltered postfit kbrr residuals           Tregoning 201113
!   gracefit                       : increased formatting size of output computing time               Tregoning 201114
!   gracesim                       : many changes to make it compatible with the Amat gracefit versn  Tregoning 210131
!   gracesim                       : add option of reading vector of diag elements for reg file       Tregoning 210131
!
!   addnorm                        : add flag to calculate and output singular values using DGESDD     McGirr   210216 
!   addnorm                        : add feature to output stacked normal equations files (snorms)    Tregoning 210217
!   output_codes2labels            : made more generic to handle snorm files with more than 1 day of ICs Tregoning 210217
!   addnorm                        : changed cons. of mass equation to multiply by msc_scl_factor rather
!                                    than the number of files                                         Tregoning 210222
!   output_writeNORM               : add IC_offset and msc_offset to call to output_codes2labels      Tregoning 210223
!   gracefit                       : added ability to read diag_only reg file                         McGirr    210415
!   gracefit                       : fixed bug in calc of postfits from vcv                           McGirr    210415
!   gracefit                       : added ability to read gracefit or addnorm style vcv              McGirr    210415
!   gracesim                       : optionally takes .kb as 12th arg and adds postfits to pre_omc    McGirr    210416
!   addnorm                        : performs SVD and outputs file containing USVt and AtWb           McGirr    210422
!   addnorm                        : comment out duplication of IC/mascon normeq into lower triang    Tregoning/McGirr 210521
!   mask_outliers                  : remove "s" from outlier.mask                                     Tregoning 211006
!   gracesim                       : added outlier mask code to gracesim                              Tregoning 211022
!   gracefit_mod, gracefit, gracesim : control mask code via entry in gracefit/sim command file       Tregoning 211206
!   addnorm_mod,read_addnorm_cmd,addnorm   : new handling of reading information from command files   Tregoning 211206
!   kb_kbsponge                    : turned off debug                                                 Tregoning 220121
!   read_addnorm_cmd               : fix bug when a normeq file was commented out                     Tregoning 220124
!   output_openAllFiles            : commented out the opening of the .svs, .rms, .resid files        Tregoning 220125
!   command_readSetup              : initialise tmp_end_neq to zero                                   Tregoning/McQueen 220125
!   gracefit                       : pass calling_prog into generate_inert_efixed                     Tregoning 220127
!   output_writeVCV                : comment out the writing of the full VCV                          Tregoning 220127
!   gracefit                       : improve the wording of the comments when forming matrices        Tregoning 220127
!   kb_cos_fft                     : added Lachlan's fixes to infilled data. Fixed bug in going beyond end_neq  Tregoning 220127
!   addnorm                        : increased dimensioning of max number of epochs                   Tregoning 220509
!   addnorm                        : fix bug when "num_files" is commented out in command file        Tregoning 220510
!   gracefit,mask_outliers         : include GPS PO and VE as obs. Mask out RA if RR masked           Tregoning 220518
!   kb_computeKBRA                 : fixed bug in dimensioning of use_obs                             Tregoning 220524
!   addnorm                        : scale by number of days in a snorm file (-99.0 in cmd file)      Tregoning 220527
!   read_msc_ts_WRMS               : changed dimensioning  to msc_const(nparam_t,nparam_t)            Tregoning 220608
!   addnorm_mod                    : added logicals output_correl and use_correl_in_reg               Tregoning 220608
!   addnorm                        : add code to control outputting the correlation matrix            Tregoning 220608
!   addnorm                        : add code to use mascon correlations when building regularisation Tregoning 220608
!   addnorm_mod                    : added logical use_adaptive_reg                                   Tregoning 220620
!   read_addnorm_cmd               : add code to read use_adaptive_reg: option                        Tregoning 220620
!   addnorm                        : add adaptive regularisation capability, controlled by cmd file   Tregoning 220620
!   mask_outliers                  : remove 3 RA obs either side of a RR data outage                  Tregoning 220624
!   kb_cos_fft                     : fix bugs in infilling data when gap goes to end of orbit         Tregoning 220624
!   addnorm                        : improved error message info when opening files                   Tregoning 220701
!   kb_cos_fft                     : fixed bug when start_neq > 1 affecting infilling                 Tregoning 220706
!   output_writeFIT                : changed format of RMS print statements (for Herb)                Tregoning 220801
!
! 10.0   :: remove all dependencies on gamit routines                                                 Tregoning 220804
!   read_addnorm_cmd               : add control for netcdf normal equations files                    Tregoning 220824
!   addnorm_netcdf, gracefit, gracesim : modifications to replace norm_netcdf_mod with soln_mod       Tregoning 221122
!   gracefit,gracesim              : write out normal equations in netcdf format                      Tregoning 221122
!   mask_outliers                  : removed debug                                                    Tregoning 221129
!   read_addnorm_cmd               : read the adapt_sigma from the addnorm command file               Tregoning 230213
!   addnorm_netcdf                 : add code for apriori mass in terms of sea level                  Tregoning 230223
!   addnorm_netcdf                 : fix bug in AtWb part of conservation of mass equation            Tregoning 230224

RETURN

END SUBROUTINE GVERSN
