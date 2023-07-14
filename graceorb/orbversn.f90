
subroutine orbversn

  call status_update('STATUS','GRACEORB','orbversn',' ','graceorb: Integrate GRACE orbits and compute mascon and IC partials',0)

! List all the modifications made to the graceorb software.
!
! File started 1 November 2012
!
!     GRACEORB: renamed orbit.f90 to GRACEORB.  Tregoning 121101
!               removed reading of ut1. pole. etc filenames. Code just reads in default gamit names   Tregoning 121101
!     INRED   : skip over the entries DEALIAS_DIR, DATA_DIR and MASCON_DIR in GRACE.input. Also, increased "MASCON" to "MASCON " to separate it from "MASCON_DIR" and the same for DEALIAS  Tregoning 121101
!     Makefile: added orbversn to makefile  Tregoning 121101
!     ADAM_MOULTON_PART:  added calls to status_update    Tregoning 121102
!     JPLEPHEM  :  fixed declarations of integer variables  Tregoning 121211
!     INRED     :  removed duplicate declaration of variable "message"  Tregoning 121211
!     GRACEORB, GRAVCALC, ADAM_MOULTON_PART  : merged code to permit mascons to be turned off  Purcell/Tregoning 121221
!     GRAVCALC  :  fixed error in computation of longitude for mascon ternary elements   Purcell 121204
!     ADAM_MOULTON_PART: fixed bug in scale partial  Purcell 130328
!     GRACEORB, ORBITCALC, HDWRT :  added GPS antenna offset parameters to input/output lines  Tregoning 130528
!     PEPtime, JPLEPHEM: moved these routines to grace/lib  Tregoning 130620
!     ACCELEROM, SCArot, ADAM_MOULTON_PART : made usage of quaternions consistent in sign Purcell/Tregoning 130710
!     GRACEORB:  added a priori mascon values to the end of the GTORB file
!     ACCELEROM: average the accelerometer obs over 5 seconds rather than just using 1st value  Tregoning 130815
!     PLANETFIELD: fixed bug in indirect J2 contribution (missing square on z/r term)   Tregoning 130830
!     INMOD_MOD, INRED, GRACEORB: bit-map which satellite IC parameters to update a priori    Tregoning 130905
!     OCEANFIELD: fix bug in normalisation of ocean tide harmonic coefficients   Tregoning 130905
!     OTMGRD2ARR: implemented openMP over lat/lon loop for ocean tide grid (8% time saving)  Tregoning 130905
!     GETGRID: define ocean tide grid file name based on model specified in input file  Tregoning 130906
!     ACCEL_MOD, MASCON_CALC, GRACEORB : change mcon_lat variable to "mcon_colat" to represent what it is   Tregoning 130918
!     HDWRT : added extra decimal places in scale/bias a priori values in header  Tregoning 131016
!     GRACEORB : fixed name of IC file in error statement if file not opened   Tregoning 131016
!     POLTID : fixed bug when we roll over to the next PeP day  Tregoning 131030
!     GRACEORB, HDWRT_BINARY (new routine), ADAM_MOULTON_PART : changed GTORB files to be binary  Tregoning 131101
!     WHICH_ICPARAM : new routine to find the index for a particular parameter in the input vcv file   Tregoning 131107
!     GRACEORB : modified to read in an incomplete list of IC parameters and assign/update a priori values only for those present  Tregoning 131107
!     HDWRT_BINARY: output the name of the mascon_primary_file onto the header line "MODELS USED"  Tregoning 131121
!     HDWRT_BINARY: increased the significant figures of a priori bias values   Tregoning 131205
!     SPHHARMFIELD: set evaluation of dealiasing file back to starting at degree 2  Tregoning 140110
!     TIDEGRID    : save variables for subsequent use when tide grid not updated    Tregoning 131202 (reinserted into code 140110)
!     GENERALREL  : fixed bug in assigning sun velocity to array velS (1 cm effect over 12 hrs)  Tregoning 140130
!     GRAVCALC, SPHHARMFIELD, GRACEORB : compute ocean tides by either spherical harmonics or grid   Tregoning 140204
!     GENERALREL  : added relativity accelerations from Earth oblateness (2mm effect over 12 hours)  Tregoning 140205
!     GRAVCALC    : passed J2 through to GENERALREL   Tregoning 140205
!     ROTATION_MOD: fixed dimensioning bug on variable sidtm   Tregoning 140218
!
! Version 2.0
!     GENERATE_INERT_EFIXED : new routine (in inert_lib.f90) to compute the inertial->efixed rotation matrices (and their rates)
!                             for all epochs of the requested orbit, plus 5 mins at the start/end
!     READ_EOP_UT1 : read the usno.finals.data file to get pole and ut1-utc information
!     IAU_INTERP, LAGINT, PMUT1_OCEANS, PM_GRAVI   : subroutines (by Ch. Bizouard) to interpolate pole/UT1-UTC to a requested
!                 epoch and to add the sub-daily signals to it
!     ORBITCALC : replaced calls to rotsnp with new routines in inert_lib       Tregoning 140318
!     ADAM_MOULTON_PART : replaced calls to rotsnp with new routines in inert_lib       Tregoning 140318
!     EFIXED_INERT : pass in the rotation and rotation rate matrices, then compute.     Tregoning 140318
!     ICPART_CALC  : moved out of adam_moulton_part into its own file. New inert-efixed code.  Tregoning 140318
!     SOLIDFIELD   : passed in the inertial->efixed rotation matrix    Tregoning 140318
!     POLTID       : updated mean pole to IERS 2010 standards   Tregoning 140318
!     GENERATE_INERT_EFIXED, ORBITCALC, ADAM_MOULTON_PART, GRAVCALC, SPHHARMFIELD, POLTID : passed the Xp, Yp 
!                    values through to correct C21/S21 for pole tide  Tregoning 140319
!     GENERATE_INERT_EFIXED :  add the computation of the rotation acceleration matrix  Tregoning 140320
!     NTERP_INERT :  add the interpolation of the rotation acceleration matrix  Tregoning 140320
!     ORBITCALC, ADAM_MOULTON_PART, GRAVCALC : include the rotation acceleration matrix  Tregoning 140320
!     GENERATE_INERT_EFIXED :  add the computation of the roll/pitch/yaw values and their rates  Tregoning 140320
!     GRACEORB, GENERATE_INERT_EFIXED : apply 3 x rotation angle adjustments from command line to the roll,pitch,yaw of
!                    the inertial-efixed rotation matrices (for testing)  Tregoning 140328
!     GENERALREL : fixed bug (introduced by Tregoning 140205) in Schwarzchild computation (vdotv, not v) McClusky/Tregoning 140302
!     OUTPUT_LIB : increased significant figures on output prefit pos/vel RMS values   Tregoning 140410
!     GRAVCALC   : add computation of partials of mascon accelerations wrt position    Tregoning 140428
!     GRAVCALC   : fixed bug when using mascons but not ocean tide grid    Tregoning 140505
!     GENERATE_INERT_EFIXED :  fixed bug in looping through too many epochs   Tregoning 140514
!     ADAM_MOULTON_PART : changed starting partials for velocity to use rot_e2i (as it should be), set partial of position wrt vel to zero Purcell/Tregoning 140519
!     GRACEORB   : add 1000 m to a priori values of ocean mascons so that gracefit can identify them  Tregoning 140527
!     GRACEORB/GRAVCALC/MASCON_CALC/SETUP_OCEAN_TIDE_GRID : changed logic in handling ocean tide from binary grid  Tregoning 140606
!     ACCEL_MOD : changed mascon variables definitions
!
! Version 3.0  add estimation of tidal amplitudes for ocean mascons
!     ACCEL_MOD, GRACEORB : added bit-mapped flag for tidal mascons  Tregoning 140612
!     INMOD_MOD,  : add code to read from command file whether to estimate tidal mascon amplitudes   Tregoning 140612
!     HDWRT_BINARY : wrote bit-mapped tidal information for each mascon to the line after the end of the GTORB header  Tregoning 140613
!     ACCEL_MOD : added counter of number of tidal mascon constituents to estimate  Tregoning 140613
!     HDWRT_BINARY : wrote number of mascons and tidal mascon constituents to header  Tregoning 140613
!     HDWRT_BINARY : write out the actual name of the GRACE.input file (not a dummy, hardwired name)  Tregoning 140616
!     GRACEORB     : write out every a priori mascon tidal amplitude to the GTORB file (for all mascons, one component per constituent per record)  Tregoning 140619
!     GRACEORB     : fixed bug in defining period of mascon tidal constituents     Tregoning 140707
!     GRAVCALC     : initialised the mascon tidal partials to zero     Tregoning 140707
!     ICPART_CALC  : fixed bug in incrementing mascon tidal partials in a0 array  Tregoning 140708
!     GRACEORB     : repositioned the counting of the number of tidal amplitudes to be modelled   Tregoning 140722
!     ADAM_MOULTON_PART : increased the dimensioning of maxnd to allow for at least 820 tidal mascon amplitudes   Tregoning 140722
!     ORBITCALC    : increased the record length of the GTORB to allow up to 820 tidal mascons ( < 61 south)      Tregoning 140722
!     GRACEORB, GRAVCALC : added satellite name to status outputs  Tregoning 140724!
!
! Version 4.0  change the way in which we generate partial derivatives of the IC parameters
!     ADAM_MOULTON_PART : integrate perturbed orbits then difference to get partials (rather than integrating the
!                         perturbation)    Tregoning 140812
!     CALC_PERTURB_ACC  : new subroutine to adjust the total inertial acceleration for the perturbed orbit locations   Tregoning 140813
!     ACCELEROM : pass in bias/scale information as arguments rather than via the include file    Tregoning 140813
!     BIASSCALE, ORBITCALC, ACCELEROM : passed the bias values through as arguments rather than via the include file   Tregoning 140813
!     ACCELEROM, GRAVCALC, CALC_PERTURB_ACC : change the way the bias perturbation is passed in so that the computation of bias partials are not affected by scale errors   Tregoning 140816
!     ADAM_MOULTON_PART : fixed bug in time variable passed to calc_perturbed_acc   Tregoning 140817
!     ADAM_MOULTON_PART : changed the size of the perturbation for scale-perturbed orbits  Tregoning 140819
!     BIAS_FROM_ACC     : fixed status_update statement so that adjusted biases were output  Tregoning 140820
!
! Version 5.0  add a twice-per-rev parameterisation to the along-track acceleration
!     ADAM_MOULTON_PART, CALC_PERTURBED_ACC : increase number of ICs by two (2pr cos, 2pr sin)  Tregoning 140820
!     GRACEORB, ORBITCALC, HDWRT_BIN : increased efic from 6 to 8 ICs (pos/vel/2pr)             Tregoning 140820
!     HDWRT_BINARY : added 2prS/C to binary header comment information     Tregoning 140821
!     ADAM_MOULTON_PART, CALC_PERTURBED_ACC, GRACEORB, ORBITCALC, HDWRT_BINARY : added 1/rev along-track acc   Tregoning 140821
!     WHICH_ICPARAM : increased max number of IC parameters to 34   Tregoning 140822
!     BIAS_FROM_ACC : increase mean bsz by +5 nm/s^2 (rather than reduce by 30 nm/^2)  Tregoning 140822
!     output_CODES2LABELS : fixed bug in 1/rev and 2/rev label conversion   Tregoning 140822
!     SETUP_OCEAN_TIDES   : fixed bug in computation of seconds of day for hours > 12  Tregoning 140829
!     ADAM_MOULTON_PART   : increased number of patials output into GTORB to include 1/rev and 2/rev  Tregoning 140904
!     ICPARTCALC          : fixed bug in the index of the mascon partials to account for 1/rev 2/rev  Tregoning 140904
!     GMSET               : changed the GM value from 4415 to 4418, which is the IERS2010 value  Tregoning 150603
!     ICPart_calc         : replaced the local definition of GM(earth) with the value in gmset.f90  Tregoning/McClusky 150605
!     GRACEORB            : fixed bug in mascon index for reading in a priori values  Tregoning 150731
!     GRAVCALC            : fixed bug in indexing for computing mascon tidal amplitude partials  Tregoning/Purcell 150821
!     ICPART_CALC         : fixed bug in indexing for computing mascon tidal amplitude partials  Tregoning  150824
!     ADAM_MOULTON_PART   : increased maxnd to allow tidal amaplitudes on 203 mascons  Tregoning 150827
!     ORBITCALC           : increased record length of GTORB file to allow for 20 extra mascons with tidal amplitude estimates  Tregoning 150827
!     COEFF_READ          : add traps for case when input static grav file doesn't exist or contains no data   Tregoning 151014
!     GRACEORB            : turned off the read of roll/pitch/yaw errors on command line
!     POLTID              : fixed a sign bug in the ocean pole tide computation   Tregoning 160121
!     SPHHARMFIELD        : add in the ocean pole tide (as computed in poltid)    Tregoning 160121
!     GENERALREL          : fixed bugs in the J2 contribution calcs               King/Tregoning 160202
!     ACCEL_MOD           : added array to store the min dist from each mascon to the satellite within an orbit   Tregoning 160205
!     GRACEORB            : fixed indexing bug when transferring a priori mascon information   Tregoning 160209
!     GRACEORB, INRED     : bit-map the use of a priori tidal amplitudes     Tregoning 160215
!     WHICHMSCTIDE_PARAM  : increased max loop counter value   Tregoning 160321
!
! Version 6.0  Revamped mascon and tide code. Moved a lot of stuff from graceorb.f90 to subroutines
!    INMOD_MOD            : added new file name for the combined mascon file   Tregoning 160905
!    GRACEORB             : moved a lot of code reading input information into subroutines in input_lib.f90  Tregoning 16090[2-6]
!    READ_GRACEORB_INFILE, READ_GRACEORB_ic, INRED_PRINT, READ_DEALIAS, READ_MASCON_FILE,     
!                           code moved from graceorb.f90 into new subroutines   Tregoning 16090[2-6]
!    ORBITCALC, INERT_LIB, USNO_MOD : moved dimensions for max and actual number of eop values into a mod file  Tregoning 161005
!    GRAVCALC, SPHHARMFIELD, POLTID, ADAM_MOULTON_PART, ORBITCALC, INERT_LIB, CALC_PERTURB_ACC : moved Xp,Yp,ut1utc from
!                           passed variables into usno_mod.f90   Tregoning 161005
!
! Version 6.1  Revamp dealias file format and add atmospheric tides
!    READ_DEALIAS         : read the header line, which gives the max degree of the AOD1B file   Tregoning 161018
!    SPHHARMFIELD, ATM_TIDE_SPH  : add code to include the atmospheric tides   Tregoning 161018
!    ATM_TIDE_SPH, OCENFIELD_SPH : added sat_mod so that we have the name of the satellite   Tregoning 161019
!    GENERATE_MASCON_VECTOR      : turned "sat_lat" into latitude (as it should be) not colatitude   Tregoning 161025
!    GRACEORB             : read mascon file header    Tregoning 161124
!    READ_PRIM_DIST_FLAGS : read in the header lines. Compare the code in this file with the unique code from the mascon file
!                           to ensure compatibility.    Tregoning 161124
!    ADAM_MOULTON_PART    : added allocation statement to dimension efacc_part   Tregoning 161128
!    GRACEORB             : updated to ellipsoidal geometry pattern for ternary mascons  Tregoning/McQueen 161201
!    GRACEORB             : rearranged reading of mascon headers and allocating mascon arrays      Tregoning 161209
!    CALC_OCEAN_ACC_PART  : added (back) in the computation of the ocean tidal amplitude partials  Tregoning 170125
!    WHICH_MSCTIDE_PARAM, GRACEORB  : updated to using mascon_mod variables related to ocean mascon file       Tregoning 170127
!    GRAVCALC,CALC_OCEAN_ACC_PART   : added jd to argument list to compute the dt_tides variable value         Tregoning 170127
!    HDWRT_BINARY         : added mascon file names and hdr codes to GTORB header, fixed some formats          Tregoning 170127
!    ACCEL_MOD            : moved max_msc_tides declaration to lib/mascon_mod.f90                              Tregoning 170127
!    ICPART_CALC          : updated variable names for tidal amplitude partials                                Tregoning 170127
!    TERNARY_TIDAL_AMPL   : new routine to separate the ternarys into those with/without tidal amplitudes      Tregoning 170130
!    CALC_OCEAN_ACC_PART  : finished adding back the ocean tidal amplitude modelling and partials              Tregoning 170131
!    GRACEORB             : fix logic to allow ocean tide mascon accelerations when temporal mascons are off   Tregoning 170208
!    GRACEORB             : fixed bug initialising mascon a priori values                                      Tregoning 170210
!    ICpart_calc          : fix bug in variable name in if statement for tidal amplitude partials              Tregoning 170228
!    ICpart_calc          : output tidal ampl partials irrespective of whether using apriori values or not     Tregoning 170303

! Version 7.  Create "grace/active/graceorb version.
!    MAKEFILE             : removed oceanfield and all references to spotl libraries and routines              Tregoning 170314
!    MAKEFILE             : moved cart_convert to lib/cart_convert                                             Tregoning 170316
!    ADAM_MOULTON_PART, ATM_TIDE_SPH, COEFF_READ, GRAVCALC, LEGCALC, LEGCALC_NORM, OCEANFIELD, ORBITCALC
!       OCEANFIELD_SPH, SOLIDFIELD, SPHHARMFIELD, STATICFIELD, TEST_CHECK_TIDE, TEST_OMP, 
!       TEST_PLM :  renamed spherhar_mod1 to spherhar_mod                                                      Tregoning 170321
!    GENERATE_INERT_EFIXD : fix bug in dimensioning of RPOM matrix                                             Tregoning 170322
!    INERT_LIB            : moved this routine to grace/lib (removed from makefile)                            Tregoning 170327
!    ACCELEROM, ORBITCALC, ADAM_MOULTON_PART, ICPART_CALC, GRAVCALC, CALC_PERTURBED_ACC
!                           moved efixed-inertial matrices into mod_subs/rotation_mod.f90                      Tregoning 170403
!    ACCELEROM, ADAM_MOULTON_PART, GRAVCALC, CALC_PERTURBED_ACC : add partials for roll/pitch/yaw              Tregoning 170405
!    INRED                : increased length of character variable that reads the GRACE.input file             Tregoning 170510
!    ADAM_MOULTON_PART    : increased array dimensions to allow for more mascon parameters (currently 7526)    Tregoning 170609
!    ORBITCALC            : increased record length of GTORB file to allow for more mascons                    Tregoning 170609
! 
! Version 8.  GTORB record length information sourced from mod_subs/gtorb_mod and written to GTORB file       
!    GRACEORB, ORBITCALC, ADAM_MOULTON_PART : include gtorb_mod (and stop passing "nrec" around)               Tregoning 170609
!    HDWRT_BINARY         : write the GTORB record length to the first line of the GTORB file                  Tregoning 170609
!    HDWRT_BINARY         : add total_ocean_prim to DATA RECORD FORMAT line, since it is required in gracefit  Tregoning 170611
!    ORBITCALC, GENERATE_INERT_EFIXED : added rotaccs, being the time derivative of rotdots                    Tregoning 170821
!    INERT_INTERP         : added rotaccs, rotacc_e2i, rotacc_i2e                                              Tregoning 170821
!    ORBITCALC, ADAM_MOULTON_PART  : added accelerations to inert_interp calls                                 Tregoning 170821
!    ADAM_MOULTON_PART, EFIX_iNERT : add acceleration transformations efixed<->inert                           Tregoning 170821
!    ORBITCALC                     : fixed bug (missing argument) in call to generate_inert_efixed             Tregoning 170821
!    IAU_INTERP                    : removed from makefie. This is already in ../lib                           Tregoning 170821
!    ICpart_calc, calc_perturb_acc : added rotation acceleration arguments to call to efixed_inert             Tregoning 170821
!    GRAVCALC             :  added trap for when radius of satellite orbit is < earth radius                   Tregoning 171003
!    ORBITCALC, GRAVCALC  : only read anc compute accelerometer/star camera data if using accelerometer obs    Tregoning 171003
!    ORBITCALC            : fatal stop if need to calculate default acc bias but didn't read acc data          Tregoning 171003
!    CALC_PERTURBED_ACC   : do not call accelerom subroutine if we are not using accelerometer obs             Tregoning 171003
!    ADAM_MOULTON_PART    : fixed a couple of spelling mistakes                                                Tregoning 180216
!    ATM_TIDE_SPH         : fixed bug in initialising the read of coefficients for 1-second integration        Tregoning 180302
!
! Version 9. read generic L1B ascii data files for ACC1B and SCA1B rather than in-house versions
!    read_graceorb_infile : read additional lines in input file (how to model non-linear accelerometer data)   Tregoning 180619
!    read_graceorb_IC     : no longer read polynomial coefficients for bias values                             Tregoning 180619
!    ORBITCALC            : moved call of read ACC1B and SCA1B data to graceorb.f90                            Tregoning 180620
!    BIAS_FROM_ACC        : changed format statement to output mean bias values                                Tregoning 180621
!    INRED                : added GRACE C and GRACE D GPS antenna offsets to list of permissible lines         Tregoning 180621
!    HDWRT_BINARY         : removed debug                                                                      Tregoning 180621
!    read_acc_sca_thr     : renamed, added the read of the THR1B information and removing a non-linear model   Tregoning 180705
!    read_acc_sca_thr     : added debug flag to call to acc_fit_quadratic                                      Tregoning 180727
!    read_acc_sca_thr     : added option of using EMD to linearise accelerometer obs. Changed counter names    Tregoning 180809
!    read_acc_sca_thr     : add number of columns of ACC1B data to be read                                     Tregoning 180810
!    read_acc_sca_thr     : pass the acc1b obs span (in seconds) to acc_EMD1                                   Tregoning 180816
!    read_acc_sca_thr     : added a call to the subroutine EMD2                                                McGirr    180914
!    acc_EMD2,acc_zero_crossing : limit data passed to acc_EMD2 if there are missing acc obs                   Tregoning/McGirr 181003
!    read_graceorb_infile : added reading the surface on which the mascons will be modelled (geoid/ell/topo)   Tregoning 181023
!    gravcalc,solidfield  : update to IERS2010 solid body tides                                                Tregoning 181029
!    read_acc_sca_thr     : added a call to the subroutine EMD3                                                McGirr    181203
!
! Version 10. : changed equation to calibrate accelerometer observations
!    calc_perturb_acc,gravcalc,accelerom : pass scale perturbation separately into accelerom                   Tregoning 190222
!    accelerom            : changed equation to calibrate accelerometer observations                           Tregoning/Purcell 190222
!    input_lib, orb_mascon_lib : fixed integer multiplication for large mascon files	 		       McGirr/Tregoning 190305
!    adam_moulton_part    : increased no of params/partials to be integrated to account for 48000 mascons      McGirr    190305
!    input_read_apriori   : updated to use mod_subs/soln_mod. Renamed "prm_input". Removed VCV_obs_local       Tregoning 190312
!    input_update_apriori_params : updated to use mod_subs/soln_mod. Removed "prm_input" to "prmnam"           Tregoning 190312
!    graceorb             : removed prm_input from calling arguments to input_update_apriori_params
!    graceorb             : moved declaration of GTORB_recl from mod_subs/gtorb_mod.f90 to graceorb.f90        Tregoning 190403
!    input_read_apriori   : updated to use revised lib/read_IC_list_v1 using mod_subs/soln_mod                 McQueen   190610
!    read_acc_sca_thr     : added a call to the subroutine fft                                                 McGirr    190701
!    graceorb,read_acc_sca_thr : passed in the position of the satellite                                       McGirr/Tregoning 190703
!    graceorb,read_acc_sca_thr : changed nfft to 2**20 for smoother spectrum                                   McGirr    190716
!    poltid               : updated mean pole tide model to the IERS 2018 update linear model                  Tregoning 190722
!    read_acc_sca_thr     : removed unnecessary print statements                                               Tregoning 190807
!    read_L1b_files       : added fatal stop if file "L1B_files.txt" doesn't exist                             Tregoning 190712
!    read_acc_sca_thr     : lowered high frequency cutoff send to acc_cos_fft                                  McGirr    190816
!    read_acc_sca_thr     : changed fhi and flo (again)                                                        McGirr/Tregoning 190819
!    bias_from_ACC        : trim 1000 epochs from beginning/end to remove fft edge effects                     McGirr/Tregoning 190820
!    graceorb,input_read_apriori   : add flag to indicate whether mascon apriori values have been read         Tregoning 190917
!    graceorb,read_GRACEinput : renamed "inred", added msc_apriori_file to GRACE.input file                    Tregoning 190917
!    INPUT_read_msc_apriori : detect the apriori mascon file type and stop if it can't be resolved             McQueen   190919
!    graceorb             : start storing mascon apriori values at index 39 of apriori and soln                McQueen   190925
!    read_acc_sca_thr     : fixed bug in high-pass filter frequency settings (from e-9 to e-3)                 Tregoning/McGirr 190927)
!    graceorb,input_read_apriori : mods to read only ICs from APRIORI_FILE, mascons from MSC_APRIORI_FILE      Tregoning 190927
!    input_read_msc_apriori      : allocate array apriori_msc, then read in the apriori and solution           Tregoning 190927
!    adam_moulton__part   :  Preprocessor flags to remove some logging information in production mode          Allgeyer  200130
!    orb_mascon_lib       :  Preprocessor flags to remove some logging information in production mode          Allgeyer  200130  
!    sphharmfield         : fixed bug in interpolating AOD1B coefficients                                      Tregoning/McGirr 200304
!    read_acc_sca_thr     : changed fhi and flo (again). Now 0.45d-4 and 0.55d-4                               McGirr    200512
!    accelerom,poltid     : undid changes that put code back to an old way of running (calib eqn, poltide)     Tregoning 200622
!    lib/msc_lib/calc_mascon_xyz : changes ocean ternary mascon CoM from ellipsoid to geoid. Not in this program
!                                  directory but affects how graceorb now works                                Tregoning 201031
!
! Version 11.0  :: remove all dependencies on gamit/globk routines
!    solidfield           : replaced "tide_angles" with "fund_arguments_iers1992"                              Tregoning 210615
!    graceorb             : add runstring information if no command arguments supplied                         Tregoning 210615
!    read_GRACEinput      : replaced all calls to trimlen with calls to "trim"                                 Tregoning 210716
!    graceorb             : add "O" option to gt_mcon, being to use mascon accels but not partials             Tregoning 210806
!    gravcalc             : add "O" option for mascon accels but no partials                                   Tregoning 210806
!    adam_moulton_part    : only declare array for mascon partials if we want to output mascon partials        Tregoning 210806
!    graceorb             : ensure soln array is allocated when we want apriori mascons but not ICs            Tregoning 210811
!    graceorb, INPUT_update_apriori_params, calc_mascon_acc_part: changed mcon_EWH to mcon_prim_EWH            Tregoning 210818
!    graceorb             : fixed logic around setting zero apriori mascon values when required                Tregoning 210823
!    INPUT_read_msc_apriori : fixed bug in assigning a priori mascon values                                    McGirr    211011
!    etides_iers2010      : fixed bugs in Doodson numbers                                                      Tregoning 211014
!    graceorb             : added call to etides_iers2010                                                      Tregoning 211014
!    atm_tide_sph         : fix bug on first read, when using 1 sec integration step                           Allgeyer/Tregoning 211221
!    read_prim_dist_flags : fixed bug in do loop counter. Also fixed name in status calls                      Tregoning/Purcell 220104
!    generate_mascon_vector : fixed indexing bug. Now reads lower-tri (no diag elements) properly              Tregoning 220106
!    orbitcalc            : pass program name into generate_inert_efixed                                       Tregoning 220127
!    makefile,sphharmfield: added Seb's code for the long-period ocean tides                                   Tregoning 220204
!    read_acc_sca_thr     : removed debug                                                                      Tregoning 220429
!    read_l1b_file        : removed debug                                                                      Tregoning 220429
!    INPUT_update_apriori_params : removed debug                                                               Tregoning 220429
!    sphharmfield         : changed name of mod file for ocean tides to omtides_mod                            Tregoning 220429
!    read_acc_sca_thr     : fix along-track thermal transfer functions for transplant data                     Tregoning 220429
!
! Version 12.0 :: removed (again) all dependencies on gamit/globk routines                                     Tregoning 220804
!    calc_mascon_part     : reduced from 20,000 to 1000 km the distance beyond which primary mscs are point sources Tregoning 221129
!    adam_moulton_part    : reduced number of partials orbits from 19 to 12 (no 1/rev, 2/rev, rpy)             Tregoning 221129
 
  end

