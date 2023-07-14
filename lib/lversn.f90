    SUBROUTINE LVERSN(VERSION)

!   Update history for GRACE common library grace/lib. (This subroutine is never called, it is just for documenting changes)

!   0.10     LVERSN:  documenting changes commenced  Tregoning 130620
!            JPLEPHEM: subroutine moved to lib from graceorb
!            GET_LAMBDA: subroutine from within gamit/arc/shadow.f extracted and placed in grace/lib
!            READ_SOLN_V2: added argument to pass out the ascii descriptors of the parameters   Tregoning 131107
!            TIDE_LIB/spotl_point_setup_v2: changed output debug to limit the amount   Tregoning 140106
!            INVERT : added simple matrix inversion routine   Tregoning 140318
!            TIDE_LIB/read_grid_header : changed to match new format of binary header file  Tregoning 140603
!            READ_SOLN_V2 : added pointers to first parameter entries for mascons and tidal mascons  Tregoning 140612
!            READ_SOLN_V3 : read subroutine for GRACEKAL, including pointers for first parameter entries   Tregoning 141125
!            TIDE_LIB/write_tide_grid_head : write out the lat/lon values of the mascons  Tregoning 150806
!            NRLMSIS_lib  : new subroutines to compute the NRLMSIS atmospheric density   Tregoning 151023
!            READ_SOLN_v3 : fixed output message related to number of tidal amplitudes   Tregoning 160429
!            CALC_WHICH_TERNARY : new routine to determine in which ternary mascon the satellite resides  Tregoning 160902
!            MASCON_MOD, READ_MSC_HDR, READ_OCEAN_MASCONS  : fixed dimensioning bug for ocean mascons     Tregoning 161209
!            ALLOCATE_MASCONS   : added allocation of ocean tidal amplitude partials array                Tregoning 170127
!            ALLOCATE_MASCONS   : increased the number of pointers in mcon_tern_ptr to 3(but didn't use the 3rd)  Tregoning 170131
!            READ_SOLN_V3       : fix bug in storing off-diagonal covariance terms of VCV matrix          Tregoning 170308
!            makefile           : moved graceorb/cart_convert.f90 to lib                                  Tregoning 170316
!            makefile           : moved graceorb/inert_lib.f90 to lib                                     Tregoning 170327
!            allocate_mascons   : added extra "geoid" column to prim/sec/tern mascon arrays               Tregoning 170530
!            read_mascon_file   : added "geoid" column to mascon arrays                                   Tregoning 170530
!            write_mascon_record: added msc_geoid to argument list. Added "geoid" to output mascon lines  Tregoning 170530
!            write_mascon_record: increased to i8 (from i7) the #mascons in prim and sec lines            Tregoning/McQueen 170531
!            read_ocean_mascons : removed unnecessary check on ternary mascon number                      Tregoning 170602
!            write_mascon_record: increased to f18.0 for area of primary and secondary mascons            Tregoning 170607
!          
!   0.20     gtorb_lib,makefile : new subroutine to determine the record length of a GTORB file           Tregoning 170609
!            READ_SOLN_V3       : increased significant figures in solution vector                        Tregoning 170713
!            NOISE_ROBUST_DERIV : moved from ../gracefit, and commented out gracefit.h90                  Tregoning 170818
!            MATRICES_LIB       : added new routine that contains transp.f90                              Tregoning 170821
!
!  
!   1.00     makefile           : renamed infill_SCA.f90 to infill_TS.f90 to make it generic              Tregoning 180620
!            ACC_READ_HDR, ACC_READ_DATA : new routines to read ACC1B data for all missions               Tregoning 180531
!            level1B_OPEN, SCA_READ_HDR, SCA_READ_DATA, THR_READ_HDR, THR_READ_DATA : L1B reading         Tregoning 180626
!            read_thr_hdr       : read the number of thrusts in GRACE FO data files                       Tregoning 180727
!            read_gnv_hdr       : read start/stop times for GRACE GNV1B. Read GRACE FO data.              Tregoning 180802
!            read_mascon_file   : made character length of mascon file name generic                       Tregoning 180808
!            acc_read_data      : make variable the number of obs/epoch read from the ACC1B file          Tregoning 180810
!            acc_EMD1,acc_inflexions : set the required number of inflexions based on obs data length     Tregoning 180816
!            gtorb_record_length: add fatal stop if GTORB file does not exist                             Tregoning 180907
!            acc_EMD2, EMD2     : added new subroutine calling the new EMD2                               McGirr    180914
!            acc_read_data      : fix bug when orbit length is shorter than ACC1B data length             Tregoning 180917
!            acc_inflexions     : added conditional to deal with extended acc data                        McGirr    180924
!            acc_zero_crossings : determine frequency content of EMD components by # zero crossings       Tregoning 180926
!            kbr_read_data      : transferred code from gracefit/input_readKBRR to new subroutine         Tregoning 181119
!            kbr_read_data      : fixed bugs in declaration of passed variables                           Tregoning 181120
!            acc_EMD3, EMD3     : subroutine calls emd_v3, stopping crietria determined by zero crossings McGirr    181203
!            deallocate_mascon_arrays : new subroutine. Name says it all ....                             Tregoning 181205
!            get_lambda         : fixed single precision pi                                               McQueen   181224
!            tide_lib           : fixed single precision pi                                               McQueen   181224
!            msc_lib            : fixed single percision pi                                               McQueen   181224
!            gtorb_record_length: increased check size on maximum record length of binary GTORB. Fixed some
!                                 minor error reporting information                                       Tregoning/McGirr 190306
!            read_soln_v1       : updated parameter dimensioning to use mod_subs/soln_mod.f90             Tregoning 190312
!            read_soln_v2       : updated parameter dimensioning to use mod_subs/soln_mod.f90             Tregoning 190312
!            write_soln         : updated parameter dimensioning to use mod_subs/soln_mod.f90             Tregoning 190312
!            write_soln_v2      : updated parameter dimensioning to use mod_subs/soln_mod.f90             Tregoning 190312
!            read_soln_v3       : updated parameter dimensioning to use mod_subs/soln_mod.f90             Tregoning/McGirr 190313
!       inert_lib/read_eop_ut1  : added fatal trap if bull_b values are missing from usno file            Tregoning/McQueen/Allgeyer 190326
!       bubble_sort_R8          : new routine, floating point bubble sort (also outputs index order)      Tregoning 190408
!       acc_lib,inert_lib	: removed legacy commas before i/o list	in write statements               McQueeen  190424
!       adjust_mascon_surface   : only makes geoid/ellipsoid correction to surface mascons (not GIA)      McGirr    190522
!       allocate_mascon_arrays  : update ocean mascon array allocations                                   Tregoning 190522
!       sca_read_hdr            : don't trust the yaml header for the number of observations!!!           Tregoning 190525
!       gnv_read_data           : changed logic to return only 5 s GNV1B data, irrespective of data in input file   Tregoning 190529
!       acc_lib, level_1B_lib   : read start/end seconds as real in read_hdr routines in case it's <10^8  McQueeen  190608
!       read_IC_list_v1         : updated to use mod_subs/soln_mod to allocate arrays                     McQueen   190610 
!       fft_lib/compute_fft     : FFT analysis, filters coefficients, FFT synthesis, spectra, freqs       McGirr    190628
!       fft_lib/cos_filt        : raised cosine high/low pass filter                                      McGirr    190628
!       fft_lib/various         : contains routines required to perform forward/backward FFT (FFTPACK5.1) McGirr    190628
!       acc_lib/acc_fft         : added subroutine which calls compute_fft                                McGirr    190628
!       acc_lib/acc_fft         : compute and write out both filtered and unfiltered power spectrums      McGirr    190710
!       fft_lib                 : added hann window function, window option and padding option            McGirr    190715
!       acc_lib/acc_fft         : added call to hann window and padding                                   McGirr    190716
!       gnv_read_data           : removed some debug                                                      Tregoning 190724
!       acc_lib                 : removed a bunch of unnecessary print statements                         Tregoning 190807
!       level1B_lib             : removed a bunch of unnecessary print statements                         Tregoning 190807
!       fft_lib/compute_fft     : simplified, moved mean and window calculation to calling program        McGirr    190815
!       fft_lib/window_function : extended version of hann_window, returns window given type and length   McGirr    190815
!       acc_lib/acc_cos_fft     : changed name, added call to window_function and mean calculation        McGirr    190815
!       read_soln_v3            : changed fatal to warning when VCV not found in .vcv file                Tregoning 190828
!       read_IC_list_v1,read_soln_v2 : removed allocate statements for prmnam                             Tregoning 190905
!       read_msc_hdr            : only output the first and last 4 header lines of a mascon file          Tregoning 190906
!       read_IC_list_v1         : include number of mascons when dimensioning "apriori", "soln", if required  Tregoning 190917
!       read_IC_list_v1         : fix satellite names for GRACE FO                                        Tregoning/McQueen 190917
!       read_IC_list_v1         : fixed bug in renaming satellite IC prmnam strings                       Tregoning 190918
!       msc_lib/calc_which_ternary : fixed bug mislocating ternaries by 1 cell in lat and lon             McQueen   191114
!       msc_lib/calc_which_ternary_ell : removed - not used                                               McQueen   191114
!       msc_lib/tern_lat_bands  : commented out - tern_lat_bands_ell should be used instead               McQueen   191114          
!       decompose_lib           : rename .f95 to f90; intel compiler are failing                          Allgeyer  191120
!       srf_2_los_quadprec      : replace quad definition kind=10 to kind=16; intel comp failing          Allgeyer  191120 
!       quat_lib_quad           : replace quad definition kind=10 to kind=16; intel comp failing          Allgeyer  191120 
!       findquat_quad           : replace quad definition kind=10 to kind=16; intel comp failing          Allgeyer  191120 
!       read_soln_addnorm       : changed format read of n_files to match changes to other programs       Tregoning 191231
!       matrices_lib            : added matmult, invert                                                   Tregoning 200117
!       read-soln_addnorm       : fixed bug when reading addnorm vcv with only one orbit soln+mascons     Tregoning 200128
!       write_mascon_record     : fixed formatting when writing our ternary region                        Tregoning 200129
!       read_msc_model          : new subroutine to read the output of get_rms_and_lowpass.py             Tregoning 200227
!       lininterp_r8            : commented out debug                                                     Tregoning 200304
!       level1b_lib/read_gnv_data : fixed bug in reading gnv1b data onto 5 second sampling                Tregoning 200624
!       SRF_2_LOS               : fixed bug in normalising position vectors (makes pitch calc wrong)      Tregoning 200724
!       msc_lib/calc_mascon_xyz : moved ocean ternary mascons CoM from ellipsoid to geoid                 Tregoning 201031
!
!  2.0  tools_lib               : new library to replace all the gamit/globk subroutines                  Tregoning 210621
!       fund_arguments_iers1992 : defined secs_in_360                                                     Tregoning 210806
!       allocate_mascon_arrays, deallocate_mascon_arrays : changed mcon_EWH to mcon_prim_EWH              Tregoning 210818
!       fund_arguments_iers1992 : fixed bug in the calculation of the mean elongation of the Moon from the Sun Tregoning 210916
!       deallocate_mascon_arrays: fixed bug in mascon EWH deallocation when input/output files different  Tregoning 220114
!       read_msc_model          : fixedbug when requested date is beyond last date in model               Tregoning 220117
!       some EMD subroutine     : changed "Status" to "STATUS"                                            Tregoning 220122
!       generate_inert_efixed,read_eop_ut1 : pass in calling_prog to use in status_update calls             Tregoning 220127
!       omtides_lib             : new subroutine containing "code" moved from mod_subs/oceantides.f90     Tregoning 220429
!
!  3.0  
!       all code                : all references to gamit routines removed                                Tregoning 220804
!       netcdf_lib              : subroutines for netcdf file handling: read_norm_netcdf, check_status_netcdf Tregoning 220824
!       read_norm_netcdf        : replace norm_netcdf_mod with soln_mod                                   Tregoning 221122
!       netcdf_lib              : turned off some debug                                                   Tregoning 230214

     END SUBROUTINE LVERSN

