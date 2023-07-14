  subroutine uversn

! started 6 January 2014 (even though files exist from much earlier than that !)
!
! GET_OCEAN_HEIGHT: add code to display the tidal consituents used   Tregoning 140106
! READ_GTORB_RECORD: added fatal stop if not enough command line arguments   Tregoning 140811
! CORREL_APR  :  modified to include inverting the covariance matrix for the msc constraints  Tregoning 140923
! MAKE_MSC_CONSTRAINT : new program to generate directly the msc constraint matrix to be added
!                       to the normal equations before inversion   Tregoning 140923
! CONFIGURE_MASCONS   : new program to generate a set of primary mascons that have no (or few) mixed land/ocean 
!                       ternary mascons within the same primary mascon  Tregoning 141104
! MASCON_lib  :  library of various tools for mascon readjustments  Tregoning 141104
! RESHAPE_PRIMARY_MASCONS : new program to reconfigure the distribution of ternary mascons within a set of 
!                           non-mixed primary mascons so that they all have roughly the same area and shape  Tregoning 141105
! VORONOI_MASCONS    : new program to reshape the primary mascons. Replaces reshape_primary_mascons that I never finished  Tregoning 141110
! CONFIGURE_MASCONS  : separate out the continental shelf ocean ternary mascons from the deeper mascons and put them all into the
!                      same primary mascon. voronoi_mascons can then break them up into different primary mascons   Tregoning 141113
! MASCON_TIDE_GRID   : pass lat/lon values to subroutine that writes out the tide grid header  Tregoning 150806
! READ_GTORB_RECORD  : increased character read of line from 100 to 200 characters  Tregoning 160404
! READ_GTORB_RECORD  : added capability to read lines of tidal amplitude information  Tregoning 160622
! TIDETOASC          : increased dimensions of the number of grid points   Tregoning 160725   
! secondary_from_primary : new routine to use the voronoi algorithm to break primary mascons into secondary. Also writes
!                          out the new mascon file format, where primary/secondary/ternary info is all contained in one file    Tregoning 160902
! secondary_from_primary : changed formatted write of line of primary macson info. Also, switched around the density/depth of ternary macons  Tregoning 160906
! tidecdftoasc   : new program to read the netcdf ocean tide file and output in ascii format
! make_daily_AOD1B  : new program to remove atm tides from the 6-hourly AOD1B data, then output a file (deg ord C S)  Tregoning 17601017
! range_accel  : new program to compute range accelerations from range rate residuals  Tregoning 161020
! WHICH_MASCON : updated to the new mascon file format   Tregoning 161201
! SEPARATE_MASCONS  : fixed bug in primary mascon codes   Tregoning 161205
! RESHAPE_MASCONS   : fixed bugs in prim-sec pointer and ternary mascon codes   Tregoning 161205
! RESHAPE_MASCONS   : modified to reshape a single, user-nominated mascon e.g. "MC0042"   Tregoning 170205
! COLOUR_MASCONS    : modified to also read a vcv file for mascon values        Tregoning 170308
!
! Version 2
! SECONDARY_MASCONS, SECONDARY_FROM_PRIMARY, MASCON_TIDE_GRID, TIDETOASC  :   removed obsolete programs    Tregoning 170316
! SPHHARM2XYZ       : fixed bugs found by new gfortran compiler    Tregoning 170316
! DIFF_LEGENDRE     : renamed spherhar_mod1 to spherhar_mod        Tregoning 170321
! REGULARIZE_MASCONS : limit sigma_polar for Antarctica to only land mascons  Tregoning 170517
! MODEL_ACC          : new program to remove a quadratic from level1B ACC files  Tregoning 170522
! ADDNORM2VCV        : new program to extract a single IC/mascon set from an addnorm VCV file   Tregoning 170530
! MASCON_transfer_ternary : added comments that "geoid" column is added to mascon arrays        Tregoning 170530
! MASCON_percent_land, MASCON_calc_CoM     : added "geoid" to the ternary array size            Tregoning 170531
! SEPARATE_MASCONS,  : updated for the addition of "geoid" to mascon files                      Tregoning 170531
! MASCON_percent_land, MASCON_transfer_ternary, MASCON_compute_CoM   : passed in nvar_tern variable  Tregoning 170531
! VORONOI_RESHAPE  : passed nvar_prim, nvar_tern variables to MASCON_calc_CoM                   Tregoning 170531
! VORONOI_RESHAPE  : update code to use nvar_tern, nvar_prim in various places                  Tregoning 170605
! RESHAPE_MASCONS  : updated to include "geoid" in mascon file                                  Tregoning 170606
! RESHAPE_MASCONS  : added rms_limit to command line arguments (after area)                     Tregoning 170608
! VORONOI_RESHAPE  : fixed (long-standing) bug causing NaN problems                             Tregoning 170608
! REGULARIZE_MASCONS : fixed bug in mcon_prim pointer (line 156)                                Tregoning 170609
! REGULARIZE_MASCONS : changed sigma_msc to be an allocatable variable                          Tregoning 170609
! READ_GTORB_DIFF    : GTORB record length defined in mod_subs/gtorb_mod                        Tregoning 170609
! COLOUR_MASCONS     : increased max mascon dimention to 8000                                   Tregoning 170611
! MODEL_ACC          : add catch for when there are no accelerometer observations               Tregoning 170615
! RESHAPE_MASCONS    : fixed name of mascon polygon file in output message. Added rms_limit
!                      to help of runstring                                                     Tregoning 170619
! REGULARIZE_MASCONS : added sigma for Greenland. Changed Hudson Bay to be treated as land.     Tregoning 170703
! REGULARIZE_MASCONS : fixed several bugs in program. Seems to work now ...                     Tregoning 170709
! MODEL_ACC          : output accelerometer obs for all epochs, not just when present           Tregoning 170724
! DIST_FLAG          : fixed bug (missing argument) in calling calc_mascon_xyz                  Tregoning 171123
! ADDNORM2VCV        : fixed indexing in reading IC parameter lines                             Tregoning 180226
! RESHAPE_MASCONS    : removed trimlen when calling read_mascon_file                            Tregoning 180307
! MODEL_ACC          : iterated the fitting of a quadratic to the ACC data                      Tregoning 180330
! REGULARIZE_MASCONS : added separate uncertainty (but not length) constraint for Alaska glaciers Tregoning 180427
! MODEL_ACC_V2       : allow a different model per accelerometer direction                      Tregoning 180705
! MODEL_ACC_V2       : added debug flag to call to acc_fit_quadratic                            Tregoning 180726
! separate_mascons   : changed character string passed to read_mascon_file                      Tregoning 180808
! RESHAPE_MASCON     : fixed error in output area of reshaped mascon header lines               Tregoning 180808
! EXTEND_L1B         : program to concatenate three L1B files together. Only ACC1B GRACE at present  Tregoning 180809
! EXTEND_L1B         : fixed naming of preceding and subsequent days                            Tregoning 180816
! MODEL_ACC_V2       : pass obs data length (in seconds) to acc_EMD1                            Tregoning 180816
! RESHAPE_MASCONS    : added additional ocean regions to code                                   Tregoning 180822
! REGULARISE_MASCONS : renamed from regularize_mascons                                          Tregoning 180831
! ADDNORM2VCV        : updated to read new and old VCV header formats                           Tregoning 180913
! MODEL_ACC_V2       : added acc_EMD2                                                           McGirr    180914
! MODEL_ACC_V2       : fixed passsing of obs data to acc_EMD1 and acc_EMD2                      McGirr    180924
! read_GTORB_record  : re-instated this program. Fixed call to GTORB_record_length              Tregoning 180926
! extract_ocean_mascons : changed length of filename passed to read_mascon_file                 Tregoning 181102
! extract_ocean_mascons : changed length of hashcode variable to be consistent with subroutine  Tregoning 181102
! MERGE_MASCONS      : increased max number of header records to 200 (from 50)                  Tregoning 181108
! DIST_FLAG          : excluded Caspian primary mascons from being nearest ocean mascon         Tregoning 181122
! RESHAPE_MASCONS    : fixed bug in numbering of secondary mascons                              Tregoning/McQueen 181123
! RESHAPE_MASCONS    : fixed this bug again                                                     Tregoning 181124
! MODEL_ACC_V2       : added acc_EMD3                                                           McGirr    181203
! DIST_FLAG          : fixed single precision pi, adding ".d0"                                  McQueen   181224
! WHICH_MASCON       : fixed single precision pi, adding ".d0"                                  McQueen   181224
! INIT_MASCONS       : fixed single precision pi, adding ".d0"                                  McQueen   181224
! CLASSIFY_MASCONS   : fixed single precision pi, adding ".d0"                                  McQueen   181224
! AUDIT_MASCONS      : added checks for close/small primaries and P-T parameter mismatch        McQueen   181231
! GSFC_to_apr_v2.4   : new program to convert GSFC mascons into values in ANU fit/vcv format    Tregoning 190111
! RESHAPE_MASCONS    : fixed array dimensioning for large numbers of mascons to be reshaped     McGirr    190204
! RESHAPE_MASCONS    : added check for empty reshaped primary mascons, produces warning         McGirr    190212
! RESHAPE_MASCONS    : added capability to reshape primary mascons to mimic ternary mascons     McGirr    190212
! RESHAPE_MASCONS    : fixed bug in indexing when writing out reshaped mascons                  McGirr/Tregoning 190226
! GSFC_to_apr_v2.4   : fixed bug in code identifying in which GSFC mascon an ANU mascon resides Tregoning 190314
! RESHAPE_MASCONS    : fixed bug to no longer try to write out reshaped primarys with 0 ternarys.
!                      Involved changing counters sent to write_mascon_record                   Tregoning 190319
! MODEL_ACC_V2       : fixed issues re passing of acc obs to acc_fit_quadratic                  McGirr    190320
! DIST_FLAG          : added coords to min_ternary messages. Added min_tern, min_sep to command line Tregoning 190321
! DIST_FLAG          : loop from mascon 1 so that info is output for all mascons                Tregoning 190406
! voronoi_reshape    : set max number of ternarys per reshaped primary. Improved OMP usage      Tregoning 190408
! voronoi_reshape    : increased the max permitted number of mascons for Antarctica to 110%     Tregoning 190410
! regularise_mascons : added separate sigmas for Arctic and Antarctic ocean mascons             Tregoning 190507
! regularise_mascons : added separate sigmas for GIA mascons					McGirr	  190521
! mascon_percent_land: sets percent land to -1 for GIA mascons based on density of ternarys     McGirr    190522
! voronoi_reshape    : sets density of GIA primary to 3300 when percent land is -1              McGirr    190522
! GIA_mascons        : creates 1 GIA prim containing all terns w/ lat <-60 from mascons_stage4  McGIrr    190522
! voronoi_mascons    : added logical flag to calling arguments (ocean mascon file or not?)      Tregoning 190523
! make_daily_AOD1B   : made RL06 compatible (3-hourly, no need to remove atm tides)             Tregoning 190526
! reshape_mascons    : added colour option for GIA ternary mascons                              McGirr    190527
! diff_ternarys      : fixed bug in allocation declarations for some solution arrays            Tregoning 190604
! regularise_mascons : updated help runstring to include sigma_GIA_mascons                      Tregoning 190607
! model_acc_v2       : added acc_fft option                                                     McGirr    190628
! model_acc_v2       : flo and fhi defined as 1/period and 1/period*1.5  (somewhat arbitrary)   McGirr    190703
! model_acc_v2       : applies hann window and padding to next power of 2                       McGirr    190716
! extend_L1B         : made it GRACE FO compatible, also for ACT data                           Tregoning 190807
! model_acc_v2       : lowered high freq cutoff sent to fft, made compatible with acc_extend    McGirr    190816
! addnorm2vcv        : changed 1st line of header information to make compatible with other programs Tregoning 190826
! regularise_mascons : change to dabs(0.3) the setting for reweighting sigmas based on .vcv     Tregoning 190826
! GSFC_to_apr_v2.4   : added option of extracting a time series of GSFC estimates on one mascon Tregoning 190830
! diff_ternarys      : resolved which piece of code allocates the "apriori" and "soln" arrays   Tregoning 190906
! reshape_mascons    : fixed bug in reading 5-digit mascon numbers to reshape                   Tregoning 190910
! regularise_mascons : improved the shape of the Greenland mask                                 Tregoning 190910
! calc_msc_correl    : fixed bug in reading 5-digit mascon numbers                              Tregoning 190912
! model_acc_v2       : changed filtering cutoff frequencies to 0.09e-3 and 0.11e-3              McGirr    190924
! apriori_mascons    : new program to create apriori vcv file with non-zero mascons             Tregoning 190927
! grab_mascons       : fixed bug in outputting just primary mascon numbers within a region      Tregoning 191017
! get_sca_data       : new routine to open, read, close a SCA1B file                            Tregoning 191104
! apriori_mascons    : changed calc_which_ternary_ell to calc_which_ternary                     McQueen   191114
! tidecdftoasc       : changed tern_lat_bands to tern_lat_bands_ell                             McQueen   191114
! diff_ternarys      : added fatal stop if soln2 file type not recognised                       Tregoning 191231
! grab_mascons       : fixed up defining region for SthnO and Antar                             Tregoning/McGirr 200114
! grab_mascons       : fix bugs in calculating the new total number of primary mascons          McGirr/Tregoning 200114
! apriori_mascons    : added total number of mascons to the second line of the file             Tregoning 200117
! diff_ternarys      : add capability to read VCV files created by apriori_mascons              Tregoning 200117
! dist_flag          : improved runstring help information                                      Tregoning 200128
! diff_ternarys      : fixed bug (in lib/read_soln_addnorm) when only 1 soln in addnorm vcv     Tregoning 200128
! separate_mascons   : set to land some ternarys in NW Africa that are below sea level          Tregoning 200129
! make_daily_AOD1B   : changed number of epochs to 145 (from 146) for AOD1B RL05 data           McQueen/Tregoning 200217
! atm_tide_lib       : fixed frequency error in time_arg calculation, changed constants to DP   McQueen   200223
! GSFC_to_apr_v2.4   : if output is vcv program writes to file instead of screen                McGirr    200226
! apriori_mascons    : modified to read a second type of input apriori mascon information       Tregoning 200227
! GSFC_to_apr_v2.4   : modified to conserve mass                                                McGirr    200302
! colour_mascons     : added extra columns to output file to match output of diff_ternarys      Tregoning 200304
! regularise_mascons : added Patagonia, Baffin as regions, moved Hudson Bay code lower          Tregoning 200310
! blur_ocean_boundaries: program to merge ocean mascons along a line of lon/lat                 McGirr    200422
! get_sca_data       : initialised "mission" to -1 so that the code will determine it from file Tregoning 200726
! regularise_mascons : use mcon_region descriptor to identify amazon and caspian mascons        Tregoning 200806
! msc_to_grid        : adapt to read fit files as well as vcv                                   Tregoning 200810
! colour_mascons     : output EWH values to 2 decimal places (rather than integers)             McQueen/Tregoning 200818
! regularise_mascons : add option to output only diagonal elements ("output_diag_vector")       Tregoning 201023
! regularise_mascons : control scaling of ocean mascons from command file (scale_ocean_sigmas)  Tregoning 201119
! regularise_mascons : control scaling of all mascons (using EWH values) from command file      Tregoning 201120
! xyz_to_llr         : new program to compute lat/lon/rad from GNV1B                            McGirr    201210
! bin_obs_to_msc     : new program to bin the number of overpasses per mascon                   McGirr    210208
! colour_mascons     : added capability to colour mascons according to sigmas given reg a file  McGirr    210706
! MASCON_nearest_msc : check that a nearby mascon actually has ternaries in it                  Tregoning 210902
! separate_mascons   : merge partial land and/or ocean primaries into neighbouring primaries
!                      if coastline has cut through more than 10% of the mixed primary          Tregoning 210902
! regularise_mascons : break W_Ant into W_Ant_coast and W_Ant_south                             Tregoning 211019
! regularise_mascons : permit name AMAZON or Amazon (as per V006 mascon file)                   Tregoning 211105
! dist_flag          : change back to writing out only lower-tri without diagonal elements      Tregoning 220106
! extend_L1B         : fixed bug in name of calling program                                     Tregoning 220114
! diff_ternarys      : changed "Status" to "STATUS"                                             Tregoning 220122
! apriori_mascons    : fixed but in setting ocean to zero when input/output mascons the same    Tregoning 220208
! reformat_TN13      : new program to reformat the degree-1 terms into format for Seb's python  Tregoning 220216
! colour_mascons     : increased print statement range to put space in output format            Tregoning 220301
! diff_mascons       : modify to read VCV files converted from hdf5 format                      Tregoning 220316
! model_acc_v2       : updated to version model_acc_v2_bec. Fixed bugs in status_update calls     Tregoning 220422
! colour_mascons     : add command line argument to decimate ternary output (default is every)  Tregoning 220608
! diff_ternarys      : add a scale factor to command line, to convert diffs into velocities     Tregoning 220621
! apriori_mascons    : fix bug when requested mascon field has few mascons than modelled field  Tregoning 220719
!
! Version 3.   :: remove dependencies on gamit code
! all code           : All reference to gamit routines removed                                  Tregoning 220804
! model_acc_v2       : removed duplicate declaration of bitmap                                  Tregoning 220804
! extend_L1B         : added ACH1B as a supported data type                                     Tregoning 221025
! apriori_mascons    : added flag to zero ocean mascons                                         McGirr    221101
! colour_mascons     : added density to ternary output line. Fixed formatting also              Tregoning 230103
! norm2netcdf        : new program, moved from test, to convert .norm to .nc                    Tregoning 230106

 return

 end