! mod_subs directory created in March 2017 to put all the mod files into the same directory
!
! ROTATION_MOD  : added efixed-inertial matrices to this mod file                                     Tregoning 170404
! MASCON_MOD    : add an extra "geoid" column to prim/sec/tern msc arrays                             Tregoning 170530
! MASCON_MOD    : increased nvar_tern to 9                                                            Tregoning 170605
! SPHERHARM_MOD : added atm_sph_first_call logical variable                                           Tregoning 170721
! TIDE_NETCDF_MOD: Read netcdf tides file                                                             Allgeyer long time ago
! GTORB_MOD     : Read Write orbits in HDF5 Format                                                    Allgeyer 170930 
! INMOD_MOD     :  added variables for linearising the accelerometer observations                     Tregoning 180619
! ACCRED_MOD    :  added sca_obs allocatable variable to store star camera quaternions                Tregoning 180619
! ACCRED_MOD    : added sca_step (the time interval between obs in SCA1B data)                        Tregoning 180821
! ACCRED_MOD    : add acc_obs allocatable variable to store accelerometer values for each satellite   Tregoning 180625
! INMOD_MOD     : added variable to describe on which surface the mascons are modelled                Tregoning 181023
! ETIDES_MOD    : new file containing the IERS 2010 solid body tide variables                         Tregoning 181029
! MASCON_MOD    : increased mascon file header lines to 200 (from 50)                                 Tregoning 181102
! ACCRED_MOD    : added array for constant accelerometer offset values for each sat (A/B/C/D)         Tregoning 190212
! GTORB_MOD     : increased the GTORB record length to permit  up to 48000 mascons                    McGirr    190305
! SOLN_MOD      : replaces write_soln_mod. Includes allocatable variables as well as maxparm          Tregoning 190312
! WRITE_SOLN_MOD: decommissioned ....                                                                 Tregoning 190312
! GTORB_mod     : moved declaration of GTORB_recl from gtorb_mod.f90 to graceorb/graceorb.f90         Tregoning 190403
! mascon_mod    : added extra variables for ocean mascons                                             Tregoning 190522
! dealias_mod   : increased to deg 180 (RL06 AOD1B is deg 180, up from 100) and 10 epochs (from 9)    Tregoning 190528
! usno_mod      : added eop_type array, being flag for Bull B (B), Bull A (A) or Bull A predicted (P) Tregoning 190701
! inmod_mod     : added msc_apriori_file variable                                                     Tregoning 190917
! soln_mod      : added apriori_msc array                                                             Tregoning 190927
! accred_mod    : added an array for a second star camera                                             Tregoning 191104
! gtorb_mod     : added gtorb_partials variable                                                       Tregoning 191126
! orbits_rw     : new subroutine to read the orbits partials directly in partials array               Allgeyer200401
! orbits_rw     : turned off debug                                                                    Tregoning 201007
! orbits_rw     : new subroutine to read the partials into a 2D part array                            Allgeyer201021
! etides_mod    : increased dimensions of legendre matrices to accommodate iers2010 solid body tide   Tregoning 210122
! orbit_rw_v2, orbits_rw : remove debug                                                               Tregoning 210610
! gm_mod        : include wgs84 ellipsoid parameter definitions                                       Tregoning 210612
! orbits_rw     : added argument as to whether to output mascon partials                              Tregoning 210802
! mascon_mod    : changed mcon_EWH to mcon_prim_EWH and added mcon_sec_EWH, mcon_tern_EWH             Tregoning 210818
! sig_process_mod : new mod_file for variables related to numerical diff + fft stuff                  Tregoning 220428
! omtides_mod   : removed "code" from oceantides.f90 to make omtides_mod with just variable defs      Tregoning 220429
!
! Version 2
! all code      : any reference to gamit routines removed
! Tregoning 220804
