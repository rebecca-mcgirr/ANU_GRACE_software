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

    END SUBROUTINE LVERSN






















