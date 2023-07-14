# Makefile for orbits on Darwin 10.4.0 by unimake at Wed Jul 28 14:51:06 EST 2010
#

all: grace_lib

# ASSIGMENTS        
### GCC V4 gfortran flags ###
FTN = gfortran
#FFLAGS = -O3 -mcmodel=medium -Wuninitialized -fno-f2c -ffast-math -fno-automatic -fno-backslash  -fbounds-check    
FFLAGS = -O3 -mcmodel=medium -Wuninitialized -fno-f2c -ffast-math -fno-automatic -fno-backslash   `nf-config --fflags` -I../mod_subs                
CC = gcc
CFLAGS = -g

GRACE_LIB = grace_lib.a
SLIB = ~/gg/gamit/lib/gamit_lib.a ~/gg/libraries/comlib/com_lib.a ~/gg/libraries/matrix/kinv_lib.a ../gracefit/gracefit_lib.a 
X11LIB = /usr/X11R6/lib
II = /usr/include/X11

grace_lib :  ~/gg/libraries/comlib/com_lib.a $(GRACE_LIB) $(SLIB) ../includes/grace.param ../includes/gracefit.h90 ../includes/write_soln.h90 ../mod_subs/libmod.a

$(GRACE_LIB) : \
	$(GRACE_LIB)(amag3_quadprec.o) \
	$(GRACE_LIB)(cart_convert.o) \
	$(GRACE_LIB)(cart_to_sph.o) \
	$(GRACE_LIB)(calc_vector_component.o) \
	$(GRACE_LIB)(decompose_lib.o) \
	$(GRACE_LIB)(get_lambda.o) \
	$(GRACE_LIB)(gtorb_lib.o) \
	$(GRACE_LIB)(inert_lib.o) \
	$(GRACE_LIB)(invert.o) \
	$(GRACE_LIB)(infill_SCA.o) \
	$(GRACE_LIB)(jplephem.o) \
	$(GRACE_LIB)(legendre_matrix.o) \
	$(GRACE_LIB)(matmult.o) \
	$(GRACE_LIB)(matmult_quad.o) \
	$(GRACE_LIB)(matrices_lib.o) \
	$(GRACE_LIB)(msc_lib.o) \
	$(GRACE_LIB)(noise_robust_deriv.o) \
	$(GRACE_LIB)(printmat.o) \
        $(GRACE_LIB)(read_IC_list_v1.o)\
        $(GRACE_LIB)(read_soln_v1.o)\
        $(GRACE_LIB)(read_soln_v2.o)\
        $(GRACE_LIB)(read_soln_v3.o)\
        $(GRACE_LIB)(running_mean.o)\
        $(GRACE_LIB)(cross.o)\
	$(GRACE_LIB)(cross_quadprec.o) \
        $(GRACE_LIB)(dot.o)\
        $(GRACE_LIB)(dot_quad.o)\
        $(GRACE_LIB)(findquat.o)\
        $(GRACE_LIB)(findquat_quad.o)\
        $(GRACE_LIB)(gsec_to_ymdhms.o)\
        $(GRACE_LIB)(lininterp_r8.o)\
        $(GRACE_LIB)(PEPtime.o)\
	$(GRACE_LIB)(pnm_calc.o)\
        $(GRACE_LIB)(quat_lib.o)\
        $(GRACE_LIB)(quat_lib_quad.o)\
	$(GRACE_LIB)(read_sphharm.o)\
	$(GRACE_LIB)(rotmat2rpy.o)\
	$(GRACE_LIB)(rpy2rotmat.o)\
        $(GRACE_LIB)(smooth.o)\
        $(GRACE_LIB)(spotl_tide_lib.o)\
        $(GRACE_LIB)(SRF_2_LOS.o)\
        $(GRACE_LIB)(srf_2_los_quadprec.o)\
        $(GRACE_LIB)(tide_lib.o)\
        $(GRACE_LIB)(utils_constants.o)\
        $(GRACE_LIB)(utils_spline.o)\
	$(GRACE_LIB)(write_soln.o)\
	$(GRACE_LIB)(write_soln_v2.o)\
	$(GRACE_LIB)(write_soln_gracekal_v1.o)\
	$(GRACE_LIB)(which_primary_msc.o) \
	$(GRACE_LIB)(ymdhms_to_gracesec.o)  \
	$(GRACE_LIB)(ymdhms_to_mjd.o)  
	ranlib grace_lib.a
 
#$(GRACE_LIB)(write_soln.o): ../includes/write_soln.h90 ../includes/grace.param
#$(GRACE_LIB)(write_soln_v2.o): ../includes/write_soln.h90 ../includes/grace.param
#$(GRACE_LIB)(read_soln_v1.o): ../includes/write_soln.h90 ../includes/grace.param
#$(GRACE_LIB)(read_soln_v2.o): ../includes/write_soln.h90 ../includes/grace.param
#$(GRACE_LIB)(running_mean.o): ../includes/grace.param
#$(GRACE_LIB)(tide_lib.o): ../mod_subs/accel_mod.mod ../mod_subs/tide_mod.mod ../mod_subs/tide_netcdf_mod.mod
#$(GRACE_LIB)(msc_lib.o): ../mod_subs/mascon_mod.mod            
#$(GRACE_LIB)(tide_netcdf_mod.o): /usr/include/netcdf.mod            
#$(GRACE_LIB)(gtorb_lib.o): ../mod_subs/gtorb_mod.mod
.c.a:
	$(CC) -c $(CFLAGS)  -I$(II) $<
	ar rv -U $@ *.o
	rm -f $*.o

.f90.a:
	$(FTN) -c $(FFLAGS) -I../mod_subs/  $<
	ar rv -U $@ *.o
	rm -f $*.o

.f95.a:
	$(FTN) -c $(FFLAGS)  -I../mod_subs $<
	ar rv -U $@ *.o
	rm -f $*.o

.f.a:
	$(FTN) -c $(FFLAGS) -I../mod_subs/ $<
	ar rv -U $@ *.o
	rm -f $*.o

.SUFFIXES: .f90 .f95
