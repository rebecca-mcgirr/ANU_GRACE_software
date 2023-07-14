# Makefile for GRACEFIT
include ../makefile.variables
all: gracefit addnorm gracesim 

gracefit : gracefit.f90 $(GRACEFIT_LIB) $(GRACE_LIB)  $(SOFA_LIB)  ../mod_subs/libmod.a gracefit_mod.mod
	$(FTN) $(FFLAGS) gracefit.f90 $(GRACEFIT_LIB) $(GRACE_LIB) $(SOFA_LIB)     -I../mod_subs ../mod_subs/libmod.a $(BLAS_LIB) -o gracefit $(H5_LIB) $(NC_LIB) 
	rm -f gracefit.o

gracesim : gracesim.f90 $(GRACEFIT_LIB) $(GRACE_LIB) $(SOFA_LIB)  ../mod_subs/libmod.a  gracefit_mod.mod
	$(FTN) $(FFLAGS) gracesim.f90 $(GRACEFIT_LIB) $(GRACE_LIB) $(SOFA_LIB)   -I../mod_subs ../mod_subs/libmod.a $(BLAS_LIB)  -o gracesim  $(H5_LIB) $(NC_LIB)
	rm -f gracesim.o

addnorm :  addnorm.f90 $(GRACEFIT_LIB) $(GRACE_LIB)  ../mod_subs/libmod.a gracefit_mod.mod
	$(FTN) $(FFLAGS)  addnorm.f90 $(GRACEFIT_LIB) $(GRACE_LIB)  -o  addnorm  -I../mod_subs ../mod_subs/libmod.a  $(BLAS_LIB) $(H5_LIB) $(NC_LIB) 
	rm -f addnorm.o

addnorm_netcdf :  addnorm_netcdf.f90 $(GRACEFIT_LIB) $(GRACE_LIB)  ../mod_subs/libmod.a gracefit_mod.mod
	$(FTN) $(FFLAGS)  addnorm_netcdf.f90 $(GRACEFIT_LIB) $(GRACE_LIB) $(NC_FFLAGS) $(NC_LIB)   -o  addnorm_netcdf  -I../mod_subs ../mod_subs/libmod.a  $(BLAS_LIB)  $(H5_LIB) $(NC_LIB)
	rm -f addnorm_netcdf.o

$(GRACEFIT_LIB) : \
	$(GRACEFIT_LIB)(addnorm_mod.o) \
	$(GRACEFIT_LIB)(add_CoM_cond.o) \
	$(GRACEFIT_LIB)(add_meanACC_cond.o) \
	$(GRACEFIT_LIB)(addnorm_lib.o) \
	$(GRACEFIT_LIB)(calc_beta_angle.o) \
	$(GRACEFIT_LIB)(command_lib.o) \
	$(GRACEFIT_LIB)(compute_aorc.o) \
	$(GRACEFIT_LIB)(debug_lib.o) \
	$(GRACEFIT_LIB)(header_lib.o) \
	$(GRACEFIT_LIB)(kb_lib.o) \
	$(GRACEFIT_LIB)(input_lib.o) \
	$(GRACEFIT_LIB)(gobs_lib.o) \
	$(GRACEFIT_LIB)(gversn.o) \
	$(GRACEFIT_LIB)(los_to_srf.o) \
	$(GRACEFIT_LIB)(LS_lib.o) \
	$(GRACEFIT_LIB)(mask_outliers.o) \
	$(GRACEFIT_LIB)(output_lib.o) \
	$(GRACEFIT_LIB)(plot_lib.o) \
	$(GRACEFIT_LIB)(rpy_AOC.o) \
	$(GRACEFIT_LIB)(shadow_lib.o) \
	$(GRACEFIT_LIB)(xyz2rac.o) \
	$(GRACEFIT_LIB)(output_netcdf.o)
	ranlib gracefit_lib.a

$(GRACEFIT_LIB)(add_CoM_cond.o):       gracefit_mod.mod
$(GRACEFIT_LIB)(add_meanACC_cond.o):   gracefit_mod.mod
$(GRACEFIT_LIB)(calc_beta_angle.o):    gracefit_mod.mod input.h90
$(GRACEFIT_LIB)(command_lib.o):        gracefit_mod.mod input.h90
$(GRACEFIT_LIB)(compute_aorc.o):       gracefit_mod.mod
$(GRACEFIT_LIB)(header_lib.o):         gracefit_mod.mod input.h90
$(GRACEFIT_LIB)(input_lib.o):          gracefit_mod.mod input.h90
$(GRACEFIT_LIB)(gobs_lib.o):           gracefit_mod.mod
$(GRACEFIT_LIB)(kb_lib.o):             gracefit_mod.mod
$(GRACEFIT_LIB)(LS_lib.o):             gracefit_mod.mod
$(GRACEFIT_LIB)(output_lib.o):         gracefit_mod.mod           output.h90
$(GRACEFIT_LIB)(plot_lib.o):           gracefit_mod.mod           output.h90
$(GRACEFIT_LIB)(rpy_AOC.o):            gracefit_mod.mod
$(GRACEFIT_LIB)(shadow_lib.o):         gracefit_mod.mod input.h90
.c.a:
	$(CC) -c $(CFLAGS) -I$(II) $<
	ar rv $@ *.o
	rm -f $*.o

.f90.a:
	$(FTN) -c $(FFLAGS) $(NC_FFLAGS) -I../mod_subs/  $<
	ar rv $@ *.o
	rm -f $*.o

.f.a:
	$(FTN) -c $(FFLAGS) $<
	ar rv $@ *.o
	rm -f $*.o

.f90.mod:
	$(FTN) -c  $< -I../mod_subs

.SUFFIXES: .f90
