# makefile for programs in graceorb
#
include ../makefile.variables

all: graceorb

# ASSIGNMENTS

#FFLAGS = -O3 -finit-local-zero -fno-automatic
# With OSX 10.6.7 and gcc-4.6.0 don't use -O this will lead to incorrect integration results.
#FFLAGS = -finit-local-zero -fno-automatic


ORBIT_LIB = orbit_lib.a
GRACE_LIB = ../lib/grace_lib.a

doy   : doy.f $(ORBIT_LIB)
	$(FTN) $(FFLAGS) doy.f $(ORBIT_LIB) $(SLIB) -o doy
	rm -f doy.o

graceorb : graceorb.f90 $(ORBIT_LIB) $(GRACE_LIB)  ../mod_subs/libmod.a
	$(FTN) $(FFLAGS) graceorb.f90 -I../mod_subs $(ORBIT_LIB) $(GRACE_LIB) ../mod_subs/libmod.a   $(NC_FFLAGS) $(SOFA_LIB) $(NC_LIB) $(H5_LIB) -o graceorb
	rm -f graceorb.o

#	$(ORBIT_LIB)(calpnm.o) \

#	$(ORBIT_LIB)(gauleg.o) \
#	$(ORBIT_LIB)(legcalc.o) \
# 	$(ORBIT_LIB)(matmult.o) \
# 	$(ORBIT_LIB)(noise_robust_deriv.o) \

$(ORBIT_LIB) : \
	$(ORBIT_LIB)(legcalc_norm.o) \
	$(ORBIT_LIB)(gauleg.o) \
	$(ORBIT_LIB)(calpnm.o) \
	$(ORBIT_LIB)(accelerom.o) \
	$(ORBIT_LIB)(adam_moulton_part.o) \
	$(ORBIT_LIB)(atm_tide_sph.o) \
	$(ORBIT_LIB)(oceanfield_sph.o) \
	$(ORBIT_LIB)(biasscale.o) \
	$(ORBIT_LIB)(bias_from_ACC.o) \
	$(ORBIT_LIB)(calc_corners.o) \
	$(ORBIT_LIB)(calc_perturbed_acc.o) \
	$(ORBIT_LIB)(cart_sph_convert.o) \
	$(ORBIT_LIB)(coeff_read.o) \
	$(ORBIT_LIB)(compute_nutation.o) \
	$(ORBIT_LIB)(efix_inert.o) \
	$(ORBIT_LIB)(etides_iers2010.o) \
	$(ORBIT_LIB)(evrtcf.o) \
	$(ORBIT_LIB)(generalrel.o) \
	$(ORBIT_LIB)(getgrid.o) \
	$(ORBIT_LIB)(gravcalc.o) \
	$(ORBIT_LIB)(ICpart_calc.o) \
	$(ORBIT_LIB)(input_lib.o) \
	$(ORBIT_LIB)(interpol.o) \
	$(ORBIT_LIB)(jd_to_ydsnew.o) \
	$(ORBIT_LIB)(jd_to_ymdhmsnew.o) \
	$(ORBIT_LIB)(hdwrt.o) \
	$(ORBIT_LIB)(hdwrt_binary.o) \
	$(ORBIT_LIB)(gmset.o) \
	$(ORBIT_LIB)(nutred_grace.o) \
	$(ORBIT_LIB)(orb_mascon_lib.o) \
	$(ORBIT_LIB)(orbitcalc.o) \
	$(ORBIT_LIB)(orbversn.o) \
	$(ORBIT_LIB)(otmgrd2arr.o) \
	$(ORBIT_LIB)(planetfield.o) \
	$(ORBIT_LIB)(pnrot_grace.o) \
	$(ORBIT_LIB)(polcoeff_read.o) \
	$(ORBIT_LIB)(polredek.o) \
	$(ORBIT_LIB)(poltid.o) \
	$(ORBIT_LIB)(quat2rotmat.o) \
	$(ORBIT_LIB)(quatcomb.o) \
	$(ORBIT_LIB)(realft.o) \
	$(ORBIT_LIB)(read_acc_sca_thr.o) \
	$(ORBIT_LIB)(read_GRACEinput.o) \
	$(ORBIT_LIB)(read_precnut.o) \
	$(ORBIT_LIB)(reconstruct.o) \
	$(ORBIT_LIB)(rotmat2quat.o) \
	$(ORBIT_LIB)(SCArot.o) \
	$(ORBIT_LIB)(shadow_grace.o) \
	$(ORBIT_LIB)(setup_ocean_tides.o) \
	$(ORBIT_LIB)(setup_ocean_tides_v2.o) \
	$(ORBIT_LIB)(solidfield.o) \
	$(ORBIT_LIB)(solidfield_IERS2010.o) \
	$(ORBIT_LIB)(sphharm.o) \
	$(ORBIT_LIB)(sphharmfield.o) \
	$(ORBIT_LIB)(tidegrid.o) \
	$(ORBIT_LIB)(trnsp.o) \
	$(ORBIT_LIB)(which_ICparam.o) \
	$(ORBIT_LIB)(which_MSCtide_param.o) \
	$(ORBIT_LIB)(ymd_to_doynew.o) \
	$(ORBIT_LIB)(ymdhms_to_jdnew.o) \
	$(ORBIT_LIB)(ymdhms_to_ydoyhms.o)
	ranlib orbit_lib.a

# $(ORBIT_LIB)(accelerom.o): ../mod_subs/accel_mod.mod ../mod_subs/accred_mod.mod ../mod_subs/bsscl_mod.mod
# $(ORBIT_LIB)(accred.o): ../mod_subs/accred_mod.mod ../mod_subs/inmod_mod.mod
# $(ORBIT_LIB)(adam_moulton_part.o): ../mod_subs/coeff_mod.mod ../mod_subs/spherhar_mod.mod ../mod_subs/gm_mod.mod ../mod_subs/rotation_mod.mod ../mod_subs/accel_mod.mod ../mod_subs/inmod_mod.mod
# $(ORBIT_LIB)(atm_tide_sph.o): ../mod_subs/gm_mod.mod ../mod_subs/gauleg_mod.mod ../mod_subs/lovenum_mod.mod ../mod_subs/coeff_mod.mod ../mod_subs/sat_mod.mod
# $(ORBIT_LIB)(biasscale.o): ../mod_subs/sat_mod.mod ../mod_subs/bsscl_mod.mod ../mod_subs/accel_mod.mod ../mod_subs/inmod_mod.mod
# $(ORBIT_LIB)(bias_from_ACC.o):  ../mod_subs/bsscl_mod.mod ../mod_subs/accel_mod.mod
# $(ORBIT_LIB)(calpnm.o): ../mod_subs/gm_mod.mod ../mod_subs/gauleg_mod.mod
# $(ORBIT_LIB)(calc_perturbed_acc.o): ../mod_subs/bsscl_mod.mod ../mod_subs/usno_mod.mod
# $(ORBIT_LIB)(coeff_read.o): ../mod_subs/coeff_mod.mod ../mod_subs/inmod_mod.mod ../mod_subs/spherhar_mod.mod
# $(ORBIT_LIB)(dealias.o): ../mod_subs/coeff_mod.mod ../mod_subs/inmod_mod.mod ../mod_subs/dealias_mod.mod
# $(ORBIT_LIB)(efix_inert.o): ../mod_subs/rotation_mod.mod
# $(ORBIT_LIB)(gauleg.o): ../mod_subs/gm_mod.mod
# $(ORBIT_LIB)(generalrel.o): ../mod_subs/gm_mod.mod
# $(ORBIT_LIB)(getgrid.o): ../mod_subs/gm_mod.mod ../mod_subs/gauleg_mod.mod
# $(ORBIT_LIB)(gmset.o): ../mod_subs/gm_mod.mod
# $(ORBIT_LIB)(gravcalc.o): ../mod_subs/dealias_mod.mod ../mod_subs/rotation_mod.mod ../mod_subs/spherhar_mod.mod ../mod_subs/coeff_mod.mod ../mod_subs/gm_mod.mod ../mod_subs/accel_mod.mod ../mod_subs/inmod_mod.mod ../mod_subs/sat_mod.mod
# $(ORBIT_LIB)(hdwrt.o): ../mod_subs/sat_mod.mod ../mod_subs/inmod_mod.mod ../mod_subs/bsscl_mod.mod ../mod_subs/accel_mod.mod
# $(ORBIT_LIB)(ICpart_calc.o): ../mod_subs/accel_mod.mod ../mod_subs/gm_mod.mod ../mod_subs/inmod_mod.mod ../mod_subs/coeff_mod.mod ../mod_subs/bsscl_mod.mod ../mod_subs/mascon_mod.mod
# $(ORBIT_LIB)(input_lib.o): ../mod_subs/accel_mod.mod ../mod_subs/mascon_mod.mod
# $(ORBIT_LIB)(inred.o): ../mod_subs/sat_mod.mod ../mod_subs/inmod_mod.mod
# $(ORBIT_LIB)(legcalc.o): ../mod_subs/spherhar_mod.mod
# $(ORBIT_LIB)(oceanfield_sph.o): ../mod_subs/gm_mod.mod ../mod_subs/gauleg_mod.mod ../mod_subs/lovenum_mod.mod ../mod_subs/coeff_mod.mod  ../mod_subs/sat_mod.mod
# $(ORBIT_LIB)(orb_mascon_lib.o): ../mod_subs/accel_mod.mod ../mod_subs/mascon_mod.mod
# $(ORBIT_LIB)(orbitcalc.o): ../mod_subs/usno_mod.mod ../mod_subs/sat_mod.mod ../mod_subs/rotation_mod.mod ../mod_subs/coeff_mod.mod ../mod_subs/spherhar_mod.mod ../mod_subs/gm_mod.mod ../mod_subs/evcoef_mod.mod ../mod_subs/timxtr_mod.mod ../mod_subs/gauleg_mod.mod ../mod_subs/lovenum_mod.mod ../mod_subs/inmod_mod.mod ../mod_subs/accred_mod.mod ../mod_subs/gpsant_mod.mod ../mod_subs/accel_mod.mod ../mod_subs/dealias_mod.mod ../mod_subs/bsscl_mod.mod
# $(ORBIT_LIB)(otmgrd2arr.o): ../mod_subs/gauleg_mod.mod
# $(ORBIT_LIB)(planetfield.o): ../mod_subs/coor_mod.mod ../mod_subs/gm_mod.mod ../mod_subs/tabsm_mod.mod
# $(ORBIT_LIB)(polcoeff.o): ../mod_subs/coeff_mod.mod ../mod_subs/inmod_mod.mod
# $(ORBIT_LIB)(polcoeff_read.o): ../mod_subs/coeff_mod.mod ../mod_subs/inmod_mod.mod
# $(ORBIT_LIB)(poltid.o): ../mod_subs/gm_mod.mod ../mod_subs/lovenum_mod.mod ../mod_subs/coeff_mod.mod
# $(ORBIT_LIB)(realft.o): ../mod_subs/gm_mod.mod
# $(ORBIT_LIB)(reconstruct.o): ../mod_subs/gm_mod.mod ../mod_subs/gauleg_mod.mod
# $(ORBIT_LIB)(SCArot.o): ../mod_subs/accred_mod.mod
# $(ORBIT_LIB)(setup_ocean_tides.o): ../mod_subs/tide_mod.mod
# $(ORBIT_LIB)(solidfield.o): ../mod_subs/coor_mod.mod ../mod_subs/rotation_mod.mod ../mod_subs/spherhar_mod.mod ../mod_subs/gm_mod.mod
# $(ORBIT_LIB)(sphharm.o): ../mod_subs/gauleg_mod.mod ../mod_subs/gm_mod.mod
# $(ORBIT_LIB)(sphharmfield.o): ../mod_subs/dealias_mod.mod ../mod_subs/gm_mod.mod ../mod_subs/spherhar_mod.mod ../mod_subs/coeff_mod.mod ../mod_subs/inmod_mod.mod
# $(ORBIT_LIB)(staticfield.o): ../mod_subs/spherhar_mod.mod ../mod_subs/coeff_mod.mod ../mod_subs/gm_mod.mod
# $(ORBIT_LIB)(tidegrid.o): ../mod_subs/gm_mod.mod ../mod_subs/gauleg_mod.mod ../mod_subs/tide_mod.mod
# $(ORBIT_LIB)(which_MSCtide_param.o): ../mod_subs/accel_mod.mod
# $(ORBIT_LIB)(write_soln_mod.o): ../includes/grace.param ../includes/write_soln.h90

.c.a:
	$(CC) $(CFLAGS) -I$(II) $<
	ar rv $@ *.o
	rm -f $*.o

.f90.a:
	$(FTN) $(FFLAGS) -I ../mod_subs/ -c   $<
	ar rv $@ *.o
	rm -f $*.o

.f.a:
	$(FTN) $(FFLAGS) -I ../mod_subs/ -c  $<
	ar rv $@ *.o
	rm -f $*.o

.SUFFIXES: .f90
