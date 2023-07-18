include ../makefile.variables
MOD_LIB_NAME = libmod.a


all: $(MOD_LIB_NAME)


MOD_OBS = \
	accred_mod.o \
	mascon_mod.o \
	accel_mod.o \
	bsscl_mod.o \
	coeff_mod.o \
	coor_mod.o  \
	dealias_mod.o  \
	etides_mod.o  \
	evcoef_mod.o \
	gm_mod.o \
	orbits_rw.o \
	gtorb_mod.o \
	gauleg_mod.o \
	GPSant_mod.o \
	inmod_mod.o \
	lovenum_mod.o \
	nrtype_mod.o \
	omtides_mod.o \
	rpy_mod.o \
	rotation_mod.o \
	sat_mod.o \
	sig_process_mod.o \
	soln_mod.o \
	spherhar_mod.o \
	tabsm_mod.o \
	tide_mod.o \
	tide_netcdf_mod.o \
	timing_mod.o \
	timxtr_mod.o \
	utils_spline.o \
	utils_constants.o \
	physics_constants.o \
	physics_msis.o \
	usno_mod.o



#---- DEPENDENCIES
utils_spline.o:  utils_constants.o
orbits_rw.o: sat_mod.o inmod_mod.o
#-----------------------------------------------------------------------
#
#  TARGETS
#
#-----------------------------------------------------------------------
# Build (but do not install) the library.
all : $(MOD_LIB_NAME)

clean:
	rm *.o *.a *.mod

# Build the library.
$(MOD_LIB_NAME): $(MOD_OBS)
	ar rv $(MOD_LIB_NAME) $?

#-----------------------------------------------------------------------
%.o: %.f90
	$(FTN) $(FFLAGS) $(NC_FFLAGS) $(H5_FFLAGS) -c $<

#-----------------------------------------------------------------------
