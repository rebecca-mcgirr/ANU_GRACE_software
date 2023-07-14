include ../makefile.variables.nogamit
# makefile for programs in my bin directory
#
all: calc_range_accel repack_vcv  addnorm2vcv \
	model_acc model_acc_v2 calc_daily_mean_acc grab_mascons extract_ocean_mascons \
	recolour_mascons reshape_mascons reshape_ocean_mascons \
	separate_mascons  range_accel  make_daily_AOD1B   dist_flag \
	colour_mascons         \
	msc_to_grid            \
	which_mascon           \
	apriori_mascons          \
	mascon_threshold          \
	apriori_model_timeseries          \
	voronoi_mascons  \
	diff_ternarys    \
	GSFC_to_apr_v2.4     \
	covar_vs_dist       \
	atm_tide_model    \
	init_mascons classify_mascons extract_mascons \
	integrate_mascons writeout_mascons \
	extend_L1B             \
	norm2netcdf   \
	regularise_mascons  mascon_groundtrack  \
	rotate_mascons          \
	reformat_TN13    \
	get_beta_angle       \
	pitch_correction    \
	fix_GRACEFO_orient    \
	dist_deric    \
	extract_ocean_mascons_JP20190807 dist_flag_seb gia_over_oceans resolution_operator norm2asc \

#cluster_mascons
#timeseries_mascons
#tidetoasc
#tide_to_accel
# mascon_tide_sph
#msc_to_grid compare_inertial_frames_v2
# srp_geometry_test2 	add_srp_panel_acc tidegrid tidemascons
# reg_msc_by_lat convert_sphcoeff tidecdftoasc  read_GTORB_record sphharm2xyz sphharm2mascon normeq_solve LS_rpy read_AOD1B
# make_msc_constraint  make_msc_constraint_v2 smooth beta_angle  get_beta_angle  decimate test_svd square_thrust_ACC1B
# pass_over_msc mascon_sizes configure_mascons
# decompose_kbrr
# ASSIGNMENTS

#FFLAGS = -O3 -finit-local-zero -fno-automatic
# With OSX 10.6.7 and gcc-4.6.0 don't use -O this will lead to incorrect integration results.
#FFLAGS = -finit-local-zero -fno-automatic
#FFLAGS =  -O3 -fno-automatic -fbounds-check -fopenmp
#FFLAGS = -O3 -mcmodel=medium -Wuninitialized -fno-f2c -ffast-math                -fno-backslash  -fbounds-check -finit-local-zero -fopenmp  -I../mod_subs/ -L../mod_subs/ 
#NFFLAGS = -L/usr/lib -lnetcdff -Wl,-Bsymbolic-functions -Wl,-z,relro -lnetcdf -lnetcdf       
#FFLAGSPSD =  -O3
#FTN = gfortran
#FC = gfortran
#FTN = f90 
#CFLAGS += -g -c
#CC = cc
#LDFLAGS =  $(NC_LIB)
UTIL_LIB = util_lib.a
SPOTL_LIB = 
#LDFLAGS = -lnetcdf -lhdf5_hl -lhdf5 -lpthread -lsz -lz -ldl -lm -lcurl
#SLIB = ~/gg/gamit/lib/gamit_lib.a ~/gg/libraries/comlib/com_lib.a  ~/gg/libraries/matrix/kinv_lib.a  ../lib/grace_lib.a ../mod_subs/libmod.a 
#BLAS_LIB = ../BLAS/blas_LINUX.a
#COMLIB = ~/gg/libraries/comlib/com_lib.a
GRACEFIT_LIB = ../gracefit/gracefit_lib.a
GRACEORB_LIB = ../graceorb/orbit_lib.a
GRACEKAL_LIB = ../gracekal/gracekal_lib.a
SOFA_LIB = ../external/libsofa.a
GAMIT_LIB = ~/gg/gamit/lib/gamit_lib.a
#COM_LIB = ~/gg/libraries/comlib/com_lib.a
GRACE_LIB = ../mod_subs/libmod.a ../lib/grace_lib.a
#KINV_LIB = ~/gg/libraries/matrix/kinv_lib.a
#CFLAGS = -O3 -I$(FES_DIR) -I$(UTILC_DIR)
# tidegrid: tidegrid.o $(FES_DIR)/$(FES_LIB) $(UTILC_DIR)/$(UTILC_LIB)
# 	$(CC) $(CFLAGS) tidegrid.o $(FES_DIR)/$(FES_LIB) $(UTILC_DIR)/$(UTILC_LIB) -o tidegrid  $(LDFLAGS)

# tidemascons: tidemascons.o $(FES_DIR)/$(FES_LIB) $(UTILC_DIR)/$(UTILC_LIB)
# 	$(CC) tidemascons.o $(FES_DIR)/$(FES_LIB) $(UTILC_DIR)/$(UTILC_LIB) -o tidemascons  $(LDFLAGS)

# smooth : smooth.f90  
# 	$(FTN) $(FFLAGS) smooth.f90  -o smooth
# 	rm -f smooth.o

# test_svd : test_svd.f90  $(UTIL_LIB)  $(SLIB)
# 	$(FTN) $(FFLAGS) test_svd.f90 $(UTIL_LIB)  $(SLIB)  -o test_svd
# 	rm -f test_svd.o

# normeq_solve : normeq_solve.f90 $(GRACEFIT_LIB) $(UTIL_LIB)  $(SLIB)
# 	$(FTN) $(FFLAGS) normeq_solve.f90 $(GRACEFIT_LIB) $(SLIB) $(UTIL_LIB)  -o normeq_solve
# 	rm -f normeq_solve.o

# decimate : decimate.f90  
# 	$(FTN) $(FFLAGS) decimate.f90  -o smooth
# 	rm -f decimate.o

# beta_angle : beta_angle.f90  $(UTIL_LIB) $(SLIB)
# 	$(FTN) $(FFLAGS) beta_angle.f90 $(SLIB) $(UTIL_LIB) -o beta_angle
# 	rm -f beta_angle.o

#get_beta_angle : get_beta_angle.f90  $(UTIL_LIB) $(GRACE_LIB)
#	$(FTN) $(FFLAGS) get_beta_angle.f90 $(UTIL_LIB) $(GRACE_LIB)  -o get_beta_angle
#	rm -f get_beta_angle.o


# correl_apr : correl_apr.f90  $(UTIL_LIB) $(SLIB)
# 	$(FTN) $(FFLAGS) correl_apr.f90 $(SLIB) $(UTIL_LIB) -o correl_apr
# 	rm -f correl_apr.o

#make_msc_constraint : make_msc_constraint.f90  $(UTIL_LIB) $(SLIB) $(GRACEKAL_LIB) $(BLAS_LIB)
#	$(FTN) $(FFLAGS) make_msc_constraint.f90 $(SLIB) $(UTIL_LIB) $(GRACEKAL_LIB) $(BLAS_LIB)  -o make_msc_constraint
#	rm -f make_msc_constraint.o

#make_msc_constraint_v2 : make_msc_constraint_v2.f90  $(UTIL_LIB) $(SLIB) $(GRACEKAL_LIB) $(BLAS_LIB)
#	$(FTN) $(FFLAGS) make_msc_constraint_v2.f90 $(SLIB) $(UTIL_LIB) $(GRACEKAL_LIB) $(BLAS_LIB)  -o make_msc_constraint_v2
#	rm -f make_msc_constraint_v2.o

#make_msc_constraint_feb2016 : make_msc_constraint_feb2016.f90  $(UTIL_LIB) $(SLIB) $(GRACEKAL_LIB) $(BLAS_LIB)
#	$(FTN) $(FFLAGS) make_msc_constraint.f90 $(SLIB) $(UTIL_LIB) $(GRACEKAL_LIB) $(BLAS_LIB)  -o make_msc_constraint_feb2016
#	rm -f make_msc_constraint_feb2016.o

find_mascon : find_mascon.f90  $(UTIL_LIB)  $(GRACE_LIB)
	$(FTN) $(FFLAGS) find_mascon.f90  $(UTIL_LIB)  $(GRACE_LIB)-o find_mascon
	rm -f find_mascon.o

sphharm2mascon : sphharm2mascon.f90  $(UTIL_LIB) $(GRACE_LIB)   ../mod_subs ../mod_subs/libmod.a 
	$(FTN) $(FFLAGS) sphharm2mascon.f90 $(GRACE_LIB) $(UTIL_LIB)  -I../mod_subs ../mod_subs/libmod.a  -o sphharm2mascon
	rm -f sphharm2mascon.o

sphharm2xyz : sphharm2xyz.f90  $(UTIL_LIB) $(GRACE_LIB)
	$(FTN) $(FFLAGS) sphharm2xyz.f90 $(GRACE_LIB) $(UTIL_LIB) -o sphharm2xyz
	rm -f sphharm2xyz.o

# test_sphharm : test_sphharm.f90  $(UTIL_LIB) $(SLIB)
# 	$(FTN) $(FFLAGS) test_sphharm.f90 $(SLIB) $(UTIL_LIB) -o test_sphharm
# 	rm -f test_sphharm.o

mascon_tide_sph : mascon_tide_sph.f90 $(UTIL_LIB) $(GRACE_LIB) $(GRACEFIT_LIB) $(SPOTL_LIB)
	$(FTN) $(FFLAGS)  mascon_tide_sph.f90  -o mascon_tide_sph  $(GRACE_LIB) $(UTIL_LIB) $(SPOTL_LIB)
	rm -f  mascon_tide_sph.o

otm_latlon : otm_latlon.f90 $(UTIL_LIB) $(GRACE_LIB) $(GRACEFIT_LIB) $(SPOTL_LIB)
	$(FTN) $(FFLAGS)  otm_latlon.f90  -o otm_latlon $(GRACE_LIB) $(UTIL_LIB) $(SPOTL_LIB)
	rm -f  otm_latlon.o

square_thrust_ACC1B : square_thrust_ACC1B.f90  $(GRACE_LIB) $(GRACEFIT_LIB)
	$(FTN) $(FFLAGS)  square_thrust_ACC1B.f90  -o square_thrust_ACC1B $(GRACE_LIB) $(GRACEFIT_LIB)
	rm -f  square_thrust_ACC1B.o

ternary_to_secondary : ternary_to_secondary.f90  $(GRACEORB_LIB) $(GRACE_LIB) $(COM_LIB)
	$(FTN) $(FFLAGS) ternary_to_secondary.f90  -o ternary_to_secondary $(GRACEORB_LIB) $(GRACE_LIB) $(COM_LIB)
	rm -f ternary_to_secondary.o

#tidetoasc  : tidetoasc.f90  $(UTIL_LIB) $(SLIB)
#	$(FTN) $(FFLAGS) tidetoasc.f90 $(SLIB) $(UTIL_LIB) $(UTIL_LIB)  $(GRACE_LIB) $(COM_LIB) -lnetcdf -lnetcdff -o tidetoasc 
#	rm -f tidetoasc.o

tide_to_accel  : tide_to_accel.f90  $(UTIL_LIB) $(GRACE_LIB)
	$(FTN) $(FFLAGS) tide_to_accel.f90  $(GRACE_LIB) $(UTIL_LIB)  -lnetcdff -lnetcdf -o tide_to_accel 
	rm -f tide_to_accel.o

decompose_kbrr : decompose_kbrr.f95  $(UTIL_LIB) $(GRACE_LIB)
	$(FTN) $(FFLAGS) decompose_kbrr.f95  -o decompose_kbrr $(UTIL_LIB) $(GRACE_LIB)
	rm -f decompose_kbrr.o

LS_rpy : LS_rpy.f90  $(UTIL_LIB)  $(GRACE_LIB)
	$(FTN) $(FFLAGS) LS_rpy.f90 $(UTIL_LIB)  $(GRACE_LIB)  -o LS_rpy  
	rm -f LS_rpy.o

pitch_correction : pitch_correction.f90  $(UTIL_LIB)  $(GRACE_LIB)  ../mod_subs ../mod_subs/libmod.a $(SOFA_LIB)
	$(FTN) $(FFLAGS) pitch_correction.f90 $(UTIL_LIB)  $(GRACE_LIB) $(SOFA_LIB)  -I../mod_subs/ -mcmodel=medium   -o pitch_correction
	rm -f pitch_correction.o

fix_GRACEFO_orient : fix_GRACEFO_orient.f90  $(UTIL_LIB)  $(GRACE_LIB)  ../mod_subs ../mod_subs/libmod.a $(SOFA_LIB)
	$(FTN) $(FFLAGS) fix_GRACEFO_orient.f90 $(UTIL_LIB)  $(GRACE_LIB) $(SOFA_LIB)  -I../mod_subs/ -mcmodel=medium   -o fix_GRACEFO_orient
	rm -f fix_GRACEFO_orient.o

read_AOD1B : read_AOD1B.f90  $(UTIL_LIB)  $(GRACE_LIB)
	$(FTN) $(FFLAGS) read_AOD1B.f90 $(UTIL_LIB)  $(GRACE_LIB)  -o read_AOD1B
	rm -f read_AOD1B.o

reformat_TN13 : reformat_TN13.f90  $(UTIL_LIB)  $(GRACE_LIB)
	$(FTN) $(FFLAGS) reformat_TN13.f90 $(UTIL_LIB)  $(GRACE_LIB)  -o reformat_TN13
	rm -f reformat_TN13.o

read_GTORB_record : read_GTORB_record.f90  $(UTIL_LIB)   $(GRACE_LIB)
	$(FTN) $(FFLAGS) read_GTORB_record.f90 $(UTIL_LIB)  $(GRACE_LIB) -I../mod_subs -o read_GTORB_record
	rm -f read_GTORB_record.o

pass_over_msc : pass_over_msc.f90  $(UTIL_LIB)  $(GRACE_LIB)
	$(FTN) $(FFLAGS) pass_over_msc.f90 $(UTIL_LIB)  $(GRACE_LIB)  -o pass_over_msc
	rm -f pass_over_msc.o

mascon_sizes : mascon_sizes.f90  $(UTIL_LIB)  $(GRACE_LIB)
	$(FTN) $(FFLAGS) mascon_sizes.f90 $(UTIL_LIB)  $(GRACE_LIB)  -o mascon_sizes
	rm -f mascon_sizes.o

mascon_threshold : mascon_threshold.f90  $(UTIL_LIB)  $(GRACE_LIB)  ../mod_subs ../mod_subs/libmod.a 
	$(FTN) $(FFLAGS) mascon_threshold.f90 $(UTIL_LIB)  $(GRACE_LIB)  -I../mod_subs/ -mcmodel=medium  -o mascon_threshold
	rm -f mascon_threshold.o

mascon_threshold_series : mascon_threshold_series.f90  $(UTIL_LIB)  $(GRACE_LIB)  ../mod_subs ../mod_subs/libmod.a 
	$(FTN) $(FFLAGS) mascon_threshold_series.f90 $(UTIL_LIB)  $(GRACE_LIB)  -I../mod_subs/ -mcmodel=medium  -o mascon_threshold_series
	rm -f mascon_threshold_series.o

configure_mascons : configure_mascons.f90  $(UTIL_LIB)  $(GRACE_LIB)
	$(FTN) $(FFLAGS) configure_mascons.f90 $(UTIL_LIB)  $(GRACE_LIB) -I../mod_subs   -o configure_mascons
	rm -f configure_mascons.o

extract_mascons : extract_mascons.f90  $(UTIL_LIB)  $(GRACE_LIB)
	$(FTN) $(FFLAGS) extract_mascons.f90 $(UTIL_LIB)  $(GRACE_LIB) -I../mod_subs   -o extract_mascons
	rm -f extract_mascons.o

fix_ternary_file : fix_ternary_file.f90  $(UTIL_LIB)  $(GRACE_LIB)
	$(FTN) $(FFLAGS) fix_ternary_file.f90 $(UTIL_LIB)  $(GRACE_LIB)  -o fix_ternary_file
	rm -f fix_ternary_file.o

voronoi_mascons : voronoi_mascons.f90  $(UTIL_LIB)  $(GRACE_LIB)
	$(FTN) $(FFLAGS) voronoi_mascons.f90 $(UTIL_LIB)  $(GRACE_LIB)  -o voronoi_mascons
	rm -f voronoi_mascons.o

# srp_geometry : srp_geometry.f90  $(UTIL_LIB)  ../lib/libsofa.a $(GAMIT_LIB) $(COM_LIB) $(GRACE_LIB) $(GRACEORB_LIB) $(GRACEFIT_LIB)
# 	$(FTN) $(FFLAGS) srp_geometry.f90 $(UTIL_LIB) $(GRACEFIT_LIB) $(GRACE_LIB) $(GRACEORB_LIB) $(SOFA_LIB) $(GAMIT_LIB) $(COM_LIB) -o srp_geometry
# 	rm -f srp_geometry.o

# srp_geometry_test : srp_geometry_test.f90  $(UTIL_LIB)   $(SOFA_LIB)  $(GRACEORB_LIB) $(GAMIT_LIB) $(COM_LIB) $(GRACE_LIB) $(GRACEFIT_LIB)
# 	$(FTN) $(FFLAGS) srp_geometry_test.f90 $(UTIL_LIB)  $(GRACEORB_LIB) $(GRACEFIT_LIB) $(GRACE_LIB) $(SOFA_LIB) $(GAMIT_LIB) $(COM_LIB) $(GRACE_LIB) -o srp_geometry_test
# 	rm -f srp_geometry.o

# srp_geometry_test2 : srp_geometry_test2.f90  $(UTIL_LIB)   $(SOFA_LIB) $(GRACEORB_LIB) $(GAMIT_LIB) $(COM_LIB) $(GRACE_LIB) $(GRACEFIT_LIB)
# 	$(FTN) $(FFLAGS) srp_geometry_test2.f90 $(UTIL_LIB) $(GRACEORB_LIB)  $(GRACEFIT_LIB) $(GRACE_LIB) $(SOFA_LIB) $(GAMIT_LIB) $(COM_LIB) $(GRACE_LIB) -o srp_geometry_test2
# 	rm -f srp_geometry.o

# add_srp_panel_acc : add_srp_panel_acc.f90  $(UTIL_LIB)   $(SOFA_LIB)  $(GRACEORB_LIB) $(GAMIT_LIB) $(COM_LIB) $(GRACE_LIB) $(GRACEFIT_LIB)
# 	$(FTN) $(FFLAGS) add_srp_panel_acc.f90 $(UTIL_LIB)  $(GRACEORB_LIB) $(GRACEFIT_LIB) $(GRACE_LIB) $(SOFA_LIB) $(GAMIT_LIB) $(COM_LIB) $(GRACE_LIB) -o add_srp_panel_acc 
# 	rm -f add_srp_panel_acc.o

# sph_to_tide_grid : sph_to_tide_grid.f90  $(UTIL_LIB)  $(GAMIT_LIB) $(COM_LIB) $(GRACE_LIB) $(GRACEFIT_LIB) $(SPOTL_LIB)
# 	$(FTN) $(FFLAGS) sph_to_tide_grid.f90 $(UTIL_LIB)  $(GRACEFIT_LIB) $(GRACE_LIB) $(SOFA_LIB) $(GAMIT_LIB) $(COM_LIB) $(GRACE_LIB) $(GRACEORB_LIB)   $(SPOTL_LIB) -o sph_to_tide_grid 
# 	rm -f sph_to_tide_grid.o

which_mascon : which_mascon.f90  $(UTIL_LIB)   $(GRACE_LIB) $(GRACEFIT_LIB) 
	$(FTN) $(FFLAGS) which_mascon.f90 $(UTIL_LIB)  $(GRACE_LIB) -I../mod_subs -o which_mascon           
	rm -f which_mascon.o

bin_obs_to_msc : bin_obs_to_msc.f90  $(UTIL_LIB) $(GRACE_LIB) $(GRACEFIT_LIB) 
	$(FTN) $(FFLAGS) bin_obs_to_msc.f90 $(UTIL_LIB)  $(GRACE_LIB) -I../mod_subs -o bin_obs_to_msc           
	rm -f bin_obs_to_msc.o

dist_flag_seb : dist_flag_seb.f90  $(UTIL_LIB)  $(GRACE_LIB) $(GRACEFIT_LIB) 
	$(FTN) $(FFLAGS) dist_flag_seb.f90 $(UTIL_LIB)  $(GRACE_LIB)  -I../mod_subs -o dist_flag_seb           

dist_flag : dist_flag.f90  $(UTIL_LIB)   $(GRACE_LIB) $(GRACEFIT_LIB) 
	$(FTN) $(FFLAGS) dist_flag.f90 $(UTIL_LIB)  $(GRACE_LIB) -I../mod_subs -o dist_flag           
	rm -f dist_flag.o

dist_flag_bin : dist_flag_bin.f90  $(UTIL_LIB)  $(GRACE_LIB) $(GRACEFIT_LIB) 
	$(FTN) $(FFLAGS) dist_flag_bin.f90 $(UTIL_LIB)  $(GRACE_LIB)  -I../mod_subs -o dist_flag_bin           
	rm -f dist_flag_bin.o

tidecdftoasc : tidecdftoasc.f90  $(UTIL_LIB)  $(GRACE_LIB) $(GRACEFIT_LIB)    
	$(FTN) $(FFLAGS) tidecdftoasc.f90 $(UTIL_LIB)  $(GRACE_LIB)  -lnetcdf -lnetcdff -I../mod_subs -o tidecdftoasc       
	rm -f tidecdftoasc.o

make_daily_AOD1B : make_daily_AOD1B.f90  $(UTIL_LIB)  $(GRACE_LIB) $(GRACEFIT_LIB) 
	$(FTN) $(FFLAGS) make_daily_AOD1B.f90 $(UTIL_LIB)  $(GRACE_LIB)   -lnetcdff -lnetcdf  -o make_daily_AOD1B           
	rm -f make_daily_AOD1B.o

convert_sphcoeff : convert_sphcoeff.f90  $(UTIL_LIB)  $(GRACE_LIB) $(GRACEFIT_LIB) 
	$(FTN) $(FFLAGS) convert_sphcoeff.f90 $(UTIL_LIB)  $(GRACE_LIB)  -o convert_sphcoeff           
	rm -f convert_sphcoeff.o

range_accel : range_accel.f90  $(UTIL_LIB)  $(GRACE_LIB) $(GRACEORB_LIB)
	$(FTN) $(FFLAGS) range_accel.f90 $(UTIL_LIB)  $(GRACE_LIB)  $(GRACEORB_LIB) -o range_accel          
	rm -f range_accel.o

# reg_msc_by_lat   : reg_msc_by_lat.f90  $(UTIL_LIB)  $(GAMIT_LIB) $(COM_LIB) $(GRACE_LIB) $(GRACEFIT_LIB) 
# 	$(FTN) $(FFLAGS) reg_msc_by_lat.f90 $(UTIL_LIB)  $(GRACE_LIB) $(COM_LIB) -o reg_msc_by_lat           
# 	rm -f reg_msc_by_lat.o

separate_mascons   : separate_mascons.f90  $(UTIL_LIB)  $(GRACE_LIB) $(GRACEFIT_LIB) 
	$(FTN) $(FFLAGS) separate_mascons.f90 $(UTIL_LIB)  $(GRACE_LIB) -I../mod_subs -o separate_mascons           
	rm -f separate_mascons.o

reshape_mascons   : reshape_mascons.f90  $(UTIL_LIB)  $(GRACE_LIB) $(GRACEFIT_LIB) 
	$(FTN) $(FFLAGS) reshape_mascons.f90 $(UTIL_LIB)  $(GRACE_LIB)  -I../mod_subs -o reshape_mascons
	rm -f reshape_mascons.o

reshape_mascons_paul   : reshape_mascons_paul.f90  $(UTIL_LIB)  $(GRACE_LIB) $(GRACEFIT_LIB) 
	$(FTN) $(FFLAGS) -mcmodel=medium   reshape_mascons_paul.f90 $(UTIL_LIB)  $(GRACE_LIB)  -I../mod_subs -o reshape_mascons_paul
	rm -f reshape_mascons.o

reshape_mascons_bec   : reshape_mascons_bec.f90  $(UTIL_LIB)  $(GRACE_LIB) $(GRACEFIT_LIB) 
	$(FTN) $(FFLAGS) reshape_mascons_bec.f90 $(UTIL_LIB)  $(GRACE_LIB)  -I../mod_subs -o reshape_mascons_bec
	rm -f reshape_mascons_bec.o
reshape_mascons_julia   : reshape_mascons_julia.f90  $(UTIL_LIB)  $(GRACE_LIB) $(GRACEFIT_LIB)
	$(FTN) $(FFLAGS) reshape_mascons_julia.f90 $(UTIL_LIB)  $(GRACE_LIB)  -I../mod_subs -o reshape_mascons_julia
	rm -f reshape_mascons_julia.o

cluster_mascons   : cluster_mascons.f90  $(UTIL_LIB)  $(GRACE_LIB) $(GRACEFIT_LIB) 
	$(FTN) $(FFLAGS) cluster_mascons.f90 $(UTIL_LIB)  $(GRACE_LIB)  -I../mod_subs -o cluster_mascons
	rm -f cluster_mascons.o

reshape_ocean_mascons   : reshape_ocean_mascons.f90  $(UTIL_LIB)   $(GRACE_LIB) $(GRACEFIT_LIB) 
	$(FTN) $(FFLAGS) reshape_ocean_mascons.f90 $(UTIL_LIB)  $(GRACE_LIB) -I../mod_subs -o reshape_ocean_mascons
	rm -f reshape_ocean_mascons.o

GIA_mascons   : GIA_mascons.f90  $(UTIL_LIB)  $(GRACE_LIB) $(GRACEFIT_LIB) 
	$(FTN) $(FFLAGS) GIA_mascons.f90 $(UTIL_LIB)  $(GRACE_LIB)  -I../mod_subs -o GIA_mascons
	rm -f GIA_mascons.o

recolour_mascons   : recolour_mascons.f90  $(UTIL_LIB)  $(GRACE_LIB) $(GRACEFIT_LIB) 
	$(FTN) $(FFLAGS) recolour_mascons.f90 $(UTIL_LIB)  $(GRACE_LIB) -I../mod_subs -o recolour_mascons
	rm -f recolour_mascons.o

colour_mascons   : colour_mascons.f90  $(UTIL_LIB)   $(GRACE_LIB) $(GRACEFIT_LIB) 
	$(FTN) $(FFLAGS) colour_mascons.f90 $(UTIL_LIB)  $(GRACE_LIB)  -I../mod_subs -o colour_mascons
	rm -f colour_mascons.o

colour_GIA_mascons   : colour_GIA_mascons.f90  $(UTIL_LIB)   $(GRACE_LIB) $(GRACEFIT_LIB) 
	$(FTN) $(FFLAGS) colour_GIA_mascons.f90 $(UTIL_LIB)  $(GRACE_LIB)  -I../mod_subs -o colour_GIA_mascons
	rm -f colour_GIA_mascons.o

regularise_mascons   : regularise_mascons.f90  $(UTIL_LIB)   $(GRACE_LIB) $(GRACEFIT_LIB) 
	$(FTN) $(FFLAGS) regularise_mascons.f90 $(UTIL_LIB)  $(GRACE_LIB)  $(GRACEFIT_LIB)  -I../mod_subs ../mod_subs/libmod.a  -o regularise_mascons
	rm -f regularise_mascons.o

regularise_mascons_bec   : regularise_mascons_bec.f90  $(UTIL_LIB)  $(GRACE_LIB) $(GRACEFIT_LIB) 
	$(FTN) $(FFLAGS) regularise_mascons_bec.f90 $(UTIL_LIB)  $(GRACE_LIB)  $(GRACEFIT_LIB)  -I../mod_subs ../mod_subs/libmod.a  -o regularise_mascons_bec
	rm -f regularise_mascons_bec.o

regularise_mascons_paul   : regularise_mascons_paul.f90  $(UTIL_LIB)  $(GRACE_LIB) $(GRACEFIT_LIB) 
	$(FTN) $(FFLAGS) regularise_mascons_paul.f90 $(UTIL_LIB)  $(GRACE_LIB)  $(GRACEFIT_LIB)  -I../mod_subs ../mod_subs/libmod.a  -o regularise_mascons_paul
	rm -f regularise_mascons_paul.o

writeout_mascons   : writeout_mascons.f90  $(UTIL_LIB)   $(GRACE_LIB) $(GRACEFIT_LIB) 
	$(FTN) $(FFLAGS) writeout_mascons.f90 $(UTIL_LIB)  $(GRACE_LIB)  $(GRACEFIT_LIB)  -I../mod_subs -o writeout_mascons
	rm -f writeout_mascons.o

covar_vs_dist  : covar_vs_dist.f90  $(UTIL_LIB)   $(GRACE_LIB) $(GRACEFIT_LIB) 
	$(FTN) $(FFLAGS) covar_vs_dist.f90 $(UTIL_LIB)  $(GRACE_LIB)  $(GRACEFIT_LIB)  -I../mod_subs -o covar_vs_dist
	rm -f covar_vs_dist.o

extract_ocean_mascons   : extract_ocean_mascons.f90  $(UTIL_LIB)   $(GRACE_LIB) $(GRACEFIT_LIB) 
	$(FTN) $(FFLAGS) extract_ocean_mascons.f90 $(UTIL_LIB)  $(GRACE_LIB)  -I../mod_subs ../mod_subs/libmod.a   -o extract_ocean_mascons
	rm -f extract_ocean_mascons.o

extract_ocean_mascons_JP20190807: extract_ocean_mascons_JP20190807.f90  $(UTIL_LIB)   $(GRACE_LIB) $(GRACEFIT_LIB)
	$(FTN) $(FFLAGS) extract_ocean_mascons_JP20190807.f90 $(UTIL_LIB)  $(GRACE_LIB)   -I../mod_subs ../mod_subs/libmod.a   -o extract_ocean_mascons_JP20190807
	rm -f extract_ocean_mascons_JP20190807.o

grab_mascons   : grab_mascons.f90  $(UTIL_LIB)   $(GRACE_LIB) $(GRACEFIT_LIB) ../mod_subs/libmod.a
	$(FTN) $(FFLAGS) grab_mascons.f90 $(UTIL_LIB)  $(GRACE_LIB)   -I../mod_subs ../mod_subs/libmod.a -o grab_mascons
	rm -f grab_mascons.o

grab_ocean_mascons   : grab_ocean_mascons.f90  $(UTIL_LIB)   $(GRACE_LIB) $(GRACEFIT_LIB) ../mod_subs/libmod.a
	$(FTN) $(FFLAGS) grab_ocean_mascons.f90 $(UTIL_LIB)  $(GRACE_LIB)   -I../mod_subs ../mod_subs/libmod.a -o grab_ocean_mascons
	rm -f grab_ocean_mascons.o

blur_ocean_boundaries   : blur_ocean_boundaries.f90  $(UTIL_LIB)  $(GRACE_LIB) $(GRACEFIT_LIB) ../mod_subs/libmod.a
	$(FTN) $(FFLAGS) blur_ocean_boundaries.f90 $(UTIL_LIB)  $(GRACE_LIB)   -I../mod_subs ../mod_subs/libmod.a -o blur_ocean_boundaries
	rm -f blur_ocean_boundaries.o

merge_ocean_mascons   : merge_ocean_mascons.f90  $(UTIL_LIB)   $(GRACE_LIB) $(GRACEFIT_LIB) ../mod_subs/libmod.a
	$(FTN) $(FFLAGS) merge_ocean_mascons.f90 $(UTIL_LIB)  $(GRACE_LIB)   -I../mod_subs ../mod_subs/libmod.a -o merge_ocean_mascons
	rm -f merge_ocean_mascons.o

init_mascons : init_mascons.f90 
	$(FTN) -g -fno-automatic -O2 init_mascons.f90 $(NFFLAGS) -o init_mascons
	rm -f init_mascons.o

classify_mascons : classify_mascons.f90 interp_2p5m.o
	$(FTN) -g -fno-automatic -O2  classify_mascons.f90 interp_2p5m.o $(NC_FFLAGS) $(NC_LIB) -mcmodel=medium  -o classify_mascons
	rm -f classify_mascons.o

audit_mascons : audit_mascons.f90
	$(FTN) -g -fno-automatic -O2 -fbounds-check audit_mascons.f90 -mcmodel=medium -o audit_mascons
	rm -f audit_mascons.o

cleanup_mascons : cleanup_mascons.f90
	$(FTN) -g -fno-automatic -O2 -fbounds-check cleanup_mascons.f90 -mcmodel=medium -o cleanup_mascons
	rm -f cleanup_mascons.o

interp_2p5m.o : interp_2p5m.f
	$(FTN) -g -fno-automatic -O2 -c interp_2p5m.f

merge_mascons : merge_mascons.f90
	$(FTN) $(FFLAGS)  -fbounds-check  merge_mascons.f90 -o merge_mascons
	rm -f merge_mascons.o

merge_ternary : merge_ternary.f90
	$(FTN) $(FFLAGS) merge_ternary.f90 -o merge_ternary
	rm -f merge_ternary.o

repair_mascons : repair_mascons.f90
	$(FTN) $(FFLAGS) mascon_setup.o repair_mascons.f90 -o repair_mascons
	rm -f repair_mascons.o

mascon_setup: mascon_setup.f90
	$(FTN) -c mascon_setup.f90

diff_ternarys  : diff_ternarys.f90  $(UTIL_LIB)   $(GRACE_LIB) $(GRACEFIT_LIB) ../mod_subs/libmod.a
	$(FTN) $(FFLAGS) diff_ternarys.f90 $(UTIL_LIB)  $(GRACE_LIB)   -I../mod_subs ../mod_subs/libmod.a -o diff_ternarys
	rm -f diff_ternarys.o

GSFC_to_apr_v2.4  : GSFC_to_apr_v2.4.f90  $(UTIL_LIB)   $(GRACE_LIB) $(GRACEFIT_LIB) ../mod_subs/libmod.a
	$(FTN) $(FFLAGS) GSFC_to_apr_v2.4.f90 $(UTIL_LIB)  $(GRACE_LIB)   -I../mod_subs ../mod_subs/libmod.a -o GSFC_to_apr_v2.4
	rm -f GSFC_to_apr_v2.4.o

integrate_mascons : integrate_mascons.f90 $(UTIL_LIB)     ../mod_subs/libmod.a ../lib/matrices_lib.a
	$(FTN) $(FFLAGS)  integrate_mascons.f90 -o integrate_mascons -I../mod_subs ../mod_subs/libmod.a   $(UTIL_LIB)    $(GRACE_LIB) ../lib/matrices_lib.a
	rm -f integrate_mascons.o

timeseries_mascons : timeseries_mascons.f90
	$(FTN) $(FFLAGS)  timeseries_mascons.f90 -o timeseries_mascons
	rm -f timeseries_mascons.o

gia_over_oceans : gia_over_oceans.f90   $(GRACE_LIB) ../mod_subs/libmod.a
	$(FTN) $(FFLAGS)  $(GRACE_LIB) ../mod_subs/libmod.a  gia_over_oceans.f90 -o gia_over_oceans ../mod_subs/libmod.a
	rm -f gia_over_oceans.o

resolution_operator : resolution_operator.f90   $(GRACE_LIB) $(BLAS_LIB) ../mod_subs/libmod.a
	$(FTN) $(FFLAGS)  $(GRACE_LIB) ../mod_subs/libmod.a  resolution_operator.f90 -o resolution_operator ../mod_subs/libmod.a $(COM_LIB) $(GRACE_LIB) $(BLAS_LIB)
	rm -f resolution_operator.o

resolution_operator_bec : resolution_operator_bec.f90   $(GRACE_LIB) $(BLAS_LIB) ../mod_subs/libmod.a
	$(FTN) $(FFLAGS)  $(GRACE_LIB) ../mod_subs/libmod.a  resolution_operator_bec.f90 -o resolution_operator_bec ../mod_subs/libmod.a $(COM_LIB) $(GRACE_LIB) $(BLAS_LIB)
	rm -f resolution_operator_bec.o

svd_normeq  : svd_normeq.f90  $(UTIL_LIB)   $(GRACE_LIB) $(GRACEFIT_LIB) $(BLAS_LIB64) ../mod_subs/libmod.a
	$(FTN) $(FFLAGS) svd_normeq.f90 $(UTIL_LIB)  $(GRACEFIT_LIB)     $(GRACE_LIB)  $(BLAS_LIB64) -I../mod_subs ../mod_subs/libmod.a -o svd_normeq
	rm -f svd_normeq.o

svd_invert  : svd_invert.f90  $(UTIL_LIB)   $(GRACE_LIB) $(GRACEFIT_LIB) $(BLAS_LIB64) ../mod_subs/libmod.a
	$(FTN) $(FFLAGS) svd_invert.f90 $(UTIL_LIB)  $(GRACEFIT_LIB)     $(GRACE_LIB)  $(BLAS_LIB64) -I../mod_subs ../mod_subs/libmod.a -o svd_invert
	rm -f svd_invert.o

norm2asc  : norm2asc.f90  $(UTIL_LIB)   $(GRACE_LIB) $(GRACEFIT_LIB) ../mod_subs/libmod.a
	$(FTN) $(FFLAGS) norm2asc.f90 $(UTIL_LIB)  $(GRACEFIT_LIB)     $(GRACE_LIB)   -I../mod_subs ../mod_subs/libmod.a -o norm2asc
	rm -f norm2asc.o


# compare_inertial_frames_v2   : compare_inertial_frames_v2.f90  $(UTIL_LIB)  $(GAMIT_LIB) $(COM_LIB) $(GRACE_LIB) $(GRACEORB_LIB) $(SOFA_LIB) $(GAMIT_LIB) ../mod_subs/libmod.a
# 	$(FTN) $(FFLAGS) compare_inertial_frames_v2.f90 $(UTIL_LIB)  $(GRACE_LIB) $(COM_LIB) $(GRACEORB_LIB) $(SOFA_LIB) $(GAMIT_LIB) ../mod_subs/libmod.a -o compare_inertial_frames_v2
# 	rm -f compare_inertial_frames_v2.o
msc_to_grid : msc_to_grid.f90  $(UTIL_LIB)    $(GRACE_LIB)  
	$(FTN) $(FFLAGS) msc_to_grid.f90 $(UTIL_LIB)  $(GRACE_LIB)  -I../mod_subs -o msc_to_grid
	rm -f msc_to_grid.o

model_acc   : model_acc.f90  $(UTIL_LIB)    $(GRACE_LIB)  
	$(FTN) $(FFLAGS) model_acc.f90 $(UTIL_LIB)    $(GRACE_LIB)  -I../mod_subs ../mod_subs/libmod.a  -o model_acc
	rm -f model_acc.o

model_acc_v2     : model_acc_v2.f90  $(UTIL_LIB)   $(GRACE_LIB)  
	$(FTN) $(FFLAGS) model_acc_v2.f90 $(UTIL_LIB)  $(GRACE_LIB)  -I../mod_subs ../mod_subs/libmod.a  -o model_acc_v2    
	rm -f model_acc_v2.o

calc_daily_mean_acc  : calc_daily_mean_acc.f90  $(UTIL_LIB)    $(GRACE_LIB) 
	$(FTN) $(FFLAGS) calc_daily_mean_acc.f90 $(UTIL_LIB)  $(GRACE_LIB)  -I../mod_subs ../mod_subs/libmod.a -o calc_daily_mean_acc 
	rm -f calc_daily_mean_acc.o

model_acc_v2_bec     : model_acc_v2_bec.f90  $(UTIL_LIB)    $(GRACE_LIB)  
	$(FTN) $(FFLAGS) model_acc_v2_bec.f90 $(UTIL_LIB)  $(GRACE_LIB) -o model_acc_v2_bec   
	rm -f model_acc_v2_bec.o

model_acc_v2_paul : model_acc_v2_paul.f90  $(UTIL_LIB)    $(GRACE_LIB) 
	$(FTN) $(FFLAGS) model_acc_v2_paul.f90 $(UTIL_LIB)  $(GRACE_LIB)  -I../mod_subs ../mod_subs/libmod.a  -o model_acc_v2_paul    
	rm -f model_acc_v2.o

rfft_acc_test : rfft_acc_test.f90  $(UTIL_LIB)    $(GRACE_LIB)  
	$(FTN) $(FFLAGS) rfft_acc_test.f90 $(UTIL_LIB)  $(GRACE_LIB) -o rfft_acc_test   
	rm -f rfft_acc_test.o

addnorm2vcv   : addnorm2vcv.f90  $(UTIL_LIB)    $(GRACE_LIB)  
	$(FTN) $(FFLAGS) addnorm2vcv.f90 $(UTIL_LIB)  $(GRACE_LIB)  -I../mod_subs ../mod_subs/libmod.a  -o addnorm2vcv
	rm -f addnorm2vcv.o

mascon_groundtrack   : mascon_groundtrack.f90  $(UTIL_LIB)    $(GRACE_LIB)  ../mod_subs
	$(FTN) $(FFLAGS) mascon_groundtrack.f90 $(UTIL_LIB)  $(GRACE_LIB)  -I../mod_subs ../mod_subs/libmod.a  -o mascon_groundtrack
	rm -f mascon_groundtrack.o

repack_vcv   : repack_vcv.f90  $(UTIL_LIB)    $(GRACE_LIB)  
	$(FTN) $(FFLAGS) repack_vcv.f90 $(UTIL_LIB)  $(GRACE_LIB) -o repack_vcv
	rm -f repack_vcv.o

calc_range_accel   : calc_range_accel.f90  $(UTIL_LIB)    $(GRACE_LIB)  
	$(FTN) $(FFLAGS) calc_range_accel.f90 $(UTIL_LIB)  $(GRACE_LIB) -I../mod_subs ../mod_subs/libmod.a  -o calc_range_accel
	rm -f calc_range_accel.o

extend_L1B   : extend_L1B.f90  $(UTIL_LIB)    $(GRACE_LIB)  
	$(FTN) $(FFLAGS) extend_L1B.f90 $(UTIL_LIB)  $(GRACE_LIB) -o extend_L1B
	rm -f extend_L1B.o

calc_msc_correl : calc_msc_correl.f90 $(UTIL_LIB)   $(GRACE_LIB)  
	$(FTN) $(FFLAGS) calc_msc_correl.f90 $(UTIL_LIB)  $(GRACE_LIB) -o calc_msc_correl
	rm -f calc_msc_correl.o

apriori_mascons  : apriori_mascons.f90  $(UTIL_LIB)   $(GRACE_LIB) $(GRACEFIT_LIB) ../mod_subs/libmod.a
	$(FTN) $(FFLAGS) apriori_mascons.f90 $(UTIL_LIB)  $(GRACE_LIB)   -I../mod_subs ../mod_subs/libmod.a -o apriori_mascons
	rm -f apriori_mascons.o

apriori_mascons_paul  : apriori_mascons.f90  $(UTIL_LIB)  $(GRACE_LIB) $(GRACEFIT_LIB) ../mod_subs/libmod.a
	$(FTN) $(FFLAGS) apriori_mascons.f90 $(UTIL_LIB)  $(GRACE_LIB)   -I../mod_subs ../mod_subs/libmod.a -o apriori_mascons_paul
	rm -f apriori_mascons.o

apriori_model_timeseries  : apriori_model_timeseries.f90  $(UTIL_LIB)   $(GRACE_LIB) $(GRACEFIT_LIB) ../mod_subs/libmod.a
	$(FTN) $(FFLAGS) apriori_model_timeseries.f90 $(UTIL_LIB)  $(GRACE_LIB)   -I../mod_subs ../mod_subs/libmod.a -o apriori_model_timeseries
	rm -f apriori_mascons.o

bec_rotate_mascons  : bec_rotate_mascons.f90  $(UTIL_LIB)   $(GRACE_LIB) $(GRACEFIT_LIB) ../mod_subs/libmod.a
	$(FTN) $(FFLAGS) bec_rotate_mascons.f90 $(UTIL_LIB)  $(GRACE_LIB)  -I../mod_subs ../mod_subs/libmod.a -o bec_rotate_mascons
	rm -f bec_rotate_mascons.o

rotate_mascons  : rotate_mascons.f90  $(UTIL_LIB)   $(GRACE_LIB) $(GRACEFIT_LIB) ../mod_subs/libmod.a
	$(FTN) $(FFLAGS) rotate_mascons.f90 $(UTIL_LIB)  $(GRACE_LIB)   -I../mod_subs ../mod_subs/libmod.a -o rotate_mascons
	rm -f rotate_mascons.o

rotate_mascons_paul  : rotate_mascons_paul.f90  $(UTIL_LIB)   $(GRACE_LIB) $(GRACEFIT_LIB) ../mod_subs/libmod.a
	$(FTN) $(FFLAGS) rotate_mascons_paul.f90 $(UTIL_LIB)  $(GRACE_LIB)   -I../mod_subs ../mod_subs/libmod.a -o rotate_mascons_paul
	rm -f rotate_mascons_paul.o

atm_tide_model  : atm_tide_model.f90  $(UTIL_LIB)   $(GRACE_LIB) $(GRACEFIT_LIB) ../mod_subs/libmod.a
	$(FTN) $(FFLAGS) atm_tide_model.f90 $(UTIL_LIB)  $(GRACE_LIB)   -I../mod_subs ../mod_subs/libmod.a -o atm_tide_model
	rm -f atm_tide_model.o

xyz_to_llr : xyz_to_llr.f90  $(UTIL_LIB)  $(GRACE_LIB)  
	$(FTN) $(FFLAGS) xyz_to_llr.f90 $(UTIL_LIB)  $(GRACE_LIB)   -o xyz_to_llr
	rm -f xyz_to_llr.o

coastline_mascons : coastline_mascons.f90  $(UTIL_LIB)  $(GRACE_LIB) ../mod_subs/libmod.a
	        $(FTN) $(FFLAGS) coastline_mascons.f90 $(UTIL_LIB)  $(GRACE_LIB)   -I../mod_subs ../mod_subs/libmod.a   -o coastline_mascons
		        rm -f coastline_mascons.o

dist_deriv : dist_deriv.f95 $(UTIL_LIB) $(GRACE_LIB) $(GRACEFIT_LIB)
	$(FTN) $(FFLAGS)  dist_deriv.f95  -o dist_deriv  $(GRACE_LIB) $(UTIL_LIB) $(GRACEFIT_LIB) -I../mod_subs ../mod_subs/libmod.a
	rm -f  dist_deriv.o

transplant  : transplant.f90  $(UTIL_LIB)   $(GRACE_LIB) $(GRACEFIT_LIB) ../mod_subs/libmod.a
	$(FTN) $(FFLAGS) transplant.f90 $(UTIL_LIB)  $(GRACEFIT_LIB) $(GRACE_LIB)  -I../mod_subs ../mod_subs/libmod.a -o transplant
	rm -f rotate_mascons_paul.o

altgrid_to_primary  : altgrid_to_primary.f90  $(UTIL_LIB)   $(GRACE_LIB) $(GRACEFIT_LIB) ../mod_subs/libmod.a
	$(FTN) $(FFLAGS) altgrid_to_primary.f90 $(UTIL_LIB)  $(GRACEFIT_LIB) $(GRACE_LIB)  -I../mod_subs ../mod_subs/libmod.a -o altgrid_to_primary 
	rm -f altgrid_to_primary.o

norm2netcdf : norm2netcdf.f90  $(GRACE_LIB) ../mod_subs/libmod.a
	$(FTN)  $(FFLAGS)  norm2netcdf.f90  $(NC_FFLAGS) $(NC_LIB) $(GRACE_LIB)  -I../mod_subs ../mod_subs/libmod.a  -mcmodel=medium  -o norm2cdf  $(NC_FFLAGS) $(NC_LIB)
	rm -f norm2netcdf.o

read_netcdf : read_netcdf.f90  $(GRACE_LIB) ../mod_subs/libmod.a
	$(FTN) -g -fno-automatic -O2  read_netcdf.f90  $(NC_FFLAGS) $(NC_LIB) $(GRACE_LIB)  -I../mod_subs ../mod_subs/libmod.a  -mcmodel=medium  -o read_netcdf  $(NC_FFLAGS) $(NC_LIB)
	rm -f read_netcdf.o


$(UTIL_LIB) : \
	$(UTIL_LIB)(analyse_lib.o) \
	$(UTIL_LIB)(atm_tide_lib.o) \
	$(UTIL_LIB)(calc_beta_angle.o) \
	$(UTIL_LIB)(calc_yrate.o) \
	$(UTIL_LIB)(gauleg.o) \
	$(UTIL_LIB)(get_ocean_height.o) \
	$(UTIL_LIB)(get_sca_data.o) \
	$(UTIL_LIB)(mascon_lib.o) \
	$(UTIL_LIB)(mean_bias.o) \
	$(UTIL_LIB)(nr.o) \
	$(UTIL_LIB)(nrutil.o) \
	$(UTIL_LIB)(pythag.o) \
	$(UTIL_LIB)(read_ACC1B.o) \
	$(UTIL_LIB)(read_level1B_files.o) \
	$(UTIL_LIB)(read_THR1B.o) \
	$(UTIL_LIB)(sort_longitude.o) \
	$(UTIL_LIB)(svbksb.o) \
	$(UTIL_LIB)(svdvar.o) \
	$(UTIL_LIB)(svdcmp.o) 
	ranlib util_lib.a

#$(UTIL_LIB)(mascon_lib.o): ../mod_subs/mascon_mod.mod
#$(UTIL_LIB)(make_msc_constraint_v2.o): ../mod_subs/mascon_mod.mod
#$(UTIL_LIB)(separate_mascons.o): ../mod_subs/mascon_mod.mod
#$(UTIL_LIB)(reshape_mascons.o): ../mod_subs/mascon_mod.mod
#$(UTIL_LIB)(regularise_mascons.o): ../mod_subs/mascon_mod.mod
#$(UTIL_LIB)(extract_ocean_mascons.o): ../mod_subs/mascon_mod.mod
#$(UTIL_LIB)(extract_mascons.o): ../mod_subs/mascon_mod.mod

.c.a:
	$(CC) $(CFLAGS)  -I$(II) $<
	ar rv -U $@ *.o
	rm -f $*.o

#.c:
#	$(CC) $(CFLAGS) $< $(LDFLAGS)

#.o:
#	$(FTN) $(FFLAGS) $< $(FLDFLAGS)
.f90.a:
	$(FTN) $(FFLAGS) -I../mod_subs/  -c $<
	ar rv -U $@ *.o
	rm -f $*.o

.f95.a:
	$(FTN) $(FFLAGS)  -c $<
	ar rv -U $@ *.o
	rm -f $*.o

.f.a:
	$(FTN) $(FFLAGS) -c $<
	ar rv -U $@ *.o
	rm -f $*.o


.SUFFIXES: .f90 .f95
