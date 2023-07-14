#####################    DEFAULT VALUES     #######################
# version of GRACE software to use

set version = "~/gt_nogamit"

# uncertainties
# 70 mm    uncertainty on all GNV1B positions
set pos_unc = 0.07        
# 70 um/s  uncertainty on all GNV1B velocities
set vel_unc = 0.07e-3     
# 0.3 um/s uncertainty on all KBR1B range rate observations
set kbrr_unc_iter0 = 1.e0
set kbrr_unc = 0.8e-6     
set kbrr_unc_iter3 = 10.e-6
# 1.0 nm/s^2 uncertainty on all KBR1B range acceleration observations
set kbra_unc_iter0 = 1.e0
set kbra_unc = 1.e-9     
# duration of orbit integration 2 minutes short of 24 hrs
set duration = 86280     
# duration of zeroth iteration
set iter0_duration = 21600
# duration of 1st iteration
set iter1_duration = 43200

# a priori mascon file (either a temporal gravity field or zero, or whatever you want)
set apr_mascon_file = ~/gg/grace/templates/mascon_zero.vcv

# a priori mascon model file (default ZAP = zero apriori mascons)
set msc_model_file = "ZAP"

# jump to iter3 gracefit ?
set process = all
set iter3_extn = ""

# mascon regularisation file
#set mascon_reg_file =  "mascons_stage4_350km_700km_0.5_0.1_0.5_0.1.reg"
#set mascon_reg_file = "mascons_0.1_0.025_0.25_0.25_0.25.reg"
#set mascon_reg_file = "mascons_stage4_V002_0.1_0.025_0.25.reg"
set mascon_reg_file = "mascons_stage5_V005_200km.reg"
#set mascon_reg_file = "mascons_stage5_V004.reg"
#set mascon_reg_file = "mascons_stage5_V005_aq.reg"

ln -s ~/gg/grace/tables/$mascon_reg_file .

# type of GTORB file
set gtorb_type = "h5"

# PT180622: which mission ? 0: GRACE, 1: GRACE FO, 2: GRACE II
set mission = 0
set RL_num = 02

# model to linearise the accelerometer data
set acc_model_AT = none
set acc_model_CT = fft
set acc_model_RD = none

##################### DECIPHER COMMAND LINE #######################
while ($#argv > 0 )
  set input = ( $argv )
  switch($input[1])
# get the day to be processed
    case -t:
      set yr     = $input[2] ; shift argv
      set month  = $input[3] ; shift argv
      set day    = $input[4] ; shift argv
    breaksw
# duration of orbit integration
     case -duration:
      set duration = $input[2]
    breaksw
# duration of zeroeth iteration orbit integration
     case -iter0_duration:
      set iter0_duration = $input[2]
    breaksw
# duration of 1st iteration orbit integration
     case -iter1_duration:
      set iter1_duration = $input[2]
    breaksw
# process ICs, mascons or all ?
     case -process:
      set process =  $input[2]
    breaksw
# iteration 3 gracefit extension on file names
     case -iter3_extn:
      set iter3_extn = $input[2]
    breaksw
# uncertainties
     case -pos_unc:
      set pos_unc = $input[2]
    breaksw
     case -vel_unc:
      set vel_unc = $input[2]
    breaksw
      case -kbrr_unc:
      set kbrr_unc = $input[2]
    breaksw
     case -kbra_unc:
      set kbra_unc = $input[2]
    breaksw
# version of software
     case -version:
      set version = $input[2]
    breaksw
# apriori mascon file
     case -apr_mascon_file:
      set apr_mascon_file = $input[2]
    breaksw
# apriori mascon model
     case -msc_apriori_model:
      set msc_model_file = $input[2]
    breaksw
#  mascon regularisation file
     case -mascon_reg_file:
      set mascon_reg_file = $input[2]
    breaksw
# type of GTORB file
     case -type:
      set gtorb_type = $input[2]
    breaksw
# mission
     case -mission:
      set mission = $input[2]
      if($mission == 0)set RL_num = 02
      if($mission == 1)set RL_num = 04
    breaksw
    case -acc_model:
      set acc_model_AT     = $input[2] ; shift argv
      set acc_model_CT     = $input[3] ; shift argv
      set acc_model_RD     = $input[4] ; shift argv
    breaksw
 endsw
  if ( $#argv > 0 ) shift argv
end
alldone:
######################################################################


######################################################################
# set the input file depending on mission
if ($mission == 0)then
  set setup_cmdfile = "GRACE.input"
  set gracefit_cmdfile = "gracefit.cmd.template_GRACE"
  set sat1 = "A"
  set sat2 = "B"


###############################
##                           ##
##  I T E R A T I O N    3   ##
##                           ##
###############################

${version}/com/sh_setup_grace_210913 -t $yr $month $day 00 00 00 $duration -mission $mission -type $gtorb_type -acc_model $acc_model_AT $acc_model_CT $acc_model_RD -mascon_file ~/gg/grace/tables/mascons_stage5_V005_200km -extend 21600 -mascon_surface topography -msc_apriori_model ~/gg/grace/tables/apriori_batch7z_Tik2.5mm_adapt_casp20mm_5yr.txt

echo
echo ## Begin Iteration 3 ##
iter3:

# set the required duration in the ICS file (default is 24hrs)           
##@ epochs = ( $duration / 5 )
##foreach FILE (ICS*)
##    sed -i "8s/.*/$epochs/" $FILE
##    end


set apr_IC = IC_${yr}_${month}_${day}_iter2.vcv

~/gt/util/apriori_mascons ~/gg/grace/tables/mascons_stage5_V005_200km /home/geod/lachland/gg/grace/tables/apriori_b7i1d_c4_2yr.txt 2016 08 01 ~/gg/grace/tables/mascons_stage5_V005_200km 3

cp ${setup_cmdfile}  ${setup_cmdfile}.iter3
sed -i 's/APRIORI_FILE      :.*/APRIORI_FILE      : '"${apr_IC}"'/' ${setup_cmdfile}.iter3
sed -i 's/USE_APRIORI_ICS   : 0/USE_APRIORI_ICS   : 15/'  ${setup_cmdfile}.iter3
sed -i 's/MASCON            : N/MASCON            : Y/' ${setup_cmdfile}.iter3

set msc_model_file = "ZAP"
if ($msc_model_file != "ZAP" ) then
  sed -i 's/MSC_APRIORI_FILE  : none/MSC_APRIORI_FILE  : msc_apriori_'"${yr}"'_'"${month}"'.vcv/' ${setup_cmdfile}.iter3
  sed -i 's/USE_APRIORI_MCS   : N/USE_APRIORI_MCS   : Y/' ${setup_cmdfile}.iter3
endif

# reintegrate the orbits
setenv OMP_NUM_THREADS 10
${version}/graceorb/graceorb ${setup_cmdfile}.iter3 in_file_${sat1}_${yr}-${month}-${day}_00 0 0 0 &
${version}/graceorb/graceorb ${setup_cmdfile}.iter3 in_file_${sat2}_${yr}-${month}-${day}_00 0 0 0 &
wait

# rename the orbits
mv GTORB_${yr}-${month}-${day}_00-00-00_${sat1}_${RL_num}.$gtorb_type  GTORB_${yr}-${month}-${day}_00-00-00_${sat1}_${RL_num}_iter3.$gtorb_type
mv GTORB_${yr}-${month}-${day}_00-00-00_${sat2}_${RL_num}.$gtorb_type  GTORB_${yr}-${month}-${day}_00-00-00_${sat2}_${RL_num}_iter3.$gtorb_type

# run gracefit, estimating pos/vel/scl/bias + tightly constrained mascons

iter3_gracefit:
# PT190725: switch on the cosine fft filter on the kbrr and kbra. Low- and high-frequency cutoffs defined by Bec McGirr from various tests.
cp ~/gg/grace/templates/${gracefit_cmdfile} gracefit_msc_iter3.cmd
sed -i 's/shadow:   0/shadow:   1/'  gracefit_msc_iter3.cmd
sed -i 's/ mission: M/ mission: '"$mission"'/' gracefit_msc_iter3.cmd
sed -i 's/PPP/'"$pos_unc"'/g' gracefit_msc_iter3.cmd
sed -i 's/VVV/'"$vel_unc"'/g' gracefit_msc_iter3.cmd
sed -i 's/KBRR_SIG/'"$kbrr_unc_iter3"'/' gracefit_msc_iter3.cmd
sed -i 's/KBRA_SIG/'"$kbra_unc"'/' gracefit_msc_iter3.cmd
sed -i 's/SCLX/1.0/' gracefit_msc_iter3.cmd
sed -i 's/SCLY/1.0/' gracefit_msc_iter3.cmd
sed -i 's/SCLZ/1.0/' gracefit_msc_iter3.cmd
sed -i 's/ kb_cos_fft: 0 / kb_cos_fft: 1 /' gracefit_msc_iter3.cmd
sed -i 's/kb_NR: 0 0/kb_NR: 0 1/' gracefit_msc_iter3.cmd
sed -i 's/mascons:            0/mascons:            1/' gracefit_msc_iter3.cmd
sed -i 's/msc_len_const:      0/msc_len_const:      1/' gracefit_msc_iter3.cmd

setenv OMP_NUM_THREADS 10
${version}/gracefit/gracefit gracefit_msc_iter3.cmd msc_${yr}_${month}_${day}_iter3.rms 4 GTORB_${yr}-${month}-${day}_00-00-00_${sat1}_${RL_num}_iter3.$gtorb_type GTORB_${yr}-${month}-${day}_00-00-00_${sat2}_${RL_num}_iter3.$gtorb_type $yr $month $day $mascon_reg_file 0 0 0 0 0 0 0 0 0 0 0 0

mv prefit.res prefit_iter3.res






#######################################################################
