#!/bin/bash --norc

export LC_ALL="C"
export SHELL=$(type -p bash)

if [ $# -lt 3 ]
then
    echo "./batch_1node.sh filename startindx stopindx"
    exit;
fi


list_of_day=$1
declare -i start_indx=$2
declare -i stop_indx=$3
declare -i ndays=0
ndays=$(( stop_indx - start_indx + 1 ))

declare -a days_arr
for (( i=0; i<=$ndays - 1 ; i++ ))
do 
    days_arr[$i]=`sed -n $((start_indx + i))p $list_of_day`
done

yr=`echo $days_arr[1] |cut -c 1-4`
mm=`echo $days_arr[1] |cut -c 6-7`
echo "Year: " $yr


# To use an apriori mascon model put filename here
export msc_model_file=~/gg/grace/tables/apriori_batch7z_Tik2.5mm_adapt_casp20mm_2yr.txt
echo " apriori masccon model = " ${msc_model_file}

if [ ${yr} -lt 2018 ] 
then
    mission=0
    sat1="A"
    sat2="B"
    setup_cmdfile="GRACE.input"
    gracefit_cmdfile="gracefit.cmd.template_GRACE"
else
    mission=1
    sat1="C"
    sat2="D"
    setup_cmdfile="GRACE_FO.input"
    gracefit_cmdfile="gracefit.cmd.template_GRACEFO"
fi

export pos_unc=0.07
export vel_unc=0.07e-3
export kbra_unc=1.e-9
export kbrr_unc_iter3=10.e-6
export mission
export sat1 sat2
export setup_cmdfile
export gracefit_cmdfile
export reg=~/gg/grace/tables/mascons_stage5_V005_200km.reg  
export msc_file=~/gg/grace/tables/mascons_stage5_V005_200km
export grace_root=/g/data/tn40/GRACE/sfw/grace/test_201014rev2
export CMD=`echo "CMD: " $0 $*`
export host=`hostname`


# echo "---"
# module list 2>&1
# echo 'PATH= ' ${PATH}
# echo "---"

 source ~/bin/setpath_GRACE
# module load parallel/20191022
 echo "iter3_1node"
 echo " running " $ndays " days"
 module list 2>&1
 echo 'LD_LIBRARY_PATH= ' ${LD_LIBRARY_PATH}
 echo 'PATH= ' ${PATH}
 echo

printf "Code will run the following days on %s\n" "${host}"
printf '  %s   %s  %s\n' "${days_arr[@]}" 

# declare -f stage_init
#export pwd=/scratch/tn40/hm3063/grace/mrun #`pwd`

stage_init()
{
    #printf " --  Stage_init on day %s %(%Y%m%d-%H:%M:%S)T\n" "$1" "-1"
    _s=`date +%s`
    yr=`date --date ${1} +%Y`
    month=`date --date ${1} +%m`
    day=`date --date ${1} +%d`

    run_dir=msc${yr}-${month}-${day}
    work_dir=${PBS_JOBFS}/${run_dir}
    init_dir=${PBS_O_WORKDIR}/${yr}
    out_dir=${init_dir}/${run_dir}
    mkdir -p $work_dir
    mkdir -p $out_dir

    cd $work_dir
#    cd $out_dir
    cp -p ${PBS_O_WORKDIR}/iter3_1node.sh .
    echo $CMD > setup${yr}${month}${day}.log
    if [ -z "${msc_model_file}" ] 
    then
	${grace_root}/com/sh_setup_grace_200602 -t $yr $month $day 00 00 00 86280 -mission ${mission} -type h5 -acc_model therm fft none -mascon_file ${msc_file} -extend 21600 >> setup${yr}${month}${day}.log
    else
	${grace_root}/com/sh_setup_grace_200602 -t $yr $month $day 00 00 00 86280 -mission ${mission} -type h5 -acc_model therm fft none -mascon_file ${msc_file} -extend 21600 -msc_apriori_model ${msc_model_file} >> setup${yr}${month}${day}.log
    fi
    _e=`date +%s`
    printf " --  Stage_init on day %s %(%Y%m%d-%H:%M:%S)T  ... DONE in %i sec \n" "$1" "-1" "$((_e - _s))"
}


 stage3_orbs()
 {
#    printf " --  Stage_orb on day %s  %(%Y%m%d-%H:%M:%S)T\n" "$1" "-1"
    _s=`date +%s`
    yr=`date --date ${1} +%Y`
    month=`date --date ${1} +%m`
    day=`date --date ${1} +%d`
    run_dir=msc${yr}-${month}-${day}
    work_dir=${PBS_JOBFS}/${run_dir}
    init_dir=${PBS_O_WORKDIR}/${yr}
    out_dir=${init_dir}/${run_dir}
    cd ${work_dir}
    apr_IC=${yr}_IC.dat
    cp ${PBS_O_WORKDIR}/ICs/${apr_IC} .
    cp ${setup_cmdfile}  ${setup_cmdfile}.iter3

    sed -i 's/APRIORI_FILE      :.*/APRIORI_FILE      : '"${apr_IC}"'/' ${setup_cmdfile}.iter3
    sed -i 's/USE_APRIORI_ICS   : 0/USE_APRIORI_ICS   : 15/'  ${setup_cmdfile}.iter3
    sed -i 's/MASCON            : N/MASCON            : Y/' ${setup_cmdfile}.iter3
    if [ ! -z "${msc_model_file}" ] 
    then
	APR_MSC="msc_apriori_${yr}_${month}.vcv"
	#echo "monthly apriori mascon file: " $APR_MSC
    	sed -i 's/MSC_APRIORI_FILE.*/MSC_APRIORI_FILE  : '"${APR_MSC}"'/' ${setup_cmdfile}.iter3
        sed -i 's/USE_APRIORI_MCS.*/USE_APRIORI_MCS   : Y /' ${setup_cmdfile}.iter3
    fi
    OMP_NUM_THREADS=2
    ${grace_root}/graceorb/graceorb ${setup_cmdfile}.iter3 in_file_${sat1}_${yr}-${month}-${day}_00 > gorb${sat1}${yr}${month}${day}.log &
    ${grace_root}/graceorb/graceorb ${setup_cmdfile}.iter3 in_file_${sat2}_${yr}-${month}-${day}_00 > gorb${sat2}${yr}${month}${day}.log &
    wait
    _e=`date +%s`
    printf " --  Stage_orb on day %s  %(%Y%m%d-%H:%M:%S)T   ... DONE in %i sec\n" "$1" "-1" "$((_e - _s))"
}


stage3_fits()
{
    #printf " --  Stage_fits on day %s  %(%Y%m%d-%H:%M:%S)T  \n" "$1" "-1"
    _s=`date +%s`
    yr=`date --date ${1} +%Y`
    month=`date --date ${1} +%m`
    day=`date --date ${1} +%d`

    run_dir=msc${yr}-${month}-${day}
    work_dir=${PBS_JOBFS}/${run_dir}
    init_dir=${PBS_O_WORKDIR}/${yr}
    out_dir=${init_dir}/${run_dir}
    cd $work_dir

    cp ~/gg/grace/templates/${gracefit_cmdfile} gracefit_msc_iter3.cmd
    sed -i 's/shadow:   0/shadow:   0/'  gracefit_msc_iter3.cmd
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
    sed -i 's/mascons: .*/mascons: 1/' gracefit_msc_iter3.cmd
    sed -i 's/msc_len_const: .*/msc_len_const: 1/' gracefit_msc_iter3.cmd

    OMP_NUM_THREADS=12
    OPENBLAS_NUM_THREADS=12
    orb1=`ls GTORB_${yr}-${month}-${day}_*_${sat1}_*.h5`
    orb2=`ls GTORB_${yr}-${month}-${day}_*_${sat2}_*.h5`
    ${grace_root}/gracefit/gracefit gracefit_msc_iter3.cmd msc_${yr}${month}${day}.rms 4 $orb1 $orb2 $yr $month $day $reg 0 0 0 0 0 0 0 0 0 0 0 0 > gfit.log &
    wait

    du -sh .
    _e=`date +%s`
    ls -l > file_stats.log
    printf " --  Stage_fits on day %s  %(%Y%m%d-%H:%M:%S)T   ... DONE in %i sec \n" "$1" "-1" "$((_e - _s))"
    rm -f $orb1 $orb2
    gzip  msc_${yr}${month}${day}.norm    
}


stage3_transfer()
{
    #printf " --  Stage_transfer all days"
    _s=`date +%s`

    cd ${PBS_JOBFS}
    
    for DAY in `ls -d msc*`
    do
	yr=`echo ${DAY} |cut -c 4-7`
	run_dir=${DAY}
	work_dir=${PBS_JOBFS}/${run_dir}
	init_dir=${PBS_O_WORKDIR}/${yr}
	out_dir=${init_dir}/${run_dir}

	#echo ${yr} "   " ${DAY}
	#echo run_dir= $run_dir
	#echo work_dir= $work_dir
	#echo init_dir= $init_dir
	echo out_dir= $out_dir
	cd $out_dir
	
	if ls ${work_dir}/*.fatal 1> /dev/null 2>&1; then
	    rsync -ltr $work_dir/* .
	    cd ..
	    mv $run_dir ${run_dir}_failed
	else
	    rsync -t $work_dir/*.fit .
	    # rsync -t $work_dir/*.vcv .
	    rsync -t $work_dir/*.norm.gz .
	    rsync -t $work_dir/*.status .
	    rsync -t $work_dir/*.warning .
	    rsync -t $work_dir/plt* .
	    rsync -t $work_dir/GRACE*.input* .
	    rsync -t $work_dir/ICS* .
	    rsync -t $work_dir/in_file* .
	    rsync -t $work_dir/gracefit*.cmd* .
	    rsync -t $work_dir/*.log .
	    rsync -t ${PBS_O_WORKDIR}/iter3_1node.sh .
	    rsync -t ${grace_root}/com/sh_setup_grace_200602 .
	    #  rsync -t $work_dir/GT*iter3.h5 .
	fi
	rm -rf ${work_dir}
    done
    _e=`date +%s`
    printf " --  Stage_transfer on node %s  %(%Y%m%d-%H:%M:%S)T   ... DONE in %i sec \n" "${host}" "-1" "$((_e - _s))"
    
}

export -f stage_init
export -f stage3_orbs
export -f stage3_fits
export -f stage3_transfer

echo "STAGE INIT"
parallel --jobs 24 stage_init ::: ${days_arr[@]}
du -sh ${PBS_JOBFS}

echo "STAGE ORBS"
parallel --jobs 24 stage3_orbs ::: ${days_arr[@]}
du -sh ${PBS_JOBFS}

echo "STAGE FIT"
parallel --jobs 4 stage3_fits ::: ${days_arr[@]}
du -sh ${PBS_JOBFS}

echo "STAGE TRANSFER"
stage3_transfer

exit

