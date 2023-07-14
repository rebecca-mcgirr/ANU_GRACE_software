#!/usr/bin/env python3
import numpy as np
import datetime as dt
import os
from fnmatch import fnmatch
import argparse
import grace_utils as gr
###############################################################################
## READ ARGUMENTS    
############################################################################### 
## Define arguments   
parser = argparse.ArgumentParser(description=' This programs aim to compute kbr and kbra residuals from monthly mascon solutions.')
parser.add_argument("datadir",default="/scratch/compute1/julia/solutions/monthly_solutions_10cm_filt/",type=str,help="complete path to the directory that contains your monthly solutions")
parser.add_argument("outputdir",default="/scratch/compute1/julia/solutions/monthly_residudals",type=str,help="complete path to the directory that will contain your daily residuals")
parser.add_argument("--masconfile",default='/scratch/compute1/julia/solutions/masconfiles/mascons_stage4_V003a',type=str,help="mascon file ")
parser.add_argument("--regfile",default='/scratch/compute1/julia/solutions/masconfiles/mascons_stage4_V003a.reg',type=str,help="regularisation file")
parser.add_argument("--start",default=0,type=int,help="first month of the timeseries (0 = Dec. 2002)")
parser.add_argument("--end",default=1,type=int,help="last month of the timeseries (if end > nb months in timeseries, the program will stop at the last month of the timeseries)")
## Read arguments
args = parser.parse_args()
## Declare arguments
outputdir = args.outputdir
datadir = args.datadir
startloop = args.start
endloop = args.end
masconfile=args.masconfile
regfile=args.regfile
## data dir

###############################################################################
# List all cmd file in subdirectories
if outputdir[-1]!='/':
    outputdir=outputdir+'/'  
if datadir[-1]!='/':
    datadir=datadir+'/' 
pattern_cmd=datadir+'????/??/addnorm.cmd'
cmd_filelist=[]
pattern_vcv=datadir+'????/??/addnorm*.vcv'
vcv_filelist=[]
for path, subdirs, files in os.walk(datadir):
    for name in files:
        if fnmatch(os.path.join(path, name), pattern_cmd): 
            cmd_filelist.append(os.path.join(path, name))
        if fnmatch(os.path.join(path, name), pattern_vcv): 
            vcv_filelist.append(os.path.join(path, name))            
cmd_filelist.sort()
vcv_filelist.sort()
nb_months=len(cmd_filelist)
if endloop>nb_months:
    endloop=nb_months
###############################################################################
# Create a directory for all residual
if os.path.isdir(outputdir):
    pass
else:
    os.system("mkdir %s"%(outputdir))
# Loop on each month       
for cpt in np.arange(startloop,endloop):
    myyear=cmd_filelist[cpt][-19:-15]
    mymonth=cmd_filelist[cpt][-14:-12]
    # mission setings depending on the year
    if int(myyear)<2018:
        mission=0
        graceorb_cmd = "GRACE.input"
        gracefit_cmd = "/scratch/compute1/julia/solutions/cmdfiles/grace_step1.cmd"
        sat1='A'
        sat2='B'
        RL_num='02'
        duration=86280
    else:
        mission=1
        graceorb_cmd = "GRACE_FO.input"
        gracefit_cmd = "/scratch/compute1/julia/solutions/cmdfiles/gracefo_step1.cmd"
        sat1='C'
        sat2='D'
        RL_num='04'
        duration=86280
    # create monthly directory and copy addnorm
    monthly_dir=outputdir+myyear+mymonth+'/'  
    if os.path.isdir(monthly_dir):
        pass
    else:
        os.system("mkdir %s"%(monthly_dir))
    vcv_filename=vcv_filelist[cpt].split('/')[-1]    
    os.system("cp %s %s"%(vcv_filelist[cpt],monthly_dir+vcv_filename))
    # get days selected in addnorm
    daylist=gr.read_addnorm_cmd(cmd_filelist[cpt])
    for myday in daylist:   
        if myday=='01':
            daily_dir=monthly_dir+myday+'/'
            print('###########################################################################################')
            print('##### ESTIMATE RESIDUALS FROM ADDNORM SOLUTIONS, DAY: %s %s %s: %s #####'%(myyear,mymonth,myday,str(dt.datetime.now())[0:19]))
            print('###########################################################################################')
            if os.path.isdir(daily_dir):
                pass
            else:
                os.system("mkdir %s"%(daily_dir))
            os.chdir(daily_dir)
            # Create setup : get level1B data & auxiliary
            print('## Get data for %s %s %s : %s '%(myday,mymonth,myyear,str(dt.datetime.now())[0:19]))
            os.system('~/ga/com/sh_setup_grace_190927 -t %s %s %s 00 00 00 %s -mission %d -type h5 -acc_model none fft none -mascon_file %s \
                      -extend 21600 -msc_apriori_model msc_apriori_quadratic.model & wait'%(myyear,mymonth,myday,duration,mission,masconfile))
            print('')
            # Get initial conditions from a priori addnorm vcv file
            print('## Create ICs from postfit addnorm values for %s %s %s : %s '%(myday,mymonth,myyear,str(dt.datetime.now())[0:19]))
            os.system('~/ga/util/addnorm2vcv ../%s IC_%s_%s_%s_iter3.vcv %s %s %s pos_E vel_E scl_E bs__E msc_E & wait'%(vcv_filename,myyear,mymonth,myday,myyear,mymonth,myday))
            print('')
            # Prepare graceorb cmd files
            print('## Write input file for %s %s %s : %s '%(myday,mymonth,myyear,str(dt.datetime.now())[0:19]))
            gr.write_setup_iter3(graceorb_cmd,myyear,mymonth,myday)
            print('## %s written'%(graceorb_cmd+'.iter3') )
            print('')
            # Integrate the orbits
            print('## Integrate orbits for %s %s %s : %s '%(myday,mymonth,myyear,str(dt.datetime.now())[0:19]))
            os.system('export OMP_NUM_THREADS=2')
            os.system('~/gt/graceorb_paul_191126/graceorb %s.iter3 in_file_%s_%s-%s-%s_00 0 0 0 & \
                       ~/gt/graceorb_paul_191126/graceorb %s.iter3 in_file_%s_%s-%s-%s_00 0 0 0 & \
                       wait'%(graceorb_cmd,sat1,myyear,mymonth,myday,graceorb_cmd,sat2,myyear,mymonth,myday))
            print('')
            # Rename the orbits
            print('## Rename the orbits for %s %s %s : %s '%(myday,mymonth,myyear,str(dt.datetime.now())[0:19]))
            os.system('mv GTORB_%s-%s-%s_00-00-00_%s_%s.h5  GTORB_%s-%s-%s_00-00-00_%s_%s_iter3.h5'%(myyear,mymonth,myday,sat1,RL_num,myyear,mymonth,myday,sat1,RL_num))
            os.system('mv GTORB_%s-%s-%s_00-00-00_%s_%s.h5  GTORB_%s-%s-%s_00-00-00_%s_%s_iter3.h5'%(myyear,mymonth,myday,sat2,RL_num,myyear,mymonth,myday,sat2,RL_num))
            print('')
            ## Compute residuals with gracefit
            print('## Estimate residuals for msc_%s_%s_%s.vcv: %s '%(myyear,mymonth,myday,str(dt.datetime.now())[0:19]))
            os.system('export OMP_NUM_THREADS=2')
            os.system('~/gt/gracefit/gracefit %s resid_%s_%s_%s.rms 4 GTORB_%s-%s-%s_00-00-00_%s_%s_iter3.h5 \
                      GTORB_%s-%s-%s_00-00-00_%s_%s_iter3.h5 %s %s %s %s 0 0 0 0 0 0 0 0 0 0 0 0 & wait'\
                      %(gracefit_cmd,myyear,mymonth,myday,myyear,mymonth,myday,sat1,RL_num,myyear,mymonth,myday,sat2,RL_num,myyear,mymonth,myday,regfile))
            os.system('mv prefit.res ../prefit_rr_%s_%s_%s.res'%(myyear,mymonth,myday)) 
            os.system('mv prefitRA_filt.res ../prefit_ra_%s_%s_%s.res'%(myyear,mymonth,myday)) 
            print('## prefit_rr_%s_%s_%s.res, prefit_ra_%s_%s_%s.res copied to %s'%(myyear,mymonth,myday,myyear,myyear,mymonth,myday,myyear,mymonth,myday,myyear,mymonth,myday,myyear,mymonth,myday,monthly_dir))
            print('')                   
                        ## Delete daily directory (everything execpt kb & kba residuals)
            print('## Delete daily directtory for %s %s %s %s '%(myday,mymonth,myyear,str(dt.datetime.now())[0:19]))
            os.chdir(monthly_dir)
            os.system('rm -rf %s'%(daily_dir))
            print('## %s deleted'%(daily_dir))
        print('')
## The End