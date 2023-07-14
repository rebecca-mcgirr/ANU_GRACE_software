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
parser.add_argument("datadir",default="/scratch/compute1/julia/solutions/monthly_solutions",type=str,help="complete path to your monthly solutions")
parser.add_argument("outputdir",default="/scratch/compute1/julia/solutions/monthly_resiudals",type=str,help="complete path to the directory that will contain your monthly solutions")
parser.add_argument("--masconfile",default='/scratch/compute1/julia/solutions/masconfiles/mascons_stage4_V003a',type=str,help="mascon file for the inversion")
parser.add_argument("--start",default=0,type=int,help="first month of the timeseries (0 = Dec. 2002)")
parser.add_argument("--end",default=1,type=int,help="last month of the timeseries (if end > nb months in timeseries, the program will stop at the last month of the timeseries)")
## Read arguments
args = parser.parse_args()
## Declare arguments
filedir = args.datadir
outputdir = args.outputdir
startloop = args.start
endloop = args.end
masconfile=args.masconfile
###############################################################################
# List all cmd file in subdirectories
if outputdir[-1]!='/':
    outputdir=outputdir+'/'
if filedir[-1]=='/':
    pattern=filedir+'????/??/addnorm.cmd'
else:
    pattern=filedir+'/????/??/addnorm.cmd'
filelist=[]
for path, subdirs, files in os.walk(filedir):
    for name in files:
        if fnmatch(os.path.join(path, name), pattern): 
            filelist.append(os.path.join(path, name))
filelist.sort()
nb_months=len(filelist)
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
    myyear=filelist[cpt][-19:-15]
    mymonth=filelist[cpt][-14:-12]
    # mission setings depending on the year
    if int(myyear)<2018:
        mission=0
        graceorb_cmd = "GRACE.input"
        gracefit_cmd = "/scratch/compute1/julia/solutions/cmdfiles/grace.cmd"
        sat1='A'
        sat2='B'
        RL_num='02'
        duration=1000
    else:
        mission=1
        graceorb_cmd = "GRACE_FO.input"
        gracefit_cmd = "/scratch/compute1/julia/solutions/cmdfiles/gracefo.cmd"
        sat1='C'
        sat2='D'
        RL_num='04'
        duration=1000
    # create monthly directory and copy addnorm
    monthly_dir=outputdir+myyear+mymonth+'/'  
    if os.path.isdir(monthly_dir):
        pass
    else:
        os.system("mkdir %s"%(monthly_dir))
    os.system("cp %s %s"%(filelist[cpt][:-3]+'vcv',monthly_dir+'addnorm.vcv'))
    # get days selected in addnorm
    daylist=gr.read_addnorm_cmd(filelist[cpt])
    for myday in daylist:        
        daily_dir=monthly_dir+myday
        print('################################################################')
        print('##### %s %s %s: %s #####'%(myyear,mymonth,myday,str(dt.datetime.now())[0:19]))
        print('################################################################')
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
        print('## Create ICs from a priori addnorm values for %s %s %s : %s '%(myday,mymonth,myyear,str(dt.datetime.now())[0:19]))
        os.system('~/ga/util/addnorm2vcv ../addnorm.vcv IC_%s_%s_%s_iter3.vcv %s %s %s pos_A vel_A scl_A bs__A msc_A & wait'%(myyear,mymonth,myday,myyear,mymonth,myday))
        print('')
        # Prepare graceorb cmd files
        print('## Write input file for %s %s %s : %s '%(myday,mymonth,myyear,str(dt.datetime.now())[0:19]))
        gr.write_setup_iter3(graceorb_cmd,myyear,mymonth,myday)
        print('## %s written'%(graceorb_cmd+'iter3') )
        print('')
        # Integrate the orbits
        print('## Integrate orbits for %s %s %s : %s '%(myday,mymonth,myyear,str(dt.datetime.now())[0:19]))
        os.system('export OMP_NUM_THREADS=10')
        os.system('~/ga/graceorb/graceorb %s.iter3 in_file_%s_%s-%s-%s_00 0 0 0 & \
                   ~/ga/graceorb/graceorb %s.iter3 in_file_%s_%s-%s-%s_00 0 0 0 & \
                   wait'%(graceorb_cmd,sat1,myyear,mymonth,myday,graceorb_cmd,sat2,myyear,mymonth,myday))
        print('')
        # Rename the orbits
        print('## Rename the orbits for %s %s %s : %s '%(myday,mymonth,myyear,str(dt.datetime.now())[0:19]))
        os.system('mv GTORB_%s-%s-%s_00-00-00_%s_%s.h5  GTORB_%s-%s-%s_00-00-00_%s_%s_iter3.h5'%(myyear,mymonth,myday,sat1,RL_num,myyear,mymonth,myday,sat1,RL_num))
        os.system('mv GTORB_%s-%s-%s_00-00-00_%s_%s.h5  GTORB_%s-%s-%s_00-00-00_%s_%s_iter3.h5'%(myyear,mymonth,myday,sat2,RL_num,myyear,mymonth,myday,sat2,RL_num))
        print('')
        ## Get daily msc solution from addnorm (daily estimates for satellites, monthly estimates for mascons)
        print('## Create daily vcv from estimated (postfit) addnorm values for %s %s %s : %s '%(myday,mymonth,myyear,str(dt.datetime.now())[0:19]))
        os.system('~/ga/util/addnorm2vcv ../addnorm.vcv msc_%s_%s_%s_iter3.vcv %s %s %s pos_E vel_E scl_E bs__E msc_E & wait'%(myyear,mymonth,myday,myyear,mymonth,myday))          
        print('')
        ## Compute residuals with gracefit
        print('## Estimate residuals from vcv file with gracefit for %s %s %s : %s '%(myday,mymonth,myyear,str(dt.datetime.now())[0:19]))
        os.system('export OMP_NUM_THREADS=10')
        os.system('~/gt/gracefit/gracefit %s resid_%s_%s_%s.rms 4 GTORB_%s-%s-%s_00-00-00_%s_%s_iter3.h5 \
                  GTORB_%s-%s-%s_00-00-00_%s_%s_iter3.h5 %s %s %s msc_%s_%s_%s_iter3.vcv 0 0 0 0 0 0 0 0 0 0 0 0 & wait'\
                  %(gracefit_cmd,myyear,mymonth,myday,myyear,mymonth,myday,sat1,RL_num,myyear,mymonth,myday,sat2,RL_num,myyear,mymonth,myday,myyear,mymonth,myday))
        print('')
        ## Move residuals in monthly directory
        print('## Copy daily residuals in monthly directory for %s %s %s %s '%(myday,mymonth,myyear,str(dt.datetime.now())[0:19]))
        os.system('mv plt_resid_%s_%s_%s.kb ../'%(myyear,mymonth,myday)) 
        os.system('mv plt_resid_%s_%s_%s.kba ../'%(myyear,mymonth,myday)) 
        print('## plt_resid_%s_%s_%s.kb and plt_resid_%s_%s_%s.kba copied to %s'%(myyear,mymonth,myday,myyear,mymonth,myday,monthly_dir) )
        print('')
        ## Delete daily directory (everything execpt kb & kba residuals)
        print('## Delete daily directtory for %s %s %s %s '%(myday,mymonth,myyear,str(dt.datetime.now())[0:19]))
        os.chdir(monthly_dir)
        os.system('rm -rf %s'%(daily_dir))
        print('## %s deleted'%(daily_dir))
        print('')
## The End