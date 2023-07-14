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
parser.add_argument("outputdir",default="/scratch/compute1/julia/solutions/monthly_residudals",type=str,help="complete path to the directory that will contain your monthly solutions")
parser.add_argument("--masconfile",default='/scratch/compute1/julia/solutions/masconfiles/mascons_stage4_V003a',type=str,help="mascon file for the inversion")
parser.add_argument("--start",default=0,type=int,help="first month of the timeseries (0 = Dec. 2002)")
parser.add_argument("--end",default=1,type=int,help="last month of the timeseries (if end > nb months in timeseries, the program will stop at the last month of the timeseries)")
## Read arguments
args = parser.parse_args()
## Declare arguments
outputdir = args.outputdir
startloop = args.start
endloop = args.end
masconfile=args.masconfile
## data dir
datadir_10cm="/scratch/compute1/julia/solutions/monthly_solutions/"
datadir_bom="/scratch/compute1/julia/solutions/monthly_solutions_bom/"
datadir_20cm="/scratch/compute1/julia/solutions/monthly_solutions_aus20cm/"
###############################################################################
# List all cmd file in subdirectories
if outputdir[-1]!='/':
    outputdir=outputdir+'/'   
pattern=datadir_10cm+'????/??/addnorm.cmd'
filelist_10cm=[]
for path, subdirs, files in os.walk(datadir_10cm):
    for name in files:
        if fnmatch(os.path.join(path, name), pattern): 
            filelist_10cm.append(os.path.join(path, name))
filelist_10cm.sort()
nb_months=len(filelist_10cm)
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
    myyear=filelist_10cm[cpt][-19:-15]
    mymonth=filelist_10cm[cpt][-14:-12]
    # mission setings depending on the year
    if int(myyear)<2018:
        mission=0
        graceorb_cmd = "GRACE.input"
        gracefit_cmd = "/scratch/compute1/julia/solutions/cmdfiles/grace.cmd"
        sat1='A'
        sat2='B'
        RL_num='02'
        duration=86280
    else:
        mission=1
        graceorb_cmd = "GRACE_FO.input"
        gracefit_cmd = "/scratch/compute1/julia/solutions/cmdfiles/gracefo.cmd"
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
    os.system("cp %s %s"%(filelist_10cm[cpt][:-3]+'vcv',monthly_dir+'addnorm_10cm.vcv'))
    # Copy addnorm bom
    vcv_bom=datadir_bom+myyear+'/'+mymonth+'/addnorm.vcv'
    if os.path.isfile(vcv_bom):
        os.system("cp %s %s"%(vcv_bom,monthly_dir+'addnorm_bom.vcv'))
    # Copy addnorm 20cm
    vcv_20cm=datadir_20cm+myyear+'/'+mymonth+'/addnorm.vcv'
    if os.path.isfile(vcv_20cm):
        os.system("cp %s %s"%(vcv_20cm,monthly_dir+'addnorm_20cm.vcv'))    
    # get days selected in addnorm
    daylist=gr.read_addnorm_cmd(filelist_10cm[cpt])
    for myday in daylist:        
        daily_dir=monthly_dir+myday+'/'
        print('######################################################################################')
        print('##### ESTIMATE RESIDUALS FROM ADDNORM SOLUTIONS, DAY: %s %s %s: %s #####'%(myyear,mymonth,myday,str(dt.datetime.now())[0:19]))
        print('######################################################################################')
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
        os.system('~/ga/util/addnorm2vcv ../addnorm_10cm.vcv IC_%s_%s_%s_iter3.vcv %s %s %s pos_A vel_A scl_A bs__A msc_A & wait'%(myyear,mymonth,myday,myyear,mymonth,myday))
        print('')
        # Prepare graceorb cmd files
        print('## Write input file for %s %s %s : %s '%(myday,mymonth,myyear,str(dt.datetime.now())[0:19]))
        gr.write_setup_iter3(graceorb_cmd,myyear,mymonth,myday)
        print('## %s written'%(graceorb_cmd+'.iter3') )
        print('')
        # Integrate the orbits
        print('## Integrate orbits for %s %s %s : %s '%(myday,mymonth,myyear,str(dt.datetime.now())[0:19]))
        os.system('export OMP_NUM_THREADS=2')
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
        print('## Create msc_%s_%s_%s_10cm.vcv from estimated addmorm values in %s : %s '%(myyear,mymonth,myday,monthly_dir+'addnorm_10cm.vcv',str(dt.datetime.now())[0:19]))
        os.system('~/ga/util/addnorm2vcv ../addnorm_10cm.vcv msc_%s_%s_%s_10cm.vcv %s %s %s pos_E vel_E scl_E bs__E msc_E & wait'%(myyear,mymonth,myday,myyear,mymonth,myday)) 
        print('')
        if os.path.isfile(monthly_dir+'addnorm_bom.vcv'): 
            print('## Create msc_%s_%s_%s_bom.vcv from estimated addmorm values in %s : %s '%(myyear,mymonth,myday,monthly_dir+'addnorm_bom.vcv',str(dt.datetime.now())[0:19]))
            os.system('~/ga/util/addnorm2vcv ../addnorm_bom.vcv msc_%s_%s_%s_bom.vcv %s %s %s pos_E vel_E scl_E bs__E msc_E & wait'%(myyear,mymonth,myday,myyear,mymonth,myday))
            print('')
        if os.path.isfile(monthly_dir+'addnorm_20cm.vcv'): 
            print('## Create msc_%s_%s_%s_20cm.vcv from estimated addmorm values in %s : %s '%(myyear,mymonth,myday,monthly_dir+'addnorm_20cm.vcv',str(dt.datetime.now())[0:19]))
            os.system('~/ga/util/addnorm2vcv ../addnorm_20cm.vcv msc_%s_%s_%s_20cm.vcv %s %s %s pos_E vel_E scl_E bs__E msc_E & wait'%(myyear,mymonth,myday,myyear,mymonth,myday))
            print('')
        ## Compute residuals with gracefit
        print('## Estimate residuals for msc_%s_%s_%s_10cm.vcv: %s '%(myyear,mymonth,myday,str(dt.datetime.now())[0:19]))
        os.system('export OMP_NUM_THREADS=2')
        os.system('~/gt/gracefit/gracefit %s resid_%s_%s_%s_10cm.rms 4 GTORB_%s-%s-%s_00-00-00_%s_%s_iter3.h5 \
                  GTORB_%s-%s-%s_00-00-00_%s_%s_iter3.h5 %s %s %s msc_%s_%s_%s_10cm.vcv 0 0 0 0 0 0 0 0 0 0 0 0 & wait'\
                  %(gracefit_cmd,myyear,mymonth,myday,myyear,mymonth,myday,sat1,RL_num,myyear,mymonth,myday,sat2,RL_num,myyear,mymonth,myday,myyear,mymonth,myday))
        os.system('mv prefit.res ../prefit_rr_10cm_%s_%s_%s.res'%(myyear,mymonth,myday)) 
        os.system('mv prefitRA_filt.res ../prefit_ra_10cm_%s_%s_%s.res'%(myyear,mymonth,myday)) 
        print('## prefit_rr_10cm_%s_%s_%s.res and prefit_ra_10cm_%s_%s_%s.res copied to %s'%(myyear,mymonth,myday,myyear,mymonth,myday,monthly_dir))
        print('')      
        if os.path.isfile(daily_dir+'msc_%s_%s_%s_bom.vcv'%(myyear,mymonth,myday)): 
            print('## Estimate residuals for msc_%s_%s_%s_bom.vcv : %s '%(myyear,mymonth,myday,str(dt.datetime.now())[0:19]))
            os.system('export OMP_NUM_THREADS=2')
            os.system('~/gt/gracefit/gracefit %s resid_%s_%s_%s_bom.rms 4 GTORB_%s-%s-%s_00-00-00_%s_%s_iter3.h5 \
                      GTORB_%s-%s-%s_00-00-00_%s_%s_iter3.h5 %s %s %s msc_%s_%s_%s_bom.vcv 0 0 0 0 0 0 0 0 0 0 0 0 & wait'\
                      %(gracefit_cmd,myyear,mymonth,myday,myyear,mymonth,myday,sat1,RL_num,myyear,mymonth,myday,sat2,RL_num,myyear,mymonth,myday,myyear,mymonth,myday))
            os.system('mv prefit.res ../prefit_rr_bom_%s_%s_%s.res'%(myyear,mymonth,myday)) 
            os.system('mv prefitRA_filt.res ../prefit_ra_bom_%s_%s_%s.res'%(myyear,mymonth,myday))
            print('## prefit_rr_bom_%s_%s_%s.res and prefit_ra_bom_%s_%s_%s.res copied to %s'%(myyear,mymonth,myday,myyear,mymonth,myday,monthly_dir))
            print('')      
        if os.path.isfile(daily_dir+'msc_%s_%s_%s_20cm.vcv'%(myyear,mymonth,myday)): 
            print('## Estimate residuals for msc_%s_%s_%s_20cm.vcv : %s '%(myyear,mymonth,myday,str(dt.datetime.now())[0:19]))
            os.system('export OMP_NUM_THREADS=2')
            os.system('~/gt/gracefit/gracefit %s resid_%s_%s_%s_20cm.rms 4 GTORB_%s-%s-%s_00-00-00_%s_%s_iter3.h5 \
                      GTORB_%s-%s-%s_00-00-00_%s_%s_iter3.h5 %s %s %s msc_%s_%s_%s_20cm.vcv 0 0 0 0 0 0 0 0 0 0 0 0 & wait'\
                      %(gracefit_cmd,myyear,mymonth,myday,myyear,mymonth,myday,sat1,RL_num,myyear,mymonth,myday,sat2,RL_num,myyear,mymonth,myday,myyear,mymonth,myday))               
            os.system('mv prefit.res ../prefit_rr_20cm_%s_%s_%s.res'%(myyear,mymonth,myday)) 
            os.system('mv prefitRA_filt.res ../prefit_ra_20cm_%s_%s_%s.res'%(myyear,mymonth,myday))
            print('## prefit_rr_20cm_%s_%s_%s.res and prefit_ra_20cm_%s_%s_%s.res copied to %s'%(myyear,mymonth,myday,myyear,mymonth,myday,monthly_dir))
            print('')                 
                    ## Delete daily directory (everything execpt kb & kba residuals)
        print('## Delete daily directtory for %s %s %s %s '%(myday,mymonth,myyear,str(dt.datetime.now())[0:19]))
        os.chdir(monthly_dir)
        os.system('rm -rf %s'%(daily_dir))
        print('## %s deleted'%(daily_dir))
        print('')
## The End