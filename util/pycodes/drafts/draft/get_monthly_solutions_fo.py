#!/usr/bin/env python3
from netCDF4 import Dataset
import numpy as np
import datetime as dt
from dateutil.relativedelta import relativedelta
from math import ceil
import sys
sys.path.append('~/gt/util/pycodes')
import grace_utils as gr
import os
import argparse
###############################################################################
## READ ARGUMENTS    
############################################################################### 
# Define arguments   
parser = argparse.ArgumentParser(description=' This programs aim to stack daily ANU mascons solutions towards monthly solutions.')
parser.add_argument("datadir",default="/Gdata/GRAV/GRACE/solutions/raijin/batch5",type=str,help="complete path to your daily solutions")
parser.add_argument("--ncfile",default='/scratch/compute1/julia/solutions/ncfiles/daily_solutions_raijin_batch5.nc',type=str,help="netcdffile compiling daily fitfiles")
parser.add_argument("outputdir",default="/scratch/compute1/julia/solutions/monthly_solutions",type=str,help="complete path to the directory that will contain your monthly solutions")
parser.add_argument("--regfile",default='~/gg/grace/tables/mascons_stage4_V003a.reg',type=str,help="regularisation file for the inversion")
parser.add_argument("--suffix",default='10cm',type=str,help="suffix to identify regulaiation in output vcv and fit file")
parser.add_argument("--start",default=0,type=int,help="first month of the timeseries (0 = Dec. 2002)")
parser.add_argument("--end",default=999,type=int,help="last month of the timeseries (if end > nb months in timeseries, the program will stop at the last month of the timeseries)")
# Read arguments
args = parser.parse_args()
# Declare arguments
datadir = args.datadir
ncfile = args.ncfile
mydir = args.outputdir
if mydir[-1]!='/':
    mydir=mydir+'/'	
regularisation_file = args.regfile
suffix_reg = args.suffix
startloop = args.start
endloop = args.end
###############################################################################
## INPUT FILES
#grace_file='/scratch/compute1/julia/solutions/cmdfiles/grace_ga.cmd'
gracefo_file='/scratch/compute1/julia/solutions/cmdfiles/gracefo_ga.cmd'
##########################
## read daily fitfiles in netcdf
nc=Dataset(ncfile,'r')
daydate_all=nc.variables['date'][:]
decyear_all=nc.variables['decyear'][:]
mascons=nc.variables['mascons'][:]
prefit_rms_all=nc.variables['prefit-stats'][:]
postfit_rms_all=nc.variables['postfit-stats'][:]
params_sat1_all=nc.variables['params_sat1'][:]
params_sat2_all=nc.variables['params_sat2'][:]
params_msc_all=nc.variables['params_msc'][:]
nc.close()
del nc, ncfile
##########################
## get dates and rms pn postfit kbra
daydate=daydate_all[postfit_rms_all[:,15]<9999]
decyear=decyear_all[postfit_rms_all[:,15]<9999]
kbra=postfit_rms_all[postfit_rms_all[:,15]<9999,15]
postfit=params_msc_all[postfit_rms_all[:,15]<9999,:,2]
nb_msc=len(mascons)
nbd=len(decyear)
rms_postfit=np.zeros(nbd)
for cpt in np.arange(nbd):
    rms_postfit[cpt]=np.sqrt(np.sum(np.square(postfit[cpt,:])))/nb_msc
######################
## Create monthly directories with a selection of norm files to run adnorm
nb_months=ceil((decyear[-1]-decyear[0])*12) 
firstdate=dt.date(2000,1,1)+dt.timedelta(days=int(daydate[0]))
if endloop >=nb_months:
    endloop=nb_months
for cpt in np.arange(startloop,endloop):
    startdate=gr.toYearFraction(dt.date(firstdate.year,firstdate.month,1)+relativedelta(months=cpt))
    enddate=gr.toYearFraction(dt.date(firstdate.year,firstdate.month,1)+relativedelta(months=cpt+1))
    if startdate>2018:
        dates_select=daydate[(decyear>=startdate)&(decyear<enddate)&(kbra<=1.5)]
        nbdays=len(dates_select)
    else:
        nbdays=0
    if nbdays>=7: # number of minimum daily solutions in a month
        mydate=dt.date(2000,1,1)+relativedelta(days=dates_select[0])
        print('%4d %02d: %d days: %s'%(mydate.year,mydate.month,nbdays,str(dt.datetime.now())[0:19]))
        monthly_dir="%s%4d/%02d"%(mydir,mydate.year,mydate.month)
        if os.path.isdir("%s/%4d"%(mydir,mydate.year)):
            os.system("mkdir %s"%(monthly_dir))
        else:
            os.system("mkdir %s/%4d"%(mydir,mydate.year))
            os.system("mkdir %s"%(monthly_dir))      
        os.system("cp -s %s %s/gracefit.cmd"%(gracefo_file,monthly_dir))
        ofile=open('%s/addnorm.cmd'%(monthly_dir), 'w')
        ofile.write('Y Y N\n')
        ofile.write('{:2.1f}\n'.format(nbdays))
        ofile.write('{:s}\n'.format(regularisation_file))
        ofile.write('{:2d}\n'.format(nbdays))
        for cpt in np.arange(nbdays):
            mydate=dt.date(2000,1,1)+relativedelta(days=dates_select[cpt])
            if datadir[-1]=='/':
                ofile.write('{:s}\n'.format(datadir+'%4d/msc%4d-%02d-%02d/msc_%4d_%02d_%02d_iter3.norm'%(mydate.year,mydate.year,mydate.month,mydate.day,mydate.year,mydate.month,mydate.day)))
            else:
                ofile.write('{:s}\n'.format(datadir+'/%4d/msc%4d-%02d-%02d/msc_%4d_%02d_%02d_iter3.norm'%(mydate.year,mydate.year,mydate.month,mydate.day,mydate.year,mydate.month,mydate.day)))
        ofile.close() 
        os.chdir(monthly_dir)
        os.system("~/gt/gracefit/addnorm_no_vcv_matrix addnorm.cmd gracefit.cmd")
        os.system("mv addnorm.fit addnorm_%s_%4d_%02d.fit"%(suffix_reg,mydate.year,mydate.month))
        os.system("mv addnorm.vcv addnorm_%s_%4d_%02d.vcv"%(suffix_reg,mydate.year,mydate.month))
######################
## The End
