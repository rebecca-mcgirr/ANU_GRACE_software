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
##########################
## INPUT FILES
grace_file='/scratch/compute1/julia/solutions/cmdfiles/grace.cmd'
gracefo_file='/scratch/compute1/julia/solutions/cmdfiles/gracefo.cmd'
regularisation_file='/scratch/compute1/julia/masconfiles/mascons_stage4_V003a.reg'
datadir='/Gdata/GRAV/GRACE/solutions/raijin/batch5'
mydir='/scratch/compute1/julia/solutions/monthly_solutions'
ncfile='/scratch/compute1/julia/solutions/ncfiles/daily_solutions_raijin_batch5.nc' 
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
for cpt in np.arange(nb_months):
    startdate=gr.toYearFraction(dt.date(firstdate.year,firstdate.month,1)+relativedelta(months=cpt))
    enddate=gr.toYearFraction(dt.date(firstdate.year,firstdate.month,1)+relativedelta(months=cpt+1))
    if startdate<2018:
        dates_select=daydate[(decyear>=startdate)&(decyear<enddate)&(kbra<=3)&(rms_postfit<0.01)] 
    else:
        dates_select=daydate[(decyear>=startdate)&(decyear<enddate)&(kbra<=2)&(rms_postfit<0.01)]
    nbdays=len(dates_select)
    if nbdays>=7: # number of minimum daily solutions in a month
        mydate=dt.date(2000,1,1)+relativedelta(days=dates_select[0])
        print('%4d %02d: %d days: %s'%(mydate.year,mydate.month,nbdays,str(dt.datetime.now())[0:19]))
        monthly_dir="%s/%4d/%02d"%(mydir,mydate.year,mydate.month)
        if os.path.isdir("%s/%4d"%(mydir,mydate.year)):
            os.system("mkdir %s"%(monthly_dir))
        else:
            os.system("mkdir %s/%4d"%(mydir,mydate.year))
            os.system("mkdir %s"%(monthly_dir))      
        if mydate.year<2018:
            os.system("ln -s %s %s/gracefit.cmd"%(grace_file,monthly_dir))
        else:
            os.system("ln -s %s %s/gracefit.cmd"%(gracefo_file,monthly_dir))
        ofile=open('%s/addnorm.cmd'%(monthly_dir), 'w')
        ofile.write('Y Y N\n')
        ofile.write('{:2.1f}\n'.format(nbdays))
        ofile.write('{:s}\n'.format(regularisation_file))
        ofile.write('{:2d}\n'.format(nbdays))
        for cpt in np.arange(nbdays):
            mydate=dt.date(2000,1,1)+relativedelta(days=dates_select[cpt])
            ofile.write('{:s}\n'.format(datadir+'/%4d/msc%4d-%02d-%02d/msc_%4d_%02d_%02d_iter3.norm'%(mydate.year,mydate.year,mydate.month,mydate.day,mydate.year,mydate.month,mydate.day)))
        ofile.close() 
        os.chdir(monthly_dir)
        #os.system("~/gt/gracefit/addnorm addnorm.cmd gracefit.cmd")
######################
## The End
