#!/usr/bin/env python3
import numpy as np
import datetime as dt
from dateutil.relativedelta import relativedelta
from math import ceil
import sys
sys.path.append('~/gt/util/pycodes')
import grace_utils as gr
import os
from fnmatch import fnmatch
import argparse
###############################################################################
## READ ARGUMENTS    
############################################################################### 
# Define arguments   
parser = argparse.ArgumentParser(description=' This programs aim to stack daily ANU mascons solutions towards monthly solutions.')
parser.add_argument("datadir",default="/Gdata/GRAV/GRACE/solutions/raijin/batch5",type=str,help="complete path to your daily solutions")
parser.add_argument("pattern",default="????/??/addnorm*.fit",type=str,help="fitfile pattern in datadir")

parser.add_argument("outputdir",default="/scratch/compute1/julia/solutions/monthly_solutions",type=str,help="complete path to the directory that will contain your monthly solutions")
parser.add_argument("--masconfile",default='/home/geod/julia/gg/grace/tables/mascons_stage4_V003a',type=str,help="mascon file for the inversion")
parser.add_argument("--regfile",default='~/gg/grace/tables/mascons_stage4_V003a.reg',type=str,help="regularisation file for the inversion")
parser.add_argument("--prefix",default='10cm',type=str,help="prefix to identify regulaiation in output vcv and fit file")
parser.add_argument("--start",default=0,type=int,help="first month of the timeseries (0 = Dec. 2002)")
parser.add_argument("--end",default=999,type=int,help="last month of the timeseries (if end > nb months in timeseries, the program will stop at the last month of the timeseries)")
# Read arguments
args = parser.parse_args()
# Declare arguments
datadir = args.datadir
mydir = args.outputdir
if mydir[-1]!='/':
    mydir=mydir+'/'	
mypattern=args.pattern
regularisation_file = args.regfile
prefix_reg = args.prefix
startloop = args.start
endloop = args.end
mascon_file = args.masconfile
###############################################################################
## INPUT FILES
grace_file='/scratch/compute1/julia/solutions/cmdfiles/grace_ga.cmd'
gracefo_file='/scratch/compute1/julia/solutions/cmdfiles/gracefo_ga.cmd'
#mascon_file='/scratch/compute1/julia/solutions/masconfiles/mascons_stage4_V003a'
##########################
### read daily fitfiles 
if datadir[-1]=='/':
    pattern=datadir+mypattern
else:
    pattern=datadir+'/'+mypattern
print('Make list of files fitting pattern %s: %s'%(pattern,str(dt.datetime.now())[0:19]))
filelist=[]
for path, subdirs, files in os.walk(datadir):
    for name in files:
        if fnmatch(os.path.join(path, name), pattern):
            filelist.append(os.path.join(path, name))
filelist.sort()
nbt=len(filelist)
# Read masconfile to get number of mascons
print('Read \'%s\' to get the number of mascons: %s'%(mascon_file,str(dt.datetime.now())[0:19]))
f=open(mascon_file)
lines=f.readlines()
f.close()
nb_msc=int(lines[0].split()[1])
del f,lines
# Create arrays to store daily data
print('Read %d fitfiles and create nd-arrays with prefit and postfit satatistics, and, satellite and mascon parameters: %s'%(nbt,str(dt.datetime.now())[0:19]))
datelist=np.zeros(nbt)
decyear=np.zeros(nbt)
prestats=9999*np.ones((nbt,16))
poststats=9999*np.ones((nbt,16))
params_sat1=9999*np.ones((nbt,12,5))
params_sat2=9999*np.ones((nbt,12,5))
params_msc=9999*np.ones((nbt,nb_msc,5))
tcount=-1
for fitfile in filelist:
    if gr.read_fitfile(fitfile)==None:
        pass
    else:
        tcount=tcount+1
        mydate,prestats[tcount,:],poststats[tcount,:],params_sat1[tcount,:,:],params_sat2[tcount,:,:],params_msc[tcount,:,:]=gr.read_fitfile(fitfile)
        delta=dt.date(int(mydate[0:4]),int(mydate[5:7]),int(mydate[8:10]))-dt.date(2000,1,1)
        datelist[tcount]=delta.days
        decyear[tcount]=gr.toYearFraction(dt.datetime.strptime(mydate,'%Y-%m-%d'))
print('Found %d non-empty fitfiles in \'%s\': %s'%(tcount,filelist[0],str(dt.datetime.now())[0:19]))
### get daily dates and daily rms for kbra postfit 
daydate=datelist[poststats[:,15]<9999]
decyear=decyear[poststats[:,15]<9999]
kbra=poststats[poststats[:,15]<9999,15]
postfit=params_msc[poststats[:,15]<9999,:,2]
######################
## Create monthly directories with a selection of norm files to run adnorm
nb_months=ceil((decyear[-1]-decyear[0])*12) 
firstdate=dt.date(2000,1,1)+dt.timedelta(days=int(daydate[0]))
if endloop >=nb_months:
    endloop=nb_months
for cpt in np.arange(startloop,endloop):
    startdate=gr.toYearFraction(dt.date(firstdate.year,firstdate.month,1)+relativedelta(months=cpt))
    enddate=gr.toYearFraction(dt.date(firstdate.year,firstdate.month,1)+relativedelta(months=cpt+1))
    if startdate<2018:
        dates_select=daydate[(decyear>=startdate)&(decyear<enddate)&(kbra<=3)] 
    else:
        dates_select=daydate[(decyear>=startdate)&(decyear<enddate)&(kbra<=2)]
    nbdays=len(dates_select)
    if nbdays>=7: # number of minimum daily solutions in a month
        mydate=dt.date(2000,1,1)+relativedelta(days=dates_select[0])
        print('%4d %02d: %d days: %s'%(mydate.year,mydate.month,nbdays,str(dt.datetime.now())[0:19]))
        monthly_dir="%s%4d/%02d"%(mydir,mydate.year,mydate.month)
        if os.path.isfile(monthly_dir+'/addnorm_%s_%4d_%02d.fit'%(prefix_reg,mydate.year,mydate.month)):
            print('fitfile already exists for %4d %02d. Skip.'%(mydate.year,mydate.month))
        else:
            if os.path.isdir("%s/%4d"%(mydir,mydate.year)):
                os.system("mkdir %s"%(monthly_dir))
            else:
                os.system("mkdir %s/%4d"%(mydir,mydate.year))
                os.system("mkdir %s"%(monthly_dir))      
            if mydate.year<2018:
                os.system("cp -s %s %s/gracefit.cmd"%(grace_file,monthly_dir))
            else:
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
            myprefix='addnorm_%s_%4d_%02d'%(prefix_reg,mydate.year,mydate.month)
            os.system("~/gt/gracefit/addnorm addnorm.cmd gracefit.cmd %s %s"%(mascon_file,myprefix))
            #os.system("mv addnorm.fit addnorm_%s_%4d_%02d.fit"%(prefix_reg,mydate.year,mydate.month))
            #os.system("mv addnorm.vcv addnorm_%s_%4d_%02d.vcv"%(prefix_reg,mydate.year,mydate.month))
######################
## The End
