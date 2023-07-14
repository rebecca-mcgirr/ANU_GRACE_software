#!/usr/bin/env python3
from fnmatch import fnmatch
from netCDF4 import Dataset
import os
import numpy as np
import datetime as dt
import grace_utils as gr
import argparse
################################################################################
### READ ARGUMENTS    
################################################################################ 
### Define arguments   
parser = argparse.ArgumentParser(description=' This program reads all daily fitfiles in a directory and convert them into one netcdf.')
parser.add_argument("dirname",default="None",type=str,help="complete path to your directory containing fitfiles (will run through subdirectories) ")
parser.add_argument("--outfile",default="daily_fits.nc",type=str,help="complete path to your output netcdf file (default='daily_fits.nc')")
### Read arguments
args = parser.parse_args()
### Declare arguments
filedir = args.dirname
ncname = args.outfile
###############################################################################
print('Start fitfiles2netcdf.py: %s'%(str(dt.datetime.now())[0:19]))
# List all fit file in subdirectories
if filedir[-1]=='/':
    pattern=filedir+'????/msc????-??-??/msc_????_??_??_iter3.fit'
else:
    pattern=filedir+'/????/msc????-??-??/msc_????_??_??_iter3.fit'
print('Make list of files fitting pattern %s: %s'%(pattern,str(dt.datetime.now())[0:19]))
filelist=[]
for path, subdirs, files in os.walk(filedir):
    for name in files:
        if fnmatch(os.path.join(path, name), pattern): 
            filelist.append(os.path.join(path, name))
filelist.sort()
nbt=len(filelist) 
# Read first file to get number of mascons
print('Read \'%s\' to get the number of mascons: %s'%(filelist[0],str(dt.datetime.now())[0:19]))
f=open(filelist[0])
lines=f.readlines()
f.close()
if lines[-1]=='\n':
    nb_msc=int(lines[-2].split()[1][2:])
else:
    nb_msc=int(lines[-1].split()[1][2:]) 
mascon_number=np.zeros(nb_msc)
pcount=0
for line in lines[72:]:
    if line=='\n':
        pass
    else:
        mascon_number[pcount]=int(line.split()[1][2:])
        pcount=pcount+1
    
del lines,f
print('Found %d mascons in \'%s\': %s'%(nb_msc,filelist[0],str(dt.datetime.now())[0:19]))
# Create arrays to store data
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
################################################################################
# Write netcdf file
print('Create netcdf file %s: %s'%(ncname,str(dt.datetime.now())[0:19]))
ncw=Dataset(ncname,'w',format='NETCDF4')
ncw.description = 'stack of all fitfiles in %s'%(filedir)
# Create dimensions
ncw.createDimension('nbt', tcount) # number of fitfiles each corresponding to a date
ncw.createDimension('nb_msc', nb_msc) # number mascons in each fitfile
ncw.createDimension('nb_stats', 16) # nummber of statistics on prefit and postfit
ncw.createDimension('nbl_sat', 12) # number of parameters per satellite
ncw.createDimension('nb_params', 5) # params of the inversion (a priori, adjust, postfit, sigma, frac)
# Create variables
# date
date_vec= ncw.createVariable('date', 'f4', ('nbt'))
date_vec.units = 'days'
date_vec.description='number of days since 1st Jan. 2000 (first Jan. 2000 = 0)'
date_vec[:]=datelist[:tcount]
# decimal year
decyear_vec= ncw.createVariable('decyear', 'f4', ('nbt'))
decyear_vec.units = 'decimal year'
decyear_vec.description='date in decimal year'
decyear_vec[:]=decyear[:tcount]
# mascon list
mascon_vec= ncw.createVariable('mascons', 'f4', ('nb_msc'))
mascon_vec.units = 'number'
mascon_vec.description='primary mascon number'
mascon_vec[:]=mascon_number
# prefit statistics
prestats_arr= ncw.createVariable('prefit-stats', 'f4', ('nbt','nb_stats'))
prestats_arr.units = '(m)(m)(m)(m/s)(m/s)(m/s)(m)(m)(m)(m/s)(m/s)(m/s)(m)(m/s)(um)(um/s^2)'
prestats_arr.description='prefit rms on position sat1 (x,y & z), velocity sat1 (x,y,z), position sat2 (x,y & z), velocity sat2 (x,y,z), position, velocity, range, acceleration'
prestats_arr[:,:]=prestats[:tcount,:]
# postfit statistics
poststats_arr= ncw.createVariable('postfit-stats', 'f4', ('nbt','nb_stats'))
poststats_arr.units = '(m)(m)(m)(m/s)(m/s)(m/s)(m)(m)(m)(m/s)(m/s)(m/s)(m)(m/s)(um)(um/s^2)'
poststats_arr.description='postfit rms on position sat1 (x,y & z), velocity sat1 (x,y,z), position sat2 (x,y & z), velocity sat2 (x,y,z), position, velocity, range, acceleration'
poststats_arr[:,:]=poststats[:tcount,:]
# parameters of the inversion for first satellite
params_sat1_arr= ncw.createVariable('params_sat1', 'f4', ('nbt','nbl_sat','nb_params'))
params_sat1_arr.units = '(mm)(mm)(mm)(mm/s)(mm/s)(mm/s)(n/a)(n/a)n/a)(um/s^2)(um/s^2)(um/s^2)'
params_sat1_arr.description='# parameters for sat1 : first dim = time, second dim = position (x,y,z), velocity (x,y,z), scale (x,y,z), bias (x,y,z), third dim = a priori, adjust,postfit,sigma and frac.'
params_sat1_arr[:,:]=params_sat1[:tcount,:,:]
# parameters of the inversion for second satellite
params_sat2_arr= ncw.createVariable('params_sat2', 'f4', ('nbt','nbl_sat','nb_params'))
params_sat2_arr.units = '(mm)(mm)(mm)(mm/s)(mm/s)(mm/s)(n/a)(n/a)n/a)(um/s^2)(um/s^2)(um/s^2)'
params_sat2_arr.description='# parameters for sat2 : first dim = time, second dim = position (x,y,z), velocity (x,y,z), scale (x,y,z), bias (x,y,z), third dim = a priori, adjust,postfit,sigma and frac.'
params_sat2_arr[:,:]=params_sat2[:tcount,:,:]
# parameters of the inversion for mascons
params_msc_arr= ncw.createVariable('params_msc', 'f4', ('nbt','nb_msc','nb_params'))
params_msc_arr.units = '(m)'
params_msc_arr.description='parameters for mascons : first dim = time, second dim = mascon number, third dim = a priori, adjust,postfit,sigma and frac'
params_msc_arr[:,:]=params_msc[:tcount,:,:]
# Close NetCF file
ncw.close() 
print('End fitfiles2netcdf.py: %s'%(str(dt.datetime.now())[0:19]))