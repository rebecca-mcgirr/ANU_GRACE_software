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
parser = argparse.ArgumentParser(description=' This program reads all daily kba files in a directory and convert them into one netcdf.')
parser.add_argument("dirname",default="None",type=str,help="complete path to your directory containing  kba files (will run through subdirectories) ")
parser.add_argument("--outfile",default="daily_kba.nc",type=str,help="complete path to your output netcdf file (default='daily_kba.nc')")
### Read arguments
args = parser.parse_args()
### Declare arguments
filedir = args.dirname
ncname = args.outfile
###############################################################################
print('Start kba2netcdf.py: %s'%(str(dt.datetime.now())[0:19]))
# List all fit file in subdirectories
if filedir[-1]=='/':
    pattern=filedir+'????/msc????-??-??/plt_msc_????_??_??_iter3.kba'
else:
    pattern=filedir+'/????/msc????-??-??/plt_msc_????_??_??_iter3.kba'
print('Make list of files fitting pattern %s: %s'%(pattern,str(dt.datetime.now())[0:19]))
filelist=[]
for path, subdirs, files in os.walk(filedir):
    for name in files:
        if fnmatch(os.path.join(path, name), pattern): 
            filelist.append(os.path.join(path, name))
filelist.sort()
nbt=len(filelist) 
# Create arrays to store data
print('Read %d fitfiles and create nd-arrays with timevec, decyear, X, Y, Z, prefit, postfit: %s'%(nbt,str(dt.datetime.now())[0:19]))
timevec=np.zeros((nbt,17257))
decyear=np.zeros((nbt,17257))
X=np.zeros((nbt,17257))
Y=np.zeros((nbt,17257))
Z=np.zeros((nbt,17257))
prefit=np.zeros((nbt,17257))
postfit=np.zeros((nbt,17257))
tcount=-1
for kbafile in filelist:
    if gr.read_kbafile(kbafile)==None:
        pass
    else:
        tcount=tcount+1
        timevec[tcount,:],decyear[tcount,:],X[tcount,:],Y[tcount,:],Z[tcount,:],prefit[tcount,:],postfit[tcount,:]=gr.read_kbafile(kbafile)
print('Found %d non-empty kba files in \'%s\': %s'%(tcount,filelist[0],str(dt.datetime.now())[0:19]))
################################################################################
# Write netcdf file
print('Create netcdf file %s: %s'%(ncname,str(dt.datetime.now())[0:19]))
ncw=Dataset(ncname,'w',format='NETCDF4')
ncw.description = 'stack of all kba files in %s'%(filedir)
# Create dimensions
ncw.createDimension('nbt', tcount) # number of fitfiles each corresponding to a date
ncw.createDimension('nb_epochs', 17257) # params of the inversion (a priori, adjust, postfit, sigma, frac)
# Create variables
# date
date_vec= ncw.createVariable('date', 'f4', ('nbt'))
date_vec.units = 'days'
date_vec.description='number of days since 1st Jan. 2000 (first Jan. 2000 = 0)'
date_vec[:]=timevec[:tcount,:]
# decimal year
decyear_vec= ncw.createVariable('decyear', 'f4', ('nbt'))
decyear_vec.units = 'decimal year'
decyear_vec.description='date in decimal year'
decyear_vec[:]=decyear[:tcount,:]
# lat
lat_arr= ncw.createVariable('lat', 'f4', ('nbt','nb_epochs'))
lat_arr.units = '(deg)'
lat_arr.description='latitude'
lat_arr[:,:]=X[:tcount,:]
# lon
lon_arr= ncw.createVariable('lon', 'f4', ('nbt','nb_epochs'))
lon_arr.units = '(deg)'
lon_arr.description='longitude'
lon_arr[:,:]=Y[:tcount,:]
# elevation
alt_arr= ncw.createVariable('altitude', 'f4', ('nbt','nb_epochs'))
alt_arr.units = '(m)'
alt_arr.description='distance from the satellite to the Earth center of mass'
alt_arr[:,:]=Z[:tcount,:]
# prefit
prefit_arr= ncw.createVariable('prefit', 'f4', ('nbt','nb_epochs'))
prefit_arr.units = '(um/s2)'
prefit_arr.description='prefit acceleration'
prefit_arr[:,:]=prefit[:tcount,:]
# postfit
postfit_arr= ncw.createVariable('postfit', 'f4', ('nbt','nb_epochs'))
postfit_arr.units = '(um/s2)'
postfit_arr.description='postfit acceleration'
postfit_arr[:,:]=postfit[:tcount,:]
# Close NetCF file
ncw.close() 
print('End kba2netcdf.py: %s'%(str(dt.datetime.now())[0:19]))