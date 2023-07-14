#!/usr/bin/env python3
import os
from fnmatch import fnmatch
import numpy as np
import datetime as dt
import grace_utils as gr
from netCDF4 import Dataset
import argparse
###############################################################################
## READ ARGUMENTS    
############################################################################### 
# Define arguments   
parser = argparse.ArgumentParser(description=' This programs aim to write netcdf from monthly fitfiles.')
parser.add_argument("datadir",default="/scratch/compute1/julia/solutions/monthly_solutions",type=str,help="complete path to your monthly solutions")
parser.add_argument("ncfile",default='/scratch/compute1/julia/solutions/ncfiles/monhly_solutions_raijin_batch5.nc',type=str,help="netcdf file compiling monthly fitfiles")
parser.add_argument("pattern",default='????/??/addnorm*.fit',type=str,help="pattern describing fitfiles")
parser.add_argument("masconfile",default="/scratch/compute1/julia/solutions/masconfiles/mascons_stage4_V003a",type=str,help="complete path to your masconfile")
# Read arguments
args = parser.parse_args()
# Declare arguments
filedir = args.datadir
ncname= args.ncfile
mypattern= args.pattern
masconfile= args.masconfile
###############################################################################
if filedir[-1]=='/':
    pattern=filedir+mypattern
else:
    pattern=filedir+'/'+mypattern
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
_,Pdata,_,_=gr.read_mascon_file(masconfile)
nb_msc=len(Pdata[:,0])
mascon_number=Pdata[:,0]
print('Found %d mascons in \'%s\': %s'%(nb_msc,filelist[0],str(dt.datetime.now())[0:19]))
# Create arrays to store data
print('Read %d fitfiles and create nd-arrays with prefit and postfit satatistics, and, satellite and mascon parameters: %s'%(nbt,str(dt.datetime.now())[0:19]))
decyear=np.zeros(nbt)
params_msc=9999*np.ones((nbt,nb_msc,4))
nbdays=np.zeros(nbt)
tcount=-1
for fitfile in filelist:
    if gr.read_addnorm_fit(fitfile,nb_msc)==None:
        pass
    else:
        tcount=tcount+1
        decyear[tcount],params_msc[tcount,:,:],nbdays[tcount],_=gr.read_addnorm_fit(fitfile,nb_msc)
print('Found %d non-empty fitfiles in \'%s\': %s'%(tcount,filelist[0],str(dt.datetime.now())[0:19]))
################################################################################
# Write netcdf file
print('Create netcdf file %s: %s'%(ncname,str(dt.datetime.now())[0:19]))
ncw=Dataset(ncname,'w',format='NETCDF4')
ncw.description = 'stack of all fitfiles in %s'%(filedir)
# Create dimensions
ncw.createDimension('nbt', tcount) # number of fitfiles each corresponding to a date
ncw.createDimension('nb_msc', nb_msc) # number mascons in each fitfile
# Create variables
# decimal year
decyear_vec= ncw.createVariable('decyear', 'f4', ('nbt'))
decyear_vec.units = 'decimal year'
decyear_vec.description='date in decimal year'
decyear_vec[:]=decyear[:tcount]
# mascon list
mascon_vec= ncw.createVariable('mascons', 'f4', ('nb_msc'))
mascon_vec.units = '-'
mascon_vec.description='primary mascon number'
mascon_vec[:]=mascon_number[:]
# nb days per solution
nbdays_vec= ncw.createVariable('nbdays', 'f4', ('nbt'))
nbdays_vec.units = '-'
nbdays_vec.description='number of days per monthly solution'
nbdays_vec[:]=nbdays[:tcount]
# a priori mascon value
apriori_arr= ncw.createVariable('apriori', 'f4', ('nbt','nb_msc'))
apriori_arr.units = 'EWH(m)'
apriori_arr.description='a_priori equivalent water height'
apriori_arr[:,:]=params_msc[:tcount,:,0]
# postfit mascon value
postfit_arr= ncw.createVariable('postfit', 'f4', ('nbt','nb_msc'))
postfit_arr.units = 'EWH(m)'
postfit_arr.description='postfit equivalent water height'
postfit_arr[:,:]=params_msc[:tcount,:,2]
# sigma mascon value
sigma_arr= ncw.createVariable('sigma', 'f4', ('nbt','nb_msc'))
sigma_arr.units = 'EWH(m)'
sigma_arr.description='misfit'
sigma_arr[:,:]=params_msc[:tcount,:,3]
# Close NetCF file
ncw.close() 
print('End fitfiles2netcdf.py: %s'%(str(dt.datetime.now())[0:19]))
