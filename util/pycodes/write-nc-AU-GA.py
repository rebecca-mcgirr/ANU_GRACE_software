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
# trim arrays for valid months only
decyear=decyear[:tcount]
mascon=mascon_number[:]
nbdays=nbdays[:tcount]
apriori=params_msc[:tcount,:,0]
postfit=params_msc[:tcount,:,2]
sigma=params_msc[:tcount,:,4]
# trim arrays from Jan 2003 to Sep. 2016
postfit=postfit[(decyear>=2003)&(decyear<=2016.8),:]
sigma=sigma[(decyear>=2003)&(decyear<=2016.8),:]
decyear=decyear[(decyear>=2003)&(decyear<=2016.8)]
nbt=len(decyear)
# select australian mascons
Pcode=[]
lat=[]
lon=[]
region=[]
aus_ewh=[]
aus_sigma=[]
area=[]
for pcpt in np.arange(nb_msc):
    if Pdata[pcpt,-1][0:3]=='AUS':
        Pcode.append(Pdata[pcpt,0])
        lat.append(Pdata[pcpt,4])
        lon.append(Pdata[pcpt,5])
        region.append(Pdata[pcpt,-1])
        area.append(Pdata[pcpt,7])
        aus_ewh.append(postfit[:,pcpt])
        aus_sigma.append(sigma[:,pcpt])
lat=np.asarray(lat)
lon=np.asarray(lon)
area=np.asarray(area)
region=np.asarray(region,dtype=object)
aus_ewh=np.vstack(aus_ewh)
aus_sigma=np.vstack(aus_ewh)
nbm_aus=len(lat)
#####
#
ncw=Dataset(ncname,'w',format='NETCDF4')
ncw.description = 'ANU mascon solutions for australian surface catchments iteration 2'
# Create dimensions
ncw.createDimension('nbt', nbt)
ncw.createDimension('nbp', nbm_aus)
# Create variables
# decimal year
decyear_vec= ncw.createVariable('decyear', 'f4', ('nbt'))
decyear_vec.units = 'decimal year'
decyear_vec.description='date in decimal year'
decyear_vec[:]=decyear[:]
# primary mascon
pcode_vec= ncw.createVariable('primary', 'f4', ('nbp'))
pcode_vec.units = '-'
pcode_vec.description='primary mascon number'
pcode_vec[:]=Pcode[:]
# primary latitude
lat_vec= ncw.createVariable('lat', 'f4', ('nbp'))
lat_vec.units = 'degrees'
lat_vec.description='latitude of the primary mascon centroid in WGS84'
lat_vec[:]=lat[:]
# primary longitude
lon_vec= ncw.createVariable('lon', 'f4', ('nbp'))
lon_vec.units = 'degrees'
lon_vec.description='longitude of the primary mascon centroid in WGS84'
lon_vec[:]=lon[:]
# drainage basin
region_vec= ncw.createVariable('drainage_basin', 'str', ('nbp'))
region_vec.units = '-'
region_vec.description='name of the drainage basin'
region_vec[:]=region[:]
# australia ewh
aewh_arr= ncw.createVariable('ewh', 'f4', (('nbp','nbt')))
aewh_arr.units = 'm'
aewh_arr.description='equivalent water height'
aewh_arr[:]=aus_ewh[:,:]
# australia sigma
asigma_arr= ncw.createVariable('sigma', 'f4', (('nbp','nbt')))
asigma_arr.units = '(m)'
asigma_arr.description='uncertainty at one sigma'
asigma_arr[:]=aus_sigma[:,:]
# Close NetCF file
ncw.close()


print('End fitfiles2netcdf.py: %s'%(str(dt.datetime.now())[0:19]))
