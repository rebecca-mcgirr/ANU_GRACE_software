#!/usr/bin/env python3
from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import griddata
import grace_utils as gr
import datetime as dt
import pyproj
import argparse
###############################################################################
## READ ARGUMENTS    
############################################################################### 
# Define arguments   
parser = argparse.ArgumentParser(description=' This programs aim to extract regrid the ANU mascons solutions towards a regular grid.')
parser.add_argument("input_file",default="None",type=str,help="complete path to your input solution file (netcdf) ")
parser.add_argument("mascon_file",default="None",type=str,help="complete path to the mascon file associated with the solution")
parser.add_argument("output_file",default="None",type=str,help="complete path to your output solution file (netcdf)")
parser.add_argument("--spacing",default=0.5,type=float,help="resolution of the global grid (degree)")
# Read arguments
args = parser.parse_args()
# Declare arguments
input_file = args.input_file
mascon_file = args.mascon_file
output_file=args.output_file
spacing=args.spacing
###############################################################################
## READ INPUT FILES    
###############################################################################
# solutions
print('Read solutions: %s' %(str(dt.datetime.now())))
nc=Dataset(input_file,'r')
decyear=nc.variables['decyear'][:]
mascons=nc.variables['mascons'][:]
postfit=nc.variables['postfit'][:]
sigma=nc.variables['sigma'][:]
nc.close()
del nc
nbt=len(decyear)
#mascon file
print('Read mascon file: %s' %(str(dt.datetime.now())))
_,Pdata,_,Tdata=gr.read_mascon_file(mascon_file)
nbp=len(Pdata[:,0])
nbtern=len(Tdata[:,0])
###############################################################################
## EXPAND PRIMARIES VALUES AT TERNARIES LOCATIONS   
###############################################################################
print('primary to ternary: %s' %(str(dt.datetime.now())))
lat=Tdata[:,2]
lon=Tdata[:,3]
prim=Tdata[:,9]
region=Tdata[:,9]
Tewh=np.zeros((nbtern,nbt))
Tsigma=np.zeros((nbtern,nbt))
tcount=0
for pcpt in np.arange(nbp):
    ntip=Pdata[pcpt,3]
    Tewh[tcount:tcount+ntip,:]=np.tile(postfit[:,pcpt],(ntip,1))
    Tsigma[tcount:tcount+ntip,:]=np.tile(sigma[:,pcpt],(ntip,1))
    tcount=tcount+ntip
###############################################################################
## REGRID FROM TERNARIES TO REGULAR GRID
##############################################################################
print('ternary to grid: %s' %(str(dt.datetime.now())))   
latgrd=np.r_[-90:90.01:spacing]
longrd=np.r_[0:360.01:spacing]
lon2d,lat2d = np.meshgrid(longrd,latgrd)
ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
x, y = ecef(lon,lat)
xi, yi = ecef(lon2d,lat2d)
grid_ewh=np.ones((len(latgrd),len(longrd),nbt))
grid_sigma=np.ones((len(latgrd),len(longrd),nbt))
for tcpt in np.arange(nbt):
    print (tcpt, dt.datetime.now())
    grid_ewh[:,:,tcpt] = griddata((x,y),Tewh[:,tcpt],(xi,yi),method='linear')
    grid_sigma[:,:,tcpt] = griddata((x,y),Tsigma[:,tcpt],(xi,yi),method='linear')
###############################################################################
## WRITE OUTPUT GRIDS IN NETCDF  
##############################################################################
# Write gridded data in netcdf
print('Create netcdf file %s: %s'%(output_file,str(dt.datetime.now())[0:19]))
ncw=Dataset(output_file,'w',format='NETCDF4')
ncw.description = 'Monthly ANU solutions of the total water storage anomalies derived from GRACE and GRACE-FO satellite missions'
# Create dimensions
ncw.createDimension('nbt', nbt) # number of fitfiles each corresponding to a date
ncw.createDimension('nblat', len(latgrd)) # number mascons in each fitfile
ncw.createDimension('nblon', len(longrd)) # number mascons in each fitfile
# Create variables
# decimal year
decyear_vec= ncw.createVariable('decyear', 'f4', ('nbt'))
decyear_vec.units = 'decimal year'
decyear_vec.description='date in decimal year'
decyear_vec[:]=decyear[:]
# latitude
lat_vec= ncw.createVariable('lat', 'f4', ('nblat'))
lat_vec.units = 'degree'
lat_vec.description='latitude'
lat_vec[:]=latgrd[:]
# longitude
lon_vec= ncw.createVariable('lon', 'f4', ('nblon'))
lon_vec.units = 'degree'
lon_vec.description='longitude'
lon_vec[:]=longrd[:]
# postfit ewh
postfit_arr= ncw.createVariable('postfit', 'f4', (('nblat','nblon','nbt')))
postfit_arr.units = 'EWH(m)'
postfit_arr.description='postfit equivalent water height'
postfit_arr[:,:,:]=grid_ewh[:,:,:]
# sigma mascon value
sigma_arr= ncw.createVariable('sigma', 'f4', (('nblat','nblon','nbt')))
sigma_arr.units = 'EWH(m)'
sigma_arr.description='misfit'
sigma_arr[:,:,:]=grid_sigma[:,:,:]
# Close NetCF file
ncw.close() 
print('End fitfiles2netcdf.py: %s'%(str(dt.datetime.now())[0:19]))   