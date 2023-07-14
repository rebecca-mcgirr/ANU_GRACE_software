#!/usr/bin/env python3
from math import sin, cos, sqrt, atan2, radians
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy
import cartopy.mpl.geoaxes
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import argparse
###############################################################################
## READ ARGUMENTS    
############################################################################### 
# Define arguments   
parser = argparse.ArgumentParser(description=' This program will plot an equivalent water height time series for the mascon closest to the selected latitude and longitude.')
parser.add_argument("datafile",default='ANU_GRACE_solutions_australia_BUGGED.nc',type=str,help="complete path to your monthly solutions (netcdf)")
parser.add_argument("lat",default=None,type=float,help="latitude")
parser.add_argument("lon",default=None,type=float,help="longitude")
parser.add_argument("--figname",default='my_grace_ts.png',type=str,help="name of your output figure")

# Read arguments
args = parser.parse_args()
# Declare arguments
ncfile = args.datafile
point_lat= args.lat
point_lon= args.lon
figname= args.figname
# Use positive longitudes
if point_lon<0:
    point_lon=360+point_lon
# Read data file
nc=Dataset(ncfile,'r')
decyear=nc.variables['decyear'][:]
lat=nc.variables['lat'][:]
lon=nc.variables['lon'][:]
ewh=nc.variables['ewh'][:]
nc.close()
del nc, ncfile
decyear=decyear.__array__()
ewh=ewh.__array__()
NT=len(decyear)
NB=len(lat)
# select closest primary
distance=99999e3*np.ones(NB)
R = 6373.0
lat1 = radians(point_lat)
lon1 = radians(point_lon)
for cpt in np.arange(NB):
    lat2 = radians(lat[cpt])
    lon2 = radians(lon[cpt])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))  
    distance[cpt] = R * c
index=np.where(distance==np.min(distance))[0][0]
mylon=lon[index]
if mylon>180:
    mylon=mylon-360
mylat=lat[index]

### plot
fig=plt.figure(figsize=(10,7))
ax = fig.subplots()
ax.plot(decyear,ewh[index,:],'k.--',linewidth=1.5,label='total ewh anomaly (m)')
ax.set_ylim(np.min(ewh[index,:]),np.max(ewh[index,:])+1/2*(np.max(ewh[index,:])-np.min(ewh[index,:])))
ax.set_xlim(2002.85,np.max(decyear)+0.25)
ax.legend(loc=2)
ax.set_ylabel ('EWH (m)')
axins = inset_axes(ax, width="30%", height="30%", loc="upper right", 
                   axes_class=cartopy.mpl.geoaxes.GeoAxes, 
                   axes_kwargs=dict(map_projection=cartopy.crs.PlateCarree()))
axins.add_feature(cartopy.feature.COASTLINE)
axins.set_extent([110, 160, -45, -7])
axins.plot(mylon,mylat, 'ro', markersize=7)

ax.set_title('latitude: %3.2f, longitude : %3.2f'%(mylat,mylon))
plt.savefig(figname)
plt.show()
