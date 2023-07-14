#!/usr/bin/env python3
import os
from fnmatch import fnmatch
import numpy as np
import scipy as sc
from scipy.signal import butter
import datetime as dt
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.basemap import Basemap
import argparse
from math import sin, cos, sqrt, atan2, radians
##
import grace_utils as gr
###############################################################################
## READ ARGUMENTS    
############################################################################### 
# Define arguments   
parser = argparse.ArgumentParser(description=' This programs will plot the data, low[ass filtered data, annual signal and residual for a prticular mascon selected by latitude and longitude.')
parser.add_argument("datadir",default="/scratch/compute1/julia/solutions/monthly_solutions",type=str,help="complete path to your monthly solutions")
parser.add_argument("pattern",default='????/??/addnorm*.fit',type=str,help="pattern describing fitfiles")
parser.add_argument("lat",default=None,type=float,help="latitude")
parser.add_argument("lon",default=None,type=float,help="longitude")
parser.add_argument("--masconfile",default='/scratch/compute1/geodynamics/software/gtgk_1064/grace/tables/mascons_stage4_V003a',type=str,help="complete path to yourmascon file (default=mascons_stage4_V003a)")
parser.add_argument("--cutoff",default=2,type=float,help="cutoff period in years (default=2)")
parser.add_argument("--figname",default='lowpass_ts.png',type=str,help="name of figure output")
# Read arguments
args = parser.parse_args()
# Declare arguments
filedir = args.datadir
mypattern= args.pattern
masconfile= args.masconfile
myfigname= args.figname
Tc=args.cutoff
point_lat= args.lat
point_lon= args.lon
if point_lon>180:
    point_lon=point_lon-360
###############################################################################
##  Read fitfile
###############################################################################
## List files
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
print('Read \'%s\' to get the number of mascons: %s'%(masconfile,str(dt.datetime.now())[0:19]))
f=open(masconfile)
lines=f.readlines()
f.close()
nb_msc=int(lines[0].split()[1])
del f, lines

print('Found %d mascons in \'%s\': %s'%(nb_msc,masconfile,str(dt.datetime.now())[0:19]))
# Create arrays to store data
print('Read %d fitfiles and create nd-arrays with prefit and postfit statistics, and, satellite and mascon parameters: %s'%(nbt,str(dt.datetime.now())[0:19]))
decyear=np.zeros(nbt)
params_msc=9999*np.ones((nbt,nb_msc,4))
nbdays=np.zeros(nbt)
## Read addnom fitfiles
tcount=-1
for fitfile in filelist:
    if gr.read_addnorm_fit(fitfile,nb_msc)==None:
        pass
    else:
        tcount=tcount+1
        decyear[tcount],params_msc[tcount,:,:],nbdays[tcount],_=gr.read_addnorm_fit(fitfile,nb_msc)
print('Found %d non-empty fitfiles in \'%s\': %s'%(tcount,filelist[0],str(dt.datetime.now())[0:19]))
# Extract time and postfit values
decyear=decyear[:tcount]
postfit=params_msc[:tcount,:,2]
NT=len(decyear) # nb times
NB=len(postfit[0,:]) # nb mascons
###############################################################################
##  Select mascon given lat, lon
###############################################################################
print('Select mascon near lat = %5.2f , lon = %5.2f: %s'%(point_lat,point_lon,str(dt.datetime.now())[0:19]))
#mascon file
_,Pdata,_,_=gr.read_mascon_file(masconfile)
nbp=len(Pdata[:,0])
lat=Pdata[:,4]
lon=Pdata[:,5]
lon[lon>180]=lon[lon>180]-360*np.ones(len(lon[lon>180]))
# select
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
#index=gr.find_closest_mascon(point_lat,point_lon,lat,lon)
data=postfit[:,index]-np.mean(postfit[:,index])*np.ones(NT)
mylat=lat[index]
mylon=lon[index]
print('Closest mascon at lat = %5.2f , lon = %5.2f: %s'%(mylat,mylon,str(dt.datetime.now())[0:19]))
###############################################################################
##  Compute annual, lowpass, HF signals for mascon closest to lat, lon
###############################################################################
print('Separate annual, lowpass and hf signals for %d mascons: %s'%(NB,str(dt.datetime.now())[0:19]))
# annual sinusoid
model=np.ones((NT,3))
model[:,0]=np.cos(decyear*2*np.pi)
model[:,1]=np.sin(decyear*2*np.pi)
# OLS inversion of an annual sinusoid
coefs=np.linalg.lstsq(model,data)
coefs=coefs[0]
annual_predictions=np.dot(model,coefs)
annual_residuals=data-annual_predictions
# linear interpolation
decyear_i=np.arange(decyear[0]-1/12,decyear[-1]+1/12,1/12)
fi=sc.interpolate.interp1d(decyear, annual_residuals,kind='linear',fill_value='extrapolate')
data_i=fi(decyear_i)
# low pass filter on interpolated data
fs = 12
fc = 1/Tc  # Cut-off frequency of the filter
w = fc / (fs / 2) # Normalize the frequency
b, a = butter(5, w, 'low')
lowpass_2yr_i = sc.signal.filtfilt(b, a, data_i)
# interp back to original sampling
fi_low=sc.interpolate.interp1d(decyear_i, lowpass_2yr_i,kind='linear',fill_value='extrapolate')
lowpass_2yr=fi_low(decyear)
###############################################################################
##  Plot
###############################################################################
print('Plot times series at lat: %5.2f, lon; %5.2f: %s'%(mylat,mylon,str(dt.datetime.now())[0:19]))
fig=plt.figure(figsize=(10,7))
ax = fig.subplots()
ax.plot(decyear,postfit[:,index],'k.--',linewidth=1.5,label='postfit values')
ax.plot(decyear,annual_predictions,'b.--',linewidth=1.5,label='annual')
ax.plot(decyear,lowpass_2yr+np.mean(postfit[:,index])*np.ones(NT),'r.--',linewidth=1.5,label='lowpass (%5.2f yr)'%(Tc))
ax.plot(decyear,data-annual_predictions-lowpass_2yr,'c.--',linewidth=1.5,label='high frequencies')
ax.set_ylim(np.min(data),np.max(data)+1/6*(np.max(data)-np.min(data)))
ax.set_xlim(2002.85,2019.85)
ax.legend(loc=2)
ax.set_ylabel ('EWH (m)')
axins = inset_axes(ax, width="30%", height="30%", loc="upper right")
m = Basemap(llcrnrlon=-180,llcrnrlat=-90,urcrnrlon=180,urcrnrlat=90, ax=axins)
plon,plat = m(mylon,mylat) # Transforms lat/lon into plotting coordinates for projection
m.drawcoastlines(zorder=3)
m.drawmapboundary(zorder=1)
m.plot(plon,plat,'ro ')
ax.set_title('latitude: %3.2f, longitude : %3.2f'%(mylat,mylon))
plt.show()
plt.savefig(myfigname)
print('%s: You\'re all done! '%(str(dt.datetime.now())[0:19]))
###############################################################################
##  End
###############################################################################

