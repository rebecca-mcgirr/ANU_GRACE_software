#!/usr/bin/env python3
from math import sin, cos, sqrt, atan2, radians
from netCDF4 import Dataset
import numpy as np

import matplotlib as mpl
mpl.use('qt5agg') # Case insensitive. Must be before backends+pyplot import
import matplotlib.rcsetup as rcsetup
import matplotlib.backend_bases as bb

import matplotlib.pyplot as plt
import cartopy
import cartopy.mpl.geoaxes
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import argparse
import os
from grace_utils import read_lowpass_file
from gracetools.io.mascons import mascons

def closest_point(point_lat,point_lon,lat,lon):

    import numpy as np
    from math import sin, cos, sqrt, atan2, radians

    # select closest primary
    distance=99999e3*np.ones(NB)
    R = 6373.0
    lat1 = radians(point_lat)
    lon1 = radians(point_lon)
    print("searching for nearest primary mascon to point",point_lat, point_lon)
    for cpt in np.arange(NB):
        lat2 = radians(lat[cpt])
        lon2 = radians(lon[cpt])
        dlon = lon2 - lon1
        dlat = lat2 - lat1
        a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
        c = 2 * atan2(sqrt(a), sqrt(1 - a))  
        distance[cpt] = R * c
    index=np.where(distance==np.min(distance))[0][0]
    print("min distance = ",distance[index])
    mylon=lon[index]
    if mylon>180:
        mylon=mylon-360
    mylat=lat[index]
    print("closest primary mascon: ",index,mylat,mylon)
    return (index,mylat,mylon)

def read_nc_coords(ncfile):

    nc=Dataset(ncfile,'r')
    decyear=nc.variables['decyear'][:]
    lat=nc.variables['lat'][:]
    lon=nc.variables['lon'][:]
    ewh=nc.variables['ewh'][:]
    sigma=nc.variables['sigma'][:]
    nc.close()
    del nc, ncfile
    decyear=decyear.__array__()
    ewh=ewh.__array__()
    sigma=sigma.__array__()
    
    return(decyear,lat,lon,ewh,sigma)


###############################################################################
## READ ARGUMENTS    
############################################################################### 
# Define arguments   
parser = argparse.ArgumentParser(description=' This program will plot an equivalent water height time series for the mascon closest to the selected latitude and longitude.')
parser.add_argument("--series1",default='none',type=str,help="netcdf file of monthly solutions")
parser.add_argument("--lat",default=None,type=float,help="latitude")
parser.add_argument("--lon",default=None,type=float,help="longitude")
parser.add_argument("--figname",default='my_grace_ts.png',type=str,help="name of your output figure")
parser.add_argument("--ymin",default='0.0',type=float,help="min y axis value")
parser.add_argument("--ymax",default='0.0',type=float,help="max y axis value")
parser.add_argument("--series2",default="none",type=str,help="add a second series to plot")
parser.add_argument("--series3",default="none",type=str,help="add a third  series to plot")
parser.add_argument("--series4",default="none",type=str,help="add a fourth series to plot")
parser.add_argument("--series5",default="none",type=str,help="add a fifth  series to plot")
parser.add_argument("--series6",default="none",type=str,help="add a sixth  series to plot")
parser.add_argument("--series7",default="none",type=str,help="add a seventh  series to plot")
parser.add_argument("--series8",default="none",type=str,help="add an eigth  series to plot")
parser.add_argument("--series9",default="none",type=str,help="add a  ninth  series to plot")
parser.add_argument("--apr_series1",default="none",type=str,help="add an apriori model file to plot the a priori mascon values")
parser.add_argument("--apr_mascon1",default="none",type=str,help="mascon file for apriori model")
parser.add_argument("--apr_series2",default="none",type=str,help="add an apriori model file to plot the a priori mascon values")
parser.add_argument("--apr_mascon2",default="none",type=str,help="mascon file for apriori model")
parser.add_argument("--large",default="yes",type=str,help="plot large single figure")


# Read arguments
args = parser.parse_args()
# Declare arguments
ncfile = args.series1
point_lat= args.lat
point_lon= args.lon
figname= args.figname
ymin = args.ymin
ymax = args.ymax
# additional series to plot
series1 = args.series1
series2 = args.series2
series3 = args.series3
series4 = args.series4
series5 = args.series5
series6 = args.series6
series7 = args.series7
series8 = args.series8
series9 = args.series9
apr_series1 = args.apr_series1
apr_mascon1 = args.apr_mascon1
apr_series2 = args.apr_series2
apr_mascon2 = args.apr_mascon2
large = args.large

print(series2, series3, series4)

# Use positive longitudes
if point_lon<0:
    point_lon=360+point_lon

# Read data file
decyear,lat,lon,ewh,sigma = read_nc_coords(ncfile)
NT=len(decyear)
NB=len(lat)

## select closest primary
index,mylat,mylon = closest_point(point_lat,point_lon,lat,lon)

print(index,mylat,mylon)

#############
### plot  ###
#############
print("Plotting time series")
legend_loc = 0  # 0: best; 2: upper-left; 3: lower left

if large == "yes" :
  fig=plt.figure(figsize=(15,10))
else:
  fig=plt.figure(figsize=(8,5))

ax = fig.subplots()

for label in (ax.get_xticklabels() + ax.get_yticklabels()):
	label.set_fontsize(16)


ax.errorbar(decyear,ewh[index,:] , yerr=sigma[index,:],ls='--',)
ax.plot(decyear,ewh[index,:],'b.--',linewidth=1.5,label=series1)
# PT220325: if ymin and ymax are both zero then set the limits
if ymin == 0. and ymax == 0:
    #ax.set_ylim(np.min(ewh[index,:])-1/2*(np.max(ewh[index,:])-np.min(ewh[index,:])),np.max(ewh[index,:])+1/2*(np.max(ewh[index,:])-np.min(ewh[index,:])))
    ax.set_ylim(np.min(ewh[index,:])-0.1,np.max(ewh[index,:])+1/2*(np.max(ewh[index,:])-np.min(ewh[index,:])))
else:
    ax.set_ylim(ymin,ymax)


#ax.set_xlim(np.min(decyear)-0.25,np.max(decyear)+0.25)
ax.set_xlim(2002.,2023.)
ax.legend(loc=legend_loc)
ax.set_ylabel ('EWH (m)',fontsize=16)
#axins = inset_axes(ax, width="30%", height="30%", loc="upper right", 
#                   axes_class=cartopy.mpl.geoaxes.GeoAxes, 
#                   axes_kwargs=dict(map_projection=cartopy.crs.PlateCarree()))
#axins.add_feature(cartopy.feature.COASTLINE)
# PT220127: fit it around the point requested
#axins.set_extent([110, 160, -45, -7])
minlon = mylon - 25
maxlon = mylon + 25
minlat = mylat - 16
maxlat = mylat + 16
if minlat < -90:
    minlat = -90
if maxlat > 90.:
    maxlat = 90.
#axins.set_extent([minlon, maxlon, minlat, maxlat])
#axins.plot(mylon,mylat, 'ro', markersize=7)

ax.set_title('Mascon: %6d,   File: %20s, latitude: %3.2f, longitude : %3.2f'%(index+1,args.series1,mylat,mylon))

#### Additional time series are read and plotted here
if series2 != "none":
    print("Read and plot time series from file: ",series2)
    decyear2,lat2,lon2,ewh2,sigma2 = read_nc_coords(series2)
    NT=len(decyear2)
    NB=len(lat2)

    ## select closest primary
    index2,mylat,mylon = closest_point(point_lat,point_lon,lat2,lon2)

    # plot it
    ax.plot(decyear2,ewh2[index2,:],linewidth=3,label=series2)
    #ax.errorbar(decyear2,ewh2[index2,:] , yerr=sigma2[index2,:])
    ax.legend(loc=legend_loc)
 
if series3 != "none":
    print("Read and plot time series from file: ",series3)
    decyear3,lat3,lon3,ewh3,sigma3 = read_nc_coords(series3)
    NT=len(decyear3)
    NB=len(lat3)

    ## select closest primary
    index3,mylat,mylon = closest_point(point_lat,point_lon,lat3,lon3)

    # plot it
    ax.plot(decyear3,ewh3[index3,:],'m',linewidth=3,label=series3)
    #ax.errorbar(decyear3,ewh3[index3,:] , yerr=sigma3[index3,:])
    ax.legend(loc=legend_loc)
   
if series4 != "none":
    print("Read and plot time series from file: ",series4)
    decyear4,lat4,lon4,ewh4,sigma4 = read_nc_coords(series4)
    NT=len(decyear4)
    NB=len(lat4)

    ## select closest primary
    index4,mylat,mylon = closest_point(point_lat,point_lon,lat4,lon4)

    # plot it
    ax.plot(decyear4,ewh4[index4,:],linewidth=3,label=series4)
    #ax.errorbar(decyear4,ewh4[index4,:] , yerr=sigma4[index4,:])
    ax.legend(loc=legend_loc)
  
if series5 != "none":
    print("Read and plot time series from file: ",series5)
    decyear5,lat5,lon5,ewh5,sigma5 = read_nc_coords(series5)
    NT=len(decyear5)
    NB=len(lat5)

    ## select closest primary
    index5,mylat,mylon = closest_point(point_lat,point_lon,lat5,lon5)

    # plot it
    ax.plot(decyear5,ewh5[index5,:],linewidth=2,label=series5)
    #ax.errorbar(decyear4,ewh4[index4,:] , yerr=sigma4[index4,:])
    ax.legend(loc=legend_loc)
  
if series6 != "none":
    print("Read and plot time series from file: ",series6)
    decyear6,lat6,lon6,ewh6,sigma6 = read_nc_coords(series6)
    NT=len(decyear6)
    NB=len(lat6)

    ## select closest primary
    index6,mylat,mylon = closest_point(point_lat,point_lon,lat6,lon6)

    # plot it
    ax.plot(decyear6,ewh6[index6,:],linewidth=3,label=series6)
    #ax.errorbar(decyear4,ewh4[index4,:] , yerr=sigma4[index4,:])
    ax.legend(loc=legend_loc)
  
if series7 != "none":
    print("Read and plot time series from file: ",series7)
    decyear7,lat7,lon7,ewh7,sigma7 = read_nc_coords(series7)
    NT=len(decyear7)
    NB=len(lat7)

    ## select closest primary
    index7,mylat,mylon = closest_point(point_lat,point_lon,lat7,lon7)

    # plot it
    ax.plot(decyear7,ewh7[index7,:],linewidth=3,label=series7)
    #ax.errorbar(decyear4,ewh4[index4,:] , yerr=sigma4[index4,:])
    ax.legend(loc=legend_loc)
  
if series8 != "none":
    print("Read and plot time series from file: ",series8)
    decyear8,lat8,lon8,ewh8,sigma8 = read_nc_coords(series8)
    NT=len(decyear8)
    NB=len(lat8)

    ## select closest primary
    index8,mylat,mylon = closest_point(point_lat,point_lon,lat8,lon8)

    # plot it
    ax.plot(decyear8,ewh8[index8,:],linewidth=3,label=series8)
    #ax.errorbar(decyear4,ewh4[index4,:] , yerr=sigma4[index4,:])
    ax.legend(loc=legend_loc)
  
if series9 != "none":
    print("Read and plot time series from file: ",series9)
    decyear9,lat9,lon9,ewh9,sigma9 = read_nc_coords(series9)
    NT=len(decyear9)
    NB=len(lat9)

    ## select closest primary
    index9,mylat,mylon = closest_point(point_lat,point_lon,lat9,lon9)

    # plot it
    ax.plot(decyear9,ewh9[index9,:],linewidth=3,label=series9)
    #ax.errorbar(decyear4,ewh4[index4,:] , yerr=sigma4[index4,:])
    ax.legend(loc=legend_loc)
  
if apr_series1 != "none":
    # get the name of the relevant mascon file
    print("Mascon file for ",apr_series1," is ",apr_mascon1)
    print("Reading mascon file:",apr_mascon1)
    msc = mascons(apr_mascon1)
    msc.read_mascons_txt()
        
    # read the coords from the mascon file
    ## select closest primary
    NB=len(msc.Plat)
    index_apr,mylat,mylon = closest_point(point_lat,point_lon,msc.Plat,msc.Plon)
    
    
    # call function to extract the apriori value for this mason
    print("Get a priori information for mascon: ",index_apr)
    apr_dec_yr,apr_density,apr_ann_ampl,apr_lowpass = read_lowpass_file(apr_series1,index_apr+1)
    ann_signal = np.zeros(len(apr_dec_yr))
    for iepoch in range(len(apr_dec_yr)):
        ann_signal[iepoch] = apr_ann_ampl[0]*np.cos(apr_dec_yr[iepoch]*2*np.pi) + apr_ann_ampl[1]*np.sin(apr_dec_yr[iepoch]*2*np.pi)

    # plot it
    ax.plot(apr_dec_yr,apr_lowpass+ann_signal,linewidth=1.5,label="apriori model: "+apr_series1)
    ax.legend(loc=legend_loc)
    
  
if apr_series2 != "none":
    # get the name of the relevant mascon file
    print("Mascon file for ",apr_series2," is ",apr_mascon2)
    print("Reading mascon file:",apr_mascon2)
    msc = mascons(apr_mascon2)
    msc.read_mascons_txt()
        
    # read the coords from the mascon file
    ## select closest primary
    NB=len(msc.Plat)
    index_apr,mylat,mylon = closest_point(point_lat,point_lon,msc.Plat,msc.Plon)
    
    
    # call function to extract the apriori value for this mason
    print("Get a priori information for mascon: ",index_apr)
    apr_dec_yr,apr_density,apr_ann_ampl,apr_lowpass = read_lowpass_file(apr_series2,index_apr+1)
    ann_signal = np.zeros(len(apr_dec_yr))
    for iepoch in range(len(apr_dec_yr)):
        ann_signal[iepoch] = apr_ann_ampl[0]*np.cos(apr_dec_yr[iepoch]*2*np.pi) + apr_ann_ampl[1]*np.sin(apr_dec_yr[iepoch]*2*np.pi)

    # plot it
    ax.plot(apr_dec_yr,apr_lowpass+ann_signal,linewidth=1.5,color='black',label="apriori model: "+apr_series2)
    ax.legend(loc=legend_loc)
    
  



#plt.savefig(figname)
plt.show()
