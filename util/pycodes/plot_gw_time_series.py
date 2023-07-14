#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 15:38:24 2020

@author: julia pfeffer
"""

from math import sin, cos, sqrt, atan2, radians
from netCDF4 import Dataset
import numpy as np
import argparse


def read_arguments():
    # Define arguments   
    parser = argparse.ArgumentParser(description=' This program will plot an equivalent water height time series for the selected water layers at the selected latitude and longitude.')
    parser.add_argument("grace_file",default='None',type=str,help="complete path to your monthly grace solutions (netcdf)")    
    parser.add_argument("sm_file",default='None',type=str,help="complete path to your monthly soil moisture estimates (netcdf)")
    parser.add_argument("lat",default=None,type=float,help="latitude")
    parser.add_argument("lon",default=None,type=float,help="longitude")
    parser.add_argument("--layer1",default='SWS',type=str,help="name of the first water layer ")
    parser.add_argument("--layer2",default='None',type=str,help="name of the second water layer to added to the previous layer")
    parser.add_argument("--layer3",default='None',type=str,help="name of the third water layer to added to the previous layers")
    parser.add_argument("--layer4",default='None',type=str,help="name of the fourth water layer to added to the previous layers")
    parser.add_argument("--layer5",default='None',type=str,help="name of the fifthth water layer to added to the previous layers")
    parser.add_argument("--layer6",default='None',type=str,help="name of the sixth water layer to added to the previous layers")
    parser.add_argument("--figname",default='my_sm_ts.png',type=str,help="name of your output figure")
    
    # Read arguments
    args = parser.parse_args()
    # Declare arguments
    gfile = args.grace_file
    smfile = args.sm_file
    point_lat= args.lat
    point_lon= args.lon
    # Use positive longitudes
    if point_lon<0:
        point_lon=360+point_lon
    layer1 = args.layer1
    layer2 = args.layer2
    layer3 = args.layer3
    layer4 = args.layer4
    layer5 = args.layer5
    layer6 = args.layer6
    figname= args.figname
    
    return gfile,smfile,point_lat,point_lon,layer1,layer2,layer3,layer4,layer5,layer6,figname

def read_grace_data(gfile):
    # Read data file
    nc=Dataset(gfile,'r')
    decyear=nc.variables['decyear'][:]
    lat=nc.variables['lat'][:]
    lon=nc.variables['lon'][:]
    ewh=nc.variables['ewh'][:]
    nc.close()
    del nc, gfile
    decyear=decyear.__array__()
    ewh=ewh.__array__()
    return lat,lon,decyear,ewh

def select_closest_primary(point_lat,point_lon,lat,lon):
        # select closest primary
    NB=len(lat)    
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
    return index

def read_smfile(smfile,layer1,layer2,layer3,layer4,layer5,layer6):
    # Read data file
    nc=Dataset(smfile,'r')
    lon=nc.variables['lon'][:]
    lat=nc.variables['lat'][:]
    decyear=nc.variables['decyear'][:]
    NT=len(decyear)
    NB=len(lat)
    if layer1=='None':
        mylayer1=np.zeros((NB,NT))
    else:        
        mylayer1=nc.variables[layer1][:,:]
        
    if layer2=='None':
        mylayer2=np.zeros((NB,NT))
    else:        
        mylayer2=nc.variables[layer2][:,:]
        
    if layer3=='None':
        mylayer3=np.zeros((NB,NT))
    else:        
        mylayer3=nc.variables[layer3][:,:]
        
    if layer4=='None':
        mylayer4=np.zeros((NB,NT))
    else:        
        mylayer4=nc.variables[layer4][:,:]
    
    if layer5=='None':
        mylayer5=np.zeros((NB,NT))
    else:        
        mylayer5=nc.variables[layer5][:,:]
    
    if layer6=='None':
        mylayer6=np.zeros((NB,NT))
    else:        
        mylayer6=nc.variables[layer6][:,:]

    return lat,lon,decyear,mylayer1,mylayer2,mylayer3,mylayer4,mylayer5,mylayer6

def from_sm_to_grace_times(decyear,decyeari,index,mylayer1i,mylayer2i,mylayer3i,mylayer4i,mylayer5i,mylayer6i):
    l1i=mylayer1i[index,:]
    mylayer1=np.interp(decyear,decyeari, l1i)
    
    l2i=mylayer2i[index,:]
    mylayer2=np.interp(decyear,decyeari, l2i)
    
    l3i=mylayer3i[index,:]
    mylayer3=np.interp(decyear,decyeari, l3i)
    
    l4i=mylayer4i[index,:]
    mylayer4=np.interp(decyear,decyeari, l4i)
    
    l5i=mylayer5i[index,:]
    mylayer5=np.interp(decyear,decyeari, l5i)
    
    l6i=mylayer6i[index,:]
    mylayer6=np.interp(decyear,decyeari, l6i)
    
    return mylayer1,mylayer2,mylayer3,mylayer4,mylayer5,mylayer6

def plot_layers(decyear,lat,lon,index,ewh,mylayer1,layer1,mylayer2,layer2,mylayer3,layer3,mylayer4,layer4,mylayer5,layer5,mylayer6,layer6,figname):
    import matplotlib
#    matplotlib.use('GTK')
    import os
    havedisplay = "DISPLAY" in os.environ
    if not havedisplay:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import cartopy
    import cartopy.mpl.geoaxes
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    # latitude, longitude of mascon centroid
    mylon=lon[index]
    if mylon>180:
        mylon=mylon-360
    mylat=lat[index]
    # total of all layers
    mysoiltotal=mylayer1+mylayer2+mylayer3+mylayer4+mylayer5+mylayer6
    
    fig=plt.figure(figsize=(9,9))
    ax = fig.subplots()
    
    ax.plot(decyear,ewh[index,:],'k.--',linewidth=1.5,label='total water storage')
    ax.plot(decyear,mysoiltotal,'b.--',linewidth=1.5,label='soil water storage')
    ax.plot(decyear,ewh[index,:]-mysoiltotal,'r.--',linewidth=1.5,label='ground water storage')
    
    minvalue=np.min([np.min(ewh[index,:]),np.min(mysoiltotal),np.min(ewh[index,:]-mysoiltotal)])    
    maxvalue=np.max([np.max(ewh[index,:]),np.max(mysoiltotal),np.max(ewh[index,:]-mysoiltotal)])   
    
    ax.set_ylim(minvalue,maxvalue+1/4*(maxvalue-minvalue))
    ax.set_xlim(np.min(decyear)-0.1,np.max(decyear)+0.1)
    ax.legend(loc=2)
    ax.set_ylabel ('EWH (m)')
    ax.grid()
    axins = inset_axes(ax, width="30%", height="30%", loc="upper right", 
                        axes_class=cartopy.mpl.geoaxes.GeoAxes, 
                        axes_kwargs=dict(map_projection=cartopy.crs.PlateCarree()))
    axins.add_feature(cartopy.feature.COASTLINE)
    axins.set_extent([110, 155, -45, -8])
    axins.plot(mylon,mylat, 'ro', markersize=7)
    
    ax.set_title('latitude: %3.2f, longitude : %3.2f'%(mylat,mylon))
    plt.savefig(figname)
    plt.show()

def main():
    #read arguments
    gfile,smfile,point_lat,point_lon,layer1,layer2,layer3,layer4,layer5,layer6,figname=read_arguments()
    # read grace file
    lat,lon,decyear,ewh=read_grace_data(gfile)
    # read soil moisture file
    lat,lon,decyeari,mylayer1i,mylayer2i,mylayer3i,mylayer4i,mylayer5i,mylayer6i=read_smfile(smfile,layer1,layer2,layer3,layer4,layer5,layer6)
    # find closest mascon to coordinates pair
    index=select_closest_primary(point_lat,point_lon,lat,lon)
    # interpolate soil moisture data at grace times
    mylayer1,mylayer2,mylayer3,mylayer4,mylayer5,mylayer6=from_sm_to_grace_times(decyear,decyeari,index,mylayer1i,mylayer2i,mylayer3i,mylayer4i,mylayer5i,mylayer6i)
    ### plot
    plot_layers(decyear,lat,lon,index,ewh,mylayer1,layer1,mylayer2,layer2,mylayer3,layer3,mylayer4,layer4,mylayer5,layer5,mylayer6,layer6,figname)
 
if __name__ == '__main__':
    main() 