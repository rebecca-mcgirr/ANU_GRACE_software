#!/usr/bin/env python3
from math import sin, cos, sqrt, atan2, radians
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import cartopy
import cartopy.mpl.geoaxes
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import argparse
import os
import h5py
import sys
from gracetools.io.mascons import mascons

###############################################################################
## READ ARGUMENTS    
############################################################################### 
# Define arguments   
parser = argparse.ArgumentParser(description=' This program will plot sea level height time series for the ocean.')
parser.add_argument("--hdf5_file",default='none',type=str,help="hdf5 file of monthly solutions")
parser.add_argument("--figname",default='my_grace_ts.png',type=str,help="name of your output figure")
parser.add_argument("--ymin",default='-20.0',type=float,help="min y axis value")
parser.add_argument("--ymax",default='40.0',type=float,help="max y axis value")
parser.add_argument("--series2",default="none",type=str,help="add a second series to plot")
parser.add_argument("--series3",default="none",type=str,help="add a third  series to plot")
parser.add_argument("--series4",default="none",type=str,help="add a fourth series to plot")
parser.add_argument("--apr_series",default="none",type=str,help="add an apriori model file to plot the a priori mascon values")
parser.add_argument("--mascon_file",default="mascons_stage5_V006_200km",type=str,help="ascii mascon file")
parser.add_argument("--plot_GSFC",default="Y",type=str,help="add GSFC ocean time series to plot")
parser.add_argument("--plot_JPL",default="Y",type=str,help="add JPL ocean time series to plot")

# Read arguments
args = parser.parse_args()
# Declare arguments
hdf5_file = args.hdf5_file
figname= args.figname
ymin = args.ymin
ymax = args.ymax
# additional series to plot
series1 = args.hdf5_file
series2 = args.series2
series3 = args.series3
series4 = args.series4
apr_series = args.apr_series
mascon_file = args.mascon_file
plot_GSFC = args.plot_GSFC
plot_JPL = args.plot_JPL

print(series2, series3, series4)


# read the ascii mascons file
print("Reading mascon file:",mascon_file)
msc = mascons(mascon_file)
msc.read_mascons_txt()
    
# Read hdf5 data file
with h5py.File(hdf5_file, "r") as h5f:
    sol = h5f['solution/ewh'][:,:]
    sigma_ewh = h5f['solution/sigma_ewh'][:,:]
    year = h5f['time/year'][:]
    decyear = h5f['time/decyear'][:]
    month = h5f['time/month'][:]
    
# get the dimensions of epochs and mascons
n_epochs = len(decyear)
n_mascons = len(sol[1,:])
print("There are ",n_epochs," epochs and ",n_mascons," mascons in file: ",hdf5_file)

# loop over the epochs
gsl = []
Ant = []
Grn = []
for iepoch in range(n_epochs):
    sum_vol = 0.
    ocean_area = 0.
    sum_Ant_vol = 0.
    Ant_area = 0.
    sum_Grn_vol = 0.
    Grn_area = 0.
    for imsc in range (n_mascons):
        #print(imsc,msc.Pdensity[imsc],msc.Parea[imsc])
        if msc.Pdensity[imsc] > 1010. :
            sum_vol += sol[iepoch,imsc] * msc.Parea[imsc]
            ocean_area += msc.Parea[imsc]
        elif msc.Pdensity[imsc] < 1010  and msc.Plat[imsc] < -60. :
            sum_Ant_vol += sol[iepoch,imsc] * msc.Parea[imsc]
            Ant_area += msc.Parea[imsc]
        elif msc.Pdensity[imsc] < 1010  and msc.Pdesc[imsc] == "Greenland" :
            sum_Grn_vol += sol[iepoch,imsc] * msc.Parea[imsc]
            Grn_area += msc.Parea[imsc]
         
    print(decyear[iepoch],sum_vol/ocean_area, sum_Ant_vol/1e9, sum_Grn_vol/1e9)
    gsl.append(1000.0*sum_vol/ocean_area)
    Ant.append(sum_Ant_vol/1e9)
    Grn.append(sum_Grn_vol/1e9)

print("Area of Antarctica:",Ant_area*1e-12," million km^2")
print("Area of the oceans:",ocean_area*1e-12," million km^2")

# plot series2
if series2 != "none" :
    # Read hdf5 data file
    with h5py.File(series2, "r") as h5f:
        sol2 = h5f['solution/ewh'][:,:]
        sigma_ewh2 = h5f['solution/sigma_ewh'][:,:]
        decyear2 = h5f['time/decyear'][:]

    n_epochs2 = len(decyear2)
    gsl2 = []
    Ant2 = []
    for iepoch in range(n_epochs2):
        sum_vol2 = 0.
        sum_Ant_vol2 = 0.
        for imsc in range (n_mascons):
            if msc.Pdensity[imsc] > 1010. :
                sum_vol2 += sol2[iepoch,imsc] * msc.Parea[imsc]

            elif msc.Plat[imsc] < -60. :
                sum_Ant_vol2 += sol2[iepoch,imsc] * msc.Parea[imsc]
            
        print(decyear2[iepoch],sum_vol2/ocean_area, sum_Ant_vol2/1e9)
        gsl2.append(1000.0*sum_vol2/ocean_area)
        Ant2.append(sum_Ant_vol2/1e9)




# plot the time series
f1 = plt.figure(1,figsize=(12,6))
plt.plot(decyear,gsl)
plt.title(hdf5_file)
plt.xlabel('Year')
plt.ylim(ymin,ymax)
plt.ylabel('GSL (mm EWH)')

# prepare the GSFC data for plotting if needed
if plot_GSFC == "Y":
    gsfc_ocean = np.loadtxt("/Mdata/GRACE/software/tables/data/GSFC_ocean_ts.dat",skiprows = 12)
    plt.plot(gsfc_ocean[:,0],gsfc_ocean[:,1]*10)

if plot_JPL == "Y":
    jpl_ocean = np.loadtxt("/Mdata/GRACE/software/tables/data/JPL_ocean.dat",skiprows = 13)
    # the JPL epoch is days since 01/01/2002. Convert to decimal year
    decyear_JPL = np.zeros(len(jpl_ocean[:,0]))
    for i in range (len(jpl_ocean[:,0])):
        decyear_JPL[i] = 2002.0 + jpl_ocean[i,0]/365.0
        
    plt.plot(decyear_JPL,jpl_ocean[:,1]*10)

# plot series2 if needed
if series2 != "none" :
    plt.plot(decyear2,gsl2)
    

plt.legend([hdf5_file,"GSFC","JPL",series2])

# Antarctic plot   
f2 = plt.figure(2,figsize=(12,6))    
plt.plot(decyear,Ant)
plt.title(hdf5_file)
plt.xlabel('Year')
plt.ylim(-2000,1200)
plt.ylabel('Antarctic (GT)')

if plot_GSFC == "Y":
    gsfc_EAnt = np.loadtxt("/Mdata/GRACE/software/tables/data/GSFC_EAnt.dat",skiprows = 12)
    gsfc_WAnt = np.loadtxt("/Mdata/GRACE/software/tables/data/GSFC_WAnt.dat",skiprows = 12)
    gsfc_AntPen = np.loadtxt("/Mdata/GRACE/software/tables/data/GSFC_AntPeninsula.dat",skiprows = 12)
    gsfc_tot = np.zeros(len(gsfc_EAnt))
    gsfc_tot = gsfc_EAnt[:,1]+gsfc_WAnt[:,1]+gsfc_AntPen[:,1] + 30.
    plt.plot(gsfc_WAnt[:,0],gsfc_tot*10. )
    plt.legend([hdf5_file,'GSFC'])

if plot_JPL == "Y":
    jpl_Ant = np.loadtxt("/Mdata/GRACE/software/tables/data/JPL_Antarctic.dat",skiprows = 13)
    plt.plot(decyear_JPL,jpl_Ant[:,1]*100. +300. )
    plt.legend([hdf5_file,'GSFC','JPL'])

# Greenland plot
f3 = plt.figure(3,figsize=(12,6))    
plt.plot(decyear,Grn)
plt.title(hdf5_file)
plt.xlabel('Year')
plt.ylim(-4000,1200)
plt.ylabel('Greenland (GT)')
plt.legend([hdf5_file])

if plot_GSFC == "Y":
    gsfc_Grn = np.loadtxt("/Mdata/GRACE/software/tables/data/GSFC_Greenland.dat",skiprows = 12)
    plt.plot(gsfc_Grn[:,0],gsfc_Grn[:,1]*10. )
    plt.legend([hdf5_file,'GSFC'])

if plot_JPL == "Y":
    jpl_Grn = np.loadtxt("/Mdata/GRACE/software/tables/data/JPL_Greenland.dat",skiprows = 13)
    plt.plot(decyear_JPL,jpl_Grn[:,1]*10. )
    plt.legend([hdf5_file,'GSFC','JPL'])



plt.show()

