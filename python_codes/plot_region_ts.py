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
from grace_utils import read_lowpass_file
from scipy.interpolate import interp1d

###############################################################################
## READ ARGUMENTS    
############################################################################### 
# Define arguments   
parser = argparse.ArgumentParser(description=' This program will plot sea level height time series for the ocean.')
parser.add_argument("--hdf5_file",default='none',type=str,help="hdf5 file of monthly solutions")
parser.add_argument("--region",default='none',type=str,help="region to plot (must match mascon file descriptor region)")
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
parser.add_argument("--lowpass_apr",default="none",type=str,help="name of apriori.txt file to include in plot")
parser.add_argument("--plotsize",default="large",type=str,help="large or small plot sizes")

# Read arguments
args = parser.parse_args()
# Declare arguments
hdf5_file = args.hdf5_file
region = args.region
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
lowpass_apr = args.lowpass_apr
plotsize = args.plotsize

print(series2, series3, series4)


# read the ascii mascons file
print("Reading mascon file:",mascon_file)
msc = mascons(mascon_file)
msc.read_mascons_txt()
    
# Read hdf5 data file
with h5py.File(hdf5_file, "r") as h5f:
    sol = h5f['solution/ewh'][:,:]
    sigma_ewh = h5f['solution/sigma_ewh'][:,:]
    apr = h5f['solution/prefit'][:,:]
    year = h5f['time/year'][:]
    decyear = h5f['time/decyear'][:]
    month = h5f['time/month'][:]
    
# get the dimensions of epochs and mascons
n_epochs = len(decyear)
n_mascons = len(sol[1,:])
print("There are ",n_epochs," epochs and ",n_mascons," mascons in file: ",hdf5_file)

# read series2 if required
if series2 != "none":
    print("reading solution from ",series2)
    with h5py.File(series2, "r") as h5f:
        sol2 = h5f['solution/ewh'][:,:]
        sigma_ewh2 = h5f['solution/sigma_ewh'][:,:]
        apr2 = h5f['solution/prefit'][:,:]
        year2 = h5f['time/year'][:]
        decyear2 = h5f['time/decyear'][:]
        month2 = h5f['time/month'][:]
  

# loop over the epochs
region_tot = []
apr_tot = []
lowpass_tot = []
model_apr = np.zeros((n_mascons,n_epochs))
for iepoch in range(n_epochs):
    sum_region_vol = 0.
    apr_region_vol = 0.
    lowpass_region_vol = 0.
    region_area = 0.
    crds = np.zeros((n_mascons,3))
    nmsc_region = 0
        
    for imsc in range (n_mascons):
        #print(imsc,msc.Pdensity[imsc],msc.Parea[imsc])

        if region == "Oceania" :
            if  msc.Pdesc[imsc] == "Australia" or msc.Pdesc[imsc] == "s_newzeala" or msc.Pdesc[imsc] == "n_newzeala" :
                sum_region_vol += sol[iepoch,imsc] * msc.Parea[imsc]
                apr_region_vol += apr[iepoch,imsc] * msc.Parea[imsc]
                region_area += msc.Parea[imsc]
                crds[nmsc_region,0] = msc.Plon[imsc]
                crds[nmsc_region,1] = msc.Plat[imsc]
                crds[nmsc_region,2] = imsc
                nmsc_region += 1

                if iepoch == 0 and lowpass_apr != "none":
                    apr_dec_yr,apr_density,apr_ann_ampl,apr_lowpass_tmp = read_lowpass_file(lowpass_apr,imsc)
                    apr_lowpass = np.interp(decyear,apr_dec_yr,apr_lowpass_tmp)
                    
                    model_apr[imsc,iepoch] = apr_lowpass[iepoch] #+ apr_ann_ampl[0]*np.cos(decyear[iepoch]*2*np.pi) + apr_ann_ampl[1]*np.sin(decyear[iepoch]*2*np.pi)
                    #model_apr[imsc,iepoch] =  apr_ann_ampl[0]*np.cos(decyear[iepoch]*2*np.pi) + apr_ann_ampl[1]*np.sin(decyear[iepoch]*2*np.pi)
                lowpass_region_vol += model_apr[imsc,iepoch] * msc.Parea[imsc]
                
                
                
        elif region == "Ocean" :
            if  msc.Pdensity[imsc] > 1010. :
                sum_region_vol += sol[iepoch,imsc] * msc.Parea[imsc]
                apr_region_vol += apr[iepoch,imsc] * msc.Parea[imsc]
                region_area += msc.Parea[imsc]

                if iepoch == 0 and lowpass_apr != "none":
                    apr_dec_yr,apr_density,apr_ann_ampl,apr_lowpass_tmp = read_lowpass_file(lowpass_apr,imsc)
                    apr_lowpass = np.interp(decyear,apr_dec_yr,apr_lowpass_tmp)
                    
                    model_apr[imsc,iepoch] = apr_lowpass[iepoch] + apr_ann_ampl[0]*np.cos(decyear[iepoch]*2*np.pi) + apr_ann_ampl[1]*np.sin(decyear[iepoch]*2*np.pi)
                lowpass_region_vol += model_apr[imsc,iepoch] * msc.Parea[imsc]
            
        elif  msc.Pdesc[imsc] == region :
            sum_region_vol += sol[iepoch,imsc] * msc.Parea[imsc]
            apr_region_vol += apr[iepoch,imsc] * msc.Parea[imsc]
            region_area += msc.Parea[imsc]

            if iepoch == 0 and lowpass_apr != "none":
                apr_dec_yr,apr_density,apr_ann_ampl,apr_lowpass_tmp = read_lowpass_file(lowpass_apr,imsc)
                apr_lowpass = np.interp(decyear,apr_dec_yr,apr_lowpass_tmp)
                    
                model_apr[imsc,iepoch] = apr_lowpass[iepoch] + apr_ann_ampl[0]*np.cos(decyear[iepoch]*2*np.pi) + apr_ann_ampl[1]*np.sin(decyear[iepoch]*2*np.pi)
            lowpass_region_vol += model_apr[imsc,iepoch] * msc.Parea[imsc]

    print(decyear[iepoch],sum_region_vol/region_area,sum_region_vol,region_area)
    if region != "Greenland" and region != "Antarctica" :
        region_tot.append(1000. * sum_region_vol/region_area)   # in EWH
        apr_tot.append(1000. * apr_region_vol/region_area)   # in EWH
        lowpass_tot.append(1000. * lowpass_region_vol/region_area)   # in EWH
        label = "EWH (mm)"
    else:
        region_tot.append(1000. * sum_region_vol*1e-12)         # in Gt
        apr_tot.append(1000. * apr_region_vol*1e-12)            # in Gt
        lowpass_tot.append(1000. * lowpass_region_vol*1e-12)            # in Gt
        label = "Gt"
print("Area of region ",region,":",region_area*1e-12," million km^2")


# save the coords of the region if it was Oceania
if region == "Oceania":
    print("Saving Oceania primary mascon coordinates into file Ocean.lonlat")
    np.savetxt("Oceania.lonlat",crds[0:nmsc_region,:])

# add 5 mm to our ocean time series
if region == "Ocean":
    offset = 0.0
    print("Adding ",offset," mm to ANU ocean time series")
    
    # make an array of epoch and ocean mass
    output = np.zeros((len(decyear),2))
    output[:,0] = decyear
    output[:,1] = np.array(region_tot) + offset
    np.savetxt("Ocean.dat",output)
else:
    offset = 0.


# integrate region for series 2
if series2 != "none" :
    n_epochs = len(decyear2)

    region_tot2 = []
    for iepoch in range(n_epochs):
        sum_region_vol2 = 0.
        nmsc_region = 0
        
        for imsc in range (n_mascons):
            if region == "Ocean" :
                if  msc.Pdensity[imsc] > 1010. :
                    sum_region_vol2 += sol2[iepoch,imsc] * msc.Parea[imsc]
            
            elif  msc.Pdesc[imsc] == region :
                sum_region_vol2 += sol2[iepoch,imsc] * msc.Parea[imsc]

        if region != "Greenland" and region != "Antarctica" :
            region_tot2.append(1000. * sum_region_vol2/region_area)   # in EWH
            label = "EWH (mm)"
        else:
            region_tot2.append(1000. * sum_region_vol2*1e-12)         # in Gt
            label = "Gt"

        print(decyear2[iepoch],sum_region_vol2/region_area,sum_region_vol2,region_area,series2)

    region_tot2 = np.array(region_tot2)

### plotting #####        
# plot the time series
region_tot = np.array(region_tot)


if plotsize == "large":
    f1 = plt.figure(1,figsize=(12,6))
elif plotsize == "small" :
    fl = plt.figure(1,figsize=(7,5))
    
plt.tick_params(axis='both', which='major', labelsize=16)

plt.plot(decyear,region_tot+offset,"black",linewidth=5.0,label=series1)
    
plt.title(region + " from " + hdf5_file)
plt.xlabel('Year',fontsize=16)
#plt.ylim(ymin,ymax)
plt.ylabel(label,fontsize=16)
plt.grid()

# add the GSFC and JPL time series
if plot_JPL == "Y":
    jpl_region = np.loadtxt("/Mdata/GRACE/software/tables/data/JPL_"+region+".dat",skiprows = 15)
    # the JPL epoch is days since 01/01/2002. Convert to decimal year
    decyear_JPL = np.zeros(len(jpl_region[:,0]))
    for i in range (len(jpl_region[:,0])):
        decyear_JPL[i] = 2002.0 + jpl_region[i,0]/365.0
    
    if region != "Greenland" and region != "Antarctica" :  
        plt.plot(decyear_JPL,jpl_region[:,1]*10,"blue",linewidth=3.0,label="JPL")
    else :
        plt.plot(decyear_JPL,jpl_region[:,1]*10*region_area*1e-12,"blue",linewidth=3.0,label="JPL")
    

if plot_GSFC == "Y" and region != "Antarctica" :
    gsfc_region = np.loadtxt("/Mdata/GRACE/software/tables/data/GSFC_"+region+".dat",skiprows = 13)        
    if region != "Greenland" :    
        plt.plot(gsfc_region[:,0],gsfc_region[:,1]*10,"red",linewidth=1.0,label="GSFC")
    else :
        plt.plot(gsfc_region[:,0],gsfc_region[:,1]*10.*region_area*1e-12,"red",linewidth=1.0,label="GSFC")
    
elif region == "Antarctica" :
    gsfc_EAnt = np.loadtxt("/Mdata/GRACE/software/tables/data/GSFC_EAnt.dat",skiprows = 12)
    gsfc_WAnt = np.loadtxt("/Mdata/GRACE/software/tables/data/GSFC_WAnt.dat",skiprows = 12)
    gsfc_AntPen = np.loadtxt("/Mdata/GRACE/software/tables/data/GSFC_AntPeninsula.dat",skiprows = 12)
    gsfc_tot = np.zeros(len(gsfc_EAnt))
    gsfc_tot = gsfc_EAnt[:,1]+gsfc_WAnt[:,1]+gsfc_AntPen[:,1]
    plt.plot(gsfc_WAnt[:,0],gsfc_tot*region_area*1e-12,"red",linewidth=1.0 )

if lowpass_apr != "none":
    print("plot the lowpass apriori model: ",lowpass_apr)
    plt.plot(decyear,lowpass_tot,"pink",linewidth=2.0)

if series2 != "none" :
    plt.plot(decyear2,region_tot2,"gray",linewidth=2.0,label=series2)

#plt.legend([hdf5_file,'apriori','JPL','GSFC','lowpass'])
#plt.legend(["ANU",'JPL','GSFC',series2])
plt.legend()

plt.show()

