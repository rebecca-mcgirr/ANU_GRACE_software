# python script to plot kbra data as groundtrack and as time series. It includes a slider
# to permit interactive movement of the dot on the map
#
# P. Tregoning
# 5 July 2022

from __future__ import print_function
from ipywidgets import interact, interactive, fixed, interact_manual
import numpy as np
import ipywidgets as widgets
import os
import sys
import matplotlib.pyplot as plt

def slider_kbra(x):
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    
    # define the colour map
    colour_map = "jet"
    
    # read the data file from this generic name
    kbra_data = np.loadtxt("tmp.data",skiprows = 6)
    
    fig = plt.figure(figsize=(10,5))
    ax = plt.axes(projection=ccrs.Robinson())
    ax.coastlines()

    ymin = -3
    ymax = 3

    # plot the data on the groundtrack
    plt.scatter(kbra_data[:,3],kbra_data[:,2],c=kbra_data[:,1]*1000,s=10,vmin=ymin,vmax=ymax,cmap=colour_map,transform=ccrs.PlateCarree())

    # and a dot at the epoch
    plt.scatter(kbra_data[x-1,3],kbra_data[x-1,2],s=500, c='black', marker='o',transform=ccrs.PlateCarree())

    bar = plt.colorbar(orientation='vertical')
    bar.set_label('KBRA residual (nm/s^2)') #check this
    ax.set_global()

    plt.show()
    
    # can I now add a time series under it?
    fig = plt.figure(figsize=(10,4))
    plt.scatter(kbra_data[:,0],kbra_data[:,1]*1000.,c=kbra_data[:,1]*1000,s=10,vmin=ymin,vmax=ymax,cmap=colour_map)

    # and a dot at the epoch
    plt.scatter(kbra_data[x-1,0],kbra_data[x-1,1]*1000.,s=500, c='black', marker='o')

    #plt.show()
    
    return x




# copy the kbra file to a generic name to be read in the function
infile = sys.argv[1]
cmd = "cp -f "+infile+" tmp.data"
os.system(cmd)

# call the function that does everything else
interact(slider_kbra,x=widgets.IntSlider(min=1, max=17257, step=1, value=17257,continuous_update=False))

plt.show()


