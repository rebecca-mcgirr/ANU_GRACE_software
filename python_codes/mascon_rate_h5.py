### Script to compute the linear rate of change through a set of mascon solutions stacked together in an hdf5 file
#
# P. Tregoning
# 29 July 2022

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import h5py
from gracetools.io.mascons import mascons
import argparse
from sklearn import datasets, linear_model

###############################################################################
## READ ARGUMENTS    
############################################################################### 
# Define arguments   
parser = argparse.ArgumentParser(description=' This program will compute the rate of change for each mascon from a set of solutions in an hdf5 file.')
parser.add_argument("--soln_h5",default='none',type=str,help="input h5 solution file")
parser.add_argument("--mascons_h5",default='none',type=str,help="input mascon file in hdf5 format ")

# Read arguments
args = parser.parse_args()
# Declare arguments
solns_hdf5 = args.soln_h5
mascons_h5= args.mascons_h5

print("Will compute rates of mascons in file ",mascons_h5," using solutions from file ",solns_hdf5)

# read in the mascon solutions from the hdf5 file
print("reading solutions from ",solns_hdf5)
with h5py.File(solns_hdf5, "r") as h5f:
    sol = h5f['solution/ewh'][:,:]
    sigma_ewh = h5f['solution/sigma_ewh'][:,:]
    year = h5f['time/year'][:]
    decyear = h5f['time/decyear'][:]
    month = h5f['time/month'][:]
    density = h5f['time/month'][:]

# get the dimensions of epochs and mascons
n_epochs = len(decyear)
n_mascons = len(sol[1,:])
print("There are ",n_epochs," epochs and ",n_mascons," mascons in file: ",solns_hdf5)


# read in the mascon file
print("Reading in the mascon hdf5 file")
msc = mascons(mascons_h5)
msc.read_mascons_h5()


#####################################################
# loop over the mascons and calculate the linear rate
msc_rates = np.zeros((n_mascons))
model = linear_model.LinearRegression()
x = decyear.reshape((-1,1))

print("Calculating the rates for each mascon ...")
for imsc in range (n_mascons):
    model.fit(x,sol[:,imsc])
    msc_rates[imsc] = model.coef_
    print("%9.4f %9.4f" % (msc.Plat[imsc],msc.Plon[imsc])," %7.2f" % (msc_rates[imsc]*1e3)," mm/yr"," mascon ",imsc)
    
#####################################################


#####################################################
#
#                 Plotting
#
#####################################################
#import tkinter
#import matplotlib
#matplotlib.use('Qt5Agg')

fig = plt.figure(figsize=(16,8))
ax = plt.axes (projection=ccrs.Robinson() )
plt.title('EWH rate (mm/yr)' )
mesh = ax.scatter (msc.Plon, msc.Plat, transform=ccrs.PlateCarree(), c=msc_rates*1000, cmap='gist_rainbow',vmin=-300,vmax=300)
bar = plt.colorbar (mesh, orientation='vertical')
ax.coastlines ()
#plt.savefig("rate.jpg")
plt.show()

    
      




