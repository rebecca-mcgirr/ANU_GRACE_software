#!/usr/bin/env python3
import os
import sys
import argparse
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from gracetools.io.mascons import *

###############################################################################
## READ ARGUMENTS    
############################################################################### 
debug = False
# Define arguments
parser = argparse.ArgumentParser(description='This program will integrate EWH ...')
parser.add_argument("hdf5_file",default=None,type=str,help="ANU hdf5 file of monthly solutions")

# Read arguments
args = parser.parse_args()

# Declare arguments
hdf5_file = args.hdf5_file

###############################################################################
## LOAD FILES    
############################################################################### 
# Load ANU solution files
sANU = ANU(fname=hdf5_file)
sANU.line = {'col':'k','width':1.75,'style':'-','label':'ANU'} 
load_mask(sANU,"/Mdata/GRACE/software/tables/data/mascons_stage5_V006_200km_region_mask.h5")

sANU2 = ANU(fname="/Mdata/GRACE/software/tables/data/b8i3_Seb_corrected.h5")
sANU2.line = {'col':'grey','width':1.5,'style':'-','label':'ANU b8i3'} 
load_mask(sANU2,"/Mdata/GRACE/software/tables/data/mascons_stage5_V006_200km_region_mask.h5")

# Load GSFC files
sGSFC = GSFC(fname="/Mdata/GRACE/software/tables/data/gsfc.glb_.200204_202112_rl06v2.0_obp-ice6gd.h5")
sGSFC.ewh = sGSFC.ewh.T
load_mask(sGSFC,"/Mdata/GRACE/software/tables/data/GSFC_mascon_region_mask.h5")
sGSFC.line = {'col':'red','width':1,'style':'--','label':'GSFC'}

# Load JPL files
sJPL = JPL(fname="/Mdata/GRACE/software/tables/data/GRCTellus.JPL.200204_202207.GLO.RL06M.MSCNv02CRI.nc")
load_mask(sJPL,"/Mdata/GRACE/software/tables/data/JPL_mascon_areas_region_mask.nc")
sJPL.line = {'col':'blue','width':1,'style':'--','label':'JPL'}

# Load CSR files 
sCSR = CSR(fname="/Mdata/GRACE/software/tables/data/CSR_GRACE_GRACE-FO_RL06_Mascons_all-corrections_v02.nc")
load_mask(sCSR,"/Mdata/GRACE/software/tables/data/CSR_mascon_areas_region_mask.nc")
sCSR.line = {'col':'green','width':1,'style':'--','label':'CSR'}

###############################################################################
## INTEGRATE AND PLOT RESULTS
###############################################################################
print("Integrating land and ocean mascons")
series = [sGSFC, sJPL, sCSR, sANU2, sANU]
#series = [sANU, sANU3, sCSR]
fig, axs = plt.subplots(3,1,figsize=(8,13))
demean = False
if np.min(sANU.date) <= 2004 and np.max(sANU.date) >= 2010: demean = True
for si in series:
    si.land = intg_ewh(si,si.land_mask,'GT',demean=demean)
    si.ocean = intg_ewh(si,si.ocean_mask,'GT',demean=demean)
    si.all = intg_ewh(si,np.ones(si.land_mask.shape),'GT',demean=demean)      
    axs[0].plot(si.date, si.land, color=si.line['col'], label=si.line['label'], linewidth=si.line['width'], linestyle=si.line['style'])
    axs[0].set_title('land')
    axs[1].plot(si.date, si.ocean, color=si.line['col'], label=si.line['label'], linewidth=si.line['width'], linestyle=si.line['style'])
    axs[1].set_title('ocean')
    axs[2].plot(si.date, si.all, color=si.line['col'], label=si.line['label'], linewidth=si.line['width'], linestyle=si.line['style'])
    axs[2].set_title('all')
for ax in axs: 
    ax.grid()
    ax.set_ylabel('GT')
handles, labels = plt.gca().get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=len(series))
fig.suptitle(hdf5_file.split("/")[-1], fontsize=14)
plt.show()
