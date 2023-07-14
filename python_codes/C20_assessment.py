#!/usr/bin/env python
# coding: utf-8

# In[3]:


from gracetools.io.fit import Fit
from gracetools.io.aod import AOD1b
from gracetools.io.mascons import mascons
import numpy as np
from scipy.interpolate import interp1d
import pyshtools as sht
import datetime
from pathlib import Path
import h5py
import concurrent.futures
# PT220307: import the library of functions to extend time series
import extend_lib as ex


# In[4]:


def toYearFraction(date):
    import time
    from datetime import datetime as dt
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction


# In[5]:


if __name__ == "__main__":

    print('python3 ~/ga/python_codes/C20_assessment.py --mascons /Mdata/GRACE/software/tables/data/mascons_stage5_V006_200km.h5 --h5 h5file.h5')    

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--mascons', type=str, help='mascons files')
    parser.add_argument('--h5', type=str, help=" hdf5 file containing solutions")
    #parser.add_argument('--love', type=str,default = "/home/geod/pault/bin/Load_Love2_CM.dat", help='love number file')
    args = parser.parse_args()


    print(f"reading {args.mascons}")
    mANU = mascons(args.mascons)
    mANU.read_mascons_h5()
    mANU.get_mapping_sht()


# In[6]:



print("input solution hdf5 file: ",args.h5)
with h5py.File(args.h5, "r") as h5f:
    apriori_ewh = h5f['solution/prefit'][:]
    # this is our EWH solution, with deg1, AOD1B, deg2 and GIA corrections all applied
    sol = h5f['solution/ewh'][:]
    sigma_ewh = h5f['solution/sigma_ewh'][:]
    # this is the C20 anomaly on each mascon wrt the TN14 mean value
    C20_mascon_stokes_anomaly = h5f['correction/C20_mascon_stokes_anomaly'][:]
    C20_mascon_EWH_anomaly = h5f['solution/C20_mascon_EWH_anomaly'][:]

    year = h5f['time/year'][:]
    decyear = h5f['time/decyear'][:]
    month = h5f['time/month'][:]
    # this "TN14_corr" is to account for the difference in the mean gravity field C20 value, plus the permanent tide added to ours
    TN14_corr = h5f['correction/C20_TN14_corr'][:]
    
n_epochs = len(decyear)
n_mascons = len(sol[1,:])
print("Number of epochs: ",n_epochs," Number of mascons: ",n_mascons)

## read in the TN14 C20 and C30 data
TN14 = np.loadtxt('/Mdata/GRACE/software/tables/data/TN-14_C30_C20_GSFC_SLR.txt', skiprows=37)
dec_yr_TN14      = TN14[:,1]
C20_anomaly_TN14 = TN14[:,3]  
C30_anomaly_TN14 = TN14[114:,5] 

# In[7]:



#degmax=180
#n, h, l, k = np.loadtxt(args.love, unpack=True)
#n = n[:degmax + 1]
#h = h[:degmax + 1]
#k = k[:degmax + 1]
#l = l[:degmax + 1]
#k[0] = 0.0
#k[1] = 0.026
## PT220908: change Earth radius from 6371.0e3 to 6378.1363e3
#Conv = 1 / 3.0 * 6378.1363e3 *5.515 * (2 * n + 1) / (1 + k)

# PT220915: convert the mascon solutions into a degree-3 spherical harmonic series expansion, using different subsets of the 
#           regions of Earth: all, ocean, Antarctica, Greenland, Amazon, all_ice  etc

#regions = (["all","TN14_corr","ocean","Ant","Grn","Amazon","Eurasia","NAmerica","Africa"])
regions = (["all","TN14_corr","Grn"])
regions = ["all","Ant","Grn"]
C20_mascon = np.zeros( (len(sol[:,0]),20) )
C20_TN14_corr = np.zeros(len(sol[:,0]))
EWH_tmp = np.zeros(n_mascons)
stokes_tmp = np.zeros(n_mascons)
msc_mask = np.zeros(n_mascons)
#for iregion in range (2):   
for iregion in range (len(regions)):   
    print("Calculating C20 for region: ",regions[iregion])
             
    for iepoch, (dec, y, m) in enumerate(zip(decyear, year, month)):
        # need to set to non-zero only the mascons within the requested region
        d = 0.
        if regions[iregion] == "all" :
            d = C20_mascon_stokes_anomaly[iepoch, :] 
            msc_sph = sht.expand.SHExpandDH(d[mANU.idx_near], sampling=2, csphase=1,norm=1, lmax_calc=2)
            C20_mascon[iepoch,10] = msc_sph[0,2,0]  
            #print(decyear[iepoch]," all",C20_mascon[iepoch,0] )               

        if regions[iregion] == "ocean" :
            stokes_ocean = np.zeros(n_mascons)
            for imsc in range (n_mascons):
                if imsc < 64 :
                    print(imsc,mANU.Pdesc_label[imsc])
                if mANU.Pdensity[imsc] > 1001. :
                    #print(iregion,iepoch,imsc,sol[iepoch,imsc],mANU.Pdensity[imsc])
                    stokes_ocean[imsc] = C20_mascon_stokes_anomaly[iepoch,imsc] 
                                        
            msc_sph = sht.expand.SHExpandDH(stokes_ocean[mANU.idx_near], sampling=2, csphase=1,norm=1, lmax_calc=2)
            C20_mascon[iepoch,11] = msc_sph[0,2,0]
            #print(decyear[iepoch]," ocean",C20_mascon[iepoch,1]    )            

        if regions[iregion] == "Ant" :
            stokes_Ant = np.zeros(n_mascons)
            for imsc in range (n_mascons):     # Antarctica is label 57
                if mANU.Pdesc[imsc] == 57 :
                    #print(iregion,iepoch,imsc,sol[iepoch,imsc],mANU.Pdensity[imsc])
                    EWH_tmp[imsc] = sol[iepoch,imsc] - TN14_corr[iepoch,imsc]
                    stokes_Ant[imsc] = C20_mascon_stokes_anomaly[iepoch,imsc] 
                                        
            msc_sph = sht.expand.SHExpandDH(stokes_Ant[mANU.idx_near], sampling=2, csphase=1,norm=1, lmax_calc=2)
            C20_mascon[iepoch,12] = np.float(msc_sph[0,2,0])
            #print(decyear[iepoch]," Ant",C20_mascon[iepoch,1]    )            

        if regions[iregion] == "Grn" :
            stokes_Grn = np.zeros(n_mascons)
            for imsc in range (n_mascons):
                #print(iregion,regions[iregion],n_mascons,imsc,mANU.Pdesc[imsc])
                if mANU.Pdesc[imsc] == 56 :   # Greenland is label 56
                    #print(iregion,iepoch,imsc,C20_mascon_stokes_anomaly[iepoch,imsc],mANU.Pdensity[imsc])
                    stokes_Grn[imsc] = C20_mascon_EWH_anomaly[iepoch,imsc]
                                        
            msc_sph = sht.expand.SHExpandDH(stokes_Grn[mANU.idx_near], sampling=2, csphase=1,norm=1, lmax_calc=2)
            C20_mascon[iepoch,13] = np.float(msc_sph[0,2,0])
            #print(decyear[iepoch]," Grn",C20_mascon[iepoch,1]    )            

        if regions[iregion] == "Amazon" :
            stokes_Amazon = np.zeros(n_mascons)
            for imsc in range (n_mascons):
                if mANU.Pdesc[imsc] == 63 :    # Amazon is label 63
                    if iepoch == 1:
                        print(iregion,iepoch,imsc,sol[iepoch,imsc],C20_mascon_stokes_anomaly[iepoch,imsc],mANU.Pdensity[imsc])
                    stokes_Amazon[imsc] = C20_mascon_stokes_anomaly[iepoch,imsc]
                                        
            msc_sph = sht.expand.SHExpandDH(stokes_Amazon[mANU.idx_near], sampling=2, csphase=1,norm=1, lmax_calc=2)
            C20_mascon[iepoch,14] = np.float(msc_sph[0,2,0])
            print(iepoch,"Grn Amazon C20:",C20_mascon[iepoch,12],C20_mascon[iepoch,13],C20_mascon[iepoch,14])
            
        if regions[iregion] == "Eurasia" :
            for imsc in range (n_mascons):
                if mANU.Pdesc[imsc] == 63 :    # Amazon is label 63
                    #print(iregion,iepoch,imsc,sol[iepoch,imsc],mANU.Pdensity[imsc])
                    EWH_tmp[imsc] = sol[iepoch,imsc] - TN14_corr[iepoch,imsc]
                    EWH_tmp2[imsc] = ourC20[iepoch,imsc]
                                        
            msc_sph = sht.expand.SHExpandDH(EWH_tmp[mANU.idx_near], sampling=2, csphase=1,norm=1, lmax_calc=3)
            C20_mascon[iepoch,5] = msc_sph[0,2,0] - TN14_corr[iepoch,imsc]
            msc_sph = sht.expand.SHExpandDH(EWH_tmp2[mANU.idx_near], sampling=2, csphase=1,norm=1, lmax_calc=3)
            C20_mascon[iepoch,15] = msc_sph[0,2,0]

        if regions[iregion] == "NAmerica" :
            stokes_NAm = np.zeros(n_mascons)
            for imsc in range (n_mascons):
                if mANU.Pdesc[imsc] == 60 :    # NAmerica is label 63
                    #print(iregion,iepoch,imsc,sol[iepoch,imsc],mANU.Pdensity[imsc])
                    stokes_NAm[imsc] = C20_mascon_stokes_anomaly[iepoch,imsc]
                                        
            msc_sph = sht.expand.SHExpandDH(stokes_NAm[mANU.idx_near], sampling=2, csphase=1,norm=1, lmax_calc=2)
            C20_mascon[iepoch,16] = msc_sph[0,2,0]

        if regions[iregion] == "Africa" :
            for imsc in range (n_mascons):
                if mANU.Pdesc[imsc] == 61 :    # Africa is 61
                    #print(iregion,iepoch,imsc,sol[iepoch,imsc],mANU.Pdensity[imsc])
                    EWH_tmp[imsc] = sol[iepoch,imsc] - TN14_corr[iepoch,imsc]
                    EWH_tmp2[imsc] = C20_mascon_stokes_anomaly[iepoch,imsc]
                                        
            msc_sph = sht.expand.SHExpandDH(EWH_tmp[mANU.idx_near], sampling=2, csphase=1,norm=1, lmax_calc=3)
            C20_mascon[iepoch,7] = msc_sph[0,2,0] - TN14_corr[iepoch,imsc]
            msc_sph = sht.expand.SHExpandDH(EWH_tmp2[mANU.idx_near], sampling=2, csphase=1,norm=1, lmax_calc=3)
            C20_mascon[iepoch,17] = msc_sph[0,2,0]

        if regions[iregion] == "TN14_corr" :
            d = TN14_corr[iepoch, :]
            msc_sph = sht.expand.SHExpandDH(d[mANU.idx_near], sampling=2, csphase=1,norm=1, lmax_calc=3)
            C20_TN14_corr[iepoch] = msc_sph[0,2,0] 

# now plot the various C20 time series
import matplotlib.pyplot as plt

f1 = plt.figure(1,figsize=(14,7))

# conversion factor for C20 from EWH to Stoke's dimensionless
Conv_C20 = 1 / 3.0 * 6378.1363e3 *5.515 * (2 * 2 + 1) / (1 + 0.026)

# TN14
plt.plot(dec_yr_TN14,C20_anomaly_TN14,"black",linewidth=3.0,label="TN14")

# the whole Earth
#plt.plot(decyear,C20_mascon[:,0]*1e10/Conv_C20,"red",linewidth=1.0)
plt.plot(decyear,C20_mascon[:,10]*1e10,"red",linewidth=3.0,label="ANU")

plt.title("C20 TN14 vs ANU " )
plt.xlabel('Year')
#plt.ylim(-0.06,0.02)
plt.ylabel("C20 (x 1e-10)")
plt.legend(["TN14","b8i2e_200v4 uncorrected","ourC20"])

#plt.show()


f2 = plt.figure(2,figsize=(14,7))
plt.title("C20 by region " )
plt.xlabel('Year')
plt.ylabel("C20 (x 1e-10)")


# Our C20 estimates, uncorrected to TN14 (but other corrections applied)
#plt.plot(decyear,C20_mascon[:,11]*1e10,"darkblue",label="ocean")
plt.plot(decyear,C20_mascon[:,12]*1e10,"plum",label="Ant")
plt.plot(decyear,C20_mascon[:,13]*1e10/Conv_C20,"darkgreen",label="Grn")
#plt.plot(decyear,C20_mascon[:,14]*1e10,"darkorange",label="Amazon")
#plt.plot(decyear,C20_mascon[:,15]*1e2,"darkgoldenrod",label="Eurasia",linewidth=3.0)
#plt.plot(decyear,C20_mascon[:,16]*1e10,"peru",label="NAm")
#plt.plot(decyear,C20_mascon[:,17]*1e2,"salmon",label="Africa")


#plt.plot(decyear,C20_mascon[:,8]*1e10/Conv_C20,"darkblue",label="ocean ourC20",linewidth=3.0)
#plt.plot(decyear,C20_mascon[:,9]*1e10/Conv_C20,"darkgreen",label="Grn ourC20",linewidth=3.0)

# plot the sum
#plt.plot(decyear,(C20_mascon[:,1]+C20_mascon[:,2]+C20_mascon[:,3]+C20_mascon[:,4]+C20_mascon[:,5]+C20_mascon[:,6]+C20_mascon[:,7])*1e10/Conv_C20,"red",label="Sum",linewidth = 6.0)

#plt.legend(["Ocean","Antarctica","Greenland","Amazon"])
plt.legend()
plt.show()


# In[ ]:




