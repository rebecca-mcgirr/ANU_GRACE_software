#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import h5py
import sys

# In[14]:


# read the information from the hdf5 file
#hdf5_file = "2021_corrected.h5"
hdf5_file = sys.argv[1]
with h5py.File(hdf5_file, "r") as h5f:
    sol = h5f['solution/ewh'][:]
    sigma_ewh = h5f['solution/sigma_ewh'][:]
    year = h5f['time/year'][:]
    decyear = h5f['time/decyear'][:]
    month = h5f['time/month'][:]


# In[49]:


# get the dimensions of epochs and mascons
n_epochs = len(decyear)
n_mascons = len(sol[1,:])
print("There are ",n_epochs," epochs and ",n_mascons," mascons in file: ",hdf5_file)


# In[51]:


# loop through each solution and write it out to an ascii file in .vcv format
for iepoch in range(n_epochs):
    # make the vcv file name based on year/month
    file = "h5_" + str(year[iepoch]) + "_" + str("%02i" % month[iepoch]) + ".vcv"
    f = open(file, "w")
    print ("Outputting file ",iepoch+1,": ",file)
    
    # write some header information
    header = "VCV file written from hdf5 file:" + hdf5_file + "\n"
    f.write(header)
    f.write(str(n_mascons))
    f.write("\n\n\n SOLUTION A PRIORI AND VECTOR:\n")
    f.write(" PARAMETER                     A PRIORI             VECTOR            SIGMA\n")
    # now, loop through the mascons for this epoch and write out the mascon solution and uncertainty
    #    25. MC00001 (m)                   -0.0348000       -0.0293975        0.0148639

    for imascon in range (n_mascons):
        msc = "%05i" % (imascon+1)
        f.write("{0:5.0f}. MC{1:5s} (m) {2:28.7f}{3:17.7f}{4:17.7f}\n".format(imascon+25,msc,sol[iepoch,imascon],sol[iepoch,imascon],sigma_ewh[imascon,iepoch]))
    f.close()
    


# In[ ]:





# In[ ]:




