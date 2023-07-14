
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib 
import sys

matplotlib.use('TkAgg')  
import matplotlib.pyplot as plt

import matplotlib.patches as patches

print( "Runstring: python ~/ga/com/plot_sca.py SCA1B_2008-08-12.asc")


# In[ ]:


# read in the data from a star camera file
data = np.genfromtxt(sys.argv[1], skip_header=24, autostrip = True, usecols=(0,2,3,4,5,6))

# create different arrays of the data
nvals = len(data[:,0])
gracesec = data[:,0]
quat = data[:,2:6]
camera = data[:,1]

# make a time index that starts with the first epoch being 1
epochs = data[:,0] - data[0,0]

quat.shape


# In[16]:


# plot the data
plt.figure(num=5, figsize=(15,10))

# plot the four elements of the quaternion
plt.subplot(511)
plt.plot(epochs,quat[:,0],'ro',markersize=1)
plt.ylabel('Quat 1',fontsize=15)
plt.xlabel('epoch',fontsize=15)

plt.subplot(512)
plt.plot(epochs,quat[:,1],'ro',markersize=1)
plt.ylabel('Quat 2',fontsize=15)
plt.xlabel('epoch',fontsize=15)

plt.subplot(513)
plt.plot(epochs,quat[:,2],'ro',markersize=1)
plt.ylabel('Quat 3',fontsize=15)
plt.xlabel('epoch',fontsize=15)

plt.subplot(514)
plt.plot(epochs,quat[:,3],'ro',markersize=1)
plt.ylabel('Quat 4',fontsize=15)
plt.xlabel('epoch',fontsize=15)







plt.subplot(515)
plt.plot(epochs,camera,'ro')
plt.ylabel('Camera info',fontsize=15)
plt.xlabel('epoch',fontsize=15)
plt.ylim(0,5)

#plt.ion()
plt.show()

