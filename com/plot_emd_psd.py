
# coding: utf-8

# In[7]:


# coding: utf-8

# In[9]:

import numpy as np
import matplotlib
import sys
import scipy.fftpack as fft 

import tkinter #python 3 syntax
root = tkinter.Tk()
root.withdraw()


#matplotlib.use('TkAgg')  
import matplotlib.pyplot as plt
#plt.ion()

print( "Runstring: python ~/ga/com/plot_emd_decomp.py <input file>)")

# PT180706: set the window size based on the screen size. Make it 80% of the height and 90% of the width
width, height = root.winfo_screenwidth(), root.winfo_screenheight()


# In[3]:


# In[13]:

# read in the GRACE A accelerometer data
print ("reading in the data")
emd_decomp = np.loadtxt(sys.argv[1])
#emd_decomp = np.loadtxt('emd.decomp2')


nepochs = len(emd_decomp[:,0])


# In[26]:


# In[ ]:

# plot the decomposed accelerometer values, one component per figure box
print("Plotting components separately..")
fig, axs = plt.subplots(nrows=13, ncols = 2, sharex = True, figsize=( int(width*0.012), int(height*0.011)) )

for i in range(13):
	axs[i,0].grid(color='grey')
	axs[i,0].plot(emd_decomp[:,0],emd_decomp[:,i+1]*1.e9,'red')
	axs[i,0].set_ylabel(' C%i'%(i+1))
axs[-1,0].set_xlabel('Seconds of day')


# periods = np.zeros((86399,13))
# power = np.zeros((86399,13))

xhat =  fft.fft(emd_decomp[:,1:], axis=0)
N = xhat.shape[0]
print(xhat.shape)
dt = 1
frequencies = fft.fftfreq(N,dt)
xhat = xhat[1:N//2] / np.sqrt(N)
frequencies = frequencies[1:N//2] 
periods = 1.0/frequencies
periods=periods[::-1]
power= np.sqrt(abs(xhat * xhat))
power=power[::-1,:]
print(power.shape)

pp = np.cumsum(power, axis=0)
print(power[:,2])
print(pp[:,2])
for i in range(9):
	axs[i,1].grid(color='grey')
	mpower = power[:,i].max()
	axs[i,1].plot(periods,power[:,i]/mpower,'k')
	mpp = pp[:,i].max()
	pp[:,i] = pp[:,i]/mpp
	axs[i,1].plot(periods,pp[:,i] ,'red')
	axs[i,1].set_ylabel(' C%i'%(i+1))
	print (i, periods[np.where(pp[:,i]>0.1)[0][0]])
axs[-1,1].set_xlabel('Period')




# for i in range(13):
#     xhat =  fft.fft(emd_decomp[:,i])
#     N = len(xhat)
#     dt = 1

#     frequencies = fft.fftfreq(N,dt)

#     # We are interested into the positive frequencies.
#     xhat = xhat[1:N//2] / np.sqrt(N)   #normalisation
#     frequencies = frequencies[1:N//2]   
#     periods[0:len(frequencies),i] = 1.0/frequencies

#     # The sign “//” in N//2 just return the integer part of the division
#     power[0:len(xhat),i] = np.sqrt(abs(xhat * xhat))

#     print('number of elements in periods and power',i,len(periods),len(power))
    





# plt.subplot(13,2,3)
# plt.grid(color='grey')
# plt.plot(emd_decomp[:,0],emd_decomp[:,2]*1.e9,'red')
# plt.ylabel('C2 ',fontsize=15)
# plt.xlabel('Seconds of day',fontsize=15)

# plt.subplot(13,2,5)
# plt.grid(color='grey')
# plt.plot(emd_decomp[:,0],emd_decomp[:,3]*1.e9,'red')
# plt.ylabel('C3 ',fontsize=15)
# plt.xlabel('Seconds of day',fontsize=15)

# plt.subplot(13,2,7)
# plt.grid(color='grey')
# plt.plot(emd_decomp[:,0],emd_decomp[:,4]*1.e9,'red')
# plt.ylabel('C4 ',fontsize=15)
# plt.xlabel('Seconds of day',fontsize=15)

# plt.subplot(13,2,9)
# plt.grid(color='grey')
# plt.plot(emd_decomp[:,0],emd_decomp[:,5]*1.e9,'red')
# plt.ylabel('C5 ',fontsize=15)
# plt.xlabel('Seconds of day',fontsize=15)

# plt.subplot(13,2,11)
# plt.grid(color='grey')
# plt.plot(emd_decomp[:,0],emd_decomp[:,6]*1.e9,'red')
# plt.ylabel('C6 ',fontsize=15)
# plt.xlabel('Seconds of day',fontsize=15)

# plt.subplot(13,2,13)
# plt.grid(color='grey')
# plt.plot(emd_decomp[:,0],emd_decomp[:,7]*1.e9,'red')
# plt.ylabel('C7 ',fontsize=15)
# plt.xlabel('Seconds of day',fontsize=15)

# plt.subplot(13,2,15)
# plt.grid(color='grey')
# plt.plot(emd_decomp[:,0],emd_decomp[:,8]*1.e9,'red')
# plt.ylabel('C8 ',fontsize=15)
# plt.xlabel('Seconds of day',fontsize=15)

# plt.subplot(13,2,17)
# plt.grid(color='grey')
# plt.plot(emd_decomp[:,0],emd_decomp[:,9]*1.e9,'red')
# plt.ylabel('C9 ',fontsize=15)
# plt.xlabel('Seconds of day',fontsize=15)

# plt.subplot(13,2,19)
# plt.grid(color='grey')
# plt.plot(emd_decomp[:,0],emd_decomp[:,10]*1.e9,'red')
# plt.ylabel('C10',fontsize=15)
# plt.xlabel('Seconds of day',fontsize=15)

# plt.subplot(13,2,21)
# plt.grid(color='grey')
# plt.plot(emd_decomp[:,0],emd_decomp[:,11]*1.e9,'red')
# plt.ylabel('C11',fontsize=15)
# plt.xlabel('Seconds of day',fontsize=15)

# plt.subplot(13,2,23)
# plt.grid(color='grey')
# plt.plot(emd_decomp[:,0],emd_decomp[:,12]*1.e9,'red')
# plt.ylabel('C12',fontsize=15)
# plt.xlabel('Seconds of day',fontsize=15)

# plt.subplot(13,2,25)
# plt.grid(color='grey')
# plt.plot(emd_decomp[:,0],emd_decomp[:,13]*1.e9,'red')
# plt.ylabel('C13',fontsize=15)
# plt.xlabel('Seconds of day',fontsize=15)

# #plt.show()


# # In[23]:


# # now calculate and plot the PSD for each component and plot it alongside the time series
# # perform the FFT for the analytical kbra
# periods = np.zeros((86399,13))
# power = np.zeros((86399,13))

# for i in range(13):
#     xhat =  fft.fft(emd_decomp[:,i])
#     N = len(xhat)
#     dt = 1

#     frequencies = fft.fftfreq(N,dt)

#     # We are interested into the positive frequencies.
#     xhat = xhat[1:N//2] / np.sqrt(N)   #normalisation
#     frequencies = frequencies[1:N//2]   
#     periods[0:len(frequencies),i] = 1.0/frequencies

#     # The sign “//” in N//2 just return the integer part of the division
#     power[0:len(xhat),i] = np.sqrt(abs(xhat * xhat))

#     print('number of elements in periods and power',i,len(periods),len(power))
    


# # In[27]:


# ######### plot the PSDs for each component
# plt.subplot(13,2,2)
# plt.grid(color='grey')
# plt.plot(periods[0:len(frequencies),0],power[0:len(frequencies),0]*1.e9,'red')
# plt.ylabel(' C1',fontsize=15)
# plt.xlabel('Period',fontsize=15)

# plt.subplot(13,2,4)
# plt.grid(color='grey')
# plt.plot(periods[0:len(frequencies),1],power[0:len(frequencies),1]*1.e9,'red')
# plt.ylabel(' C2',fontsize=15)
# plt.xlabel('Period',fontsize=15)

# plt.subplot(13,2,6)
# plt.grid(color='grey')
# plt.plot(periods[0:len(frequencies),2],power[0:len(frequencies),2]*1.e9,'red')
# plt.plot(periods[:,2],power[:,2]*1.e9,'red')
# plt.ylabel(' C3',fontsize=15)
# plt.xlabel('Period',fontsize=15)

# plt.subplot(13,2,8)
# plt.grid(color='grey')
# plt.plot(periods[0:len(frequencies),3],power[0:len(frequencies),3]*1.e9,'red')
# plt.ylabel(' C4',fontsize=15)
# plt.xlabel('Period',fontsize=15)

# plt.subplot(13,2,10)
# plt.grid(color='grey')
# plt.plot(periods[0:len(frequencies),4],power[0:len(frequencies),4]*1.e9,'red')
# plt.ylabel(' C5',fontsize=15)
# plt.xlabel('Period',fontsize=15)

# plt.subplot(13,2,12)
# plt.grid(color='grey')
# plt.plot(periods[0:len(frequencies),5],power[0:len(frequencies),5]*1.e9,'red')
# plt.ylabel(' C6',fontsize=15)
# plt.xlabel('Period',fontsize=15)

# plt.subplot(13,2,14)
# plt.grid(color='grey')
# plt.plot(periods[0:len(frequencies),6],power[0:len(frequencies),6]*1.e9,'red')
# plt.ylabel(' C7',fontsize=15)
# plt.xlabel('Period',fontsize=15)

# plt.subplot(13,2,16)
# plt.grid(color='grey')
# plt.plot(periods[0:len(frequencies),7],power[0:len(frequencies),7]*1.e9,'red')
# plt.ylabel(' C8',fontsize=15)
# plt.xlabel('Period',fontsize=15)

# plt.subplot(13,2,18)
# plt.grid(color='grey')
# plt.plot(periods[0:len(frequencies),8],power[0:len(frequencies),8]*1.e9,'red')
# plt.ylabel(' C9',fontsize=15)
# plt.xlabel('Period',fontsize=15)

# plt.subplot(13,2,20)
# plt.grid(color='grey')
# plt.plot(periods[0:len(frequencies),9],power[0:len(frequencies),9]*1.e9,'red')
# plt.ylabel(' C10',fontsize=15)
# plt.xlabel('Period',fontsize=15)

# plt.subplot(13,2,22)
# plt.grid(color='grey')
# plt.plot(periods[0:len(frequencies),10],power[0:len(frequencies),10]*1.e9,'red')
# plt.ylabel(' C11',fontsize=15)
# plt.xlabel('Period',fontsize=15)

# plt.subplot(13,2,24)
# plt.grid(color='grey')
# plt.plot(periods[0:len(frequencies),11],power[0:len(frequencies),11]*1.e9,'red')
# plt.ylabel(' C12',fontsize=15)
# plt.xlabel('Period',fontsize=15)

# plt.subplot(13,2,26)
# plt.grid(color='grey')
# plt.plot(periods[0:len(frequencies),12],power[0:len(frequencies),12]*1.e9,'red')
# plt.ylabel(' C13',fontsize=15)
# plt.xlabel('Period',fontsize=15)

plt.show()

