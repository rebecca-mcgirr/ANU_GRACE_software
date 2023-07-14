
# coding: utf-8

# In[9]:

import numpy as np
import matplotlib
import sys
import tkinter #python 3 syntax
root = tkinter.Tk()
root.withdraw()


#matplotlib.use('TkAgg')  
import matplotlib.pyplot as plt
#plt.ion()

print( "Runstring: python ~/ga/com/plot_emd_decomp.py <input file>)")

# PT180706: set the window size based on the screen size. Make it 80% of the height and 90% of the width
width, height = root.winfo_screenwidth(), root.winfo_screenheight()

# In[13]:

# read in the GRACE A accelerometer data
print ("reading in the data")
emd_decomp = np.loadtxt(sys.argv[1])


nepochs = len(emd_decomp[:,0])


# In[ ]:

# plot the decomposed accelerometer values, one component per figure box
print("Plotting components separately..")
plt.figure(26,figsize=( int(width*0.012), int(height*0.011)) )

plt.subplot(13,2,1)
plt.grid(color='grey')
plt.plot(emd_decomp[:,0],emd_decomp[:,1]*1.e9,'red')
plt.ylabel(' C1 (nm/s$^2$)',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,3)
plt.grid(color='grey')
plt.plot(emd_decomp[:,0],emd_decomp[:,2]*1.e9,'red')
plt.ylabel('C2 ',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,5)
plt.grid(color='grey')
plt.plot(emd_decomp[:,0],emd_decomp[:,3]*1.e9,'red')
plt.ylabel('C3 ',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,7)
plt.grid(color='grey')
plt.plot(emd_decomp[:,0],emd_decomp[:,4]*1.e9,'red')
plt.ylabel('C4 ',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,9)
plt.grid(color='grey')
plt.plot(emd_decomp[:,0],emd_decomp[:,5]*1.e9,'red')
plt.ylabel('C5 ',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,11)
plt.grid(color='grey')
plt.plot(emd_decomp[:,0],emd_decomp[:,6]*1.e9,'red')
plt.ylabel('C6 ',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,13)
plt.grid(color='grey')
plt.plot(emd_decomp[:,0],emd_decomp[:,7]*1.e9,'red')
plt.ylabel('C7 ',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,15)
plt.grid(color='grey')
plt.plot(emd_decomp[:,0],emd_decomp[:,8]*1.e9,'red')
plt.ylabel('C8 ',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,17)
plt.grid(color='grey')
plt.plot(emd_decomp[:,0],emd_decomp[:,9]*1.e9,'red')
plt.ylabel('C9 ',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,19)
plt.grid(color='grey')
plt.plot(emd_decomp[:,0],emd_decomp[:,10]*1.e9,'red')
plt.ylabel('C10',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,21)
plt.grid(color='grey')
plt.plot(emd_decomp[:,0],emd_decomp[:,11]*1.e9,'red')
plt.ylabel('C11',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,23)
plt.grid(color='grey')
plt.plot(emd_decomp[:,0],emd_decomp[:,12]*1.e9,'red')
plt.ylabel('C12',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,25)
plt.grid(color='grey')
plt.plot(emd_decomp[:,0],emd_decomp[:,13]*1.e9,'red')
plt.ylabel('C13',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)

######### add them up, from high-frequency to low-frequency
plt.subplot(13,2,2)
plt.grid(color='grey')
summed = emd_decomp[:,1]
plt.plot(emd_decomp[:,0],summed*1.e9,'red')
plt.ylabel(' to C1 (nm/s$^2$)',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,4)
plt.grid(color='grey')
summed = summed + emd_decomp[:,2]
plt.plot(emd_decomp[:,0],summed*1.e9,'red')
plt.ylabel('to C2 ',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,6)
plt.grid(color='grey')
summed = summed + emd_decomp[:,3]
plt.plot(emd_decomp[:,0],summed*1.e9,'red')
plt.ylabel('to C3 ',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,8)
plt.grid(color='grey')
summed = summed + emd_decomp[:,4]
plt.plot(emd_decomp[:,0],summed*1.e9,'red')
plt.ylabel('to C4 ',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,10)
plt.grid(color='grey')
summed = summed + emd_decomp[:,5]
plt.plot(emd_decomp[:,0],summed*1.e9,'red')
plt.ylabel('to C5 ',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,12)
plt.grid(color='grey')
summed = summed + emd_decomp[:,6]
plt.plot(emd_decomp[:,0],summed*1.e9,'red')
plt.ylabel('to C6 ',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,14)
plt.grid(color='grey')
summed = summed + emd_decomp[:,7]
plt.plot(emd_decomp[:,0],summed*1.e9,'red')
plt.ylabel('to C7 ',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,16)
plt.grid(color='grey')
summed = summed + emd_decomp[:,8]
plt.plot(emd_decomp[:,0],summed*1.e9,'red')
plt.ylabel('to C8 ',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,18)
plt.grid(color='grey')
summed = summed + emd_decomp[:,9]
plt.plot(emd_decomp[:,0],summed*1.e9,'red')
plt.ylabel('to C9 ',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,20)
plt.grid(color='grey')
summed = summed + emd_decomp[:,10]
plt.plot(emd_decomp[:,0],summed*1.e9,'red')
plt.ylabel('to C10',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,22)
plt.grid(color='grey')
summed = summed + emd_decomp[:,11]
plt.plot(emd_decomp[:,0],summed*1.e9,'red')
plt.ylabel('to C11',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,24)
plt.grid(color='grey')
summed = summed + emd_decomp[:,12]
plt.plot(emd_decomp[:,0],summed*1.e9,'red')
plt.ylabel('to C12',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,26)
plt.grid(color='grey')
summed = summed + emd_decomp[:,13]
plt.plot(emd_decomp[:,0],summed*1.e9,'red')
plt.ylabel('to C13',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)


# add them up, from lowest frequency through to highest, stopping the addition at component 6
plt.subplot(13,2,25)
plt.grid(color='grey')
summed = emd_decomp[:,13]
plt.plot(emd_decomp[:,0],summed*1.e9,'blue')
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,23)
plt.grid(color='grey')
summed = summed + emd_decomp[:,12]
plt.plot(emd_decomp[:,0],summed*1.e9,'blue')
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,21)
plt.grid(color='grey')
summed = summed + emd_decomp[:,11]
plt.plot(emd_decomp[:,0],summed*1.e9,'blue')
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,19)
plt.grid(color='grey')
summed = summed + emd_decomp[:,10]
plt.plot(emd_decomp[:,0],summed*1.e9,'blue')
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,17)
plt.grid(color='grey')
summed = summed + emd_decomp[:,9]
plt.plot(emd_decomp[:,0],summed*1.e9,'blue')
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,15)
plt.grid(color='grey')
summed = summed + emd_decomp[:,8]
plt.plot(emd_decomp[:,0],summed*1.e9,'blue')
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,13)
plt.grid(color='grey')
summed = summed + emd_decomp[:,7]
plt.plot(emd_decomp[:,0],summed*1.e9,'blue')
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(13,2,11)
plt.grid(color='grey')
summed = summed + emd_decomp[:,6]
plt.plot(emd_decomp[:,0],summed*1.e9,'blue')
plt.xlabel('Seconds of day',fontsize=15)



plt.show()

