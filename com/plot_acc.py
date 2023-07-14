
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

print( "Runstring: python ~/ga/com/plot_acc.py GRACE_file_2 GRACE_file_1 1/0 (0: do not swap sign of cross/radial on 2nd file, 1: swap sign)")

# PT180706: set the window size based on the screen size. Make it 80% of the height and 90% of the width
width, height = root.winfo_screenwidth(), root.winfo_screenheight()

# In[13]:

# read in the GRACE A accelerometer data
print ("reading in the data")
acc_A = np.loadtxt(sys.argv[1])
acc_B = np.loadtxt(sys.argv[2])
swap  = sys.argv[3]

#acc_A = np.loadtxt('junk.A')
#acc_B = np.loadtxt('junk.B')

# reverse the sign of GRACE B along-track and cross-track to account for the different orientation
if (swap == 1):
	print ("Swap sign of along-track and cross-track for second satellite")
	acc_B[:,1:3] = -acc_B[:,1:3]

nepochs = len(acc_A[:,0])


# In[ ]:

# plot the accelerometer values for both satellites
print("Plotting ..")
plt.figure(3,figsize=( int(width*0.012), int(height*0.010)) )

plt.subplot(311)
plt.grid(color='grey')
plt.plot(acc_A[:,0]-acc_A[0,0],acc_A[:,1]*1.e9,'red')
plt.plot(acc_B[:,0]-acc_A[0,0],acc_B[:,1]*1.e9,'blue')
plt.plot(acc_B[:,0]-acc_A[0,0],((acc_B[:,1])+np.mean(acc_B[:,1])-acc_A[:,1])*1.e9,'green')
plt.ylabel('Along-track (nm/s$^2$)',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)

plt.subplot(312)
plt.grid(color='grey')
plt.plot(acc_A[:,0]-acc_A[0,0],acc_A[:,2]*1.e9,'red')
plt.plot(acc_B[:,0]-acc_A[0,0],acc_B[:,2]*1.e9,'blue')
plt.plot(acc_B[:,0]-acc_A[0,0],((acc_A[:,2]-acc_B[:,2])+np.mean(acc_B[:,2]))*1.e9,'green')
plt.ylabel('Cross-track (nm/s$^2$)',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)
print("Min and max of removed drift = ",np.min(acc_A[:,2]-acc_B[:,2])*1.e9,np.max(acc_A[:,2]-acc_B[:,2])*1.e9)

plt.subplot(313)
plt.grid(color='grey')
plt.plot(acc_A[:,0]-acc_A[0,0],acc_A[:,3]*1.e9,'red',label=sys.argv[1])
plt.plot(acc_B[:,0]-acc_A[0,0],acc_B[:,3]*1.e9,'blue',label=sys.argv[2])
#plt.plot(acc_B[:,0]-acc_A[0,0],((acc_A[:,3]-acc_B[:,3])+np.mean(acc_B[:,3]))*1.e9,'green',label="difference")
plt.ylabel('Radial (nm/s$^2$)',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)
plt.legend(loc=4)

plt.show()

