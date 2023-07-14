
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

print( "Runstring: python ~/ga/com/plot_AHK_temps.py GRACE_file_1")

# PT180706: set the window size based on the screen size. Make it 80% of the height and 90% of the width
width, height = root.winfo_screenwidth(), root.winfo_screenheight()

# In[13]:

# read in the GRACE temperature data
print ("reading in the data")
ahk_A = np.loadtxt(sys.argv[1])


nepochs = len(ahk_A[:,0])


# In[ ]:

# plot the temperature values for a single satellites
print("Plotting ..")
plt.figure(1,figsize=( int(width*0.012), int(height*0.010)) )

plt.subplot(111)
plt.grid(color='grey')
plt.plot(ahk_A[:,0]-ahk_A[0,0],ahk_A[:,1],'red',  label="sensor unit")
plt.plot(ahk_A[:,0]-ahk_A[0,0],ahk_A[:,2],'blue', label="internal core")
plt.plot(ahk_A[:,0]-ahk_A[0,0],ahk_A[:,3],'green',label="ICU power board")
plt.plot(ahk_A[:,0]-ahk_A[0,0],ahk_A[:,4],'black',label="ICU A/D convertor")
plt.ylabel('Temperature (degrees C)',fontsize=15)
plt.xlabel('Seconds of day',fontsize=15)
plt.legend(loc=4)





plt.show()

