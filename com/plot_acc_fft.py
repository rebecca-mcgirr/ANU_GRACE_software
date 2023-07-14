import numpy as np
import matplotlib
import sys
import tkinter #python 3 syntax
root = tkinter.Tk()
root.withdraw()


#matplotlib.use('TkAgg')  
import matplotlib.pyplot as plt

print( "Runstring: python ~/ga/com/plot_acc_fft.py <input file>")

# PT180706: set the window size based on the screen size. Make it 80% of the height and 90% of the width
width, height = root.winfo_screenwidth(), root.winfo_screenheight()

with open(sys.argv[1]) as f:
    line = f.readline()

cutoffs = line.split()

# read in the GRACE A accelerometer data
print ("reading in the data")
fft_spec = np.loadtxt(sys.argv[1],skiprows=1)
flo = float(cutoffs[0])
fhi = float(cutoffs[1])

nfft = len(fft_spec[1:,0])

# plot the decomposed accelerometer values, one component per figure box
print("Plotting power spectrum")
plt.figure(2,figsize=( int(width*0.012), int(height*0.011)) )

plt.subplot(311)
plt.grid(color='grey')
plt.plot(fft_spec[:,0],fft_spec[:,1])
ymin, ymax = plt.ylim()
xshift = 0.000005
plt.ylabel('Power spectrum (unfiltered)',fontsize=15)
plt.xlabel('Frequency (cycles/sec)',fontsize=15)
plt.plot([fhi,fhi],[0,ymax],'k--')
plt.text(fhi - xshift,0.9*ymax,'High frequency cutoff',rotation=90,fontsize=10,va='top')
plt.plot([flo,flo],[0,ymax],'k--')
plt.text(flo - xshift,0.9*ymax,'Low frequency cutoff',rotation=90,fontsize=10,va='top')
plt.xlim(0,0.0005)
plt.ylim(ymin,ymax)

plt.subplot(312)
plt.grid(color='grey')
plt.plot(fft_spec[:,0],fft_spec[:,2])
ymin, ymax = plt.ylim()
xshift = 0.000005
plt.ylabel('Power spectrum (filtered)',fontsize=15)
plt.xlabel('Frequency (cycles/sec)',fontsize=15)
plt.plot([fhi,fhi],[0,ymax],'k--')
plt.text(fhi - xshift,0.9*ymax,'High frequency cutoff',rotation=90,fontsize=10,va='top')
plt.plot([flo,flo],[0,ymax],'k--')
plt.text(flo - xshift,0.9*ymax,'Low frequency cutoff',rotation=90,fontsize=10,va='top')
plt.xlim(0,0.0005)
plt.ylim(ymin,ymax)

plt.subplot(313)
plt.grid(color='grey')
plt.plot(fft_spec[:,0],fft_spec[:,3])
ymin, ymax = plt.ylim()
xshift = 0.000005
plt.ylabel('Filter',fontsize=15)
plt.xlabel('Frequency (cycles/sec)',fontsize=15)
plt.plot([fhi,fhi],[0,ymax],'k--')
plt.text(fhi - xshift,0.9*ymax,'High frequency cutoff',rotation=90,fontsize=10,va='top')
plt.plot([flo,flo],[0,ymax],'k--')
plt.text(flo - xshift,0.9*ymax,'Low frequency cutoff',rotation=90,fontsize=10,va='top')
plt.xlim(0,0.0005)
plt.ylim(ymin,ymax)
plt.show()
