#!/usr/bin/env python

from h5py  import File
import sys
import matplotlib.pyplot as plt
import numpy as np
orbA = File(sys.argv[1])
tA = orbA['time'][:]
coordsA = orbA['coords'][:]


orbB = File(sys.argv[2])
tB = orbB['time'][:]
coordsB = orbB['coords'][:]

nt = min(len(tA),len(tB))

fig, ax = plt.subplots(nrows=6, sharex=True)
#print(tA-tA[0])
#print(coordsA[:,0]-coordsB[:,0])
ax[0].plot(tA[:nt]-tA[0], (coordsA[:nt,0]-coordsB[:nt,0])*1e3,linewidth=2,alpha=0.8)
ax[0].set_ylabel('DeltaX [mm]')
ax[1].plot(tA[:nt]-tA[0], (coordsA[:nt,1]-coordsB[:nt,1])*1e3,linewidth=2,alpha=0.8)
ax[1].set_ylabel('DeltaY [mm]')
ax[2].plot(tA[:nt]-tA[0], (coordsA[:nt,2]-coordsB[:nt,2])*1e3,linewidth=2,alpha=0.8)
ax[2].set_ylabel('DeltaZ [mm]')
ax[3].plot(tA[:nt]-tA[0], (coordsA[:nt,3]-coordsB[:nt,3])*1e3,linewidth=2,alpha=0.8)
ax[3].set_ylabel('DeltaVx [mm/s]')
ax[4].plot(tA[:nt]-tA[0], (coordsA[:nt,4]-coordsB[:nt,4])*1e3,linewidth=2,alpha=0.8)
ax[4].set_ylabel('DeltaVy [mm/s]')
ax[5].plot(tA[:nt]-tA[0], (coordsA[:nt,5]-coordsB[:nt,5])*1e3,linewidth=2,alpha=0.8)
ax[5].set_ylabel('DeltaX [mm/s]')


if len(sys.argv) == 4:
    prefit = np.loadtxt(sys.argv[3], skiprows = 1)
    for i in range(6):
        if i < 3:
            scale = 1e3
        else:
            scale = 1
        ax[i].plot(tA[:nt]-tA[0], prefit[:,i+4]*scale,linewidth=1,alpha=0.8)
    print(prefit)


plt.show()
