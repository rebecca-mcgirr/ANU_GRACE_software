# Python script that decomposes acc1B data created by Sebastien, modified by Rebecca
from gracetools.io.level1b import  grace1bdata
import matplotlib.pyplot as plt
import datetime
# from gracetools.math.emd import emd
import numpy as np
import sys
from pyeemd import emd
import matplotlib.dates as mdates
from netCDF4 import num2date
plt.style.use(['seaborn-ticks', 'seaborn-paper'])
plt.ion()

print( "Runstring: python ~/gt/com/ACC_EMD.py year month day")

yr = sys.argv[1]
month = sys.argv[2]
day = sys.argv[3]
d = int(sys.argv[3])

print(yr, month, day)

datadir = '/scratch/compute1/geodynamics/grace_data/RL02/%s/' % (yr)
#datadir = '/Users/rebecca/GRACE/gracedata/'

figs = []
axs = []

i=d-1
previous_day = '0%d' % (i) if i < 10 else i 

grace_seconds = int(
    (datetime.datetime(int(yr), int(month), i, 0, 0, 0) - datetime.datetime(2000, 1, 1, 12, 0, 0)).total_seconds())
data = grace1bdata(datadir+'grace_1B_%s-%s-%s_02.tar.gz' % (yr, month, previous_day), start=grace_seconds)
print(data.acc1b)
print(data.acc1b[0]['lin_accl_y'])
print( data.acc1b[0]['gps_time'])
t = data.acc1b[0]['gps_time'][:]
obs2 = data.acc1b[0]['lin_accl_y'][:]

for i in range(d,d+1):
    print(i)
    grace_seconds = int(
        (datetime.datetime(int(yr), int(month), i, 0, 0, 0) - datetime.datetime(2000, 1, 1, 12, 0, 0)).total_seconds())
    data = grace1bdata(datadir+'grace_1B_%s-%s-%s_02.tar.gz' % (yr, month, day), start=grace_seconds)#, cat=True)
    t = np.hstack([t,data.acc1b[0]['gps_time']])
    print(data.acc1b[0]['gps_time'])
    obs2 = np.hstack([obs2,data.acc1b[0]['lin_accl_y']])

#obs = data.acc1b[0]['lin_accl_y']
fig1,ax= plt.subplots()
obs2 -= np.mean(obs2)
t2 = num2date(t, "seconds since 2000-01-01 12:00:00")
plt.plot(t2, obs2)
print('press enter')
raw_input()

# decompose with EMD
imf = emd(obs2)
figt, axt = plt.subplots(nrows=len(imf) + 1, ncols=2, sharex=True)
# figs.append(figt)
# axs.append(axt)
d = np.zeros_like(imf[0])
for i in range(len(imf)):
    axt[i, 0].plot(imf[i])
    d += imf[i]
    if i > 0:
        axt[i, 1].plot(d)

axt[-1, 0].plot(obs2)
ax.plot(t2, np.sum(imf[0:10],axis=0),'k')
print('press enter')
raw_input()
plt.show()
import sys
sys.exit()

# number of components/IMFs
maximf = 8
for i in range(int(d),int(d)+1):
    grace_seconds = int(
        (datetime.datetime(int(yr), int(month), i, 0, 0, 0) - datetime.datetime(2000, 1, 1, 12, 0, 0)).total_seconds())
    data = grace1bdata(datadir+'grace_1B_%s-%s-%s_02.tar.gz' % (yr, month, day), start=grace_seconds)
    obs = data.acc1b[0]['lin_accl_y']
    # plt.figure()
    obs -= np.mean(obs)
    # plt.plot(obs)

    imf = emd(obs)#,nIMF=maximf)
    print(len(imf))
    figt ,axt = plt.subplots(nrows=len(imf)+1,ncols=2,sharex=True)
    figs.append(figt)
    axs.append(axt)
    d = np.zeros_like(imf[0])
    for i in range(len(imf)):
        axs[-1][i,0].plot(imf[i])
        d += imf[i]
        if i>0:
            axs[-1][i,1].plot(d)

    axs[-1][-1,0].plot(obs)

#plt.plot(data.acc1b[1]['lin_accl_y'])
plt.draw()

print('press enter')
raw_input()
