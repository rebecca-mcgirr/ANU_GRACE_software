import sys
import numpy as np
import matplotlib.pyplot as plt

#Import Arguements
yr = str(sys.argv[1])
month = str(sys.argv[2])
day = str(sys.argv[3])

RR_Masked = np.loadtxt('plt_msc_{yr}_{month}_{day}_masked.kb'.format(yr = yr, month = month, day = day), skiprows=6, usecols = [0,1,11])
RR_Nomask = np.loadtxt('plt_msc_{yr}_{month}_{day}_nomask.kb'.format(yr = yr, month = month, day = day), skiprows=6, usecols = [0,1,11])

RA_Masked = np.loadtxt('plt_msc_{yr}_{month}_{day}_masked.kb'.format(yr = yr, month = month, day = day), skiprows=6, usecols = [0,1,11])
RA_Nomask = np.loadtxt('plt_msc_{yr}_{month}_{day}_nomask.kb'.format(yr = yr, month = month, day = day), skiprows=6, usecols = [0,1,11])

'''
RR_Masked = np.loadtxt('plt_msc_2017_01_15_masked.kb', skiprows=6, usecols = [0,1,11])
#RR_Masked = np.loadtxt('plt_msc_2016_12_07_masked3.kb', skiprows=6, usecols = [0,1,11])
#RR_Masked = np.loadtxt('plt_msc_2016_12_07_masked_nofft.kb', skiprows=6, usecols = [0,1,11])
#RR_Masked = np.loadtxt('plt_msc_2016_12_07_masked2_nofft.kb', skiprows=6, usecols = [0,1,11])

RR_Nomask = np.loadtxt('plt_msc_2017_01_15_nomask.kb', skiprows=6, usecols = [0,1,11])


RA_Masked = np.loadtxt('plt_msc_2017_01_15_masked.kba', skiprows=6, usecols = [0,1,11])
#RA_Masked = np.loadtxt('plt_msc_2016_12_07_masked3.kba', skiprows=6, usecols = [0,1,11])
#RA_Masked = np.loadtxt('plt_msc_2016_12_07_masked_nofft.kba', skiprows=6, usecols = [0,1,11])
#RA_Masked = np.loadtxt('plt_msc_2016_12_07_masked2_nofft.kba', skiprows=6, usecols = [0,1,11])

RA_Nomask = np.loadtxt('plt_msc_2017_01_15_nomask.kba', skiprows=6, usecols = [0,1,11])'''
Outliers = np.loadtxt('outlier.mask', skiprows=2, usecols = [5,6])

RR_Diff = RR_Masked - RR_Nomask
RA_Diff = RA_Masked - RA_Nomask


if (len(Outliers)) > 2:
	plt.subplot(2, 2, 1)
	plt.title("RANGE RATE PREFITS")
	y = [-0.05,0.05]
	for i in range(len(Outliers)):
		A = [Outliers[i,0] + Outliers[0,1], Outliers[i,0] + Outliers[0,1]]
		plt.plot(A, y, '-k')
		B = [Outliers[i,1] + Outliers[0,1], Outliers[i,1] + Outliers[0,1]]
		plt.plot(B, y, '-k')
	#plt.plot(RR_Masked[:,0],RR_Diff[:,1])
	plt.plot(RR_Masked[:,2], '-r')
	plt.plot(RR_Nomask[:,0],RR_Nomask[:,2], '-g')

	plt.subplot(2, 2, 2)
	plt.title("RANGE RATE POSTFITS")
	y = [-0.05,0.05]
	for i in range(len(Outliers)):
		A = [Outliers[i,0] + Outliers[0,1], Outliers[i,0] + Outliers[0,1]]
		plt.plot(A, y, '-k')
		B = [Outliers[i,1] + Outliers[0,1], Outliers[i,1] + Outliers[0,1]]
		plt.plot(B, y, '-k')
	#plt.plot(RR_Masked[:,0],RR_Diff[:,2])
	plt.plot(RR_Masked[:,0],RR_Masked[:,1], '-r')
	plt.plot(RR_Nomask[:,0],RR_Nomask[:,1], '-g')

	plt.subplot(2, 2, 3)
	plt.title("RANGE ACCELERATION PREFITS")
	y = [-0.1,0.1]
	for i in range(len(Outliers)):
		A = [Outliers[i,0] + Outliers[0,1], Outliers[i,0] + Outliers[0,1]]
		plt.plot(A, y, '-k')
		B = [Outliers[i,1] + Outliers[0,1], Outliers[i,1] + Outliers[0,1]]
		plt.plot(B, y, '-k')
	#plt.plot(RA_Masked[:,0],RA_Diff[:,1])
	plt.plot(RA_Masked[:,0],RA_Masked[:,2], '-r')
	plt.plot(RA_Nomask[:,0],RA_Nomask[:,2], '-g')

	plt.subplot(2, 2, 4)
	plt.title("RANGE ACCELERATION POSTFITS")
	y = [-0.1,0.1]
	for i in range(len(Outliers)):
		A = [Outliers[i,0] + Outliers[0,1], Outliers[i,0] + Outliers[0,1]]
		plt.plot(A, y, '-k')
		B = [Outliers[i,1] + Outliers[0,1], Outliers[i,1] + Outliers[0,1]]
		plt.plot(B, y, '-k')
	#plt.plot(RA_Masked[:,0],RA_Diff[:,2])
	plt.plot(RA_Masked[:,0],RA_Masked[:,1], '-r')
	plt.plot(RA_Nomask[:,0],RA_Nomask[:,1], '-g')
else:
	plt.subplot(2, 2, 1)
	plt.title("RANGE RATE PREFITS")
	y = [-0.05,0.05]
	for i in range(len(Outliers)):
		A = [Outliers[0], Outliers[0]]
		plt.plot(A, y, '-k')
		B = [Outliers[1], Outliers[1]]
		plt.plot(B, y, '-k')
	#plt.plot(RR_Masked[:,0],RR_Diff[:,1])
	plt.plot(RR_Masked[:,0],RR_Masked[:,2], '-r')
	plt.plot(RR_Nomask[:,0],RR_Nomask[:,2], '-g')

	plt.subplot(2, 2, 2)
	plt.title("RANGE RATE POSTFITS")
	y = [-0.05,0.05]
	for i in range(len(Outliers)):
		A = [Outliers[0], Outliers[0]]
		plt.plot(A, y, '-k')
		B = [Outliers[1], Outliers[1]]
		plt.plot(B, y, '-k')
	#plt.plot(RR_Masked[:,0],RR_Diff[:,2])
	plt.plot(RR_Masked[:,0],RR_Masked[:,1], '-r')
	plt.plot(RR_Nomask[:,0],RR_Nomask[:,1], '-g')

	plt.subplot(2, 2, 3)
	plt.title("RANGE ACCELERATION PREFITS")
	y = [-0.1,0.1]
	for i in range(len(Outliers)):
		A = [Outliers[0], Outliers[0]]
		plt.plot(A, y, '-k')
		B = [Outliers[1], Outliers[1]]
		plt.plot(B, y, '-k')
	#plt.plot(RA_Masked[:,0],RA_Diff[:,1])
	plt.plot(RA_Masked[:,0],RA_Masked[:,2], '-r')
	plt.plot(RA_Nomask[:,0],RA_Nomask[:,2], '-g')

	plt.subplot(2, 2, 4)
	plt.title("RANGE ACCELERATION POSTFITS")
	y = [-0.1,0.1]
	for i in range(len(Outliers)):
		A = [Outliers[0], Outliers[0]]
		plt.plot(A, y, '-k')
		B = [Outliers[1], Outliers[1]]
		plt.plot(B, y, '-k')
	#plt.plot(RA_Masked[:,0],RA_Diff[:,2])
	plt.plot(RA_Masked[:,0],RA_Masked[:,1], '-r')
	plt.plot(RA_Nomask[:,0],RA_Nomask[:,1], '-g')

plt.show()

