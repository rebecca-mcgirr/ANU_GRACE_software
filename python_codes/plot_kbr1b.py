import numpy as np
import sys
import matplotlib.pyplot as plt

# read in the kbr1b file
kbr1b_file = sys.argv[1]
data = np.loadtxt(kbr1b_file,skiprows=163)
gracesec = data[:,0]
kbrr = data[:,1]



# plot
f1 = plt.figure(figsize=(15,10))

plt.plot(gracesec,kbrr,'o')


plt.show()


