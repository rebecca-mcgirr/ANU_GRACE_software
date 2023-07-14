def read_lowpass_file(apr_lowpass_file,msc_number):

    # read a lowpass filter file and extract the model values for a particular mascon (annual signal, lowpass values)\
    # P. Tregoning
    # 20 June 2022

    import numpy as np

    # get the list of decimal years of the epochs within the file
    decimal_year = np.array(np.genfromtxt(apr_lowpass_file ,skip_header=3,max_rows=1))
    #print(decimal_year)

    # get the model values of a particular mascon ...
    msc_model = np.array(np.genfromtxt(apr_lowpass_file,skip_header=msc_number+4,max_rows=1))
    #print(msc_model)
    
    
    return (decimal_year,msc_model[2],msc_model[3:5],msc_model[10:])
    

####### Main Python Script starts here ######   
import numpy as np
import matplotlib.pyplot as plt
import sys

# read apriori low-pass file and mascon number from command line
apr_lowpass_file = sys.argv[1]
msc_number = int(sys.argv[2])
print("Will extract a priori values for mascon ",msc_number," from file: ",apr_lowpass_file)

#msc_number = 2000

apr_ann_ampl = np.zeros((2))
apr_dec_yr,apr_density,apr_ann_ampl[:],apr_lowpass = read_lowpass_file(apr_lowpass_file,msc_number)

#print("Decimal years of epochs in lowpass file:",apr_dec_yr)
#print("Annual cosine/sine amplitudes for mascon",msc_number,": ",apr_ann_ampl)
#print("Low-pass filter values:",apr_lowpass)

# calculate the annual signal as a time series
ann_signal = np.zeros(len(apr_dec_yr))
#print("dimensions of apr_dec_yr and ann_signal:",np.shape(apr_dec_yr),np.shape(ann_signal))

for iepoch in range(len(apr_dec_yr)):
    ann_signal[iepoch] = apr_ann_ampl[0]*np.cos(apr_dec_yr[iepoch]*2*np.pi) + apr_ann_ampl[1]*np.sin(apr_dec_yr[iepoch]*2*np.pi)

#print("dimensions of apr_dec_yr and ann_signal:",np.shape(apr_dec_yr),np.shape(ann_signal))

# plot it just to check what it is doing
fig = plt.figure(figsize=(10,5))

plt.plot(apr_dec_yr,apr_lowpass,label="lowpass")
plt.plot(apr_dec_yr,ann_signal,label="annual")
plt.plot(apr_dec_yr,apr_lowpass+ann_signal,label="lowpass+annual")
plt.legend(loc=3)
plt.title("A priori file:" + apr_lowpass_file + "     Mascon: " + str(msc_number))

plt.show()

