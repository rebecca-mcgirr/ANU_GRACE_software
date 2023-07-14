import numpy as np
import scipy as sc
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.signal import butter
################################################
## Python functions for extending time series ##
##
## calc_annual_signal : LS inversion to solve for annual amplitudes
## calc_lowpass       : uses Julia's low-pass filter to extract the low filter component of a time series of data
## extend_timeseries  : extends a time series, given lowpass filter and annual amplitude inputs

def calc_annual_signal(dec_yr,data):
    # calculate amplitudes of in-phase and out-of-phase annual variation through the original-lowpass_series. Simple LS nversion.
    # P. Tregoning
    # 17 February 2022
    
    n_epochs = len(data)
    A = np.zeros((n_epochs,2))
    B = np.zeros(n_epochs)

    for i in range (n_epochs):
    	A[i,0] = np.cos(dec_yr[i]*np.pi*2)
    	A[i,1] = np.sin(dec_yr[i]*np.pi*2)
    
    	# use a priori values of zero for each amplitude. The OMC is then just the data itself
    	B[i] = data[i]
    
    # perform the least squares solution
    At = np.transpose(A)
    AtA = np.dot(At,A)
    VCV = np.linalg.inv(AtA)
    AtB = np.dot(At,B)
    ann_amplitudes = np.dot(VCV,AtB)
    
    return ann_amplitudes

def calc_lowpass (dec_yr,data,nvals,cutoff):
    # Use Julia's low-pass filter to extract the low filter component of a time series of data
    # P. Tregoning
    # 17 February 2022
    #
    # inputs:   data: 2D array of decimal year and data points
    # output:   lowpass_series: time series of just lowpass components, at the epochs of the input data

    # try to take code from Julia's get_lowpass_values.py ....
    #dec_yr = data[:,0]
    fi=sc.interpolate.interp1d(dec_yr,data,kind='linear',fill_value='extrapolate') 
    data_i=fi(dec_yr)

    # the butterworth filter part
    fs = nvals
    Tc = cutoff # 2-year cutoff on the filter (arbitrarily chosen!)
    fc = 1/Tc
    w = fc / (fs / 2) # Normalize the frequency
    b, a = butter(5, w, 'low')
    lowpass_2yr_i = sc.signal.filtfilt(b, a, data_i)
    # interp back to original sampling
    fi_low=sc.interpolate.interp1d(dec_yr,lowpass_2yr_i,kind='linear' \
    ,fill_value='extrapolate')
    lowpass_series=fi_low(dec_yr)

    return lowpass_series


def extend_timeseries (n_extend,time_inc, epochs, data, lowpass_series, ann_amplitudes):
    # extend a time series of data by extrapolating what is passed in, plus adding an annual variation into the extended bit
    # P. Tregoning
    # 17 February 2022
    #
    # input:
    #  n_extend  : number of epochs to extend the data after the last actual epoch
    #  time_inc. : time step (in decimal year) between extended epochs
    #  epochs.   : vector of actual epochs (in decimal year)
    #  data.     : vector of raw data
    #  lowpass_series : as it says
    # ann_amplitudes : amplitude of cosine/sine for an annual variation

    n_epochs = len(epochs)
    lowpass_plus_annual = np.zeros(n_epochs+n_extend)
    
    # backfill the extended arrays with the original data first
    epochs_extended = np.zeros(n_epochs+n_extend)
    epochs_extended[0:n_epochs] = epochs
    data_extended = np.zeros(n_epochs+n_extend)
    data_extended[0:n_epochs] = data

    # fill in the values for the extended epochs
    for i in range (n_extend):
        epochs_extended[n_epochs+i] = epochs[-1] + (i+1)*1/12

    # extrapolate the lowpass series, using quadratic extension 
    s = InterpolatedUnivariateSpline(epochs[-20:], lowpass_series[-20:], k=2)
    data_extended[-n_extend:] = s(epochs_extended[n_epochs:])

    
    # now add the annual variation
    for i in range (n_extend):
        data_extended[n_epochs+i] = data_extended[n_epochs+i] + \
        ann_amplitudes[0]*np.cos(epochs_extended[n_epochs+i]*2*np.pi) + \
        ann_amplitudes[1]*np.sin(epochs_extended[n_epochs+i]*2*np.pi)

        
    #print("epochs_extended",epochs_extended)
    #print("data_extended",data_extended)
    
    return epochs_extended, data_extended


