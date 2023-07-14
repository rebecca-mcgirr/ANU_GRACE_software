#!/usr/bin/env python3
## standard packages
import os
from fnmatch import fnmatch
import numpy as np
import datetime as dt
import argparse
import scipy as sc
from scipy.signal import butter
## subroutines
import grace_utils as gr
###############################################################################
## READ ARGUMENTS
###############################################################################
# Define arguments
parser = argparse.ArgumentParser(description=' This programs aims to extract annual, low-frequency and high-frequency signals from postfit values. It will write the output in an asci file, containing latitude, longitude, annual amplitude, rms of total signal, rms of annual signal, rms of lowpass filtered series, rms of high frequencies, rms of annual + high frequencies and the monthly values of the lowpass filtered series.')
parser.add_argument("datadir",default="/scratch/compute1/julia/solutions/monthly_solutions",type=str,help="complete path to your monthly solutions")
parser.add_argument("pattern",default='????/??/addnorm*.fit',type=str,help="pattern describing fitfiles")
parser.add_argument("output_file",default='my_output_file.txt',type=str,help="complete path to your output file")
parser.add_argument("--masconfile",default='/scratch/compute1/geodynamics/software/gtgk_1064/grace/tables/mascons_stage4_V003a',type=str,help="complete path to yourmascon file (default=mascons_stage4_V003a)")
parser.add_argument("--cutoff",default=0.5,type=float,help="cutoff period in years (default=0.5 i.e. 6 months)")
parser.add_argument("--min_days",default=15,type=int,help="min number of days in monthly solution")
parser.add_argument("--calc_annual",default='yes',type=str,help="yes/no to calculate annual periodic signal")

# Read arguments
args = parser.parse_args()
# Declare arguments
filedir = args.datadir
mypattern= args.pattern
outfile= args.output_file
masconfile= args.masconfile
Tc=args.cutoff
min_days = args.min_days
calc_annual = args.calc_annual

# PT220607: strip off any directory information from the mascon file name
found = 0
start_pos = 0
while found > -1 and start_pos < len(masconfile):
    found = masconfile.find('/',start_pos+1,-1)
    if found > -1:
        start_pos = found
trimmed_masconfile = masconfile[start_pos+1:]
###############################################################################
##  READ FITFILES
###############################################################################
if filedir[-1]=='/':
    pattern=filedir+mypattern
else:
    pattern=filedir+'/'+mypattern
print('Make list of files fitting pattern %s: %s'%(pattern,str(dt.datetime.now())[0:19]))
filelist=[]
for path, subdirs, files in os.walk(filedir):
    for name in files:
        if fnmatch(os.path.join(path, name), pattern):
            #PT220702: find out whether there are sufficient days in the monthly solution file
            ndays=int(np.array( np.genfromtxt(os.path.join(path, name), skip_header=1,max_rows=1))[11])
            if ndays >= min_days:
                filelist.append(os.path.join(path, name))
            else:
                print("File ",os.path.join(path, name)," has only ",ndays," days. Do not use.")
                
filelist.sort()
nbt=len(filelist)
# Read first file to get number of mascons
print('Read \'%s\' to get the number of mascons: %s'%(masconfile,str(dt.datetime.now())[0:19]))
f=open(masconfile)
lines=f.readlines()
f.close()
nb_msc=int(lines[0].split()[1])
del f,lines
mascon_number=np.zeros(nb_msc)
print('Found %d mascons in \'%s\': %s'%(nb_msc,filelist[0],str(dt.datetime.now())[0:19]))
# Create arrays to store data
print('Read %d fitfiles and create nd-arrays with prefit and postfit statistics, and, satellite and mascon parameters: %s'%(nbt,str(dt.datetime.now())[0:19]))
decyear=np.zeros(nbt)
sigmas=np.zeros((nbt,nb_msc))
params_msc=9999*np.ones((nbt,nb_msc,4))
nbdays=np.zeros(nbt)
tcount=-1
#print(filelist)
for fitfile in filelist:
    print(fitfile)
    if gr.read_addnorm_fit(fitfile,nb_msc)==None:
        pass
    else:
        tcount=tcount+1
        decyear[tcount],params_msc[tcount,:,:],nbdays[tcount],sigmas[tcount,:]=gr.read_addnorm_fit(fitfile,nb_msc)
print('Found %d non-empty fitfiles in \'%s\': %s'%(tcount,filelist[0],str(dt.datetime.now())[0:19]))
decyear=decyear[:tcount+1]
postfit=params_msc[:tcount+1,:,2]
NT=len(decyear) # nb times
NB=len(postfit[0,:]) # nb mascons
###############################################################################
## COMPUTE RMS AND LOWPASS
###############################################################################
print('Separate annual, lowpass and hf signals for %d mascons: %s'%(NB,str(dt.datetime.now())[0:19]))
#arrays
annual_coefs=np.zeros((NB,2))
annual_series=np.zeros((NT,NB))
lowpass_series=np.zeros((NT,NB))
lowpass_series_2nd=np.zeros((NT,NB))
hf_series=np.zeros((NT,NB))
for cpt in np.arange(NB):
    mean_mascon=np.mean(postfit[:,cpt])
    data=postfit[:,cpt]-mean_mascon*np.ones(NT)
    # linear interpolation
    decyear_i=np.arange(decyear[0]-1/12,decyear[-1]+1/12,1/12)
    fi=sc.interpolate.interp1d(decyear, data,kind='linear',fill_value='extrapolate')
    data_i=fi(decyear_i)
    # low pass filter on interpolated data
    fs = 12
    fc = 1/Tc  # Cut-off frequency of the filter
    w = fc / (fs / 2) # Normalize the frequency
    b, a = butter(5, w, 'low')
    lowpass_2yr_i = sc.signal.filtfilt(b, a, data_i)
    # interp back to original sampling
    fi_low=sc.interpolate.interp1d(decyear_i, lowpass_2yr_i,kind='linear',fill_value='extrapolate')
    lowpass_2yr=fi_low(decyear)
    lowpass_series[:,cpt]=lowpass_2yr+mean_mascon*np.ones(NT)
    hf_series[:,cpt]=data-lowpass_2yr

    ######## Annual Sinusoid model
    # PT220620: this was missing completely from this python script, although it had been in earler versions
    #           Original approach was to fit an annual sinusoid to the original data. Modify it here th
    #           remove first the low-pass filter then fit to the residuals
    # annual sinusoid
    if calc_annual == "yes" :
        model=np.ones((NT,3))
        model[:,0]=np.cos(decyear*2*np.pi)
        model[:,1]=np.sin(decyear*2*np.pi)
        # OLS inversion of an annual sinusoid. Note, we remove the lowpass_series, whereas Julia did not.
        coefs=np.linalg.lstsq(model,data-lowpass_series[:,cpt])
        coefs=coefs[0]
        annual_predictions=np.dot(model,coefs)
        annual_residuals=data-annual_predictions

        # PT220620: save correctly the annual sinusoid model values
        annual_series[:,cpt]=annual_predictions
        annual_coefs[cpt,0]=coefs[0]
        annual_coefs[cpt,1]=coefs[1]
    else :
        annual_residuals = np.zeros_like(data)
    
    # PT/HM 220620: now remove the annual from the original series, then low-pass the residual
    fi_ann_resid=sc.interpolate.interp1d(decyear, annual_residuals,kind='linear',fill_value='extrapolate')
    ann_resid_i=fi(decyear_i)

    lowpass_2nd_i = sc.signal.filtfilt(b, a, ann_resid_i)
    #print("lowpass_2nd_i:",lowpass_2nd_i)
    # interp back to original sampling
    fi_low=sc.interpolate.interp1d(decyear_i, lowpass_2nd_i,kind='linear',fill_value='extrapolate')
    lowpass_2nd=fi_low(decyear)
    lowpass_series_2nd[:,cpt]=lowpass_2nd+mean_mascon*np.ones(NT)
    hf_series[:,cpt]=data-lowpass_2nd

    # PT220620: DEBUG DEBUG. set lowpass_series_2nd to be the orig minus annual
    #lowpass_series_2nd[:,cpt] = annual_residuals
    
        
## Compute RMS
print('Compute RMS for total, annual, lowpass and hf signals: %s'%(str(dt.datetime.now())[0:19]))
rms_annual=np.sqrt(np.sum(np.square(annual_series),axis=0)/NT)
rms_lowpass=np.sqrt(np.sum(np.square(lowpass_series_2nd),axis=0)/NT)
rms_highfreq=np.sqrt(np.sum(np.square(hf_series),axis=0)/NT)
rms_total=np.sqrt(np.sum(np.square(postfit),axis=0)/NT)
rms_hfan=np.sqrt(np.sum(np.square(hf_series+annual_series),axis=0)/NT)
###############################################################################
## WRITE RESULTS IN ASCII FILE
###############################################################################
print('Read %s to get latitudes and longitudes of primary mascons: %s'%(masconfile,str(dt.datetime.now())[0:19]))
#mascon file
_,Pdata,_,_=gr.read_mascon_file(masconfile)
lat=Pdata[:,4]
lon=Pdata[:,5]
density=Pdata[:,10]
print('Write RMS and lowpass monthly values in %s: %s '%(outfile,str(dt.datetime.now())[0:19]))
ofile=open(outfile, 'w')
ofile.write('# n_msc n_epochs n_columns mascon_file cutoff_freq(years)\n')
ofile.write('{:5d}\t{:5d}\t{:5d}\t{:s}\t{:f}\n'.format(NB,NT,NT+10,trimmed_masconfile,Tc))
ofile.write('# Epochs in decimal year \n')
for tpt in np.arange(NT):
        ofile.write('{:10.5f}\t'.format(decyear[tpt]))
ofile.write('\n')
ofile.write('# Latitude \t Longitude \t Density \t annual_cos \t annual_sin \t total_rms \t annual_rms \t lowpass_rms \t HF_rms  \t annual+HF_rms  \t lowpass_0 ... lowpass_%d\n'%(NT))
for cpt in np.arange(NB):
    ofile.write('{:10.4f}\t{:10.4f}\t{:10.4f}\t{:10.4f}\t{:10.4f}\t{:10.4f}\t{:10.4f}\t{:10.4f}\t{:10.4f}\t{:10.4f}'.format(lat[cpt],lon[cpt],density[cpt],annual_coefs[cpt,0],annual_coefs[cpt,1],rms_total[cpt],rms_annual[cpt],rms_lowpass[cpt],rms_highfreq[cpt],rms_hfan[cpt]))
    for tpt in np.arange(NT):
        ofile.write('{:10.4f}\t'.format(lowpass_series_2nd[tpt,cpt]))
    ofile.write('\n')

ofile.close()
print('%s: You\'re all done! '%(str(dt.datetime.now())[0:19]))
###############################################################################
## END
###############################################################################
