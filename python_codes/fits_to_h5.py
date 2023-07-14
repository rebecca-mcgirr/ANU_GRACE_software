#import h5py
#import numpy as np
#import re
import logging
#from pathlib import Path
#import msc
logging.basicConfig(level=logging.NOTSET)


def cmdline():
    import argparse
    parser = argparse.ArgumentParser(description='Convert a set of fit files into a hdf5 file \n example: python3 fits_to_h5.py -yr_start 2003 -yr_end 2022 -o batch7_i2_loose.h5 "/scratch/compute1/pault/ANU_mascons/b8i3/<YYYY>/<MM>/msc_<YYYY>-<MM>_b7i2.fit"  ')
    parser.add_argument('-o', '--output', metavar='output', type=str, 
                    help='output file')
    parser.add_argument('path', metavar='path', type=str, help='path of the fit files.')
    parser.add_argument('-yr_start', metavar='yr_start', type=str, help='start year of the fit files.')
    parser.add_argument('-yr_end', metavar='yr_end', type=str, help='end year of the fit files.')
    parser.add_argument('-mascon_file',metavar='mascon_file',type=str,help='mascon ascii file')
    return parser.parse_args()

def merge_fit_into_h5(yearmin, yearmax, pattern, mascon_file, outfile):
    import os
    from pathlib import Path
    import re
    from h5py import File
    from numpy import zeros, asarray, float32
    #yearmin=2021
    #yearmax=2022

    num_epochs=0
    for yr in range(int(yearmin), int(yearmax)+1):
        for month in range(1,12+1):
            fn = pattern.replace('<YYYY>','%04i'%(yr)).replace('<MM>','%02i'%(month))
            path=Path(fn)
            if path.exists():
                num_epochs+=1

    logging.info(f'there are {num_epochs} input solution files')
    
    # PT220604: we need to know how many mascons in the mascon file. How to do in python?
    with open(mascon_file) as f:
        firstline = f.readline().rstrip()

    print(firstline)
    num_mascons = int(firstline[15:21:1])
    print("num_mascons=",int(num_mascons))
    
    postfit_ewh = zeros((num_epochs, int(num_mascons)))
    prefit_ewh = zeros((num_epochs, int(num_mascons)))
    sigmas_ewh = zeros((num_epochs,int(num_mascons)))
    dec_yr = zeros(num_epochs)
    years = zeros(num_epochs)
    months = zeros(num_epochs)
    years = years.astype(int)
    months = months.astype(int)
    
    print("Reading solution files")
    iepoch=0
    for yr in range(int(yearmin), int(yearmax)+1):
        for month in range(1,12+1):
            fn = pattern.replace('<YYYY>','%04i'%(yr)).replace('<MM>','%02i'%(month))
            path=Path(fn)
            if path.exists():
                logging.info('reading %s'%(fn))
                a = [re.findall(r'^.*MC.*$', line) for line in open(fn)]
                c = [re.findall("[-+]?\d+[\.]?\d*", s[0]) for s in a if s != []]
                cc = asarray(list(map(float32, c)))
                postfit_ewh[iepoch,:] = cc[:, 4]
                prefit_ewh[iepoch,:] = cc[:,2]
                sigmas_ewh[iepoch,:] = cc[:,5]
                dec_yr[iepoch]=yr+float(month-0.5)/12.
                years[iepoch]=int(yr)
                months[iepoch]=int(month)
                iepoch+=1

    print("Writing hdf5 output file")        
    with File(outfile, 'w') as hf:
        hf.attrs['fit_files'] = pattern
        pst = hf.create_dataset('solution/ewh', data=postfit_ewh)
        pst = hf.create_dataset('solution/prefit', data=prefit_ewh)
        pst = hf.create_dataset('solution/sigma_ewh', data=sigmas_ewh)
        tt = hf.create_dataset('time/year', data=years)
        tt = hf.create_dataset('time/decyear', data=dec_yr)
        tt = hf.create_dataset('time/month', data=months)

 
    
    # PT220201: try writing it out as a netcdf as well
    from netCDF4 import Dataset
    import grace_utils as gr
    import numpy as np

    ncname = outfile[:-3]+".nc"
    ncw=Dataset(ncname,'w',format='NETCDF4')
    ncw.description = 'ANU mascon solutions for 2003-2016 (Allgeyer et al, 2022)'

    # PT220201: read the mascon file
    _,Pdata,_,_=gr.read_mascon_file(mascon_file)
    num_msc=len(Pdata[:,0])
    mascon_number=Pdata[:,0]
    
    Pcode=[]
    lat=[]
    lon=[]
    region=[]
    area=[]

    print("Getting mascon coordinate information")
    for pcpt in np.arange(num_msc):
        lat.append(Pdata[pcpt,4])
        lon.append(Pdata[pcpt,5])
        region.append(Pdata[pcpt,-1])
        area.append(Pdata[pcpt,7])
        #print(pcpt,lat[pcpt],lon[pcpt])
        
    # Create dimensions
    ncw.createDimension('num_epochs', num_epochs)
    ncw.createDimension('num_msc', num_msc)

    # PT220201: for the netcdf file the ewh and sigma arrays have the columns reversed
    ewh = zeros((num_msc,num_epochs))
    ewh = np.transpose(postfit_ewh)
    np.shape(ewh)
    sigmas = zeros((num_msc,num_epochs))
    sigmas = np.transpose(sigmas_ewh)
    # PT220615: include the a priori mascons in the netcdf file
    apriori_ewh = zeros((num_msc,num_epochs))
    apriori_ewh = np.transpose(prefit_ewh)
    np.shape(apriori_ewh)
    
    
    print("Writing netcdf file: ",ncname)
    # Create variables
    # decimal year
    decyear_vec= ncw.createVariable('decyear', 'f4', ('num_epochs'))
    decyear_vec.units = 'decimal year'
    decyear_vec.description='Epoch of monthly estimate in decimal year'
    decyear_vec[:]=dec_yr[:]
    # primary mascon
    pcode_vec= ncw.createVariable('mascon_number', 'f4', ('num_msc'))
    pcode_vec.units = '-'
    pcode_vec.description='primary mascon number'
    pcode_vec[:]=mascon_number[:]
    # primary latitude
    lat_vec= ncw.createVariable('lat', 'f4', ('num_msc'))
    lat_vec.units = 'degrees'
    lat_vec.description='latitude of the primary mascon centroid in WGS84'
    lat_vec[:]=lat[:]
    ## primary longitude
    lon_vec= ncw.createVariable('lon', 'f4', ('num_msc'))
    lon_vec.units = 'degrees'
    lon_vec.description='longitude of the primary mascon centroid in WGS84'
    lon_vec[:]=lon[:]
    # drainage basin
    region_vec= ncw.createVariable('drainage_basin', 'S1', ('num_msc'))
    region_vec.units = '-'
    region_vec.description='name of the drainage basin'
    region_vec[:]=region[:]
    #  ewh
    ewh_arr= ncw.createVariable('ewh', 'f4', (('num_msc','num_epochs')))
    ewh_arr.units = 'm'
    ewh_arr.description='estimated mascon equivalent water height anomaly'
    ewh = np.vstack(ewh)
    ewh_arr[:]=ewh[:,:]
    #  apriori_ewh
    apriori_ewh_arr= ncw.createVariable('apriori_ewh', 'f4', (('num_msc','num_epochs')))
    apriori_ewh_arr.units = 'm'
    apriori_ewh_arr.description='apriori mascon equivalent water height anomaly'
    apriori_ewh = np.vstack(apriori_ewh)
    apriori_ewh_arr[:]=apriori_ewh[:,:]
    #  sigma
    sigma_arr= ncw.createVariable('sigma', 'f4', (('num_msc','num_epochs')))
    sigma_arr.units = '(m)'
    sigma_arr.description='formal uncertainty of mascon EWH'
    sigmas = np.vstack(sigmas)
    sigma_arr[:]=sigmas[:,:]
    # Close NetCF file
    ncw.close()



if __name__ == "__main__":
    cmd = cmdline()
    logging.info("Will run the code on the fit files : %s, output will be %s, mascon_file %s", cmd.path, cmd.output, cmd.mascon_file)
    merge_fit_into_h5(cmd.yr_start, cmd.yr_end, cmd.path, cmd.mascon_file, cmd.output)


import sys
sys.exit()



