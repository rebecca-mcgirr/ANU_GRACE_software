#!/usr/bin/env python3
import os, re, argparse, logging
import numpy as np
from pathlib import Path
from h5py import File
from gracetools.io.mascons import mascons
from grace_utils import *
from netCDF4 import Dataset

if __name__ == "__main__":

    logging.basicConfig(level=logging.NOTSET)

    # RM221122: removed unecessary functions
    runstr = 'Convert a set of fit files into a hdf5 file \n eg: python3 fitlist_to_h5.py -list_file fit-10day.lst -o batch7_i2_loose.h5 -mascon_file mascons_stage5_V006_200km'
    
    parser = argparse.ArgumentParser(description=runstr)
    parser.add_argument('-o', '--output', metavar='output', type=str,help='output file')
    #parser.add_argument('path', metavar='path', type=str, help='path of the fit files.')
    #parser.add_argument('-yr_start', metavar='yr_start', type=str, help='start year of the fit files.')
    #parser.add_argument('-yr_end', metavar='yr_end', type=str, help='end year of the fit files.')
    parser.add_argument('-mascon_file',metavar='mascon_file',type=str,help='mascon ascii file')
    parser.add_argument('-list_file',metavar='list_file',type=str,help='file containing list of fit files')
    args = parser.parse_args()
    logging.info("Will run the code on the fit files : %s, output will be %s, mascon_file %s", args.list_file, args.output, args.mascon_file)

    list_file = args.list_file
    mascon_file = args.mascon_file
    outfile = args.output

    num_epochs=0
    # PT220804: read all the names of fit files from an input file
    file_list = np.asarray(open(list_file, 'r').read().splitlines())
    num_epochs = len(file_list)
    
    logging.info(f'there are {num_epochs} input solution files')
    print(file_list[:])
    # PT220604: we need to know how many mascons in the mascon file. How to do in python?
    with open(mascon_file) as f:
        firstline = f.readline().rstrip()

    print(firstline)
    num_mascons = int(firstline[15:21:1])
    print("num_mascons=",int(num_mascons))
    
    postfit_ewh = np.zeros((num_epochs, int(num_mascons)))
    prefit_ewh = np.zeros((num_epochs, int(num_mascons)))
    sigmas_ewh = np.zeros((num_epochs,int(num_mascons)))
    dec_yr = np.zeros(num_epochs)
    years = np.zeros(num_epochs)
    months = np.zeros(num_epochs)
    years = years.astype(int)
    months = months.astype(int)
    
    print("Reading solution files")
    iepoch=0
    for iepoch in range (num_epochs):
        fn = str(file_list[iepoch])
        path=Path(fn)
        if path.exists():
            logging.info('reading %s'%(fn))
            a = [re.findall(r'^.*MC.*$', line) for line in open(fn)]
            c = [re.findall("[-+]?\d+[\.]?\d*", s[0]) for s in a if s != []]           
            cc = np.asarray(list(map(np.float32, c)))           
            postfit_ewh[iepoch,:] = cc[:, 4]
            prefit_ewh[iepoch,:] = cc[:,2]
            sigmas_ewh[iepoch,:] = cc[:,5]
            
            vcv = open(fn)
            # PT221012: get the number of solution files, within this fit file, from the second line
            for i in range (2):
                nextline = vcv.readline().rstrip()
                
            nsolns = int(nextline[75:82])
            print("There are ",nsolns," days in this file")
            
            # PT220805: get the year/month for the first solution from the 5th line of the addnorm file
            for i in range (2):
                nextline = vcv.readline().rstrip()

            yr = int(nextline[19:23])
            month = int(nextline[23:28])
            day = int(nextline[29:33])

            # RM221122: fixed bug in dec_yr calculation
            dec_yr[iepoch]=toYearFraction(datetime.date(yr, month, day))
            years[iepoch]=yr
            months[iepoch]=month

            # PT221012: get the year/month for the last solution of the addnorm file
            if nsolns > 1:
                for i in range (nsolns-1):
                    nextline = vcv.readline().rstrip()
                
                yr_end = int(nextline[19:23])
                month_end = int(nextline[23:28])
                day_end = int(nextline[29:33])
                # RM221122: fixed bug in dec_yr calculation
                dec_yr_end=toYearFraction(datetime.date(yr_end, month_end, day_end))
                dec_yr[iepoch] = 0.5 * (dec_yr[iepoch] + dec_yr_end)
                
        else:
            print("File ",fn," does not exist")


    print(f"reading {mascon_file}.h5")
    mANU = mascons(mascon_file+".h5")
    mANU.read_mascons_h5()
    mANU.get_mapping_sht()

    print("Writing hdf5 output file")        
    with File(outfile, 'w') as hf:
        hf.attrs['fit_files'] = list_file
        pst = hf.create_dataset('solution/ewh', data=postfit_ewh)
        pst = hf.create_dataset('solution/prefit', data=prefit_ewh)
        pst = hf.create_dataset('solution/sigma_ewh', data=sigmas_ewh)
        tt = hf.create_dataset('time/year', data=years)
        tt = hf.create_dataset('time/decyear', data=dec_yr)
        tt = hf.create_dataset('time/month', data=months)

        hf.create_dataset('mascon/Plat', data=mANU.Plat)
        hf.create_dataset('mascon/Plon', data=mANU.Plon)
        hf.create_dataset('mascon/Pdensity', data=mANU.Pdensity)
        hf.create_dataset('mascon/Parea', data=mANU.Parea)

 
    
    # PT220201: try writing it out as a netcdf as well
    ncname = outfile[:-3]+".nc"
    ds=Dataset(ncname,'w',format='NETCDF4')
    ds.description = 'ANU mascon solutions for 2003-2016 (Allgeyer et al, 2022)'

    # PT220201: read the mascon file
    _,Pdata,_,_=read_mascon_file(mascon_file)
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
    ds.createDimension('num_epochs', num_epochs)
    ds.createDimension('num_msc', num_msc)

    # PT220201: for the netcdf file the ewh and sigma arrays have the columns reversed
    ewh = np.zeros((num_msc,num_epochs))
    ewh = np.transpose(postfit_ewh)
    np.shape(ewh)
    sigmas = np.zeros((num_msc,num_epochs))
    sigmas = np.transpose(sigmas_ewh)
    # PT220615: include the a priori mascons in the netcdf file
    apriori_ewh = np.zeros((num_msc,num_epochs))
    apriori_ewh = np.transpose(prefit_ewh)
    np.shape(apriori_ewh)
    
    
    print("Writing netcdf file: ",ncname)
    # Create variables
    # decimal year
    decyear_vec= ds.createVariable('decyear', 'f4', ('num_epochs'))
    decyear_vec.units = 'decimal year'
    decyear_vec.description='Epoch of monthly estimate in decimal year'
    decyear_vec[:]=dec_yr[:]
    # primary mascon
    pcode_vec= ds.createVariable('mascon_number', 'f4', ('num_msc'))
    pcode_vec.units = '-'
    pcode_vec.description='primary mascon number'
    pcode_vec[:]=mascon_number[:]
    # primary latitude
    lat_vec= ds.createVariable('lat', 'f4', ('num_msc'))
    lat_vec.units = 'degrees'
    lat_vec.description='latitude of the primary mascon centroid in WGS84'
    lat_vec[:]=lat[:]
    ## primary longitude
    lon_vec= ds.createVariable('lon', 'f4', ('num_msc'))
    lon_vec.units = 'degrees'
    lon_vec.description='longitude of the primary mascon centroid in WGS84'
    lon_vec[:]=lon[:]
    # drainage basin
    region_vec= ds.createVariable('drainage_basin', 'S1', ('num_msc'))
    region_vec.units = '-'
    region_vec.description='name of the drainage basin'
    region_vec[:]=region[:]
    #  ewh
    ewh_arr= ds.createVariable('ewh', 'f4', (('num_msc','num_epochs')))
    ewh_arr.units = 'm'
    ewh_arr.description='estimated mascon equivalent water height anomaly'
    ewh = np.vstack(ewh)
    ewh_arr[:]=ewh[:,:]
    #  apriori_ewh
    apriori_ewh_arr= ds.createVariable('apriori_ewh', 'f4', (('num_msc','num_epochs')))
    apriori_ewh_arr.units = 'm'
    apriori_ewh_arr.description='apriori mascon equivalent water height anomaly'
    apriori_ewh = np.vstack(apriori_ewh)
    apriori_ewh_arr[:]=apriori_ewh[:,:]
    #  sigma
    sigma_arr = ds.createVariable('sigma', 'f4', (('num_msc','num_epochs')))
    sigma_arr.units = '(m)'
    sigma_arr.description='formal uncertainty of mascon EWH'
    sigmas = np.vstack(sigmas)
    sigma_arr[:]=sigmas[:,:]

    # Close NetCDF file
    ds.close()

    ncname = outfile[:-3]+"_AuScope.nc"
    print("Writing old-style netcdf file for AuScope Portal: ",ncname)

    ncw=Dataset(ncname,'w',format='NETCDF4')
    ncw.description = 'ANU mascon solutions for AuScope Portal'
    # Create dimensions
    ncw.createDimension('nbt', num_epochs)
    ncw.createDimension('nbp', num_msc)
    # Create variables
    # decimal year
    decyear_vec= ncw.createVariable('decyear', 'f4', ('nbt'))
    decyear_vec.units = 'decimal year'
    decyear_vec.description='date in decimal year'
    decyear_vec[:]=dec_yr[:]
    # primary mascon
    pcode_vec= ncw.createVariable('primary', 'f4', ('nbp'))
    pcode_vec.units = '-'
    pcode_vec.description='primary mascon number'
    pcode_vec[:]=mascon_number[:]
    # primary latitude
    lat_vec= ncw.createVariable('lat', 'f4', ('nbp'))
    lat_vec.units = 'degrees'
    lat_vec.description='latitude of the primary mascon centroid in WGS84'
    lat_vec[:]=lat[:]
    # primary longitude
    lon_vec= ncw.createVariable('lon', 'f4', ('nbp'))
    lon_vec.units = 'degrees'
    lon_vec.description='longitude of the primary mascon centroid in WGS84'
    lon_vec[:]=lon[:]
    # drainage basin
    region_vec= ncw.createVariable('drainage_basin', 'S1', ('nbp'))
    region_vec.units = '-'
    region_vec.description='name of the drainage basin'
    region_vec[:]=region[:]
    # australia ewh
    ewh_arr= ncw.createVariable('ewh', 'f4', (('nbp','nbt')))
    ewh_arr.units = 'm'
    ewh_arr.description='equivalent water height'
    ewh = np.vstack(ewh)
    ewh_arr[:]=ewh[:,:]
    # australia sigma
    sigma_arr= ncw.createVariable('sigma', 'f4', (('nbp','nbt')))
    sigma_arr.units = '(m)'
    sigma_arr.description='uncertainty at one sigma'
    sigmas = np.vstack(sigmas)
    sigma_arr[:]=sigmas[:,:]
    # Close NetCF file
    ncw.close()







