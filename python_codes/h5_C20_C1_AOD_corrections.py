#!/usr/bin/env python3
# RM221122: moved all imports to beginning of script
import numpy as np
import pyshtools as sht
import datetime
import h5py
import argparse
import concurrent.futures
# PT220307: import the library of functions to extend time series
import extend_lib as ex
from gracetools.io.fit import Fit
from gracetools.io.aod import AOD1b
from gracetools.io.mascons import mascons
from netCDF4 import Dataset
from scipy.interpolate import interp1d
from pathlib import Path
# RM221122: import grace_utils library
from grace_utils import *

if __name__ == "__main__":

    runstr = 'python3 ~/ga/python_codes/h5_C20_C1_AOD_corrections.py --mascons /Mdata/GRACE/software/tables/data/mascons_stage5_V006_200km.h5 --tn14 /Mdata/GRACE/software/tables/data/TN13_SLR.txt --h5 b8i3_uncorrected.h5 --gia /Mdata/GRACE/software/tables/data/ICE6G_D_Peltier.hs --love /Mdata/GRACE/software/tables/data/Load_Love2_CM.dat --aod "/Mdata/GRACE/L1B/AOD1B/RL06/monthly_averaged_AOD1B/AOD1B_<YYYY>-<MM>_oba.sph.gz" --outfile b8i3_corrected'    

    parser = argparse.ArgumentParser(description=runstr)
    parser.add_argument('h5', default=None, type=str, help="combined hdf5 file")
    parser.add_argument('--mascons', type=str, default='/Mdata/GRACE/software/tables/data/mascons_stage5_V006_200km.h5', help='mascon files')
    #parser.add_argument('--tn14', type=str, default='/Mdata/GRACE/software/tables/data/TN-14_C30_C20_GSFC_SLR.txt', help='SLR time series')
    parser.add_argument('--aod', type=str, default='/Mdata/GRACE/software/tables/data/mascons_stage5_V006_200km_AOD1B_oba.h5', help='path to AOD1B files or output h5 file from mascon_monthly_avg_AOD1B.py')
    parser.add_argument('--gia', type=str, default='/Mdata/GRACE/software/tables/data/ICE6G_D_Peltier.hs', help='gia coeff file')
    parser.add_argument('--love', type=str, default='/Mdata/GRACE/software/tables/data/Load_Love2_CM.dat', help='love number file')
    parser.add_argument('--outfile', type=str,help='output stem')
    args = parser.parse_args()

    print(f"reading {args.mascons}")
    mANU = mascons(args.mascons)
    mANU.read_mascons_h5()
    mANU.get_mapping_sht()
    
    print("input solution hdf5 file: ",args.h5)
    with h5py.File(args.h5, "r") as h5f:
        apriori_ewh = h5f['solution/prefit'][:]
        sol = h5f['solution/ewh'][:]
        sigma_ewh = h5f['solution/sigma_ewh'][:]
        year = h5f['time/year'][:]
        decyear = h5f['time/decyear'][:]
        month = h5f['time/month'][:]

    degmax=180
    n, h, l, k = np.loadtxt(args.love, unpack=True)
    n = n[:degmax + 1]
    h = h[:degmax + 1]
    k = k[:degmax + 1]
    l = l[:degmax + 1]
    k[0] = 0.0
    k[1] = 0.026
    # PT220908: change Earth radius from 6371.0e3 to 6378.1363e3
    Conv = 1 / 3.0 * 6378.1363e3 * 5.515 * (2 * n + 1) / (1 + k)

    print('Calculate corrections:')
    print('1. Low degree harmonics')
    tn14 = np.loadtxt('/Mdata/GRACE/software/tables/data/TN-14_C30_C20_GSFC_SLR.txt', skiprows=37)
    # jpl = np.loadtxt('files/TN-11_C20_SLR.txt')
    # PT220216: -4.841652961880e-04 is the value of C20 in GOCE.deg200, -4.841694573200E-04 is the mean value in GSFC C20_C30 file.
    #            The value -4.179e-9 is the permanent tide contribution
    # PT220217: the C20_static value we should use is our C20 plus the permanent tide contribution. That is, it should be
    #           -4.841652961880e-04 +  -4.1376e-9  = -4.84169433788e-4
    DC20 = -4.841652961880e-04+4.841694573200E-04-4.179e-9
    # this was Seb's value. Don't know where it came from
    #C20_static = -4.8416938905481E-04
    C20_GOCE_deg200 = -4.841652961880e-04
    permanent_tide_contribution = -4.1376e-9
    C20_static = C20_GOCE_deg200 + permanent_tide_contribution
    
    # PT220307: extrapolate the time series before interpolating. We need to do this for GRACE-FO solutions that
    #           are available prior to the TN14 becoming available
    #  First, derive the lowpass filter series. Do this on the C20 anomalies
    C20_anomaly = (tn14[:,2] - np.mean(tn14[:,2])) 
    C30_anomaly = (tn14[114:,5] - np.mean(tn14[114:,5]))
    lowpass_series = ex.calc_lowpass(tn14[:,1],C20_anomaly,60,0.5) + np.mean(tn14[:,2])
    #  Next, derive a sinusoidal model from all the data
    ann_amplitudes = ex.calc_annual_signal(tn14[:,1],C20_anomaly)
    print("  Amplitudes of cosine/sine annual C20 signal:",ann_amplitudes)
    # extend the time series, using the lowpass_series and the annual amplitude
    dec_yr_extended_C20, C20_series_extended = ex.extend_timeseries(5,1/12,tn14[:,1],tn14[:,2],lowpass_series,ann_amplitudes)
    print("  TN14 extrapolated from [%.2f,%.2f] to [%.2f,%.2f]"%(tn14[0,-1],tn14[-1,-1],dec_yr_extended_C20[0],dec_yr_extended_C20[-1]))
    # RM221122: HOW FAR IS THIS BEING INTERPOLATED!!???
    # now, interpolate it
    tn14_f  = interp1d(dec_yr_extended_C20, C20_series_extended, fill_value=0.0, bounds_error=False)
    
    # PT220217: add a C30 correction for the GRACE-FO solutions
    #tn14_C30_f = interp1d(tn14[:, 1], tn14[:, 5], fill_value=0.0, bounds_error=False) 
    
    #deg 1
    tn13 = True
    if tn13:
        tn13_data = np.loadtxt('/Mdata/GRACE/software/tables/data/TN13_CSR')
        deg1_year = np.zeros(tn13_data.shape[0])
        C10 = np.zeros(tn13_data.shape[0])
        C11 = np.zeros(tn13_data.shape[0])
        S11 = np.zeros(tn13_data.shape[0])
        dec_yr_deg1 = np.zeros(tn13_data.shape[0])
        
        for i, l in enumerate(tn13_data):
            st = str(l[0])
            y_ = np.int_(st[:4])
            m_ = np.int_(st[4:6])
            d_ = np.int_(st[6:8])
            st2 = str(l[1])
            y_2 = np.int_(st2[:4])
            m_2 = np.int_(st2[4:6])
            d_2 = np.int_(st2[6:8])
            deg1_year[i] = 0.5* (toYearFraction(datetime.date(y_2, m_2, d_2)) + toYearFraction(datetime.date(y_, m_, d_)))
            C10[i] = l[2]
            C11[i] = l[3]
            S11[i] = l[4]
            
    else:
        data_deg1 = np.loadtxt('/Mdata/GRACE/software/tables/data/GCN_L1_L2_30d_CF-CM.txt', skiprows=1)
        rad_e = 6378136360
        deg1_year = data_deg1[:,0]
        C10 = data_deg1[:,3]/(rad_e*np.sqrt(3.0))
        C11 = data_deg1[:,1]/(rad_e*np.sqrt(3.0))
        S11 = data_deg1[:,2]/(rad_e*np.sqrt(3.0))
        # k1 here is the k^'_{1) of the CF frame, which has a value of 0.021 (Blewitt, 2003)
        k1 = 0.021
        for i in range(len(deg1_year)):
            print(deg1_year[i], C10[i], C11[i], S11[i])
        coef_d1 = 1/3.0*6378136.46*5.513*(2*1+1)/(1+k1)
    C10 *= Conv[1]
    C11 *= Conv[1]
    S11 *= Conv[1]

    # PT220307: extrapolate the degree-1 time series so that we can use them for recent GRACE-FO data
    # now, extrapolate for each of the coefficients C10, C11, S11
    lowpass_series_C10 = ex.calc_lowpass(deg1_year, C10 - np.mean(tn13_data[:,2]),60,0.5) + np.mean(tn13_data[:,2])
    lowpass_series_C11 = ex.calc_lowpass(deg1_year, C11 - np.mean(tn13_data[:,3]),60,0.5) + np.mean(tn13_data[:,3])
    lowpass_series_S11 = ex.calc_lowpass(deg1_year, S11 - np.mean(tn13_data[:,4]),60,0.5) + np.mean(tn13_data[:,4])

    ann_amplitudes_C10 = ex.calc_annual_signal(deg1_year,C10 - np.mean(tn13_data[:,2]))
    ann_amplitudes_C11 = ex.calc_annual_signal(deg1_year,C11 - np.mean(tn13_data[:,3]))
    ann_amplitudes_S11 = ex.calc_annual_signal(deg1_year,S11 - np.mean(tn13_data[:,4]))

    deg1_yr_extended, C10_series_extended = ex.extend_timeseries(5,1/12,deg1_year,C10,lowpass_series_C10,ann_amplitudes_C10)
    deg1_yr_extended, C11_series_extended = ex.extend_timeseries(5,1/12,deg1_year,C11,lowpass_series_C11,ann_amplitudes_C11)
    deg1_yr_extended, S11_series_extended = ex.extend_timeseries(5,1/12,deg1_year,S11,lowpass_series_S11,ann_amplitudes_S11)

    # now interpolate
    # Seb's code
    #C10_interpolator = interp1d(deg1_year, C10, fill_value=0.0, bounds_error=False)
    #C11_interpolator = interp1d(deg1_year, C11, fill_value=0.0, bounds_error=False)
    #S11_interpolator = interp1d(deg1_year, S11, fill_value=0.0, bounds_error=False)
    # using interpolated vectors
    C10_interpolator = interp1d(deg1_yr_extended, C10_series_extended, fill_value=0.0, bounds_error=False)
    C11_interpolator = interp1d(deg1_yr_extended, C11_series_extended, fill_value=0.0, bounds_error=False)
    S11_interpolator = interp1d(deg1_yr_extended, S11_series_extended, fill_value=0.0, bounds_error=False)

    print("  Rate of ICE6G_D C20 = ",1.38173003E-11 *Conv[2], Conv[2])

    # PT220906: break the corrections into a C20_mascon (to be subtracted), C20_SLR (to be added), deg1_correction (to be added)
    C20_mascon = np.zeros(len(sol[:,0]))
    C20_mascon_ewh = np.zeros_like(sol)
    C20_TN14_corr = np.zeros_like(sol)
    
    # PT221017: a whole pile more variables, just to confuse everyone later on
    C20_mascon_EWH = np.zeros_like(sol)
    C20_mascon_stokes_anomaly = np.zeros_like(sol)
    C20_mascon_EWH_anomaly = np.zeros_like(sol)


    deg1_correction = np.zeros_like(sol)
    for i, (dec, y, m) in enumerate(zip(decyear, year, month)):
        d = sol[i, :]
        rem_corr = sht.expand.SHExpandDH(d[mANU.idx_near], sampling=2, csphase=1,norm=1, lmax_calc=2)
        corr_C20 = np.zeros_like(rem_corr)
        corr_C20[0][2][0] = (tn14_f(dec)-C20_static)*Conv[2]

        corr_deg1 = np.zeros_like(rem_corr)
        corr_deg1[0][2][0] = 0.0
        corr_deg1[0][1][0] = C10_interpolator(dec) #- C.coeffs[0,1,0] * coef_d1
        corr_deg1[0][1][1] = C11_interpolator(dec) #- C.coeffs[0,1,1] * coef_d1
        corr_deg1[1][1][1] = S11_interpolator(dec) #- C.coeffs[1,1,1] * coef_d1

        c20 = rem_corr[0,2,0] #   - 0.004  # PT220908: subtracting 0.004 aligns it with the TN14 time series. Not sure why it is necessary!
        rem_corr[:,:,:]  = 0.0
        rem_corr[0,2,0] = c20
        #####print(i," c20",c20,Conv[2])
        C20_mascon[i] = c20
        # PT220906: calculate the C20 per mascon using the C20 estimated from our mascon fields
        C20_mascon_ewh[i, :] = sht.expand.MakeGridPoint(rem_corr, mANU.Plat, mANU.Plon, csphase=1, norm=1)
        # PT220906: calculate the C20 per mascon using the TN14 C20 coefficients
        C20_TN14_corr[i, :]   = sht.expand.MakeGridPoint(corr_C20, mANU.Plat, mANU.Plon, csphase=1, norm=1)
        # PT220906: calculate the C10, C11, S11 per mascon from the TN13 coefficient values
        deg1_correction[i, :] += sht.expand.MakeGridPoint(corr_deg1, mANU.Plat, mANU.Plon, csphase=1, norm=1)


        #############
        ## PT221017: start again!!!!
        #            delta_C20 wrt TN14 mean =  delta_C20^mascons + C20^GOCE - C20^TN14
        
        # save the C20 coefficient, converted to a dimensionless Stoke's coefficeint plus the static GOCE.deg200 C20 value plus the permanent tide MINUS the TN14 mean value
        mascon_stokes_C20 = np.zeros_like(rem_corr)
        mascon_stokes_C20[0,2,0] = rem_corr[0,2,0] /Conv[2] + C20_static - (-4.841694573200E-04) 
        #####print("mascon_stokes_C20: ",mascon_stokes_C20[0,2,0])
        # calculate what signal this coefficient puts on each mascon 
        C20_mascon_stokes_anomaly[i, :] = sht.expand.MakeGridPoint(mascon_stokes_C20, mANU.Plat, mANU.Plon, csphase=1, norm=1) 
        C20_mascon_EWH_anomaly[i, :] = sht.expand.MakeGridPoint(mascon_stokes_C20*Conv[2], mANU.Plat, mANU.Plon, csphase=1, norm=1) 
        #####print("evaluated C20 for mascon 100:",C20_mascon_stokes_anomaly[i, 100])
        
    with h5py.File('d12_correction.h5', 'w') as hf:
        hf.create_dataset('correction/C20_mascon', data=C20_mascon)
        hf.create_dataset('correction/C20_mascon_ewh', data=C20_mascon_ewh)
        hf.create_dataset('correction/C20_TN14_corr', data=C20_TN14_corr)
        hf.create_dataset('correction/deg1', data=deg1_correction)
        hf.create_dataset('time/decyear', data=decyear)
        hf.create_dataset('time/year', data=year)
        hf.create_dataset('time/month', data=month)

    # PT220307: Seb's original code just computed GIA correction from the first epoch. Change it to always do it from
    #           a fixed epoch. Let's make it from the epoch of the static gravity field, which is ~2008.0
    print('2. GIA')
    GIA_t0 = 2008.0
    print('  GIA (t_0 = ',GIA_t0,')')
    mANU.compute_GIA(fname=args.gia, love_fn=args.love)
    #gia_correct =  mANU.P_GIA[np.newaxis, :] * (decyear[:, np.newaxis]-decyear[0])
    gia_correct =  mANU.P_GIA[np.newaxis, :] * (decyear[:, np.newaxis]-GIA_t0)

    with h5py.File('GIA_correction.h5', 'w') as hf:
        pst = hf.create_dataset('correction/gia', data=gia_correct)
        tt = hf.create_dataset('time/decyear', data=decyear)
        ty = hf.create_dataset('time/year', data=year)
        tm = hf.create_dataset('time/month', data=month)
        hf.create_dataset('lat', data=mANU.Plat)
        hf.create_dataset('lon', data=mANU.Plon)
        hf.create_dataset('density', data=mANU.Pdensity)
        hf.create_dataset('area', data=mANU.Parea)

    print('3. Non-tidal ocean (OBA or GAD)')
    correction_aod = np.zeros_like(sol)
    ocean = mANU.Pdensity > 1001.0
    s = np.sum(mANU.Parea[ocean])
    fp_out = open("aod_oce.txt","w")
    if args.aod[-7:] == ".sph.gz" :
        for i, (y, m) in enumerate(zip(year, month)):
            fname = args.aod.replace('<YYYY>', "%04i"%y).replace('<MM>',"%02i"%m).replace('<TT>', 'oba')
            print(fname)
            C = sht.SHCoeffs.from_file(fname, header=True, header2=True)
            #print("what are these C coeffs?",C.coeffs[0,0:5,0:5])
            C.coeffs = C.coeffs * Conv[None, :, None]
            #print("shape of C.coeffs:",np.shape(C.coeffs))
            correction_aod[i, :] = sht.expand.MakeGridPoint(C.coeffs, mANU.Plat, mANU.Plon, norm=1, csphase=1)
            x = np.sum(correction_aod[i, ocean] * mANU.Parea[ocean]) / s
            #print(f"  {y:04d}-{m:02d}  => {x: .4f} m")
            fp_out.write(f" {y:04d}-{m:02d}-15  {x: .4f} \n")
            print(i, y, m, fname)
        fp_out.close()
        with h5py.File('AOD_correction.h5', 'w') as hf:
            pst = hf.create_dataset('correction/aod', data=correction_aod)
            tt = hf.create_dataset('time/decyear', data=decyear)
            ty = hf.create_dataset('time/year', data=year)
            tm = hf.create_dataset('time/month', data=month)
            hf.create_dataset('lat', data=mANU.Plat)
            hf.create_dataset('lon', data=mANU.Plon)
            hf.create_dataset('density', data=mANU.Pdensity)
            hf.create_dataset('area', data=mANU.Parea)
    # RM221121: modified code to read .h5 and find correct dates
    elif args.aod[-3:] == ".h5" :
        h5f = h5py.File(args.aod, "r")
        if len(mANU.Pn) != len(h5f['lat']): 
    	    raise Exception("AOD file incompatible with solution file. Were they created with different mascon geometries?")
        year_aod = h5f['time/year'][:]
        month_aod = h5f['time/month'][:]
        j = 0
        for i, (y, m) in enumerate(zip(year_aod, month_aod)):
            if y == year[j] and m == month[j]:
                correction_aod[j,:] = h5f['correction/aod'][i,:]
                x = np.sum(correction_aod[j, ocean] * mANU.Parea[ocean]) / s
                #print(f" {y:04d}-{m:02d}  => {x: .4f} m")
                fp_out.write(f" {y:04d}-{m:02d}-15  {x: .4f} \n")
                j += 1
            if j == len(month): break
        fp_out.close()
        print('  Timeseries of AOD integrated over ocean mascons written to "aod_oce.txt"')
        print("Have read AOD1B corrections from file %s"%args.aod)
               
    output_h5 = args.outfile + ".h5"
    output_nc = args.outfile + ".nc"

    with h5py.File(output_h5, "w") as hf:
        hf.create_dataset('solution/prefit', data = apriori_ewh + correction_aod - gia_correct - C20_mascon_ewh + C20_TN14_corr + deg1_correction) 
        hf.create_dataset('solution/ewh', data = sol + correction_aod - gia_correct - C20_mascon_ewh + C20_TN14_corr + deg1_correction) 
        hf.create_dataset('solution/sigma_ewh',data = np.vstack(np.transpose(sigma_ewh)) )

        # C20 signal calculated on each mascon
        hf.create_dataset('solution/C20_mascon_EWH_anomaly', data = C20_mascon_EWH_anomaly + correction_aod)
        # C20 values estimated from mascons
        hf.create_dataset('solution/C20_mascon', data=C20_mascon)

        hf.create_dataset('time/year', data=year)
        hf.create_dataset('time/month', data=month)
        hf.create_dataset('time/decyear', data=decyear)
                
        print("Write corrections (C20, GIA, OBA) to output file: ",output_h5)
        # PT220908: also write all corrections out to one, single hdf5 file
        # C20 signal calculated on each mascon, in dimensionless Stoke's values
        hf.create_dataset('correction/C20_mascon_stokes_anomaly', data=C20_mascon_stokes_anomaly)
        # C20 signal from TN14, calculated on each mascon
        hf.create_dataset('correction/C20_TN14_corr', data=C20_TN14_corr)
        # deg1 TN13_CSR calculated on each mascon
        hf.create_dataset('correction/deg1', data=deg1_correction)
        # non-tidal ocean (OBA, or "GAD") on each mascon. Should have land = zero
        hf.create_dataset('correction/aod', data=correction_aod)
        # GIA correction. Mean date of 2008.0 used
        hf.create_dataset('correction/gia', data=gia_correct)

        print("Add some mascon information as well (lat, lon, density, area)")
        hf.create_dataset('mascon/Plat', data=mANU.Plat)
        hf.create_dataset('mascon/Plon', data=mANU.Plon)
        hf.create_dataset('mascon/Pdensity', data=mANU.Pdensity)
        hf.create_dataset('mascon/Parea', data=mANU.Parea)
            
    # PT220201: also write out the fully corrected version in netcdf format
    ds=Dataset(output_nc,'w',format='NETCDF4')
    ds.description = 'ANU mascon solutions for 2003-2016 (Allgeyer et al, 2022)'
        
    # Create dimensions
    num_msc=len(mANU.Pn[:])
    num_epochs = len(decyear)
    ds.createDimension('num_epochs', num_epochs)
    ds.createDimension('num_msc', num_msc)

    # PT220201: for the netcdf file the ewh and sigma arrays have the columns reversed
    ewh = np.zeros((num_msc,num_epochs))
    ewh = np.transpose(sol + correction_aod - gia_correct - C20_mascon_ewh + C20_TN14_corr + deg1_correction)
    apr_ewh = np.zeros((num_msc,num_epochs))
    apr_ewh = np.transpose(apriori_ewh + correction_aod - gia_correct - C20_mascon_ewh + C20_TN14_corr + deg1_correction)
    sigmas = np.zeros((num_msc,num_epochs))
    sigmas = np.transpose(sigma_ewh)
    
    print("Writing netcdf file")
    # Create variables
 
    # decimal year
    decyear_vec= ds.createVariable('decyear', 'f4', ('num_epochs'))
    decyear_vec.units = 'decimal year'
    decyear_vec.description='Epoch of monthly estimate in decimal year'
    decyear_vec[:]=decyear[:]
 
    # primary mascon
    pcode_vec= ds.createVariable('mascon_number', 'f4', ('num_msc'))
    pcode_vec.units = '-'
    pcode_vec.description='primary mascon number'
    pcode_vec[:]=mANU.Pn[:]
 
    # primary latitude
    lat_vec= ds.createVariable('lat', 'f4', ('num_msc'))
    lat_vec.units = 'degrees'
    lat_vec.description='latitude of the primary mascon centroid in WGS84'
    lat_vec[:]=mANU.Plat[:]
 
    # primary longitude
    lon_vec= ds.createVariable('lon', 'f4', ('num_msc'))
    lon_vec.units = 'degrees'
    lon_vec.description='longitude of the primary mascon centroid in WGS84'
    lon_vec[:]=mANU.Plon[:]
 
    # drainage basin
    region_vec= ds.createVariable('drainage_basin', 'S1', ('num_msc'))
    region_vec.units = '-'
    region_vec.description='name of the drainage basin'
    region_vec[:]=mANU.Pdesc[:]
 
    #  ewh
    ewh_arr= ds.createVariable('ewh', 'f4', (('num_msc','num_epochs')))
    ewh_arr.units = 'm'
    ewh_arr.description='estimated mascon equivalent water height anomaly'
    ewh = np.vstack(ewh)
    ewh_arr[:]=ewh[:,:]
 
    #  apriori ewh
    apr_ewh_arr= ds.createVariable('apriori_ewh', 'f4', (('num_msc','num_epochs')))
    apr_ewh_arr.units = 'm'
    apr_ewh_arr.description='apriori mascon equivalent water height anomaly'
    apr_ewh = np.vstack(apr_ewh)
    apr_ewh_arr[:]=apr_ewh[:,:]
 
    #  sigma
    sigma_arr= ds.createVariable('sigma', 'f4', (('num_msc','num_epochs')))
    sigma_arr.units = '(m)'
    sigma_arr.description='formal uncertainty of mascon EWH'
    sigma_ewh = np.vstack(sigmas)
    sigma_arr[:]=sigma_ewh[:,:]
 
    # Close NetCF file
    ds.close()
    
