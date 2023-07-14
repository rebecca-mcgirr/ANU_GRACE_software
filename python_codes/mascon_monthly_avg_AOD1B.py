#!/usr/bin/env python3
import numpy as np
import h5py
import pyshtools as sht
from datetime import datetime
from gracetools.io.mascons import mascons
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--mascons', default='/Mdata/GRACE/software/tables/mascons_stage5_V006_200km.h5', type=str, help='mascons files')
parser.add_argument('--aod', default='/Mdata/GRACE/L1B/AOD1B/RL06/monthly_averaged_AOD1B/AOD1B_<YYYY>-<MM>_oba.sph.gz', type=str, help='AOD1B files')
parser.add_argument('--start_date', type=str, help='YYYY-MM')
parser.add_argument('--end_date', type=str, help='YYYY-MM')
parser.add_argument('--love', default='/Mdata/GRACE/software/tables/data/Load_Love2_CM.dat', type=str, help='love number file')
parser.add_argument('--outfile', type=str, help='output h5 file')
args = parser.parse_args()

print('python3 ~/ga/python_codes/mascon_monthly_avg_AOD1B.py --mascons /Mdata/GRACE/software/tables/mascons_stage5_V006_200km.h5 --start_date 2002-04 --end_date 2022-09 --aod "/Mdata/GRACE/L1B/AOD1B/RL06/monthly_averaged_AOD1B/AOD1B_<YYYY>-<MM>_oba.sph.gz" --love /Mdata/GRACE/software/tables/data/Load_Love2_CM.dat --outfile /Mdata/GRACE/software/tables/data/mascons_stage5_V006_200km_AOD1B_oba.h5') 

print(f"reading {args.mascons}")
mANU = mascons(args.mascons)
mANU.read_mascons_h5()
mANU.get_mapping_sht()

# convert start and end date to list of years and months
sd = datetime.strptime(args.start_date, "%Y-%m") 
ed = datetime.strptime(args.end_date, "%Y-%m") 

dates = [datetime.strptime('%2.2d-%2.2d' % (y, m), '%Y-%m').strftime('%Y %m') \
       for y in np.arange(sd.year, ed.year+1) \
       for m in np.arange(sd.month if y==sd.year else 1, ed.month+1 if y == ed.year else 13)]

year = [int(d.split(' ')[0]) for d in dates]
month = [int(d.split(' ')[1]) for d in dates]
decyear = [year[i] + int(datetime.strptime('%2.2d-%2.2d-15'%(year[i], month[i]), '%Y-%m-%d').strftime('%j'))/ \
           int(datetime.strptime('%2.2d-12-31'%year[i], '%Y-%m-%d').strftime('%j')) \
           for i in range(len(dates))]

degmax=180
n, h, l, k = np.loadtxt(args.love, unpack=True)
n = n[:degmax + 1]
h = h[:degmax + 1]
k = k[:degmax + 1]
l = l[:degmax + 1]
k[0] = 0.0
k[1] = 0.026
# PT220908: change Earth radius from 6371.0e3 to 6378.1363e3
Conv = 1 / 3.0 * 6378.1363e3 *5.515 * (2 * n + 1) / (1 + k)


print("Mascon monthly AOD1B corrections written to %s"%args.outfile)

correction_aod = np.zeros((len(year),len(mANU.Pn)))
ocean = mANU.Pdensity > 1001.0 #why only ocean?
s = np.sum(mANU.Parea[ocean])
for i, (y, m) in enumerate(zip(year, month)):
    fname = args.aod.replace('<YYYY>', "%04i"%y).replace('<MM>',"%02i"%m).replace('<TT>', 'oba')
    print(fname)
    C = sht.SHCoeffs.from_file(fname, header=True, header2=True)
    #print("what are these C coeffs?",C.coeffs[0,0:5,0:5])
    C.coeffs = C.coeffs * Conv[None, :, None]
    #print("shape of C.coeffs:",np.shape(C.coeffs))
    correction_aod[i, :] = sht.expand.MakeGridPoint(C.coeffs, mANU.Plat, mANU.Plon, norm=1, csphase=1)
    x = np.sum(correction_aod[i, ocean] * mANU.Parea[ocean]) / s
    print(f" {y:04d}-{m:02d}  => {x: .4f} m")

# build h5 file
with h5py.File(args.outfile, 'w') as hf:
    pst = hf.create_dataset('correction/aod', data=correction_aod)
    tt = hf.create_dataset('time/decyear', data=decyear)
    ty = hf.create_dataset('time/year', data=year)
    tm = hf.create_dataset('time/month', data=month)
    hf.create_dataset('lat', data=mANU.Plat)
    hf.create_dataset('lon', data=mANU.Plon)
    hf.create_dataset('density', data=mANU.Pdensity)
    hf.create_dataset('area', data=mANU.Parea)

print("Mascon monthly AOD1B corrections written to %s"%args.outfile)

