from collections import defaultdict
from itertools import count
from os.path import splitext
import numpy as np
from h5py import File
from netCDF4 import Dataset
from pyshtools.expand import MakeGridPoint
from scipy.spatial import cKDTree


def lon_lat_to_cartesian(lon, lat, R=1):
    """
    calculates lon, lat coordinates of a point on a sphere with
    radius R
    """
    lon_r = np.radians(lon)
    lat_r = np.radians(lat)

    x = R * np.cos(lat_r) * np.cos(lon_r)
    y = R * np.cos(lat_r) * np.sin(lon_r)
    z = R * np.sin(lat_r)
    return x, y, z

def intg_ewh(self,mask,units='EWH',demean=True):
    '''
    units: EWH/GT/ESL
    '''
    intg_series = np.zeros(self.nepoch)
    for iepoch in range(self.nepoch):
        if units == 'EWH':
            intg_series[iepoch] = np.sum(self.ewh[iepoch] * self.area * mask) / np.sum(self.area * mask)
        elif units == 'GT':
            intg_series[iepoch] = np.sum(self.ewh[iepoch] * self.area * mask) * 1e-9
        elif units == 'ESL':
            raise Exception("Units not coded. Exiting script.")
        else:
            raise Exception("Units not coded. Exiting script.")
    if demean:
        return intg_series - np.mean(intg_series[(self.date > 2004) * (self.date < 2011)])
    else:
        return intg_series
        
def load_mask(self,fname):
    '''
    Load the mask files that have been created 
    for each of the mascon solutions.

    ftype is 'nc' for CSR and JPL solutions or 
    'h5' for ANU and GSFC solutions.
    '''
    if fname[-2:] == 'nc':
        with Dataset(fname) as f:
            self.land_mask = f['land_mask'][:]
            self.ocean_mask = f['ocean_mask'][:]
            self.AIS_mask = f['antarctic_mask'][:]
            self.GIS_mask = f['greenland_mask'][:]
            self.MDB_mask = f['murrayDB_mask'][:]
            self.caspian_mask = f['caspian_mask'][:]
            self.amazon_mask = f['amazon_mask'][:]
            self.area = f['area'][:]
    elif fname[-2:] == 'h5':
        with File(fname,'r') as f:
            self.land_mask = f['mask/land'][:]
            self.ocean_mask = f['mask/ocean'][:]
            self.AIS_mask = f['mask/antarctic'][:]
            self.GIS_mask = f['mask/greenland'][:]
            self.MDB_mask = f['mask/murrayDB'][:]
            self.caspian_mask = f['mask/caspian'][:]
            self.amazon_mask = f['mask/amazon'][:]
    else:
        raise Exception("Error: file type '%s' not recognised"%fname[-2:])
        
class mascons:
    def __init__(self, fname):
        self._fname = fname
        self.code = None
        self.ntnmsc = None
        self.lat = None
        self.lon = None
        self.pnum = None
        self.lat_grd = None
        self.lon_grd = None
        self.idx_near = None
        self.P_GIA = None

    def read_mascons(self):
        if splitext(self._fname)[-1] == '.h5':
            self.read_mascons_h5()
        else:
            self.read_mascons_txt()

    def read_mascons_h5(self):
        with File(self._fname, 'r') as hf:
            self.Pn = hf['mascons/P/idx'][:]
            self.Plat = hf['mascons/P/lat'][:]
            self.Plon = hf['mascons/P/lon'][:]
            self.Pradius = hf['mascons/P/R'][:]
            self.Parea = hf['mascons/P/area'][:]
            self.Palt = hf['mascons/P/alt'][:]
            self.Pgeod = hf['mascons/P/geod'][:]
            self.Pdensity = hf['mascons/P/density'][:]
            # 'mascons/P/description_label', data=np.string_(list(Pdesclabel.keys())))
            self.Pdesc = hf['mascons/P/desc'][:]
            self.Pdesc_label = hf['mascons/P/description_label'][:]
            # hf.create_dataset('mascons/P/type_label', data=np.string_(list(Ptypelabel.keys())))
            # hf.create_dataset('mascons/P/type', data=Ptypemap - 1)

            self.Tnum = hf['mascons/T/idx'][:]
            self.Tlat = hf['mascons/T/lat'][:]
            self.Tlon = hf['mascons/T/lon'][:]
            self.Tradius = hf['mascons/T/R'][:]
            self.Tarea = hf['mascons/T/area'][:]
            self.Talt = hf['mascons/T/alt'][:]
            self.Tgeod = hf['mascons/T/geod'][:]
            self.Tdens = hf['mascons/T/density'][:]
            self.pnum = hf['mascons/T/primary'][:]
            self.desc = hf['mascons/T/desc'][:]
            # 'mascons/T/description_label', data=np.string_(list(Tdesclabel.keys())))
            # hf.create_dataset('mascons/T/desc', data=Tdescmap - 1)
            # hf.create_dataset('mascons/T/type_label', data=np.string_(list(Ttypelabel.keys())))
            # hf.create_dataset('mascons/T/type', data=Ttypemap - 1)

    def read_mascons_txt(self):
        with open(self._fname, 'r') as fp:
            ii = -1
            ip = -1
            for i, line in enumerate(fp):
                if i == 0:
                    a = line.split()
                    self.code = a[0][1:]
                    self.Tnum = np.empty(int(a[3]), dtype='i')
                    self.Ttype = np.empty(int(a[3]), dtype='U10')
                    self.Tlat = np.empty(int(a[3]))
                    self.Tlon = np.empty(int(a[3]))
                    self.Tradius = np.empty(int(a[3]))
                    self.Tarea = np.empty(int(a[3]))
                    self.Talt = np.empty(int(a[3]))
                    self.Tgeod = np.empty(int(a[3]))
                    self.Tdens = np.empty(int(a[3]))
                    self.pnum = np.empty(int(a[3]), dtype='i')
                    self.desc = np.empty(int(a[3]), dtype='U10')

                    self.Plat = np.empty(int(a[1]))
                    self.Plon = np.empty(int(a[1]))
                    self.Pn = np.empty(int(a[1]), dtype='i')
                    self.Pradius = np.empty(int(a[1]))
                    self.Palt = np.empty(int(a[1]))
                    self.Pgeod = np.empty(int(a[1]))
                    self.Parea = np.empty(int(a[1]))
                    self.Pdensity = np.empty(int(a[1]))
                    self.Pdesc = np.empty(int(a[1]), dtype='U10')
                    self.Ptype = np.empty(int(a[1]), dtype='U10')
                else:
                    if line[0] != '#':
                        a = line.split()
                        if a[1][0] == 'T':
                            ii = ii + 1
                            self.Tnum[ii] = int(a[0])
                            self.Ttype[ii] = a[1][:]
                            self.Tlat[ii] = float(a[2])
                            self.Tlon[ii] = float(a[3])
                            self.Tradius[ii] = float(a[4])
                            self.Tarea[ii] = float(a[5])
                            self.Talt[ii] = float(a[6])
                            self.Tgeod[ii] = float(a[7])
                            self.Tdens[ii] = float(a[8])
                            self.pnum[ii] = int(a[9])
                            self.desc[ii] = a[13][:]
                        if a[1][0] == 'P':
                            ip = ip + 1
                            self.Pn[ip] = int(a[0])
                            self.Ptype[ip] = a[1][:]
                            self.Plat[ip] = float(a[4])
                            self.Plon[ip] = float(a[5])
                            self.Pradius[ip] = float(a[6])
                            self.Parea[ip] = float(a[7])
                            self.Palt[ip] = float(a[8])
                            self.Pgeod[ip] = float(a[9])
                            self.Pdensity[ip] = float(a[10])
                            self.Pdesc[ip] = a[13][:]

    def get_mapping(self, spacing=30.0 / 60.0):
        self.lat_grd = np.r_[-90.0:90.0:spacing]
        self.lon_grd = np.r_[0:360:spacing]
        lon2d, lat2d = np.meshgrid(self.lon_grd, self.lat_grd)
        xG, yG, zG = lon_lat_to_cartesian(lon2d.flatten(), lat2d.flatten())
        xT, yT, zT = lon_lat_to_cartesian(self.Tlon, self.Tlat)
        tree = cKDTree(list(zip(xT, yT, zT)))
        d, inds = tree.query(list(zip(xG, yG, zG)), k=1)
        self.idx_near = np.int_(self.pnum[inds].reshape(lon2d.shape) - 1)

    def get_mapping_sht(self, spacing=0.5):
        self.lat_grd = np.r_[90.0:-90.0:-1 * spacing]
        self.lon_grd = np.r_[0:360:spacing]
        lon2d, lat2d = np.meshgrid(self.lon_grd, self.lat_grd)
        xG, yG, zG = lon_lat_to_cartesian(lon2d.flatten(), lat2d.flatten())
        xT, yT, zT = lon_lat_to_cartesian(self.Tlon, self.Tlat)
        tree = cKDTree(list(zip(xT, yT, zT)))
        d, inds = tree.query(list(zip(xG, yG, zG)), k=1)
        self.idx_near = np.int_(self.pnum[inds].reshape(lon2d.shape) - 1)

    def compute_GIA(self, fname='../input/ICE-6G_High_Res_Stockes_trend.txt',
                    maxdeg=180, love_fn='../input/Load_Love2_CM.dat'):
        coeff = np.zeros((2, maxdeg + 1, maxdeg + 1))
        with open(fname, 'r') as fp:
            for i, line in enumerate(fp):
                d = line.split()
                if int(d[0]) > maxdeg:
                    break
                coeff[0, int(d[0]), int(d[1])] = float(d[2])
                coeff[1, int(d[0]), int(d[1])] = float(d[3])

        n, h, l, k = np.loadtxt(love_fn, unpack=True)
        coef = np.zeros(maxdeg + 1)
        coef[2:] = 1 / 3.0 * 6378136.46 * 5.513 * (2 * n[2:maxdeg + 1] + 1) / (1 + k[2:maxdeg + 1])
        coeff[0, 2:, :] *= coef[2:, None]
        coeff[1, 2:, :] *= coef[2:, None]
        self.P_GIA = np.zeros(len(self.Plat))
        # for i in range(len(self.Plat)):
        self.P_GIA = MakeGridPoint(coeff, self.Plat, self.Plon, norm=1, csphase=1)

    def to_h5(self, fname):

        with File(fname, 'w') as hf:
            Pdesclabel = defaultdict(count(1).__next__)
            Pdescmap = np.array([Pdesclabel[k] for k in self.Pdesc], dtype=np.int_)

            Ptypelabel = defaultdict(count(1).__next__)
            Ptypemap = np.array([Ptypelabel[k] for k in self.Ptype], dtype=np.int_)

            Tdesclabel = defaultdict(count(1).__next__)
            Tdescmap = np.array([Tdesclabel[k] for k in self.desc], dtype=np.int_)

            Ttypelabel = defaultdict(count(1).__next__)
            Ttypemap = np.array([Ttypelabel[k] for k in self.Ttype], dtype=np.int_)

            hf.create_dataset('mascons/P/idx', data=self.Pn)
            hf.create_dataset('mascons/P/lat', data=self.Plat, dtype='f')
            hf.create_dataset('mascons/P/lon', data=self.Plon, dtype='f')
            hf.create_dataset('mascons/P/R', data=self.Pradius, dtype='f')
            hf.create_dataset('mascons/P/area', data=self.Parea, dtype='f')
            hf.create_dataset('mascons/P/alt', data=self.Palt, dtype='f')
            hf.create_dataset('mascons/P/geod', data=self.Pgeod, dtype='f')
            hf.create_dataset('mascons/P/density', data=self.Pdensity, dtype='f')
            hf.create_dataset('mascons/P/description_label', data=np.string_(list(Pdesclabel.keys())))
            hf.create_dataset('mascons/P/desc', data=Pdescmap - 1)
            hf.create_dataset('mascons/P/type_label', data=np.string_(list(Ptypelabel.keys())))
            hf.create_dataset('mascons/P/type', data=Ptypemap - 1)

            hf.create_dataset('mascons/T/idx', data=self.Tnum)
            hf.create_dataset('mascons/T/lat', data=self.Tlat, dtype='f')
            hf.create_dataset('mascons/T/lon', data=self.Tlon, dtype='f')
            hf.create_dataset('mascons/T/R', data=self.Tradius, dtype='f')
            hf.create_dataset('mascons/T/area', data=self.Tarea, dtype='f')
            hf.create_dataset('mascons/T/alt', data=self.Talt, dtype='f')
            hf.create_dataset('mascons/T/geod', data=self.Tgeod, dtype='f')
            hf.create_dataset('mascons/T/density', data=self.Tdens, dtype='f')
            hf.create_dataset('mascons/T/primary', data=self.pnum)
            hf.create_dataset('mascons/T/description_label', data=np.string_(list(Tdesclabel.keys())))
            hf.create_dataset('mascons/T/desc', data=Tdescmap - 1)
            hf.create_dataset('mascons/T/type_label', data=np.string_(list(Ttypelabel.keys())))
            hf.create_dataset('mascons/T/type', data=Ttypemap - 1)

class JPL:
    def __init__(self, fname=None):
        if fname is None:
            print("File name needed")
        else:
            self.fname = fname
            self.read()

    def read(self):
        with Dataset(self.fname) as ds:
            self.lon = ds['lon'][:]
            self.lat = ds['lat'][:]
            self.date = 2002 + ds['time'][:]/365.25 # convert days since 1 Jan 2002 to dec yr
            self.ewh = ds["lwe_thickness"][:]*1e-2 # convert cm to m
            self.nepoch = len(self.date)
            
class CSR:
    def __init__(self, fname=None):
        if fname is None:
            print("File name needed")
        else:
            self.fname = fname
            self.read()

    def read(self):
        with Dataset(self.fname) as ds:
            self.lon = ds['lon'][:]
            self.lat = ds['lat'][:]
            self.date = 2002 + ds['time'][:]/365.25 # convert days since 1 Jan 2002 to dec yr
            self.ewh = ds["lwe_thickness"][:]*1e-2 # convert cm to m
            self.nepoch = len(self.date)

class GSFC:
    def __init__(self, fname=None):
        self.lat = None
        self.lon = None
        self.date = None
        self.lat_grd = None
        self.lon_grd = None
        if fname is None:
            print("File name needed")
        else:
            self.fname = fname
            self.read()

    def read(self):
        with File(self.fname, 'r') as f:
            self.lat = f['mascon/lat_center'][0, :]
            self.lon = f['mascon/lon_center'][0, :]
            self.date = f['time/yyyy_doy_yrplot_middle'][2,:]
            self.ewh = f['solution/cmwe'][:] * 1e-2
            self.location = f['mascon/location'][0, :]
            self.basin = f['mascon/basin'][0, :]
            self.area = f['mascon/area_km2'][0, :] * 1e6
            self.nmsc = len(self.lat)
            self.nepoch = len(self.date)

    def compute_GIA(self, fname='ICE-6G_High_Res_Stockes_trend.txt', load_love_fname='Load_Love2_CM',
                    maxdeg=180):
        from pyshtools.expand import MakeGridPoint
        coeff = np.zeros((2, maxdeg + 1, maxdeg + 1))
        with open(fname, 'r') as fp:
            for i, line in enumerate(fp):
                d = line.split()
                if int(d[0]) > maxdeg:
                    break
                coeff[0, int(d[0]), int(d[1])] = float(d[2])
                coeff[1, int(d[0]), int(d[1])] = float(d[3])
        n, h, l, k = np.loadtxt(load_love_fname, unpack=True)
        coef = np.zeros(maxdeg + 1)
        coef[2:] = 1 / 3.0 * 6378136.46 * 5.513 * (2 * n[2:maxdeg + 1] + 1) / (1 + k[2:maxdeg + 1])
        coeff[0, 2:, :] *= coef[2:, None]
        coeff[1, 2:, :] *= coef[2:, None]
        self.P_GIA = np.zeros(len(self.lat))
        for i in range(len(self.lat)):
            self.P_GIA[i] = MakeGridPoint(coeff, self.lat[i], self.lon[i], norm=1, csphase=1)

    def select_time(self, idx):
        self.ewh = self.ewh  # - self.ewh.mean(axis=1)[:,np.newaxis]
        self.ewh = self.ewh[:, idx] - self.ewh.mean(axis=1)
        print('TIME: ', self.date[:, idx])

    def get_mapping(self, spacing=5 / 60.0):
        self.lat_grd = np.r_[-90:90.0:spacing]
        self.lon_grd = np.r_[0:360.:spacing]
        lon2d, lat2d = np.meshgrid(self.lon_grd, self.lat_grd)
        xG, yG, zG = lon_lat_to_cartesian(lon2d.flatten(), lat2d.flatten())
        xT, yT, zT = lon_lat_to_cartesian(self.lon, self.lat)
        self.tree = cKDTree(np.vstack([xT, yT, zT]).transpose())
        d, inds = self.tree.query(np.vstack([xG, yG, zG]).transpose(), k=1)
        return self.ewh[inds].reshape(lon2d.shape), self.lat_grd, self.lon_grd

    def project(self, lat, lon):
        xP, yP, zP = lon_lat_to_cartesian(lon, lat)
        d, inds = self.tree.query(np.vstack([xP, yP, zP]).transpose(), k=1)
        return self.ewh[inds]  # .reshape(lon2d.shape)

class ANU:
    def __init__(self, fname, corrected=False, apriori=False):
        self.fname = fname
        if fname[-2:] == 'h5':
            with File(self.fname, 'r') as f:
                self.lat = f['mascon/Plat'][:]
                self.lon = f['mascon/Plon'][:]
                self.date = f['time/decyear'][:]
                self.year = f['time/year'][:]
                self.month = f['time/month'][:]
                self.area = f['mascon/Parea'][:]
                self.density = f['mascon/Pdensity'][:]
                self.ewh = f['solution/ewh'][:,:]
                if apriori: self.ewh = f['solution/prefit'][:,:]
                if corrected: 
                    self.C20_ewh = f['solution/C20_mascon'][:]
                    self.gia = f['correction/gia'][:,:]
                    self.deg1 = f['correction/deg1'][:,:]
                    self.aod = f['correction/aod'][:,:]
                self.nmsc = len(self.lat)
                self.nepoch = len(self.date)
                
if __name__ == '__main__':
    mANU = mascons(fname='/scratch/paper_results/mascons_file/stage5_V005_200km/mascons_V005_200km.h5')
    mANU.read_mascons()
    mGSFC = GSFC(fname='/scratch/paper_results/model/GSFC.glb.200301_201607_v02.4.hdf')
    mGSFC.select_time(91)
    mGSFC.get_mapping()
    res = mGSFC.project(mANU.lat, mANU.lon)
    ewh = np.zeros((int(mANU.pnum.max())))
    nmsc = len(mANU.Pn)
    for i in range(1, nmsc + 1):
        ewh[i - 1] = res[mANU.pnum == i].mean()
    print(res)
    print(res.min(), res.max())
    with open(f'aprGRGS091.vcv', 'w') as fp:
        fp.write("V2 APRIORI  VCV  created for 2010-09 from GSFC mascon solution\n")
        fp.write("%i mascons in file\n" % (nmsc))
        fp.write("Output primary mascon geometry file: mascons_stage5_V004\n")
        fp.write("An incomplete VCV file, containing only mascon values to be read\n")
        fp.write("SOLUTION A PRIORI AND VECTOR:\n")
        fp.write("PARAMETER                     A PRIORI             VECTOR            SIGMA\n")
        for i in range(1, nmsc + 1):
            fp.write(f'{i + 24: 6d}.MC{i:05d} (m)           {0: 17.7f} {ewh[i - 1]: 19.9f}{0: 19.9f}\n')
