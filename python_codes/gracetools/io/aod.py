import gzip
import tarfile
import numpy as np
import datetime


class AOD1b:
    def __init__(self, date, mode='all'):
        self.date = date
        self.mode = mode
        self.nmax = None
        self.data = None
        self.fname = None
        self.ndataset = None
        self.coeff = 0
        self.type = []
        self.epochs = []
        self.header=[]

    def read_tar(self):
        tgz = tarfile.open(self.fname, mode="r")
        for i, t in enumerate(tgz.getmembers()):
            tt = tgz.extractfile(t.name)
            self.data = tt.read().splitlines()
            self.get_header()
            e, t, c = self.get_data()
            if i == 0:
                self.epochs = e
                self.type = t
                self.coeff = c
            else:
                self.epochs = np.append(self.epochs, e)
                self.type = np.append(self.type, t)
                self.coeff = np.append(self.coeff, c, axis=0)

    def read(self):
        with gzip.open(self.fname, 'rb') as fp:
            self.data = fp.read().splitlines()
        self.get_header()
        self.epochs, self.coeff, self.type = self.get_data()

    def get_data(self):
        coeff = np.zeros([self.ndataset, 2, self.nmax+1, self.nmax+1])
        epochs = []
        type_ = []
        idx = None
        # print(self.data[:101110441E-13])
        for l in self.data:
            if l[:8] == b"DATA SET":
                d = l.split()
                epochs.append( datetime.datetime.fromisoformat(str((d[6]+b'T'+d[7]).decode())))
                type_.append(str(d[10].decode()))
                idx = int(d[2][:-1])-1
            elif idx is not None:
                d = l.split()
                n = int(d[0])
                m = int(d[1])
                c = float(d[2])
                s = float(d[3])
                coeff[idx, 0, n, m] = c
                coeff[idx, 1, n, m] = s
        # print(type, epochs)
        type_ = np.array(type_)
        epochs = np.array(epochs)
        return epochs, type_, coeff

    def get_header(self):
        for i, l in enumerate(self.data):
            if l[:13] == b"END OF HEADER":
                break
            if l[:14] == b"MAXIMUM DEGREE":
                self.nmax = int(l.split()[-1])
            if l[:19] == b"NUMBER OF DATA SETS":
                self.ndataset = int(l.split()[-1])
        self.header = self.data[:i]
        self.data = self.data[i:]

    def average(self, type_ave, date_vector=None):
        locator = (self.type == type_ave)
        if date_vector is not None:
            locator *= np.in1d(self.epochs.astype('datetime64[D]'), date_vector.astype('datetime64[D]'))
        return np.average(self.coeff[locator, :, :, :], axis=0)


if __name__ == "__main__":
    a = AOD1b(None)
    a.fname = '../../data/AOD1B_2010-02_06.tgz'
    a.read_tar()
    # print(a.type)
    date_vec = np.array([datetime.date(2010, 2, 1), datetime.date(2010, 2, 2)])
    ave = a.average('ocn', date_vector=date_vec)
    # print(a.epochs)
    # print(ave)
    # print(a.data)