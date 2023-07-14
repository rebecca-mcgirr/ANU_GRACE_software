import numpy as np
import re

class Fit:
    def __init__(self, fname=None):
        self.addnorm = False
        self.date=None
        self.fname=fname
        self.n_mascons = 0
        self.length = 0
        self.msc_apriori = None
        self.msc_adjustment = None
        self.msc_posteriori = None

    def read(self):
        with open(self.fname) as fp:
            data = fp.readlines()
        if data[0][:10] == 'V2 ADDNORM':
            self.addnorm = True
            self.parse_addnorm(data)

    def parse_addnorm(self, data):
        line = data[1].split()
        self.n_mascons = int(line[-1])
        self.length = int(line[-5])
        self.date = np.empty((self.length), dtype='datetime64[D]')
        for i, line in enumerate(data[3:3+self.length]):
            yy = int(line.split()[3])
            mm = int(line.split()[4])
            dd = int(line.split()[5])
            self.date[i] = np.datetime64(f"{yy:04d}-{mm:02d}-{dd:02d}")
        select = [re.findall(r'^.*MC.*$', line) for line in data]
        split = [re.findall("[-+]?\d+[\.]?\d*", s[0]) for s in select if s != []]
        info = np.asarray(list(map(np.float32, split)))
        assert info.shape[0] == self.n_mascons, f"mascons information incompatible in {self.fname} between header and data"
        self.msc_apriori = info[:, 2]
        self.msc_adjustment = info[:, 3]
        self.msc_posteriori = info[:, 4]



if __name__ == "__main__":
    f = Fit()
    f.fname = "../../data/msc_2010-02_b7i0d.fit"
    f.read()