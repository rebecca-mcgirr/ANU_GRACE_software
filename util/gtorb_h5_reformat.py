#!/apps/python3/3.7.4/bin/python3
import h5py
import sys

if (len(sys.argv) != 3):
    print("Input arguments are old_GTORB.h5 new_GTORB.h5")
    sys.exit(2)

fname = sys.argv[1]
fname2 = sys.argv[2]

file1 = h5py.File(fname,'r')
file2 = h5py.File(fname2, 'w')

for name,value in file1.attrs.items():
    file2.attrs[name]=value 

#file2['coords'] = file1['coords']
dset = file2.create_dataset("coords", data=file1['coords'][:])
dset = file2.create_dataset("time", data=file1['time'][:])
dset = file2.create_dataset("msc_apr", data=file1['msc_apr'][:])
dset = file2.create_dataset("bitmap_ocean", data=file1['bitmap_ocean'][:])
dset = file2.create_dataset("tmsc_part", data=file1['tmsc_part'][:])
dset = file2.create_dataset("quat", data=file1['quat'][:])

d1 = file1['partials'][:]
print(d1.shape)
dset = file2.create_dataset("partials", data=d1.reshape(d1.shape[0],d1.shape[1]//6, 6))


d1 = file1['msc_part'][:]
print(d1.shape)
dset = file2.create_dataset("msc_part", data=d1.reshape(d1.shape[0],d1.shape[1]//6, 6))


file1.close()
file2.close()

