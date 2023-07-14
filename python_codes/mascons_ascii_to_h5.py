import h5py
from gracetools.io.mascons import mascons
import argparse


###############################################################################
## READ ARGUMENTS    
############################################################################### 
# Define arguments   
parser = argparse.ArgumentParser(description=' This program will convert an ascii mascons file into a hdf5 version.')
parser.add_argument("--mascons_ascii",default="mascons_stage5_V006_200km",type=str,help="input ascii mascon file")
parser.add_argument("--mascons_hdf5",default='none',type=str,help="output hdf5 file ")

print("mascons_ascii_to_h5 --mascons_ascii mascons_stage5_V006_200km --mascons_hdf5 mascons_stage5_V006_200km.h5")


# Read arguments
args = parser.parse_args()
# Declare arguments
mascons_hdf5 = args.mascons_hdf5
mascons_ascii= args.mascons_ascii


# read the ascii mascons file
print("Reading ascii mascon file:",mascons_ascii)
msc_ascii = mascons(mascons_ascii)
msc_ascii.read_mascons_txt()


# write it out in hdf5 format
print("Creating file",mascons_hdf5)
mascons.to_h5(msc_ascii, mascons_hdf5)


