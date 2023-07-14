#!/usr/bin/env python3
import sys
sys.path.append('~/gt/util/pycodes')
import os
# arguments
infile='/scratch/compute1/julia/solutions/newcontinents/mascons_stage4_JP20'
tempfile='/scratch/compute1/julia/solutions/newcontinents/mascons_stage4_JP20_temp'
finalfile='/scratch/compute1/julia/solutions/newcontinents/mascons_stage4_JP21'
###################
# Indian/ Atlantic
os.system("~/gt/util/blur_ocean_boundaries %s %s 6. -55 20 -35 20 1"%(infile,tempfile))
os.system("~/gt/util/pycodes/pyshape_constrained_size_mascons_RM.py %s 90000000000 --overflow=10 --precision=7500 --maxsize=10000000000 --out_masconfile=%s --region=OceanBlur "%(tempfile,tempfile))
os.system("~/gt/util/pycodes/replace_region_name.py %s %s OceanBlur Indian"%(tempfile,tempfile))  
# Pacific 
#2nd quarter
os.system("~/gt/util/blur_ocean_boundaries %s %s 6. 0 200 55 200 1"%(tempfile,tempfile))
os.system("~/gt/util/pycodes/pyshape_constrained_size_mascons_RM.py %s 90000000000 --overflow=10 --precision=7500 --maxsize=10000000000 --out_masconfile=%s --region=OceanBlur "%(tempfile,tempfile))
os.system("~/gt/util/pycodes/replace_region_name.py %s %s OceanBlur Pacific"%(tempfile,tempfile))  
#4th quarter
os.system("~/gt/util/blur_ocean_boundaries %s %s 6. -55 200 0 200 1"%(tempfile,tempfile))
os.system("~/gt/util/pycodes/pyshape_constrained_size_mascons_RM.py %s 90000000000 --overflow=10 --precision=7500 --maxsize=10000000000 --out_masconfile=%s --region=OceanBlur "%(tempfile,tempfile))
os.system("~/gt/util/pycodes/replace_region_name.py %s %s OceanBlur Pacific"%(tempfile,tempfile))  
#rename
