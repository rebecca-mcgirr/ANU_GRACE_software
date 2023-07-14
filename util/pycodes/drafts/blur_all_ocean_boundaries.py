#!/usr/bin/env python3
import sys
sys.path.append('~/gt/util/pycodes')
import os
# arguments
infile='/scratch/compute1/julia/solutions/newcontinents/mascons_stage4_JP19'
tempfile='/scratch/compute1/julia/solutions/newcontinents/mascons_stage4_JP19_temp'
finalfile='/scratch/compute1/julia/solutions/newcontinents/mascons_stage4_JP20'
###################
###WIndian / Eindian
os.system("~/gt/util/blur_ocean_boundaries %s %s 6. -55.0 77.5 8.0 77.5 1"%(infile,tempfile))
os.system("~/gt/util/pycodes/pyshape_constrained_size_mascons_RM.py %s 90000000000 --overflow=10 --precision=7500 --maxsize=10000000000 --out_masconfile=%s --region=OceanBlur "%(tempfile,tempfile))
os.system("~/gt/util/pycodes/replace_region_name.py %s %s OceanBlur Indian"%(tempfile,tempfile))  
# rename 
os.system("~/gt/util/pycodes/replace_region_name.py %s %s EIndian Indian"%(tempfile,tempfile))  
os.system("~/gt/util/pycodes/replace_region_name.py %s %s WIndian Indian"%(tempfile,tempfile))  
# Indian/ Atlantic
os.system("~/gt/util/blur_ocean_boundaries %s %s 6. -35 20 -55 20 1"%(tempfile,tempfile))
os.system("~/gt/util/pycodes/pyshape_constrained_size_mascons_RM.py %s 90000000000 --overflow=10 --precision=7500 --maxsize=10000000000 --out_masconfile=%s --region=OceanBlur "%(tempfile,tempfile))
os.system("~/gt/util/pycodes/replace_region_name.py %s %s OceanBlur Indian"%(tempfile,tempfile))  
# Atlantic North/south
#1st half
os.system("~/gt/util/blur_ocean_boundaries %s %s 6. 0 0 0 10 1"%(tempfile,tempfile))
os.system("~/gt/util/pycodes/pyshape_constrained_size_mascons_RM.py %s 90000000000 --overflow=10 --precision=7500 --maxsize=10000000000 --out_masconfile=%s --region=OceanBlur "%(tempfile,tempfile))
os.system("~/gt/util/pycodes/replace_region_name.py %s %s OceanBlur Atlantic"%(tempfile,tempfile))  
#2nd half
os.system("~/gt/util/blur_ocean_boundaries %s %s 6. 0 310 0 360 1"%(tempfile,tempfile))
os.system("~/gt/util/pycodes/pyshape_constrained_size_mascons_RM.py %s 90000000000 --overflow=10 --precision=7500 --maxsize=10000000000 --out_masconfile=%s --region=OceanBlur "%(tempfile,tempfile))
os.system("~/gt/util/pycodes/replace_region_name.py %s %s OceanBlur Atlantic"%(tempfile,tempfile))  
#rename
os.system("~/gt/util/pycodes/replace_region_name.py %s %s NAtlantic Atlantic"%(tempfile,tempfile))  
os.system("~/gt/util/pycodes/replace_region_name.py %s %s SAtlantic Atlantic"%(tempfile,tempfile))  
# Pacific 
#1st quarter
os.system("~/gt/util/blur_ocean_boundaries %s %s 6. 0 200 0 280 1"%(tempfile,tempfile))
os.system("~/gt/util/pycodes/pyshape_constrained_size_mascons_RM.py %s 90000000000 --overflow=10 --precision=7500 --maxsize=10000000000 --out_masconfile=%s --region=OceanBlur "%(tempfile,tempfile))
os.system("~/gt/util/pycodes/replace_region_name.py %s %s OceanBlur Pacific"%(tempfile,tempfile))  
#2nd quarter
os.system("~/gt/util/blur_ocean_boundaries %s %s 6. 55 200 0 200 1"%(tempfile,tempfile))
os.system("~/gt/util/pycodes/pyshape_constrained_size_mascons_RM.py %s 90000000000 --overflow=10 --precision=7500 --maxsize=10000000000 --out_masconfile=%s --region=OceanBlur "%(tempfile,tempfile))
os.system("~/gt/util/pycodes/replace_region_name.py %s %s OceanBlur Pacific"%(tempfile,tempfile))  
#3rd quarter
os.system("~/gt/util/blur_ocean_boundaries %s %s 6. 0 135 0 200 1"%(tempfile,tempfile))
os.system("~/gt/util/pycodes/pyshape_constrained_size_mascons_RM.py %s 90000000000 --overflow=10 --precision=7500 --maxsize=10000000000 --out_masconfile=%s --region=OceanBlur "%(tempfile,tempfile))
os.system("~/gt/util/pycodes/replace_region_name.py %s %s OceanBlur Pacific"%(tempfile,tempfile))  
#4th quarter
os.system("~/gt/util/blur_ocean_boundaries %s %s 6. 0 200 -55 200 1"%(tempfile,tempfile))
os.system("~/gt/util/pycodes/pyshape_constrained_size_mascons_RM.py %s 90000000000 --overflow=10 --precision=7500 --maxsize=10000000000 --out_masconfile=%s --region=OceanBlur "%(tempfile,tempfile))
os.system("~/gt/util/pycodes/replace_region_name.py %s %s OceanBlur Pacific"%(tempfile,tempfile))  
#rename
os.system("~/gt/util/pycodes/replace_region_name.py %s %s NEPacific Pacific"%(tempfile,tempfile))  
os.system("~/gt/util/pycodes/replace_region_name.py %s %s NWPacific Pacific"%(tempfile,tempfile)) 
os.system("~/gt/util/pycodes/replace_region_name.py %s %s SEPacific Pacific"%(tempfile,tempfile))  
os.system("~/gt/util/pycodes/replace_region_name.py %s %s SWPacific Pacific"%(tempfile,tempfile))  
# Southern Ocean
#1st quarter
os.system("~/gt/util/blur_ocean_boundaries %s %s 6. -55 300 -55 360 1"%(tempfile,tempfile))
os.system("~/gt/util/pycodes/pyshape_constrained_size_mascons_RM.py %s 90000000000 --overflow=10 --precision=7500 --maxsize=10000000000 --out_masconfile=%s --region=OceanBlur "%(tempfile,tempfile))
os.system("~/gt/util/pycodes/replace_region_name.py %s %s OceanBlur Southern"%(tempfile,tempfile))  
#2nd quarter
os.system("~/gt/util/blur_ocean_boundaries %s %s 6. -55 180 -55 280 1"%(tempfile,tempfile))
os.system("~/gt/util/pycodes/pyshape_constrained_size_mascons_RM.py %s 90000000000 --overflow=10 --precision=7500 --maxsize=10000000000 --out_masconfile=%s --region=OceanBlur "%(tempfile,tempfile))
os.system("~/gt/util/pycodes/replace_region_name.py %s %s OceanBlur Southern"%(tempfile,tempfile))  
#3rd quarter
os.system("~/gt/util/blur_ocean_boundaries %s %s 6. -55 0 -55 90 1"%(tempfile,tempfile))
os.system("~/gt/util/pycodes/pyshape_constrained_size_mascons_RM.py %s 90000000000 --overflow=10 --precision=7500 --maxsize=10000000000 --out_masconfile=%s --region=OceanBlur "%(tempfile,tempfile))
os.system("~/gt/util/pycodes/replace_region_name.py %s %s OceanBlur Southern"%(tempfile,tempfile))  
#4th quarter
os.system("~/gt/util/blur_ocean_boundaries %s %s 6. -55 90 -55 180 1"%(tempfile,tempfile))
os.system("~/gt/util/pycodes/pyshape_constrained_size_mascons_RM.py %s 90000000000 --overflow=10 --precision=7500 --maxsize=10000000000 --out_masconfile=%s --region=OceanBlur "%(tempfile,tempfile))
os.system("~/gt/util/pycodes/replace_region_name.py %s %s OceanBlur Southern"%(tempfile,tempfile))  
# Arctic
os.system("~/gt/util/blur_ocean_boundaries %s %s 6. 62 300 62 360 1"%(tempfile,tempfile))
os.system("~/gt/util/pycodes/pyshape_constrained_size_mascons_RM.py %s 90000000000 --overflow=10 --precision=7500 --maxsize=10000000000 --out_masconfile=%s --region=OceanBlur "%(tempfile,tempfile))
os.system("~/gt/util/pycodes/replace_region_name.py %s %s OceanBlur ArcticO"%(tempfile,tempfile))  
#
os.system("cp %s %s"%(tempfile,finalfile))