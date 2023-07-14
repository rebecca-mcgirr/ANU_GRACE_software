#!/usr/bin/env python3
import numpy as np
import sys
sys.path.append('~/gt/util/pycodes')
import grace_utils as gr
import os
#import argparse
################################################################################
### READ ARGUMENTS    
################################################################################ 
## Define arguments   
#parser = argparse.ArgumentParser(description=' This programs aim to grab ternary mascons inside polygons.')
#parser.add_argument("in_masconfile",default="/scratch/compute1/julia/solutions/masconfile/mascons_stage4",type=str,help="complete path to your mascon file")
#parser.add_argument("area",default=90000000000,type=int,help="resolution of primary mascon in m2")
#parser.add_argument("outputfile",default="/scratch/compute1/julia/solutions/masconfile/mascons_stage4_islands",type=str,help="complete path to your output mascon file")
#parser.add_argument("overflow",default=10,type=int,help="maximum diiference between mascons")
#parser.add_argument("maxsize",default=1000000000,type=int,help="maximum size of the MCF problem")
#parser.add_argument("precision",default=7500,type=int,help="maximum distance displacement")
## Read arguments
#args = parser.parse_args()
# Declare arguments
infile = '/scratch/compute1/julia/solutions/newcontinents/mascons_stage4_JP16'
area = 90000000000
outfile = '/scratch/compute1/julia/solutions/newcontinents/mascons_stage4_JP17'
maxsize = 1000000000000
overflow = 20
precision = 7500
### Read polygon file
_,Pdata,_,_=gr.read_mascon_file(infile)
tempfile=infile+'_temp'
os.system("cp %s %s"%(infile,tempfile))
regions,coords,nbv=gr.read_polygons_ascii('/scratch/compute1/julia/solutions/newcontinents/islands_to_reshape.polygons','all')
nbr=len(regions)
for cpt in np.arange(nbr):
    os.system("~/gt/util/pycodes/pyshape_constrained_size_mascons_RM.py %s %d --overflow=%d --precision=%d --maxsize=%d --region=%s --out_masconfile=%s"%(tempfile,area,overflow,precision,maxsize,regions[cpt],tempfile))

os.system("cp %s %s"%(tempfile,outfile))    
