#!/usr/bin/env python3
import numpy as np
import sys
sys.path.append('~/gt/util/pycodes')
import grace_utils as gr
import os
import argparse
###############################################################################
## READ ARGUMENTS    
############################################################################### 
# Define arguments   
parser = argparse.ArgumentParser(description=' This programs aim to grab ternary mascons inside polygons.')
parser.add_argument("in_masconfile",default="/scratch/compute1/julia/solutions/masconfile/mascons_stage4",type=str,help="complete path to your mascon file")
parser.add_argument("polygonfile",default='/scratch/compute1/julia/solutions/masconfile/islands.polygons',type=str,help="path to your polygon file")
parser.add_argument("outputfile",default="/scratch/compute1/julia/solutions/masconfile/mascons_stage4_islands",type=str,help="complete path to your output mascon file")
parser.add_argument("density",default=1000,type=int,help="default = 1000 for land (oceans = 1029)")
# Read arguments
args = parser.parse_args()
# Declare arguments
infile = args.in_masconfile
polfile = args.polygonfile
outfile = args.outputfile
mydensity = args.density
### Read polygon file
regions,coords,nbv=gr.read_polygons_ascii(polfile,'all')
nbpol=len(coords)
## Run pygrab mascons sequentially
tempfile=infile+'_temp'
os.system("cp %s %s"%(infile,tempfile))
for cpt in np.arange(nbpol):
    os.system("~/gt/util/pycodes/pygrab_mascons.py %s %s --region=%s --density=%d --out_masconfile=%s"%(tempfile,polfile,regions[cpt],mydensity,tempfile))
os.system("cp %s %s"%(tempfile,outfile))    
