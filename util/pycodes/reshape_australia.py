#!/usr/bin/env python3
import numpy as np
import grace_utils as gr
import datetime as dt
import os
import argparse
################################################################################
### READ ARGUMENTS    
################################################################################ 
## Define arguments   
parser = argparse.ArgumentParser(description=' This programs aim to reshape australian basins at requited resolution.')
parser.add_argument("input_file",default="None",type=str,help="complete path to your input mascon file (netcdf) ")
#parser.add_argument("output_file",default="None",type=str,help="complete path to your output mascon file (netcdf)")
parser.add_argument("--resolution",default=40000000000,type=float,help="resolution of the mascon file (m2)")
## Read arguments
args = parser.parse_args()
## Declare arguments
input_file = args.input_file
#output_file=args.output_file
parea=args.resolution
###############################################################################
## READ INPUT FILES    
###############################################################################
#mascon file
print('Read mascon file: %s' %(str(dt.datetime.now())))
headers,Pdata,Sdata,Tdata=gr.read_mascon_file(input_file)
nbp=len(Pdata[:,0])
nbtern=len(Tdata[:,0])
#
plist=[]
regions=[]
for cpt in np.arange(len(Pdata)):
    if Pdata[cpt,-1][0:3]=='AUS':
        #plist.append(cpt)
        regions.append(Pdata[cpt,-1])
for myregion in regions: 
    headers,Pdata,Sdata,Tdata=gr.read_mascon_file(input_file)
    pnumber=Pdata[Pdata[:,-1]==myregion,0]
    nt_ideal=parea/np.mean(Tdata[Tdata[:,9]==pnumber,5])   
    threshold=1.75*nt_ideal
    if Pdata[Pdata[:,-1]==myregion,3]>threshold:
        os.system('~/gt/util/pycodes/pyshape_constrained_size_mascons_JP.py %s %f --region=%s --out_masconfile=%s'%(input_file,parea,myregion,input_file))
###############################################################################
