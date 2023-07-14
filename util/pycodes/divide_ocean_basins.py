#!/usr/bin/env python3
import numpy as np
import sys
sys.path.append('~/gt/util/pycodes')
import grace_utils as gr
import argparse
###############################################################################
## READ ARGUMENTS    
############################################################################### 
## Define arguments   
parser = argparse.ArgumentParser(description=' This programs aim to divide ocean basins in subbasins of smaller area.')
parser.add_argument("in_masconfile",default="/scratch/compute1/julia/solutions/masconfiles/mascons_stage5_V004",type=str,help="complete path to your mascon file")
parser.add_argument("outputfile",default="/scratch/compute1/julia/solutions/masconfile/mascons_stage4_islands",type=str,help="complete path to your output mascon file")
## Read arguments
args = parser.parse_args()
## Declare arguments
infile = args.in_masconfile
outfile = args.outputfile
### Read mascon file
headers,Pdata,Sdata,Tdata=gr.read_mascon_file(infile)
nbp=len(Pdata[:,0])
nbtern=len(Tdata[:,0])
### Read all primaries an rename regions where appropriate
tcount=0
for cpt in np.arange(nbp):
    if Pdata[cpt,-1]=='Pacific':
        if (Pdata[cpt,4]>0)&(Pdata[cpt,5]>200):
            myregion='NEPacific'
        elif (Pdata[cpt,4]>0)&(Pdata[cpt,5]<=200):
            myregion='NWPacific'
        elif (Pdata[cpt,4]<=0)&(Pdata[cpt,5]>200):
            myregion='SEPacific'    
        else:
            myregion='SWPacific' 
    elif Pdata[cpt,-1]=='Atlantic':
        if Pdata[cpt,4]>0:
            myregion='NAtlantic'  
        else:
            myregion='SAtlantic'  
    elif Pdata[cpt,-1]=='Indian':
        if Pdata[cpt,5]>77.5:
            myregion='EIndian'  
        else:
            myregion='WIndian'                    
    else:
        myregion=Pdata[cpt,-1]    
    Pdata[cpt,-1]=myregion
    Sdata[cpt,-1]=myregion
    nbt=Pdata[cpt,3]
    #Read each ternary in this primary and rename region where appropriate
    for tpt in np.arange(nbt):
        Tdata[tcount,-1]=myregion
        Tdata[tcount,1]='T'+Pdata[cpt,1][1:]
        tcount=tcount+1  
gr.write_mascon_file(outfile,headers,Pdata,Tdata)
