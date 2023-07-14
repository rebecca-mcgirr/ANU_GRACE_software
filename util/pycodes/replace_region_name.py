#!/usr/bin/env python3
import numpy as np
import grace_utils as gr
import argparse
###############################################################################
## READ ARGUMENTS    
############################################################################### 
# Define arguments   
parser = argparse.ArgumentParser(description=' This programs aim to rename a region in primary, secondary and ternary mascons.')
parser.add_argument("in_masconfile",default="/scratch/compute1/julia/solutions/masconfile/mascons_stage4",type=str,help="complete path to your input mascon file")
parser.add_argument("out_masconfile",default="/scratch/compute1/julia/solutions/masconfile/mascons_stage4",type=str,help="complete path to your output mascon file")
parser.add_argument("in_region",default="myregion1",type=str,help="name of the region to change")
parser.add_argument("out_region",default="myregion2",type=str,help="name of the region in the output mascon file")
# Read arguments
args = parser.parse_args()
# Declare arguments
inputfile = args.in_masconfile
outputfile = args.out_masconfile
inregion= args.in_region
outregion = args.out_region
###############################################################################
## READ INPUT FILES    
###############################################################################
#mascon file
print('')
print('Replace region %s by %s in %s'%(inregion,outregion,inputfile))
headers,Pdata,Sdata,Tdata=gr.read_mascon_file(inputfile)
nbp=len(Pdata[:,0])
nbtern=len(Tdata[:,0])
tcount=0
for cpt in np.arange(nbp):
    if Pdata[cpt,-1]==inregion:
        myregion=outregion
    else:
        myregion=Pdata[cpt,-1]    
    Pdata[cpt,-1]=myregion
    Sdata[cpt,-1]=myregion
    nbt=Pdata[cpt,3]
    for tpt in np.arange(nbt):
        Tdata[tcount,-1]=myregion
        Tdata[tcount,1]='T'+Pdata[cpt,1][1:]
        tcount=tcount+1  
gr.write_mascon_file(outputfile,headers,Pdata,Tdata)
print('New region name written in %s'%(outputfile))
print('')
