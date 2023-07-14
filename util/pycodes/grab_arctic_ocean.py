#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 13:30:30 2020

@author: julia pfeffer
"""
import numpy as np
import grace_utils as gr
import argparse
import datetime

###############################################################################
## READ ARGUMENTS    
############################################################################### 
# Define arguments   
parser = argparse.ArgumentParser(description=' This programs aims to create an Arctic Ocean region (one primary mascon) including all ocean mascons above 58 degree north')
parser.add_argument("in_masconfile",default="None",type=str,help="complete path to your input mascon file")
parser.add_argument("--out_masconfile",default="None",type=str,help="complete path to your output mascon file")
# Read arguments
args = parser.parse_args()
# Declare arguments
inputfile = args.in_masconfile
outputfile = args.out_masconfile
if outputfile=="None":
     outputfile=inputfile+'_grab_ao'
     
print('')
print('############################################################################################')
print('######################### GRAB ARCTIC OCEAN:',str(datetime.datetime.now())[0:19],'#########################')
print('##########################################################################################')      

## read mascon file
print('')
print('Read input files:',str(datetime.datetime.now())[0:19])    
headers,Pdata,_,Tdata=gr.read_mascon_file(inputfile)
nbt=len(Tdata)

## select ocean ternaries below 60S as Southern Ocean and land ternaries below 60S as Antarctica
print('')
print('Select seawater mascons above 58 degree North:',str(datetime.datetime.now())[0:19])
inTdata=Tdata[(Tdata[:,2]>=58)&(Tdata[:,8]==1029)] # ternaries inside polygons
outTdata=np.vstack((Tdata[(Tdata[:,2]<58)],Tdata[(Tdata[:,2]>=58)&(Tdata[:,8]==1000)]))

## create new primary and ternary arrays
print('')
print('Create new primary and ternary arrays. This can take a few minutes. Relax and finish your tea :-)',str(datetime.datetime.now())[0:19])

##### Count the primaries 
ip_out=np.unique(outTdata[:,9])#  list of unique primaries outside southern ocean
nbp_out=len(ip_out) # number of primaries outside southern ocean
new_nbp=nbp_out+1 # total number of primaries

##### Build empty primary and ternary arrays
new_Pdata=np.array(np.zeros((new_nbp,14)),dtype=object)
new_Tdata=np.array(np.zeros((nbt,14)),dtype=object)

### Fill arrays for primaries and ternaries above 60 degree South and in Antarctica
tcount=-1
for pcpt in np.arange(nbp_out):# primaries
    ternaries=outTdata[outTdata[:,9]==ip_out[pcpt]] # select ternaries outside polygons
    new_tip=np.shape(ternaries)[0] # number of ternary inside one primary 
    pcland=100*len(ternaries[:,8][ternaries[:,8]==1000])/new_tip # percentage of land in primary 
    new_Pdata[pcpt,0]=pcpt+1 # renumber the primary in ascending order
    new_Pdata[pcpt,1]=Pdata[:,1][Pdata[:,0]==ip_out[pcpt]][0] # type of primary unchanged
    new_Pdata[pcpt,2]=1
    new_Pdata[pcpt,3]=new_tip # number of ternary in primary
    new_Pdata[pcpt,4]=Pdata[:,4][Pdata[:,0]==ip_out[pcpt]][0]# latitude unchanged
    new_Pdata[pcpt,5]=Pdata[:,5][Pdata[:,0]==ip_out[pcpt]][0]# longitude unchanged
    new_Pdata[pcpt,6]=Pdata[:,6][Pdata[:,0]==ip_out[pcpt]][0]# radius unchanged
    new_Pdata[pcpt,7]=Pdata[:,7][Pdata[:,0]==ip_out[pcpt]][0] # area unchanged
    new_Pdata[pcpt,8]=Pdata[:,8][Pdata[:,0]==ip_out[pcpt]][0]# altitude unchanged 
    new_Pdata[pcpt,9]=Pdata[:,9][Pdata[:,0]==ip_out[pcpt]][0]# geoid hei unchangedght
    new_Pdata[pcpt,10]=Pdata[:,10][Pdata[:,0]==ip_out[pcpt]][0] # density unchanged
    new_Pdata[pcpt,11]=pcland
    new_Pdata[pcpt,12]=Pdata[:,12][Pdata[:,0]==ip_out[pcpt]][0] # tidal mask unchanged
    new_Pdata[pcpt,13]=Pdata[:,13][Pdata[:,0]==ip_out[pcpt]][0] # region unchanged
    for tcpt in np.arange(new_tip):#ternaries
        tcount=tcount+1
        new_Tdata[tcount,0:9]=ternaries[tcpt,0:9] # ternary number, type, lat, lon, radius, area, alt., geoid height, density unchanged
        new_Tdata[tcount,9]=pcpt+1 # renumber associated primary
        new_Tdata[tcount,10]=pcpt+1 # renumber associated secondary
        new_Tdata[tcount,11:14]=ternaries[tcpt,11:14] # pcolour, scolour and region unchanged
    del ternaries, new_tip, pcland
del outTdata, ip_out   

### Fill arrays for primaries and ternaries in Southern Ocean
ternaries=inTdata# select ternaries for each primary inside polygon
new_tip=np.shape(ternaries)[0] # number of ternaries inside primary
pcland=100*len(ternaries[:,8][ternaries[:,8]==1000])/new_tip # percentage of land in primary
new_Pdata[pcpt+1,0]=nbp_out+1 # renumber primary
new_Pdata[pcpt+1,1]='PDeep'
new_Pdata[pcpt+1,2]=1 #number of secondary in primary
new_Pdata[pcpt+1,3]=new_tip #number of ternary in primary
new_Pdata[pcpt+1,4]=np.mean(ternaries[:,2])# new latitude
new_Pdata[pcpt+1,5]=np.mean(ternaries[:,3])# new longitude
new_Pdata[pcpt+1,6]=np.mean(ternaries[:,4])# new radius
new_Pdata[pcpt+1,7]=np.sum(ternaries[:,5]) # new area
new_Pdata[pcpt+1,8]=np.mean(ternaries[:,6])# new altitude 
new_Pdata[pcpt+1,9]=np.mean(ternaries[:,7])# new geoid height
new_Pdata[pcpt+1,10]=1029#  density defined as a selection criteria
new_Pdata[pcpt+1,11]=pcland
new_Pdata[pcpt+1,12]=0 # tidal mask (default = 0)
new_Pdata[pcpt+1,13]='Arctic'# new region name defined in polygon file
for tcpt in np.arange(new_tip):
    tcount=tcount+1
    new_Tdata[tcount,0:9]=ternaries[tcpt,0:9] # number of ternary unchanged
    new_Tdata[tcount,9]=nbp_out+1  # renumber associated primary
    new_Tdata[tcount,10]=nbp_out+1  # renumber associated secondary
    new_Tdata[tcount,11:13]=ternaries[tcpt,11:13] # pcolour and scolour unchanged
    new_Tdata[tcount,13]='Arctic' # new region name defined in polygon file
del ternaries, new_tip, pcland    

## Create new headers
print('')
print('Create headers:',str(datetime.datetime.now())[0:19])
npex=int(new_nbp)
ntex=int(nbt)
maxtip=int(np.max(new_Pdata[:,3]))
lastcomment='# grab_arctic_ocean.py : define arctic ocean as all ocean mascons above 58 degree north.  Timetag: %s \n'%(str(datetime.datetime.now())[0:19])
new_headers=gr.define_new_mascon_headers(npex,ntex,maxtip,headers,lastcomment)

## Write mascon file
print('Write mascon in file %s:'%(outputfile),str(datetime.datetime.now())[0:19])
gr.write_mascon_file(outputfile,new_headers,new_Pdata,new_Tdata)
del new_Pdata,new_Tdata,new_headers 