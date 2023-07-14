#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 13:54:01 2020

@author: julia pfeffer
"""

import numpy as np
import grace_utils as gr
import argparse
import datetime

def compute_distance(plat,plon,latvec,lonvec):
    nb=len(latvec)
    R = 6373.0
    latvec=latvec.astype(float)
    lonvec=lonvec.astype(float)
    lonvec[lonvec>180]==lonvec[lonvec>180]-360*np.ones(len(lonvec[lonvec>180]))
    lat1 = np.radians(plat)*np.ones(nb)
    lon1 = np.radians(plon)*np.ones(nb)
    lat2 = np.radians(latvec)
    lon2 = np.radians(lonvec)
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.square(np.sin(dlat / 2)) + np.cos(lat1) * np.cos(lat2) * np.square(np.sin(dlon / 2))
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))  
    distance = R * c
    return distance

###############################################################################
## READ ARGUMENTS    
############################################################################### 
# Define arguments   
parser = argparse.ArgumentParser(description=' This programs aims to merge ocean mascons left alone with the closest appropriate mascon')
parser.add_argument("in_masconfile",default="None",type=str,help="complete path to your input mascon file")
parser.add_argument("--out_masconfile",default="None",type=str,help="complete path to your output mascon file")
# Read arguments
args = parser.parse_args()
# Declare arguments
inputfile = args.in_masconfile
outputfile = args.out_masconfile
if outputfile=="None":
      outputfile=inputfile+'_ofixed'
 
print('')
print('Merge solo ocean mascons with closest primary:',str(datetime.datetime.now())[0:19])

print('')
print('Read mascon file: %s'%(str(datetime.datetime.now())[0:19]))   
headers,Pdata,_,Tdata=gr.read_mascon_file(inputfile)
nbt=len(Tdata)

print('')
print('Select solo ocean ternaries: %s'%(str(datetime.datetime.now())[0:19])) 
outTdata=Tdata[Tdata[:,13]=='Ocean'] # ternaries inside polygons
nbt_out=len(outTdata[:,0])
print('Found %d ocean ternary mascons to merge with other primaries: %s'%(nbt_out,str(datetime.datetime.now())[0:19])) 

if nbt_out>0:
    print('')
    print('Look for closest land/oceans mascon: %s'%(str(datetime.datetime.now())[0:19])) 
    
    ### select all ternaries with regions other than Land
    inTdata=Tdata[Tdata[:,13]!='Ocean']
    
    ## create a ternary-primary-distance array 
    
    tpd=np.zeros((nbt_out,3))
    for tpt in np.arange(nbt_out):
        tselect=inTdata[inTdata[:,8]==1029]
        dist=compute_distance(outTdata[tpt,2],outTdata[tpt,3],tselect[:,2],tselect[:,3])
        if np.min(dist)>1000:
            del tselect, dist
            tselect=inTdata[inTdata[:,8]==1000]
            dist=compute_distance(outTdata[tpt,2],outTdata[tpt,3],tselect[:,2],tselect[:,3])    
         
        index=np.where(dist==np.nanmin(dist))[0][0]
        tpd[tpt,0]=int(outTdata[tpt,0]) # solo ternary
        tpd[tpt,1]=int(tselect[index,9]) # closest primary
        tpd[tpt,2]=np.nanmin(dist) # distance  
        print('Ocean ternary %d (%6.2f,%6.2f) merged with primary %d in %s. Timetag: %s'%(tpd[tpt,0],outTdata[tpt,2],outTdata[tpt,3],tpd[tpt,1],tselect[index,13],str(datetime.datetime.now())[0:19]))  
        del tselect,dist,index
    
    
    print('')
    print('Create new primary and ternary arrays. This can take a few minutes. Relax and finish your tea :-)',str(datetime.datetime.now())[0:19])        
    
    ## associate new primary to solo land mascons
    all_Tdata=np.vstack((inTdata,outTdata))
    all_Tdata[len(inTdata[:,0]):,9]=tpd[:,1]
    all_Tdata[len(inTdata[:,0]):,10]=tpd[:,1]
    
    ##### Get list and number of unique primary mascons
    ip_in=np.unique(inTdata[:,9])
    nbp_in=len(ip_in)
    
    ##### Build empty primary and ternary arrays
    new_Pdata=np.array(np.zeros((nbp_in,14)),dtype=object)
    new_Tdata=np.array(np.zeros((nbt,14)),dtype=object)   
    
    ### Fill arrays for primaries and ternaries above 60 degree South and in Antarctica
    tcount=-1
    for pcpt in np.arange(nbp_in):# primaries
        ternaries=all_Tdata[all_Tdata[:,9]==ip_in[pcpt]] # select ternaries outside polygons
        new_tip=np.shape(ternaries)[0] # number of ternary inside one primary 
        pcland=100*len(ternaries[:,8][ternaries[:,8]==1000])/new_tip # percentage of land in primary 
        new_Pdata[pcpt,0]=pcpt+1 # renumber the primary in ascending order
        new_Pdata[pcpt,1]=Pdata[:,1][Pdata[:,0]==ip_in[pcpt]][0] # type of primary unchanged
        new_Pdata[pcpt,2]=1
        new_Pdata[pcpt,3]=new_tip # number of ternary in primary
        new_Pdata[pcpt,4]=np.mean(ternaries[:,2])# new latitude
        new_Pdata[pcpt,5]=np.mean(ternaries[:,3])# new longitude
        new_Pdata[pcpt,6]=np.mean(ternaries[:,4])# new radius
        new_Pdata[pcpt,7]=np.sum(ternaries[:,5]) # new area
        new_Pdata[pcpt,8]=np.mean(ternaries[:,6])# new altitude 
        new_Pdata[pcpt,9]=np.mean(ternaries[:,7])# new geoid height
        if pcland>50:
            new_Pdata[pcpt,10]=1000
        else: 
            new_Pdata[pcpt,10]=1029
        new_Pdata[pcpt,11]=pcland
        new_Pdata[pcpt,12]=Pdata[:,12][Pdata[:,0]==ip_in[pcpt]][0] # tidal mask unchanged
        new_Pdata[pcpt,13]=Pdata[:,13][Pdata[:,0]==ip_in[pcpt]][0] # region unchanged
        for tcpt in np.arange(new_tip):#ternaries
            tcount=tcount+1
            new_Tdata[tcount,0:9]=ternaries[tcpt,0:9] # ternary number, type, lat, lon, radius, area, alt., geoid height, density unchanged
            new_Tdata[tcount,9]=pcpt+1 # renumber associated primary
            new_Tdata[tcount,10]=pcpt+1 # renumber associated secondary
            new_Tdata[tcount,11:13]=ternaries[tcpt,11:13] # pcolour, scolour and region unchanged
            new_Tdata[tcount,13]=Pdata[:,13][Pdata[:,0]==ip_in[pcpt]][0]
        del ternaries, new_tip, pcland
    
    
    ## Create new headers
    print('')
    print('Create headers:',str(datetime.datetime.now())[0:19])
    npex=int(nbp_in)
    ntex=int(nbt)
    maxtip=int(np.max(new_Pdata[:,3]))
    lastcomment='# fix_solo_ocean_mascons.py : merge solo ocean mascons with appropriate primary. Timetag: %s \n'%(str(datetime.datetime.now())[0:19])
    new_headers=gr.define_new_mascon_headers(npex,ntex,maxtip,headers,lastcomment)
    
    ## Write mascon file
    print('Write mascon in file %s:'%(outputfile),str(datetime.datetime.now())[0:19])
    gr.write_mascon_file(outputfile,new_headers,new_Pdata,new_Tdata)
    del new_Pdata,new_Tdata,new_headers  

else:
    print('All ocean ternary mascons in %s are already associated with a region. No output mascon file created. Timetag: %s.'%(inputfile),str(datetime.datetime.now())[0:19])