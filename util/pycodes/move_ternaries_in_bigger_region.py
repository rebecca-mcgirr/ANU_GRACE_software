#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 14:52:42 2020

@author: julia
"""


import os
import argparse
import grace_utils as gr
import numpy as np
import datetime as datetime

def read_arguments():
    parser = argparse.ArgumentParser(description=' This programs aims to merge all primary mascons in a given region')
    parser.add_argument("in_masconfile",default="None",type=str,help="complete path to your input mascon file")
    parser.add_argument("--out_masconfile",default="None",type=str,help="complete path to your output mascon file")
    parser.add_argument("--minsize",default=25,type=int,help="size under which a region is considered too small")
    parser.add_argument("--maxdist",default=1000,type=int,help="maximum distance between regions to be merged (km)")    
    parser.add_argument("--ftype",default=0,type=int,help="set to 1 to allow to merge regions of different type")
  
    # Read arguments
    args = parser.parse_args()
    # Declare arguments
    infile = args.in_masconfile
    outfile = args.out_masconfile
    if outfile=="None":
        outfile=infile+'_merged'
    minsize = args.minsize
    maxdist = args.maxdist
    ftype = args.ftype
    return infile,outfile,minsize,maxdist,ftype

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

def main():
    #Print program identification
    print('~/gt/util/pycodes/move_ternaries_in_bigger_region.py started at %s'%(str(datetime.datetime.now())[0:19]))
    
    # read arguments
    infile,outfile,minsize,maxdist,ftype=read_arguments()
    
    #create temporary file
    tempfile=outfile+'_temp'
    
    # Create temporary mascon file to write in
    tempfile=outfile+'_temp'
    os.system("cp %s %s"%(infile,tempfile))
    
    # Read input masconfile
    print('')
    print('Read mascon file: %s'%(str(datetime.datetime.now())[0:19]))   
    headers,Pdata,_,Tdata=gr.read_mascon_file(infile)
    nbt=len(Tdata)
    
    #Get number of ternary per regions
    regions=np.unique(Tdata[:,13])  
    nbr=len(regions)
    ntir=np.zeros(nbr)
    Tsmall=[]   
    Tlarge=[]
    for rct in np.arange(nbr):
        ternaries=Tdata[Tdata[:,13]==regions[rct],:]
        ntir[rct]=len(ternaries[:,0])
        if ntir[rct]<minsize:
            Tsmall.append(ternaries)
        else:
            Tlarge.append(ternaries)
    print('Found %d ternaries in small regions. Timetag: %s.'%(len(Tsmall),str(datetime.datetime.now())[0:19]))
    if len(Tsmall)>0:            
        Tsmall=np.vstack(Tsmall) 
        Tlarge=np.vstack(Tlarge) 
        
        # Find primary where to move ternaries of small regions       
        tpd=np.zeros((len(Tsmall),3))
        for tpt in np.arange(len(Tsmall)):            
            tpool=Tdata[(Tdata[:,13]!=Tsmall[tpt,13])&(Tdata[:,8]==Tsmall[tpt,8]),:] # ternary can be merged to a different primary with the same density
            dist=compute_distance(Tsmall[tpt,2],Tsmall[tpt,3],tpool[:,2],tpool[:,3])
            index=np.where(dist==np.nanmin(dist))[0][0]
            if np.min(dist)<maxdist:
                tpd[tpt,0]=int(Tsmall[tpt,0]) # small ternary
                tpd[tpt,1]=int(tpool[index,9]) # closest primary
                tpd[tpt,2]=np.nanmin(dist) # distance 
                del tpool,dist,index
            else:
                if ftype==1:
                    tpool=Tdata[(Tdata[:,13]!=Tsmall[tpt,13])&(Tdata[:,8]!=Tsmall[tpt,8]),:] # ternary can be merged to a different primary with a different density
                    dist=compute_distance(Tsmall[tpt,2],Tsmall[tpt,3],tpool[:,2],tpool[:,3]) 
                    index=np.where(dist==np.nanmin(dist))[0][0]
                    tpd[tpt,0]=int(Tsmall[tpt,0]) # small ternary
                    tpd[tpt,1]=int(tpool[index,9]) # closest primary
                    tpd[tpt,2]=np.nanmin(dist) # distance 
                    del tpool,dist,index
                else:
                    tpd[tpt,0]=int(Tsmall[tpt,0])
                    tpd[tpt,1]=int(Tsmall[tpt,9]) # ternary not moved, same primary
                    tpd[tpt,2]=0
            newregion=Pdata[Pdata[:,0]==tpd[tpt,1],13]        
            print('ternary %d in region %s moved to primary %d in %s. Distance= %8.2f km. Timetag: %s'%(tpd[tpt,0],Tsmall[tpt,13],tpd[tpt,1],newregion,tpd[tpt,2],str(datetime.datetime.now())[0:19]))  
  
            
        print('')
        print('Create new primary and ternary arrays. This can take a few minutes. Relax and finish your tea :-)',str(datetime.datetime.now())[0:19])        
        
        ## associate new primary to solo land mascons
        all_Tdata=np.vstack((Tlarge,Tsmall))
        all_Tdata[len(Tlarge[:,0]):,9]=tpd[:,1]
        all_Tdata[len(Tlarge[:,0]):,10]=tpd[:,1]
        
        ##### Get list and number of unique primary mascons
        ip=np.unique(all_Tdata[:,9])
        nbp=len(ip)
        
        ##### Build empty primary and ternary arrays
        new_Pdata=np.array(np.zeros((nbp,14)),dtype=object)
        new_Tdata=np.array(np.zeros((nbt,14)),dtype=object)   
        
        ### Fill arrays for primaries and ternaries above 60 degree South and in Antarctica
        tcount=-1
        for pcpt in np.arange(nbp):# primaries
            ternaries=all_Tdata[all_Tdata[:,9]==ip[pcpt]] # select ternaries outside polygons
            new_tip=np.shape(ternaries)[0] # number of ternary inside one primary 
            pcland=100*len(ternaries[:,8][ternaries[:,8]==1000])/new_tip # percentage of land in primary 
            new_Pdata[pcpt,0]=pcpt+1 # renumber the primary in ascending order
            new_Pdata[pcpt,1]=Pdata[:,1][Pdata[:,0]==ip[pcpt]][0] # type of primary unchanged
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
            new_Pdata[pcpt,12]=Pdata[:,12][Pdata[:,0]==ip[pcpt]][0] # tidal mask unchanged
            new_Pdata[pcpt,13]=Pdata[:,13][Pdata[:,0]==ip[pcpt]][0] # region unchanged
            for tcpt in np.arange(new_tip):#ternaries
                tcount=tcount+1
                new_Tdata[tcount,0:9]=ternaries[tcpt,0:9] # ternary number, type, lat, lon, radius, area, alt., geoid height, density unchanged
                new_Tdata[tcount,9]=pcpt+1 # renumber associated primary
                new_Tdata[tcount,10]=pcpt+1 # renumber associated secondary
                new_Tdata[tcount,11:13]=ternaries[tcpt,11:13] # pcolour, scolour and region unchanged
                new_Tdata[tcount,13]=Pdata[:,13][Pdata[:,0]==ip[pcpt]][0]
            del ternaries, new_tip, pcland
        
        
        ## Create new headers
        print('')
        print('Create headers:',str(datetime.datetime.now())[0:19])
        npex=int(nbp)
        ntex=int(nbt)
        maxtip=int(np.max(new_Pdata[:,3]))
        if ftype==1:
            answer='yes'
        else:
            answer='no'
        lastcomment='# merge_small_regions.py : move ternaries of small regions (<%d). Maximum distance allowed between ternary and new primary = %d km. Merged with primary of different type: %s. Timetag: %s \n'%(minsize,maxdist,answer,str(datetime.datetime.now())[0:19])
        new_headers=gr.define_new_mascon_headers(npex,ntex,maxtip,headers,lastcomment)
        
        ## Write mascon file
        print('Write mascon in file %s:'%(outfile),str(datetime.datetime.now())[0:19])
        gr.write_mascon_file(outfile,new_headers,new_Pdata,new_Tdata)
        del new_Pdata,new_Tdata,new_headers       
    
 
if __name__ == '__main__':
    main()    
