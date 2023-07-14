#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 12:17:51 2020

@author: julia
"""
import argparse
import grace_utils as gr
import numpy as np
import datetime as datetime

def read_arguments():
    parser = argparse.ArgumentParser(description=' This programs aims to merge all primary mascons in a given region')
    parser.add_argument("in_masconfile",default="None",type=str,help="complete path to your input mascon file")
    parser.add_argument("--outfile",default="None",type=str,help="complete path to your output mascon file")
    parser.add_argument("--minsize",default=25,type=int,help="size under which a region is considered too small")
    parser.add_argument("--maxdist",default=1000,type=int,help="maximum distance between regions to be merged (km)")    
    parser.add_argument("--ftype",default=0,type=int,help="set to 1 to allow to merge regions of different type")
  
    # Read arguments
    args = parser.parse_args()
    # Declare arguments
    infile = args.in_masconfile
    outfile = args.outfile
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
    print('')
    print('~/gt/util/pycodes/merge_small_primaries.py started at %s'%(str(datetime.datetime.now())[0:19]))
    
    # read arguments
    infile,outfile,minsize,maxdist,ftype=read_arguments()   
    
    # Read input masconfile
    print('')
    print('Read mascon file: %s'%(str(datetime.datetime.now())[0:19]))   
    headers,Pdata,_,Tdata=gr.read_mascon_file(infile)
    nbt=len(Tdata)
    
    #Find small primaries
    Psmall=Pdata[Pdata[:,3]<minsize]        
    print('Found %d primaries with less than %d ternaries. Timetag: %s.'%(len(Psmall),minsize,str(datetime.datetime.now())[0:19]))
    if len(Psmall)>0:

        ppd=np.zeros((len(Psmall),3))
        #Find closest primary
        for pct in np.arange(len(Psmall)):
            #print(Psmall[pct,:])
            # Find closest primary of same density
            Ppool=Pdata[(Pdata[:,0]!=Psmall[pct,0])&(Pdata[:,10]==Psmall[pct,10]),:] # ternary can be merged to a different primary with the same density
            dist=compute_distance(Psmall[pct,4],Psmall[pct,5],Ppool[:,4],Ppool[:,5])
            try:
                index=np.where(dist==np.nanmin(dist))[0][0]
            except:
                print('Could not find closest primary for primary number %d in %s'%(Psmall[pct,0],Psmall[pct,13]))
                ppd[pct,0]=int(Psmall[pct,0])
                ppd[pct,1]=int(Psmall[pct,0]) # ternary not moved, same primary
                ppd[pct,2]=0  
            else:    
                if np.min(dist)<maxdist:
                    ppd[pct,0]=int(Psmall[pct,0]) # small ternary
                    ppd[pct,1]=int(Ppool[index,0]) # closest primary
                    ppd[pct,2]=np.nanmin(dist) # distance 
                    del Ppool,dist,index
                else:
                    if ftype==1:
                        Ppool=Pdata[(Pdata[:,0]!=Psmall[pct,0])&(Pdata[:,10]!=Psmall[pct,10]),:] # ternary can be merged to a different primary with the same density
                        dist=compute_distance(Psmall[pct,4],Psmall[pct,5],Ppool[:,4],Ppool[:,5])
                        index=np.where(dist==np.nanmin(dist))[0][0]
                        if np.min(dist)<maxdist:
                            ppd[pct,0]=int(Psmall[pct,0]) # small ternary
                            ppd[pct,1]=int(Ppool[index,0]) # closest primary
                            ppd[pct,2]=np.nanmin(dist) # distance 
                            del Ppool,dist,index
                        else:
                            ppd[pct,0]=int(Psmall[pct,0])
                            ppd[pct,1]=int(Psmall[pct,0]) # ternary not moved, same primary
                            ppd[pct,2]=0                                       
                    else:
                        ppd[pct,0]=int(Psmall[pct,0])
                        ppd[pct,1]=int(Psmall[pct,0]) # ternary not moved, same primary
                        ppd[pct,2]=0
      
        #Apply new primary number to ternaries
        all_Tdata=Tdata       
        for pct in np.arange(len(ppd[:,0])):
            jpt=np.where(ppd[:pct,1]==ppd[pct,0])[0]
            if len(jpt)==1:
                if ppd[jpt[0],1]==ppd[pct,0]:
                    ppd[pct,1]=ppd[pct,0]
                    ppd[pct,2]=0
            old_region=Pdata[Pdata[:,0]==int(ppd[pct,0]),13][0]
            new_region=Pdata[Pdata[:,0]==int(ppd[pct,1]),13][0]          
            print('Primary %d in %s moved to primary %d in %s. Distance= %8.2f km. Timetag: %s'%(ppd[pct,0],old_region,ppd[pct,1],new_region,ppd[pct,2],str(datetime.datetime.now())[0:19]))                 
            all_Tdata[all_Tdata[:,9]==ppd[pct,0],9]=int(ppd[pct,1])*np.ones(len(all_Tdata[all_Tdata[:,9]==ppd[pct,0],9]))
    
        # count new primaries
        new_ip=np.unique(all_Tdata[:,9])
        new_nbp=len(new_ip)
        
        #Create new primary and ternary arrays
        print('')
        print('Create new primary and ternary arrays. This can take a few minutes. Relax and finish your tea :-)',str(datetime.datetime.now())[0:19])        
            
        new_Pdata=np.array(np.zeros((new_nbp,14)),dtype=object)
        new_Tdata=np.array(np.zeros((nbt,14)),dtype=object)   
        
        ### Fill arrays for primaries and ternaries above 60 degree South and in Antarctica
        tcount=-1
        for pcpt in np.arange(new_nbp):# primaries
            ternaries=all_Tdata[all_Tdata[:,9]==new_ip[pcpt]] # select ternaries outside polygons
            new_tip=np.shape(ternaries)[0] # number of ternary inside one primary 
            pcland=100*len(ternaries[:,8][ternaries[:,8]==1000])/new_tip # percentage of land in primary 
            new_Pdata[pcpt,0]=pcpt+1 # renumber the primary in ascending order
            new_Pdata[pcpt,1]=Pdata[:,1][Pdata[:,0]==new_ip[pcpt]][0] # type of primary unchanged
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
            new_Pdata[pcpt,12]=Pdata[:,12][Pdata[:,0]==new_ip[pcpt]][0] # tidal mask unchanged
            new_Pdata[pcpt,13]=Pdata[:,13][Pdata[:,0]==new_ip[pcpt]][0] # region unchanged
            for tcpt in np.arange(new_tip):#ternaries
                tcount=tcount+1
                new_Tdata[tcount,0:9]=ternaries[tcpt,0:9] # ternary number, type, lat, lon, radius, area, alt., geoid height, density unchanged
                new_Tdata[tcount,9]=pcpt+1 # renumber associated primary
                new_Tdata[tcount,10]=pcpt+1 # renumber associated secondary
                new_Tdata[tcount,11:13]=ternaries[tcpt,11:13] # pcolour, scolour and region unchanged
                new_Tdata[tcount,13]=Pdata[:,13][Pdata[:,0]==new_ip[pcpt]][0]
            del ternaries, new_tip, pcland
            
        ## Create new headers
        print('')
        print('Create headers:',str(datetime.datetime.now())[0:19])
        npex=int(new_nbp)
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
    