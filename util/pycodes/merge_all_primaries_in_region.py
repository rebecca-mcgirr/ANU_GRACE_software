#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 12:03:10 2020

@author: julia pfeffer
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
    parser.add_argument("--region",default="all",type=str,help="name of the region in which to merge primaries")
    
    # Read arguments
    args = parser.parse_args()
    # Declare arguments
    infile = args.in_masconfile
    outfile = args.out_masconfile
    if outfile=="None":
        outfile=infile+'_rmerged'
    region = args.region   
    
    return infile,outfile,region

def main():
    #start program
    print('~/gt/util/pycodes/merge_all_primaries_in_region.py started at %s'%(str(datetime.datetime.now())[0:19]))

    # read arguments
    infile,outfile,region=read_arguments()
    
    #create temporary file
    tempfile=outfile+'_temp'
    
    # Create temporary mascon file to write in
    tempfile=outfile+'_temp'
    os.system("cp %s %s"%(infile,tempfile))
    
    # Read input masconfile
    print('')
    print('Read mascon file: %s'%(str(datetime.datetime.now())[0:19]))   
    headers,Pdata,_,Tdata=gr.read_mascon_file(infile)
    nbt=len(Tdata[:,0])
    
    # Merge all regions
    if region=='all':
        
        ## get all regions
        regions=np.unique(Pdata[:,13])
        nbr=len(regions)
        
        ##### Build empty primary and ternary arrays
        new_Pdata=np.array(np.zeros((nbr,14)),dtype=object)
        new_Tdata=np.array(np.zeros((nbt,14)),dtype=object)
        
        ## fill primary and ternary arrays
        print('')
        print('Create new primary and ternary arrays. This can take a few minutes. Relax and finish your tea :-)',str(datetime.datetime.now())[0:19])

        tcount=-1
        for pct in np.arange(nbr):# primaries
            ternaries=Tdata[Tdata[:,13]==regions[pct]] # select ternaries outside polygons
            new_tip=np.shape(ternaries)[0] # number of ternary inside one primary 
            pcland=100*len(ternaries[:,8][ternaries[:,8]==1000])/new_tip # percentage of land in primary 
            new_Pdata[pct,0]=pct+1 # renumber the primary in ascending order
            new_Pdata[pct,1]=Pdata[:,1][Pdata[:,13]==regions[pct]][0] # type of primary unchanged
            new_Pdata[pct,2]=1
            new_Pdata[pct,3]=new_tip # number of ternary in primary
            new_Pdata[pct,4]=Pdata[:,4][Pdata[:,13]==regions[pct]][0]# latitude unchanged
            new_Pdata[pct,5]=Pdata[:,5][Pdata[:,13]==regions[pct]][0]# longitude unchanged
            new_Pdata[pct,6]=Pdata[:,6][Pdata[:,13]==regions[pct]][0]# radius unchanged
            new_Pdata[pct,7]=Pdata[:,7][Pdata[:,13]==regions[pct]][0] # area unchanged
            new_Pdata[pct,8]=Pdata[:,8][Pdata[:,13]==regions[pct]][0]# altitude unchanged 
            new_Pdata[pct,9]=Pdata[:,9][Pdata[:,13]==regions[pct]][0]# geoid hei unchangedght
            new_Pdata[pct,10]=Pdata[:,10][Pdata[:,13]==regions[pct]][0] # density unchanged
            new_Pdata[pct,11]=pcland
            new_Pdata[pct,12]=Pdata[:,12][Pdata[:,13]==regions[pct]][0] # tidal mask unchanged
            new_Pdata[pct,13]=Pdata[:,13][Pdata[:,13]==regions[pct]][0] # region unchanged
            for tcpt in np.arange(new_tip):#ternaries
                tcount=tcount+1
                new_Tdata[tcount,0:9]=ternaries[tcpt,0:9] # ternary number, type, lat, lon, radius, area, alt., geoid height, density unchanged
                new_Tdata[tcount,9]=pct+1 # renumber associated primary
                new_Tdata[tcount,10]=pct+1 # renumber associated secondary
                new_Tdata[tcount,11:14]=ternaries[tcpt,11:14] # pcolour, scolour and region unchanged
            del ternaries, new_tip, pcland
    else:
        # if only one region to process
        inTdata=Tdata[Tdata[:,13]==region] # select ternaries inside region
        outTdata=Tdata[Tdata[:,13]!=region] # select ternaries outside region
        
        ##### Count the primaries 
        ip_out=np.unique(outTdata[:,9])#  list of unique primaries outside southern ocean
        nbp_out=len(ip_out) # number of primaries outside southern ocean
        new_nbp=nbp_out+1 # total number of primaries

        ##### Build empty primary and ternary arrays
        new_Pdata=np.array(np.zeros((new_nbp,14)),dtype=object)
        new_Tdata=np.array(np.zeros((nbt,14)),dtype=object)
        
        ### Fill arrays for primaries and ternaries outside region
        print('')
        print('Create new primary and ternary arrays. This can take a few minutes. Relax and finish your tea :-)',str(datetime.datetime.now())[0:19])

        tcount=-1
        for pct in np.arange(nbp_out):# primaries
            ternaries=outTdata[outTdata[:,9]==ip_out[pct]] # select ternaries outside polygons
            new_tip=np.shape(ternaries)[0] # number of ternary inside one primary 
            pcland=100*len(ternaries[:,8][ternaries[:,8]==1000])/new_tip # percentage of land in primary 
            new_Pdata[pct,0]=pct+1 # renumber the primary in ascending order
            new_Pdata[pct,1]=Pdata[:,1][Pdata[:,0]==ip_out[pct]][0] # type of primary unchanged
            new_Pdata[pct,2]=1
            new_Pdata[pct,3]=new_tip # number of ternary in primary
            new_Pdata[pct,4]=Pdata[:,4][Pdata[:,0]==ip_out[pct]][0]# latitude unchanged
            new_Pdata[pct,5]=Pdata[:,5][Pdata[:,0]==ip_out[pct]][0]# longitude unchanged
            new_Pdata[pct,6]=Pdata[:,6][Pdata[:,0]==ip_out[pct]][0]# radius unchanged
            new_Pdata[pct,7]=Pdata[:,7][Pdata[:,0]==ip_out[pct]][0] # area unchanged
            new_Pdata[pct,8]=Pdata[:,8][Pdata[:,0]==ip_out[pct]][0]# altitude unchanged 
            new_Pdata[pct,9]=Pdata[:,9][Pdata[:,0]==ip_out[pct]][0]# geoid hei unchangedght
            new_Pdata[pct,10]=Pdata[:,10][Pdata[:,0]==ip_out[pct]][0] # density unchanged
            new_Pdata[pct,11]=pcland
            new_Pdata[pct,12]=Pdata[:,12][Pdata[:,0]==ip_out[pct]][0] # tidal mask unchanged
            new_Pdata[pct,13]=Pdata[:,13][Pdata[:,0]==ip_out[pct]][0] # region unchanged
            for tcpt in np.arange(new_tip):#ternaries
                tcount=tcount+1
                new_Tdata[tcount,0:9]=ternaries[tcpt,0:9] # ternary number, type, lat, lon, radius, area, alt., geoid height, density unchanged
                new_Tdata[tcount,9]=pct+1 # renumber associated primary
                new_Tdata[tcount,10]=pct+1 # renumber associated secondary
                new_Tdata[tcount,11:14]=ternaries[tcpt,11:14] # pcolour, scolour and region unchanged
            del ternaries, new_tip, pcland
        del outTdata, ip_out  
        
        ### Fill arrays for primaries and ternaries in Southern Ocean
        ternaries=inTdata# select ternaries for each primary inside polygon
        new_tip=np.shape(ternaries)[0] # number of ternaries inside primary
        pcland=100*len(ternaries[:,8][ternaries[:,8]==1000])/new_tip # percentage of land in primary
        new_Pdata[pct+1,0]=nbp_out+1 # renumber primary
        new_Pdata[pct+1,1]=Pdata[:,1][Pdata[:,13]==region][0] # type of primary unchanged
        new_Pdata[pct+1,2]=1 #number of secondary in primary
        new_Pdata[pct+1,3]=new_tip #number of ternary in primary
        new_Pdata[pct+1,4]=np.mean(ternaries[:,2])# new latitude
        new_Pdata[pct+1,5]=np.mean(ternaries[:,3])# new longitude
        new_Pdata[pct+1,6]=np.mean(ternaries[:,4])# new radius
        new_Pdata[pct+1,7]=np.sum(ternaries[:,5]) # new area
        new_Pdata[pct+1,8]=np.mean(ternaries[:,6])# new altitude 
        new_Pdata[pct+1,9]=np.mean(ternaries[:,7])# new geoid height
        new_Pdata[pct+1,10]=Pdata[:,10][Pdata[:,13]==region][0]#  density unchanged
        new_Pdata[pct+1,11]=pcland
        new_Pdata[pct+1,12]=Pdata[:,12][Pdata[:,13]==region][0]#  density unchanged
        new_Pdata[pct+1,13]=region # new region name
        for tcpt in np.arange(new_tip):
            tcount=tcount+1
            new_Tdata[tcount,0:9]=ternaries[tcpt,0:9] # number of ternary unchanged
            new_Tdata[tcount,9]=nbp_out+1  # renumber associated primary
            new_Tdata[tcount,10]=nbp_out+1  # renumber associated secondary
            new_Tdata[tcount,11:13]=ternaries[tcpt,11:13] # pcolour and scolour unchanged
            new_Tdata[tcount,13]=region # new region name defined in polygon file
        del ternaries, new_tip, pcland 
    
    ## Create new headers
    print('')
    print('Create headers:',str(datetime.datetime.now())[0:19])
    npex=int(len(new_Pdata))
    ntex=int(nbt)
    maxtip=int(np.max(new_Pdata[:,3]))
    lastcomment='# merge_all_primaries_in_region.py : all primaries of %s region merged per region.  Timetag: %s \n'%(region,str(datetime.datetime.now())[0:19])
    new_headers=gr.define_new_mascon_headers(npex,ntex,maxtip,headers,lastcomment)

    ## Write mascon file
    print('Write mascon in file %s:'%(outfile),str(datetime.datetime.now())[0:19])
    gr.write_mascon_file(outfile,new_headers,new_Pdata,new_Tdata)
    del new_Pdata,new_Tdata,new_headers     
  
        

if __name__ == '__main__':
    main() 
