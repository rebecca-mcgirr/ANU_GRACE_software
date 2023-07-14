#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thursday 18 Jun 2020 14:55:24 

@author: julia pfeffer
"""

import os
import argparse
import grace_utils as gr
import numpy as np
import datetime


def read_arguments():
    
    # Define arguments
    parser = argparse.ArgumentParser('Divide mascon file in regions')
    parser.add_argument('infile', default="None",  type=str, help='complete pah to your input masconfile')
    parser.add_argument('--outfile', default="None",  type=str, help='complete pah to your output masconfile')
    parser.add_argument('--cpol', default="/scratch/compute1/geodynamics/software/grace/tables/continents.polygons",  type=str, help='complete pah to your continents polygon file (default = /scratch/compute1/geodynamics/software/grace/tables/continents.polygons')
    parser.add_argument('--opol', default="/scratch/compute1/geodynamics/software/grace/tables/oceans.polygons",  type=str, help='complete pah to your oceans polygon file (default = /scratch/compute1/geodynamics/software/grace/tables/oceans.polygons')
    parser.add_argument('--ipol', default="/scratch/compute1/geodynamics/software/grace/tables/islands_peninsulas.polygons",  type=str, help='complete pah to your islands and peninsulas polygon file (default = /scratch/compute1/geodynamics/software/grace/tables/islands_peninsulas.polygons')
    parser.add_argument('--spol', default="/scratch/compute1/geodynamics/software/grace/tables/seas_gulfs.polygons",  type=str, help='complete pah to your enclosed seas and gulfs polygon file (default = /scratch/compute1/geodynamics/software/grace/tables/seas_gulfs.polygons')
    parser.add_argument('--apol', default="/scratch/compute1/geodynamics/software/grace/tables/australian_catchments.polygons",  type=str, help='complete pah to your australian catchments polygon file (default = /scratch/compute1/geodynamics/software/grace/tables/australian_catchments_L1.polygons')
    parser.add_argument('--rbpol', default="/scratch/compute1/geodynamics/software/grace/tables/riverbasins_grdc.polygons",  type=str, help='complete pah to your river basin polygon file (default = /scratch/compute1/geodynamics/software/grace/tables/riverbasins_grdc.polygons')
    parser.add_argument('--lpol', default="/scratch/compute1/geodynamics/software/grace/tables/largestlakes.polygons",  type=str, help='complete pah to your river basin polygon file (default = /scratch/compute1/geodynamics/software/grace/tables/largestlakes.polygons')
    parser.add_argument('--fc', default=1,  type=bool, help='flag to extract continents (default=1) or not (set to zero)')
    parser.add_argument('--fo', default=1,  type=bool, help='flag to extract oceans (default=1) or not (set to zero)')
    parser.add_argument('--fi', default=1,  type=bool, help='flag to extract islands (default=1) or not (set to zero)')
    parser.add_argument('--fs', default=1,  type=bool, help='flag to extract seas (default=1) or not (set to zero)')
    parser.add_argument('--fa', default=1,  type=bool, help='flag to extract australia (default=1) or not (set to zero)')
    parser.add_argument('--frb', default=1,  type=bool, help='flag to extract river basins (default=1) or not (set to zero)')
    parser.add_argument('--fl', default=1,  type=bool,help='flag to extract lakes (default=1) or not (set to zero)')

    ## Read arguments
    args = parser.parse_args()
    
    ## Declare arguments
    infile = args.infile
    outfile = args.outfile
    if outfile=='None':
        outfile=infile+'_stage4'
    continents_polfile = args.cpol
    oceans_polfile = args.opol
    islands_polfile = args.ipol
    seas_polfile = args.spol
    australian_polfile = args.apol
    riverbasins_polfile = args.rbpol
    lakes_polfile = args.lpol
    fc = args.fc
    fo = args.fo
    fi = args.fi
    fs = args.fs
    fa = args.fa
    frb = args.frb
    fl = args.fl
    
    return infile,outfile,continents_polfile,oceans_polfile,islands_polfile,seas_polfile,australian_polfile,riverbasins_polfile,lakes_polfile,fc,fo,fi,fs,fa,frb,fl


def grab_antarctica(masconfile):
    print('')
    print('############################################################################################')
    print('######################### GRAB Antarctica:',str(datetime.datetime.now())[0:19],'#########################')
    print('##########################################################################################')      
    
    print('')
    print('Read mascon file:',str(datetime.datetime.now())[0:19])    
    headers,Pdata,_,Tdata=gr.read_mascon_file(masconfile)
    nbt=len(Tdata)
    
    print('')
    print('Select land ternaries below 60S:',str(datetime.datetime.now())[0:19])
    inTdata=Tdata[(Tdata[:,2]<=-60)&(Tdata[:,8]==1000)] # ternaries inside polygons
    outTdata=np.vstack((Tdata[(Tdata[:,2]>-60)],Tdata[(Tdata[:,2]<=-60)&(Tdata[:,8]==1029)]))

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
    new_Pdata[pcpt+1,2]=nbp_out+1 #number of secondary in primary
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
    new_Pdata[pcpt+1,13]='SOcean'# new region name defined in polygon file
    for tcpt in np.arange(new_tip):
        tcount=tcount+1
        new_Tdata[tcount,0:9]=ternaries[tcpt,0:9] # number of ternary unchanged
        new_Tdata[tcount,9]=nbp_out+1  # renumber associated primary
        new_Tdata[tcount,10]=nbp_out+1  # renumber associated secondary
        new_Tdata[tcount,11:13]=ternaries[tcpt,11:13] # pcolour and scolour unchanged
        new_Tdata[tcount,13]='SOcean' # new region name defined in polygon file
    del ternaries, new_tip, pcland    

    ## Create new headers
    print('')
    print('Create headers:',str(datetime.datetime.now())[0:19])
    npex=int(new_nbp)
    ntex=int(nbt)
    maxtip=int(np.max(new_Pdata[:,3]))
    lastcomment='# Define Southern Ocean as seawater mascons below 55 degree South. Timetag: %s \n'%(str(datetime.datetime.now())[0:19])
    new_headers=gr.define_new_mascon_headers(npex,ntex,maxtip,headers,lastcomment)

    ## Write mascon file
    print('Write mascon in file %s:'%(masconfile),str(datetime.datetime.now())[0:19])
    gr.write_mascon_file(masconfile,new_headers,new_Pdata,new_Tdata)
    del new_Pdata,new_Tdata,new_headers    


def grab_southern_ocean(masconfile):
    print('')
    print('############################################################################################')
    print('######################### GRAB SOUTHERN OCEAN:',str(datetime.datetime.now())[0:19],'#########################')
    print('##########################################################################################')      
    
    ## read mascon file
    print('')
    print('Read input files:',str(datetime.datetime.now())[0:19])    
    headers,Pdata,_,Tdata=gr.read_mascon_file(masconfile)
    nbt=len(Tdata)
    
    ## select ocean ternaries below 55S as Southern Ocean
    print('')
    print('Select ocean ternaries below 55S:',str(datetime.datetime.now())[0:19])
    inTdata=Tdata[(Tdata[:,2]<=-55)&(Tdata[:,8]==1029)] # ternaries inside polygons
    outTdata=np.vstack((Tdata[(Tdata[:,2]>-55)],Tdata[(Tdata[:,2]<=-55)&(Tdata[:,8]==1000)]))

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
    new_Pdata[pcpt+1,1]='PLand'
    new_Pdata[pcpt+1,2]=nbp_out+1 #number of secondary in primary
    new_Pdata[pcpt+1,3]=new_tip #number of ternary in primary
    new_Pdata[pcpt+1,4]=np.mean(ternaries[:,2])# new latitude
    new_Pdata[pcpt+1,5]=np.mean(ternaries[:,3])# new longitude
    new_Pdata[pcpt+1,6]=np.mean(ternaries[:,4])# new radius
    new_Pdata[pcpt+1,7]=np.sum(ternaries[:,5]) # new area
    new_Pdata[pcpt+1,8]=np.mean(ternaries[:,6])# new altitude 
    new_Pdata[pcpt+1,9]=np.mean(ternaries[:,7])# new geoid height
    new_Pdata[pcpt+1,10]=1000#  density defined as a selection criteria
    new_Pdata[pcpt+1,11]=pcland
    new_Pdata[pcpt+1,12]=0 # tidal mask (default = 0)
    new_Pdata[pcpt+1,13]='Antarctica'# new region name defined in polygon file
    for tcpt in np.arange(new_tip):
        tcount=tcount+1
        new_Tdata[tcount,0:9]=ternaries[tcpt,0:9] # number of ternary unchanged
        new_Tdata[tcount,9]=nbp_out+1  # renumber associated primary
        new_Tdata[tcount,10]=nbp_out+1  # renumber associated secondary
        new_Tdata[tcount,11:13]=ternaries[tcpt,11:13] # pcolour and scolour unchanged
        new_Tdata[tcount,13]='Antarctica' # new region name defined in polygon file
    del ternaries, new_tip, pcland    

    ## Create new headers
    print('')
    print('Create headers:',str(datetime.datetime.now())[0:19])
    npex=int(new_nbp)
    ntex=int(nbt)
    maxtip=int(np.max(new_Pdata[:,3]))
    lastcomment='# Define Antarctica as freshwater mascons below 60 degree South. Timetag: %s \n'%(str(datetime.datetime.now())[0:19])
    new_headers=gr.define_new_mascon_headers(npex,ntex,maxtip,headers,lastcomment)

    ## Write mascon file
    print('Write mascon in file %s:'%(masconfile),str(datetime.datetime.now())[0:19])
    gr.write_mascon_file(masconfile,new_headers,new_Pdata,new_Tdata)
    del new_Pdata,new_Tdata,new_headers 

def grab_arctic_ocean(masconfile):
    print('')
    print('############################################################################################')
    print('######################### GRAB ARCTIC OCEAN:',str(datetime.datetime.now())[0:19],'#########################')
    print('##########################################################################################')      
    
    ## read mascon file
    print('')
    print('Read input files:',str(datetime.datetime.now())[0:19])    
    headers,Pdata,_,Tdata=gr.read_mascon_file(masconfile)
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
    new_Pdata[pcpt+1,2]=nbp_out+1 #number of secondary in primary
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
    lastcomment='# Define Arctic as seawater mascons above 58 degree North. Timetag: %s \n'%(str(datetime.datetime.now())[0:19])
    new_headers=gr.define_new_mascon_headers(npex,ntex,maxtip,headers,lastcomment)

    ## Write mascon file
    print('Write mascon in file %s:'%(masconfile),str(datetime.datetime.now())[0:19])
    gr.write_mascon_file(masconfile,new_headers,new_Pdata,new_Tdata)
    del new_Pdata,new_Tdata,new_headers     
    
def main():
    # read arguments
    infile,outfile,continents_polfile,oceans_polfile,islands_polfile,seas_polfile,australian_polfile,riverbasins_polfile,lakes_polfile,fc,fo,fi,fs,fa,frb,fl=read_arguments()
    
    #create temporary file
    tempfile=outfile+'_temp'
    
    # Create temporary mascon file to write in
    tempfile=outfile+'_temp'
    os.system("cp %s %s"%(infile,tempfile))
        
    #grab oceans
    if fo==1:
        oregions,ocoords,_=gr.read_polygons_ascii(oceans_polfile,'all')
        nbpol_o=len(ocoords)
        for cpt in np.arange(nbpol_o):
           os.system("~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1029 --region=%s"%(tempfile,oceans_polfile,tempfile,oregions[cpt]))    
        del oregions,ocoords,nbpol_o
        
        #grab Southern ocean 
        grab_southern_ocean(tempfile)    
        
        #grab Arctic ocean
        grab_arctic_ocean(tempfile)
    
    if fs==1:
        #grab seas and gulfs
        sregions,scoords,_=gr.read_polygons_ascii(seas_polfile,'all')
        nbpol_s=len(scoords)
        for cpt in np.arange(nbpol_s):
            os.system("~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1029 --region=%s"%(tempfile,seas_polfile,tempfile,sregions[cpt]))    
        del sregions,scoords,nbpol_s    
        

    if fc==1:
        #grab continents
        cregions,ccoords,_=gr.read_polygons_ascii(continents_polfile,'all')
        nbpol_c=len(ccoords)
        for cpt in np.arange(nbpol_c):
            os.system("~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1000 --region=%s"%(tempfile,continents_polfile,tempfile,cregions[cpt]))    
        del cregions,ccoords,nbpol_c 
        
        #grab antarctica
        grab_antarctica(tempfile)
    
    if fi==1:
        #grab islands and peninsulas  
        iregions,icoords,_=gr.read_polygons_ascii(islands_polfile,'all')
        nbpol_i=len(icoords)
        for cpt in np.arange(nbpol_i):
            os.system("~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1000 --region=%s"%(tempfile,islands_polfile,tempfile,iregions[cpt]))    
        del iregions,icoords,nbpol_i
       
    if frb==1:
        #grab river basins
        rbregions,rbcoords,_=gr.read_polygons_ascii(riverbasins_polfile,'all')
        nbpol_rb=len(rbcoords)
        for cpt in np.arange(nbpol_rb):
            os.system("~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1000 --region=%s"%(tempfile,riverbasins_polfile,tempfile,rbregions[cpt]))    
        del rbregions,rbcoords,nbpol_rb
    
    if fl==1:
        #grab large lakes
        lregions,lcoords,_=gr.read_polygons_ascii(lakes_polfile,'all')
        nbpol_l=len(lcoords)
        for cpt in np.arange(nbpol_l):
            os.system("~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1000 --region=%s"%(tempfile,lakes_polfile,tempfile,lregions[cpt]))    
        del lregions,lcoords,nbpol_l
    
    if fa==1:
        #grab australian catchments
        os.system("~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1000 --bufferdist=2 --threshold=50"%(tempfile,australian_polfile,tempfile))    
    
    #repair masconfile
    os.system("~/gt/util/repair %s "%(tempfile))
    
    #move repaired masconfile as final output file
    os.system("mv %s %s "%(tempfile+'_repaired',outfile))
    
if __name__ == '__main__':
    main()  
