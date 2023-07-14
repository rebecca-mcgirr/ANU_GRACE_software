#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 14:55:24 2020

@author: julia
"""
import numpy as np
import datetime as dt
import grace_utils as gr
import os
import argparse


def read_arguments():
    
    # Define arguments
    parser = argparse.ArgumentParser('Separate mascon file in ocean and continents mascons at requested resolution')
    parser.add_argument('infile', default="None",  type=str, help='complete pah to your input masconfile')
    parser.add_argument('--outfile', default="None",  type=str, help='complete pah to your output masconfile (default = input file name + stage3)')
    parser.add_argument('--area', default=90000000000,  type=int, help='mascon area in m2 (default = 90 000 000 000 m2)')
    parser.add_argument('--region', default="all",  type=str, help='region to reshape : Land, Ocean or all (default=all)')
    parser.add_argument('--precision',default=10000,type=int,help='maximum mascon displacement allowed to stop reshaping process (default = 10 000 m)')
    parser.add_argument('--maxsize',default=9000000,type=int,help='maximum size of the minimum cost flow problem (default = 9 000 000 arcs)')
    parser.add_argument('--flag',default=0,type=int,help='set flag value to 1 to avoid reshaping process (default=0, allowing reshaping process)')
    
    ## Read arguments
    args = parser.parse_args()
    
    ## Declare arguments
    infile = args.infile
    outfile = args.outfile
    if outfile=='None':
        outfile=infile+'stage3'
    area = args.area
    precision = args.precision
    maxsize = args.maxsize
    flag = args.flag
    region = args.region
    
    return infile,outfile,area,region,precision,maxsize,flag

def separate_land_from_oceans(infile,outfile,flag):
    
    #read input masconfile
    headers,Pdata,Sdata,Tdata=gr.read_mascon_file(infile)
    nbp=len(Pdata)
    tcount=0
    
    # Separate ocean and continents based on mascon density
    for cpt in np.arange(nbp):
        if Pdata[cpt,10]==1000:
            Pdata[cpt,-1]='Land' #region
            Pdata[cpt,1]='PLand' #type
        else:
            Pdata[cpt,-1]='Ocean'#region
            Pdata[cpt,1]='PDeep' #type
        ntip=Pdata[cpt,3]
        for tpt in np.arange(ntip):
            if Tdata[cpt,8]==1000:
                Tdata[tcount,-1]='Land'#region
                Tdata[tcount,1]='TLand'#type
            else:
                Tdata[tcount,-1]='Ocean'#region
                Tdata[tcount,1]='TDeep' #type             
            tcount=tcount+1
    
    #Write headers in output masconfile
    npex=nbp
    ntex=len(Tdata)
    maxtip=int(np.max(Pdata[:,3]))
    headers=['#','# mascon file stage 3 generated from %s %s\n'%(infile,str(dt.datetime.now())[0:19])]
    lastcomment='# %s separated in ocean and continents based on density values \n'%(str(dt.datetime.now())[0:19])
    new_headers=gr.define_new_mascon_headers(npex,ntex,maxtip,headers,lastcomment)
    
    if flag==0:
        tempfile=infile+'temp'
    else:
        tempfile=outfile
    #Write mascon lines
    gr.write_mascon_file(tempfile,new_headers,Pdata,Tdata)
    
    #return name of temporary output file
    return tempfile
   

def main():
    # read arguments
    infile,outfile,area,region,precision,maxsize,flag=read_arguments()
    
    #separate land from oceans
    tempfile=separate_land_from_oceans(infile,outfile,flag)
    
    #reshape land/ocean mascon file
    if flag==0:
        if region=='all':
            os.system("~/gt/util/pycodes/pyshape_constrained_size_mascons.py %s %d --precision=%d --maxsize=%d --region=Land --out_masconfile=%s"%(tempfile,area,precision,maxsize,tempfile))
            os.system("~/gt/util/pycodes/pyshape_constrained_size_mascons.py %s %d --precision=%d --maxsize=%d --region=Ocean --out_masconfile=%s"%(tempfile,area,precision,maxsize,outfile))
        else:
            os.system("~/gt/util/pycodes/pyshape_constrained_size_mascons.py %s %d --precision=%d --maxsize=%d --region=%s --out_masconfile=%s"%(tempfile,area,precision,maxsize,region,outfile))
        
    #repair land/ocean mascon file    
    os.system("~/gt/util/repair_mascons %s"%(outfile))
    
    #rename repaired file
    os.system("mv %s %s"%(outfile+'_repaired',outfile))
    
if __name__ == '__main__':
    main()  
