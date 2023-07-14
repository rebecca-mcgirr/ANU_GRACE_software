#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 15:39:42 2020

@author: julia
"""
import argparse
import grace_utils as gr
import numpy as np
import datetime 
import os

def read_arguments():
    parser = argparse.ArgumentParser(description=' This programs aims to merge all primary mascons in a given region')
    parser.add_argument("in_masconfile",default="None",type=str,help="complete path to your input mascon file")
    parser.add_argument("--outfile",default="None",type=str,help="complete path to your output mascon file")
    parser.add_argument("--maxarea",default=105000000000,type=int,help="size under which a region is considered too small")
    parser.add_argument("--minarea",default=75000000000,type=int,help="maximum distance between regions to be merged (km)")    
    parser.add_argument("--area",default=90000000000,type=int,help="maximum distance between regions to be merged (km)")    

    # Read arguments
    args = parser.parse_args()
    # Declare arguments
    infile = args.in_masconfile
    outfile = args.outfile
    if outfile=="None":
        outfile=infile+'_merged'
    area = args.area    
    minarea = args.minarea
    maxarea = args.maxarea
    return infile,outfile,area,minarea,maxarea

def main():
    
    #read arguments
    infile,outfile,area,minarea,maxarea=read_arguments()
    
    # Read input masconfile
    print('')
    print('Read mascon file: %s'%(str(datetime.datetime.now())[0:19]))   
    headers,Pdata,_,Tdata=gr.read_mascon_file(infile)
    
    #Find large primaries
    Plarge=Pdata[Pdata[:,7]>maxarea]    
    lregions=np.unique(Plarge[:,13])
    nbr=len(lregions)
    print('Found %d regions with primaries over %f km2. Timetag: %s.'%(nbr,maxarea/1000000,str(datetime.datetime.now())[0:19]))
    
    #Reshape regions with large primaries with one more primary
    for rpt in np.arange(nbr):
        print('')
        ternaries=Tdata[Tdata[:,13]==lregions[rpt]]
        nbp_init=len(np.unique(ternaries[:,9]))
        new_nbp=nbp_init+1
        new_area=np.floor((np.sum(ternaries[:,5])/new_nbp))
        if new_area<minarea:
            print('%s cannot be reshaped to fit an area between %d and %d m2. To fit one more primary mascon, you would need to allow a minimum area of %d m2.'%(lregions[rpt],minarea,maxarea,int(np.floor(new_area))))
        else:
            print('%s  reshaped to fit an area of %d.'%(lregions[rpt],new_area))
            os.system("~/gt/util/pycodes/pyshape_constrained_size_mascons.py %s %d --region=%s --out_masconfile=%s"%(infile,new_area,lregions[rpt],infile))
    

if __name__ == '__main__':
    main()     

