#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 11:36:22 2020

@author: julia
"""
import os
import argparse
import grace_utils as gr
import numpy as np
import datetime

def read_arguments():
    # Define arguments   
    parser = argparse.ArgumentParser(description=' This programs aim to grab ternary mascons inside polygons.')
    parser.add_argument("in_masconfile",default="/scratch/compute1/julia/masconfiles/mascons_stage4",type=str,help="complete path to your mascon file")
    parser.add_argument("--area",default=90000000000,type=int,help="resolution of primary mascon in m2")
    parser.add_argument("--outputfile",default="None",type=str,help="complete path to your output mascon file")
    parser.add_argument("--overflow",default=10,type=int,help="maximum diiference between mascons")
    parser.add_argument("--maxsize",default=3000000,type=int,help="maximum size of the MCF problem")
    parser.add_argument("--precision",default=7500,type=int,help="maximum distance displacement")
    # Read arguments
    args = parser.parse_args()
    # Declare arguments
    infile = args.in_masconfile
    area = args.area
    outfile = args.outputfile
    if outfile=='None':
        outfile=infile+'_stage5'
    maxsize = args.maxsize
    overflow = args.overflow
    precision = args.precision
    return infile,outfile,area,maxsize,overflow,precision

def main():
    
    #get current directory
    mydir=os.getcwd()
    if mydir[-1]=='/':
        pass
    else:
        mydir=mydir+'/'
    
    #create status file    
    status_file=mydir+'mascons_stage5.status'
    sfile=open(status_file, 'w')
    sfile.write('started create_mascons_stage5_with_status_file.py. Timetag: %s. \n'%(str(datetime.datetime.now())[0:19]))
 
    
    #Read arguments
    try:
        infile,outfile,area,maxsize,overflow,precision=read_arguments()
    except:
        sfile.write('Error in reading arguments. Timetag: %s. \n'%(str(datetime.datetime.now())[0:19]))
    else:
        sfile.write('read arguments: infile=%s, outfile=%s, area=%d,maxsize=%d, overflow=%d, precision=%d. Timeflag=%s. \n'%(infile,outfile,area,maxsize,overflow,precision,str(datetime.datetime.now())[0:19]))    

    
    #create temporary file
    try:
        tempfile=outfile+'_temp'
        os.system("cp %s %s"%(infile,tempfile))
    except:
        sfile.write('Error in creating temporary file. Timetag: %s. \n'%(str(datetime.datetime.now())[0:19]))
    else:
        sfile.write('temporary file %s created at %s\n'%(tempfile,str(datetime.datetime.now())[0:19]))
     
    
    #Read input masconfile
    try:
        _,Pdata,_,_=gr.read_mascon_file(tempfile)
        regions=np.unique(Pdata[:,13])
        nbr=len(regions)
        del Pdata
    except:
        sfile.write('Error. Cannot read %s to extract regions. Timetag: %s'%(tempfile,str(datetime.datetime.now())[0:19]))
    else:
        sfile.write('Read %s to extract regions. Timetag: %s'%(tempfile,str(datetime.datetime.now())[0:19]))
    
    ## run pyshape mascons sequentially on all regions
    for cpt in np.arange(nbr):
        try:
            #print(regions[cpt])
            os.system("~/gt/util/pycodes/pyshape_constrained_size_mascons.py %s %d --overflow=%d --precision=%d --maxsize=%d --region=%s --out_masconfile=%s"%(tempfile,area,overflow,precision,maxsize,regions[cpt],tempfile))
        except:
            sfile.write("Error in ~/gt/util/pycodes/pyshape_constrained_size_mascons.py %s %d --overflow=%d --precision=%d --maxsize=%d --region=%s --out_masconfile=%s. Timetag: %s. \n"%(tempfile,area,overflow,precision,maxsize,regions[cpt],tempfile,str(datetime.datetime.now())[0:19]))
        else:
            sfile.write("~/gt/util/pycodes/pyshape_constrained_size_mascons.py %s %d --overflow=%d --precision=%d --maxsize=%d --region=%s --out_masconfile=%s. Timetag: %s. \n"%(tempfile,area,overflow,precision,maxsize,regions[cpt],tempfile,str(datetime.datetime.now())[0:19]))
    
    ## merge small islands to closest landmass if closer than 750 km, otherwise merge with ocean
    try:
        os.system("~/gt/util/pycodes/merge_small_primaries.py %s --outfile=%s --minsize=25 --maxdist=750 --ftype=1"%(tempfile,tempfile))
    except:
        sfile.write("Error in ~/gt/util/pycodes/merge_small_primaries.py %s --outfile=%s --minsize=25 --maxdist=750 --ftype=1. Timetag: %s. \n"%(tempfile,tempfile,str(datetime.datetime.now())[0:19]))
    else:
        sfile.write("~/gt/util/pycodes/merge_small_primaries.py %s --outfile=%s --minsize=25 --maxdist=750 --ftype=1. Timetag: %s. \n"%(tempfile,tempfile,str(datetime.datetime.now())[0:19]))
    
    ## merge small primary mascons to closest mascon of same type
    try:
        _,Pdata,_,_=gr.read_mascon_file(tempfile)
        minsize=np.median(Pdata[:,3])/2
        maxdist=2.5*np.sqrt(area)/1000
        os.system("~/gt/util/pycodes/merge_small_primaries.py %s --outfile=%s --minsize=%d --maxdist=%d --ftype=0"%(tempfile,tempfile,minsize,maxdist))
    except:
        sfile.write("Error in ~/gt/util/pycodes/merge_small_primaries.py %s --outfile=%s --minsize=%d --maxdist=%d --ftype=0. Timetag: %s. \n"%(tempfile,tempfile,minsize,maxdist,str(datetime.datetime.now())[0:19]))
    else:
        sfile.write("~/gt/util/pycodes/merge_small_primaries.py %s --outfile=%s --minsize=%d --maxdist=%d --ftype=0. Timetag: %s. \n"%(tempfile,tempfile,minsize,maxdist,str(datetime.datetime.now())[0:19]))
        del Pdata,minsize,maxdist          
    
    ## reshape regions with really large primaries
    try:
        os.system("~/gt/util/pycodes/reshape_large_primaries.py %s --outfile=%s --minarea=%d --maxarea=%d --area=%d"%(tempfile,tempfile,int(0.75*area),int(1.25*area),area))
    except:
        sfile.write("Error in ~/gt/util/pycodes/reshape_large_primaries.py %s --outfile=%s --minarea=%d --maxarea=%d --area=%d. Timetag: %s. \n"%(tempfile,tempfile,int(0.75*area),int(1.25*area),area,str(datetime.datetime.now())[0:19]))
    else:
        sfile.write("~/gt/util/pycodes/reshape_large_primaries.py %s --outfile=%s --minarea=%d --maxarea=%d --area=%d. Timetag: %s. \n"%(tempfile,tempfile,int(0.75*area),int(1.25*area),area,str(datetime.datetime.now())[0:19]))
        
    #repair masconfile
    try:
        lakefile='/scratch/compute1/geodynamics/software/grace/tables/lake_polygons.dat'
        os.system("cp %s %s "%(lakefile,mydir))
        os.system("~/gt/util/repair_mascons %s -lakefix"%(tempfile))
    except:
        sfile.write("Error in ~/gt/util/repair_mascons %s -lakefix. Timetag: %s. \n"%(tempfile,str(datetime.datetime.now())))
    else:
        sfile.write("~/gt/util/repair_mascons %s -lakefix. Timetag: %s. \n"%(tempfile,str(datetime.datetime.now())))        
        
    #move repaired masconfile as final output file
    os.system("mv %s %s "%(tempfile+'_repaired',outfile))

if __name__ == '__main__':
    main()     
