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
    parser.add_argument('--apol', default="/scratch/compute1/geodynamics/software/grace/tables/australian_catchments_L1.polygons",  type=str, help='complete pah to your australian catchments polygon file (default = /scratch/compute1/geodynamics/software/grace/tables/australian_catchments_L1.polygons')
    parser.add_argument('--rbpol', default="/scratch/compute1/geodynamics/software/grace/tables/riverbasins_grdc.polygons",  type=str, help='complete pah to your river basin polygon file (default = /scratch/compute1/geodynamics/software/grace/tables/riverbasins_grdc.polygons')
    parser.add_argument('--lpol', default="/scratch/compute1/geodynamics/software/grace/tables/largestlakes.polygons",  type=str, help='complete pah to your river basin polygon file (default = /scratch/compute1/geodynamics/software/grace/tables/largestlakes.polygons')
    parser.add_argument('--fc', default=1,  type=int, help='flag to extract continents (default=1) or not (set to zero)')
    parser.add_argument('--fo', default=1,  type=int, help='flag to extract oceans (default=1) or not (set to zero)')
    parser.add_argument('--fi', default=1,  type=int, help='flag to extract islands (default=1) or not (set to zero)')
    parser.add_argument('--fs', default=1,  type=int, help='flag to extract seas (default=1) or not (set to zero)')
    parser.add_argument('--fa', default=1,  type=int, help='flag to extract australia (default=1) or not (set to zero)')
    parser.add_argument('--frb', default=1,  type=int, help='flag to extract river basins (default=1) or not (set to zero)')
    parser.add_argument('--fl', default=1,  type=int,help='flag to extract lakes (default=1) or not (set to zero)')

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
    
def main():
    while True:
        try:
            #code
            #get current directory
            mydir=os.getcwd()
            if mydir[-1]=='/':
                pass
            else:
                mydir=mydir+'/'
            
            #create status file    
            status_file=mydir+'mascons_stage4.status'
            sfile=open(status_file, 'w')
            sfile.write('started create_mascons_stage4.py. Timetag: %s. \n'%(str(datetime.datetime.now())[0:19]))
            
            # read arguments
            try:
                infile,outfile,continents_polfile,oceans_polfile,islands_polfile,seas_polfile,australian_polfile,riverbasins_polfile,lakes_polfile,fc,fo,fi,fs,fa,frb,fl=read_arguments()        
            except:
                sfile.write('Error in reading arguments. Timetag: %s. \n'%(str(datetime.datetime.now())[0:19]))
            else:
                sfile.write('read arguments: infile=%s, outfile=%s, continents_polfile=%s,oceans_polfile=%s, islands_polfile=%s, seas_polfile=%s, australian_polfile=%s, riverbasins_polfile=%s, lakes_polfile=%s, fc=%d, fo=%d, fi=%d,fs=%d, fa=%d, frb=%d, fl=%d. Timeflag=%s\n'%(infile,outfile,continents_polfile,oceans_polfile,islands_polfile,seas_polfile,australian_polfile,riverbasins_polfile,lakes_polfile,fc,fo,fi,fs,fa,frb,fl,str(datetime.datetime.now())[0:19]))    

            # Create temporary mascon file to write in
            tempfile=outfile+'_temp'
            try:
                os.system("cp %s %s"%(infile,tempfile))            
            except:
                sfile.write('Error in creating temporary file. Timetag: %s\n'%(str(datetime.datetime.now())[0:19]))
            else:
                sfile.write('temporary file %s created at %s\n'%(tempfile,str(datetime.datetime.now())[0:19]))

                
            #grab oceans
            if fo==1:
                oregions,ocoords,_=gr.read_polygons_ascii(oceans_polfile,'all')
                nbpol_o=len(ocoords)
                for cpt in np.arange(nbpol_o):
                    try:
                        os.system("~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1029 --region=%s"%(tempfile,oceans_polfile,tempfile,oregions[cpt]))    
       
                    except:
                        sfile.write("Error in ~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1029 --region=%s. Timetag: %s. \n"%(tempfile,oceans_polfile,tempfile,oregions[cpt],str(datetime.datetime.now())[0:19]))                
                    else:
                        sfile.write("~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1029 --region=%s. Timetag: %s. \n"%(tempfile,oceans_polfile,tempfile,oregions[cpt],str(datetime.datetime.now())[0:19]))
  
                del oregions,ocoords,nbpol_o
                
                #grab Southern ocean 
                try:
                    os.system("~/gt/util/pycodes/grab_southern_ocean.py %s --out_masconfile=%s "%(tempfile,tempfile))        
                except:
                    sfile.write("Error in ~/gt/util/pycodes/grab_southern_ocean.py %s --out_masconfile=%s. Timetag: %s. \n"%(tempfile,tempfile,str(datetime.datetime.now())))
                else:
                    sfile.write("~/gt/util/pycodes/grab_southern_ocean.py %s --out_masconfile=%s. Timetag: %s. \n"%(tempfile,tempfile,str(datetime.datetime.now())))
  
                #grab Arctic ocean
                try:
                    os.system("~/gt/util/pycodes/grab_arctic_ocean.py %s --out_masconfile=%s "%(tempfile,tempfile))  
                except:
                    sfile.write("Error in ~/gt/util/pycodes/grab_arctic_ocean.py %s --out_masconfile=%s. Timetag: %s. \n "%(tempfile,tempfile,str(datetime.datetime.now())))
                else:
                    sfile.write("~/gt/util/pycodes/grab_arctic_ocean.py %s --out_masconfile=%s. Timetag: %s. Timetag: %s. \n"%(tempfile,tempfile,str(datetime.datetime.now())))
            
            if fs==1:
                #grab seas and gulfs
                sregions,scoords,_=gr.read_polygons_ascii(seas_polfile,'all')
                nbpol_s=len(scoords)
                for cpt in np.arange(nbpol_s):
                    try:
                        os.system("~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1029 --region=%s"%(tempfile,seas_polfile,tempfile,sregions[cpt]))            
                    except:
                        sfile.write("Error in ~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1029 --region=%s. Timetag: %s. \n"%(tempfile,seas_polfile,tempfile,sregions[cpt],str(datetime.datetime.now())))
                    else:
                        sfile.write("~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1029 --region=%s. Timetag: %s. \n"%(tempfile,seas_polfile,tempfile,sregions[cpt],str(datetime.datetime.now())))

                del sregions,scoords,nbpol_s    
            
            if (fo==0)&(fs==0):
                pass
            else:
                #merge ocean ternaries that have not been associated with any ocean basin or sea
                try:
                    os.system("~/gt/util/pycodes/fix_solo_ocean_mascons.py %s --out_masconfile=%s "%(tempfile,tempfile))          
                except:
                    sfile.write("Error in ~/gt/util/pycodes/fix_solo_ocean_mascons.py %s --out_masconfile=%s. Timetag: %s. \n"%(tempfile,tempfile,str(datetime.datetime.now())))
                else:
                    sfile.write("~/gt/util/pycodes/fix_solo_ocean_mascons.py %s --out_masconfile=%s. Timetag: %s. \n"%(tempfile,tempfile,str(datetime.datetime.now())))       

            if fc==1:        
                #grab continents
                cregions,ccoords,_=gr.read_polygons_ascii(continents_polfile,'all')
                nbpol_c=len(ccoords)
                for cpt in np.arange(nbpol_c):
                    try:
                        os.system("~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1000 --region=%s"%(tempfile,continents_polfile,tempfile,cregions[cpt]))            
                    except:
                        sfile.write("Error in ~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1000 --region=%s. Timetag: %s. \n"%(tempfile,continents_polfile,tempfile,cregions[cpt],str(datetime.datetime.now())))  
                    else:
                        sfile.write("~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1000 --region=%s. Timetag: %s. \n"%(tempfile,continents_polfile,tempfile,cregions[cpt],str(datetime.datetime.now())))            

                del cregions,ccoords,nbpol_c 
                
                # Rename Eurasia_a as Eurasia
                try:
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'Eurasia_a','Eurasia'))          
                except:
                    sfile.write("Error in ~/gt/util/pycodes/replace_region_name.py %s %s %s %s. Timetag: %s. \n"%(tempfile,tempfile,'Eurasia_a','Eurasia',str(datetime.datetime.now())))
                else:
                    sfile.write("~/gt/util/pycodes/replace_region_name.py %s %s %s %s. Timetag: %s. \n"%(tempfile,tempfile,'Eurasia_a','Eurasia',str(datetime.datetime.now())))

                # Rename Eurasia_b as Eurasia
                try:
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'Eurasia_b','Eurasia'))          
                except:
                    sfile.write("Error in ~/gt/util/pycodes/replace_region_name.py %s %s %s %s. Timetag: %s. \n"%(tempfile,tempfile,'Eurasia_b','Eurasia',str(datetime.datetime.now())))
                else:
                    sfile.write("~/gt/util/pycodes/replace_region_name.py %s %s %s %s. Timetag: %s. \n"%(tempfile,tempfile,'Eurasia_b','Eurasia',str(datetime.datetime.now())))

                #grab antarctica
                try:
                    os.system("~/gt/util/pycodes/grab_antarctica.py %s --out_masconfile=%s "%(tempfile,tempfile))          
                except:
                    sfile.write("Error in ~/gt/util/pycodes/grab_antarctica.py %s --out_masconfile=%s. Timetag: %s. \n "%(tempfile,tempfile,str(datetime.datetime.now())))        
                else:
                   sfile.write("~/gt/util/pycodes/grab_antarctica.py %s --out_masconfile=%s. Timetag: %s. \n "%(tempfile,tempfile,str(datetime.datetime.now())))       

            if fi==1:
                #grab islands and peninsulas  
                iregions,icoords,_=gr.read_polygons_ascii(islands_polfile,'all')
                nbpol_i=len(icoords)
                for cpt in np.arange(nbpol_i):
                    try:
                        os.system("~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1000 --region=%s"%(tempfile,islands_polfile,tempfile,iregions[cpt]))            
                    except:
                        sfile.write("Error in ~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1000 --region=%s. Timetag: %s. \n"%(tempfile,islands_polfile,tempfile,iregions[cpt],str(datetime.datetime.now())))                          
                    else:
                        sfile.write("~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1000 --region=%s. Timetag: %s. \n"%(tempfile,islands_polfile,tempfile,iregions[cpt],str(datetime.datetime.now())))

                del iregions,icoords,nbpol_i
                
                try:
                    # Rename Fiji_E as Fiji
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'Fiji_E','Fiji'))                  
                    # Rename Fiji_W as Fiji
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'Fiji_W','Fiji'))           
                    # Rename EBeringIs as BeringIs
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'EBeringIs','BeringIs'))                  
                    # Rename WBeringIs as BeringIs
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'WBeringIs','BeringIs'))                          
                except:
                    sfile.write("Error in renaming Fiji and Bering islands. Timetag: %s. \n"%(str(datetime.datetime.now())))
                else:
                    sfile.write("Renaming Fiji_E and Fiji_W as Fiji. Renaming EBering_Is and W_Bering_is as Bering islands. Timetag: %s. \n"%(str(datetime.datetime.now())))
                

            
            if (fc==0)&(fi==0):
                pass
            else:
                #fix land primaries that have not been associated with any continent or islands
                try:
                    os.system("~/gt/util/pycodes/fix_solo_land_mascons.py %s --out_masconfile=%s "%(tempfile,tempfile))          
                except:
                    sfile.write("Error in ~/gt/util/pycodes/fix_solo_land_mascons.py %s --out_masconfile=%s. Timetag: %s. \n"%(tempfile,tempfile,str(datetime.datetime.now())))  
                else:
                    sfile.write("~/gt/util/pycodes/fix_solo_land_mascons.py %s --out_masconfile=%s. Timetag: %s. \n"%(tempfile,tempfile,str(datetime.datetime.now())))  

            if frb==1:
                #grab river basins
                rbregions,rbcoords,_=gr.read_polygons_ascii(riverbasins_polfile,'all')
                nbpol_rb=len(rbcoords)
                for cpt in np.arange(nbpol_rb):
                    try:
                        os.system("~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1000 --region=%s"%(tempfile,riverbasins_polfile,tempfile,rbregions[cpt]))            
                    except:
                        sfile.write("Error in ~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1000 --region=%s. Timetag: %s. \n"%(tempfile,riverbasins_polfile,tempfile,rbregions[cpt],str(datetime.datetime.now())))
                    else:
                        sfile.write("~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1000 --region=%s. Timetag: %s. \n"%(tempfile,riverbasins_polfile,tempfile,rbregions[cpt],str(datetime.datetime.now())))
                del rbregions,rbcoords,nbpol_rb
            
            if fl==1:
                #grab large lakes
                lregions,lcoords,_=gr.read_polygons_ascii(lakes_polfile,'all')
                nbpol_l=len(lcoords)
                for cpt in np.arange(nbpol_l):
                    try:
                        os.system("~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1000 --region=%s"%(tempfile,lakes_polfile,tempfile,lregions[cpt]))         
                    except:
                        sfile.write("Error in ~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1000 --region=%s. Timetag: %s. \n"%(tempfile,lakes_polfile,tempfile,lregions[cpt],str(datetime.datetime.now())))          
                    else:
                        sfile.write("~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1000 --region=%s. Timetag: %s. \n"%(tempfile,lakes_polfile,tempfile,lregions[cpt],str(datetime.datetime.now())))
                del lregions,lcoords,nbpol_l
            
            if fa==1:
                #grab australian catchments
                try:
                    os.system("~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1000 --bufferdist=0.5 --threshold=75"%(tempfile,australian_polfile,tempfile))    
                    # Rename AUS01a,b and c as AUS01
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'AUS01a','AUS01'))  
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'AUS01b','AUS01'))       
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'AUS01c','AUS01'))
                    # Rename AUS02a,b,c and d as AUS02
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'AUS02a','AUS02'))  
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'AUS02b','AUS02'))       
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'AUS02c','AUS02'))
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'AUS02d','AUS02'))
                    # Rename AUS05a and b as AUS05
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'AUS05a','AUS05'))  
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'AUS05b','AUS05')) 
                    # Rename AUS07a,b and c as AUS07
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'AUS07a','AUS07'))  
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'AUS07b','AUS07'))       
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'AUS07c','AUS07'))
                    # Rename AUS08a and b as AUS08
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'AUS08a','AUS08'))  
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'AUS08b','AUS08')) 
                    # Rename AUS10a and b as AUS10
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'AUS10a','AUS10'))  
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'AUS10b','AUS10'))        
                    # Rename AUS11a,b,c and d as AUS11
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'AUS11a','AUS11'))  
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'AUS11b','AUS11'))       
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'AUS11c','AUS11'))
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'AUS11d','AUS11'))
                    # Rename AUS12a and b as AUS12
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'AUS12a','AUS12'))  
                    os.system("~/gt/util/pycodes/replace_region_name.py %s %s %s %s "%(tempfile,tempfile,'AUS12b','AUS12'))            
                except:
                    sfile.write("~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1000 --bufferdist=0.5 --threshold=75. Timetag:%s."%(tempfile,australian_polfile,tempfile,str(datetime.datetime.now())))
                else:
                    sfile.write("~/gt/util/pycodes/pygrab_mascons.py %s %s --out_masconfile=%s --density=1000 --bufferdist=0.5 --threshold=75. Timetag:%s."%(tempfile,australian_polfile,tempfile,str(datetime.datetime.now())))
                                 
            #make one primary per region
            try:
                os.system("~/gt/util/pycodes/merge_all_primaries_in_region.py %s --out_masconfile=%s "%(tempfile,tempfile))               
            except:
                sfile.write("Error in ~/gt/util/pycodes/merge_all_primaries_in_region.py %s --out_masconfile=%s. Timetag: %s. \n "%(tempfile,tempfile,str(datetime.datetime.now())))
            else:
                sfile.write("~/gt/util/pycodes/merge_all_primaries_in_region.py %s --out_masconfile=%s. Timetag: %s. \n "%(tempfile,tempfile,str(datetime.datetime.now())))
            
            #move ternaries of small regions to different primaries
            try:
                os.system("~/gt/util/pycodes/move_ternaries_in_bigger_region.py %s --out_masconfile=%s --minsize=10 --maxdist=1000 --ftype=1"%(tempfile,tempfile))              
            except:
                sfile.write("Error in ~/gt/util/pycodes/move_ternaries_in_bigger_region.py %s --out_masconfile=%s --minsize=10 --maxdist=1000 --ftype=1. Timetag: %s. \n"%(tempfile,tempfile,str(datetime.datetime.now())))
            else:       
                sfile.write("~/gt/util/pycodes/move_ternaries_in_bigger_region.py %s --out_masconfile=%s --minsize=10 --maxdist=1000 --ftype=1. Timetag: %s. \n"%(tempfile,tempfile,str(datetime.datetime.now())))
 
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
            sfile.close()

        except KeyboardInterrupt:
            break

        
if __name__ == '__main__':
    main()  
