#!/usr/bin/env python3
import numpy as np
import datetime
from sklearn.cluster import KMeans
import pyproj 
import sys
import argparse
from random import randint
###############################################################################
import grace_utils as gr
###############################################################################
## READ ARGUMENTS    
############################################################################### 
## Define arguments   
parser = argparse.ArgumentParser(description=' This program aims to reshape primary mascons to fit the designated area. Mandatory arguments include an input mascon file and an area in square meters. Three different options are available to select the primary mascons to be reshaped : (1) select all primaries inside the polygons of a polygon file (you can choose which polygon, default =\'all\'), (2) select all primaries with one given region name, (3) select one primary by number. You can choose an output filename, the default is your input filename + \'_reshaped\'. Author : Julia Pfeffer, 22 October 2019 ')
parser.add_argument("in_masconfile",default="None",type=str,help="complete path to your input mascon file (string) ")
parser.add_argument("area",default=4e10,type=float,help="target area for reshaped mascons (m^2)")
parser.add_argument("--out_masconfile",default="None",type=str,help="complete path to your output mascon file (default = in_masconfile + \"_reshape\")")
parser.add_argument("--polygon_file",default="None",type=str,help="complete path to your ascii polygon file (default = \"None\")")
parser.add_argument("--polygon",default="all",type=str,help="select one polygon from polygon file (default = \"all\")")
parser.add_argument("--region",default="None",type=str,help="select one region from mascon file (default = \"None\")")
parser.add_argument("--primary",default=0,type=int,help="select one primary from mascon file per number (default = 0, no selection)")
parser.add_argument("--density",default=0,type=int,help="only select primary mascons with given density. Only useful when performing a geographic selection (select primaries inside polygons, default = 0, no selection)")
## Read arguments
args = parser.parse_args()
## Declare arguments
input_mascon_file = args.in_masconfile
myarea = args.area
if args.out_masconfile== "None":
    out_masconfile=input_mascon_file+'_reshaped'
else:    
    out_masconfile = args.out_masconfile
mydensity = args.density
polygon_file = args.polygon_file
mypolygon=args.polygon
myregion = args.region
mypnumber = args.primary
############################################################################### 
## READ INPUT FILES
###############################################################################
print('')
print('#######################################################################################################')
print('############################# START PYSHAPE_MASCONS:',str(datetime.datetime.now())[0:19],'##############################')
print('#######################################################################################################')      
print('')
print('Read input files:',str(datetime.datetime.now())[0:19])
iheaders,iPdata,_,iTdata=gr.read_mascon_file(input_mascon_file)# Read mascon file
nbt=len(iTdata)
############################################################################### 
## SELECT PRIMARY MASCONS TO RESHAPE
###############################################################################
print('Select primary mascons to reshape:',str(datetime.datetime.now())[0:19])
if polygon_file=='None':
    if myregion=='None':
        if mypnumber==0:
            print('No primary has been selected to be reshaped. Please provide a valid polygon file, region name or primary mascon number.')
            sys.exit(1)
        else:
            inPdata=iPdata[iPdata[:,0]==mypnumber]
            nbp_in=len(inPdata)
            outPdata=iPdata[iPdata[:,0]!=mypnumber]
            nbp_out=len(outPdata)
            if nbp_in==0:
                print('No primary mascon associated with number %d in %s'%(mypnumber,input_mascon_file))
                sys.exit(1)
            else:
                print('Found %d primaries corresponding to primary number %d:'%(int(nbp_in),mypnumber),str(datetime.datetime.now())[0:19])
    else:
        inPdata=iPdata[iPdata[:,13]==myregion]
        nbp_in=len(inPdata)  
        outPdata=iPdata[iPdata[:,13]!=myregion]
        nbp_out=len(outPdata)
        if nbp_in==0:
            print('No primary mascon associated with region %d in %s'%(myregion,input_mascon_file))
            sys.exit(1)
        else:
            print('Found %d primaries corresponding to region %s in %s:'%(int(nbp_in),myregion,input_mascon_file),str(datetime.datetime.now())[0:19])
else:    
    polygons,coords,nbv=gr.read_polygons_ascii(polygon_file,mypolygon)
    nbpol=len(coords)
    print('Define bounding box around polygons:',str(datetime.datetime.now())[0:19])
    minlat,maxlat,minlon,maxlon=gr.boundingbox(coords) 
    if minlon<0: # if bounding box crosses the Greenwitch meridian, convert positive longitudes to negatives longitudes
        iPdata[:,5][iPdata[:,5]>180]=iPdata[:,5][iPdata[:,5]>180]-360*np.ones(len(iPdata[:,5][iPdata[:,5]>180]))
        iTdata[:,3][iTdata[:,3]>180]=iTdata[:,3][iTdata[:,3]>180]-360*np.ones(len(iTdata[:,3][iTdata[:,3]>180]))
        coords=gr.convert_longitudes_in_tuple(coords,flag='to_negative')
    if mydensity ==0:
        Pmask=(iPdata[:,4]>=minlat)&(iPdata[:,4]<=maxlat)&(iPdata[:,5]<=maxlon)&(iPdata[:,5]>=minlon) # defime mask with bounding box
    else:
        Pmask=(iPdata[:,4]>=minlat)&(iPdata[:,4]<=maxlat)&(iPdata[:,5]<=maxlon)&(iPdata[:,5]>=minlon) &(iPdata[:,10]==mydensity)# defime mask with bounding box & density
    boxPdata=iPdata[Pmask==True] # ternaries inside bounding box with density = mydensity
    nbpbox=len(boxPdata[:,0]) # number of ternaries inside the box
    print('Found %d primaries in bounding box:'%(nbpbox),str(datetime.datetime.now())[0:19])
    del minlat,maxlat,minlon,maxlon
    print('Find primaries inside polygons:',str(datetime.datetime.now())[0:19])
    id_primary,id_polygon=gr.select_primaries_in_polygons(boxPdata,polygons,coords,print_output='no')# Select primaries in polygons
    inout = np.zeros(nbpbox)
    inout[id_primary] = 1 # Define mask where 1 = inside polygon and 0 = outside polygon            
    inPdata=boxPdata[inout==1] # Define array of primaries inside polygon
    nbp_in=len(inPdata[:,0])
    outPdata=np.concatenate((iPdata[Pmask==False],boxPdata[inout==0]))# Define array of primaries outside polygons 
    nbp_out=len(outPdata[:,0])
    del boxPdata,nbpbox,inout,Pmask,iPdata,polygons,coords,id_primary,id_polygon
    if nbp_in==0:
        print('No primary mascon associated with polygon %d from %s in %s'%(mypolygon,polygon_file,input_mascon_file))
        sys.exit(1)
    else:
        print('Found %d primaries inside polygons:'%(int(nbp_in)),str(datetime.datetime.now())[0:19]) 
################################################################################
## RESHAPE SELECTED PRIMARIES TO FIT INPUT AREA  
################################################################################
print('')
print('Reshape primaries inside polygons to fit an area of about %d km^2:'%(int(myarea/1e6)),str(datetime.datetime.now())[0:19])
print('')
nbt_in=np.sum(inPdata[:,3])
inTdata=np.array(np.zeros((nbt_in,14)),dtype=object)
new_Plon=[]
new_Plat=[]
new_tap=np.array(np.zeros(nbt_in),dtype=object)
tcount=0
pcount=nbp_out+1
wgs84=pyproj.Proj(init='epsg:4326') 
for pcpt in np.arange(nbp_in):  
    ip=inPdata[pcpt,0]
    ntip=inPdata[pcpt,3]
    inTdata[tcount:tcount+ntip,:]=iTdata[iTdata[:,9]==ip]
    # project point in azimuthal equidistant projection centered on point so that distances are preserved
    pproj = pyproj.Proj(proj="aeqd", lat_0=inPdata[pcpt,4], lon_0=inPdata[pcpt,5], datum="WGS84", units="m")
    X,Y=pyproj.transform(wgs84,pproj,iTdata[:,3][iTdata[:,9]==ip],iTdata[:,2][iTdata[:,9]==ip])
    tcoords=np.ones((ntip,2))
    tcoords[:,0]=X
    tcoords[:,1]=Y
    # Compute number of clusters
    ternary_area=np.mean(iTdata[:,5][iTdata[:,9]==ip]) # average area of a ternary inside that region
    ntip_ideal=int(round(myarea/ternary_area)) # number of ternary inside a primary of ideal size
    new_nbp=int(round(ntip/ntip_ideal)) # number of reshaped primaries 
    if new_nbp==0:
        new_nbp=1     
    #print('number of primaries in that region %d'%(nb_prim))
    if new_nbp>1:
        kmeans = KMeans(n_clusters=new_nbp, init='k-means++')
        new_prim = kmeans.fit_predict(tcoords)  
        for ccpt in np.arange(new_nbp):
            Plon,Plat=pyproj.transform(pproj,wgs84,kmeans.cluster_centers_[:, 0][ccpt],kmeans.cluster_centers_[:, 1][ccpt])
            new_Plon.append(Plon)
            new_Plat.append(Plat)
        del Plon,Plat,kmeans
    elif new_nbp==1:
        new_prim=np.zeros(ntip)
        new_Plon.append(inPdata[pcpt,5])
        new_Plat.append(inPdata[pcpt,4])  
    new_tap[tcount:tcount+ntip]=new_prim+pcount*np.ones(ntip)
    pcount=pcount+new_nbp
    tcount=tcount+ntip
    print('Primary number %d (area = %d km^2) reshaped into %d primary mascons '%(ip,int(inPdata[:,7][inPdata[:,0]==ip]/1e6),new_nbp))  
    del ip,ntip,pproj,X,Y,tcoords,ternary_area,ntip_ideal,new_nbp,new_prim
del wgs84    
######## Count the new primaries 
new_ip=np.unique(new_tap).astype(int)
new_nbp_in=len(new_ip)
################################################################################
## CREATE TERNARY AND PRIMARY ARRAYS
################################################################################
print('')
print('Create new primary and ternary arrays. This can take a few minutes:',str(datetime.datetime.now())[0:19])
######## Build empty primary and ternary arrays
new_Pdata=np.array(np.zeros((new_nbp_in+nbp_out,14)),dtype=object)
new_Tdata=np.array(np.zeros((nbt,14)),dtype=object)
### Fill arrays for primaries and ternaies outside polygons
tcount=-1
for pcpt in np.arange(nbp_out):# primaries
    ternaries=iTdata[iTdata[:,9]==outPdata[pcpt,0]] # select ternaries outside polygons
    out_tip=outPdata[pcpt,3] # number of ternary inside one primary 
    new_Pdata[pcpt,0]=pcpt+1 # renumber the primary in ascending order
    new_Pdata[pcpt,1:14]=outPdata[pcpt,1:14] # rest of primary data unchanged
    for tcpt in np.arange(out_tip):#ternaries
        tcount=tcount+1
        new_Tdata[tcount,0:9]=ternaries[tcpt,0:9] # ternary number, type, lat, lon, radius, area, alt., geoid height, density unchanged
        new_Tdata[tcount,9]=pcpt+1 # renumber associated primary
        new_Tdata[tcount,10]=pcpt+1 # renumber associated secondary
        new_Tdata[tcount,11:14]=ternaries[tcpt,11:14] # pcolour, scolour and region unchanged
    del ternaries, out_tip
del outPdata,iTdata       
for pcpt in np.arange(new_nbp_in):#primaries
    ternaries=inTdata[new_tap==new_ip[pcpt]]# select ternaries for each primary inside polygon
    in_tip=len(ternaries[:,0]) # number of ternaries inside primary
    pcland=100*len(ternaries[:,8][ternaries[:,8]==1000])/in_tip # percentage of land in primary
    new_Pdata[pcpt+nbp_out,0]=int(new_ip[pcpt]) # renumber primary
    new_Pdata[pcpt+nbp_out,1]=inPdata[:,1][inPdata[:,0]==ternaries[0,9]][0] # new primary typee
    new_Pdata[pcpt+nbp_out,2]=1 #number of secondary in primary
    new_Pdata[pcpt+nbp_out,3]=in_tip #number of ternary in primary
    new_Pdata[pcpt+nbp_out,4]=new_Plat[pcpt]# new latitude
    new_Pdata[pcpt+nbp_out,5]=new_Plon[pcpt]# new longitude
    new_Pdata[pcpt+nbp_out,6]=np.mean(ternaries[:,4])# new radius
    new_Pdata[pcpt+nbp_out,7]=np.sum(ternaries[:,5]) # new area
    new_Pdata[pcpt+nbp_out,8]=np.mean(ternaries[:,6])# new altitude 
    new_Pdata[pcpt+nbp_out,9]=np.mean(ternaries[:,7])# new geoid height
    new_Pdata[pcpt+nbp_out,10]=inPdata[:,10][inPdata[:,0]==ternaries[0,9]][0]  #  density defined as a selection criteria
    new_Pdata[pcpt+nbp_out,11]=pcland
    new_Pdata[pcpt+nbp_out,12]=inPdata[:,12][inPdata[:,0]==ternaries[0,9]][0] # tidalmask defined in input (default = 0)
    new_Pdata[pcpt+nbp_out,13]=inPdata[:,13][inPdata[:,0]==ternaries[0,9]][0] # new region name defined in polygon file
    # RM200114
    if new_Pdata[pcpt+nbp_out,10] < 1010:
        col_range=[750,1500]
        prim_colour=col_range[0] + randint(0,col_range[1]-col_range[0]) 
    elif new_Pdata[pcpt+nbp_out,10] > 1010:
        col_range=[0,450]
        prim_colour=col_range[0] + randint(0,col_range[1]-col_range[0]) 
    for tcpt in np.arange(in_tip):
        tcount=tcount+1
        new_Tdata[tcount,0:9]=ternaries[tcpt,0:9] # number of ternary, type, lat, lon, radius, area, altitude, height, density unchanged
        new_Tdata[tcount,9]=int(new_ip[pcpt]) # renumber associated primary
        new_Tdata[tcount,10]=int(new_ip[pcpt]) # renumber associated secondary
        new_Tdata[tcount,11:13]=[prim_colour,prim_colour]  # RM200114: randomly assign new pcolour and scolour
        new_Tdata[tcount,13]=ternaries[tcpt,13] # region unchanged
    del ternaries, in_tip, pcland    
del inTdata, inPdata, tcount,pcpt,tcpt  
new_Pdata[:,5][new_Pdata[:,5]<0]=new_Pdata[:,5][new_Pdata[:,5]<0]+360*np.ones(len(new_Pdata[:,5][new_Pdata[:,5]<0]))
new_Tdata[:,3][new_Tdata[:,3]<0]=new_Tdata[:,3][new_Tdata[:,3]<0]+360*np.ones(len(new_Tdata[:,3][new_Tdata[:,3]<0]))  
##################################################################################
##### Define new headers
# first line
npex=nbp_out+new_nbp_in
ntex=nbt
maxtip=int(np.max(new_Pdata[:,3]))
#last line
if polygon_file=='None':
    if myregion=='None':
        if mypnumber==0:
            print('No primary has been selected to be reshaped. Please provide a valid polygon file, region name or primary mascon number.')
        else:
            newline='# pyshape_mascons.py, reshape mascon number %s in  %s to fit an area of %d km^2. Timetag: %s \n'%(mypnumber,input_mascon_file,int(myarea/1e6),str(datetime.datetime.now())[0:19])       
    else:
        newline='# pyshape_mascons.py, reshape mascons associated with region %s in  %s to fit an area of %d km^2. Timetag: %s \n'%(myregion,input_mascon_file,int(myarea/1e6),str(datetime.datetime.now())[0:19])
else:
    newline='# pyshape_mascons.py, reshape mascons associated with %s polygon in  %s to fit an area of %d km^2. Timetag: %s \n'%(mypolygon,polygon_file,int(myarea/1e6),str(datetime.datetime.now())[0:19]) 
#all lines
new_headers=gr.define_new_mascon_headers(npex,ntex,maxtip,iheaders,newline)                   
del maxtip,iheaders, newline,npex,ntex,nbp_out,nbp_in,nbt
###################################################################################
# Write headers
gr.write_mascon_file(out_masconfile,new_headers,new_Pdata,new_Tdata)
del new_headers, new_Pdata,new_Tdata
################################################################################## 
## Map
print('')
if polygon_file=='None':
    if myregion=='None':
        if mypnumber==0:
            print('No primary has been selected to be reshaped. Please provide a valid polygon file, region name or primary mascon number.')
        else:
            print('Create map of former pimary %d (map_pyshape_mascons_%d.png):'%(mypnumber,mypnumber),str(datetime.datetime.now())[0:19])
            gr.map_mascons_per_number(out_masconfile,new_ip,figname='pyshape_mascons_%d.png'%(mypnumber))
                
    else:
        print('Create map of new primaries (map_pyshape_mascons_%s.png):'%(myregion),str(datetime.datetime.now())[0:19])
        gr.map_mascons_in_region(out_masconfile,myregion,figname='pyshape_mascons_%s.png'%(myregion))
else:
    print('Create map of new primaries (map_pyshape_mascons_%s.png):'%(mypolygon),str(datetime.datetime.now())[0:19])
    gr.map_mascons_and_polygons(out_masconfile,polygon_file,mypolygon=mypolygon,figname='pyshape_mascons_%s.png'%(mypolygon))    
###################################################################################
print('')
print('#######################################################################################################')
print('############## SUGGESTION: run repair_mascons in ~/gt/util on %s #################'%(out_masconfile))
print('#######################################################################################################')         
print('## This will generate a valid hashcode and recompute the latitude, longitude, radius, area, altitude ##')
print('######## and geoid height of each primary and secondary mascon to be consistent with ternaries ########')
print('#######################################################################################################')
print('')
print('#######################################################################################################')
print('############################# END PYSHAPE_MASCONS:',str(datetime.datetime.now())[0:19],'################################')
print('#######################################################################################################')      
#################################################################################    
  
