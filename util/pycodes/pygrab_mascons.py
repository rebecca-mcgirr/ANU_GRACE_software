#!/usr/bin/env python3
import datetime
import numpy as np
import pyproj
import shapely.geometry as geometry
import shapely.ops as ops
from functools import partial
import argparse
import os
###############################################################################
## This programs aim to extract mascons within regions defined by polygons.
## It can either identified the primaries inside polygons (--primflag=1) or form new primaries to match polygons (--primflag=0, default).
## Mandatory input contains an input mascon file and a polygon file.
## You can choose a custom output mascon file name (option --out_masconfile, default = input mascon file + '_grab_' + region) 
## Tou can change the default of differemt options such as density, buffer distance, threshold on the size of primaries and tide mask to be applied to new primaries.
## To get more information about arguments and syntax, type \"pygrab_mascon -h\". 
## Author : Julia Pfeffer, 18 October 2019 
###############################################################################
## DEFINE FUNCTIONS TO READ INPUT FILES
###############################################################################
import grace_utils as gr       
###############################################################################
## READ ARGUMENTS    
############################################################################### 
# Define arguments   
parser = argparse.ArgumentParser(description=' This programs aim to extract mascons within regions defined by polygons. It can run in default mode (--primflag=0) or display mode (--primflag=1). In default mode, ternaries are rearranged into new primaries matching the geometry of polygons. In display mode, information about primaries inside polygons (number and area) is displayed on screen (primaries are not reshaped and you  don\'t generate a new mascon file) . Mandatory arguments include an input mascon file and a polygon file. You can choose an output mascon filename, the default  is the input mascon filename + \"_grab_\" + region (default=\'all\'). You can change defaults options including the region, density, buffer distance, threshold and/or tidal mask.  Author : Julia Pfeffer, 18 October 2019 ')
parser.add_argument("in_masconfile",default="None",type=str,help="complete path to your input mascon file (string) ")
parser.add_argument("polygon_file",default="None",type=str,help="complete path to your ascii polygon file (string)")
parser.add_argument("--out_masconfile",default="None",type=str,help="complete path to your output mascon file (default = input mascon filename + \"_grab_\" + region)")
parser.add_argument("--region",default="all",type=str,help="select one region from polygon file (default = all)")
parser.add_argument("--density",default=1000,type=int,help="default = 1000 for land (oceans = 1029)")
parser.add_argument("--bufferdist",default=0.25,type=float,help="margin distance around the bounding box wrapping polygons, default = 0.25 degrees")
parser.add_argument("--threshold",default=10,type=int,help="number of ternaries in the smallest primary outside all polygons, but inside the bounding box. Primaries smaller than that will be merged with the closest primary of identical density.  Default = 10")
parser.add_argument("--tidemask",default=0,type=int,help="bitmap integer defining the tidal correction to be be applied on new primaries, default = 0 (no tide correction)")
parser.add_argument("--primflag",default=0,type=int,help="Find primaries whose centroids are inside the polygons. When activated (--primflag==1), primaries inside polygons are displayed on screen, but we don't generate any  masconfile. Default = 0.")
# Read arguments
args = parser.parse_args()
# Declare arguments
input_mascon_file = args.in_masconfile
myregion = args.region
if args.out_masconfile== "None":
    out_masconfile=input_mascon_file+'_grab_'+myregion
else:    
    out_masconfile = args.out_masconfile
polygon_file = args.polygon_file

mydensity = args.density
mybufferdist = args.bufferdist
mythreshold = args.threshold
mytidemask = args.tidemask
myflag=args.primflag
############################################################################### 
## READ INPUT FILES
###############################################################################
print('')
print('#######################################################################################################')
print('############################## START PYGRAB_MASCONS:',str(datetime.datetime.now())[0:19],'##############################')
print('#######################################################################################################')      
print('')
print('Read input files:',str(datetime.datetime.now())[0:19])
# Read mascon file
iheaders,iPdata,iSdata,iTdata=gr.read_mascon_file(input_mascon_file)
nbt=len(iTdata)
# Read polygon file
regions,coords,nbv=gr.read_polygons_ascii(polygon_file,myregion)
nbpol=len(coords)
############################################################################## 
## DEFINE A BOUNDING BOX AROUND POLYGONS 
##############################################################################
print('Define bounding box around polygons:',str(datetime.datetime.now())[0:19])
minlat,maxlat,minlon,maxlon=gr.boundingbox(coords,mybufferdist)
if minlon<0: # if bounding box crosses the Greenwitch meridian, convert positive longitudes to negatives longitudes
    iPdata[:,5][iPdata[:,5]>180]=iPdata[:,5][iPdata[:,5]>180]-360*np.ones(len(iPdata[:,5][iPdata[:,5]>180]))
    iTdata[:,3][iTdata[:,3]>180]=iTdata[:,3][iTdata[:,3]>180]-360*np.ones(len(iTdata[:,3][iTdata[:,3]>180]))
    coords=gr.convert_longitudes_in_tuple(coords,flag='to_negative')  
if myflag==1: # display mode only
    Pmask=(iPdata[:,4]>=minlat)&(iPdata[:,4]<=maxlat)&(iPdata[:,5]<=maxlon)&(iPdata[:,5]>=minlon)&(iPdata[:,10]==mydensity) # defime mask with bounding box 
    boxPdata=iPdata[Pmask==True] # ternaries inside bounding box with density = mydensity
    nbpbox=len(boxPdata[:,0]) # number of ternaries inside the box
    print('Found %d primaries in bounding box:'%(nbpbox),str(datetime.datetime.now())[0:19])
############################################################################### 
## GRAB PRIMARY MASCON INSIDE POLYGONS
############################################################################### 
    print('')
    print('Find primaries inside polygons:',str(datetime.datetime.now())[0:19])
    id_primary,id_polygon=gr.select_primaries_in_polygons(boxPdata,regions,coords,print_output='yes') # return indices of  primaries whose centroid is inside polygon and display information on screen
    ## select mascons with region name in polygon file
    id_primary=np.asarray(id_primary)
    inPdata=boxPdata[id_primary]
    nbpin=len(id_primary)
    # select ternaries associated with primaries
    ternaries=[]
    for pcpt in np.arange(nbpin):
        ternaries.append(iTdata[iTdata[:,9]==inPdata[pcpt,0]])
    inTdata=np.asarray(ternaries)
    nbtin=len(inTdata)
    del boxPdata,Pmask,nbpbox, iheaders, iPdata,iSdata,iTdata,ternaries,nbt,id_primary,id_polygon
################################################################################
## MAP NEW MASCONS 
#################################################################################        
    print('')
    print('Create map of primaries inside polygons(map_pygrab_mascons_%s.png):'%(myregion),str(datetime.datetime.now())[0:19])
    import matplotlib
    havedisplay = "DISPLAY" in os.environ
    if not havedisplay:
        matplotlib.use("Agg")
    from mpl_toolkits.basemap import Basemap
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Polygon
    import matplotlib.pyplot as plt
    import math
    #### colorbar
    colors=[]
    cm20=plt.cm.Set1
    for cpt in np.arange(math.ceil(nbpol)/cm20.N):
        for cpt2 in np.arange(cm20.N):
            colors.append(cm20(cpt2))
    cm=matplotlib.colors.ListedColormap(colors)     
    ### figure 
    fig = plt.figure(figsize=(15,15))
    plt.cla()
    plt.clf()
    ax=fig.gca();
    m = Basemap(llcrnrlon=minlon,llcrnrlat=minlat, urcrnrlon=maxlon, urcrnrlat=maxlat, resolution = 'l')
    Px, Py = m(inPdata[:,5],inPdata[:,4])
    patches = []
    for pcpt in np.arange(len(coords)):
        patches.append(Polygon(coords[pcpt]))
    ax.add_collection(PatchCollection(patches, facecolor='none', edgecolor='k', linewidths=0.75,zorder=1))
    m.drawcoastlines(zorder=1)
    for cpt in np.arange(len(inPdata)):
        Tx, Ty = m(inTdata[cpt][:,3],inTdata[cpt][:,2])
        plt.scatter(Tx,Ty,5,inTdata[cpt][:,9],cmap=cm,vmin=np.min(inPdata[:,0]),vmax=np.max(inPdata[:,0]),edgecolors='face',zorder=2)
        t=plt.text(Px[cpt]-0.75,Py[cpt],inPdata[cpt,0],color='k',fontweight='bold',fontsize=10)
        t.set_bbox(dict(facecolor='white', alpha=0.6, edgecolor='none'))
    plt.savefig('map_pygrab_mascons_%s.png'%(myregion))        
    del minlat,maxlat,minlon,maxlon, coords, regions,nbv, patches
    del inPdata,inTdata,nbpin,nbtin,colors,cm20,Px,Py,Tx,Ty                 
else:
############################################################################## 
## SELECT TERNARIES INSIDE BOUNDING BOX TO SAVE TIME
##############################################################################
    Tmask=(iTdata[:,2]>=minlat)&(iTdata[:,2]<=maxlat)&(iTdata[:,3]<=maxlon)&(iTdata[:,3]>=minlon)&(iTdata[:,8]==mydensity) # defime mask with bounding box and density
    boxTdata=iTdata[Tmask==True] # ternaries inside bounding box with density = mydensity
    nbtbox=len(boxTdata) # number of ternaries inside the box
    print('Found %d ternaries in bounding box:'%(nbtbox),str(datetime.datetime.now())[0:19])
############################################################################### 
## GRAB TERNARY MASCON INSIDE POLYGONS
############################################################################### 
    print('')
    print('Find which polygons contain which ternaries. This can take a few minutes, have tea! ',str(datetime.datetime.now())[0:19])
    inout=np.zeros(nbtbox) # logical array defining if ternary is inside at least one polygon
    dist2pol=-99999*np.ones((nbtbox,nbpol)) # distance from ternary to edge of polygon
    id_polygon=-9999*np.ones(nbtbox) # number of polygon containing ternary or closest polygon
    wgs84=pyproj.Proj(init='epsg:4326') 
    ## For each ternary find polygon
    for tcpt in np.arange(nbtbox):
        plat=boxTdata[tcpt,2]
        plon=boxTdata[tcpt,3]
        # project point in azimuthal equidistant projection centered on point so that distances are preserved
        pproj = pyproj.Proj(proj="aeqd", lat_0=plat, lon_0=plon, datum="WGS84", units="m")
        project = partial(pyproj.transform,pproj,wgs84)
        mypoint= geometry.Point((plon,plat))
        mypoint = ops.transform(project, mypoint)
        allpol=np.zeros(nbpol)  # logical array defining if one polygon contains a ternary
        for pcpt in np.arange(nbpol):# check if point is contained in one of the polygons
            # project polygon in azimuthal equidistant projection
            mypoly=geometry.Polygon(coords[pcpt])
            mypoly = ops.transform(project, mypoly)
            dist2pol[tcpt,pcpt]=mypoly.exterior.distance(mypoint)# compute distance ternary to edge of polygon
            if mypoly.contains(mypoint):
                allpol[pcpt]=1 # set logical to 1 if point inside polygon
        inout[tcpt]=np.sum(allpol)
        if inout[tcpt]>0:  # find if point is inside one polygon     
            id_polygon[tcpt]=np.where(allpol==1)[0][0] # find the polygon containing the point
        elif inout[tcpt]==0: # if point is outside polygon
            id_polygon[tcpt]=np.where(dist2pol[tcpt,:]==np.min(dist2pol[tcpt,:]))[0][0]# find the closest polygon
    del mypoly, mypoint, dist2pol,allpol, tcpt, pcpt,plat,plon,pproj,wgs84 ,project     
    print('Found %d ternaries rearranged in %d polygons'%(len(boxTdata[inout>0]),nbpol))
###############################################################################
### FIND SMALL PRIMARIES IN BOUNDING BOX AND MERGE WITH CLOSEST PRIMARY
###############################################################################
    print('')
    print('Find small primaries with all ternaries outside polygons and at least one ternary inside bounding box :',str(datetime.datetime.now())[0:19])
    # Put together all ternaries the polygons
    outpol_Tdata=np.concatenate((iTdata[Tmask==False],boxTdata[inout==0]),axis=0) 
    outpol_Tprim=outpol_Tdata[:,9] # associated primary number for ternaries outside the polygons      
    ## Find small primaries outside polygons in bounding box
    small_prim=0
    for tcpt in np.arange(nbtbox):
        if inout[tcpt]==0: # ternary outside polygons inside bounding box
            tip_out=len(outpol_Tprim[outpol_Tprim==boxTdata[tcpt,9]]) # number of ternary in associated primary
            if tip_out<mythreshold: # if the number of ternary in associated primary < threshold (defined as input, default = 10) 
                inout[tcpt]=1 # put ternary inside closest polygon
                small_prim=small_prim+1
    print('Found %d small primaries (less than %d ternaries per primary) merged with the closest primary of identical density'%(small_prim,mythreshold))
    #del outpol_Tdata, outpol_Tprim, tip_out, small_prim
    ## Divide ternaries inside and outside polygons
    in_Tdata=boxTdata[inout>0] # ternaries inside polygons
    in_Tprim=id_polygon[inout>0] # new associated primary number (=temporary) = number of polygon in polygon file
    out_Tdata=np.concatenate((iTdata[Tmask==False],boxTdata[inout==0])) # ternaries outside polygons
    del Tmask, boxTdata, tcpt, nbtbox, id_polygon, inout, mythreshold
###############################################################################
## CREATE NEW PRIMARY AND TERNARY ARRAYS
#################################################################################
    print('')
    print('Create new primary and ternary arrays. This can take a few minutes. Relax and finish your tea :-)',str(datetime.datetime.now())[0:19])
    ##### Count the primaries 
    ip_out=np.unique(out_Tdata[:,9])#  list of unique primaries outside polygons
    nbp_out=len(ip_out) # number of primaries outside polygons
    ip_in=np.unique(in_Tprim) #  list of unique primaries inside polygons
    nbp_in=len(ip_in) # number of primaries inside polygons
    new_nbp=nbp_out+nbp_in # total number of primaries
    ##### Build empty primary and ternary arrays
    new_Pdata=np.array(np.zeros((new_nbp,14)),dtype=object)
    new_Tdata=np.array(np.zeros((nbt,14)),dtype=object)
    ### Fill arrays for primaries and ternaies outside polygons
    tcount=-1
    for pcpt in np.arange(nbp_out):# primaries
        ternaries=out_Tdata[out_Tdata[:,9]==ip_out[pcpt]] # select ternaries outside polygons
        new_tip=np.shape(ternaries)[0] # number of ternary inside one primary 
        pcland=100*len(ternaries[:,8][ternaries[:,8]==1000])/new_tip # percentage of land in primary 
        new_Pdata[pcpt,0]=pcpt+1 # renumber the primary in ascending order
        new_Pdata[pcpt,1]=iPdata[:,1][iPdata[:,0]==ip_out[pcpt]][0] # type of primary unchanged
        new_Pdata[pcpt,2]=1
        new_Pdata[pcpt,3]=new_tip # number of ternary in primary
        new_Pdata[pcpt,4]=iPdata[:,4][iPdata[:,0]==ip_out[pcpt]][0]# latitude unchanged
        new_Pdata[pcpt,5]=iPdata[:,5][iPdata[:,0]==ip_out[pcpt]][0]# longitude unchanged
        new_Pdata[pcpt,6]=iPdata[:,6][iPdata[:,0]==ip_out[pcpt]][0]# radius unchanged
        new_Pdata[pcpt,7]=iPdata[:,7][iPdata[:,0]==ip_out[pcpt]][0] # area unchanged
        new_Pdata[pcpt,8]=iPdata[:,8][iPdata[:,0]==ip_out[pcpt]][0]# altitude unchanged 
        new_Pdata[pcpt,9]=iPdata[:,9][iPdata[:,0]==ip_out[pcpt]][0]# geoid hei unchangedght
        new_Pdata[pcpt,10]=iPdata[:,10][iPdata[:,0]==ip_out[pcpt]][0] # density unchanged
        new_Pdata[pcpt,11]=pcland
        new_Pdata[pcpt,12]=iPdata[:,12][iPdata[:,0]==ip_out[pcpt]][0] # tidal mask unchanged
        new_Pdata[pcpt,13]=iPdata[:,13][iPdata[:,0]==ip_out[pcpt]][0] # region unchanged
        for tcpt in np.arange(new_tip):#ternaries
            tcount=tcount+1
            new_Tdata[tcount,0:9]=ternaries[tcpt,0:9] # ternary number, type, lat, lon, radius, area, alt., geoid height, density unchanged
            new_Tdata[tcount,9]=pcpt+1 # renumber associated primary
            new_Tdata[tcount,10]=pcpt+1 # renumber associated secondary
            new_Tdata[tcount,11:14]=ternaries[tcpt,11:14] # pcolour, scolour and region unchanged
        del ternaries, new_tip, pcland
    del out_Tdata, ip_out   
    ### Fill arrays for primaries and ternaies inside polygons
    for pcpt in np.arange(nbp_out,new_nbp):#primaries
        ternaries=in_Tdata[in_Tprim==ip_in[pcpt-nbp_out]]# select ternaries for each primary inside polygon
        new_tip=np.shape(ternaries)[0] # number of ternaries inside primary
        pcland=100*len(ternaries[:,8][ternaries[:,8]==1000])/new_tip # percentage of land in primary
        new_Pdata[pcpt,0]=pcpt+1 # renumber primary
        if mydensity == 1000:
            new_Pdata[pcpt,1]='PLand' # new primary typee
        else:
            new_Pdata[pcpt,1]='PDeep'
        new_Pdata[pcpt,2]=1 #number of secondary in primary
        new_Pdata[pcpt,3]=new_tip #number of ternary in primary
        new_Pdata[pcpt,4]=np.mean(ternaries[:,2])# new latitude
        new_Pdata[pcpt,5]=np.mean(ternaries[:,3])# new longitude
        new_Pdata[pcpt,6]=np.mean(ternaries[:,4])# new radius
        new_Pdata[pcpt,7]=np.sum(ternaries[:,5]) # new area
        new_Pdata[pcpt,8]=np.mean(ternaries[:,6])# new altitude 
        new_Pdata[pcpt,9]=np.mean(ternaries[:,7])# new geoid height
        new_Pdata[pcpt,10]=mydensity #  density defined as a selection criteria
        new_Pdata[pcpt,11]=pcland
        new_Pdata[pcpt,12]=mytidemask # tidalmask defined in input (default = 0)
        new_Pdata[pcpt,13]=regions[int(ip_in[pcpt-nbp_out])]# new region name defined in polygon file
        for tcpt in np.arange(new_tip):
            tcount=tcount+1
            new_Tdata[tcount,0]=ternaries[tcpt,0] # number of ternary unchanged
            if mydensity == 1000:
                new_Tdata[tcount,1]='TLand' # ternary type
            else:
                new_Tdata[tcount,1]='TDeep'
            new_Tdata[tcount,2:9]=ternaries[tcpt,2:9] # lat, lon, radius, area, alt., geoid height and density unchanged
            new_Tdata[tcount,9]=pcpt+1 # renumber associated primary
            new_Tdata[tcount,10]=pcpt+1 # renumber associated secondary
            new_Tdata[tcount,11:13]=ternaries[tcpt,11:13] # pcolour and scolour unchanged
            new_Tdata[tcount,13]=regions[int(ip_in[pcpt-nbp_out])] # new region name defined in polygon file
        del ternaries, new_tip, pcland    
    del iPdata,iSdata,iTdata    
    del in_Tdata, in_Tprim, ip_in, nbp_out, nbp_in, tcount, pcpt,tcpt, mydensity,mytidemask
    new_Pdata[:,5][new_Pdata[:,5]<0]=new_Pdata[:,5][new_Pdata[:,5]<0]+360*np.ones(len(new_Pdata[:,5][new_Pdata[:,5]<0]))
    new_Tdata[:,3][new_Tdata[:,3]<0]=new_Tdata[:,3][new_Tdata[:,3]<0]+360*np.ones(len(new_Tdata[:,3][new_Tdata[:,3]<0]))
###############################################################################
## NEW HEADERS
#################################################################################
    print('')
    print('Create headers:',str(datetime.datetime.now())[0:19])
    npex=int(new_nbp)
    ntex=int(nbt)
    maxtip=int(np.max(new_Pdata[:,3]))
    lastcomment='# pygrab_mascons.py, extract %s in  %s. Timetag: %s \n'%(myregion,polygon_file,str(datetime.datetime.now())[0:19])
    new_headers=gr.define_new_mascon_headers(npex,ntex,maxtip,iheaders,lastcomment)
    del iheaders, ntex,maxtip,npex,lastcomment
#################################################################################
### WRITE NEW MASCON FILE
#################################################################################
    print('Write mascon in file %s:'%(out_masconfile),str(datetime.datetime.now())[0:19])
    gr.write_mascon_file(out_masconfile,new_headers,new_Pdata,new_Tdata)
    del new_Pdata,new_Tdata,new_headers      
    print('')
    print('#######################################################################################################')
    print('######## SUGGESTION: run repair_mascons in ~/gt/util on %s ########'%(out_masconfile))
    print('#######################################################################################################')         
    print('## This will generate a valid hashcode and recompute the latitude, longitude, radius, area, altitude ##')
    print('######## and geoid height of each primary and secondary mascon to be consistent with ternaries ########')
    print('#######################################################################################################')
print('')
print('#######################################################################################################')
print('############################ FINISH PYGRAB_MASCONS:',str(datetime.datetime.now())[0:19],'###############################')
print('#######################################################################################################')      
#################################################################################
