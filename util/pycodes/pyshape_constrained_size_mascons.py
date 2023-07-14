#!/scratch/compute1/geodynamics/pyenv/py27/bin/python2.7
import numpy as np
import datetime
import pyproj 
import sys
import argparse
from random import randint
###############################################################################
import utils_py27 as ut
###############################################################################
## READ ARGUMENTS    
############################################################################### 
## Define arguments   
parser = argparse.ArgumentParser(description=' This program aims to reshape mascons at a certain size. \
                                 The algorithm consists of regular K-Means implementation where cluster assignment is done by solving Minimum Cost Flow (MCF) problem, \
                                 as opposed to just using clusters with smallest distance measure to the center. \
                                 The MCF is formulated in such a way that data nodes (ternaries) are nodes with unitary supplies of flow. \
                                 The cluster centers (primaries) are nodes with a minimum flow demand equal to the minimum number of ternaries needed to match the resolution defined by user.\
                                 One additional node contains demand balancing the global sum of supply and demands to zero. The latter condition is necessary for feasibility of MCF solution. \
                                 The data nodes are fully connected to center nodes with edge costs corresponding to Euclidean distances from data points to respective cluster centers. \
                                 The data nodes are projected to X,Y coordinates with, an azimutal equidistant projection to preserve distances on the sphere. \
                                 The center nodes are fully connected to the balance node with zero cost, and a maximum capicity (overflow) corresponding to the difference in size allowed between clusters.\
                                 The maximum capacity (or overflow) needs to be sufficient to transport all flow exceeding the minimum demand at each cluster node to the virtual node.\
                                 If the overflow specified by the user is too small, it will be changed to match the minimum flow required. \
                                 Mandatory arguments include an input mascon file and an area in square meters defining the resolution. \
                                 Three different options are available to select the primary mascons to be reshaped : \
                                 (1) select all primaries inside the polygons of a polygon file (you can choose which polygon, default =\'all\'), \
                                 (2) select all primaries with one given region name, \
                                 (3) select one primary by number. You can choose an output filename, the default is your input filename + \'_reshaped\'. \
                                 The user can also change the overflow (maximum difference in size between primaries: 10 by default), the precision (maximum displacement of clusters defining convergence: 5000 m by default) \
                                 and the maxsize, defining the maximum number of arcs between ternaries and primaries that the algorithm can solve. \
                                 If the number of arcs exceed maxsize, the problem will de divided in several subregions. \
                                 The user can increase the maxsize, to avoid divinding the problem in several subregions. \
                                 Reference: Bradley PS, Bennett KP, Demiriz A (2000) Constrained K-Means Clustering. Microsoft Research. Available: http://research.microsoft.com/pubs/69796/tr-2000-65.pdf \
                                 Author : Julia Pfeffer, 17 February 2020 ')
parser.add_argument("in_masconfile",default="None",type=str,help="complete path to your input mascon file (string) ")
parser.add_argument("area",default=4e10,type=float,help="target area for reshaped mascons (m^2)")
parser.add_argument("--overflow",default=10,type=int,help="maximum difference in size (in terms of number of ternaries) between the primaries. the overflow has to be a positive integer. If will be adjusted to the minimum required to solve the problem, if the value indicated is too small.")
parser.add_argument("--precision",default=5000,type=int,help="maximum displacement (in meters) authorized for a cluster center. If all clusters move less than precion, your problem has converged and the program will stop.")
parser.add_argument("--maxsize",default=3e6,type=int,help="maximum number of arcs ((nb. of ternaries + 1) x (nb. primaries)) that the program will solve at once. If the dimension of the problem is greater than maxsize, the problem will first be divided in subregions, that will be solved one by one. ")
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
overflow=args.overflow
stop=args.precision
arclimit=args.maxsize
############################################################################### 
## READ INPUT FILES
###############################################################################
print('')
print('#######################################################################################################')
print('############################# START PYSHAPE_MASCONS:'+str(datetime.datetime.now())[0:19]+' ##############################')
print('#######################################################################################################')      
print('')
print('Read input files:'+str(datetime.datetime.now())[0:19])
iheaders,iPdata,_,iTdata=ut.read_mascon_file(input_mascon_file)# Read mascon file
inbt=len(iTdata)
############################################################################### 
## SELECT PRIMARY MASCONS TO RESHAPE
###############################################################################
print('Select primary mascons to reshape: '+str(datetime.datetime.now())[0:19])
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
                print('No primary mascon associated with number '+str(mypnumber)+' in '+input_mascon_file)
                sys.exit(1)
            else:
                print('Found '+str(int(nbp_in))+' primaries corresponding to primary number '+str(int(mypnumber))+':'+str(datetime.datetime.now())[0:19])
    else:
        inPdata=iPdata[iPdata[:,13]==myregion]
        nbp_in=len(inPdata)  
        outPdata=iPdata[iPdata[:,13]!=myregion]
        nbp_out=len(outPdata)
        if nbp_in==0:
            print('No primary mascon associated with region '+myregion+' in '+input_mascon_file)
            sys.exit(1)
        else:
            print('Found '+str(int(nbp_in))+' primaries corresponding to region '+myregion+' in '+input_mascon_file+':'+str(datetime.datetime.now())[0:19])
else:    
    polygons,coords,nbv=ut.read_polygons_ascii(polygon_file,mypolygon)
    nbpol=len(coords)
    print('Define bounding box around polygons: '+str(datetime.datetime.now())[0:19])
    minlat,maxlat,minlon,maxlon=ut.boundingbox(coords) 
    if minlon<0: # if bounding box crosses the Greenwitch meridian, convert positive longitudes to negatives longitudes
        iPdata[:,5][iPdata[:,5]>180]=iPdata[:,5][iPdata[:,5]>180]-360*np.ones(len(iPdata[:,5][iPdata[:,5]>180]))
        iTdata[:,3][iTdata[:,3]>180]=iTdata[:,3][iTdata[:,3]>180]-360*np.ones(len(iTdata[:,3][iTdata[:,3]>180]))
        coords=ut.convert_longitudes_in_tuple(coords,flag='to_negative')
    if mydensity ==0:
        Pmask=(iPdata[:,4]>=minlat)&(iPdata[:,4]<=maxlat)&(iPdata[:,5]<=maxlon)&(iPdata[:,5]>=minlon) # defime mask with bounding box
    else:
        Pmask=(iPdata[:,4]>=minlat)&(iPdata[:,4]<=maxlat)&(iPdata[:,5]<=maxlon)&(iPdata[:,5]>=minlon) &(iPdata[:,10]==mydensity)# defime mask with bounding box & density
    boxPdata=iPdata[Pmask==True] # ternaries inside bounding box with density = mydensity
    nbpbox=len(boxPdata[:,0]) # number of ternaries inside the box
    print('Found '+str(int(nbpbox))+' primaries in bounding box: '+str(datetime.datetime.now())[0:19])
    del minlat,maxlat,minlon,maxlon
    print('Find primaries inside polygons: '+str(datetime.datetime.now())[0:19])
    id_primary,id_polygon=ut.select_primaries_in_polygons(boxPdata,polygons,coords,print_output='no')# Select primaries in polygons
    inout = np.zeros(nbpbox)
    inout[id_primary] = 1 # Define mask where 1 = inside polygon and 0 = outside polygon            
    inPdata=boxPdata[inout==1] # Define array of primaries inside polygon
    nbp_in=len(inPdata[:,0])
    outPdata=np.concatenate((iPdata[Pmask==False],boxPdata[inout==0]))# Define array of primaries outside polygons 
    nbp_out=len(outPdata[:,0])
    del boxPdata,nbpbox,inout,Pmask,iPdata,polygons,coords,id_primary,id_polygon
    if nbp_in==0:
        print('No primary mascon associated with polygon '+mypolygon+' from '+polygon_file+' in '+input_mascon_file)
        sys.exit(1)
    else:
        print('Found '+str(int(nbp_in))+' primaries inside polygons: '+str(datetime.datetime.now())[0:19]) 
################################################################################
## RESHAPE TERNARIES IN SELECTED PRIMARIES  
################################################################################
# select ternaries to reshape
nbt=np.sum(inPdata[:,3])
inTdata=np.array(np.zeros((nbt,14)),dtype=object)
tcount=0
for pcpt in np.arange(nbp_in):
    ntip=inPdata[pcpt,3]
    inTdata[tcount:tcount+ntip]=iTdata[iTdata[:,9]==inPdata[pcpt,0],:]
    tcount=tcount+ntip
# Project ternaries coordinates with an azimutal equidistent projection (center = region center)
lon_proj,lat_proj=ut.get_center(inTdata[:,3], inTdata[:,2]) # get region center
wgs84=pyproj.Proj(init='epsg:4326') # define wgs84 (lat,lon datum)
pproj= pyproj.Proj(proj="aeqd", lat_0=lat_proj, lon_0=lon_proj, datum="WGS84", units="m") # define azimutal equidistqnt projection
TX,TY=pyproj.transform(wgs84,pproj,inTdata[:,3],inTdata[:,2]) # project lon,lat to X,Y
coords=np.hstack((TX.reshape((nbt,1)),TY.reshape((nbt,1)))) # stack coordinates 
# Define number of primaries with reshaped size
ternary_area=np.mean(inTdata[:,5]) # average area of a ternary inside my region
target_ntip=int(np.floor(myarea/ternary_area)) # number of ternaries inside an ideal primary (size=resolution)
nbp=int(np.floor(nbt/target_ntip)) # number of primaries inside region to fit the requited resolution
if nbp==0:
    nbp=1
if (nbt+1)*nbp>arclimit:
    # Define number of primaries and number of trenaries that should be inside primaries    
    nbsub=int(np.ceil(np.sqrt((nbp*(nbt+1))/arclimit)))
    print('#################################################################################')
    print(' Problem too large : (nb ternaries + 1) x number of primaries  > arc_limit ('+str(int((nbt+1)*nbp))+' > '+str(int(arclimit))+')')
    print(' Divide '+myregion+' in '+str(int(nbsub))+' subregions without constraints, to be reshaped at requited resolution ('+str(int(myarea))+' m2).')
    print('#################################################################################')
    subclusters,subcentroids=ut.get_unconstrained_clusters(coords,nbsub,stop*2)
    inTdata=inTdata[np.argsort(subcentroids)] 
    inTdata[:,9]=subcentroids[np.argsort(subcentroids)]
    inTdata[:,10]=subcentroids[np.argsort(subcentroids)]
    print('Will now reshape each subregions at requited size.')
    clusters=[]
    centroids=[]
    for pcpt in np.arange(nbsub): 
        Tselect=inTdata[inTdata[:,9]==pcpt,:]
        nbt_select=len(Tselect)
        select_area=np.mean(Tselect[:,5]) # average area of a ternary inside that region
        select_target=int(np.floor(myarea/select_area)) # number of ternary inside a primary of ideal size
        nbp_select=int(np.floor(nbt_select/select_target))
        print('######################################################################################')
        print('Reshape subregion '+str(int(pcpt+1))+' in '+str(int(nbp_select))+' primary mascons with at least '+str(int(select_target))+' and at most '+str(int(select_target+overflow))+' ternaries.')
        print('######################################################################################')
        min_overflow=np.ceil((nbt_select-(nbp_select*select_target))/float(nbp_select))
        if overflow<min_overflow:
            print('Overflow is too small. Problem not feasible.')
            print('Overflow set at '+str(int(min_overflow))+'(minimum value).')
            overflow=min_overflow
        lon_sproj,lat_sproj=ut.get_center(Tselect[:,3], Tselect[:,2])
        sproj = pyproj.Proj(proj="aeqd", lat_0=lat_sproj, lon_0=lon_sproj, datum="WGS84", units="m")
        TX_select,TY_select=pyproj.transform(wgs84,sproj,Tselect[:,3],Tselect[:,2])
        coords_select=np.hstack((TX_select.reshape((nbt_select,1)),TY_select.reshape((nbt_select,1))))
        # get clusters with constrained size
        clusters_select,centroids_select=ut.get_constrained_size_clusters(coords_select,nbp_select,select_target,select_target+overflow,stop)
        Plon_select,Plat_select=pyproj.transform(sproj,wgs84,clusters_select[:,0],clusters_select[:,1])
        ccoords_select=np.zeros((nbp_select,2))
        ccoords_select[:,0]=Plon_select
        ccoords_select[:,1]=Plat_select
        clusters.extend(ccoords_select)
        if len(centroids)>0:
            max_pid=1+max(centroids)
            centroids.extend(centroids_select+max_pid*np.ones(len(centroids_select)))
        else:
            centroids.extend(centroids_select)
        del clusters_select,centroids_select,ccoords_select,Tselect,nbt_select,select_area,select_target,nbp_select,lon_sproj,lat_sproj,sproj,TX_select,TY_select
           ## write reshaped data in array    
    cid=np.asarray([int(i) for i in centroids])
    inTdata=inTdata[np.argsort(cid)] 
    inTdata[:,9]=cid[np.argsort(cid)]
    inTdata[:,10]=cid[np.argsort(cid)]
    inPdata=np.array(np.zeros((np.max(cid)+1,14)),dtype=object)  
    clusters=np.vstack(clusters)
    tcount=-1    
    for pcpt in np.arange(np.max(cid)+1):
        ternaries=inTdata[inTdata[:,9]==pcpt,:]# select ternaries for each primary inside polygon
        ntip=len(ternaries[:,0]) # number of ternaries inside primary
        pcland=100*len(ternaries[:,8][ternaries[:,8]==1000])/ntip # percentage of land in primary
        inPdata[pcpt,0]=pcpt
        inPdata[pcpt,1]='P'+inTdata[inTdata[:,9]==pcpt,1][0][1:]
        inPdata[pcpt,2]=1
        inPdata[pcpt,3]=ntip
        inPdata[pcpt,4]=clusters[pcpt,1]
        inPdata[pcpt,5]=clusters[pcpt,0]
        inPdata[pcpt,6]=np.mean(inTdata[inTdata[:,9]==pcpt,4])
        inPdata[pcpt,7]=np.sum(inTdata[inTdata[:,9]==pcpt,5])
        inPdata[pcpt,8]=np.mean(inTdata[inTdata[:,9]==pcpt,6])
        inPdata[pcpt,9]=np.mean(inTdata[inTdata[:,9]==pcpt,7])
        inPdata[pcpt,10]=inTdata[inTdata[:,9]==pcpt,8][0]
        inPdata[pcpt,11]=pcland  
        inPdata[pcpt,12]=0
        inPdata[pcpt,13]=inTdata[inTdata[:,9]==pcpt,13][0]

        # RM200217: randomly assign colour to primary
        if inPdata[pcpt,10] < 1010:
            col_range=[750,1500]
            prim_colour=col_range[0] + randint(0,col_range[1]-col_range[0])
        elif inPdata[pcpt,10] > 1010:
            col_range=[0,450]
            prim_colour=col_range[0] + randint(0,col_range[1]-col_range[0])

        # RM200217: assign prim_col to each ternary in new primary
        for tcpt in np.arange(ntip):
            tcount=tcount+1
	    inTdata[tcount,11:13]=[prim_colour,prim_colour]
else:    
    print('######################################################################################')
    print('Reshape region in '+str(int(nbp))+' primary mascons with at least '+str(int(target_ntip))+' and at most '+str(int(target_ntip+overflow))+' ternaries.')
    print('######################################################################################')
    # get clusters with minimum size
    min_overflow=np.ceil((nbt-(nbp*target_ntip))/float(nbp))
    if overflow<min_overflow:
        print('Overflow is too small. Problem not feasible.')
        print('Overflow set at '+str(int(min_overflow))+'(minimum value).')
        overflow=min_overflow
    clusters,centroids=ut.get_constrained_size_clusters(coords,nbp,target_ntip,target_ntip+overflow,stop)
    ## write reshaped data in array    
    inTdata=inTdata[np.argsort(centroids)] 
    inTdata[:,9]=centroids[np.argsort(centroids)]
    inTdata[:,10]=centroids[np.argsort(centroids)]
    inPdata=np.array(np.zeros((nbp,14)),dtype=object)  
    inPlon,inPlat=pyproj.transform(pproj,wgs84,clusters[:,0],clusters[:,1])
    tcount=-1
    for pcpt in np.arange(nbp):
        ternaries=inTdata[inTdata[:,9]==pcpt,:]# select ternaries for each primary inside polygon
        ntip=len(ternaries[:,0]) # number of ternaries inside primary
        pcland=100*len(ternaries[:,8][ternaries[:,8]==1000])/ntip # percentage of land in primary
        inPdata[pcpt,0]=pcpt
        inPdata[pcpt,1]='P'+inTdata[inTdata[:,9]==pcpt,1][0][1:]
        inPdata[pcpt,2]=1
        inPdata[pcpt,3]=ntip
        inPdata[pcpt,4]=inPlat[pcpt]
        inPdata[pcpt,5]=inPlon[pcpt]
        inPdata[pcpt,6]=np.mean(inTdata[inTdata[:,9]==pcpt,4])
        inPdata[pcpt,7]=np.sum(inTdata[inTdata[:,9]==pcpt,5])
        inPdata[pcpt,8]=np.mean(inTdata[inTdata[:,9]==pcpt,6])
        inPdata[pcpt,9]=np.mean(inTdata[inTdata[:,9]==pcpt,7])
        inPdata[pcpt,10]=inTdata[inTdata[:,9]==pcpt,8][0]
        inPdata[pcpt,11]=pcland
        inPdata[pcpt,12]=0
        inPdata[pcpt,13]=inTdata[inTdata[:,9]==pcpt,13][0]

        # RM200217: randomly assign colour to primary
        if inPdata[pcpt,10] < 1010:
            col_range=[750,1500]
            prim_colour=col_range[0] + randint(0,col_range[1]-col_range[0])
        elif inPdata[pcpt,10] > 1010:
            col_range=[0,450]
            prim_colour=col_range[0] + randint(0,col_range[1]-col_range[0])

        # RM200217: assign prim_col to each ternary in new primary
        for tcpt in np.arange(ntip):
            tcount=tcount+1
	    inTdata[tcount,11:13]=[prim_colour,prim_colour]

###############################################################################
## CREATE TERNARY AND PRIMARY ARRAYS
################################################################################
print('')
print('Create full primary and ternary arrays including reshaped and unchanged mascons. '+str(datetime.datetime.now())[0:19])
######## Build empty primary and ternary arrays
new_nbp_in=len(inPdata[:,0])
nbt_out=inbt-nbt
new_Pdata=np.array(np.zeros((new_nbp_in+nbp_out,14)),dtype=object)
new_Tdata=np.array(np.zeros((inbt,14)),dtype=object)
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
new_Pdata[nbp_out:,:]=inPdata 
new_Pdata[nbp_out:,0]=new_Pdata[nbp_out:,0]+(nbp_out+1)*np.ones(new_nbp_in)
new_Tdata[nbt_out:,:]=inTdata
new_Tdata[nbt_out:,9]=new_Tdata[nbt_out:,9]+(nbp_out+1)*np.ones(nbt)
new_Tdata[nbt_out:,10]=new_Tdata[nbt_out:,10]+(nbp_out+1)*np.ones(nbt)
new_Pdata[:,5][new_Pdata[:,5]<0]=new_Pdata[:,5][new_Pdata[:,5]<0]+360*np.ones(len(new_Pdata[:,5][new_Pdata[:,5]<0]))
new_Tdata[:,3][new_Tdata[:,3]<0]=new_Tdata[:,3][new_Tdata[:,3]<0]+360*np.ones(len(new_Tdata[:,3][new_Tdata[:,3]<0]))  
##################################################################################
##### Define new headers
# first line
npex=nbp_out+new_nbp_in
ntex=inbt
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
new_headers=ut.define_new_mascon_headers(npex,ntex,maxtip,iheaders,newline)                   
del maxtip,iheaders, newline,npex,ntex,nbp_out,nbp_in,nbt
###################################################################################
# Write mascon file
ut.write_mascon_file(out_masconfile,new_headers,new_Pdata,new_Tdata)
del new_headers, new_Pdata,new_Tdata
################################################################################## 
print('')
print('#######################################################################################################')
print('############## SUGGESTION: run repair_mascons in ~/gt/util on '+out_masconfile+' #################')
print('#######################################################################################################')         
print('## This will generate a valid hashcode and recompute the latitude, longitude, radius, area, altitude ##')
print('######## and geoid height of each primary and secondary mascon to be consistent with ternaries ########')
print('#######################################################################################################')
print('')
print('#######################################################################################################')
print('############################# END PYSHAPE_MASCONS:'+str(datetime.datetime.now())[0:19]+'################################')
print('#######################################################################################################')      
#################################################################################    
  
