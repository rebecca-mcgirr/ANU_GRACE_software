#!/scratch/compute1/geodynamics/pyenv/py27/bin/python2.7
# -*- coding: utf-8 -*-
import numpy as np
import pyMCFSimplex as pyMCF
import pyproj
import shapely.geometry as geometry
import shapely.ops as ops
from functools import partial
import os
import datetime as dt
import sys
"""
Created on Wed Feb 12 14:46:02 2020

@author: julia
"""
###############################################################################
def get_center(lon, lat, R=1):
    """
    calculates lon, lat coordinates of a point on a sphere with
    radius R
    """
    lon_r = np.radians(lon.astype(float))
    lat_r = np.radians(lat.astype(float))

    x = R * np.cos(lat_r) * np.cos(lon_r)
    y = R * np.cos(lat_r) * np.sin(lon_r)
    z = R * np.sin(lat_r)
    
    mx=np.sum(x)/len(x)
    my=np.sum(y)/len(y)
    mz=np.sum(z)/len(z)
    
    mlon = 180/np.pi * np.arctan2(my, mx);
    mlat = 180/np.pi * np.arctan2(mz, np.sqrt(mx * mx + my * my));
    
    return mlon, mlat
###############################################################################
def get_costs(coords,clusters):
    m=len(coords[:,0])
    k=len(clusters[:,0])
    distances=np.zeros((m,k))
    for cpt in np.arange(k):
        distances[:,cpt]=np.sqrt(np.square(coords[:,0]-clusters[cpt,0]*np.ones(m))+np.square(coords[:,1]-clusters[cpt,1]*np.ones(m)))
    costs=np.zeros(m*k+k)   
    costs[0:m*k]=distances.ravel()
    return costs
###############################################################################
def get_balance(m,k,min_nbt):
    node_balance=np.zeros(m+k+1)
    node_balance[0:m]=np.ones(m)
    node_balance[m:m+k]=-min_nbt*np.ones(k)
    node_balance[-1]=-m+min_nbt*k
    return node_balance
###############################################################################
def get_arcs(m,k):
    arcs=np.zeros((m*k+k,2))
    count=0 
    for pcpt in np.arange(m):
        for kcpt in np.arange(k):
            arcs[count,0]=pcpt+1
            arcs[count,1]=m+kcpt+1
            count=count+1
    for kcpt in np.arange(k):
        arcs[count,0]=m+kcpt+1  
        arcs[count,1]=m+k+1 
        count=count+1
    return arcs    
###############################################################################
def get_flow(m,k):
    flow=np.zeros(((m*k+k,2)))
    flow[:m*k,1]=np.ones(m*k)
    flow[m*k:,1]=np.ones(k)*m
    return flow
###############################################################################
def get_min_max_flow(m,k,min_cluster_size,max_cluster_size):
    flow=np.zeros(((m*k+k,2)))
    flow[:m*k,1]=np.ones(m*k)
    flow[m*k:,1]=np.ones(k)*(max_cluster_size-min_cluster_size)
    return flow
###############################################################################
def write_dimacs(filename,nb_nodes,nb_arcs,balance,arcs,flow,costs):
    fw=open(filename,'w')
    fw.write('p min %d %d\n'%(nb_nodes,nb_arcs))  
    for cpt in np.arange(nb_nodes):
        fw.write('n %d %d\n'%(cpt+1,balance[cpt]))
    for cpt in np.arange(nb_arcs):
        fw.write('a %d %d 0 %d %d\n'%(arcs[cpt,0],arcs[cpt,1],flow[cpt,1],int(costs[cpt]*100)))    
    fw.close()   
###############################################################################
def get_constrained_size_clusters(coords,nb_clusters,min_cluster_size,max_cluster_size,stop):
    # read function arguments
    m=len(coords[:,0])
    k=nb_clusters
    min_nbt=min_cluster_size
    epsilon=stop+999
    # choos randomly k clusters among m ternaries
    init=np.random.choice(np.arange(m),k)
    clusters=coords[init,:]
    # design constant parameters to solve the minimum-cost-flow problem
    nb_nodes=m+k+1
    nb_arcs=m*k+k
    mybalance=get_balance(m,k,min_nbt) # supply (positive) and demand (negative) at each node
    myarcs=get_arcs(m,k) # authorized arcs between nodes
    myflow=get_min_max_flow(m,k,min_cluster_size,max_cluster_size) # authorized flow between nodes (minimum and maximum capacities)
    while epsilon>stop:
        ## set variable network parameters (cost is a function of distance between the ternaries and the cluster centers)
        mycosts=get_costs(coords,clusters)
        ### write network file
        dirpath=os.getcwd()
        filename=dirpath+'/network.mx'
        write_dimacs(filename,nb_nodes,nb_arcs,mybalance,myarcs,myflow,mycosts)
        ### build network
        mcf = pyMCF.MCFSimplex()
        f = open(filename,'r')
        mxdata = f.read()
        f.close()
        mcf.LoadDMX(mxdata)
        ## solve minimum cost flow problem
        mcf.SolveMCF()
        ## get arc results
        length = mcf.MCFm() # arc length
        flow = pyMCF.new_darray(length)
        length = mcf.MCFn() # node length
        nms = pyMCF.new_uiarray(length)
        mcf.MCFGetX(flow,nms)
        t2p=np.zeros((m,2)) # arcs between ternaries and primaries solving the minimum cost flow problem
        count=0
        i=-1
        while count<m:
            i=i+1
            if (pyMCF.darray_get(flow,i)==1)&(pyMCF.uiarray_get(nms,i)<=m*k+1):
                t2p[count,:]=myarcs[pyMCF.uiarray_get(nms,i),:]
                count=count+1
        # get indices for cluster assignement for each point in coords
        t2p[:,0]=  t2p[:,0]-np.ones(m) #numbering of arcs starts at 1, python idexing starts at 0
        t2p[:,1]=  t2p[:,1]-(m+1)*np.ones(m) # primaries number = node number - m (the first m nodes are ternaries)
        indices=(t2p[:,1]) # indicate wgich primary solves the minimum cost flow problem for each ternary
        # get new clusters coordinates
        new_clusters=np.zeros((k,2))
        for pcpt in np.arange(k):  
            new_clusters[pcpt,0]=np.mean(coords[indices==pcpt,0]) 
            new_clusters[pcpt,1]=np.mean(coords[indices==pcpt,1]) 
        # compute clusters displacement
        pdist=np.sqrt(np.square(new_clusters[:,0]-clusters[:,0])+np.square(new_clusters[:,1]-clusters[:,1])) 
        clusters=new_clusters
        # stop if all cluster displacements are less than stop
        epsilon=np.max(pdist)
        # display epsilon and time to estimate if the algorithm is going to converge
        print('Maximum cluster displacement = '+str(int(np.ceil(epsilon/1000.0)))+ ' km. Stop at '+str(int(stop/1000.0))+' km. '+str(dt.datetime.now())[0:19])
        os.remove(filename)
        del mcf,length,flow,nms,t2p
    return clusters,indices
###############################################################################
def get_unconstrained_clusters(coords,nb_clusters,stop):
    TX=coords[:,0]
    TY=coords[:,1]
    nbt=len(TX)
    nbp=nb_clusters
    # initialization
    init=np.random.choice(np.arange(nbt),nbp)
    PX=TX[init]
    PY=TY[init]
    epsilon=stop+999999
    while epsilon>stop:
        # assignement
        tdist=np.zeros((nbt,nbp))
        for pcpt in np.arange(nbp):
            tdist[:,pcpt]=np.sqrt(np.square(TX-PX[pcpt]*np.ones(nbt)) + np.square(TY-PY[pcpt]*np.ones(nbt)))
        centroids=np.argmin(tdist,axis=1)
        ## update
        new_PX=np.zeros(nbp)
        new_PY=np.zeros(nbp)
        pdist=np.zeros(nbp)
        for pcpt in np.arange(nbp):
            new_PX[pcpt]=np.mean(TX[centroids==pcpt]) 
            new_PY[pcpt]=np.mean(TY[centroids==pcpt])    
            pdist[pcpt]=np.sqrt(np.square(new_PX[pcpt]-PX[pcpt])+np.square(new_PY[pcpt]-PY[pcpt]))
        PX=new_PX
        PY=new_PY
        epsilon=np.max(pdist)
        print('Maximum cluster displacement = '+str(int(np.ceil(epsilon/1000.0)))+ ' km. Stop at '+str(int(stop/1000.0))+' km. '+str(dt.datetime.now())[0:19])    
    clusters=np.zeros((nbp,2)) 
    clusters[:,0]=PX
    clusters[:,0]=PY   
    return clusters, centroids
###############################################################################
def read_mascon_file(fname):
    ifile=open(fname,'r')
    ilines=ifile.readlines()
    ifile.close()
    del ifile
    # Get number of primary,secondary and ternary mascons
    nbp=int(ilines[0].split()[1])
    nbs=int(ilines[0].split()[2])
    nbt=int(ilines[0].split()[3])
    #
    headers=[]
    Pdata=np.array(np.zeros((nbp,14)),dtype=object)
    Sdata=np.array(np.zeros((nbs,13)),dtype=object)
    Tdata=np.array(np.zeros((nbt,14)),dtype=object)
    #
    pcount=-1
    scount=-1
    tcount=-1
    for lcpt in np.arange(len(ilines)):
        if ilines[lcpt][0]=='#': #  header lines (comments starting with #)
            headers.append(ilines[lcpt])
        else:
            msc_data = ilines[lcpt].split()
            if msc_data[1][0] == 'P': # primaries
                pcount=pcount+1
                Pdata[pcount,0]=int(msc_data[0]) # primary mascon number
                Pdata[pcount,1]=msc_data[1] # mascon type (land/shelf/ocean)
                Pdata[pcount,2]=int(msc_data[2]) # number of ternaries in primary
                Pdata[pcount,3]=int(msc_data[3]) # primary mascon centre latitude
                Pdata[pcount,4]=float(msc_data[4])
                Pdata[pcount,5]=float(msc_data[5])
                Pdata[pcount,6]=float(msc_data[6])
                Pdata[pcount,7]=float(msc_data[7])
                Pdata[pcount,8]=float(msc_data[8])
                Pdata[pcount,9]=float(msc_data[9])
                Pdata[pcount,10]=float(msc_data[10])
                Pdata[pcount,11]=float(msc_data[11])
                Pdata[pcount,12]=int(msc_data[12])
                Pdata[pcount,13]=msc_data[13]
            elif msc_data[1][0] == 'S':  # secondaries
                scount=scount+1
                Sdata[scount,0]=int(msc_data[0])
                Sdata[scount,1]=msc_data[1]
                Sdata[scount,2]=int(msc_data[2])
                Sdata[scount,3]=float(msc_data[3])
                Sdata[scount,4]=float(msc_data[4])
                Sdata[scount,5]=float(msc_data[5])
                Sdata[scount,6]=float(msc_data[6])
                Sdata[scount,7]=float(msc_data[7])
                Sdata[scount,8]=float(msc_data[8])
                Sdata[scount,9]=float(msc_data[9])
                Sdata[scount,10]=int(msc_data[10])
                Sdata[scount,11]=int(msc_data[11])
                Sdata[scount,12]=msc_data[12]             
            elif msc_data[1][0] == 'T':  # ternaries
                tcount=tcount+1
                Tdata[tcount,0]=int(msc_data[0])
                Tdata[tcount,1]=msc_data[1]
                Tdata[tcount,2]=float(msc_data[2])
                Tdata[tcount,3]=float(msc_data[3])
                Tdata[tcount,4]=float(msc_data[4])
                Tdata[tcount,5]=float(msc_data[5])
                Tdata[tcount,6]=float(msc_data[6])
                Tdata[tcount,7]=float(msc_data[7])
                Tdata[tcount,8]=float(msc_data[8])
                Tdata[tcount,9]=int(msc_data[9])
                Tdata[tcount,10]=int(msc_data[10])
                Tdata[tcount,11]=int(msc_data[11])
                Tdata[tcount,12]=int(msc_data[12])
                Tdata[tcount,13]=msc_data[13]
    return headers, Pdata,Sdata,Tdata
###############################################################################
def read_polygons_ascii(fname,fregion='all'):
    f=open(fname)
    lines=f.readlines()
    f.close()
    list_vertices=[]
    polygons=[]
    startpol=[]
    coords=[]
    for cpt in np.arange(len(lines)):
        if len(lines[cpt].split())==3:
            startpol.append(cpt)
            nb_vertices=int(lines[cpt].split()[0])
            list_vertices.append(nb_vertices)
            polygons.append(lines[cpt].split()[1])
            coords_temp=[]
            for cpt2 in np.arange(cpt+1,cpt+1+nb_vertices):
                vlon=float(lines[cpt2].split()[1])
                if vlon<0:
                    vlon=vlon+360
                vlat=float(lines[cpt2].split()[0])
                coords_temp.append((vlon,vlat))
            coords.append(coords_temp)
            del coords_temp  
    if fregion=="all":
        print('The program will run for %d regions defined in %s.'%(len(polygons),fname))
    else:
        rindex=-1 # Define index to identify myregion in regions
        for polcpt in np.arange(len(polygons)):
            if polygons[polcpt]==fregion:
                rindex=polcpt
                break
        if rindex==-1:
            print('%s was not found in %s.'%(fregion,fname))
            sys.exit()
        else:
            polygons=[polygons[rindex]]
            coords=[coords[rindex]]
            list_vertices=[list_vertices[rindex]]
            del rindex
            print('The program will run for %d region named %s.'%(len(polygons),polygons[:]))         
    return polygons,coords,list_vertices  
###############################################################################
def boundingbox(coords,bufferdist=0):
    lat=[]
    lon=[]
    for pcpt in np.arange(len(coords)): # loop on polygons 
        for vcpt in np.arange(len(coords[pcpt])):# loop on vertices 
            lat.append(coords[pcpt][vcpt][1])
            lon.append(coords[pcpt][vcpt][0])
    lon=np.asarray(lon)
    lat=np.asarray(lat)
    if np.max(lon > 180) and np.min(lon) < 180:    
        lon[lon>180]=lon[lon>180]-360*np.ones(len(lon[lon>180]))      
    minlat=np.min(lat)
    maxlat=np.max(lat)
    minlon=np.min(lon)
    maxlon=np.max(lon) 
    return minlat,maxlat,minlon,maxlon
###############################################################################
def convert_longitudes_in_tuple(coords,flag='to_negative'):
    for polcpt in np.arange(len(coords)):
        for vcpt in np.arange(len(coords[polcpt])):
            clist=list(coords[polcpt][vcpt])
            if flag=='to_negative':
                if clist[0]>180:
                    clist[0]=clist[0]-360
                    coords[polcpt][vcpt]=tuple(clist)
            if flag=='to_positive':
                if clist[0]<0:
                    clist[0]=clist[0]+360
                    coords[polcpt][vcpt]=tuple(clist)               
    return coords
###############################################################################
def select_primaries_in_polygons(Pdata,polname,polcoords,print_output='no'):
    nbp=len(Pdata)
    nbpol=len(polname)
    wgs84=pyproj.Proj(init='epsg:4326') 
    ## Find primaries in each polygon
    id_primary=[]
    id_polygon=[]
    for polcpt in np.arange(nbpol):# check if point is contained in one of the polygons
        mypoly=geometry.Polygon(polcoords[polcpt])
        polproj = pyproj.Proj(proj="aeqd", lat_0=mypoly.centroid.coords[0][1], lon_0=mypoly.centroid.coords[0][0], datum="WGS84", units="m")
        project = partial(pyproj.transform,polproj,wgs84)
        mypoly=ops.transform(project,mypoly)# project polygon in azimuthal equidistant projection
        for pcpt in np.arange(nbp):
             mypoint= ops.transform(project,geometry.Point(Pdata[pcpt,5],Pdata[pcpt,4]))
             if mypoly.contains(mypoint):
                 if print_output=='yes':
                     print('Primary number = %d, area = %f, region = '%(Pdata[pcpt,0],Pdata[pcpt,7])+polname[polcpt])
                 id_primary.append(pcpt) 
                 id_polygon.append(polcpt) 
    return id_primary,id_polygon 
###############################################################################    
###############################################################################    
def write_mascon_file(masconfile,headers,Pdata,Tdata):    
    print('Write headers in '+masconfile+' : '+str(dt.datetime.now())[0:19])
    ofile=open(masconfile, 'w')
    for hcpt in np.arange(len(headers)):
        ofile.write(headers[hcpt])
    del headers    
    # Write mascon data    
    print('Write primary, secondary and ternary data in '+masconfile+' : '+str(dt.datetime.now())[0:19])
    tcount=0
    for pcpt in np.arange(len(Pdata[:,0])):
        ofile.write('{:7d}  {:6s}{:8d}{:8d}{:10.4f}{:10.4f}{:11.1f}{:17.0f}{:1s}{:9.1f}{:9.1f}{:5.0f}{:1s}{:6.1f}{:6d}  {:14s}\n'.format(int(Pdata[pcpt,0]),str(Pdata[pcpt,1]),int(Pdata[pcpt,2]),int(Pdata[pcpt,3]),float(Pdata[pcpt,4]),float(Pdata[pcpt,5]),float(Pdata[pcpt,6]),float(Pdata[pcpt,7]),'.',float(Pdata[pcpt,8]),float(Pdata[pcpt,9]),float(Pdata[pcpt,10]),'.',float(Pdata[pcpt,11]),int(Pdata[pcpt,12]),str(Pdata[pcpt,13]).rjust(14)))
        ofile.write('{:7d}  {:6s}{:8d}        {:10.4f}{:10.4f}{:11.1f}{:17.0f}{:1s}{:9.1f}{:9.1f}{:5.0f}{:1s}{:6d}{:6d}  {:14s}\n'.format(int(Pdata[pcpt,0]),'S'+str(Pdata[pcpt,1][1:]),int(Pdata[pcpt,3]),float(Pdata[pcpt,4]),float(Pdata[pcpt,5]),float(Pdata[pcpt,6]),float(Pdata[pcpt,7]),'.',float(Pdata[pcpt,8]),float(Pdata[pcpt,9]),float(Pdata[pcpt,10]),'.',int(Pdata[pcpt,0]),int(Tdata[tcount,12]),str(Pdata[pcpt,13]).rjust(14)))
        for tcpt in np.arange(Pdata[pcpt,3]):
            ofile.write('{:7d}  {:6s}{:10.4f}{:10.4f}{:11.1f}{:14.0f}{:1s}{:9.1f}{:9.1f}{:5.0f}{:1s}{:6d}{:6d}{:6d}{:6d}         {:6s}\n'.format(int(Tdata[tcount,0]),str(Tdata[tcount,1]),float(Tdata[tcount,2]),float(Tdata[tcount,3]),float(Tdata[tcount,4]),float(Tdata[tcount,5]),'.',float(Tdata[tcount,6]),float(Tdata[tcount,7]),float(Tdata[tcount,8]),'.',int(Tdata[tcount,9]),int(Tdata[tcount,10]),int(Tdata[tcount,11]),int(Tdata[tcount,12]),str(Tdata[tcount,13]).rjust(6)))
            tcount=tcount+1
    ofile.close() 
###############################################################################    
def define_new_mascon_headers(npex,ntex,maxtip,file_history,newline):
    hashcode="#HASHCOD"
    nsex=npex
    maxsip=int(1)
    maxtis=maxtip
    ##### Define new headers
    headers=[]
    ## First line
    headers.append('{:8s}{:11d}{:8d}{:10d}{:7d}{:7d}{:7d}\n'.format(hashcode,npex,nsex,ntex,maxsip,maxtip,maxtis))
    ## Former header lines
    for hcpt in np.arange(1,len(file_history)):
        headers.append(file_history[hcpt])
    # Write new header line indicating changes made    
    headers.append(newline)   
    return headers  
###############################################################################
def map_mascons_and_polygons(masconfile,polygonfile,mypolygon='all',figname='map_mascons.png'):
    import matplotlib
    havedisplay = "DISPLAY" in os.environ
    if not havedisplay:
        matplotlib.use("Agg")
    from mpl_toolkits.basemap import Basemap
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Polygon
    import matplotlib.pyplot as plt
    import math
    # read final mascon file
    _,Pdata,_,Tdata=read_mascon_file(masconfile)
    regions,coords,nbv=read_polygons_ascii(polygonfile,mypolygon)
    mybufferdist=0
    minlat,maxlat,minlon,maxlon=boundingbox(coords,mybufferdist)
    if minlon<0: # if bounding box crosses the Greenwitch meridian, convert positive longitudes to negatives longitudes
        Pdata[:,5][Pdata[:,5]>180]=Pdata[:,5][Pdata[:,5]>180]-360*np.ones(len(Pdata[:,5][Pdata[:,5]>180]))
        Tdata[:,3][Tdata[:,3]>180]=Tdata[:,3][Tdata[:,3]>180]-360*np.ones(len(Tdata[:,3][Tdata[:,3]>180]))
        coords=convert_longitudes_in_tuple(coords,flag='to_negative') 
    del mybufferdist
    ## select mascons in polygon file
    Pmask=(Pdata[:,4]>=minlat)&(Pdata[:,4]<=maxlat)&(Pdata[:,5]<=maxlon)&(Pdata[:,5]>=minlon)# defime mask with bounding box & density
    boxPdata=Pdata[Pmask==True]
    id_primary,id_polygon=select_primaries_in_polygons(boxPdata,regions,coords,print_output='no')
    id_primary=np.asarray(id_primary)
    inPdata=boxPdata[id_primary]
    nbpin=len(id_primary)
    # select ternaries associated with primaries
    ternaries=[]
    for pcpt in np.arange(nbpin):
        ternaries.append(Tdata[Tdata[:,9]==inPdata[pcpt,0]])
    inTdata=np.vstack(np.asarray(ternaries))
    #### colorbar
    colors=[]
    cm20=plt.cm.Set1
    for cpt in np.arange(math.ceil(nbpin)/cm20.N):
        for cpt2 in np.arange(cm20.N):
            colors.append(cm20(cpt2))
    cm=matplotlib.colors.ListedColormap(colors)     
    ### figure 
    fig = plt.figure(figsize=(15,15))
    plt.cla()
    plt.clf()
    ax=fig.gca();
    m = Basemap(llcrnrlon=minlon,llcrnrlat=minlat, urcrnrlon=maxlon, urcrnrlat=maxlat, resolution = 'l')
    Tx, Ty = m(inTdata[:,3],inTdata[:,2])
    Px, Py = m(inPdata[:,5],inPdata[:,4])
    patches = []
    for pcpt in np.arange(len(coords)):
        patches.append(Polygon(coords[pcpt]))
    ax.add_collection(PatchCollection(patches, facecolor='none', edgecolor='k', linewidths=0.75,zorder=1))
    m.drawcoastlines(zorder=1)
    plt.scatter(Tx,Ty,5,inTdata[:,9],cmap=cm,vmin=np.min(inTdata[:,9]),vmax=np.max(inTdata[:,9]),edgecolors='face',zorder=2)
    for cpt in np.arange(len(inPdata)):
        t=plt.text(Px[cpt]-0.75,Py[cpt],inPdata[cpt,0],color='k',fontweight='bold',fontsize=10)
        t.set_bbox(dict(facecolor='white', alpha=0.6, edgecolor='none'))
    plt.savefig(figname)        
    del minlat,maxlat,minlon,maxlon, coords, regions,nbv
    del Pdata,Tdata,inPdata,inTdata,colors,cm20
###############################################################################
def map_mascons_in_region(masconfile,myregion,figname='map_mascons.png'):
    import matplotlib
    havedisplay = "DISPLAY" in os.environ
    if not havedisplay:
        matplotlib.use("Agg")
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt
    import math
    # read final mascon file
    _,Pdata,_,Tdata=read_mascon_file(masconfile)
    inPdata=Pdata[Pdata[:,13]==myregion]
    nbpin=len(inPdata)
    inTdata=Tdata[Tdata[:,13]==myregion]
    minlon=np.min(inTdata[:,3])
    maxlon=np.max(inTdata[:,3])
    minlat=np.min(inTdata[:,2])
    maxlat=np.max(inTdata[:,2])
    #### colorbar
    colors=[]
    cm20=plt.cm.Set1
    for cpt in np.arange(math.ceil(nbpin)/cm20.N):
        for cpt2 in np.arange(cm20.N):
            colors.append(cm20(cpt2))
    cm=matplotlib.colors.ListedColormap(colors)     
    ### figure 
    plt.figure(figsize=(15,15))
    plt.cla()
    plt.clf()
    m = Basemap(llcrnrlon=minlon,llcrnrlat=minlat, urcrnrlon=maxlon, urcrnrlat=maxlat, resolution = 'l')
    Tx, Ty = m(inTdata[:,3],inTdata[:,2])
    Px, Py = m(inPdata[:,5],inPdata[:,4])
    m.drawcoastlines(zorder=1)
    plt.scatter(Tx,Ty,5,inTdata[:,9],cmap=cm,vmin=np.min(inTdata[:,9]),vmax=np.max(inTdata[:,9]),edgecolors='face',zorder=2)
    for cpt in np.arange(nbpin):
        t=plt.text(Px[cpt]-0.75,Py[cpt],inPdata[cpt,0],color='k',fontweight='bold',fontsize=10)
        t.set_bbox(dict(facecolor='white', alpha=0.6, edgecolor='none'))
    plt.savefig(figname)        
    del minlat,maxlat,minlon,maxlon
    del Pdata,Tdata,inPdata,inTdata,nbpin,colors,cm20           
###############################################################################
def map_mascons_per_number(masconfile,ip_number,figname='map_mascons.png'):
    import matplotlib, os
    havedisplay = "DISPLAY" in os.environ
    if not havedisplay:
        matplotlib.use("Agg")
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt
    import math
    # read final mascon file
    _,Pdata,_,Tdata=read_mascon_file(masconfile)
    inPdata=[]
    for pcpt in np.arange(len(ip_number)):
        inPdata.append(Pdata[Pdata[:,0]==ip_number[pcpt]])
    inPdata=np.vstack(inPdata)
    nbpin=len(inPdata)
    nbtin=np.sum(inPdata[:,3])
    tcount=0
    inTdata=np.array(np.zeros((nbtin,14)),dtype=object)
    for pcpt in np.arange(nbpin):
        ntip=inPdata[pcpt,3]
        inTdata[tcount:tcount+ntip]=Tdata[Tdata[:,9]==ip_number[pcpt]]
        tcount=tcount+ntip
    minlon=np.min(inTdata[:,3])
    maxlon=np.max(inTdata[:,3])
    minlat=np.min(inTdata[:,2])
    maxlat=np.max(inTdata[:,2])
    #### colorbar
    colors=[]
    cm20=plt.cm.Set1
    for cpt in np.arange(math.ceil(nbpin)/cm20.N):
        for cpt2 in np.arange(cm20.N):
            colors.append(cm20(cpt2))
    cm=matplotlib.colors.ListedColormap(colors)     
    ### figure 
    plt.figure(figsize=(15,15))
    plt.cla()
    plt.clf()
    m = Basemap(llcrnrlon=minlon,llcrnrlat=minlat, urcrnrlon=maxlon, urcrnrlat=maxlat, resolution = 'l')
    Tx, Ty = m(inTdata[:,3],inTdata[:,2])
    Px, Py = m(inPdata[:,5],inPdata[:,4])
    m.drawcoastlines(zorder=1)
    plt.scatter(Tx,Ty,5,inTdata[:,9],cmap=cm,vmin=np.min(inTdata[:,9]),vmax=np.max(inTdata[:,9]),edgecolors='face',zorder=2)
    for cpt in np.arange(nbpin):
        t=plt.text(Px[cpt]-0.75,Py[cpt],inPdata[cpt,0],color='k',fontweight='bold',fontsize=10)
        t.set_bbox(dict(facecolor='white', alpha=0.6, edgecolor='none'))
    plt.savefig(figname)        
    del minlat,maxlat,minlon,maxlon
    del Pdata,Tdata,inPdata,inTdata,nbpin,colors,cm20
