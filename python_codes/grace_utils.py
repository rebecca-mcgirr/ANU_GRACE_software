#!/usr/bin/env python3
import datetime
from scipy.interpolate import griddata
from dateutil.relativedelta import relativedelta
import numpy as np
import pyproj
import shapely.geometry as geometry
import shapely.ops as ops
from functools import partial
import os, sys
import time
import re
#from geopy.distance import geodesic
from math import ceil
#import statsmodels.api as sm
###############################################################################
## DEFINE FUNCTIONS TO READ INPUT FILES
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
###############################################################################    
def write_mascon_file(masconfile,headers,Pdata,Tdata):    
    print('Write headers in %s:'%(masconfile),str(datetime.datetime.now())[0:19])
    ofile=open(masconfile, 'w')
    for hcpt in np.arange(len(headers)):
        ofile.write(headers[hcpt])
    del headers    
    # Write mascon data    
    print('Write primary, secondary and ternary data in %s:'%(masconfile),str(datetime.datetime.now())[0:19])
    tcount=0
    for pcpt in np.arange(len(Pdata[:,0])):
        ofile.write('{:7d}  {:6s}{:8d}{:8d}{:10.4f}{:10.4f}{:11.1f}{:17.0f}{:1s}{:9.1f}{:9.1f}{:5.0f}{:1s}{:6.1f}{:6d}  {:14s}\n'.format(Pdata[pcpt,0],Pdata[pcpt,1],Pdata[pcpt,2],Pdata[pcpt,3],Pdata[pcpt,4],Pdata[pcpt,5],Pdata[pcpt,6],Pdata[pcpt,7],'.',Pdata[pcpt,8],Pdata[pcpt,9],Pdata[pcpt,10],'.',Pdata[pcpt,11],Pdata[pcpt,12],str(Pdata[pcpt,13]).rjust(14)))
        ofile.write('{:7d}  {:6s}{:8d}        {:10.4f}{:10.4f}{:11.1f}{:17.0f}{:1s}{:9.1f}{:9.1f}{:5.0f}{:1s}{:6d}{:6d}  {:14s}\n'.format(Pdata[pcpt,0],'S'+Pdata[pcpt,1][1:],Pdata[pcpt,3],Pdata[pcpt,4],Pdata[pcpt,5],Pdata[pcpt,6],Pdata[pcpt,7],'.',Pdata[pcpt,8],Pdata[pcpt,9],Pdata[pcpt,10],'.',Pdata[pcpt,0],Tdata[tcount,12],str(Pdata[pcpt,13]).rjust(14)))
        for tcpt in np.arange(Pdata[pcpt,3]):
            ofile.write('{:7d}  {:6s}{:10.4f}{:10.4f}{:11.1f}{:14.0f}{:1s}{:9.1f}{:9.1f}{:5.0f}{:1s}{:6d}{:6d}{:6d}{:6d}         {:6s}\n'.format(Tdata[tcount,0],Tdata[tcount,1],Tdata[tcount,2],Tdata[tcount,3],Tdata[tcount,4],Tdata[tcount,5],'.',Tdata[tcount,6],Tdata[tcount,7],Tdata[tcount,8],'.',Tdata[tcount,9],Tdata[tcount,10],Tdata[tcount,11],Tdata[tcount,12],str(Tdata[tcount,13]).rjust(6)))
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
def project_coordinates_in_tuple(coords,inProj,outProj):
    for polcpt in np.arange(len(coords)):
        for vcpt in np.arange(len(coords[polcpt])):
            xin=coords[polcpt][vcpt][0]
            yin=coords[polcpt][vcpt][1]
            xout,yout = pyproj.transform(inProj,outProj,xin,yin)
            coords[polcpt][vcpt]=(xout,yout)                        
    return coords
###############################################################################
def toYearFraction_bad(date):
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = datetime.datetime(year=year, month=1, day=1)
    startOfNextYear = datetime.datetime(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction
#############################################################################
def toYearFraction(date):
    return date.year + int(date.strftime("%j"))/int(datetime.date(date.year,12,31).strftime("%j"))
#############################################################################
def decyear2nbdays(decvec):
    NT=len(decvec)
    daytime=np.zeros(NT)
    for cpt in np.arange(NT):
        year = int(decvec[cpt])
        remain = decvec[cpt] - year
        nbdays=(datetime.datetime(year+1,1,1)-datetime.datetime(year,1,1)).days
        nbyears=(datetime.datetime(year,1,1)-datetime.datetime(2000,1,1)).days
        daytime[cpt]=nbdays*remain+nbyears
    return daytime
##############################################################################
def read_fitfile(filename):
    f=open(filename)
    lines=f.readlines()
    f.close()
    if len(lines)<1:
        print('Empty file: %s'%(filename))
        return        
    mydate=re.search('\d{4}-\d{2}-\d{2}',lines[2]).group()
    # prefit statistics
    prestats=np.zeros(16)
    for cpt in np.arange(9,25):
        prestats[cpt-9]=str2float(lines[cpt].split()[-1])
    #postfit statistics
    poststats=np.zeros(16)
    for cpt in np.arange(27,43):
        poststats[cpt-27]=str2float(lines[cpt].split()[-1])
    # parameters for sat1 
    params_sat1=np.zeros((12,5))
    for cpt in np.arange(46,58):
        params_sat1[cpt-46,0]=str2float(lines[cpt][25:47])
        params_sat1[cpt-46,1]=str2float(lines[cpt][47:59])
        params_sat1[cpt-46,2]=str2float(lines[cpt][59:78])
        params_sat1[cpt-46,3]=str2float(lines[cpt][78:90])
        params_sat1[cpt-46,4]=str2float(lines[cpt][90:])
    # parameters for sat1 : colums= a priori, adjust,postfit,sigma and frac. Lines = position (x,y,z), velocity (x,y,z), scale (x,y,z), bias (x,y,z)
    params_sat2=np.zeros((12,5))
    for cpt in np.arange(59,71):
        params_sat2[cpt-59,0]=str2float(lines[cpt][25:47])
        params_sat2[cpt-59,1]=str2float(lines[cpt][47:59])
        params_sat2[cpt-59,2]=str2float(lines[cpt][59:78])
        params_sat2[cpt-59,3]=str2float(lines[cpt][78:90])
        params_sat2[cpt-59,4]=str2float(lines[cpt][90:])
    # parameters for mascons = = a priori, adjust,postfit,sigma and frac for all mascons
    if lines[-1]=='\n':
        nb_msc=int(lines[-2].sfplit()[1][2:])
    else:
        nb_msc=int(lines[-1].split()[1][2:])   
    params_msc=np.zeros((nb_msc,5))
    pcount=0
    for cpt in np.arange(72,len(lines)):
        if lines[cpt]=='\n':
            pass
        else:
            params_msc[pcount,0]=str2float(lines[cpt][25:47])
            params_msc[pcount,1]=str2float(lines[cpt][47:59])
            params_msc[pcount,2]=str2float(lines[cpt][59:78])
            params_msc[pcount,3]=str2float(lines[cpt][78:90])
            params_msc[pcount,4]=str2float(lines[cpt][90:])
            pcount=pcount+1
    return mydate,prestats,poststats,params_sat1,params_sat2,params_msc
###############################################################################
def str2float(mystring):
    try:
        float(mystring)
    except ValueError:
        value=99999
    else:
        value=float(mystring)
    return value
###############################################################################
def read_addnorm_fit(filename,nb_msc):
    f=open(filename)
    lines=f.readlines()
    f.close()
    if len(lines)<1:
        print('Empty file: %s'%(filename))
        return   
    nbf=int(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", lines[1])[1])

    #### PT201205: these three lines work for files like: addnorm_2003_10_02*.norm
    #mydate=re.search('\d{4}_\d{2}_\d{2}',lines[4][-45:]).group()
    #startdate=toYearFraction(datetime.date(int(mydate[0:4]),int(mydate[5:7]),1))
    #enddate=toYearFraction(datetime.date(int(mydate[0:4]),int(mydate[5:7]),1)+relativedelta(months=1))

    # PT201205: works for files like: msc_20031025*.norm
    #mydate=re.search('\d{4}\d{2}\d{2}',lines[3][-45:]).group()
    # PT220117: works for Seb's b8i3 files: msc_2008-11_b8i3.fit
    #mydate=re.search('\d{4}\d{2}\d{2}',lines[3][-16:]).group()

    # PT220128: change this to read the yyyy mm dd from the start time rather than the file name:
    # File:   1  Start:  2021   12    1    0    0    0 Duration:   23.96806 ../2021/2021-12-01/msc_2021_12_01_iter3_i1_00_12hr.norm
    mydate=lines[3][19:33]
    startdate=toYearFraction(datetime.date(int(mydate[0:4]),int(mydate[7:10]),1))

    mydate=lines[3+nbf-1][19:33]
    enddate  =toYearFraction(datetime.date(int(mydate[0:4]),int(mydate[7:10]),1))

    
    #print("startdate",startdate,mydate,nb_msc)

    monthly_date=startdate+(enddate-startdate)/2
    ####
    params_msc=99999*np.ones((nb_msc,4))
    pcount=0
    # PT220620: the "+4" seems to need to be "+3". Don't even know what it represents!!
    lcount=24*nbf+3
    while pcount<nb_msc:   
        #print(lcount,np.shape(lines))
        #print("XX",lines[lcount][0:24],"YY",lcount,pcount,lines[lcount][6:9])
        if lines[lcount]=='\n' or lines[lcount][6:9] != " MC" or lines[lcount][0:24]=="MASCON CONSTRAINT SIGMAS":
            #print("blah blah",lcount,pcount,nb_msc)
            lcount=lcount+1
            pass
        else:
            #print("mascon line",lcount,lines[lcount][0:23],pcount+1)
            params_msc[pcount,0]=float(lines[lcount].split()[-4])
            params_msc[pcount,1]=float(lines[lcount].split()[-3])
            params_msc[pcount,2]=float(lines[lcount].split()[-2])
            params_msc[pcount,3]=float(lines[lcount].split()[-1])        
            lcount=lcount+1
            pcount=pcount+1
            
       # return mydate,prestats,poststats,params_sat1,params_sat2,params_msc
    sigmas=np.zeros(nb_msc)
    if lcount<len(lines):
        #print(lcount,lines[lcount])
        if lines[lcount]=='MASCON CONSTRAINT SIGMAS\n': 
            pcount=0
            while pcount<nb_msc:
                #print(len(lines),nb_msc,lcount+1,lines[lcount+1])
                sigmas[pcount]=float(lines[lcount+1].split()[1])
                lcount=lcount+1
                pcount=pcount+1
    return monthly_date,params_msc,nbf,sigmas
############################################################################## 
def from_primary_to_ternary(Pdata,Tdata,Pvalue):
     nbtern=len(Tdata[:,0])
     nbp=len(Pdata[:,0])
     Tlat=Tdata[:,2]
     Tlon=Tdata[:,3]
     Tvalue=np.zeros((nbtern))
     tcount=0
     for pcpt in np.arange(nbp):
         ntip=Pdata[pcpt,3]
         Tvalue[tcount:tcount+ntip]=Pvalue[pcpt]*np.ones(ntip)
         tcount=tcount+ntip
     return Tlat,Tlon,Tvalue
############################################################################## 
def scattered_points_to_grid(Tlat,Tlon,Tvalue,lati,loni):
    lon2d,lat2d = np.meshgrid(loni,lati)
    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    x, y = ecef(Tlon,Tlat)
    xi, yi = ecef(lon2d,lat2d)
    grid_value= griddata((x,y),Tvalue,(xi,yi),method='linear')
    return grid_value
###############################################################################
def find_closest_mascon(Plat,Plon,latvec,lonvec):
    distancekm=99999*np.ones(len(latvec))
    for cpt in np.arange(len(latvec)):
        coord1=(Plat,Plon) # coordinates of the selected point
        coord2=(latvec[cpt],lonvec[cpt]) #coordinates of one subgrid element
        distancekm[cpt]=geodesic(coord1,coord2, ellipsoid='WGS-84').km 
    index=np.where(distancekm==np.min(distancekm))[0][0]  
    return index
###############################################################################
def get_annual_climatology(timevec,data):
    nbp=len(data[0,:])
    nbt=len(timevec)
    nby=ceil(timevec[-1]-timevec[0])+1
    temp=99999*np.ma.ones((nbp,nby,12))
    coefs=np.polynomial.polynomial.polyfit(timevec,data,2,rcond=1e-9)
    ffit = np.polynomial.polynomial.polyval(timevec, coefs) 
    ffit=np.transpose(ffit)
    for cpt in np.arange(nbt):
        rem = timevec[cpt] - int(timevec[cpt])
        delta=datetime.datetime(int(timevec[cpt])+1, 1, 1)-datetime.datetime(int(timevec[cpt]), 1, 1)
        my_date = datetime.datetime(int(timevec[cpt]), 1, 1)+datetime.timedelta(days=rem*delta.days)
        temp[:,int(my_date.year-2002),int(my_date.month-1)]=data[cpt,:]-ffit[cpt,:]
    temp=np.ma.masked_where(temp>999,temp)
    annual_climatology=np.ma.mean(temp,axis=1)
    return annual_climatology
###############################################################################
def read_kbafile(filename):
    f=open(filename)
    lines=f.readlines()
    f.close()
    if len(lines)<=6:
        print('Empty file: %s'%(filename))
        return
    mydate=re.search('\d{4}_\d{2}_\d{2}',filename).group()
    delta_2000=datetime.date(int(mydate[0:4]),int(mydate[5:7]),int(mydate[8:10]))-datetime.date(2000,1,1)
    delta_1Jan=datetime.date(int(mydate[0:4]),int(mydate[5:7]),int(mydate[8:10]))-datetime.date(int(mydate[0:4]),1,1)
    diny=(datetime.date(int(mydate[0:4])+1,1,1)-datetime.date(int(mydate[0:4]),1,1)).days
    nb_epochs=17257
    X=np.zeros(nb_epochs)
    Y=np.zeros(nb_epochs)
    Z=np.zeros(nb_epochs)
    postfit=np.zeros(nb_epochs)
    prefit=np.zeros(nb_epochs)
    timevec=np.zeros(nb_epochs)
    decyear=np.zeros(nb_epochs)
    for cpt in np.arange(len(lines)-6):
        rem=(str2float(lines[cpt+6][0:169].split()[0])*5)/86400
        timevec[cpt]=delta_2000.days+rem
        decyear[cpt]=int(mydate[0:4])+(delta_1Jan.days+rem)/diny
        X[cpt]=str2float(lines[cpt+6][0:169].split()[2])
        Y[cpt]=str2float(lines[cpt+6][0:169].split()[3])
        Z[cpt]=str2float(lines[cpt+6][0:169].split()[4])
        prefit[cpt]=str2float(lines[cpt+6][0:169].split()[11])
        postfit[cpt]=str2float(lines[cpt+6][0:169].split()[1])
    return timevec,decyear,X,Y,Z,prefit,postfit  
###############################################################################
def invert_for_tides(timevec,data,constituents='all',regularisation_weight='None',regularisation_type=0):
    " LSE of tidal components in data. \
      If tidal_constituents='all'(default), invert for M2, S2, K1 and O1 (in this order), you can choose a different selection of this 4 constituents \
      If regularisation weight = 'None' (default), performs an ordinary least square inversion. If regularisation weight > 0, permors a regularised least square inversion.\
      If regularisation type = 0 (default), performs a ridge regression with a penalty = regularisation weights \
      If regularisation type = 1, performs a lasso regression with a penalty = regularisation weights \
      If regularisation type > 0 & < 1, performs an elastic net regression with a penalty = regularisation weights\
      \
      Usage: coefs, predictions = invert_for_tides (timevec, data, constituents=['M2','K1'], regularisation_weight=5, regularisation_type=0)"
    NT=len(timevec)   
    M2=12.4206/24
    S2=12
    K1=23.9344/24
    O1=25.8193/24
    if constituents=='all':
        NP=9
        tides=np.ones((NT,NP))
        tides[:,1]=np.cos(timevec*2*np.pi/M2)
        tides[:,2]=np.sin(timevec*2*np.pi/M2)    
        tides[:,3]=np.cos(timevec*2*np.pi/S2)
        tides[:,4]=np.sin(timevec*2*np.pi/S2)
        tides[:,5]=np.cos(timevec*2*np.pi/K1)
        tides[:,6]=np.sin(timevec*2*np.pi/K1)    #
        tides[:,7]=np.cos(timevec*2*np.pi/O1)
        tides[:,8]=np.sin(timevec*2*np.pi/O1)
    else:
        NP=len(constituents)*2+1
        tides=np.ones((NT,NP))
        count=1
        for tc in constituents:
            tides[:,count]=np.cos(timevec*2*np.pi/eval(tc))
            tides[:,count+1]=np.sin(timevec*2*np.pi/eval(tc))
            count=count+2
    # Adjustment of a tidal model by OLS
    if regularisation_weight=='None':
        system = sm.OLS(data, tides)
        results = system.fit()
        coefs=results.params[:] 
        predictions=results.fittedvalues

    # Adjustment of a tidal model by weighted LS
    else:
        system = sm.OLS(data, tides)
        results = system.fit_regularized(method='elastic_net',alpha=regularisation_weight,L1_wt=regularisation_type)
        coefs=results.params[:] 
        predictions=results.fittedvalues
    return coefs,predictions   
##############################################################################
def read_addnorm_cmd(filename):
    f=open(filename)
    lines=f.readlines()
    f.close()
    nb_days=int(float(lines[1].split()[0]))
    mydays=[]
    for cpt in np.arange(nb_days):
        mydate=re.search('\d{4}_\d{2}_\d{2}',lines[4+cpt]).group()
        mydays.append(mydate[-2:])
    return mydays
##############################################################################
def write_setup_iter3(filename,myyear,mymonth,myday):   
    f=open(filename)
    lines=f.readlines()
    f.close()
    ofile=open(filename+'.iter3', 'w')
    for line in lines:
        if line=='MASCON            : N\n':    
            line='MASCON            : Y\n'
        if line[0:19]=='APRIORI_FILE      :':
            line='APRIORI_FILE      : IC_%s_%s_%s_iter3.vcv\n'%(myyear,mymonth,myday)
        if line=='USE_APRIORI_ICS   : 0  \n':    
            line='USE_APRIORI_ICS   : 15  \n'
        if line=='USE_APRIORI_MCS   : N\n':    
            line='USE_APRIORI_MCS   : Y\n'  
        if line=='MSC_APRIORI_FILE  : none\n':    
            line='MSC_APRIORI_FILE  : msc_apriori_%s_%s.vcv\n'%(myyear,mymonth) 
        ofile.write(line)
        
        
##############################################################################
def read_lowpass_file(apr_lowpass_file,msc_number):

    # read a lowpass filter file and extract the model values for a particular mascon (annual signal, lowpass values)\
    # P. Tregoning
    # 20 June 2022

    # get the list of decimal years of the epochs within the file
    decimal_year = np.array(np.genfromtxt(apr_lowpass_file ,skip_header=3,max_rows=1))
    #print(decimal_year)

    # get the model values of a particular mascon ...
    msc_model = np.array(np.genfromtxt(apr_lowpass_file,skip_header=msc_number+4,max_rows=1))
    #print(msc_model)
    
    return (decimal_year,msc_model[2],msc_model[3:5],msc_model[10:])
    


