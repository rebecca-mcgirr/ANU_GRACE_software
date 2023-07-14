#!/usr/bin/env python3
import numpy as np
import grace_utils as gr
import datetime as dt
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib
import math
##
def compute_distance(plat,plon,latvec,lonvec):
    nb=len(latvec)
    R = 6373.0
    lat1 = np.radians(plat)*np.ones(nb)
    lon1 = np.radians(plon)*np.ones(nb)
    lat2 = np.radians(latvec)
    lon2 = np.radians(lonvec)
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.square(np.sin(dlat / 2)) + np.cos(lat1) * np.cos(lat2) * np.square(np.sin(dlon / 2))
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))  
    distance = R * c
    return distance
##
inputfile='/Users/julia/GRACE/mascons/mascons_stage4_JP29'
headers,Pdata,Sdata,Tdata=gr.read_mascon_file(inputfile)
nbp=len(Pdata[:,0])
nbt=len(Tdata[:,0])
##
Pselect=Pdata[Pdata[:,3]<95]
Pout=Pdata[Pdata[:,3]>=95]
nbp_select=len(Pselect[:,0])
colors=[]
cm20=plt.cm.tab20
for cpt in np.arange(math.ceil(nbp)/cm20.N):
    for cpt2 in np.arange(cm20.N):
        colors.append(cm20(cpt2))
cm=matplotlib.colors.ListedColormap(colors) 
###
mindist=[]
closest_primary=[]
small_primary=[]
for pcpt in np.arange(nbp_select):
    Pnot=Pdata[(Pdata[:,0]!=Pselect[pcpt,0])&(Pdata[:,10]==Pselect[pcpt,10])]
    dist=compute_distance(Pselect[pcpt,4],Pselect[pcpt,5],Pnot[:,4].astype(float),Pnot[:,5].astype(float))
    mindist.append(np.min(dist))
    small_primary.append(Pselect[pcpt,0])
    closest_primary.append(Pnot[dist==np.min(dist),0])
    #### colorbar
    small_ternarys=Tdata[Tdata[:,9]==small_primary[pcpt],:]
    close_ternarys=Tdata[Tdata[:,9]==closest_primary[pcpt],:]
    two_ternarys=np.vstack((small_ternarys,close_ternarys))
    minlon=np.min(two_ternarys[:,3])-5
    maxlon=np.max(two_ternarys[:,3])+5
    minlat=np.min(two_ternarys[:,2])-5
    maxlat=np.max(two_ternarys[:,2])+5
    inTdata=Tdata[(Tdata[:,3]>=minlon)&(Tdata[:,3]<=maxlon)&(Tdata[:,2]>=minlat)&(Tdata[:,2]<=maxlat)]   
    prim=np.unique(inTdata[:,9])
    nbprim=len(prim)
    inPdata=np.array(np.zeros((nbprim,14)),dtype=object)
    for cpt in np.arange(nbprim):
        inPdata[cpt,:]=Pdata[Pdata[:,0]==prim[cpt]]
    plt.figure(figsize=(10,10))
    plt.cla()
    plt.clf()
    m = Basemap(projection='cyl',llcrnrlon=minlon,llcrnrlat=minlat, urcrnrlon=maxlon, urcrnrlat=maxlat, resolution = 'l')
    Tx, Ty = m(inTdata[:,3],inTdata[:,2])
    Txsmall, Tysmall = m(small_ternarys[:,3],small_ternarys[:,2])
    Txclose, Tyclose = m(close_ternarys[:,3],close_ternarys[:,2])
    Px, Py = m(inPdata[:,5],inPdata[:,4])
    Pxsmall, Pysmall = m(Pdata[Pdata[:,0]==small_primary[pcpt],5],Pdata[Pdata[:,0]==small_primary[pcpt],4])
    Pxclose, Pyclose = m(Pdata[Pdata[:,0]==closest_primary[pcpt],5],Pdata[Pdata[:,0]==closest_primary[pcpt],4])
    m.drawcoastlines(zorder=2)
    plt.scatter(Tx,Ty,2,inTdata[:,9],cmap=cm,vmin=np.min(inTdata[:,9]),vmax=np.max(inTdata[:,9]),edgecolors='face',zorder=1)
    plt.plot(Txsmall,Tysmall,'ro') 
    #t=plt.text(Pxsmall,Pysmall,small_primary[pcpt],color='k',fontweight='bold',fontsize=8)
    #t.set_bbox(dict(facecolor='white', alpha=0.6, edgecolor='none'))
    #t=plt.text(Pxclose,Pyclose,closest_primary[pcpt][0],color='k',fontweight='bold',fontsize=8)
    #t.set_bbox(dict(facecolor='white', alpha=0.6, edgecolor='none'))
    plt.plot(Txclose,Tyclose,'bo') 
    for cpt in np.arange(nbprim):
        if (inPdata[cpt,4]<maxlat)&(inPdata[cpt,4]>minlat)&(inPdata[cpt,5]<maxlon)&(inPdata[cpt,5]>minlon):
            t=plt.text(Px[cpt],Py[cpt],inPdata[cpt,0],color='k',fontweight='bold',fontsize=8)
            t.set_bbox(dict(facecolor='white', alpha=0.6, edgecolor='none'))
 
    plt.title('small prim: %s %d tern, close prim: %s %d tern, distance = %d km'%(Pdata[Pdata[:,0]==small_primary[pcpt],13][0],Pdata[Pdata[:,0]==small_primary[pcpt],3],Pdata[Pdata[:,0]==closest_primary[pcpt],13][0],Pdata[Pdata[:,0]==closest_primary[pcpt],3],int(mindist[pcpt])))
    plt.savefig('small_primarys_%d.png'%(pcpt))  
        
    
#dt2p=np.zeros(nbt)
#tap=np.zeros((nbt,2))
#tcount=0
#for cpt in np.arange(nbp):
#    if cpt%100==0:
#        print(cpt, dt.datetime.now())
#    ntip=Pdata[cpt,3]
#    prim=Pdata[cpt,0]
#    dist=compute_distance(Pdata[cpt,4],Pdata[cpt,5],Tdata[Tdata[:,9]==prim,2].astype(float),Tdata[Tdata[:,9]==prim,3].astype(float))
#    dt2p[tcount:tcount+ntip]=dist
#    tap[tcount:tcount+ntip,0]=prim*np.ones(ntip)
#    tap[tcount:tcount+ntip,1]=Tdata[Tdata[:,9]==prim,0]
#    tcount=tcount+ntip
####
#dt2p_select=dt2p[dt2p>500]    
#tap_select=tap[dt2p>500] 
#nb_select=len(tap_select)
###
#closest_primary=np.zeros(nb_select)
#new_distance=np.zeros(nb_select)
#for cpt in np.arange(nb_select):
#    mytern=tap_select[cpt,1]
#    mydensity=Tdata[Tdata[:,0]==mytern,8]
#    tlat=Tdata[Tdata[:,0]==mytern,2][0]
#    tlon=Tdata[Tdata[:,0]==mytern,3][0]
#    pnumbers=Pdata[Pdata[:,10]==mydensity,0]
#    plats=Pdata[Pdata[:,10]==mydensity,4].astype(float)
#    plons=Pdata[Pdata[:,10]==mydensity,5].astype(float)
#    mydist=compute_distance(tlat,tlon,plats,plons)
#    closest_primary[cpt]=pnumbers[mydist==np.min(mydist)]
#    new_distance[cpt]=mydist[mydist==np.min(mydist)]
#tern2merge=[] 
#for cpt in np.arange(nb_select):
#    if tap_select[cpt,0]!=closest_primary[cpt]:
#        if dt2p_select[cpt]-new_distance[cpt]>100:
#            tern2merge.append([tap_select[cpt,1],closest_primary[cpt]])
#tern2merge=np.vstack(tern2merge)     
#outfile='/Users/julia/GRACE/mascons/tern2merge_400km'
#ofile=open(outfile, 'w')
#for tpt in np.arange(len(tern2merge)):
#    ofile.write('{d} {d}\n'.format(int(tern2merge[tpt,0]),int(tern2merge[tpt,1])))
#ofile.close()      
##outfile='/Users/julia/GRACE/mascons/t2p_distances_400km'
##ofile=open(outfile, 'w')
##ofile.write('# primary number, ternary number, distance \n')
##for tpt in np.arange(len(dt2p_select)):
##        ofile.write('{:10d}\t{:10d}\t{:10.5f}\n'.format(int(tap_select[tpt,0]),int(tap_select[tpt,1]),dt2p_select[tpt]))
##ofile.write('\n')
##ofile.close() 
##outfile='/Users/julia/GRACE/mascons/t2p_distances'
##ofile=open(outfile, 'w')
##ofile.write('# primary number, ternary number, distance \n')
##for tpt in np.arange(nbt):
##        ofile.write('{:10d}\t{:10d}\t{:10.5f}\n'.format(int(tap[tpt,0]),int(tap[tpt,1]),dt2p[tpt]))
##ofile.write('\n')
##ofile.close()    
#
