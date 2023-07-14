import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import grace_utils as gr
import datetime as dt
###
def get_center(lon, lat, R = 6371e3):
    lon_r = np.radians(lon.astype(float))
    lat_r = np.radians(lat.astype(float))
    x =  R * np.cos(lat_r) * np.cos(lon_r)
    y = R * np.cos(lat_r) * np.sin(lon_r)
    z = R * np.sin(lat_r)
    R_center=np.sqrt(np.square(np.mean(x))+np.square(np.mean(y))+np.square(np.mean(z)))
    lon_center=(180/np.pi)*np.arctan(np.divide(np.mean(y),np.mean(x)))
    lat_center=90-(180/np.pi)*np.arccos(np.divide(np.mean(z),R))
    return lon_center,lat_center,R_center
headers,Pdata,_,Tdata=gr.read_mascon_file('/scratch/compute1/rebecca/test/V004_mascons/mascons_stage4_seas')
nbp=len(Pdata[:,0])
#oceans=['Ocean_NH', 'Ocean_SH', 'Arctic', 'SthnO', 'GrndOa', 'A400km1']
oceans=['Ocean', 'Shelf']
pselect=[]
pout=Pdata[Pdata[:,0]>=7,:]
print('select ocean ternaries',dt.datetime.now())
for cpt in np.arange(len(oceans)):
    pselect.append(Pdata[Pdata[:,13]==oceans[cpt],:])
pselect=np.vstack(pselect)
tselect=[]  
for cpt in np.arange(len(oceans)):
    tselect.append(Tdata[Tdata[:,9]==pselect[cpt,0],:])
tselect=np.vstack(tselect)  
tout=[]  
for cpt in np.arange(len(pout)):
    tout.append(Tdata[Tdata[:,9]==pout[cpt,0],:])
tout=np.vstack(tout)  
##arctic
print('select arctic ternaries',dt.datetime.now())
tarctic=tselect[tselect[:,2]>=62]
iarctic=np.where(tselect[:,2]>=62)[0]
nbt_arc=len(tarctic[:,0])
print(nbt_arc)
for tpt in np.arange(nbt_arc):
    if tpt%5000==0:
        print(tpt, dt.datetime.now())
    tselect[iarctic[tpt],1]='TOcean'
    tselect[iarctic[tpt],9]=int(1)
    tselect[iarctic[tpt],10]=int(1)
    tselect[iarctic[tpt],13]='ArcticO'
plon_arc,plat_arc,pr_arc=get_center(tarctic[:,3], tarctic[:,2], R = 6371e3)
parctic=np.array(np.zeros((1,14)),dtype=object)
parctic[0,0]=int(1)
parctic[0,1]='POcean'
parctic[0,2]=int(1)
parctic[0,3]=int(nbt_arc)
parctic[0,4]=plat_arc
parctic[0,5]=plon_arc
parctic[0,6]=pr_arc
parctic[0,7]=np.sum(tarctic[:,5])
parctic[0,8]=np.mean(tarctic[:,6])
parctic[0,9]=np.mean(tarctic[:,7])
parctic[0,10]=1029
parctic[0,11]=0
parctic[0,12]=0
parctic[0,13]='ArcticO'  
##antarctic
print('select antarctic ternaries',dt.datetime.now())
tsocean=tselect[tselect[:,2]<=-55]
isocean=np.where(tselect[:,2]<=-55)[0]
nbt_so=len(tsocean[:,0])
print(nbt_so)
for tpt in np.arange(nbt_so):
    if tpt%5000==0:
        print(tpt, dt.datetime.now())
    tselect[isocean[tpt],1]='TOcean'
    tselect[isocean[tpt],9]=int(2)
    tselect[isocean[tpt],10]=int(2)
    tselect[isocean[tpt],13]='Southern'
plon_so,plat_so,pr_so=get_center(tsocean[:,3], tsocean[:,2], R = 6371e3)
psocean=np.array(np.zeros((1,14)),dtype=object)
psocean[0,0]=int(2)
psocean[0,1]='POcean'
psocean[0,2]=int(1)
psocean[0,3]=int(nbt_so)
psocean[0,4]=plat_so
psocean[0,5]=plon_so
psocean[0,6]=pr_so
psocean[0,7]=np.sum(tsocean[:,5])
psocean[0,8]=np.mean(tsocean[:,6])
psocean[0,9]=np.mean(tsocean[:,7])
psocean[0,10]=1029
psocean[0,11]=0
psocean[0,12]=0
psocean[0,13]='Southern'  
##indian
print('select indian ternaries',dt.datetime.now())
tind1=tselect[(tselect[:,2]<=30)&(tselect[:,2]>-55)&(tselect[:,3]>=20)&(tselect[:,3]<120)]
tind2=tselect[(tselect[:,2]<=-7)&(tselect[:,2]>-23)&(tselect[:,3]>=120)&(tselect[:,3]<137)]
tind3=tselect[(tselect[:,2]<=-30)&(tselect[:,2]>-55)&(tselect[:,3]>=120)&(tselect[:,3]<146.5)]
tindian=np.vstack((tind1,tind2,tind3))
iind1=np.where((tselect[:,2]<=30)&(tselect[:,2]>-55)&(tselect[:,3]>=20)&(tselect[:,3]<120))[0]
iind2=np.where((tselect[:,2]<=-7)&(tselect[:,2]>-23)&(tselect[:,3]>=120)&(tselect[:,3]<137))[0]
iind3=np.where((tselect[:,2]<=-30)&(tselect[:,2]>-55)&(tselect[:,3]>=120)&(tselect[:,3]<146.5))[0]

nbt_io=len(tindian[:,0])
iindian=np.zeros(nbt_io)
iindian[0:len(iind1)]=iind1
iindian[len(iind1):len(iind1)+len(iind2)]=iind2
iindian[len(iind1)+len(iind2):len(iind1)+len(iind2)+len(iind3)]=iind3
print(nbt_io)
for tpt in np.arange(nbt_io):
    if tpt%5000==0:
        print(tpt, dt.datetime.now())    
    tselect[int(iindian[tpt]),1]='TOcean'
    tselect[int(iindian[tpt]),9]=int(3)
    tselect[int(iindian[tpt]),10]=int(3)
    tselect[int(iindian[tpt]),13]='Indian'
plon_io,plat_io,pr_io=get_center(tindian[:,3], tindian[:,2], R = 6371e3)
pindian=np.array(np.zeros((1,14)),dtype=object)
pindian[0,0]=int(3)
pindian[0,1]='POcean'
pindian[0,2]=int(1)
pindian[0,3]=int(nbt_io)
pindian[0,4]=plat_io
pindian[0,5]=plon_io
pindian[0,6]=pr_io
pindian[0,7]=np.sum(tindian[:,5])
pindian[0,8]=np.mean(tindian[:,6])
pindian[0,9]=np.mean(tindian[:,7])
pindian[0,10]=1029
pindian[0,11]=0
pindian[0,12]=0
pindian[0,13]='Indian'  
##pacific
print('select pacific ternaries',dt.datetime.now())
tpac1=tselect[(tselect[:,2]<62)&(tselect[:,2]>10)&(tselect[:,3]>=120)&(tselect[:,3]<275)]
tpac2=tselect[(tselect[:,2]<=10)&(tselect[:,2]>-6)&(tselect[:,3]>=120)&(tselect[:,3]<285)]
tpac3=tselect[(tselect[:,2]<=-6)&(tselect[:,2]>-35)&(tselect[:,3]>=140)&(tselect[:,3]<295)]
tpac4=tselect[(tselect[:,2]<=-35)&(tselect[:,2]>-55)&(tselect[:,3]>=146.5)&(tselect[:,3]<290)]
tpacific=np.vstack((tpac1,tpac2,tpac3,tpac4))
ipac1=np.where((tselect[:,2]<62)&(tselect[:,2]>10)&(tselect[:,3]>=120)&(tselect[:,3]<275))[0]
ipac2=np.where((tselect[:,2]<=10)&(tselect[:,2]>-6)&(tselect[:,3]>=120)&(tselect[:,3]<285))[0]
ipac3=np.where((tselect[:,2]<=-6)&(tselect[:,2]>-35)&(tselect[:,3]>=140)&(tselect[:,3]<295))[0]
ipac4=np.where((tselect[:,2]<=-35)&(tselect[:,2]>-55)&(tselect[:,3]>=146.5)&(tselect[:,3]<290))[0]


ipacific=np.vstack((np.reshape(ipac1,(len(ipac1),1)),np.reshape(ipac2,(len(ipac2),1)),np.reshape(ipac3,(len(ipac3),1)),np.reshape(ipac4,(len(ipac4),1))))
nbt_pac=len(tpacific[:,0])
print(nbt_pac)
for tpt in np.arange(nbt_pac):
    if tpt%5000==0:
        print(tpt, dt.datetime.now())
    tselect[ipacific[tpt],1]='TOcean'
    tselect[ipacific[tpt],9]=int(4)
    tselect[ipacific[tpt],10]=int(4)
    tselect[ipacific[tpt],13]='Pacific'
plon_pac,plat_pac,pr_pac=get_center(tpacific[:,3], tpacific[:,2], R = 6371e3)
ppacific=np.array(np.zeros((1,14)),dtype=object)
ppacific[0,0]=int(4)
ppacific[0,1]='POcean'
ppacific[0,2]=int(1)
ppacific[0,3]=int(nbt_pac)
ppacific[0,4]=plat_pac
ppacific[0,5]=plon_pac
ppacific[0,6]=pr_pac
ppacific[0,7]=np.sum(tpacific[:,5])
ppacific[0,8]=np.mean(tpacific[:,6])
ppacific[0,9]=np.mean(tpacific[:,7])
ppacific[0,10]=1029
ppacific[0,11]=0
ppacific[0,12]=0
ppacific[0,13]='Pacific'  
print('select atlantic ternaries',dt.datetime.now())
tatlantic=[]  
iatlantic=[]
oceans=['Ocean_NH', 'Ocean_SH', 'Arctic', 'SthnO', 'GrndOa', 'A400km1','Ocean','Africa']
for cpt in np.arange(len(oceans)):
    tatlantic.append(tselect[tselect[:,13]==oceans[cpt],:])
    nb_temp=len(tselect[tselect[:,13]==oceans[cpt],:])
    itemp=np.where(tselect[:,13]==oceans[cpt])[0]
    iatlantic.append(np.reshape(itemp,(nb_temp,1)))
    del itemp,nb_temp
iatlantic=np.vstack(iatlantic)
tatlantic=np.vstack(tatlantic) 
nbt_atl=len(tatlantic[:,0])
print(nbt_atl)
for tpt in np.arange(nbt_atl):
    if tpt%5000==0:
        print(tpt, dt.datetime.now())
    tselect[iatlantic[tpt],1]='TOcean'
    tselect[iatlantic[tpt],9]=int(5)
    tselect[iatlantic[tpt],10]=int(5)
    tselect[iatlantic[tpt],13]='Atlantic'
plon_atl,plat_atl,pr_atl=get_center(tatlantic[:,3], tatlantic[:,2], R = 6371e3)
patlantic=np.array(np.zeros((1,14)),dtype=object)
patlantic[0,0]=int(5)
patlantic[0,1]='POcean'
patlantic[0,2]=int(1)
patlantic[0,3]=int(nbt_atl)
patlantic[0,4]=plat_atl
patlantic[0,5]=plon_atl
patlantic[0,6]=pr_atl
patlantic[0,7]=np.sum(tatlantic[:,5])
patlantic[0,8]=np.mean(tatlantic[:,6])
patlantic[0,9]=np.mean(tatlantic[:,7])
patlantic[0,10]=1029
patlantic[0,11]=0
patlantic[0,12]=0
patlantic[0,13]='Atlantic' 

##pacific
#tpac1=tselect[(tselect[:,2]<62)&(tselect[:,2]>10)&(tselect[:,3]>=120)&(tselect[:,3]<275)]
#tpac2=tselect[(tselect[:,2]<=10)&(tselect[:,2]>-6)&(tselect[:,3]>=120)&(tselect[:,3]<285)]
#tpac3=tselect[(tselect[:,2]<=-6)&(tselect[:,2]>-35)&(tselect[:,3]>=140)&(tselect[:,3]<295)]
#tpac4=tselect[(tselect[:,2]<=-35)&(tselect[:,2]>-55)&(tselect[:,3]>=146.5)&(tselect[:,3]<290)]
#tpacific=np.vstack((tpac1,tpac2,tpac3,tpac4))
plt.figure(figsize=(10,10))
plt.cla()
plt.clf()
m = Basemap(llcrnrlon=0,llcrnrlat=-90, urcrnrlon=360, urcrnrlat=90, resolution = 'l')
Tx_arc, Ty_arc = m(tarctic[:,3],tarctic[:,2])
Tx_so, Ty_so = m(tsocean[:,3],tsocean[:,2])
Tx_ind, Ty_ind = m(tindian[:,3],tindian[:,2])
Tx_pac, Ty_pac = m(tpacific[:,3],tpacific[:,2])
Tx_atl, Ty_atl = m(tatlantic[:,3],tatlantic[:,2])
m.drawcoastlines(zorder=2)
#plt.scatter(Tx,Ty,5,inTdata[:,9],cmap=plt.cm.Set1,edgecolors='face',zorder=1)
plt.plot(Tx_arc,Ty_arc,'c. ')
plt.plot(Tx_so,Ty_so,'b. ')
plt.plot(Tx_ind,Ty_ind,'r. ')
plt.plot(Tx_pac,Ty_pac,'g. ')
plt.plot(Tx_atl,Ty_atl,'m. ')
#for cpt in np.arange(nbprim):
#    t=plt.text(Px[cpt]-0.75,Py[cpt],inPdata[cpt,0],color='k',fontweight='bold',fontsize=8)
#    t.set_bbox(dict(facecolor='white', alpha=0.6, edgecolor='none'))

plt.title('ocean')
plt.savefig('testoceans.png')     

# Change selected attributes
nbpnew=5+len(pout)
newPdata=np.array(np.zeros((nbpnew,14)),dtype=object)
newPdata[0,:]=parctic
newPdata[1,:]=psocean
newPdata[2,:]=pindian
newPdata[3,:]=ppacific
newPdata[4,:]=patlantic
nbtnew=len(tselect)+len(tout)
newTdata=np.array(np.zeros((nbtnew,14)),dtype=object)

ind = np.argsort( tselect[:,9] ); tselect = tselect[ind]
for tcpt in np.arange(len(tselect)):
    newTdata[tcpt,:]=tselect[tcpt,:]

tcount=len(tselect)
for cpt in np.arange(len(pout)):
    newPdata[5+cpt,:]=pout[cpt]  
    newPdata[5+cpt,0]=int(6+cpt)
    newPdata[5+cpt,0]=int(6+cpt)
    ttemp=tout[tout[:,9]==pout[cpt,0],:]
    for tcpt in np.arange(len(ttemp[:,0])):
       newTdata[tcount,:]=ttemp[tcpt,:]
       newTdata[tcount,9]=int(6+cpt)
       tcount=tcount+1
out_masconfile='/scratch/compute1/rebecca/test/V004_mascons/mascons_stage4_oceans'
newline='# Define Pacific, Atlantic, Indian, Socean and Arctic regions\n'
npex=len(newPdata)
ntex=len(newTdata)
maxtip=np.max(newPdata[:,3])
new_headers=gr.define_new_mascon_headers(npex,ntex,maxtip,headers,newline)
gr.write_mascon_file(out_masconfile,new_headers,newPdata,newTdata)    
