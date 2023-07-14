import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib
import grace_utils as gr
import math
###

headers,Pdata,_,Tdata=gr.read_mascon_file('/Users/julia/GRACE/mascons/mascons_stage4_islands')
nbp=len(Pdata[:,0])
oceans=['Ocean_NH', 'Ocean_SH', 'Arctic', 'SthnO', 'GrndOa', 'A400km1']
pselect=[]
pout=Pdata[Pdata[:,0]>=7,:]
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
# Change selected attributes
nbpnew=1+len(pout)
newPdata=np.array(np.zeros((nbpnew,14)),dtype=object)
newPdata[0,0]=int(1)
newPdata[0,1]='POcean'
newPdata[0,2]=int(1)
newPdata[0,3]=len(tselect)
newPdata[0,4:10]=np.zeros(6)
newPdata[0,10]=100
newPdata[0,11]=0
newPdata[0,12]=0
newPdata[0,13]='Ocean'
nbtnew=len(tselect)+len(tout)
newTdata=np.array(np.zeros((nbtnew,14)),dtype=object)
for tcpt in np.arange(len(tselect)):
    newTdata[tcpt,:]=tselect[tcpt,:]
    newTdata[tcpt,1]='TOcean'
    newTdata[tcpt,9]=int(1)
    newTdata[tcpt,10]=int(1)
    newTdata[tcpt,11]=int(1)
    newTdata[tcpt,12]=int(1)
    newTdata[tcpt,13]='Ocean'
tcount=len(tselect)
for cpt in np.arange(len(pout)):
    newPdata[1+cpt,:]=pout[cpt]  
    newPdata[1+cpt,0]=int(2+cpt)
    newPdata[1+cpt,0]=int(2+cpt)
    ttemp=tout[tout[:,9]==pout[cpt,0],:]
    for tcpt in np.arange(len(ttemp[:,0])):
       newTdata[tcount,:]=ttemp[tcpt,:]
       newTdata[tcount,9]=int(2+cpt)
       tcount=tcount+1
out_masconfile='/Users/julia/GRACE/mascons/mascons_stage4_JP1'
newline='# merge Ocean_NH, Ocean_SH, Arctic, SthnO, GrndOa, and A400km1 \n'
npex=len(newPdata)
ntex=len(newTdata)
maxtip=np.max(newPdata[:,3])
new_headers=gr.define_new_mascon_headers(npex,ntex,maxtip,headers,newline)
gr.write_mascon_file(out_masconfile,new_headers,newPdata,newTdata)    