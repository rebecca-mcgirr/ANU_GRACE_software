#!/usr/bin/env python3  
import numpy as np
import grace_utils as gr 
import datetime
import argparse
##########################################################################################
# Define arguments   
parser = argparse.ArgumentParser(description=' This programs aim to extract mascons in Antarctica. Mandatory arguments are the input mascon file and output mascon file in this order. You can choose the maximum latitude and density as options. Author : Julia Pfeffer, 13 January 2020 ')
parser.add_argument("in_masconfile",default='/Users/julia/GRACE/mascons/mascons_stage4_V003a',type=str,help="complete path to your input mascon file (string) ")
parser.add_argument("out_masconfile",default='/Users/julia/GRACE/mascons/mascons_stage4_V003a_Ant',type=str,help="complete path to your output mascon file  (string) ")
#parser.add_argument("--density",default=1000,type=float,help="grab all mascons with this density")
parser.add_argument("--region",default='Antarctica',type=str,help="grab all mascons in region Antarctica or Ocean")
#parser.add_argument("--latitude",default=-60,type=float,help="grab all mascons below this latitude")
# Read arguments
args = parser.parse_args()
# Declare arguments
in_masconfile = args.in_masconfile
out_masconfile = args.out_masconfile
#mydensity = args.density
myregion = args.region
############################################################################################
## Read input mascon file
iheaders,iPdata,_,iTdata=gr.read_mascon_file(in_masconfile)
## Select primaries with given density and latitudes
if myregion[:5]=='Antar':
    mydensity=1000
    mylatitude=-60
    Pmask=(iPdata[:,4]<mylatitude)&(iPdata[:,10]==mydensity)
elif myregion[:5]=='Ocean': # Grab entire ocean
    mydensity=1029
    mylatitude=90
    Pmask=(iPdata[:,4]<mylatitude)&(iPdata[:,10]==mydensity)
elif myregion[:5]=='SthnO':
    mydensity=1029
    mylatitude=-60
    Pmask=(iPdata[:,4]<mylatitude)&(iPdata[:,10]==mydensity)
elif myregion[:5]=='Arcti':
    mydensity=1029
    mylatitude=85
    Pmask=(iPdata[:,4]>mylatitude)&(iPdata[:,10]==mydensity)
else:
    print('ERROR')
inPdata=iPdata[(Pmask==True)]
outPdata=iPdata[(Pmask==False)]
print('len outPdata',len(outPdata))
## Find ternaries associated with selected primaries
nbt=len(iTdata)
ternaries=[]
for pcpt in np.arange(len(inPdata[:,0])):
    ternaries.append(iTdata[iTdata[:,9]==inPdata[pcpt,0]])
inTdata=np.asarray(np.vstack(ternaries))
print('Found ternaries',len(inTdata))
##########################################################################################
print('')
print('Create new primary and ternary arrays. This can take a few minutes. Relax and finish your tea :-)',str(datetime.datetime.now())[0:19])
##### Count the primaries 
nbp_out=len(outPdata[:,0]) # number of primaries outside polygons
#nbp_out=len(outPdata[:,0])-1 # number of primaries outside polygons
new_nbp=nbp_out+1 # total number of primaries
##### Build empty primary and ternary arrays
new_Pdata=np.array(np.zeros((new_nbp,14)),dtype=object)
new_Tdata=np.array(np.zeros((nbt,14)),dtype=object)
### Fill arrays for primaries and ternaies outside polygons
tcount=-1
for pcpt in np.arange(nbp_out):# primaries
    ternaries=iTdata[iTdata[:,9]==outPdata[pcpt,0]] # select ternaries outside polygons
    new_tip=np.shape(ternaries)[0] # number of ternary inside one primary 
    if new_tip==0:print(outPdata[pcpt,:])
    pcland=100*len(ternaries[:,8][ternaries[:,8]==1000])/new_tip # percentage of land in primary 
    new_Pdata[pcpt,0]=pcpt+1 # renumber the primary in ascending order
    new_Pdata[pcpt,1:13]=outPdata[pcpt,1:13]
    for tcpt in np.arange(new_tip):#ternaries
        tcount=tcount+1
        new_Tdata[tcount,0:9]=ternaries[tcpt,0:9] # ternary number, type, lat, lon, radius, area, alt., geoid height, density unchanged
        new_Tdata[tcount,9]=pcpt+1 # renumber associated primary
        new_Tdata[tcount,10]=pcpt+1 # renumber associated secondary
        new_Tdata[tcount,11:14]=ternaries[tcpt,11:14] # pcolour, scolour and region unchanged
    del ternaries, new_tip, pcland
### Fill arrays for primaries and ternaies inside polygons
for pcpt in np.arange(nbp_out,new_nbp):#primaries
    ternaries=inTdata# select ternaries for each primary inside polygon
    new_tip=np.shape(ternaries)[0] # number of ternaries inside primary
    pcland=100*len(ternaries[:,8][ternaries[:,8]==1000])/new_tip # percentage of land in primary
    new_Pdata[pcpt,0]=pcpt+1 # renumber primary
    new_Pdata[pcpt,1]='PLand'
    new_Pdata[pcpt,2]=1 #number of secondary in primary
    new_Pdata[pcpt,3]=new_tip #number of ternary in primary
    new_Pdata[pcpt,4]=np.mean(ternaries[:,2])# new latitude
    new_Pdata[pcpt,5]=np.mean(ternaries[:,3])# new longitude
    new_Pdata[pcpt,6]=np.mean(ternaries[:,4])# new radius
    new_Pdata[pcpt,7]=np.sum(ternaries[:,5]) # new area
    new_Pdata[pcpt,8]=np.mean(ternaries[:,6])# new altitude 
    new_Pdata[pcpt,9]=np.mean(ternaries[:,7])# new geoid height
    new_Pdata[pcpt,10]=1000 #  density defined as a selection criteria
    new_Pdata[pcpt,11]=pcland
    new_Pdata[pcpt,12]=0 # tidalmask defined in input (default = 0)
    new_Pdata[pcpt,13]='Antarctica'# new region name defined in polygon file
    for tcpt in np.arange(new_tip):
        tcount=tcount+1
        new_Tdata[tcount,0]=ternaries[tcpt,0] # number of ternary unchanged
        if mydensity==1000:new_Tdata[tcount,1]=str('TLand')
        if mydensity>1000:new_Tdata[tcount,1]=str('TDeep')
        new_Tdata[tcount,2:9]=ternaries[tcpt,2:9] # lat, lon, radius, area, alt., geoid height and density unchanged
        new_Tdata[tcount,9]=pcpt+1 # renumber associated primary
        new_Tdata[tcount,10]=pcpt+1 # renumber associated secondary
        #new_Tdata[tcount,11:13]=ternaries[tcpt,11:13] # pcolour and scolour unchanged
        if mydensity==1000:new_Tdata[tcount,11:13]=1600 # pcolour and scolour unchanged
        if mydensity>1000:new_Tdata[tcount,11:13]=0 # pcolour and scolour unchanged
        new_Tdata[tcount,13]=myregion # new region name defined in polygon file
    del ternaries, new_tip, pcland    
del iPdata,_,iTdata    
del inTdata, inPdata,outPdata
#################################################################################
print('')
print('Create headers:',str(datetime.datetime.now())[0:19])
npex=int(new_nbp)
ntex=int(nbt)
maxtip=int(np.max(new_Pdata[:,3]))
lastcomment='# Grab Antarctica. Timetag: %s \n'%(str(datetime.datetime.now())[0:19])
new_headers=gr.define_new_mascon_headers(npex,ntex,maxtip,iheaders,lastcomment)
del iheaders, ntex,maxtip,npex,lastcomment
#################################################################################
### WRITE NEW MASCON FILE
Pdata=new_Pdata
Tdata=new_Tdata
#################################################################################
print('Write mascon in file %s:'%(out_masconfile),str(datetime.datetime.now())[0:19])
ofile=open(out_masconfile, 'w')
for hcpt in np.arange(len(new_headers)):
    ofile.write(new_headers[hcpt])
# Write mascon data    
print('Write primary, secondary and ternary data in %s:'%(out_masconfile),str(datetime.datetime.now())[0:19])
tcount=0
for pcpt in np.arange(len(Pdata[:,0])):
    ofile.write('{:7d}  {:6s}{:8d}{:8d}{:10.4f}{:10.4f}{:11.1f}{:17.0f}{:1s}{:9.1f}{:9.1f}{:5.0f}{:1s}{:6.1f}{:6d}  {:14s}\n'.format(Pdata[pcpt,0],Pdata[pcpt,1],Pdata[pcpt,2],Pdata[pcpt,3],Pdata[pcpt,4],Pdata[pcpt,5],Pdata[pcpt,6],Pdata[pcpt,7],'.',Pdata[pcpt,8],Pdata[pcpt,9],Pdata[pcpt,10],'.',Pdata[pcpt,11],Pdata[pcpt,12],str(Pdata[pcpt,13]).rjust(14)))
    ofile.write('{:7d}  {:6s}{:8d}        {:10.4f}{:10.4f}{:11.1f}{:17.0f}{:1s}{:9.1f}{:9.1f}{:5.0f}{:1s}{:6d}{:6d}  {:14s}\n'.format(Pdata[pcpt,0],'S'+Pdata[pcpt,1][1:],Pdata[pcpt,3],Pdata[pcpt,4],Pdata[pcpt,5],Pdata[pcpt,6],Pdata[pcpt,7],'.',Pdata[pcpt,8],Pdata[pcpt,9],Pdata[pcpt,10],'.',Pdata[pcpt,0],Tdata[tcount,12],str([pcpt,13]).rjust(14)))
    for tcpt in np.arange(Pdata[pcpt,3]):
        ofile.write('{:7d}  {:6s}{:10.4f}{:10.4f}{:11.1f}{:14.0f}{:1s}{:9.1f}{:9.1f}{:5.0f}{:1s}{:6d}{:6d}{:6d}{:6d}         {:6s}\n'.format(Tdata[tcount,0],Tdata[tcount,1],Tdata[tcount,2],Tdata[tcount,3],Tdata[tcount,4],Tdata[tcount,5],'.',Tdata[tcount,6],Tdata[tcount,7],Tdata[tcount,8],'.',Tdata[tcount,9],Tdata[tcount,10],Tdata[tcount,11],Tdata[tcount,12],str(Tdata[tcount,13]).rjust(6)))
        tcount=tcount+1
ofile.close() 
del new_Pdata,new_Tdata,new_headers  
