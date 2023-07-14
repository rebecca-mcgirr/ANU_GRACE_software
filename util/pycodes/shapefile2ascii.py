#!/usr/bin/env python3
import fiona
import numpy as np
import datetime
import pyproj
import argparse
###############################################################################
import grace_utils as gr
################################################################################
### READ ARGUMENTS    
################################################################################ 
### Define arguments   
parser = argparse.ArgumentParser(description=' This program reads shapefiles, extract polygons and write the polygons coordinates in an ascii file legible by the GRACE software. Polygons are named with a prefix (default=\'POL\') and ascending numbers.')
parser.add_argument("shp_filename",default="None",type=str,help="complete path to your shapefile (string) ")
parser.add_argument("--outfilename",default="my_ascii_file.polygons",type=str,help="complete path to your ascii polygon file (default = my_ascii_file.polygons)")
parser.add_argument("--prefix",default="POL",type=str,help="Prefix used to name your polygons (default = \'POL\')")
### Read arguments
args = parser.parse_args()
### Declare arguments
shp_filename = args.shp_filename
#in_epsg = args.in_epsg
outfilename = args.outfilename
prefix = args.prefix
############################################################################### 
print('Start shapefile2ascii.py: %s'%(str(datetime.datetime.now())[0:19]))
# read shapefile
print('Read \'%s\': %s'%(shp_filename,str(datetime.datetime.now())[0:19]))
sf= fiona.open(shp_filename)
nbpol=len(sf)
lenpol=len(str(nbpol))
## Define output values
polygons=[]
coords=[]
nbv=[]
legend=[]
for cpt in range(nbpol):
    if sf[cpt]['geometry']['type']=='Polygon':
        polygons.append(prefix+str(cpt+1).zfill(lenpol))
        coords.append(sf[cpt]['geometry']['coordinates'][0])
        nbv.append(len(coords[cpt]))    
    elif sf[cpt]['geometry']['type']=='MultiPolygon': 
        multi=sf[cpt]['geometry']['coordinates']
        mergedpol=[]
        for poly in multi:
        	for pcoord in poly[0]:
        		mergedpol.append((pcoord[0],pcoord[1]))
        polygons.append(prefix+str(cpt+1).zfill(lenpol))
        coords.append(mergedpol)
        nbv.append(len(coords[cpt]))    
        
## Define positive longitudes in WGS84
if bool(sf.crs):
    inProj = pyproj.Proj(sf.crs)
    outProj = pyproj.Proj(init='epsg:4326')
    if inProj.srs==outProj.srs:
        print('Geographical coordinates already in WGS84. Convert negative longitudes to positive longitudes if any: %s.'%(str(datetime.datetime.now())[0:19]))
        coords=gr.convert_longitudes_in_tuple(coords,flag='to_positive')
    else:
        print('Shapefile crs = %s'%(str(sf.crs)))
        print('Geographical coordinates converted to WGS84 (positive longitudes only): %s.'%(str(datetime.datetime.now())[0:19]))
        coords=gr.project_coordinates_in_tuple(coords,inProj,outProj)
        coords=gr.convert_longitudes_in_tuple(coords,flag='to_positive')
else:
	print('CRS not found in shapefile. Polygons coordinates cannot be converted to WGS84. %s'%(str(datetime.datetime.now())[0:19]))        
### Write coordinates in a ascii polygon file similar to mascon.polygons
print('Write polygons headers and coordinates in \'%s\': %s.'%(outfilename,str(datetime.datetime.now())[0:19])) 
ofile=open(outfilename, 'w')
for cpt in range(nbpol):
    ofile.write('{:5d}   {:10s}  {:10.0f}{:1s}\n'.format(int(nbv[cpt]),polygons[cpt].rjust(10),1000,'.'))
    for cpt2 in np.arange(nbv[cpt]):
        ofile.write('{:15.8f} {:15.8f}\n'.format(coords[cpt][cpt2][1],coords[cpt][cpt2][0]))
ofile.close() 
del ofile
### Legend
print('Write legend about polygons in \'ascii_polygons.legend\': %s'%(str(datetime.datetime.now())[0:19]))
for cpt in range(nbpol):
	ofile=open('ascii_polygons.legend', 'w')
	mykeys=sf[0]['properties'].keys()
	ofile.write('{:s}\t'.format('#polname'))
	for mykey in mykeys:
		ofile.write('{:s}\t'.format(mykey))
	ofile.write('\n')	
	for cpt in range(nbpol):
		ofile.write('{:s}\t'.format(prefix+str(cpt+1).zfill(lenpol)))
		myproperties=sf[cpt]['properties']
		for mykey in mykeys:
			ofile.write('{:s}\t'.format(str(myproperties[mykey])))
		ofile.write('\n')				
ofile.close() 
del ofile
print('End shapefile2ascii.py: %s'%(str(datetime.datetime.now())[0:19]))
################################################################################