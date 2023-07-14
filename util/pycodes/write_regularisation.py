#!/usr/bin/env python3
import numpy as np
import grace_utils as gr
import argparse
##

parser = argparse.ArgumentParser('Write regularisation file for a given mascon file. ')
parser.add_argument('mcfile', default="None",  type=str, help='complete path to input masconfile')
parser.add_argument('regfile', default="None",  type=str, help='complete path to output regularisation file')
parser.add_argument('--IceSheet', default=0.100,  type=float, help='constraint value for ice seets in meters')
parser.add_argument('--Casp',default=0.100,type=float,help='constraint value for the Caspian Sea in meters')
parser.add_argument('--Amazon',default=0.050,type=float,help='constraint value for the Amazon in meters')
parser.add_argument('--Land',default=0.020,type=float,help='constraint value for land in meters')
parser.add_argument('--Ocean',default=0.020,type=float,help='constraint value for oceans in meters')
   
## Read arguments
args = parser.parse_args()

## Declare arguments
mcfile = args.mcfile
regfile = args.regfile
val_IS = args.IceSheet
val_C = args.Casp
val_A = args.Amazon
val_L = args.Land
val_O = args.Ocean


headers,Pdata,_,_=gr.read_mascon_file(mcfile)
nbp=len(Pdata[:,0])
   
newline='# Sigmas : ocean= %5.3f m, land= %5.3f m, amazon= %5.3f m, caspian= %5.3f, ice-sheets= %5.3f m \n'%(val_O,val_L,val_A,val_C,val_IS)
ofile=open(regfile, 'w')
ofile.write(headers[0])
ofile.write(newline)
for pcpt in np.arange(nbp):
    for cpt in np.arange(nbp):
        if cpt==pcpt:
            if Pdata[pcpt,10]==1029:
                ofile.write('{:24.15f}  '.format(1/(np.square(val_O))))
            else:
                if Pdata[pcpt,13]=='Greenland':
                    ofile.write('{:24.15f}  '.format(1/(np.square(val_IS))))
                elif Pdata[pcpt,13]=='Antarctica':
                    ofile.write('{:24.15f}  '.format(1/(np.square(val_IS)))) 
                elif Pdata[pcpt,13]=='Thwaites':
                    ofile.write('{:24.15f}  '.format(1/(np.square(val_IS))))   
                elif Pdata[pcpt,13]=='Pine_Is':
                    ofile.write('{:24.15f}  '.format(1/(np.square(val_IS)))) 
                elif Pdata[pcpt,13]=='AMAZONAS':
                    ofile.write('{:24.15f}  '.format(1/(np.square(val_A))))     
                elif Pdata[pcpt,13]=='Casp':
                    ofile.write('{:24.15f}  '.format(1/(np.square(val_C))))                
                else:
                    ofile.write('{:24.15f}  '.format(1/(np.square(val_L))))
        else:
            ofile.write('{:24.15f}  '.format(0))
    ofile.write('\n')
ofile.close()
