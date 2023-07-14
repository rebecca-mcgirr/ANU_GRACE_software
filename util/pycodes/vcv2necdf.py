#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 11:40:26 2019
@author: julia
"""
# import packages
import numpy as np
import os 
from fnmatch import fnmatch
import datetime as dt
from netCDF4 import Dataset
import grace_utils as gr
import argparse
## Define arguments   
parser = argparse.ArgumentParser(description=' This programs stacks monthly vcv file into one vcv file.')
parser.add_argument("datadir",default='/scratch/compute1/julia/solutions/monthly_solutions/',type=str,help="path to directory containing all vcv files to stack")
parser.add_argument("ncname",default='/scratch/compute1/julia/solutions/ncfiles/monthly_vcv_10cm.nc',type=str,help="name of output netcdf file")
parser.add_argument("--pattern",default='*.vcv',type=str,help="unix pattern to list vcv files")
## Read arguments
args = parser.parse_args()
## Declare arguments
ncname = args.ncname
datadir = args.datadir
pattern = args.pattern
# make a list of vcv file to stack
filelist=[]
for path, subdirs, files in os.walk(datadir):
    for name in files:
        if fnmatch(os.path.join(path, name), pattern): 
            filelist.append(os.path.join(path, name))            
filelist.sort()
nb_months=len(filelist)

# read first file to get number of msc
_,_,init_msc=gr.read_addnorm_vcv(filelist[0])
nbmsc=len(init_msc[:,0])
# create empy arrays to store values
nb_days_vec=99999*np.ma.ones(nb_months)
decyear=99999*np.ma.ones(nb_months)
daily_sat=99999*np.ma.ones((nb_months,31,24,3))
daily_dates=99999*np.ma.ones((nb_months,31,3))
monthly_msc=99999*np.ma.ones((nb_months,nbmsc,3))
count=0
# real all vcv files in directory
for file in filelist:
    dates,sat_params,msc_values=gr.read_addnorm_vcv(file)
    nbdays=len(dates[:,0])
    daily_dates[count,:nbdays,:,:]=dates
    daily_sat[count,:nbdays,:,:]=sat_params
    monthly_msc[count,:,:]=msc_values
    decyear[count]=(gr.toYearFraction(dt.date(dates[0,0],dates[0,1],dates[0,2]))+gr.toYearFraction(dt.date(dates[-1,0],dates[-1,1],dates[-1,2])))/2
    nb_days_vec[count]=nbdays
    count=count+1
    del dates,sat_params,msc_values,nbdays
decyear=np.ma.masked_where(decyear>9999,decyear)
nb_days_vec=np.ma.masked_where(nb_days_vec>9999,nb_days_vec)
daily_sat=np.ma.masked_where(daily_sat>9999,daily_sat)
daily_dates=np.ma.masked_where(daily_dates>9999,daily_dates)
monthly_msc=np.ma.masked_where(monthly_msc>9999,monthly_msc)
# Write netcdf file
ncw=Dataset(ncname,'w',format='NETCDF4')
ncw.description = 'stack of all vcv files in %s'%(datadir)
# Create dimensions
ncw.createDimension('nb_months', nb_months) # number of fitfiles each corresponding to a date
ncw.createDimension('nb_days', 31)
ncw.createDimension('nb_msc', nbmsc) # number mascons in each fitfile
ncw.createDimension('nb_sat_params', 24)
ncw.createDimension('nb_estimates', 3)
# Create variables
# decimal year
decyear_vec= ncw.createVariable('decyear', 'f4', ('nb_months'))
decyear_vec.units = 'decimal year'
decyear_vec.description='date of the month in decimal year'
decyear_vec[:]=decyear
# mascon list
mascon_vec= ncw.createVariable('mascons', 'f4', ('nb_msc'))
mascon_vec.units = '-'
mascon_vec.description='primary mascon number'
mascon_vec[:]=np.arange(1,1+nbmsc)
# nb days per solution
day_vec= ncw.createVariable('nbdays', 'f4', ('nbt'))
day_vec.units = '-'
day_vec.description='number of days per monthly solution'
day_vec[:]=nb_days_vec
# dates of the days used in the monthly solution
dates_arr= ncw.createVariable('daily_dates', 'f4', ('nb_months','nb_days','nb_estimates'))
dates_arr.units = '-'
dates_arr.description='dates of the days used in the monthly solution'
dates_arr[:,:,:]=daily_dates
# satellites parameters
sat_arr= ncw.createVariable('sat_params', 'f4', ('nb_months','nb_days','nb_sat_params','nb_estimates'))
sat_arr.units = '-'
sat_arr.description='24 parameters ([position(X,Y,Z), velocity(X,Y,Z), scale(X,Y,Z), bias(X,Y,Z)] * 2 satellites) for each day: a priori, postfit, sigma'
sat_arr[:,:,:,:]=daily_sat
# msc values
msc_arr= ncw.createVariable('msc_values', 'f4', ('nb_months','nb_msc','nb_estimates'))
msc_arr.units = 'm'
msc_arr.description='mascon values for each month: a priori, postfit,sigma'
msc_arr[:,:,:]=daily_dates
# Close NetCF file
ncw.close() 