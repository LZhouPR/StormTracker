'''
Author: Alex Crawford
Date Created: 10 Mar 2015
Date Modified: 18 Apr 2016; 10 Jul 2019 (update for Python 3);
                10 Sep 2020 (switch from geotiff to netcdf), switch to uniform_filter from scipy.ndimage
                30 Sep 2020 (switch back to slower custom smoother because of what scipy does to NaNs)
                18 Feb 2021 (edited seasonal caluclations to work directly from months, not monthly climatology,
                             allowing for cross-annual averaging)
                09 Sep 2021: If a pre-existing file exists, this script will append new results
                            instead of overwriting for all years. Climatologies no longer in this script.
                01 Nov 2021: Added the possibility of appending prior years as will as subsequent years.
                23 Jan 2023: Adapted to version 13
Purpose: Calculate aggergate track density (Eulerian-Lagrangian hybrid) for either
cyclone tracks or system tracks.

User inputs:
    Path Variables, including the reanalysis (ERA, MERRA, CFSR)
    Track Type (typ): Cyclone or System
    Bounding Box ID (bboxnum): 2-digit character string
    Time Variables: when to start, end, the time step of the data
    Aggregation Parameters (minls, mintl, kSizekm)

Note: Units for track density are tracks/month/gridcell
'''

'''********************
Import Modules
********************'''
# Import clock:
from time import perf_counter as clock
start = clock()
import warnings
warnings.filterwarnings("ignore")

# print("Loading modules.")
import os
import pandas as pd
from scipy import ndimage
from scipy import interpolate
import numpy as np
import netCDF4 as nc
# import pickle5
import CycloneModule_13_2 as md
from settings import *

current_file = os.path.basename(__file__)
print(f"Start {current_file}")

'''*******************************************
Set up Environment
*******************************************'''
# print("Setting up environment.")
# bboxfull_27 = "BBox" + bboxnum_27
# bboxfull_27 = "BBox27" # use "" if performing on all cyclones; or BBox##
# typ = "System"
# verd = "13_2"

# path = "/Users/simontin/Desktop/cyclonetracking/testdata/Cressida"
# inpath = path+"/CycloneTracking/tracking"+version
# outpath = inpath+"/"+bboxfull_27

'''*******************************************
Define Variables
*******************************************'''
# Variables
vName = "trkden"
vunits = 'count'

'''*******************************************
Main Analysis
*******************************************'''
# print("Main Analysis")
# Ensure that folders exist to store outputs
try:
   os.chdir(outpath_agg+"/Aggregation"+typ)
except:
    os.mkdir(outpath_agg+"/Aggregation"+typ)
    os.chdir(outpath_agg+"/Aggregation"+typ)
try:
    os.chdir(outpath_agg+"/Aggregation"+typ+"/"+str(kSizekm)+"km")
except:
    os.mkdir(outpath_agg+"/Aggregation"+typ+"/"+str(kSizekm)+"km")
    os.chdir(outpath_agg+"/Aggregation"+typ+"/"+str(kSizekm)+"km")
priorfiles = os.listdir()

print("Step 1. Load Files and References")
# Read in attributes of reference files
params = pd.read_pickle(inpath_agg+"/cycloneparams.pkl")
# params = pickle5.load(open(inpath+"/cycloneparams.pkl",'rb'))
# timestep_list = params['timestep']
# spres = params['spres']

proj = nc.Dataset(suppath_reproj+"/EASE2_N0_"+str(spres)+"km_Projection.nc")
lats = proj['lat'][:]

kSize = int(kSizekm/spres) # This needs to be the full width ('diameter'), not the half width ('radius') for ndimage filters

print("Step 2. Aggregation requested for " + str(starttime[0]) + "-" + str(endtime_nextmonth[0]-1))
name = version+"_AggregationFields_Monthly_"+vName+".nc"
if name in priorfiles:
    prior = nc.Dataset(name)
    nextyear = int(np.ceil(prior['time'][:].max()))
    firstyear = int(np.floor(prior['time'][:].min()))
    if starttime[0] < firstyear: # If the desired time range starts before the prior years...
        if endtime_nextmonth[0] >= firstyear:
            startyear, endyear = starttime[0], firstyear
            print("Years " + str(firstyear) + "-"+str(nextyear-1) + " were already aggregated.\nAggregating for " + str(startyear) + "-" + str(endyear-1) + ".")
        else:
            raise Exception("There is a gap between the ending year requested ("+str(endtime_nextmonth[0]-1)+") and the first year already aggregated ("+str(firstyear)+"). Either increase the ending year or choose a different destination folder.")
    elif endtime_nextmonth[0] > nextyear: # If the desired range ends after the prior years...
        if starttime[0] <= nextyear:
            startyear, endyear = nextyear, endtime_nextmonth[0]
            print("Years " + str(firstyear) + "-"+str(nextyear-1) + " were already aggregated.\nAggregating for " + str(startyear) + "-" + str(endyear-1) + ".")
        else:
            raise Exception("There is a gap between the last year already aggregated ("+str(nextyear-1)+") and the starting year requested ("+str(starttime[0])+"). Either decrease the starting year or choose a different destination folder.")
    else:
        raise Exception("All requested years are already aggregated.")
else:
    startyear, endyear, firstyear, nextyear = starttime[0], endtime_nextmonth[0], starttime[0], endtime_nextmonth[0]

# Start at the earliest necessary time for ALL variables of interest
newstarttime = [startyear,1,1,0,0,0]
newendtime = [endyear,1,1,0,0,0]

vlists = []

mt = newstarttime
while mt != newendtime:
    # Extract date
    Y = str(mt[0])
    MM = months[mt[1]-1]
    M = mons[mt[1]-1]
    if MM == "Jan":
        print(" " + Y)

    ### LOAD TRACKS ###
    # Load Cyclone/System Tracks
    # cs = pickle5.load(open(inpath+"/"+bboxnum+"/"+typ+"Tracks/"+Y+"/"+bboxnum+typ.lower()+"tracks"+Y+M+".pkl",'rb'))
    cs = pd.read_pickle(inpath_agg+"/"+bboxfull_27+"/"+typ+"Tracks/"+Y+"/"+bboxfull_27+typ.lower()+"tracks"+Y+M+".pkl")

    ### LIMIT TRACKS & IDS ###
    # Limit to tracks that satisfy minimum lifespan and track length
    trs = [c for c in cs if ((c.lifespan() > minlifespan) and (c.trackLength() >= mintracklength))]

    ### CALCULATE FIELDS ###
    trk_field = np.zeros_like(lats)
    for tr in trs:
        # Extract time and location
        xs = np.array(tr.data.x)
        ys = np.array(tr.data.y)
        hours = np.array(tr.data.time*24)

        # Interpolate to hourly
        f = interpolate.interp1d(hours,xs)
        xs2 = f(np.arange(hours[0],hours[-1])).astype(int)
        f = interpolate.interp1d(hours,ys)
        ys2 = f(np.arange(hours[0],hours[-1])).astype(int)

        # Zip together ys and xs and find unique values
        yxs2 = np.transpose(np.vstack( (ys2,xs2) ))
        yxs3 = np.unique(yxs2,axis=0)

        # Record Existance of Track at each unique point
        for i in range(yxs3.shape[0]):
            x = yxs3[i,1]
            y = yxs3[i,0]

            trk_field[y,x] += 1

    ### SMOOTH FIELDS ###
    varFieldsm = ndimage.uniform_filter(trk_field,kSize,mode="nearest") # --> This cannot handle NaNs
    vlists.append(varFieldsm) # append to list

    # Increment Month
    mt = md.timeAdd(mt,monthstep,lys=1)
    mt[2] = 1

### SAVE FILE ###
print("Step 3. Write to NetCDF")
mnc = nc.Dataset(version+"_AggregationFields_Monthly_"+vName+"_NEW.nc",'w')
mnc.createDimension('y', lats.shape[0])
mnc.createDimension('x', lats.shape[1])
mnc.createDimension('time', (max(nextyear,newendtime[0])-min(firstyear,newstarttime[0]))*12)
mnc.description = 'Aggregation of cyclone track characteristics on monthly time scale.'

ncy = mnc.createVariable('y', np.float32, ('y',))
ncx = mnc.createVariable('x', np.float32, ('x',))
ncy.units, ncx.units = 'm', 'm'
ncy[:] = np.arange(proj['lat'].shape[0]*spres*1000/-2 + (spres*1000/2),proj['lat'].shape[0]*spres*1000/2, spres*1000)
ncx[:] = np.arange(proj['lat'].shape[1]*spres*1000/-2 + (spres*1000/2),proj['lat'].shape[1]*spres*1000/2, spres*1000)

# Add times, lats, and lons
nctime = mnc.createVariable('time', np.float32, ('time',))
nctime.units = 'years'
nctime[:] = np.arange(min(firstyear,newstarttime[0]),max(nextyear,newendtime[0]),1/12)

nclon = mnc.createVariable('lon', np.float32, ('y','x'))
nclon.units = 'degrees'
nclon[:] = proj['lon'][:]

nclat = mnc.createVariable('lat', np.float32, ('y','x'))
nclat.units = 'degrees'
nclat[:] = proj['lat'][:]

ncvar_height = mnc.createVariable(vName, np.float64, ('time','y','x'))
ncvar_height.units = vunits + ' -- Smoothing:' + str(kSizekm) + ' km'
vout = np.array(vlists)

name = version+"_AggregationFields_Monthly_"+vName+".nc"
if name in priorfiles: # Append data if prior data existed...
    if vout.shape[0] > 0: # ...and there is new data to be added
        prior = nc.Dataset(name)

        if starttime[0] < firstyear:
            ncvar_height[:] = np.concatenate( ( np.where(vout == 0,np.nan,vout) , prior[vName][:].data ) )
        else:
            ncvar_height[:] = np.concatenate( ( prior[vName][:].data , np.where(vout == 0,np.nan,vout) ) )

        mnc.close()

        os.remove(name) # Remove old file
        os.rename(version+"_AggregationFields_Monthly_"+vName+"_NEW.nc", name) # rename new file to standard name

else: # Create new data if no prior data existed
    ncvar_height[:] = np.where(vout == 0,np.nan,vout)
    mnc.close()
    os.rename(version+"_AggregationFields_Monthly_"+vName+"_NEW.nc", name) # rename new file to standard name

if (nextyear < endtime_nextmonth[0]) & (firstyear > starttime[0]):
    print("Completed aggregating " + str(startyear) + "-" + str(endyear-1)+".\nRe-run this script to aggregate " + str(nextyear) + "-" + str(endtime_nextmonth[0]-1) + ".")
else:
    print("Completed aggregating " + str(startyear) + "-" + str(endyear-1)+".")

print(f"{current_file} has finished")