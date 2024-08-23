'''
Author: Alex Crawford
Date Created: 28 Jul 2015
Date Modified: 12 Jun 2019 --> Modified for Python 3
                18 May 2020 --> Modified for using netcdf files instead of geotiffs
                19 Jan 2021 --> Added pickles as acceptable file input for masks
                11 Jun 2021 --> Added option for minimum displacement
                17 Apr 2023 --> modified for version 13
Purpose: Identify tracks that spend any point of their lifetime within a
bounding box defined by a list of (lon,lat) ordered pairs in a csv file.

Inputs: User must define the...
    Type of Tracks (typ) -- Cyclone Centers ("Cyclone") or System Centers ("System")
    Bounding Box Number (bboxnum) -- An ID for organizing directories
    Bounding Box Mask (bboxName) -- pathway for the mask to be used
    Versions of Module and Algorithm Run (e.g. 11.1, 12.4)
    Dates of interest
'''

'''********************
Import Modules
********************'''
# Import clock:
from time import perf_counter as clock
# Start script stopwatch. The clock starts running when time is imported
start = clock()

import pandas as pd
import numpy as np
import CycloneModule_13_2 as md
import os
import netCDF4 as nc
from scipy import interpolate
from settings import *
# import pickle5

current_file = os.path.basename(__file__)
print(f"Start {current_file}")

def maxDistFromGenPnt(data):
    return  np.max([md.haversine(data.lat[0],data.lat[i],data.lon[0],data.lon[i]) for i in range(len(data.lon))]) / 1000


'''*******************************************
Main Analysis
*******************************************'''
# print("Main Analysis")

# Load Mask
if suppath_detect.endswith(".nc"):
    bboxnc = nc.Dataset(suppath_detect)
    bbox = bboxnc[ncvar_height][:].data
elif suppath_detect.endswith(".pkl"):
    bbox = pd.read_pickle(suppath_detect)

if bboxmin == None:
    mask0 = bbox
else:
    mask0 = np.where((bbox >= bboxmin) & (bbox <= bboxmax), 1, 0)

mask = np.isin(mask0,values).reshape(mask0.shape)

# Set up output paths
try:
    os.chdir(inpath_system+"/BBox"+bboxnum_10)
except:
    os.mkdir(inpath_system+"/BBox"+bboxnum_10)
    os.chdir(inpath_system+"/BBox"+bboxnum_10)

try:
    os.chdir(inpath_system+"/BBox"+bboxnum_10+"/"+typ+"Tracks")
except:
    os.mkdir(inpath_system+"/BBox"+bboxnum_10+"/"+typ+"Tracks")
    os.chdir(inpath_system+"/BBox"+bboxnum_10+"/"+typ+"Tracks")

# Main Loop
mt = starttime
# while mt != endtime:
while mt[0] <= endtime[0] and mt[1] <= endtime[1]: # Simon changed
    # Extract date
    Y = str(mt[0])
    MM = months[mt[1]-1]
    M = mons[mt[1]-1]
    print(" " + Y + " - " + MM)

    # Load Tracks
    cs = pd.read_pickle(inpath_system+"/"+bboxmain+"/"+typ+"Tracks/"+Y+"/"+bboxmain+typ.lower()+"tracks"+Y+M+".pkl")
    # cs = pickle5.load(open(inpath_system+"/"+bboxmain+"/"+typ+"Tracks/"+Y+"/"+bboxmain+typ.lower()+"tracks"+Y+M+".pkl",'rb'))
    cs = [tr for tr in cs if ((tr.lifespan() >= minlifespan) and (tr.trackLength() >= mintracklength)) and (maxDistFromGenPnt(tr.data) >= mindisplacement)]

    trs = []
    for tr in cs: # For each track
        # Extract time and location
        xs = np.array(tr.data.x)
        ys = np.array(tr.data.y)
        hours = np.array(tr.data.time*24)

        # Interpolate to hourly
        f = interpolate.interp1d(hours,xs)
        xs2 = f(np.arange(hours[0],hours[-1])).astype(int)
        f = interpolate.interp1d(hours,ys)
        ys2 = f(np.arange(hours[0],hours[-1])).astype(int)
        # Test if at least one point is within the mask
        if mask[ys2,xs2].sum() > 0:
            trs.append(tr)

    # Save Tracks
    try:
        os.chdir(inpath_system+"/BBox"+bboxnum_10+"/"+typ+"Tracks/"+Y)
    except:
        os.mkdir(inpath_system+"/BBox"+bboxnum_10+"/"+typ+"Tracks/"+Y)
        os.chdir(inpath_system+"/BBox"+bboxnum_10+"/"+typ+"Tracks/"+Y)

    pd.to_pickle(trs,"BBox"+bboxnum_10+typ.lower()+"tracks"+Y+M+".pkl")

    # Increment Month
    mt = md.timeAdd(mt,monthstep)
# Print elapsed time
# print('Elapsed time:',round(clock()-start,2),'seconds')
# print("Complete")
print(f"{current_file} has finished")