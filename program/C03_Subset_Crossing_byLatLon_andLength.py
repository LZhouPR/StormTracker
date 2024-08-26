'''
Author: Alex Crawford
Date Created: 28 Jul 2015
Date Modified: 12 Jun 2019 --> Modified for Python 3
                18 May 2020 --> Modified for using netcdf files instead of geotiffs
                23 Jan 2023 --> Adapted to version 13
Purpose: Identify tracks that spend any point of their lifetime within a
bounding box defined by a list of (long,lat) ordered pairs in a csv file.

Inputs: User must define the...
    Type of Tracks (typ) -- Cyclone Centers ("Cyclone") or System Centers ("System")
    Bounding Box Number (bboxnum_27) -- An ID for organizing directories
    Bounding Box Mask (bboxName) -- pathway for the mask to be used
    Versions of Module and Algorithm Run (e.g. 7.8, 9.5)
    Spatial Resolution
    Dates of interest
    Minimum track length and lifespan
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
from settings import *

current_file = os.path.basename(__file__)
print(f"Start {current_file}")

'''*******************************************
Set up Environment
*******************************************'''
# path = "/Users/simontin/Desktop/cyclonetracking/testdata/Cressida"
# inpath_system = path+"/CycloneTracking/tracking13_2"
# outpath = inpath_system

'''*******************************************
Main Analysis
*******************************************'''
# print("Main Analysis")

# Load Mask
bboxnc = nc.Dataset(suppath_detect)
lats = bboxnc['lat'][:].data
lons = bboxnc['lon'][:].data

if bbox[-2] > bbox[-1]:
    mask = ((lats >= bbox[0]) & (lats <= bbox[1]) & ( (lons >= bbox[2]) | (lons <= bbox[3])) )
else:
    mask = ((lats >= bbox[0]) & (lats <= bbox[1]) & (lons >= bbox[2]) & (lons <= bbox[3]))

# Set up output paths
try:
    os.chdir(inpath_system+"/BBox"+bboxnum_27)
except:
    os.mkdir(inpath_system+"/BBox"+bboxnum_27)
    os.chdir(inpath_system+"/BBox"+bboxnum_27)

try:
    os.chdir(inpath_system+"/BBox"+bboxnum_27+"/"+typ+"Tracks")
except:
    os.mkdir(inpath_system+"/BBox"+bboxnum_27+"/"+typ+"Tracks")
    os.chdir(inpath_system+"/BBox"+bboxnum_27+"/"+typ+"Tracks")

# Main Loop
mt = starttime
while mt != endtime_nextmonth:
# while mt[0] <= endtime[0] and mt[1] <= endtime[1]: # Simon changed
    # Extract date
    Y = str(mt[0])
    MM = months[mt[1]-1]
    M = mons[mt[1]-1]
    print(" " + Y + " - " + MM)

    # Load Tracks
    cs = pd.read_pickle(inpath_system+"/"+bboxmain+typ+"Tracks/"+Y+"/"+bboxmain+typ+"tracks"+Y+M+".pkl")
    cs = [tr for tr in cs if ((tr.lifespan() > minlifespan) and (tr.trackLength() >= mintracklength))]

    trs = []
    for tr in cs: # For each track
        # Collect lats and longs
        xs = list(tr.data.x)
        ys = list(tr.data.y)

        # Prep while loop
        test = 0
        i = 0
        while test == 0 and i < len(xs):
            # If at any point the cyclone enters the bbox, keep it
            if mask[int(ys[i]),int(xs[i])] == 1:
                trs.append(tr)
                test = 1
            else:
                i = i+1

    # Save Tracks
    try:
        os.chdir(inpath_system+"/BBox"+bboxnum_27+"/"+typ+"Tracks/"+Y)
    except:
        os.mkdir(inpath_system+"/BBox"+bboxnum_27+"/"+typ+"Tracks/"+Y)
        os.chdir(inpath_system+"/BBox"+bboxnum_27+"/"+typ+"Tracks/"+Y)

    # pd.to_pickle(trs,"BBox"+bboxnum_27+typ.lower()+"tracks"+Y+M+".pkl")
    pd.to_pickle(trs, typ.lower()+"tracks"+Y+M+".pkl")

    # Increment Month
    mt = md.timeAdd(mt,monthstep)
    mt[2] = 1

# # Print elapsed time
# print('Elapsed time:',round(clock()-start,2),'seconds')
# print("Complete")

print(f"{current_file} has finished")
