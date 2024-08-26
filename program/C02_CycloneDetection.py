'''
Author: Alex Crawford
Date Created: 20 Jan 2015
Date Modified: 10 Sep 2020 -> Branch from 12_1 --> kernel size is now based on km instead of cells
                16 Dec 2020 --> updated comments
                13 Jan 2021 --> changed  when masking for elevation happens -- after minima are detected, not before
                02 Mar 2022 --> transferred some constants to the module file to save space
                            --> changed the output for cyclone fields from a
                            single field per pickled file to to a list of fields for an
                            entire month in each pickled file
                14 Nov 2022 --> Added an if statement to make this more functional for the Southern Hemisphere
                23 Jan 2023 --> Branch from 12_4 --> added kSizekm to cycloneparams.pkl and switched from "surf" to "field"
                03 Apr 2023 --> Improved flexibility running the code using prior cyclone detection for a new temporal
                            resolution while tracking -- there were some bugs if the original run was at a finer temporal resolution
                            and a later run was at a coarser resolution

Purpose: Given a series of sea level pressure fields in netcdf files, this
    script performs several steps:
    1) Identify closed low pressure centers at each time step
    2) Store information to characterize these centers at each time step
    3) Identify multi-center and single-center cyclone systems
  Steps in the tracking part:
    4) Associate each cyclone center with a corresponding center in the
        previous time step (when applicable)
    5) Combine the timeseries of cyclone center charcteristics into a data frame for each track
    6) Record cyclone life cycle events (genesis, lysis, splits, merges, secondary genesis)

User Inputs: paths for inputs, desired projection info, various detection/tracking parameters
'''

'''********************
Import Modules
********************'''
# Import clock:
import time
# Start script stopwatch. The clock starts running when time is imported
start = time.perf_counter()
import os
import copy
import numpy as np
import netCDF4 as nc
import CycloneModule_13_2 as md
import warnings
from settings import *
from datetime import datetime

current_file = os.path.basename(__file__)
print(f"Start {current_file}")

np.seterr(all='ignore') # This mutes warnings from numpy
warnings.filterwarnings('ignore',category=DeprecationWarning)

'''*******************************************
Main Analysis
*******************************************'''
##### Ensure that folders exist to store outputs #####
detpath = outpath_detect+"/detection"+version
trkpath = outpath_detect+"/tracking"+version
try:
    os.chdir(outpath_detect)
except:
    os.mkdir(outpath_detect)
    os.chdir(outpath_detect)
try:
    os.chdir(detpath)
except:
    os.mkdir(detpath)
    os.chdir(detpath)
    os.mkdir("CycloneFields")
try:
    os.chdir(trkpath)
except:
    os.mkdir(trkpath)
    os.chdir(trkpath)
    os.mkdir("CycloneTracks")
    os.mkdir("ActiveTracks")
    os.mkdir("SystemTracks")

for y in range(starttime[0],endtime_nextmonth[0]+1):
    Y = str(y)

    # Cyclone Tracks
    try:
        os.chdir(trkpath+"/CycloneTracks/"+Y)
    except:
        os.mkdir(trkpath+"/CycloneTracks/"+Y)

    # Active Tracks
    try:
        os.chdir(trkpath+"/ActiveTracks/"+Y)
    except:
        os.mkdir(trkpath+"/ActiveTracks/"+Y)

    # System Tracks
    try:
        os.chdir(trkpath+"/SystemTracks/"+Y)
    except:
        os.mkdir(trkpath+"/SystemTracks/"+Y)

##### Read in attributes of reference files #####
projnc = nc.Dataset(suppath_detect)

lats = projnc['lat'][:].data
lons = projnc['lon'][:].data
yDist = projnc['yDistance'][:].data
xDist = projnc['xDistance'][:].data
elev = projnc['z'][:]

# Generate mask based on latitude and elevation
if minlat >= 0:
    mask = np.where((elev > maxelev) | (lats < minlat),np.nan,0)
else:
    mask = np.where((elev > maxelev) | (lats > minlat),np.nan,0)

# Convert kernel size to grid cells
kSize = int(2*kSizekm/spres)+1

# Convert max speed to max distance
maxdist = maxspeed*1000*timestep_list[3]

# Save Parameters
params = dict({"path":trkpath,"timestep":timestep_list, "dateref":dateref, "minfield":minfield,
    "maxfield":maxfield,"kSize":kSize, "kSizekm":kSizekm,"nanThresh":nanThresh, "d_slp":d_slp, \
    "d_dist":d_dist, "maxelev":maxelev, "minlat":minlat, "contint":contint,
    "mcctol":mcctol, "mccdist":mccdist, "maxspeed":maxspeed, "red":red, "spres":spres})
pd.to_pickle(params,trkpath+"/cycloneparams.pkl")

# ##### The actual detection and tracking #####
# print("Cyclone Detection & Tracking")
# # Print elapsed time
# print(' Elapsed time:',round(time.perf_counter()-start,2),'seconds -- Starting first month')

# Load netcdf for initial time
ncf = nc.Dataset(inpath_detect+"/"+ra+"_EASE2_N0_"+str(spres)+"km_"+var+"_Hourly_"+str(starttime[0])+md.dd[starttime[1]-1]+".nc")
tlist = ncf['time'][:].data

# Try loading cyclone field objects for the initial month -- if none exist, make an empty list
try:
    cflist = pd.read_pickle(detpath+"/CycloneFields/CF"+str(starttime[0])+md.dd[starttime[1]-1]+".pkl")
    cftimes = np.array([cf.time for cf in cflist])
except:
    cflistnew = []

t = copy.deepcopy(starttime)
last_day_check = days[-1] + timedelta(days=1)
while t != endtime:
    # Extract date
    Y = str(t[0])
    MM = md.mmm[t[1]-1]
    M = md.dd[t[1]-1]
    date = Y+M+md.dd[t[2]-1]+"_"+md.hhmm[t[3]]
    days = md.daysBetweenDates(dateref,t)
    if round(days) == days:
        current_time = datetime.now().time()
        print(current_time, pd.to_datetime(date, format = '%Y%m%d_%H%M').strftime('%Y-%m-%d'))

    # Load field
    try: # If the cyclone field has already been calculated, no need to repeat
        cf = cflist[np.where(cftimes == days)[0][0]]
    except:
        field = ncf[ncvar_p][np.where(tlist == md.daysBetweenDates(dateref,t)*24)[0][0],:,:]
        field = np.where((field < minfield) | (field > maxfield), np.nan, field)

        # Create a cyclone field object
        cf = md.cyclonefield(days)

        # Identify cyclone centers
        cf.findCenters(field, mask, kSize, nanThresh, d_slp, d_dist, yDist, xDist, lats, lons) # Identify Cyclone Centers

        # Calculate cyclone areas (and MCCs)
        cf.findAreas(field+mask, contint, mcctol, mccdist, lats, lons, kSize) # Calculate Cyclone Areas

        # Append to lists
        cflistnew.append( cf )

    # Track Cyclones
    if t == starttime: # If this is the first time step, must initiate tracking
        if prior == 0: #If this is the first time step and there are no prior months
            ct, cf.cyclones = md.startTracks(cf.cyclones)

        else: #If this is the first time step but there is a prior month
            # Identify date/time of prior timestep
            tp = md.timeAdd(t,[-i for i in timestep_list])

            # Load cyclone tracks from prior month
            ct = pd.read_pickle(trkpath+"/ActiveTracks/"+str(tp[0])+"/activetracks"+str(tp[0])+md.dd[tp[1]-1]+".pkl")

           # Load cyclone field from prior time step
            cfs = pd.read_pickle(detpath+"/CycloneFields/CF"+str(tp[0])+md.dd[tp[1]-1]+".pkl")
            cfi = np.where( np.array( [c.time for c in cfs] ) == md.daysBetweenDates(dateref,tp) )[0][0]
            cf1 = pd.read_pickle(detpath+"/CycloneFields/CF"+str(tp[0])+md.dd[tp[1]-1]+".pkl")[cfi]

            # Reset TIDs for active cyclone tracks
            md.realignPriorTID(ct,cf1)

            # Start normal tracking
            ct, cf = md.trackCyclones(cf1,cf,ct,maxdist,red,timestep_list[3])

    else: #If this isn't the first time step, just keep tracking
        ct, cf = md.trackCyclones(cf1,cf,ct,maxdist,red,timestep_list[3])

    # Increment time step indicator
    t = md.timeAdd(t,timestep_list)
    cf1 = copy.deepcopy(cf)

    # Save Tracks & Fields (at the end of each month)
    # if t[2] == 1 and t[3] == 0: # If the next timestep is the 0th hour of the 1st day of a month,
    current_day = pd.to_datetime(f"{t[0]}-{t[1]}-{t[2]}")
    if t[3] == 0 and (current_day == last_day_check or t[2] == 1): # Simon: update such that the last day does not necessarily to be the end of a month
        print("  Exporting Tracks for " + Y + " " + MM + ' -- Elapsed Time: ' + str(round(time.perf_counter()-start,2)) + ' seconds')
        start = time.perf_counter() # Reset clock
        ct, ct_inactive = md.splitActiveTracks(ct, cf1)

        # Export inactive tracks
        pd.to_pickle(ct_inactive,trkpath+"/CycloneTracks/"+Y+"/cyclonetracks"+Y+M+".pkl")
        pd.to_pickle(ct,trkpath+"/ActiveTracks/"+Y+"/activetracks"+Y+M+".pkl")

        # Export Cyclone Fields (if new)
        try:
            pd.to_pickle(cflistnew,detpath+"/CycloneFields/CF"+Y+M+".pkl")
        except:
            cflist = []

        if t != endtime_nextmonth:
            # Load netcdf for next month
            ncf = nc.Dataset(inpath_detect+"/"+ra+"_EASE2_N0_"+str(spres)+"km_"+var+"_Hourly_"+str(t[0])+md.dd[t[1]-1]+".nc")
            tlist = ncf['time'][:].data

            # Load Cyclone Fields (if they exist)
            try:
                cflist = pd.read_pickle(detpath+"/CycloneFields/CF"+str(t[0])+md.dd[t[1]-1]+".pkl")
                cftimes = np.array([cf.time for cf in cflist])
            except:
                cflistnew = []

print(f"{current_file} has finished")
