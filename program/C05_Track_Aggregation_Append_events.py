'''
Author: Alex Crawford
Date Created: 10 Mar 2015
Date Modified:  18 Apr 2016; 10 Jul 2019 (update for Python 3);
                10 Sep 2020 (switch from geotiff to netcdf), switch to uniform_filter from scipy.ndimage
                30 Sep 2020 (switch back to slower custom smoother because of what scipy does to NaNs)
                18 Feb 2021 (edited seasonal caluclations to work directly from months, not monthly climatology,
                             allowing for cross-annual averaging)
                13 Sep 2021: If a pre-existing file exists, this script will append new results
                            instead of overwriting for all years. Climatologies no longer in this script.
                01 Nov 2021: Added the possibility of appending prior years as will as subsequent years.
                23 Jan 2023: Adapted to version 13

Purpose: Calculate aggergate statistics (Eulerian and Lagrangian) for either
cyclone tracks or system tracks.

User inputs:
    Path Variables, including the reanalysis (ERA, MERRA, CFSR)
    Track Type (typ): Cyclone or System
    Bounding Box ID (bboxnum): 2-digit character string
    Time Variables: when to start, end, the time step of the data
    Aggregation Parameters (minls, mintl, kSizekm)

Other notes:
    Units for track density are tracks/month/gridcell
    Units for event counts are raw counts (#/month/gridcell)
    Units for counts relative to cyclone obs are ratios (%/gridcell/100)
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
import numpy as np
import netCDF4 as nc
import CycloneModule_13_2 as md
from settings import *
# import pickle5

current_file = os.path.basename(__file__)
print(f"Start {current_file}")

'''*******************************************
Main Analysis
*******************************************'''
# print("Main Analysis")
# Ensure that folders exist to store outputs
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
startyears, endyears = [starttime[0] for i in vNames], [endtime_nextmonth[0] for i in vNames]
firstyears, nextyears = [starttime[0] for i in vNames], [endtime_nextmonth[0] for i in vNames]
for v in varsi:
    name = version+"_AggregationFields_Monthly_"+vNames[v]+".nc"
    if name in priorfiles:
        prior = nc.Dataset(name)
        print(prior)

        nextyears[v] = int(np.ceil(prior['time'][:].max()))
        firstyears[v] = int(np.floor(prior['time'][:].min()))

        if starttime[0] < firstyears[v]: # If the desired time range starts before the prior years...
            if endtime_nextmonth[0] >= firstyears[v]:
                startyears[v], endyears[v] = starttime[0], firstyears[v]
            else:
                raise Exception("There is a gap between the ending year requested ("+str(endtime_nextmonth[0]-1)+") and the first year already aggregated ("+str(firstyears[v])+"). Either increase the ending year or choose a different destination folder.")
        elif endtime_nextmonth[0] > nextyears[v]: # If the desired range ends after the prior years...
            if starttime[0] <= nextyears[v]:
                startyears[v], endyears[v] = nextyears[v], endtime_nextmonth[0]
            else:
                raise Exception("There is a gap between the last year already aggregated ("+str(nextyears[v]-1)+") and the starting year requested ("+str(starttime[0])+"). Either decrease the starting year or choose a different destination folder.")
        else:
            raise Exception("All requested years are already aggregated.")
    else:
        startyears[0], endyears[0] = starttime[0], endtime_nextmonth[0]

# Start at the earliest necessary time for ALL variables of interest
# newstarttime = [np.min(np.array(startyears)[varsi]),1,1,0,0,0]
# newendtime = [np.max(np.array(endyears)[varsi]),1,1,0,0,0]
newstarttime = [ystart,mstart,1,0,0,0]
newendtime = [yend,mend,1,0,0,0]

print("Some years may have already been aggregated.\nAggregating for " + str(newstarttime[0]) + "-" + str(newendtime[0]-1) + ".")

vlists = [ [] for v in vNames]

mt = newstarttime
while mt != newendtime:
    print(mt)
    # Extract date
    Y = str(mt[0])
    MM = months[mt[1]-1]
    M = mons[mt[1]-1]
    print(" " + Y + " - " + MM)

    mtdays = md.daysBetweenDates(dateref,mt,lys=1) # Convert date to days since [1900,1,1,0,0,0]
    mt0 = md.timeAdd(mt,[-i for i in monthstep],lys=1) # Identify time for the previous month

    # Define number of valid times for making %s from counting stats
    if MM == "Feb" and md.leapyearBoolean(mt)[0] == 1:
        n = 29*(24/timestep_list[3])
    else:
        n = dpm[mt[1]-1]*(24/timestep_list[3])

    ### LOAD TRACKS ###
    # Load Cyclone/System Tracks
    cs = pd.read_pickle(inpath_agg+"/"+bboxfull_27+"/"+typ+"Tracks/"+Y+"/"+typ.lower()+"tracks"+Y+M+".pkl")
    # cs = pickle5.load(open(inpath+"/"+bboxnum+"/"+typ+"Tracks/"+Y+"/"+bboxnum+typ.lower()+"tracks"+Y+M+".pkl",'rb'))
    try:
        cs0 = pd.read_pickle(inpath_agg+"/"+bboxfull_27+"/"+typ+"Tracks/"+str(mt0[0])+"/"+bboxfull_27+typ.lower()+"tracks"+str(mt0[0])+mons[mt0[1]-1]+".pkl")
        # cs0 = pickle5.load(open(inpath+"/"+bboxnum+"/"+typ+"Tracks/"+str(mt0[0])+"/"+bboxnum+typ.lower()+"tracks"+str(mt0[0])+mons[mt0[1]-1]+".pkl",'rb'))
    except:
        cs0 = []
    # Load Active tracks
    ct2 = pd.read_pickle(inpath_agg+"/ActiveTracks/"+Y+"/activetracks"+Y+M+".pkl")
    # ct2 = pickle5.load(open(inpath+"/ActiveTracks/"+Y+"/activetracks"+Y+M+".pkl",'rb'))
    if typ == "Cyclone":
        cs2 = ct2
    else:
        try: # Convert active tracks to systems as well
            cs2 = md.cTrack2sTrack(ct2,[],dateref,1)[0]
        except:
            cs2 = []

    ### LIMIT TRACKS & IDS ###
    # Limit to tracks that satisfy minimum lifespan and track length
    trs = [c for c in cs if ((c.lifespan() > minlifespan) and (c.trackLength() >= mintracklength))]
    trs2 = [c for c in cs2 if ((c.lifespan() > minlifespan) and (c.trackLength() >= mintracklength))]
    trs0 = [c for c in cs0 if ((c.lifespan() > minlifespan) and (c.trackLength() >= mintracklength))]

    ### CALCULATE FIELDS ###
    fields0 = [np.nan]
    fields1 = md.aggregateEvents([trs,trs0,trs2],typ,mtdays,lats.shape)

    fields = fields0 + fields1

    ### SMOOTH FIELDS ###
    for v in varsi:
        # varFieldsm = md.smoothField(fields[v],kSize) # Smooth
        varFieldsm = ndimage.uniform_filter(fields[v],kSize,mode="nearest") # --> This cannot handle NaNs
        vlists[v].append(varFieldsm) # append to list

    # Increment Month
    mt[2] = 1
    mt = md.timeAdd(mt,monthstep,lys=1)
    

### SAVE FILE ###
print("Step 3. Write to NetCDF")
for v in varsi:
    print(vNames[v])
    mnc = nc.Dataset(version+"_AggregationFields_Monthly_"+vNames[v]+"_NEW.nc",'w')
    mnc.createDimension('y', lats.shape[0])
    mnc.createDimension('x', lats.shape[1])
    mnc.createDimension('time', (max(nextyears[v],newendtime[0])-min(firstyears[v],newstarttime[0]))*12)
    mnc.description = 'Aggregation of cyclone track ' + vNames[v] + ' on monthly time scale.'

    ncy = mnc.createVariable('y', np.float32, ('y',))
    ncx = mnc.createVariable('x', np.float32, ('x',))
    ncy.units, ncx.units = 'm', 'm'
    ncy[:] = np.arange(proj['lat'].shape[0]*spres*1000/-2 + (spres*1000/2),proj['lat'].shape[0]*spres*1000/2, spres*1000)
    ncx[:] = np.arange(proj['lat'].shape[1]*spres*1000/-2 + (spres*1000/2),proj['lat'].shape[1]*spres*1000/2, spres*1000)

    # Add times, lats, and lons
    nctime = mnc.createVariable('time', np.float32, ('time',))
    nctime.units = 'years'
    nctime[:] = np.arange(min(firstyears[v],newstarttime[0]),max(nextyears[v],newendtime[0]),1/12)

    nclon = mnc.createVariable('lon', np.float32, ('y','x'))
    nclon.units = 'degrees'
    nclon[:] = proj['lon'][:]

    nclat = mnc.createVariable('lat', np.float32, ('y','x'))
    nclat.units = 'degrees'
    nclat[:] = proj['lat'][:]

    vout = np.array(vlists[v])
    vout = np.where(vout == 0,np.nan,vout)
    ncvar_height = mnc.createVariable(vNames[v], np.float64, ('time','y','x'))
    ncvar_height.units = vunits[v] + ' -- Smoothing:' + str(kSizekm) + ' km'

    name = version+"_AggregationFields_Monthly_"+vNames[v]+".nc"
    if name in priorfiles: # Append data if prior data existed...
        if vout.shape[0] > 0: # ...and there is new data to be added
            prior = nc.Dataset(name)

            if (newstarttime[0] <= firstyears[v]) and (newendtime[0] >= nextyears[v]): # If the new data starts before and ends after prior data
                ncvar_height[:] = vout

            elif (newstarttime[0] > firstyears[v]) and (newendtime[0] < nextyears[v]): # If the new data starts after and ends before prior data
                ncvar_height[:] = np.concatenate( ( prior[vNames[v]][prior['time'][:].data < newstarttime[0],:,:].data , vout , prior[vNames[v]][prior['time'][:].data >= newendtime[0],:,:].data ) )

            elif (newendtime[0] <= firstyears[v]): # If the new data starts and ends before the prior data
                ncvar_height[:] = np.concatenate( ( vout , prior[vNames[v]][prior['time'][:].data >= newendtime[0],:,:].data ) )

            elif (newstarttime[0] >= nextyears[v]): # If the new data starts and ends after the prior data
                ncvar_height[:] = np.concatenate( ( prior[vNames[v]][prior['time'][:].data < newstarttime[0],:,:].data , vout ) )

            else:
                mnc.close()
                raise Exception("Times are misaligned. Requested Time Range: " + str(starttime) + "-" + str(endtime_nextmonth) + ". Processed Time Range: " + str(newstarttime) + "-" + str(newendtime) + ".")

            prior.close(), mnc.close()
            os.remove(name) # Remove old file
            os.rename(version+"_AggregationFields_Monthly_"+vNames[v]+"_NEW.nc", name) # rename new file to standard name

    else: # Create new data if no prior data existed
        ncvar_height[:] = vout
        mnc.close()
        os.rename(version+"_AggregationFields_Monthly_"+vNames[v]+"_NEW.nc", name) # rename new file to standard name

if (newendtime[0] < endtime_nextmonth[0]) & (max(nextyears) < endtime_nextmonth[0]):
    print("Completed aggregating " + str(newstarttime[0]) + "-" + str(newendtime[0]-1)+".\nRe-run this script to aggregate any time after " + str(max(nextyears[v],newendtime[0])-1) + ".")
else:
    print("Completed aggregating " + str(newstarttime[0]) + "-" + str(newendtime[0]-1)+".")

print(f"{current_file} has finished")
