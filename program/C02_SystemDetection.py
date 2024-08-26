'''
Author: Alex Crawford
Date Created: 11 Jan 2016
Date Modified: 8 Dec 2017, 4 Jun 2019 (Python 3), 13 Jun 2019 (warning added)
Purpose: Convert a series of center tracks to system tracks. Warning: If a) you
wish to re-run this process on some of the data and b) you are using rg = 1
(allowing regeneration), you need to re-run from the reftime or accept that
some active storms at the re-start point will get truncated.

User inputs:
    Path Variables
    Bounding Box ID (subsetnum): 2-digit character string
    Time Variables: when to start, end, the time step of the data
    Regenesis Paramter: 0 or 1, depending on whether regenesis continues tracks
'''

'''********************
Import Modules
********************'''

# Import clock:
from time import perf_counter
# Start script stopwatch. The clock starts running when time is imported
start = perf_counter()
import pandas as pd
import CycloneModule_13_2 as md
from settings import *
import os 

current_file = os.path.basename(__file__)
print(f"Start {current_file}")
'''*******************************************
Set up Environment
*******************************************'''
# print("Setting up environment.")

'''*******************************************
Main Analysis
*******************************************'''
# print("Main Analysis")

mt = starttime
# current_day = pd.to_datetime(f"{mt[0]}-{mt[1]}-{mt[2]}")
while mt != endtime_nextmonth:

# while mt[0] <= endtime[0] and mt[1] <= endtime[1]: # Simon changed
    # Extract date
    Y = str(mt[0])
    M = mons[mt[1]-1]
    print (" " + Y + " - " + M)

    # Load Cyclone Tracks
    ct = pd.read_pickle(inpath_system+"/CycloneTracks/"+Y+"/"+"cyclonetracks"+Y+M+".pkl")
    # Create System Tracks
    if mt == reftime:
        cs, cs0 = md.cTrack2sTrack(ct,[],dateref,rg)
        pd.to_pickle(cs,inpath_system+"/SystemTracks/"+Y+"/"+"systemtracks"+Y+M+".pkl")

    else:
        # Extract date for previous month
        mt0 = md.timeAdd(mt,[-d for d in monthstep])
        Y0 = str(mt0[0])
        M0 = mons[mt0[1]-1]

        # Load previous month's system tracks
        cs0 = pd.read_pickle(inpath_system+"/SystemTracks/"+Y0+"/"+"systemtracks"+Y0+M0+".pkl")

        # Create system tracks
        cs, cs0 = md.cTrack2sTrack(ct,cs0,dateref,rg)
        pd.to_pickle(cs,inpath_system+"/SystemTracks/"+Y+"/"+"systemtracks"+Y+M+".pkl")
        pd.to_pickle(cs0,inpath_system+"/SystemTracks/"+Y0+"/"+"systemtracks"+Y0+M0+".pkl")

    # Increment Time Step
    mt = md.timeAdd(mt,monthstep)
    mt[2] = 1 # Simon: set that the day of the following month to be the last day 


# print('Elapsed time:',round(perf_counter()-start,2),'seconds')
print(f"{current_file} has finished")
