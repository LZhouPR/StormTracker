'''
Author: Alex Crawford
Date Created: 10/28/16
Date Modified: 5/6/20 (Updated for Python 3)
    18 Jan 2021 --> Built in removal of DpDr and precip columns if present
    11 Apr 2022 --> Fixed a bug in the removal of DpDr and precip columns and unit converdion
    23 Jan 2023 --> Changed order of columns, added year/month/day/hour columns, changed units to be
        hPa instead of Pa
Purpose: Export Cyclone objects' data table as CSV files with year, month, and TID in the filename.
'''

'''********************
Import Modules
********************'''
print("Loading modules.")
import pandas as pd
import numpy as np
import CycloneModule_13_2 as md
import os
from settings import *
# import pickle5

current_file = os.path.basename(__file__)
print(f"Start {current_file}")

'''*******************************************
Define Variables
*******************************************'''
# print("Defining variables")

version = "13_2"
bbox = "" # Use BBox## for subsets
# typ = "System"
# path = "E:/Simon/tracking/alexcrawford0927/final_version/Cressida/CycloneTracking/tracking"+version

# # Time Variables
# years = range(1999,1999+1)
# mos = range(1,12+1)
# dateref = [1900,1,1,0,0,0]

# # Cyclone Parameters
# # Aggregation Parameters
# minls = 1 # minimum lifespan (in days) for a track to be considered
# mintl = 1000 # in km for version 11 and beyond
# minlat = 0 # minimum latitude that must be acheived at some point in the track

'''*******************************************
Main Analysis
*******************************************'''
# print("Main Analysis")

# Load parameters
params = pd.read_pickle(inpath_system+"/cycloneparams.pkl")
# params = pickle5.load(open(path+"/cycloneparams.pkl",'rb'))
# try:
#     spres = params['spres']
# except:
#     spres = 100

# Set Up Output Directories
try:
    os.chdir(inpath_system+"/"+bbox+"/CSV"+typ)
except:
    os.mkdir(inpath_system+"/"+bbox+"/CSV"+typ)

for y in years:
    Y = str(y)
    #Create directories
    try:
        os.chdir(inpath_system+"/"+bbox+"/CSV"+typ+"/"+Y)
    except:
        os.mkdir(inpath_system+"/"+bbox+"/CSV"+typ+"/"+Y)
        os.chdir(inpath_system+"/"+bbox+"/CSV"+typ+"/"+Y)

    for m in mos:
        M = md.dd[m-1]

        #Create directories
        try:
            os.chdir(inpath_system+"/"+bbox+"/CSV"+typ+"/"+Y+"/"+M)
        except:
            os.mkdir(inpath_system+"/"+bbox+"/CSV"+typ+"/"+Y+"/"+M)
            os.chdir(inpath_system+"/"+bbox+"/CSV"+typ+"/"+Y+"/"+M)

# Write CSV for Systems
if typ == "System":
    for y in years:
        Y = str(y)
        print(Y)
        for m in mos:
            M = md.dd[m-1]

            # Load data
            trs = pd.read_pickle(inpath_system+"/"+bbox+"/"+typ+"Tracks"+"/"+Y+"/"+bbox+typ.lower()+"tracks"+Y+M+".pkl")
            # trs = pickle5.load(open(path+"/"+bbox+"/"+typ+"Tracks"+"/"+Y+"/"+bbox+typ.lower()+"tracks"+Y+M+".pkl",'rb'))

            # Write CSVs
            for tr in trs:
                if tr.lifespan() >= minlifespan and tr.trackLength() >= mintracklength and np.max(tr.data.lat) >= minlat:
                    trdata = tr.data
                    if len(np.intersect1d(trdata.columns,['precip','precipArea','DpDr'])) == 3:
                        trdata = trdata.drop(columns=['precip','precipArea','DpDr'])

                    dates = [md.timeAdd(dateref,[0,0,t,0,0,0]) for t in trdata['time']]

                    trdata['year'] = np.array([date[0] for date in dates])
                    trdata['month'] = np.array([date[1] for date in dates])
                    trdata['day'] = np.array([date[2] for date in dates])
                    trdata['hour'] = np.array([date[3] for date in dates])
                    trdata['p_cent'] = trdata['p_cent'] / 100 # Pa --> hPa
                    trdata['p_edge'] = trdata['p_edge'] / 100 # Pa --> hPa
                    trdata['depth'] = trdata['depth'] / 100 # Pa --> hPa
                    trdata['radius'] = trdata['radius'] * spres # to units of km
                    trdata['area'] = trdata['area'] * spres * spres # to units of km^2
                    trdata['DsqP'] = trdata['DsqP'] / spres / spres * 100 * 100 / 100 # to units of hPa/[100 km]^2
                    trdata['p_grad'] = trdata['p_grad'] / 100 * 1000 * 1000 # to units of hPa / [1000 km]
                    trdata['Dp'] = trdata['Dp'] / 100 # Pa --> hPa
                    trdata['DpDt'] = trdata['DpDt'] / 100 # Pa/day --> hPa/day

                    trdata.to_csv(inpath_system+"/"+bbox+"/CSV"+typ+"/"+Y+"/"+M+"/"+typ+version+bbox+"_"+Y+M+"_"+str(tr.sid)+".csv",index=False,columns=list(trdata.columns[-4:])+list(trdata.columns[:-4]))

else: print("Please review and change the coding below !!!")
# # Write CSV for Cyclones
# else:
#     for y in years:
#         Y = str(y)
#         print(Y)
#         for m in mos:
#             M = md.dd[m-1]
#             # trs = pickle5.load(open(path+"/"+bbox+"/"+typ+"Tracks"+"/"+Y+"/"+bbox+typ.lower()+"tracks"+Y+M+".pkl",'rb'))
#             trs = pd.read_pickle(path+"/"+bbox+"/"+typ+"Tracks"+"/"+Y+"/"+bbox+typ.lower()+"tracks"+Y+M+".pkl")

#             for tr in trs:
#                 if tr.lifespan() >= minls and tr.trackLength() >= mintl and np.max(tr.data.lat) >= minlat:
#                     trdata= tr.data
#                     if len(np.intersect1d(trdata.columns,['precip','precipArea','DpDr'])) == 3:
#                         trdata = trdata.drop(columns=['precip','precipArea','DpDr'])
                    
#                     trdata['year'] = np.array([date[0] for date in dates])
#                     trdata['month'] = np.array([date[1] for date in dates])
#                     trdata['day'] = np.array([date[2] for date in dates])
#                     trdata['hour'] = np.array([date[3] for date in dates])
#                     trdata['p_cent'] = trdata['p_cent'] / 100 # Pa --> hPa
#                     trdata['p_edge'] = trdata['p_edge'] / 100 # Pa --> hPa
#                     trdata['depth'] = trdata['depth'] / 100 # Pa --> hPa
#                     trdata['radius'] = trdata['radius'] * spres # to units of km
#                     trdata['area'] = trdata['area'] * spres * spres # to units of km^2
#                     trdata['DsqP'] = trdata['DsqP'] / spres / spres * 100 * 100 / 100 # to units of hPa/[100 km]^2
#                     trdata['p_grad'] = trdata['p_grad'] / 100 * 1000 * 1000 # to units of hPa / [1000 km]
#                     trdata['Dp'] = trdata['Dp'] / 100 # Pa --> hPa
#                     trdata['DpDt'] = trdata['DpDt'] / 100 # Pa/day --> hPa/day

#                     trdata.to_csv(path+"/"+bbox+"/CSV"+typ+"/"+Y+"/"+M+"/"+typ+version+bbox+"_"+Y+M+"_"+str(tr.tid)+".csv",index=False,columns=list(trdata.columns))

print(f"{current_file} has finished")