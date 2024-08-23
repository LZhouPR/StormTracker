'''
Author: Alex Crawford
Date Created: 10 Mar 2019
Date Modified: 22 Aug 2019 -- Update for Python 3
               01 Apr 2020 -- Switched output to netCDF instead of GeoTIFF;
                               no longer dependent on gdal module (start V5)
               19 Oct 2020 -- pulled the map creation out of the for loop
               06 Oct 2021 -- added a wrap-around for inputs that prevents
                              empty cells from forming along either 180° or
                              360° longitude (start V6)
               15 Nov 2022 -- replaced "np.int" with "int"

Purpose: Reads in netcdf files & reprojects to the NSIDC EASE2 Grid North.
'''

'''********************
Import Modules
********************'''
import os
import numpy as np
from netCDF4 import Dataset
import xesmf as xe
import CycloneModule_13_2 as md
from settings import *
import pandas as pd

current_file = os.path.basename(__file__)
print(f"Start {current_file}")

'''*******************************************
Main Analysis
*******************************************'''
# Obtain list of nc files:
try:
    os.chdir(outpath_reproj)
except:
    os.mkdir(outpath_reproj)
    os.chdir(outpath_reproj)
fileList = os.listdir(inpath_reproj)
fileList = [f for f in fileList if (f.endswith(ncext) & f.startswith(ra))]

# Identify the time steps:
ref_netcdf = Dataset(inpath_reproj+"/"+fileList[-1])

# Create latitude and longitude arrays:
lons = ref_netcdf.variables['longitude'][:]
lats = ref_netcdf.variables['latitude'][:]

outprjnc = Dataset(suppath_reproj+'/EASE2_N0_'+str(int(xsize/1000))+'km_Projection.nc')
outlat = outprjnc['lat'][:].data
outlon = outprjnc['lon'][:].data

# Close reference netcdf:
ref_netcdf.close()

# Define Grids as Dictionaries
grid_in = {'lon': np.r_[lons,lons[0]], 'lat': lats}
grid_out = {'lon': outlon, 'lat': outlat}

# Create Regridder
regridder = xe.Regridder(grid_in, grid_out, 'bilinear')

print("Step 2. Set up dates of analysis")
years = range(ystart,yend+1)
mos = range(mstart,mend+1)
hrs = [h*timestep for h in range(int(24/timestep))]

ly = md.leapyearBoolean(years) # annual boolean for leap year or not leap year


# Start the reprojection loop
print("Step 3. Load, Reproject, and Save")
for y in years:
    Y = str(y)

    for m in mos:
        M = mons[m-1]

        mlist, hlist = [], []
        ncList = [f for f in fileList if Y+M in f]

        if len(ncList) > 1:
            print("Multiple files with the date " + Y+M + " -- skipping.")
            continue
        if len(ncList) == 0:
            print("No files with the date " + Y+M + " -- skipping.")
        else:
            nc = Dataset(inpath_reproj+"/"+ncList[0])
            tlist = nc.variables[nctvar][:]

            # # Restrict days to those that exist:
            # if m == 2 and ly[y-ymin] == 1 and dmax > dpm[m-1]:
            #     dmax1 = 29
            # elif dmax > dpm[m-1]:
            #     dmax1 = dpm[m-1]
            # else:
            #     dmax1 = dmax

#             # For days that DO exist:
#             for d in range(dmin,dmax1+1):
            days_month = days[(days.month == m) & (days.year == y)]

            for day in days_month:
                y, m, d = day.year, day.month, day.day

                timeD = md.daysBetweenDates(startdate,[y,m,d])*24

                print(" " + Y + " " + M + " " + str(d))
                for h in hrs:

                    # Establish Time
                    timeH = timeD + h
                    # Read from netcdf array
                    inArr = nc.variables[ncvar_p][np.where(tlist == timeH)[0][0],:,:]
                    masked_values = np.ma.masked_equal(nc.variables[ncvar_p][np.where(tlist == timeH)[0][0],:,:], 1e+20)
                    inArr = np.ma.filled(masked_values, np.nan)

                    # Transform data
                    outArr = regridder(np.c_[inArr,inArr[:,0]])
                    outArr[outlat < 0] = np.nan # Limits to Northern Hemisphere

                    # Add to list
                    mlist.append(outArr)
                    hlist.append(timeH)

        # Write monthly data to netcdf file
        ncf = Dataset(ra+"_EASE2_N0_"+str(int(xsize/1000))+"km_"+var+"_Hourly_"+Y+M+".nc", 'w')
        ncf.description = 'Mean sea-level pressure from ERA5. Projection specifications\
        for the EASE2 projection (Lambert Azimuthal Equal Area;\
        lat-origin = 90°N, lon-origin=0°, # cols = ' + str(nx) + ',\
        # rows = ' + str(ny) + ', dx = ' + str(xsize) + ', dy = ' + str(ysize) + ', units = meters'
        ncf.source = 'netCDF4 python module'

        ncf.createDimension('time', len(mlist))
        ncf.createDimension('x', nx)
        ncf.createDimension('y', ny)
        ncft = ncf.createVariable('time', int, ('time',))
        ncfx = ncf.createVariable('x', np.float64, ('x',))
        ncfy = ncf.createVariable('y', np.float64, ('y',))
        ncfArr = ncf.createVariable(ncvar_p, np.float64, ('time','y','x'))

        try:
            ncft.units = nc.variables[nctvar].units
        except:
            ncft.units = f'hours since 1900-01-01 00:00:00.0'

        ncfx.units = 'm'
        ncfy.units = 'm'
        ncfArr.units = 'Pa'

        # For x and y, note that the upper left point is the edge of the grid cell, but
        ## for this we really want the center of the grid cell, hence dividing by 2.
        ncft[:] = np.array(hlist)
        ncfx[:] = np.arange(-xsize*(nx-1)/2, xsize*(nx-1)/2+xsize, xsize)
        ncfy[:] = np.arange(-ysize*(ny-1)/2, ysize*(ny-1)/2+ysize, ysize)
        ncfArr[:] = np.array(mlist)

        ncf.close()

print(f"{current_file} has finished")
