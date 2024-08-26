'''
Author: Simon Tin
Date Created: 26 August 2024
Date Modified: 

Purpose: Preprocessing for the nc file downloaded
'''

'''********************
Import Modules
********************'''
import os 
import xarray as xr
import pandas as pd
import numpy as np
from settings import *

current_file = os.path.basename(__file__)
print(f"Start {current_file}")


os.chdir(inpath)
ds = xr.open_dataset(ncname + ".nc")
try:
    os.chdir(ra)
except:
    os.mkdir(ra)
    os.chdir(ra)
try: os.chdir(var)
except:
    os.mkdir(var)
    os.chdir(var)

'''*******************************************
Preprocessing
*******************************************'''

ds.longitude.values[ds.longitude > 180] -= 360 # rearrange the longitude coordinate fom [0, 360] to [-180, 180] for era5 file
ds = ds.sortby('longitude')
ds = ds.sel(longitude = np.logical_and(ds.longitude >= bbox[2], ds.longitude < bbox[3]))
ds = ds.sel(latitude = np.logical_and(ds.latitude >= bbox[0], ds.latitude < bbox[1]))

enddate = date(yend, mend, dend) + timedelta(days=1)
# ds = ds.isel(time = (pd.to_datetime("2023-01-01") <= pd.Series(ds.time.values)) & (pd.Series(ds.time.values) < pd.to_datetime("2023-02-01")))
ds = ds.isel(time = ("{}-{}-{}".format(ystart, mstart, dstart) <= pd.Series(ds.time.values)) & (pd.Series(ds.time.values) < ("{}-{}-{}".format(enddate.year, enddate.month, enddate.day))))

for month, ds_month in ds.groupby('time.month'):
    timestamps = pd.to_datetime(ds_month.time.values)
    years = np.unique(timestamps.year)[0]
    ds_month.to_netcdf(f'{ra}_{years}{str(month).zfill(2)}.nc') # the filename should starts with variable ra in settings.py