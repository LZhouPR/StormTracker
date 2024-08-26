'''
Author: Simon Tin
Date Created: 26 August 2024
Date Modified: 

Purpose: Plot the tracking results
'''

'''********************
Import Modules
********************'''
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr 
from glob import glob
import numpy as np
import cartopy.crs as ccrs
import imageio
from settings import * 
import os 

os.chdir(inpath)
if not os.path.exists(plot_directory):
    os.mkdir(plot_directory)

filename = ncname + '.nc'
tracking_path = f'{inpath_agg}/CSVSystem/*/*/*.csv'

projection =ccrs.AzimuthalEquidistant(central_longitude=(bbox[2]+bbox[3])/2, central_latitude=(bbox[0]+bbox[1])/2)
ds = xr.open_dataset(filename)

'''*******************************************
CSV processing
*******************************************'''
print("Combining CSV files")
file_list = glob(tracking_path)
dfs = [pd.read_csv(file) for file in file_list]

dfs = []
for id, file in enumerate(file_list):
    df = pd.read_csv(file)
    df['storm_id'] = id + 1
    dfs.append(df)

df_tracks = pd.concat(dfs, ignore_index=True)
df_tracks['time'] = df_tracks[['year', 'month', 'day', 'hour']].apply(lambda x: '{}-{:02d}-{:02d} {:02d}'.format(x[0], x[1], x[2], x[3]), axis=1)
df_tracks['time'] = pd.to_datetime(df_tracks['time'], format = '%Y-%m-%d %H')
df_tracks.to_csv(f'{inpath_agg}/CSVSystem/combined_results.csv', index=False)


'''*******************************************
Plot results
*******************************************'''
print("Plotting results")
# ds_era_fil = ds_era.sel(time = ds_era['time.month']== 9)
ds_fil = ds.sel(longitude = (ds.longitude >= bbox[2]) & (ds.longitude < bbox[3]), latitude = (ds.latitude >= bbox[0]) & (ds.latitude < bbox[1]))
ds_fil.longitude.values[ds_fil.longitude > 180] -= 360

lon = ds_fil.longitude.values
lat = ds_fil.latitude.values
lon, lat = np.meshgrid(lon, lat)
extent = bbox[2:] + bbox[:2]

max_p, min_p = np.round(ds_fil[ncvar_p].values.max(), 2), np.round(ds_fil[ncvar_p].values.min(), 2)
intv = 100
try: time_range = pd.date_range(plot_time_range[0], plot_time_range[1], freq=plot_interval)
except: time_range = pd.date_range(f"{ystart}-{mstart}-{dstart}", f"{yend}-{mend}-{dend}", freq=str(timestep) + "h", inclusive='left')
for time in time_range:
# for time in pd.date_range('2023-08-25', '2023-09-05', freq='12H'):
    print(time)
    df_current_point = df_tracks.query("time == @time")
    current_id = df_current_point['storm_id'].unique()
    fig, ax = plt.subplots(1, 1, figsize=(10, 6), subplot_kw=dict(projection=projection), layout='constrained')
    ax.coastlines()
    ax.set_extent(extent)
    contourf = ax.contourf(lon, lat, ds_fil[ncvar_p].sel(time=time), transform=ccrs.PlateCarree(), cmap='jet_r', levels = np.arange(min_p, max_p+intv, intv))
    plt.colorbar(contourf, label="Surface pressure (Pa)", extend='both')
    ax.scatter(df_current_point.lon.values, df_current_point.lat.values, c=df_current_point.p_cent.values, transform=ccrs.PlateCarree(), cmap = 'jet_r')
    for i, row in df_current_point.iterrows():
        lon_p = row['lon']
        lat_p = row['lat']
        p_cent = row['p_cent']
        if ~np.isnan(p_cent):
            ax.text(lon_p+0.5, lat_p+0.2, round(p_cent/100), transform=ccrs.PlateCarree())
    time_str = pd.to_datetime(time).strftime("%Y-%m-%d %H")
    ax.set_title(time_str)
    # time_str = time.strftime('%Y-%m-%d %H')
    for id in current_id:
        df_current_track = df_tracks.query("storm_id == @id and time <= @time")
        lon_t, lat_t = df_current_track.lon.values, df_current_track.lat.values
        ax.plot(lon_t, lat_t, transform=ccrs.PlateCarree())
    fig.savefig(f"fig/{time_str}.jpg")
    plt.close(fig)  # Close the figure to free up memory

plots = glob(f'{plot_directory}/*.jpg')
# Create a GIF from the list of images
images = [imageio.imread(plot) for plot in plots]
imageio.mimsave(f'{plot_directory}/animation.gif', images, duration=500, loop=0)

# # validate ERA5 dataset
# era5_path = 'E:/Simon/tracking/alexcrawford0927/final_version/ERA5_2023.nc'
# extent = [100, 150, 5, 40]
# ds_era = xr.open_dataset(era5_path)
# ds_era['longitude'] = ds_era.longitude - 180
# ds_era = ds_era.sortby('latitude')
# lon = ds_era.longitude.values
# lat = ds_era.latitude.values
# lon, lat = np.meshgrid(lon, lat)
# ds_era = ds_era.sel(time = ds_era['time.month']== 9)

# for time in ds_era.time.values[:1]:
#     print(time)
#     df_current_point = df_tracks.query("time == @time")
#     current_id = df_current_point['storm_id'].unique()
#     fig, ax = plt.subplots(1, 1, figsize=(10, 6), subplot_kw=dict(projection=projection), layout='constrained')
#     ax.coastlines()
#     ax.set_extent(extent)
#     contourf = ax.contourf(lon, lat, ds_era.msl.sel(time=time), transform=ccrs.PlateCarree(), cmap='jet_r', levels = np.arange(98000, 102000+100, 100))
#     plt.colorbar(contourf, label="Surface pressure (Pa)", extend='both')
#     # ax.scatter(df_current_point.lon.values, df_current_point.lat.values, c=df_current_point.p_cent.values, transform=ccrs.PlateCarree(), cmap = 'jet_r')
#     # for i, row in df_current_point.iterrows():
#     #     lon_p = row['lon']
#     #     lat_p = row['lat']
#     #     p_cent = row['p_cent']
#     #     if ~np.isnan(p_cent):
#     #         ax.text(lon_p+0.5, lat_p+0.2, round(p_cent/100), transform=ccrs.PlateCarree())
#     time_str = pd.to_datetime(time).strftime("%Y-%m-%d %H")
#     ax.set_title(time_str)