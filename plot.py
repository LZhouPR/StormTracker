import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr 
from glob import glob
import numpy as np
import cartopy.crs as ccrs
import imageio

year = 2023
era5_path = 'E:/Simon/tracking/alexcrawford0927/final_version/ERA5_2023.nc'
tracking_path = f'E:/Simon/tracking/alexcrawford0927/final_version/Cressida/Cyclonetracking/tracking13_2/CSVSystem/{year}/*/*.csv'
extent = [100, 150, 5, 40]

projection =ccrs.AzimuthalEquidistant(central_longitude=125, central_latitude=20)
# ds_era = xr.open_dataset('E:/Simon/spatial analysis/tracking/era5/ERA5_1999_europe.nc')
ds_era = xr.open_dataset(era5_path)
ds_era['longitude'] = ds_era.longitude - 180

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

# ds_era_fil = ds_era.sel(time = ds_era['time.month']== 9)
ds_era_fil = ds_era.sel(longitude = (ds_era.longitude >= extent[0]) & (ds_era.longitude < extent[1]), latitude = (ds_era.latitude >= extent[2]) & (ds_era.latitude < extent[3]))
# ds_era_fil = ds_era.copy()

lon = ds_era_fil.longitude.values
lat = ds_era_fil.latitude.values
lon, lat = np.meshgrid(lon, lat)


# for time in ds_era_fil.time.values[:1]:
for time in pd.date_range('2023-08-25', '2023-09-05', freq='12H'):
    print(time)
    df_current_point = df_tracks.query("time == @time")
    current_id = df_current_point['storm_id'].unique()
    fig, ax = plt.subplots(1, 1, figsize=(10, 6), subplot_kw=dict(projection=projection), layout='constrained')
    ax.coastlines()
    ax.set_extent(extent)
    contourf = ax.contourf(lon, lat, ds_era_fil.msl.sel(time=time), transform=ccrs.PlateCarree(), cmap='jet_r', levels = np.arange(98000, 102000+100, 100))
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
    time_str = time.strftime('%Y-%m-%d %H')
    for id in current_id:
        df_current_track = df_tracks.query("storm_id == @id and time <= @time")
        lon_t, lat_t = df_current_track.lon.values, df_current_track.lat.values
        ax.plot(lon_t, lat_t, transform=ccrs.PlateCarree())
    fig.savefig(f"fig/{time_str}.jpg")
    plt.close(fig)  # Close the figure to free up memory


# plots = glob('fig/*.jpg')
# # Create a GIF from the list of images
# images = [imageio.imread(plot) for plot in plots]
# imageio.mimsave('animation.gif', images, duration=500, loop=0)

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