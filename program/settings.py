from datetime import date, timedelta
import numpy as np
import pandas as pd
'''*******************************************
nc file variables you need to change
*******************************************'''
ncname = "era5_2023"
ncvar_p = "msl" # pressure variable in the nc file
nctvar = "time" # time dimension
ncext = '.nc' # file format

timestep = 3 # in hours

# Time Variables
ystart, yend = 2023, 2023
mstart, mend = 7, 10
dstart, dend = 1, 31

# spatial resolution of the gridded data
# xsize, ysize = 25000, -25000 # in meters

inpath = "E:/Simon/tracking/alexcrawford0927/final_version/Cressida"

ra = "ERA5" # put your file within this folder
# ras = [f"ERA5_em{i}" for i in range(1, 51+1)] # for ensemble members
# ra = ras[0]

var = "SLP"

# mask the area of interest in C03
# not necessary to change
bbox = [0, 60, 100, 170] # minlat, maxlat (exclusive), minlon, maxlon (exclusive)

# set the time range and interval for plotting the results
# if not set manually, all the results will be plotted
plot_time_range = ['2023-01-01', '2023-01-15']
plot_interval = "24H"
'''*******************************************
Please adjust the parameter settings for storm detections and trackings
*******************************************'''

"""
Detection Parameters
"""
minfield = 80000 # minimum reasonable value in field array
maxfield = 200000 # maximum reasonable value in field array

# Size of kernel (km) used to determine if a grid cell is a local minimum
##  Starting in Version 12.1, users enter a distance in km instead of the kernel
## size. This way, the kSize can adapt to different spatial resolutions. The
## equivalent of a 3 by 3 kernel with 100 km resolution would be 100
## i.e., kSize = (2*kSizekm/spres)+1
kSizekm = 200

# Maximum fraction of neighboring grid cells with no data (Not a Number) allowed
### for a grid cell to be considered during minimum detection
nanThresh = 0.4

# minimum slp gradient for identifying (and eliminating) weak minima:
d_slp = 750 # slp difference in Pa (use 0 to turn off)
d_dist = 1000000 # distance in m (units that match units of cellsize)

# maximum elevation for masking out high elevation minima
maxelev = 1500. # elevation in m (use 10000 to turn off)

# minimum latitude for masking out the Equator (default should be 5 for NH and -5 for SH)
minlat = 5

# Contour interval (Pa; determines the interval needed to identify closed
### contours,and therefore cyclone area)
contint = 200

# Multi-center cyclone (mcc) tolerance is the maximum ratio permitted between the
### number of unshared and total contours in a multi-centered cyclone. "Unshared"
### contours are only used by the primary center. "Shared" contours are used
### by both the primary and secondary centers.
mcctol = 0.5 # (use 0 to turn off mcc's; higher makes mcc's more likely)
# Multi-center cyclone (mcc) distance is the maximum distance (in m) two minima can
### lie apart and still be considered part of the same cyclone system
mccdist = 1200000

# Tracking Parameters #
# Maximum speed is the fastest that a cyclone center is allowed to travel; given
### in units of km/h. To be realistic, the number should be between 100 and 200.
### and probably over 125 (based on Rudeva et al. 2014). To turn off, set to
### np.inf. Also, note that instabilities occur at temporal resolution of 1-hr.
### Tracking at 6-hr and a maxspeed of 125 km/hr is more comprable to tracking
### at 1-hr and a maxspeed of 300 km/hr (assuming spatial resolution of 50 km).
maxspeed = 150 # constant value
# maxspeed = 150*(3*math.log(timestep[3],6)+2)/timestep[3] # One example of scaling by temporal resolution

# The reduction parameter is a scaling of cyclone speed.  When tracking, the
### algorithm uses the cyclone speed and direction from the last interval to
### estimate a "best guess" position. This parameter modifies that speed, making
### it slower. This reflects how cyclones tend to slow down as they mature. To
### turn off, set to 1.
red = 0.75

# Regenesis Paramater
rg = 1
# 0 = regenesis starts a new system track;
# 1 = regenesis continues previous system track with new ptid

"""
Track criteria
"""
minlifespan = 1 # minimum lifespan (in  days) for a track to be considered
mintracklength = 1000 # minimum track length (in km)
mindisplacement = 500 # in km
mindisp = 0 # minimum surface pressure
minlat = 5 # minimum latitude

"""
Aggregation
"""
bboxmin = -1*np.inf #None #    Set to None if using values instead #
bboxmax = 500 #None #   Set to None if using values instead #
values = [1] #[26] # [10,11,12,13,15] #  If using a min and max, set to 1, otherwise, these are the
## acccepted values within the mask raster (e.g., the regions of interest)

'''*******************************************
Variables that you do not need to change
*******************************************'''
# defined by Simon
days = pd.date_range(f"{ystart}-{mstart}-{dstart}", f"{yend}-{mend}-{dend}")

nx, ny = 720, 720 # number of grid cells in Projections/EASE2_N0_{res}km_GenesisRegions.nc
spres = 25 # Spatial resolution (in km), assume 1 degree equals 100 km for both latitude and longitude
xsize, ysize = spres * 1000, spres * 1000

# Time Variables
starttime = [ystart, mstart, dstart, 0, 0, 0] # Format: [Y,M,D,H,M,S]
enddate = date(yend, mend, dend) + timedelta(days=1)
yend, mend, dend = enddate.year, enddate.month, enddate.day
endtime = [yend, mend, dend, 0, 0, 0]# stop BEFORE this time (exclusive)
timestep_list = [0,0,0, timestep,0,0] # Time step in [Y,M,D,H,M,S]

months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep","Oct","Nov","Dec"]
mons = ["01","02","03","04","05","06","07","08","09","10","11","12"]
seasons = np.array([1,2,3,4,5,6,7,8,9,10,11,12]) # Ending month for three-month seasons (e.g., 2 = DJF)
dpm = [31,28,31,30,31,30,31,31,30,31,30,31] # days per month (non leap year)
monthstep = [0,1,0,0,0,0] # A Time step that increases by 1 month [Y,M,D,H,M,S]
startdate = [1900,1,1] # The starting date for the reanalysis time steps
# C17_ExportToCSV
years = np.arange(ystart, yend+1)
mos = np.arange(mstart, mend)

bboxnum_10 = "10"
bboxmain = "" # The main bbox your subsetting from; usually "" for "all cyclones", otherwise BBox##
bboxnum_27 = "27"

bboxfull_27 = "BBox" + bboxnum_27

"""
Path Variables
"""
inpath_reproj = inpath+"/"+ra+"/"+var #
outpath_reproj = inpath + "/" +ra+"/"+var+"_EASE2_N0_"+str(int(xsize/1000))+"km" # file path for reprojection nc file
suppath_reproj = inpath+"/Projections"
version = "13_2" # Detection Version

# vert = '' # Original: Ptest # Tracking Version
inpath_detect = inpath + "/" + ra + "/SLP_EASE2_N0_"+str(spres) + "km"
outpath_detect = inpath + "/CycloneTracking_" + ra
suppath_detect = inpath + "/Projections/EASE2_N0_" + str(spres) + "km_Projection.nc"
# system detection
# subset = "" # use "" if performing on all cyclones
regpath = inpath+"/Projections/EASE2_N0_25km_GenesisRegions.nc"
inpath_system = f"{outpath_detect}/tracking13_2"
inpath_agg = inpath_system 
# outpath_agg = inpath_agg + "/" + bboxfull_27
outpath_agg = inpath_agg

dateref = [startdate[0],startdate[1],startdate[2],0,0,0]  # [Y,M,D,H,M,S]
reftime = starttime
prior = 0 # 1 = a cyclone track object exists for a prior month; 0 = otherwise

# C04 
typ = "System" # Cyclone, System, or Active
ncvar_height = 'z' # 'reg' # 'z' # Set to None unless file is a netcdf

# C05_CycloneStatSummary_AllStorms
ext = ".tif"
kind1 = "System" # System, Cyclone
kind = kind1 + "Tracks" # Can be AFZ, Arctic, or other region (or no region), followed by System or Cyclone

rg = 1 # Whether regenesis of a cyclone counts as a track split (0) or track continuation (1)

V = "_GenReg" # An optional version name; suggested to start with "_" or "-" to separate from years in file title

# C05 Aggregation
# Variables
vNames = ["countA","gen","lys","spl","mrg"]
varsi = range(1,len(vNames)) # range(0,1) #
vunits = ['ratio','count','count','count','count']
agg = [-1,-1,-1,-1,-1]

# C06B
# subset2 = '' # '_DeepeningDsqP' + ''

# C99
plot_directory = 'fig'



