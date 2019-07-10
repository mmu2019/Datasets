import numpy as np
from netCDF4 import Dataset
import pylab as plt
import os
import time
from units import unit
from units.predefined import define_units
import datetime
from mpl_toolkits.basemap import Basemap
import netCDF4 as nc
import math
from urllib.request import urlretrieve
from scipy.interpolate import griddata
from subroutines import *

# set up Data directory
DataDir = "/Users/mingquan/newDATA"

# set up initial and final years of data
start_yr = 1979
end_yr   = 2018

# Set general information for the data source
remote_source = "https://www.esrl.noaa.gov/psd/data/gridded/data.cmap.html#detail"
gist_source = "https://github.com/mmu2019/Datasets/blob/master/read-pr-cmap.py"
local_source = 'precip.mon.mean.nc'
stamp1 = '2019-04-30'

datestr = str(datetime.datetime.now())
TmpStr  = datestr.split(' ')
stamp2  = TmpStr[0]

# set up institutions where created the original dataset
sourceID = "CMAP.v1904"
instit1  = "NOAA/OAR/ESRL PSD, Boulder, Colorado, USA"

# Create temporal dimension
nyears = end_yr - start_yr + 1
nmonth = 12
smonth = np.asarray(['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'])
month_bnd = np.asarray([0,31,59,90,120,151,181,212,243,273,304,334,365],dtype=float)
tbnd  = np.asarray([((np.arange(nyears)*365)[:,np.newaxis]+month_bnd[:-1]).flatten(),
                    ((np.arange(nyears)*365)[:,np.newaxis]+month_bnd[+1:]).flatten()]).T
tbnd += (start_yr - 1850)*365
tbnd.shape
t     = tbnd.mean(axis=1)

# set up the temporal and spatial resolutions for the original and final dataset
period = str(start_yr) + "-01 through " + str(end_yr) + "-12"
origtr = "monthly"
origsr = "2.5 degree"
origut = "mm/day"
finltr = "monthly"
finlsr = "0.5 degree"
finlut = "kg/m2/s"

# Create new spatial dimension
res    = 0.5
latbnd = np.asarray([np.arange(- 90    , 90     ,res),
                     np.arange(- 90+res, 90+0.01,res)]).T
lonbnd = np.asarray([np.arange(-180    ,180     ,res),
                     np.arange(-180+res,180+0.01,res)]).T
lat    = latbnd.mean(axis=1)
lon    = lonbnd.mean(axis=1)

# Create some fake data
data   = np.ma.masked_array(np.random.rand(t.size,lat.size,lon.size))

data[:,:,:] = 0.0

ntim = t.size
nlat = lat.size
nlon = lon.size

# read single netCDF file
filename = DataDir + '/' + sourceID + '/' + local_source
print(filename)
gpcp=Dataset(filename,'r',format='NETCDF3')
print(gpcp) 
print(gpcp.variables) 

time1= gpcp.variables['time']
lat1 = gpcp.variables['lat']
lon1 = gpcp.variables['lon']
pr1  = gpcp.variables['precip']
long_name = pr1.long_name
original_unit = pr1.units

ntim1 = time1.size
nlat1 = lat1.size
nlon1 = lon1.size

# convert lon1 from 0-360 to 180W-180E
if lon1[0]>0:
   nlon12=int(nlon1/2)
   lon0   = np.ma.masked_array(np.random.rand(nlon1), type=lon1.dtype)
   pr0    = np.ma.masked_array(np.random.rand(ntim1,nlat1,nlon1), dtype=pr1.dtype)

   print(pr1.dtype)
   print(pr0.dtype)
   for i in range(nlon12):
       lon0[i]           = lon1[i] - 180.0
       lon0[i+nlon12]    = lon1[i+nlon12] - 180.0
       pr0[:,:,i]        = pr1[:,:,i+nlon12]
       pr0[:,:,i+nlon12] = pr1[:,:,i]

   del pr1
   del lon1
   lon1                  = lon0
   pr1                   = pr0
   del pr0
   del lon0

# cut data in the required period
pr   = np.ma.masked_array(np.random.rand(t.size,lat1.size,lon1.size))
pr[:,:,:]  = pr1[0:t.size,:,:]

# convert unit from mm/day to kg/m2/s
pr[:,:,:] = pr[:,:,:]/(24.0*3600.0)

ij = 0
for i in range(nyears):
    year    = i + start_yr
    for j in range(nmonth):
        temp         = pr[ij,:,:]
        data[ij,:,:] = NearestNeighborInterpolation(lat1,lon1,temp,lat,lon)
        ij = ij + 1

data_min = data.min()
data_max = data.max()

with Dataset(DataDir + "/pr.nc", mode="w") as dset:

    # Create netCDF dimensions
    dset.createDimension("time",size=  t.size)
    dset.createDimension("lat" ,size=lat.size)
    dset.createDimension("lon" ,size=lon.size)
    dset.createDimension("nb"  ,size=2       )

    # Create netCDF variables
    T  = dset.createVariable("time"       ,t.dtype   ,("time"     ))
    TB = dset.createVariable("time_bounds",t.dtype   ,("time","nb"))
    X  = dset.createVariable("lat"        ,lat.dtype ,("lat"      ))
    XB = dset.createVariable("lat_bounds" ,lat.dtype ,("lat","nb" ))
    Y  = dset.createVariable("lon"        ,lon.dtype ,("lon"      ))
    YB = dset.createVariable("lon_bounds" ,lon.dtype ,("lon","nb" ))
    D  = dset.createVariable("pr"         ,data.dtype,("time","lat","lon"), fill_value = -999.)

    print(D.shape)

    # Load data and encode attributes
    # time
    T [...]    = t
    T.units    = "days since 1850-01-01"
    T.calendar = "noleap"
    T.bounds   = "time_bounds"
    TB[...]    = tbnd
    T.standard_name = "time"
    T.long_name     = "time"

    # lat
    X [...]    = lat
    X.units    = "degrees_north"
    XB[...]    = latbnd
    X.standard_name = "latitude"
    X.long_name     = "latitude"

    # lon
    Y [...]    = lon
    Y.units    = "degrees_east"
    YB[...]    = lonbnd
    Y.standard_name = "longitude"
    Y.long_name     = "longitude"

    # data
    D[...] = data
    D.units = "kg m-2 s-1"
    D.standard_name = "precipitation"
    D.long_name     = long_name
    D.actual_range = np.asarray([data_min,data_max])

    dset.title = "CPC Merged Analysis of Precipitation (excludes NCEP Reanalysis)"
    dset.version = "1904"
    dset.institutions = "%s" % (instit1)
    dset.source = " Merged Precipitation with gauge observations, a variety of satellite observations and the NCEP–NCAR reanalysis"
    dset.history = """
%s: downloaded source from %s;
%s: converted to netCDF with %s""" % (stamp1, remote_source, stamp2, gist_source)
    dset.references  = """
@ARTICLE{Xie1997,
  author = {Xie, P., and P.A. Arkin},
  title = {Global precipitation: A 17-year monthly analysis based on gauge observations, satellite estimates, and numerical model outputs},
  journal = {Bull. Amer. Meteor. Soc.},
  year = {1997},
  number = {78},
  page = {2539-2558},
  doi = {https://doi.org/10.1175/1520-0477(1997)078<2539:GPAYMA>2.0.CO;2}
}"""
    dset.comments = """
time_period: %s; original_temporal_resolution: %s; original_spatial_resolution: %s; original_units: %s; final_temporal_resolution: %s; final_spatial_resolution: %s; final_units: %s""" % (period, origtr, origsr, origut, finltr, finlsr, finlut)
    dset.convention = "CF-1.7"
