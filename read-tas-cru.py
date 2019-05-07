import numpy as np
from netCDF4 import Dataset
import pylab as plt
import os
import time
import datetime
from mpl_toolkits.basemap import Basemap
import netCDF4 as nc
import math
from urllib.request import urlretrieve
from scipy.interpolate import griddata

# Set general information for the data source
remote_source = "https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.02/"
gist_source = "https://gist.github.com/nocollier/d73585731756fa472731065389af45dc"
local_source = 'CRU/v4.02/cru_ts4.02.1901.2017.tmp.dat.nc'
stamp1 = '2019-05-06'

datestr = str(datetime.datetime.now())
TmpStr  = datestr.split(' ')
stamp2  = TmpStr[0]

print(datestr)
print(stamp2)

instit = "Climatic Research Unit, School of Environmental Sciences, University of East Anglia, UK"

period = "1960-01 through 2017-12"
origtr = "monthly"
origsr = "0.5 degree"
origut = "C"
finltr = "monthly"
finlsr = "0.5 degree"
finlut = "K"

# Create temporal dimension
nyears    = 58
nmonth    = 12
smonth = np.asarray(['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'])
month_bnd = np.asarray([0,31,59,90,120,151,181,212,243,273,304,334,365],dtype=float)
tbnd  = np.asarray([((np.arange(nyears)*365)[:,np.newaxis]+month_bnd[:-1]).flatten(),
                    ((np.arange(nyears)*365)[:,np.newaxis]+month_bnd[+1:]).flatten()]).T
tbnd += (1960-1850)*365
tbnd.shape
t     = tbnd.mean(axis=1)

# Create spatial dimension
res    = 0.5
latbnd = np.asarray([np.arange(- 90    , 90     ,res),
                     np.arange(- 90+res, 90+0.01,res)]).T
lonbnd = np.asarray([np.arange(-180    ,180     ,res),
                     np.arange(-180+res,180+0.01,res)]).T
lat    = latbnd.mean(axis=1)
lon    = lonbnd.mean(axis=1)

ntot = (2017-1900)*12
ntim = t.size
nlat = lat.size
nlon = lon.size

# read single netCDF file
filename = 'CRU/v4.02/cru_ts4.02.1901.2017.tmp.dat.nc'
print(filename)
cru=Dataset(filename,'r',format='NETCDF3')
print(cru) 
print(cru.variables) 

tas  = cru.variables['tmp']
long_name = tas.long_name
data = tas[ntot-ntim:ntot,:,:]
    
# convert unit from C to K
data[:,:,:] = data[:,:,:] + 273.14

print(t.shape)
print(tas.shape)
print(data.shape)

data_min = data.min()
data_max = data.max()

with Dataset("tas.nc", mode="w") as dset:

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
    D  = dset.createVariable("tas"        ,data.dtype,("time","lat","lon"), fill_value = -999.)

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
    D.units = "K"
    D.standard_name = "2-meter air temperature"
    D.long_name     = long_name
    D.actual_range = np.asarray([data_min,data_max])
    
    dset.title = "CRU time series (TS) high-resolution gridded datasets"
    dset.version = "4.02"
    dset.institutions = instit
    dset.source = "monthly observations at meteorological stations across the world land areas"
    dset.history = """
%s: downloaded source from %s;
%s: converted to netCDF with %s""" % (stamp1, remote_source, stamp2, gist_source)
    dset.references  = """
@ARTICLE{Harris2014,
  author = {Harris, I., P.D. Jones, T.J. Osborn and D.H. Lister},
  title = {Updated high-resolution grids of monthly climatic observations - the CRU TS3.10 Dataset},
  journal = {International Journal of Climatology},
  year = {2014},
  number = {34(3)},
  page = {623-642},
  doi = { https://doi.org/10.1002/joc.3711}
}"""
    dset.comments = """
time_period: %s; original_temporal_resolution: %s; original_spatial_resolution: %s; original_units: %s; final_temporal_resolution: %s; final_spatial_resolution: %s; final_units: %s""" % (period, origtr, origsr, origut, finltr, finlsr, finlut)
    dset.convention = "CF-1.7"
