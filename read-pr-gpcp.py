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
remote_source = "http://www.esrl.noaa.gov/psd/data/gridded/data.gpcp.html"
gist_source = "https://gist.github.com/nocollier/d73585731756fa472731065389af45dc"
local_source = 'GPCP2.3/precip.mon.mean.nc'
stamp1 = '2019-04-18'

datestr = str(datetime.datetime.now())
TmpStr  = datestr.split(' ')
stamp2  = TmpStr[0]

print(datestr)
print(stamp2)

instit1 = "NOAA/OAR/ESRL PSD, Boulder, Colorado, USA"
instit2 = "ESSIC/CICS, University of Maryland, USA"

period = "1979-01 through 2018-12"
origtr = "monthly"
origsr = "2.5 degree"
origut = "mm/day"
finltr = "monthly"
finlsr = "0.5 degree"
finlut = "kg/m2/s"

# Create temporal dimension
nyears    = 40
nmonth    = 12
smonth = np.asarray(['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'])
month_bnd = np.asarray([0,31,59,90,120,151,181,212,243,273,304,334,365],dtype=float)
tbnd  = np.asarray([((np.arange(nyears)*365)[:,np.newaxis]+month_bnd[:-1]).flatten(),
                    ((np.arange(nyears)*365)[:,np.newaxis]+month_bnd[+1:]).flatten()]).T
tbnd += (1979-1850)*365
tbnd.shape
t     = tbnd.mean(axis=1)

# Create old spatial dimension
res    = 2.5
lat1bnd = np.asarray([np.arange(- 90    , 90     ,res),
                     np.arange(- 90+res, 90+0.01,res)]).T
lon1bnd = np.asarray([np.arange(-180    ,180     ,res),
                     np.arange(-180+res,180+0.01,res)]).T
lat1    = lat1bnd.mean(axis=1)
lon1    = lon1bnd.mean(axis=1)

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
area   = np.ma.masked_array(np.random.rand(lat.size,lon.size))

data[:,:,:] = 0.0
area[:,:]   = 0.0

ntim = t.size
nlat = lat.size
nlon = lon.size

nlat1 = lat1.size
nlon1 = lon1.size

def rebin(arr, new_shape):
    """Rebin 2D array arr to shape new_shape by averaging."""
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)

# read single netCDF file
filename = 'GPCP2.3/precip.mon.mean.nc'
print(filename)
gpcp=Dataset(filename,'r',format='NETCDF3')
print(gpcp) 
print(gpcp.variables) 

pr1 = gpcp.variables['precip']
long_name = pr1.long_name

print(t.shape)
print(pr1.shape)
print(data.shape)

pr   = np.ma.masked_array(np.random.rand(t.size,lat1.size,lon1.size))
    
# convert unit from mm/day to kg/m2/s
pr[:,:,0:72]   = pr1[0:ntim,:,72:144]/(24.0*3600.0)
pr[:,:,72:144] = pr1[0:ntim,:,0:72]/(24.0*3600.0)

print("shape of pr")
print(pr.shape)
print("size of pr")
print(pr.size)

ij = 0
for i in range(nyears):

    year    = i + 1979

    print(year)

    for j in range(nmonth):

        temp         = pr[ij,:,:]
        print(temp.shape)
        data[ij,:,:] = rebin(temp, (nlat, nlon))
        area[:,:]    = rebin(area2, (nlat, nlon)) 
        area[:,:]    = area[:,:]*4

        print(ij)
        ij = ij + 1

data_min = data.min()
data_max = data.max()

with Dataset("pr.nc", mode="w") as dset:

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
    D  = dset.createVariable("pr"        ,data.dtype,("time","lat","lon"), fill_value = -999.)

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

    dset.title = "GPCP Version 2.3 Combined Precipitation Dataset (Final)"
    dset.version = "2.3"
    dset.institutions = "%s; %s" % (instit1, instit2)
    dset.source = "GPCC gauge analysis to bias correct satellite estimates over land and merge with satellite based on sampling"
    dset.history = """
%s: downloaded source from %s;
%s: converted to netCDF with %s""" % (stamp1, remote_source, stamp2, gist_source)
    dset.references  = """
@ARTICLE{Adler2003,
  author = {Adler, R.F., G.J. Huffman, A. Chang, R. Ferraro, P. Xie, J. Janowiak, B. Rudolf, U. Schneider, S. Curtis, D. Bolvin, A. Gruber, J. Susskind, P. Arkin},
  title = {The Version 2 Global Precipitation Climatology Project (GPCP) Monthly Precipitation Analysis (1979 - Present)},
  journal = {J. Hydrometeor.},
  year = {2003},
  number = {4(6)},
  page = {1147-1167},
  doi = {https://doi.org/10.1175/1525-7541(2003)004<1147:TVGPCP>2.0.CO;2}
}"""
    dset.comments = """
time_period: %s; original_temporal_resolution: %s; original_spatial_resolution: %s; original_units: %s; final_temporal_resolution: %s; final_spatial_resolution: %s; final_units: %s""" % (period, origtr, origsr, origut, finltr, finlsr, finlut)
    dset.convention = "CF-1.7"
