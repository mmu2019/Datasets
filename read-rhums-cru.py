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

# set up Data directory
DataDir = "/Users/mingquan/projects"

# Set general information for the data source
remote_source = "https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.02/"
gist_source = "https://github.com/mmu2019/Datasets/blob/master/read-rhums-cru.py"
local_source = DataDir + '/CRU/v4.02/cru_ts4.02.1901.2017.vap.dat.nc'
stamp1 = '2019-05-06'

datestr = str(datetime.datetime.now())
TmpStr  = datestr.split(' ')
stamp2  = TmpStr[0]

print(datestr)
print(stamp2)

instit = "Climatic Research Unit, School of Environmental Sciences, University of East Anglia, UK"

period = "1980-01 through 2017-12"
origtr = "monthly"
origsr = "0.5 degree"
origut = "hPa"
finltr = "monthly"
finlsr = "0.5 degree"
finlut = "%"

# Create temporal dimension
nyears    = 38
nmonth    = 12
smonth = np.asarray(['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'])
month_bnd = np.asarray([0,31,59,90,120,151,181,212,243,273,304,334,365],dtype=float)
tbnd  = np.asarray([((np.arange(nyears)*365)[:,np.newaxis]+month_bnd[:-1]).flatten(),
                    ((np.arange(nyears)*365)[:,np.newaxis]+month_bnd[+1:]).flatten()]).T
tbnd += (1980-1850)*365
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
filename1 = DataDir + '/CRU/v4.02/cru_ts4.02.1901.2017.tmp.dat.nc'
filename2 = DataDir + '/CRU/v4.02/cru_ts4.02.1901.2017.vap.dat.nc'
print(filename1)
print(filename2)
cru1=Dataset(filename1,'r',format='NETCDF3')
cru2=Dataset(filename2,'r',format='NETCDF3')
print(cru1) 
print(cru2) 
print(cru1.variables) 
print(cru2.variables) 

tas  = cru1.variables['tmp']
vap  = cru2.variables['vap']
#long_name = tas.long_name
data1 = tas[ntot-ntim:ntot,:,:]
data2 = vap[ntot-ntim:ntot,:,:]
    
# estimate saturated vapor pressure by using air temperature
Ew = 6.112*np.exp(17.62*data1[:,:,:]/(243.12+data1[:,:,:]))

# estimate relative humidity by vapor pressure and saturated vapor pressure
data = data2[:,:,:]/Ew[:,:,:]

# convert unit from fractions to percentage
data[:,:,:] = data[:,:,:]*100

# set data at all ridcells greater than 100% to be 95%
data[:,:,:] = np.where(data[:,:,:]>=100, 95, data[:,:,:])

print(t.shape)
print(tas.shape)
print(data.shape)

data_min = data.min()
data_max = data.max()

with Dataset(DataDir + "/rhums.nc", mode="w") as dset:

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
    D  = dset.createVariable("rhums"      ,data.dtype,("time","lat","lon"), fill_value = -999.)

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
    D.units = "%"
    D.standard_name = "2-meter relative humidity"
    D.long_name     = D.standard_name
    D.actual_range = np.asarray([data_min,data_max])
    
    dset.title   = "CRU time series (TS) high-resolution gridded datasets"
    dset.version = "4.02"
    dset.institutions = instit
    dset.source  = "monthly observations at meteorological stations across the world land areas"
    dset.history = """
%s: downloaded source from %s;
%s: estimate saturated vapor pressure using air temperature;
%s: calculate relative humidity using saturated vapor pressure and vapor pressure;
%s: converted to ILAMB required netCDF with %s""" % (stamp1, remote_source, stamp2, stamp2, stamp2, gist_source)
    dset.references  = """
@ARTICLE{Harris2014,
  author = {Harris, I., P.D. Jones, T.J. Osborn and D.H. Lister},
  title = {Updated high-resolution grids of monthly climatic observations - the CRU TS3.10 Dataset},
  journal = {International Journal of Climatology},
  year = {2014},
  number = {34(3)},
  page = {623-642},
  doi = {https://doi.org/10.1002/joc.3711}
}
@ARTICLE{WMO2008,
  author = {WMO CIMO Guide},
  title = {ANNEX 4.B. FORMULAE FOR THE COMPUTATION OF MEASURES OF HUMIDITY},
  journal = {MEASUREMENT OF HUMIDITY},
  year = {2008},
  number = {Chapter 4},
  page = {163},
  doi = {https://library.wmo.int/doc_num.php?explnum_id=3151}
}"""
    dset.comments = """
time_period: %s; original_temporal_resolution: %s; original_spatial_resolution: %s; original_units: %s; final_temporal_resolution: %s; final_spatial_resolution: %s; final_units: %s""" % (period, origtr, origsr, origut, finltr, finlsr, finlut)
    dset.convention = "CF-1.7"
