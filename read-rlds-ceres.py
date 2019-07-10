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
start_yr = 2000
end_yr   = 2018

VarID         = "rlds"
RawVarID      = "sfc_lw_down_all_mon"
standard_name = "surface downwelling longwave radiation"

# Set general information for the data source
remote_source = "https://ceres.larc.nasa.gov/products-info.php?product=EBAF"
gist_source   = "https://github.com/mmu2019/Datasets/blob/master/read-llds-ceres.py"
local_source  = 'CERES_EBAF_Ed4.1_Subset_200003-201809.nc'
stamp1        = '2019-06-25'

datestr = str(datetime.datetime.now())
TmpStr  = datestr.split(' ')
stamp2  = TmpStr[0]

# set up institutions where created the original dataset
sourceID = "CERES.ed4.1"
instit1  = "NASA Langley Research Center"

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
origsr = "1 degree"
origut = "W m-2"
finltr = "monthly"
finlsr = "0.5 degree"
finlut = "W m-2"

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

data[:,:,:] = -999

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
pr1  = gpcp.variables[RawVarID]
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
pr[:,:,:] = -999.
pr[2:2+ntim1,:,:] = pr1[:,:,:]

ij = 0
for i in range(nyears):
    year    = i + start_yr
    for j in range(nmonth):
        temp         = pr[ij,:,:]
        data[ij,:,:] = NearestNeighborInterpolation(lat1,lon1,temp,lat,lon)
        ij = ij + 1

data_min = data.min()
data_max = data.max()

with Dataset(DataDir + "/rlds.nc", mode="w") as dset:

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
    D  = dset.createVariable(VarID         ,data.dtype,("time","lat","lon"), fill_value = -999.)

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
    D.units = "W m-2"
    D.standard_name = standard_name
    D.long_name     = long_name
    D.actual_range = np.asarray([data_min,data_max])

    dset.title = "CERES EBAF TOA and Surface Fluxes"
    dset.version = "Ed4.1"
    dset.institutions = "%s" % (instit1)
    dset.source = "Monhtly mean surface fluxes calculated by a radiative transfer model and constrained by the combined Terra and Aqua SSF1deg measurements"
    dset.history = """
%s: downloaded source from %s;
%s: converted to netCDF with %s""" % (stamp1, remote_source, stamp2, gist_source)
    dset.references  = """
@ARTICLE{Loeb2018,
  author = {Loeb, N.G., D.R. Doelling, H. Wang, W. Su, C. Nguyen, J.G. Corbett, L. Liang, C. Mitrescu, F.G. Rose, and S. Kato},
  title = {Clouds and the Earth's Radiant Energy System (CERES) Energy Balanced and Filled (EBAF) Top-of-Atmosphere (TOA) Edition-4.0 Data Product},
  journal = {Journal of Climate},
  year = {2018},
  number = {31(2)},
  page = {895-918},
  doi = {https://doi.org/10.1175/JCLI-D-17-0208.1}
}
@ARTICLE{Kato2018,
  author = {Kato, S., F. G. Rose, D. A. Rutan, T. E. Thorsen, N. G. Loeb, D. R. Doelling, X. Huang, W. L. Smith, W. Su, and S.-H. Ham},
  title = {Surface irradiances of Edition 4.0 Clouds and the Earth's Radiant Energy System (CERES) Energy Balanced and Filled (EBAF) data product},
  journal = {Journal of Climate},
  year = {2018},
  number = {31},
  page = {4501-4527},
  doi = {https://doi.org/10.1175/JCLI-D-17-0523.1}
}"""
    dset.comments = """
time_period: %s; original_temporal_resolution: %s; original_spatial_resolution: %s; original_units: %s; final_temporal_resolution: %s; final_spatial_resolution: %s; final_units: %s""" % (period, origtr, origsr, origut, finltr, finlsr, finlut)
    dset.convention = "CF-1.7"
