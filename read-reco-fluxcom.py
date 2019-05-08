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
remote_source = "https://doi.org/doi:10.17871/FLUXCOM_RS_METEO_CRUNCEPv6_1980_2013_v1"
gist_source = "https://github.com/mmu2019/Datasets/blob/master/read-reco-fluxcom.py"
local_source = DataDir + '/FluxCom/reco/TER.ANN.CRUNCEPv6.monthly.YYYY.nc'
stamp1 = '2019-05-07'

datestr = str(datetime.datetime.now())
TmpStr  = datestr.split(' ')
stamp2  = TmpStr[0]

print(datestr)
print(stamp2)

instit = "Department Biogeochemical Integration, Max Planck Institute for Biogeochemistry, Germany"

period = "1980-01 through 2013-12"
origtr = "monthly"
origsr = "0.5 degree"
origut = "g/m2/day"
finltr = "monthly"
finlsr = "0.5 degree"
finlut = "kg/m2/s"

# Create temporal dimension
nyears    = 34
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

ntim = t.size
nlat = lat.size
nlon = lon.size

# Create some fake data
data   = np.ma.masked_array(np.random.rand(t.size,lat.size,lon.size))
area   = np.ma.masked_array(np.random.rand(lat.size,lon.size))

data[:,:,:] = 0.0
area[:,:]   = 0.0

nlat = lat.size
nlon = lon.size

ij = 0
for i in range(nyears):

    year    = i + 1980

    print(year)

    # read single netCDF file
    filename = DataDir + '/FluxCom/reco/TER.ANN.CRUNCEPv6.monthly.' + str(year) + '.nc'
    print(filename)
    flx=Dataset(filename,'r',format='NETCDF3')
    data0 = flx.variables['TER']
    lats  = flx.variables['lat']

    #long_name = data0.long_name
    long_name = "terrestrial ecosystem respiration"

    #data1 = np.where(data0[:,:,:]<=-999, 0, data0[:,:,:])

    #latrange = data0.latitude_range

    #data2 = np.float_(data1[:,:,:])*data0.DataScaleFactor + data0.DataOffsetValue

    #ilat1 = np.where(lat==latrange[0])
    #ilat2 = np.where(lat==latrange[1])

    #j1 = int(ilat1[0])
    #j2 = int(ilat2[0]) + 1

    for j in range(nmonth):

        data[ij,:,:] = data0[j,::-1,:]

        ij = ij + 1

# convert unit from g/m2/day to kg/m2/s
data[:,:,:] = data[:,:,:]/(24*3600*1000)

print(t.shape)
print(data0.shape)
print(data.shape)

data_min = data.min()
data_max = data.max()

with Dataset(DataDir + "/reco.nc", mode="w") as dset:

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
    D  = dset.createVariable("reco"        ,data.dtype,("time","lat","lon"), fill_value = -999.)

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
    D.units = "kg/m2/s"
    D.standard_name = long_name
    D.long_name     = long_name
    D.actual_range = np.asarray([data_min,data_max])
    
    dset.title = "FLUXCOM (RS+METEO) Global Land Carbon Fluxes using CRUNCEP climate data"
    dset.version = "1"
    dset.institutions = instit
    dset.source = "Data generated by Artificial Neural Networks and forced with CRUNCEPv6 meteorological data and MODIS (RS+METEO)"
    dset.history = """
%s: downloaded source from %s;
%s: converted to netCDF with %s""" % (stamp1, remote_source, stamp2, gist_source)
    dset.references  = """
@ARTICLE{Jung2017,
  author = {Jung, M., M. Reichstein, C.R. Schwalm, C. Huntingford, S. Sitch, A. Ahlstrom, A. Arneth, G. Camps-Valls, P. Ciais, P. Friedlingstein, F. Gans, K. Ichii, A.K. Jain, E. Kato, D. Papale, B. Poulter, B. Raduly, C. Rodenbeck, G. Tramontana, N. Viovy, Y.P. Wang, U. Weber, S. Zaehle and N. Zeng},
  title = {Compensatory water effects link yearly global land CO2 sink changes to temperature},
  journal = {Nature},
  year = {2017},
  number = {541},
  page = {516-520},
  doi = {https://doi.org/10.1038/nature20780}
}
@ARTICLE{Tramontana2016,
  author = {Tramontana, G., M. Jung, C.R. Schwalm, K. Ichii, G. Camps-Valls, B. Raduly, M. Reichstein, M.A. Arain, A. Cescatti, G. Kiely, L. Merbold, P. Serrano-Ortiz, S. Sickert, S. Wolf, and D. Papale},
  title = {Predicting carbon dioxide and energy fluxes across global FLUXNET sites with regression algorithms},
  journal = {Biogeosciences},
  year = {2016},
  number = {13},
  page = {4291-4313},
  doi = {https://doi.org/10.5194/bg-13-4291-2016}
}"""
    dset.comments = """
time_period: %s; original_temporal_resolution: %s; original_spatial_resolution: %s; original_units: %s; final_temporal_resolution: %s; final_spatial_resolution: %s; final_units: %s""" % (period, origtr, origsr, origut, finltr, finlsr, finlut)
    dset.convention = "CF-1.7"
