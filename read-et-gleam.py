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
remote_source = "https://www.gleam.eu/#datasets"
gist_source = "https://gist.github.com/nocollier/d73585731756fa472731065389af45dc"
local_source = 'GLEAM3.2a/E_yyyy_GLEAM_v3.2a.nc'
stamp1 = '2019-05-01'

datestr = str(datetime.datetime.now())
TmpStr  = datestr.split(' ')
stamp2  = TmpStr[0]

print(datestr)
print(stamp2)

instit1 = "Laboratory of Hydrology and Water Management, Ghent University, Belgium"
instit2 = "Department of Earth Sciences, Vrije Universiteit Amsterdam, the Netherlands"

period = "1980-01 through 2017-12"
origtr = "daily"
origsr = "0.25 degree"
origut = "mm/day"
finltr = "monthly"
finlsr = "0.5 degree"
finlut = "kg/m2/s"

# Read desert data from MODIS 
q=Dataset('/Users/mingquan/DATA/biomes/MODIS/derived/desert_0.5x0.5.nc','r',format='NETCDF4')
q.variables 

desert = q.variables['desert']
desert.shape

#july_temp = q.variables['ts']
#jan_july = np.concatenate((may_temp, jun_temp), axis=0)
#jan_july.shape

#aver_temp = np.mean(jan_temp, axis=0)# average temperature

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
t.shape

# Create old spatial dimension
res    = 0.25
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

nlat = lat.size
nlon = lon.size

#create mesh
#Yold, Xold   = np.meshgrid(lat1, lon1)
#Ynew, Xnew   = np.meshgrid(lat,lon)

def rebin(arr, new_shape):
    """Rebin 2D array arr to shape new_shape by averaging."""
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)

ij = 0
for i in range(nyears):

    year    = i + 1980

    print(year)

    # read single netCDF file
    filename = 'GLEAM3.2a/E_'  + str(year) + '_GLEAM_v3.2a.nc'
    print(filename)
    gleam=Dataset(filename,'r',format='NETCDF4')
    print(gleam) 
    print(gleam.variables) 

    et1 = gleam.variables['E']
    et = et1[:,::-1,:]
    
    # convert unit from mm/day to kg/m2/s
    et[:,:,:] = et[:,:,:]/(24.0*3600.0)

    for j in range(nmonth):

        # Variable from multiple files.
        print(smonth[j])

        long_name = gleam.variables['ET'].long_name
        ;long_name = 'evaportranspiration'
        print(long_name)

        data[ij,:,:] = rebin(abgb, (nlat, nlon)) 

        print(ij)
        ij = ij + 1
        #print(i1)
        #print(j1)

        print('total area burned final')
        print(np.sum(area[:,:]))
        print(np.sum(data[ij-1,:,:]*area[:,:]))

    print(j)
print(i)

data_min = data.min()
data_max = data.max()

# Calculate climatology of burned area
mdata     = data.mean(axis=0)

# Set all grid water, snow, ice or desert cells where have never burned to missing (-999)
#for i in range(len(t)):
#    desert      = np.where(mdata>0, 0, desert)
#    data[i,:,:] = np.where(desert[:,:]>=70, -999, data[i,:,:])

with Dataset("et.nc", mode="w") as dset:

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
    D  = dset.createVariable("burntArea"        ,data.dtype,("time","lat","lon"), fill_value = -999.)

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
    D.standard_name = "evaportranspiration"
    D.long_name     = long_name
    D.actual_range = np.asarray([data_min,data_max])
    
    dset.title = "GFED version 4.1 burned area fraction with small fires (GFED4.1s)"
    dset.version = "4.1s"
    dset.institutions = "%s; %s; %s" % (instit1, instit2, instit3)
    dset.source = "Satellite (MODIS Terra and Aqua) derived product"
    #dset.source = "Global Fire Emissions Database, Version 4.1 (GFED4.1s) Monthly and daily fire emissions 1997-present"
    dset.history = """
%s: downloaded source from %s;
%s: converted to netCDF with %s""" % (stamp1, remote_source, stamp2, gist_source)
    dset.references  = """
@ARTICLE{Martens2017,
  author = {Martens, B., D.G. Miralles, H. Lievens, R. van der Schalie, R.A.M. de Jeu, D. Fernández-Prieto, H.E. Beck, W.A. Dorigo, and N.E.C. Verhoest},
  title = {GLEAM v3: satellite-based land evaporation and root-zone soil moisture},
  journal = {Geosci. Model Dev.},
  year = {2017},
  number = {10},
  page = {1903-1925},
  doi = {https://doi.org/10.5194/gmd-10-1903-2017}
}
@ARTICLE{Miralles2011,
  author = {Miralles, D.G., T.T.H. Holmes, R.A.M. De Jeu, J.H. Gash, A.G.C.A. Meesters, and A.J. Dolman},
  title = {Global land-surface evaporation estimated from satellite-based observations},
  journal = {Hydrol. Earth Syst. Sci.},
  year = {2011},
  number = {15},
  page = {453-469},
  doi = {https://doi.org/10.5194/hess-15-453-2011}
}"""
    dset.comments = """
time_period: %s; original_temporal_resolution: %s; original_spatial_resolution: %s; original_units: %s; final_temporal_resolution: %s; final_spatial_resolution: %s; final_units: %s""" % (period, origtr, origsr, origut, finltr, finlsr, finlut)
    dset.convention = "CF-1.7"
