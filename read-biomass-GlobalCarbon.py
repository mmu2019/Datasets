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

remote_source = "https://science.jpl.nasa.gov/people/Saatchi"
gist_source = "https://gist.github.com/nocollier/d73585731756fa472731065389af45dc"
local_source = DataDir + '/GlobalCarbon_0.25x0.25/maxent_agb_mean_global_v5.1_masked_latlon_0.25deg.nc'
stamp1 = '2017-10-31'
stamp2 = '2019-03-14'

datestr = str(datetime.datetime.now())
TmpStr  = datestr.split(' ')
stamp3  = TmpStr[0]

period = "climatology"
origtr = "climatology"
origsr = "0.25 degree"
origut = "kg/m2 biomass"
finltr = "climatology"
finlsr = "0.5 degree"
finlut = "kg/m2 carbon"

# Create temporal dimension
nyears    = 1
nmonth    = 12
smonth = np.asarray(['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'])
month_bnd = np.asarray([0,31,59,90,120,151,181,212,243,273,304,334,365],dtype=float)
tbnd1  = np.asarray([((np.arange(nyears)*365)[:,np.newaxis]+month_bnd[:-1]).flatten(),
                    ((np.arange(nyears)*365)[:,np.newaxis]+month_bnd[+1:]).flatten()]).T
tbnd1 += (2000-1850)*365
#t1     = tbnd1.mean(axis=1)

tbnd    = tbnd1[0,:]
tbnd[0] = tbnd1[0,0]
tbnd[1] = tbnd1[11,1]
t       = tbnd.mean(axis=0)

print(tbnd)
print(t)

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

    year    = i + 2000

    print(year)

    # read single netCDF file
    filename1 = DataDir + '/GlobalCarbon_0.25x0.25/maxent_agb_mean_global_v5.1_masked_latlon_0.25deg.nc'
    filename2 = DataDir + '/GlobalCarbon_0.25x0.25/maxent_bgb_mean_global_v5.1_masked_latlon_0.25deg.nc'
    filename3 = DataDir + '/GlobalCarbon_0.25x0.25/maxent_landratio_global_v5.1_masked_latlon_0.25deg.nc'
    print(filename1)
    print(filename2)
    print(filename3)
    ag=Dataset(filename1,'r',format='NETCDF4')
    bg=Dataset(filename2,'r',format='NETCDF4')
    ld=Dataset(filename3,'r',format='NETCDF4')
    print(ag) 
    print(bg) 
    print(ag.variables) 
    print(bg.variables) 

    for j in range(1):

        # Variable from multiple files.
        print(smonth[j])

        #print(ba.variables.keys())
        agb = ag.variables['Band1']
        bgb = bg.variables['Band1']
        lnd = ld.variables['Band1']

        agb1 = agb[:,:]*lnd[:,:]
        bgb1 = bgb[:,:]*lnd[:,:]

        abgb = agb1[:,:] + bgb1[:,:]
        abgb[:,:] = abgb[:,:]/20.

        print(agb1.shape)
        print(bgb1.shape)
        print(abgb.shape)

        #print(agb.variables['Band1'].long_name)
        #long_name = agb.variables['Band1'].long_name
        long_name = 'global above and below ground live biomass carbon'

        data[ij,:,:] = rebin(abgb, (nlat, nlon)) 
        #area[:,:]    = rebin(area2, (nlat, nlon)) 
        #area[:,:]    = area[:,:]*4

        #for j1 in range(nlat):
        #    for i1 in range(nlon):
        #        #interp
        #        data[ij,j1,i1] = np.mean(agb1[j1:j1+1,i1:i1+1]+bgb1[j1:j1+1,i1:i1+1])/2
        #
        ij = ij + 1

with Dataset(DataDir + "/biomass.nc", mode="w") as dset:

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
    D  = dset.createVariable("biomass"        ,data.dtype,("time","lat","lon"), fill_value = -999.)

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
    D.units = "kg m-2"
    D.standard_name = "global live biomass carbon"
    D.long_name     = long_name
    D.actual_range = np.asarray([data.min(),data.max()])
    
    dset.title       = "Global forest live biomass carbon version 5.1"
    dset.version     = "5.1"
    dset.institutions= "Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA, USA"
    dset.source      = "A combination of data in situ inventory plots and satellite light detection and ranging (Lidar) samples"
    dset.history     = """
%s: downloaded source from %s;
%s: convert original tif files to netcdf data using gdal_translate;
%s: converted to ILAMB required netCDF with %s""" % (stamp1, remote_source, stamp2, stamp3, gist_source)
    dset.references  = """
@ARTICLE{Saatchi2011,
  author = {Saatchi, S.S., N.L. Harris, S. Brown, M. Lefsky, E.T. Mitchard, W. Salas, B.R. Zutta, W. Buermann, S.L. Lewis, S. Hagen, S. Petrova, L. White, M. Silman, and A. Morel},
  title = {Benchmark map of forest carbon stocks in tropical regions across three continents},
  journal = {Proceedings of the National Academy of Sciences},
  year = {2011},
  number = {108(24)},
  page = {9899-9904},
  doi = {https://doi.org/10.1073/pnas.1019576108}
}"""
    dset.comments = """
time_period: %s; original_temporal_resolution: %s; original_spatial_resolution: %s; original_units: %s; final_temporal_resolution: %s; final_spatial_resolution: %s; final_units: %s""" % (period, origtr, origsr, origut, finltr, finlsr, finlut)
    dset.convention = "CF-1.7"
