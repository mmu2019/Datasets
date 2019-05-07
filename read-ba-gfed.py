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
remote_source = "https://www.geo.vu.nl/~gwerf/GFED/GFED4/"
gist_source = "https://gist.github.com/nocollier/d73585731756fa472731065389af45dc"
local_source = DataDir + '/GFED4S/GFED4.1s_yyyy.hdf5'
stamp1 = '2019-03-13'

datestr = str(datetime.datetime.now())
TmpStr  = datestr.split(' ')
stamp2  = TmpStr[0]

print(datestr)
print(stamp2)

instit1 = "VU University Amsterdam, Faculty of Earth and Life Sciences, Earth and Climate Cluster, Netherlands"
instit2 = "Goddard Space Flight Center, USA"
instit3 = "University of California Irvine, Department of Earth System Science, USA"

period = "1997-01 through 2016-12"
origtr = "monthly"
origsr = "0.25 degree"
origut = "fractions"
finltr = "monthly"
finlsr = "0.5 degree"
finltr = "monthly"
finlut = "%"

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
nyears    = 20
nmonth    = 12
smonth = np.asarray(['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'])
month_bnd = np.asarray([0,31,59,90,120,151,181,212,243,273,304,334,365],dtype=float)
tbnd  = np.asarray([((np.arange(nyears)*365)[:,np.newaxis]+month_bnd[:-1]).flatten(),
                    ((np.arange(nyears)*365)[:,np.newaxis]+month_bnd[+1:]).flatten()]).T
tbnd += (1997-1850)*365
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

    year    = i + 1997

    print(year)

    # read single netCDF file
    filename = DataDir + '/GFED4S/GFED4.1s_' + str(year) + '.hdf5'
    print(filename)
    gfed4s=Dataset(filename,'r',format='NETCDF4')
    print(gfed4s) 
    print(gfed4s.variables) 
    print(gfed4s.groups) 

    for j in range(nmonth):

        # Variable from multiple files.
        print(smonth[j])
        ancill = gfed4s.groups['ancill']
        burntArea = gfed4s.groups['burned_area']
        ba = burntArea.groups[smonth[j]]

        #print(ba.variables.keys())
        area1 = ancill.variables['grid_cell_area']
        ba1 = ba.variables['burned_fraction']

        print('total area burned original')
        print(np.sum(area1[:,:]))
        print(np.sum(ba1[:,:]*area1[:,:]))

        ba2 = ba1[::-1,:]
        area2 = area1[::-1,:]

        ba2[:,:] = ba2[:,:]*100.

        print('total area burned original')
        print(np.sum(area2[:,:]))
        print(np.sum(ba2[:,:]*area2[:,:]))

        print(ba1.shape)

        print(ba.variables['burned_fraction'].long_name)
        long_name = ba.variables['burned_fraction'].long_name

        data[ij,:,:] = rebin(ba2, (nlat, nlon)) 
        area[:,:]    = rebin(area2, (nlat, nlon)) 
        area[:,:]    = area[:,:]*4

        #for j1 in range(nlat):
        #    for i1 in range(nlon):
        #        #interp
        #        i2 = i1 + 1
        #        j2 = j1 + 1
        #        data[ij,j1,i1] = np.mean(ba2[j1:j2,i1:i2])
        #        area[j1,i1]    = np.sum(area2[j1:j2,i1:i2])
        #        print(i1)
        #        print(i2)
        #        print(j1)
        #        print(j2)
        #
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
for i in range(len(t)):
    desert      = np.where(mdata>0, 0, desert)
    data[i,:,:] = np.where(desert[:,:]>=70, -999, data[i,:,:])

with Dataset(DataDir + "/burntArea.nc", mode="w") as dset:

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
    D.standard_name = "burned area fraction"
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
@ARTICLE{vanderWerf2017,
  author = {van der Werf, G.R., J.T. Randerson, L. Giglio, T.T. van Leeuwen, Y. Chen, B.M. Rogers, M. Mu, M.J.E. van Marle, D.C. Morton, G.J. Collatz, R.J. Yokelson, and P.S. Kasibhatla},
  title = {Global fire emissions estimates during 1997-2016},
  journal = {Earth Syst. Sci. Data},
  year = {2017},
  number = {9},
  page = {697-720},
  doi = {https://doi.org/10.5194/essd-9-697-2017}
}
@ARTICLE{Giglio2013,
  author = {Giglio, L., J.T. Randerson and G.R. van der Werf},
  title = {Analysis of daily, monthly, and annual burned area using the fourth-generation global fire emissions database (GFED4)},
  journal = {J. of Geophys. Res.-Biogeosciences},
  year = {2013},
  number = {118(1)},
  page = {317-328},
  doi = {https://doi.org/10.1002/jgrg.20042}
}
@ARTICLE{Randerson2012,
  author = {Randerson, J.T., Y. Chen, G.R. van der Werf, B.M. Rogers and D.C. Morton},
  title = {Global burned area and biomass burning emissions from small fires},
  journal = {J. of Geophys. Res.-Biogeosciences},
  year = {2012},
  number = {117(G4)},
  page = {G04012},
  doi = {https://doi.org/10.1029/2012JG002128}
}"""
    dset.comments = """
time_period: %s; original_temporal_resolution: %s; original_spatial_resolution: %s; original_units: %s; final_temporal_resolution: %s; final_spatial_resolution: %s; final_units: %s""" % (period, origtr, origsr, origut, finltr, finlsr, finlut)
    dset.convention = "CF-1.7"
