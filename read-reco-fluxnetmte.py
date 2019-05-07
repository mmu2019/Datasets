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
remote_source = "https://www.bgc-jena.mpg.de/bgi/index.php/Services/Overview"
gist_source = "https://gist.github.com/nocollier/d73585731756fa472731065389af45dc"
local_source = 'FLUXNET.MTE/reco/YYYY/EnsembleTER_MR_May09_YYYY.nc'
stamp1 = '2013-12-02'

datestr = str(datetime.datetime.now())
TmpStr  = datestr.split(' ')
stamp2  = TmpStr[0]

print(datestr)
print(stamp2)

instit = "Department Biogeochemical Integration, Max Planck Institute for Biogeochemistry, Germany"

period = "1982-01 through 2008-12"
origtr = "monthly"
origsr = "0.5 degree"
origut = "g/m2/day"
finltr = "monthly"
finlsr = "0.5 degree"
finlut = "kg/m2/s"

# Create temporal dimension
nyears    = 27
nmonth    = 12
smonth = np.asarray(['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'])
month_bnd = np.asarray([0,31,59,90,120,151,181,212,243,273,304,334,365],dtype=float)
tbnd  = np.asarray([((np.arange(nyears)*365)[:,np.newaxis]+month_bnd[:-1]).flatten(),
                    ((np.arange(nyears)*365)[:,np.newaxis]+month_bnd[+1:]).flatten()]).T
tbnd += (1982-1850)*365
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

    year    = i + 1982

    print(year)

    # read single netCDF file
    filename = 'FLUXNET.MTE/reco/' + str(year) + '/EnsembleTER_MR_May09_' + str(year) + '.nc'
    print(filename)
    mte=Dataset(filename,'r',format='NETCDF3')
    data0 = mte.variables['EnsembleTER_MR_May09']
    lats  = mte.variables['latitude']

    long_name = data0.NameDescription

    data1 = np.where(data0[:,:,:]<=-9999, 0, data0[:,:,:])

    latrange = data0.latitude_range

    data2 = np.float_(data1[:,:,:])*data0.scale_factor + data0.add_offset

    ilat1 = np.where(lat==latrange[0])
    ilat2 = np.where(lat==latrange[1])

    j1 = int(ilat1[0])
    j2 = int(ilat2[0]) + 1

    for j in range(nmonth):

        data[ij,j1:j2,:] = data2[j,::-1,:]

        ij = ij + 1

# convert unit from g/m2/day to kg/m2/s
data[:,:,:] = data[:,:,:]/(24*3600*1000)

print(t.shape)
print(data0.shape)
print(data.shape)

data_min = data.min()
data_max = data.max()

with Dataset("reco.nc", mode="w") as dset:

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
    
    dset.title = "Fluxnet multi-tree ensemble (MTE)"
    dset.version = "May 2009"
    dset.institutions = instit
    dset.source = "Global, spatially and temporally explicit estimates of carbon and water fluxes derived from empirical up-scaling eddy covariance measurements "
    dset.history = """
%s: downloaded source from %s;
%s: converted to netCDF with %s""" % (stamp1, remote_source, stamp2, gist_source)
    dset.references  = """
@ARTICLE{Jung2011,
  author = {Jung,M., M. Reichstein, H.A. Margolis, A. Cescatti, A.D. Richardson, M.A. Arain, A. Arneth, C. Bernhofer, D. Bonal, J. Chen, D. Gianelle, N. Gobron, G. Kiely, W. Kutsch, G. Lasslop, B.E. Law, A. Lindroth, L. Merbold, L. Montagnani, E.J. Moors, D. Papale, M. Sottocornola, F. Vaccari, C. Williams},
  title = {Global patterns of land-atmosphere fluxes of carbon dioxide, latent heat, and sensible heatderived from eddy covariance, satellite, and meteorological observations},
  journal = {J. Geophys. Res.},
  year = {2011},
  number = {116},
  page = {G00J07},
  doi = {https://doi.org/10.1029/2010JG001566}
}
@ARTICLE{Beer2010,
  author = {Beer, C., M. Reichstein, E. Tomelleri, P. Ciais,M. Jung, N. Carvalhais, C. Rodenbeck, M.A. Arain, D. Baldocchi, G.B. Bonan, A. Bondeau, A. Cescatti, G. Lasslop, A. Lindroth, M. Lomas, S. Luyssaert, H. Margolis, K.W. Oleson, O. Roupsard, E. Veenendaal, N. Viovy, C. Williams, F.I. Woodward, D. Papale},
  title = {Terrestrial gross carbon dioxide uptake: Global distribution and covariation with climate},
  journal = {Science},
  year = {2010},
  number = {329},
  page = {834-838},
  doi = {https://doi.org/10.1126/science.1184984}
}
@ARTICLE{Jung2009,
  author = {Jung, M., M. Reichstein, and A. Bondeau},
  title = {Towards global empirical upscaling of FLUXNET eddy covariance observations: validation of a model tree ensemble approach using a biosphere model},
  journal = {Biogeosciences},
  year = {2009},
  number = {6},
  page = {2001-2013},
  doi = {https://doi.org/10.5194/bg-6-2001-2009}
}"""
    dset.comments = """
time_period: %s; original_temporal_resolution: %s; original_spatial_resolution: %s; original_units: %s; final_temporal_resolution: %s; final_spatial_resolution: %s; final_units: %s""" % (period, origtr, origsr, origut, finltr, finlsr, finlut)
    dset.convention = "CF-1.7"
