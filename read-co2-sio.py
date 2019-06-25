import numpy as np
from netCDF4 import Dataset
import pylab as plt
import os
import time
from urllib.request import urlretrieve

remote_source = "http://scrippsco2.ucsd.edu/assets/data/atmospheric/stations/in_situ_co2/monthly/monthly_in_situ_co2_mlo.csv"
gist_source = "https://github.com/mmu2019/Datasets/blob/master/read-co2-sio.py"
local_source = os.path.basename(remote_source)
if not os.path.isfile(local_source):
    urlretrieve(remote_source, local_source)
stamp = time.strftime('%Y-%m-%d', time.localtime(os.path.getmtime(local_source)))

print(remote_source)
print(gist_source)
print(local_source)
print(stamp)

# parse CSV file, we use the 9th column for the CO2 data
rec = np.genfromtxt("monthly_in_situ_co2_mlo.csv", delimiter=",", skip_header=57)
year = rec[:, 0]
month = rec[:, 1].astype(int)
co2 = np.ma.masked_values(rec[:, 8], -99.99).reshape((-1, 1)).astype(np.float32)

print(year)
print(month)

# create time and its bounds in days since 1850-1-1
bnd_months = np.asarray([0., 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365])
tb = ((year-1850)*365)[:, np.newaxis] + np.asarray([bnd_months[month-1], bnd_months[month]]).T
t = tb.mean(axis=1)

# location information pulled from webpage
lat = np.asarray([19.5]).astype(np.float32)
lon = np.asarray([-155.6]).astype(np.float32)
elevation = np.asarray([3397.]).astype(np.float32)

start_yr  = year[0]
end_yr    = year[t.size-1]

start_mon = month[0]
end_mon   = month[t.size-1]

print(start_yr)
print(end_yr)

print(start_mon)
print(end_mon)

period = str(start_yr) + "-" + str(start_mon) + " through " + str(end_yr) + "-" + str(end_mon)
origtr = "monthly"
origsr = "site"
origut = "ppm"
finltr = "monthly"
finlsr = "site"
finlut = "ppm"

with Dataset("co2.nc", mode="w") as dset:

    # dimensions
    dset.createDimension("time", size=t.size)
    dset.createDimension("data", size=1)
    dset.createDimension("nb", size=2)

    # time
    T = dset.createVariable("time", t.dtype, ("time"))
    T[...] = t
    T.units = "days since 1850-01-01 00:00:00"
    T.calendar = "noleap"
    T.bounds = "time_bounds"

    # time bounds
    TB = dset.createVariable("time_bounds", t.dtype, ("time", "nb"))
    TB[...] = tb

    # latitude
    X = dset.createVariable("lat", lat.dtype, ("data"))
    X[...] = lat
    X.standard_name = "latitude"
    X.long_name = "site latitude"
    X.units = "degrees_north"

    # longitude
    Y = dset.createVariable("lon", lon.dtype, ("data"))
    Y[...] = lon
    Y.standard_name = "longitude"
    Y.long_name = "site longitude"
    Y.units = "degrees_east"

    # elevation
    Z = dset.createVariable("elevation", elevation.dtype, ("data"))
    Z[...] = elevation
    Z.units = "m"
    Z.positive = "up"
    
    # data
    D = dset.createVariable("co2", co2.dtype, ("time", "data"), fill_value = -99.99)
    D[...] = co2
    D.units = "ppm"
    D.standard_name = "atmosphere_moles_of_carbon_dioxide"
    D.long_name = "CO2 concentration"
    D.actual_range = np.asarray([co2.min(),co2.max()])
    
    dset.title = "Primary Mauna Loa CO2 Record (1958 - present)"
    dset.version = "2017"
    dset.institution = "Scripps Institution of Oceanography"
    dset.source = "in situ air measurements"
    dset.history = """
%s: downloaded source from %s
%s: converted to netCDF with %s""" % (stamp, remote_source, stamp, gist_source)
    dset.references = """
@ARTICLE{Keeling2001,
  author = {Keeling, C. D., S. C. Piper, R. B. Bacastow, M. Wahlen, T. P. Whorf, M. Heimann, and H. A. Meijer},
  title = {Exchanges of atmospheric CO2 and 13CO2 with the terrestrial biosphere and oceans from 1978 to 2000},
  journal = {I. Global aspects, SIO Reference Series},
  year = {2001},
  number = {01-06},
  pages = {88 pages},
  doi = {Scripps Institution of Oceanography, San Diego}
}
@ARTICLE{Keeling2005,
  author = {Keeling, C. D., S. C. Piper, R. B. Bacastow, M. Wahlen, T. P. Whorf, M. Heimann, and H. A. Meijer},
  title = {A History of Atmospheric {CO2} and its effects on Plants, Animals, and Ecosystems},
  Journal = {Atmospheric CO2 and 13CO2 exchange with the terrestrial biosphere and oceans from 1978 to 2000: observations and carbon cycle implications, edited by J. R. Ehleringer and  T. E. Cerling and M. D. Dearing},
  year = {2005},
  pages = {83-113},
  doi = {Springer Verlag, New York}
}"""
    dset.comments = """
time_period: %s; original_temporal_resolution: %s; original_spatial_resolution: %s; original_units: %s; final_temporal_resolution: %s; final_spatial_resolution: %s; final_units: %s""" % (period, origtr, origsr, origut, finltr, finlsr, finlut)
    dset.convention = "CF-1.7"
