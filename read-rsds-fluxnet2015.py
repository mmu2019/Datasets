import numpy as np
from netCDF4 import Dataset
import pylab as plt
import os
import time
import datetime
import netCDF4 as nc
import math
from urllib.request import urlretrieve

# set up Data directory
DataDir = "/Users/mingquan/newDATA"

# set up initial and final years of data
start_yr = 1991
end_yr   = 2014

VarID       = "rsds"
RawVarID    = "SW_IN_F"
RawVarID_QC = RawVarID + "_QC"
long_name   = "surface downward shortwave radiation"

# Set general information for the data source
remote_source = "https://fluxnet.fluxdata.org/data/fluxnet2015-dataset/"
gist_source = "https://github.com/mmu2019/Datasets/blob/master/read-rsds-fluxnet2015.py"
local_source = DataDir + '/FLUXNET2015/TIER1/FULLSET/MM/FLX_STATIONNAME_FLUXNET2015_FULLSET_MM_YEAR.csv'
stamp1 = '2019-06-20'

datestr = str(datetime.datetime.now())
TmpStr  = datestr.split(' ')
stamp2  = TmpStr[0]

instit1 = "FluxNet, AmeriFlux, AfriFlux, AsiaFlux, ChinaFlux, Fluxnet-Canada, and KoFlux"
instit2 = "CarboAfrica, CarboEuropeIP, CarboItaly, CarboMont, GreenGrass, and OzFlux-TERN"
instit3 = "LBA, NECC, ICOS, TCOS-Siberia, and USCCC"

period = period = str(start_yr) + "-01 through " + str(end_yr) + "-12"
origtr = "monthly"
origsr = "site"
origut = "W/m2"
finltr = "monthly"
finlsr = "site"
finlut = "W/m2"

# Create temporal dimension
nyears = end_yr - start_yr + 1
nmonth = 12
smonth = np.asarray(['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'])
month_bnd = np.asarray([0,31,59,90,120,151,181,212,243,273,304,334,365],dtype=float)
tbnd  = np.asarray([((np.arange(nyears)*365)[:,np.newaxis]+month_bnd[:-1]).flatten(),
                    ((np.arange(nyears)*365)[:,np.newaxis]+month_bnd[+1:]).flatten()]).T
tbnd += (1991-1850)*365
tbnd.shape
t     = tbnd.mean(axis=1)
t.shape

# Input all site info
SiteInfFile = DataDir + '/FLUXNET2015/README/FluxNet2015_sites_info.txt'

#infrec = np.genfromtxt(SiteInfFile, delimiter="\t", deletechars="#", invalid_raise = False, dtype=str)
infrec = np.genfromtxt(SiteInfFile, delimiter=" ", deletechars="#", dtype=str)
nsite = len(infrec)
headerstr = infrec[0].split('\t')

indx1   = headerstr.index("SITE_ID")
indx2   = headerstr.index("SITE_NAME")
indx3   = headerstr.index("LOCATION_LAT")
indx4   = headerstr.index("LOCATION_LONG")
indx5   = headerstr.index("LOCATION_ELEV")
indx6   = headerstr.index("IGBP")

AllSiteLats  = np.ma.masked_array(np.random.rand(nsite-1))
AllSiteLons  = np.ma.masked_array(np.random.rand(nsite-1))
AllSiteElev  = np.ma.masked_array(np.random.rand(nsite-1))

AllSiteIDs   = []
AllSiteNames = []
AllSiteIGBP  = []

for i in range(nsite-1):
    datastr      = infrec[i+1].split('\t')
    AllSiteIDs.append(datastr[indx1])
    AllSiteNames.append(datastr[indx2])
    AllSiteIGBP.append(datastr[indx6])
    AllSiteLats[i]  = float(datastr[indx3])
    AllSiteLons[i]  = float(datastr[indx4])
    tempstr = datastr[indx5]
    if tempstr=='':
       tempstr = "-999"

    AllSiteElev[i] = float(tempstr)

# Set data directory
ThisDir = DataDir + '/FLUXNET2015/TIER1/FULLSET/MM/'

# list all data files in directory ThisDir
AllFileNames = os.listdir(ThisDir)

nfile = len(AllFileNames)

lon    = np.ma.masked_array(np.random.rand(nfile))
lat    = np.ma.masked_array(np.random.rand(nfile))
lev    = np.ma.masked_array(np.random.rand(nfile))
siteID = np.ma.masked_array(np.random.rand(nfile))
data   = np.ma.masked_array(np.random.rand(t.size, nfile), fill_value=-999)
data2D = np.ma.masked_array(np.random.rand(nyears,nmonth), fill_value=-999)

lat[:]      = -999
lon[:]      = -999
data[:,:]   = -999

site_id   = [] 
site_name = [] 
site_igbp = [] 

ij = 0
for FileName in AllFileNames:

    print(FileName)

    data2D[:,:] = -999

    tempstr = FileName.split("_")
    SiteID  = tempstr[1]
    site_id.append(SiteID)

    #indx0 = np.where(AllSiteIDs==SiteID)
    indx0 = AllSiteIDs.index(SiteID)

    lat[ij] = AllSiteLats[indx0]
    lon[ij] = AllSiteLons[indx0]
    lev[ij] = AllSiteElev[indx0]

    site_name.append(AllSiteNames[indx0])
    site_igbp.append(AllSiteIGBP[indx0])

    datarec = np.genfromtxt(ThisDir+FileName, delimiter="\t", deletechars="#", dtype=str)
    ndata = len(datarec)
    datastr = datarec[0].split(',')

    if RawVarID in datastr:
       indx    = datastr.index(RawVarID)
       indx_qc = datastr.index(RawVarID_QC)
    else:
       indx = -1

    if indx != -1:
       for i in range(ndata-1):
           datastr = datarec[i+1].split(',')

           YYMM = int(datastr[0])

           YY = int(YYMM/100)
           MM = int(YYMM - YY*100)

           tempdata    = float(datastr[indx])
           tempdata_qc = float(datastr[indx_qc])

           if tempdata<=-990 or tempdata_qc<0.5 and tempdata_qc>=0:
              # reset missing value
              tempdata = -999.

           # only data in the period from start_yr till end_yr are chosen.
           if YY>=start_yr and YY<=end_yr:
              iy = YY - start_yr
              im = MM - 1
              data2D[iy,im] = tempdata

    ijk = 0
    for iy in range(nyears):
        for im in range(nmonth):
            data[ijk,ij] = data2D[iy,im]
            ijk = ijk + 1

    siteID[ij] = ij + 1

    ij = ij + 1

data_min = data.min()
data_max = data.max()

# Calculate climatology of burned area
mdata     = data.mean(axis=0)

with Dataset(DataDir + "/rsds.nc", mode="w") as dset:

    # dimensions
    dset.createDimension("time", size=t.size)
    dset.createDimension("data", size=nfile)
    dset.createDimension("nb", size=2)

    # time
    T = dset.createVariable("time", t.dtype, ("time"))
    T[...]          = t
    T.units         = "days since 1850-01-01 00:00:00"
    T.calendar      = "noleap"
    T.bounds        = "time_bounds"
    T.standard_name = "time"
    T.long_name     = "time"

    # time bounds
    TB      = dset.createVariable("time_bounds", t.dtype, ("time", "nb"))
    TB[...] = tbnd

    # latitude
    X               = dset.createVariable("lat", lat.dtype, ("data"))
    X[...]          = lat
    X.standard_name = "latitude"
    X.long_name     = "site latitude"
    X.units         = "degrees_north"

    # longitude
    Y               = dset.createVariable("lon", lon.dtype, ("data"))
    Y[...]          = lon
    Y.standard_name = "longitude"
    Y.long_name     = "site longitude"
    Y.units         = "degrees_east"

    # elevation
    Z               = dset.createVariable("elevation", lev.dtype, ("data"))
    Z[...]          = lev
    Z.units         = "m"
    Z.positive      = "up"
    
    # data
    D = dset.createVariable(VarID, data.dtype, ("time", "data"), fill_value = -999)
    D[...]          = data
    D.units         = "W/m2"
    D.standard_name = long_name
    D.long_name     = long_name
    D.actual_range  = np.asarray([data_min,data_max])

    # site_info
    S = dset.createVariable("site_info", int, ("data"))
    S[...]          = siteID
    S.site_id       = site_id
    S.site_name     = site_name
    S.IGBP_class    = site_igbp
    
    dset.title = "FluxNet Tower eddy covariance measurements TIER1"
    dset.version = "2015"
    dset.institutions = "%s; %s; %s" % (instit1, instit2, instit3)
    dset.source = "Shortwave radiation, incoming consolidated from SW_IN_F_MDS and SW_IN_ERA (negative values set to zero)"
    dset.history = """
%s: downloaded source from %s;
%s: converted to netCDF with %s""" % (stamp1, remote_source, stamp2, gist_source)
    dset.references  = """
@ARTICLE{Reichstein2007,
  author = {Reichstein, M., D. Papale, R. Valentini, M. Aubinet, C. Bernhofer, A. Knohl, T. Laurila, A. Lindroth, E. Moors, K. Pilegaard, and G. Seufert},
  title = {Determinants of terrestrialecosystem carbon balance inferred from European eddy covarianceflux sites},
  journal = {Geophys. Res. Lett.},
  year = {2007},
  number = {34},
  page = {L01402},
  doi = {https://doi.org/doi:10.1029/2006GL027880}
}
@ARTICLE{Lasslop2010,
  author = {Lasslop, G., M. Reichstein, D. Papale, A.D. Richardson, A. Arneth, A. Barr, P. Stoy, and G. Wohlfahrt},
  title = {Separation of net ecosystem exchange into assimilation and respiration using a light response curve approach: critical issues and global evaluation},
  journal = {Global Change Biology},
  year = {2010},
  number = {16},
  page = {187-208},
  doi = {https://doi.org/10.1111/j.1365-2486.2009.02041.x}
}
@ARTICLE{Knauer2018,
  author = {Knauer, J., S. Zaehle, B.E. Medlyn, M. Reichstein, C.A. Williams, M. Migliavacca, M.G. De Kauwe, C. Werner, C. Keitel, P. Kolari, J.-M. Limousin, and M.-L. Linderson},
  title = {Towards physiologically meaningful water use efficiency estimates from eddy covariance data},
  journal = {Global Change Biology},
  year = {2018},
  number = {24(2)},
  page = {694-710},
  doi = {https://doi.org/10.1111/gcb.13893}
}"""
    dset.comments = """
time_period: %s; original_temporal_resolution: %s; original_spatial_resolution: %s; original_units: %s; final_temporal_resolution: %s; final_spatial_resolution: %s; final_units: %s""" % (period, origtr, origsr, origut, finltr, finlsr, finlut)
    dset.convention = "CF-1.7"
