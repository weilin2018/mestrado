# add some description here

import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import pandas as pd
import os
import pickle
from scipy.interpolate import griddata
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import dates
import datetime
import cmocean as cmo

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def load_ghrsst(fname):

    # load netcdf file
    ncdata = xr.open_dataset(fname)
    # extract data, converting from K to oC
    sst = ncdata['analysed_sst'].values - 273.15
    time = ncdata['time'].values

    return sst,time

def load_erainterim(fnames):

    data = []
    time = []

    for i,fname in enumerate(fnames):
        ncdata = xr.open_dataset(fname)

        # if fname is the first file, then load lon,lat
        if i == 0:
            lon = ncdata['longitude'].values - 360
            lat = ncdata['latitude'].values

            lon,lat = np.meshgrid(lon,lat)

        for t in range(len(ncdata['time'])):
            time.append(ncdata['time'].values[t])
            data.append(ncdata['ishf'].values[t,:,:])

    # convert lists into np.ndarray
    data = np.asarray(data)
    time = np.asarray(time)

    return lon,lat,data,time



def load_csfv2(fnames,datakey):

    # cada arquivo com 4 passos de tempo [00+1, 06+1, 12+1, 18+1]
    # SHTFL_L1 is the variable for heat flux

    data = []
    time = []


    for i,fname in enumerate(fnames):
        ncdata = xr.open_dataset(fname)

        # if fname is the first file, then load lon,lat
        if i == 0:
            lon = ncdata['lon'].values - 360
            lat = ncdata['lat'].values

            lon,lat = np.meshgrid(lon,lat)

        for t in range(len(ncdata['time'])):
            time.append(ncdata['time'].values[t])
            data.append(ncdata[datakey].values[t,:,:])

    # convert lists into np.ndarray
    data = np.asarray(data)
    time = np.asarray(time)

    return lon,lat,data,time


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

BASE_DIR = oceano.make_dir()
if BASE_DIR.split("/")[2] == 'tparente':
    DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output_modelo/exp03_variables/')

else:
    DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
    DATA_DIR = '/home/danilo/Dropbox/mestrado/data/data2model/JF2014/'

#####-------------------------- loading

# defining directories with every dataset files.
THFLX_DIR = DATA_DIR + 'hflx/'                                          # Total Downward Heat Flux (CFSv2)
SHFLX_DIR = DATA_DIR + 'sensible_hflx_cfsv2/'                           # Sensible Heat Flux (CFSv2)
IHFLX_DIR = DATA_DIR + 'instantaneous_sensible_heatflux_erainterim/'    # Instantaneous Sensible Heat Flux (ERA INTERIM)

# loading files from each dataset
thflx_fnames = glob.glob(THFLX_DIR+'*.nc')
thflx_fnames.sort()
shflx_fnames = glob.glob(SHFLX_DIR+'*.nc')
shflx_fnames.sort()
ihflx_fnames = glob.glob(IHFLX_DIR+'*.nc')
ihflx_fnames.sort()

# loading data from each dataset
lon_thflx,lat_thflx,thflx,time_thflx = load_csfv2(thflx_fnames,datakey='THFLX_L1_Avg_1')
lon_shflx,lat_shflx,shflx,time_shflx = load_csfv2(shflx_fnames,datakey='SHTFL_L1')

lon_ihflx,lat_ihflx,ihflx,time_ihflx = load_erainterim(ihflx_fnames)
# cutting data from ERA INTERIM
ihflx = ihflx[:127,:,:]
time_ihflx = time_ihflx[:127]

#####-------------------------- plotting


# create structure
fig,ax = plt.subplots(ncols=3)

contour_levels = np.arange(-800,800,1600./1000)

m1 = oceano.make_map(ax[0],resolution='i')
m2 = oceano.make_map(ax[1],resolution='i')
m3 = oceano.make_map(ax[2],resolution='i')

m1.contourf(lon_thflx,lat_thflx,thflx[-1,:,:],contour_levels,latlon=True,cmap=cmo.cm.thermal)
m2.contourf(lon_shflx,lat_shflx,shflx[-1,:,:],contour_levels,latlon=True,cmap=cmo.cm.thermal)
m3.contourf(lon_ihflx,lat_ihflx,ihflx[-1,:,:],contour_levels,latlon=True,cmap=cmo.cm.thermal)
