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

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

import scipy.ndimage

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
# beginnig of the main code
BASE_DIR = oceano.make_dir()
SAVE_DIR = BASE_DIR + 'masterThesis_analysis/figures/experiments_outputs/temperature/'
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
plt.ion()

# select which experiment you want to plot:
fname = DATA_DIR + 'EA1.cdf'
ncin  = xr.open_dataset(fname)

# extract grid and other general variables
# important: lat and lon already gridded
lon   = ncin.lon.values
lat   = ncin.lat.values
lon[lon == 0.] = np.nan
lat[lat == 0.] = np.nan
depth = ncin.depth.values
sigma = ncin.sigma.values

# import 3D variable, timestep 0 just for example
temp  = ncin.temp[0,:,:,:]

# Need to create an 1D array with all coordinates pairs for the cross section
ind = 99 # latitude reference for ubatuba
xtransect = lon[ind,:]
ytransect = lat[ind,:]
ttransect = temp[:,ind,:]

# create array with coordinates pairs
line = []
for x,y in zip(xtransect,ytransect):
    line.append((x,y))
# converting geographical coordinate in pixel coordinate
row,col = np.array(zip(*line))
row = row[16:-2]
col = col[16:-2]

# creating transect with 1000 points between each data
num = 1000
row, col = [np.linspace(item[0], item[1], num) for item in [row, col]]

zi = scipy.ndimage.map_coordinates(ttransect.values[:,16:-2],np.vstack((row,col)))


#################
nan_map = np.zeros_like(ttransect)
nan_map[np.isnan(ttransect.values)] = 1

filled_z = ttransect.values.copy()
filled_z[np.isnan(ttransect.values)] = 0

f = interp2d(x, y, filled_z, kind='linear')
f_nan = interp2d(x, y, nan_map, kind='linear')
