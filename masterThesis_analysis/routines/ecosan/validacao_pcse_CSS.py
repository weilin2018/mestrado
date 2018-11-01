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

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
MODEL_DIR = '/media/danilo/Danilo/mestrado/ventopcse/output/control2006_Qnet.cdf'
DATA_DIR  = '/media/danilo/Danilo/mestrado/ventopcse/data/LajedeSantos_ECOSAN_2006.nc'

observ = xr.open_dataset(DATA_DIR)
model  = xr.open_dataset(MODEL_DIR)

# closest position between observed and grid data
ilon = 57
ilat = 10

# extract data from gcmplt.cdf
lon = model.lon.values
lat = model.lat.values
dep = model.depth[ilon,ilat].values
sig = model.sigma.values
lon[lon==0.]=np.nan
lat[lat==0.]=np.nan

temp = model.temp[:,:,ilon,ilat]

# create new depth based on sigma * local depth
depth = sig * dep
idep = 3 # closest depth to 5m


# creating dataframe to facilitate comparation
dtRange_observ = pd.DatetimeIndex(observ.index.values)
dtRange  = pd.date_range(start='2006-01-01 01:30',end='2006-02-28 22:30',freq='3H')

modeledTemp  = pd.DataFrame({'temp':temp[:,idep]},index=pd.DatetimeIndex(model.time.values))
observedTemp = pd.DataFrame({'temp':observ.temperature.values},index=dtRange)



fig,ax = plt.subplots(ncols=2)

ax[0].plot(model.time.values,temp[:,idep])
