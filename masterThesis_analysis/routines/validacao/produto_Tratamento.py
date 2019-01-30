# Le o arquivo de saida do modelo, faz tratamento e deixa pronto para validacao

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

import decomp as dp

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
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/output')
SAVE_DIR = BASE_DIR.replace('github','ventopcse')
experiment = 'EA2.cdf'
ncin = xr.open_dataset(DATA_DIR + experiment)

# extracting data
k5 = 6
k15= 18

u5 = ncin.u[:,k5,56,16].values
v5 = ncin.v[:,k5,56,16].values
u15= ncin.u[:,k15,56,16].values
v15= ncin.u[:,k15,56,16].values

dct = {
    'u5' : u5,
    'v5' : v5,
    'u15' : u15,
    'v15' : v15
}

dfModel = pd.DataFrame(dct,index=pd.DatetimeIndex(ncin.time.values))

# rotation of 51 degrees
alpha = np.deg2rad(51)
# for 5m depth first
INT,DIR = dp.uv2intdir(u5,v5,0,0)
ur,vr = dp.intdir2uv(INT,DIR,0,alpha)

df5m = pd.DataFrame({'along 5m':vr,'cross 5m':ur},index=dfModel.index)

# now for 15m depth
INT,DIR = dp.uv2intdir(u15,v15,0,0)
ur,vr = dp.intdir2uv(INT,DIR,0,alpha)

df15m = pd.DataFrame({'along 15m':vr,'cross 15m':ur},index=dfModel.index)


# saving dataframes
df5m = df5m.to_xarray()
df5m.to_netcdf(SAVE_DIR + 'df5m.nc')
df15m = df15m.to_xarray()
df15m.to_netcdf(SAVE_DIR + 'df15m.nc')

# for temperature
T5 = ncin.temp[:,k5,56,16].values
T15= ncin.temp[:,k15,56,16].values


dfTemp = pd.DataFrame({'T5':T5,'T15':T15},index=pd.DatetimeIndex(ncin.time.values))
dfTemp = dfTemp.to_xarray()
dfTemp.to_netcdf(SAVE_DIR + 'dfTemp.nc')
