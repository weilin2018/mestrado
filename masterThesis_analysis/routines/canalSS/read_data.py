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
import scipy.io as sio

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def filtering_data(data):
    tres = 30/60.0 # Dado a cada 30 minutos! 30 min/60 min - a cada 0.5 horas.
    cutperiod = 40	# (h)
    npoints = cutperiod/tres
    x1 = np.linspace(-1, 1, npoints)
    lancfilt1 = np.sinc(x1)
    lancfilt1 = lancfilt1/np.sum(lancfilt1)

    return np.convolve(lancfilt1,data,mode='same')

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/data/IOUSP')
SAVE_DIR = BASE_DIR.replace('github','ventopcse')

# importing file .mat and desired variables
fname = DATA_DIR + 'CSS_DJF2014.mat'
ncin  = sio.loadmat(fname,
        variable_names=['TupI','TdwI','uup','vup','udw','vdw','time'])

# converting time fromordinal to readable format
time = np.squeeze(ncin['time'])
dt = []
for t in time:
    dt.append(datetime.datetime.fromordinal(int(t)) + datetime.timedelta(days = t%1) - datetime.timedelta(days = 366))

# creating dataframe
dct = {
    'T5m'  : np.squeeze(ncin['TupI']),
    'T15m' : np.squeeze(ncin['TdwI']),
    'u5m'  : np.squeeze(ncin['uup'])/100,
    'v5m'  : np.squeeze(ncin['vup'])/100,
    'u15m' : np.squeeze(ncin['udw'])/100,
    'v15m' : np.squeeze(ncin['vdw'])/100
}

df = pd.DataFrame(dct,index=pd.DatetimeIndex(dt))

# filtering with lanczos, window=40h
# removing low frequency signal (associated to tides oscillation)
dct = {
    'u5m': filtering_data(df.u5m.values),
    'v5m': filtering_data(df.v5m.values),
    'u15m': filtering_data(df.u15m.values),
    'v15m': filtering_data(df.v15m.values),
    'T5m': df.T5m.values,
    'T15m': df.T15m.values
}

dfFilt = pd.DataFrame(dct,index=pd.DatetimeIndex(dt))

# showing data
fig,axes = plt.subplots(nrows=6,sharex=True)

df.u5m.plot(ax=axes[0])
dfFilt.u5m.plot(ax=axes[0],label='u 5m')

df.v5m.plot(ax=axes[1])
dfFilt.v5m.plot(ax=axes[1],label='v 5m')

df.u15m.plot(ax=axes[2])
dfFilt.u15m.plot(ax=axes[2],label='u 15m')

df.v15m.plot(ax=axes[3])
dfFilt.v15m.plot(ax=axes[3],label='v 15m')

df.T5m.plot(ax=axes[4],label='T 5m')
df.T15m.plot(ax=axes[5],label='T 15m')

# saving all data into a netcdf file
save = dfFilt.to_xarray()
save.to_netcdf(SAVE_DIR + 'CSS_filtered.nc')
