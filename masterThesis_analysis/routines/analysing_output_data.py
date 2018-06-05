# add some description here
"""
19,16
"""

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
def extract_data(ncdata,list_vars):
    """Short summary.

    Parameters
    ----------
    ncdata : xarray.core.dataset.Dataset
        NetCDF output file from ECOM.
    list_vars : list
        List of variables we want to extract from ncdata

    Returns
    -------
    output : dictionary
        Dictionary with all variables and data.
    """
    output = {}

    return output

def plot_timeserie(index,u,v,ilon=32,ilat=2):

    # create dataframe
    data = pd.DataFrame({'u':u[:,2,ilon,ilat], 'v': v[:,2,ilon,ilat]},index=index)

    data.plot(subplots=True)
    plt.show()


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
FIGU_DIR = BASE_DIR + 'masterThesis_analysis/figures/experiments_outputs/'
OUT_FILE = 'exp01.cdf'

# reading netcdf output file from ECOM
ncdata = xr.open_dataset(DATA_DIR+OUT_FILE)

# extracting data from ncdata
# list of variables we want to extract from output file
LIST_VAR = ['time', 'lon', 'lat', 'wu', 'wv', 'elev', 'u', 'v', 'w']

for var in LIST_VAR:
    exec("%s_data = ncdata['%s'].values" % (var,var))

# organizing date - only in case run_data was set with a wrong date
dt_range = pd.date_range(start='2014-01-01 01:30', end='2014-02-25 01:30', freq='3H')

# lon and lat null data
lon_data[lon_data == 0.] = np.nan
lat_data[lat_data == 0.] = np.nan

# plotting data as animation
fig, ax = plt.subplots(figsize=(16,8))

for i in np.arange(0,len(time)-430):
    plt.gca().clear()

    m = oceano.make_map(ax, ulon=np.nanmax(lon_data))

    x,y = m(lon_data,lat_data)
    u,v = wu_data[i,:,:], wv_data[i,:,:]
    eta = elev_data[i,:,:]
    spd = np.sqrt(u**2 + v**2)
    u = u/spd
    v = v/spd

    c = m.contourf(x,y,spd,extend='max')
    q = m.quiver(x[::3,::3],y[::3,::3],u[::3,::3]*2,v[::3,::3]*2,alpha=.7,scale=150,width=0.005,minshaft=2)

    m.ax.set_title(str(dt_range[i]))

    # cb = plt.colorbar(c,orientation='horizontal',fraction=0.057,pad=.06)
    # cb.set_label('Surface Current')

    plt.pause(1)
