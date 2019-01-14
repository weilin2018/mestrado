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


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/data/GHRSST/14Fev_2013-2016')
plt.ion()

years = np.arange(2013,2018)

fig,ax = plt.subplots()


for year in years:
    ax.clear()
    m = oceano.make_map(ax)

    # extrair ano
    fname = DATA_DIR + str(year) + '0214.nc'
    ncin = xr.open_dataset(fname)
    lon  = ncin.lon.values
    lat  = ncin.lat.values
    sst  = np.squeeze(ncin.analysed_sst)

    # griddando coordenadas
    lon,lat = np.meshgrid(lon,lat)
    # plotando
    m.contourf(lon,lat,sst,latlon=True)
    plt.title(str(year),fontsize=20)

    plt.pause(2)
