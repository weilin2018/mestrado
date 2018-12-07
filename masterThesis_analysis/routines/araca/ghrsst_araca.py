# extracao de dados do ghrsst proximo ao ponto da baia do araca

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
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/data/GHRSST')

fname = DATA_DIR + 'ghrsst_summer2014.nc'
ncin  = xr.open_dataset(fname)

# localizando ponto de grade mais proximo do Araca
lon = ncin.lon.values
lat = ncin.lat.values

ilon = -45.403225
ilat = -23.816398

i = np.where(lon == oceano.find_nearest_1D(lon,ilon))[0][0]
j = np.where(lat == oceano.find_nearest_1D(lat,ilat))[0][0]

sst = pd.DataFrame({'sst':ncin.analysed_sst[:,j,i].values - 273.15},index=pd.DatetimeIndex(ncin.time.values))
sst.plot()
plt.show()

# ideia: completar a lacuna de dados de série temporal com os dados do ghrsst ???
