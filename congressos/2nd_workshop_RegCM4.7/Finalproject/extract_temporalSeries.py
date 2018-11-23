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
# beginnig of the main code
FILE_DIR = '/media/danilo/Danilo/RegCM/2nd_Workshop/Project_output/output_precipitation/output/'
SAVE_DIR = '/media/danilo/Danilo/RegCM/2nd_Workshop/Project_output/output_precipitation/figures/'
nfiles = glob.glob(FILE_DIR+"summer80.nc")

# Quilombo River location
ilon = -46.32
ilat = -23.83

# extracting coordinates
ncin = xr.open_dataset(nfiles[0])
lon = ncin.xlon
lat = ncin.xlat

i,j = 19,22

# extracting data
pr = ncin.pr[:,i,j] * 86400
pr /= 4
ti = ncin.time

df = pd.DataFrame({'Daily Precipitation [mm day-1]':pr.values},index=pd.DatetimeIndex(ti.values))

####################### SANTANDA!!!!!

### observed data
ilon = -46.61
ilat = -23.5
i,j = 20,21


# extracting data
pr = ncin.pr[:,i,j] * 86400
ti = ncin.time

df = pd.DataFrame({'Modeled Precipitation':pr.values},index=pd.DatetimeIndex(ti.values))
#resampling for daily frequency
df = df.resample('1D')

dt = pd.date_range(start='1979-11-02',end='1980-04-01',freq='1D')
data = pd.read_table('/media/danilo/Danilo/mestrado/github/congressos/2nd_workshop_RegCM4.7/Finalproject/santana_data.csv',delimiter=';',usecols=[3])
data.index = dt
data.columns = ['Observed Precipitation']

df_all = pd.DataFrame({'modeled':np.squeeze(df.values[1:]),'observed':np.squeeze(data.values)},index=df.index[1:])


fig,ax = plt.subplots(figsize=(15,5))
df_all.plot(ax=ax)
plt.suptitle(u'Daily Precipitation Observed (blue) and RegCM4.7 Output (red) for \nSão Paulo - Mirante de Santana [NDJFM 1979/80]',fontsize=18)

plt.savefig(SAVE_DIR + 'santana_plot.png',dpi=200)
