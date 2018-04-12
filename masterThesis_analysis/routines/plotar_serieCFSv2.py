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

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

def readFiles_CFSv2(nfiles):
    '''
        Read file by file, extracting u,v and time from netcdf files.

        args
            nfiles    (list): list with all files

        returns
            wu,wv and time (arrays)
    '''
    # define some variable
    wu   = []
    wv   = []
    time = []

    # read files
    for f in nfiles:
        ncdata = xr.open_dataset(f)
        u     = ncdata['U_GRD_L103'].values
        v     = ncdata['V_GRD_L103'].values
        t     = ncdata['time'].values

        for k in [0,1,2,3]:
            wu.append(u[k,:,:])
            wv.append(v[k,:,:])
            time.append(t[k])

    wu   = np.asarray(np.squeeze(wu))
    wv   = np.asarray(np.squeeze(wv))
    time = np.asarray(time)

    return wu,wv,time

BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/', '/ventopcse/data/serie_cfsv2/LajeDeSantos/')

# select files
nfiles = glob.glob(DATA_DIR + '*.nc')
nfiles.sort()

wu,wv,time = readFiles_CFSv2(nfiles)

# convert time to pd.DateTimeIndex
i = pd.DatetimeIndex(time)

# create dataframe
cfsv2 = pd.DataFrame({'wu':wu,'wv':wv}, index=i)


cfsv2.plot(subplots=True, title='NCEP nearest point to Laje de Santos Stations for 2014')
plt.show()
