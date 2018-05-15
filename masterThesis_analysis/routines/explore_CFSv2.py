"""
    Explore the data from 2011-01 to 2017-12, looking for
    periods with a positive along shore velocity component persisting
    for more than 5 days, much more above the normal conditions of
    frontal systems (cold front).

    CLIMATE FORECAST SYSTEM VERSION 2 (CFSv2)
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

def read_cfsv2(nfiles):
    """Read and extract data from CFSv2 files.

    Parameters
    ----------
    nfiles : list
        Sorted by name list.

    Returns
    -------
    wu, wv, time: arrays


    """
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
            wu.append(np.squeeze(u[k,:,:]))
            wv.append(np.squeeze(v[k,:,:]))
            time.append(np.squeeze(t[k]))

    wu   = np.asarray(np.squeeze(wu))
    wv   = np.asarray(np.squeeze(wv))
    time = np.asarray(time)

    return wu,wv,time

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

os.system('clear')
# define some constants
BASE_DIR = oceano.make_dir()
NCEP_DIR = BASE_DIR.replace('github/', 'ventopcse/data/CFSv2/data/')

# select files, creating a list with all files to be read
nfiles = glob.glob(NCEP_DIR+'cdas1*')
nfiles.sort()

### extract timeseries
# os.system('Extracting data from a singlepoint')
wu,wv,time = read_cfsv2(nfiles)

# put data into pandas.DataFrame to better visualization
os.system('Converting data into DataFrame')
cfsv2 = pd.DataFrame({'wu':wu,'wv':wv},index=pd.DatetimeIndex(time))
