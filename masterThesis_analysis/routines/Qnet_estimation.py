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
import glob

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def load_data(fname,var):
    ncin = xr.open_dataset(fname)
    data = ncin[var].values

    return data

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

nfiles = glob.glob('/media/danilo/Danilo/mestrado/ventopcse/Qnet/cdas1*')
nfiles.sort()

for fname in nfiles:
    ncin = xr.open_dataset(fname)

    qsw = ncin['DSWRF_L1'].values
    qlw = ncin['DLWRF_L1'].values
    lat = ncin['LHTFL_L1'].values
    sen = ncin['SHTFL_L1'].values

nfiles = glob.glob('/home/danilo/Dropbox/mestrado/data/data2model/JF2014/hflx/cdas1*')
nfiles.sort()

ncin = xr.open_dataset(nfiles[0])
thfx = ncin['THFLX_L1_Avg_1'].values
