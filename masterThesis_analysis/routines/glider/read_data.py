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
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/data/USP03_Missao01')

# listing all netcdf files
nfile = DATA_DIR + 'merge_all_ascii_lvl1.nc'
# nfiles = glob.glob(DATA_DIR + "*.nc")
# nfiles.sort() # and sorting them by name

# load merged file
ncin = xr.open_dataset(nfiles[0])
# extract time axis, converting into string readable
time = np.asarray([datetime.datetime.fromtimestamp(t) for t in ncin.m_present_time.values])
depth= gsw.z_from_p(ncin.sci_water_pressure.values,lat=-24.2)
# creating dataframe, just because
df = pd.DataFrame({'depth': depth,'chlor-a':ncin.sci_flbbcd_chlor_units.values,'cdom':ncin.sci_flbbcd_cdom_units.values},index=pd.DatetimeIndex(time))

dfCut = df['2014-11-16':'2014-12-25'].copy()

dfCut.plot(subplots=True)
