"""
Program to read and plot wind field extracted from CFSv2 for
southwestern atlantic ocean.

Objective: visualize filtered periods of atmospherics system forcing this region,
to compare with synoptic charts and define which system was dominating
south brazil bight.

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



##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

os.system('clear')
BASE_DIR = oceano.make_dir()
SAVE_DIR = BASE_DIR + 'masterThesis_analysis/routines/pickles/'
FIGU_DIR = BASE_DIR + 'masterThesis_analysis/figures/plotting_windField/'
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/data/CFSv2/selected_periods/')
DROP_DIR = '/home/'+BASE_DIR.split("/")[2]+'/Dropbox/mestrado/selected_periods/'

# defining all selected periods in a list
periods = '20120327_20120404 20121001_20121017 20121117_20121120 20121127_20121205 20130117_20130126 20130304_20130430 20130821_20130829 20130916_20131008 20131224_20140101 20140101_20140406 20140805_20140821 20141013_20141027 20141105_20141115 20141127_20141202 20141206_20150107 20150123_20150129 20150204_20150221 20150307_20150312 20140307_20150312 20150316_20150321 20160101_20160110 20160114_20160125 20160208_20160217 20160417_20160105'.split(" ")

filt_periods = '20120327_20120404 20130117_20130126 20130304_20130430 20160114_20160125 20140101_20140406 20131224_20140101'.split(" ")


filt_periods = '20130109_20130122 20130413_20130422 20131221_20131226 20140101_20140228 20160112_20160125 20170326_20170402'.split(" ")

# for each period we have to read all files and plot. So multiples functions is needed:
def readFiles(period):

    nfiles = glob.glob(period+"/*.nc")

    return nfiles

for p in filt_periods:
    # read files from DATA_DIR+p folder
    nfiles = readFiles(DROP_DIR+p)
    nfiles.sort()

    # create a new directory to store all figures generated
    NEW_DIR = FIGU_DIR + p + '/'

    # try:
    #     os.system('mkdir %s' % NEW_DIR)
    # except:
    #     print('Not possible to create the %s directory' % NEW_DIR)

    for fname in nfiles:

        ncdata = xr.open_dataset(fname)

        lon = ncdata['lon'].values - 360
        lat = ncdata['lat'].values

        lon,lat = np.meshgrid(lon,lat)

        wu = ncdata['U_GRD_L103']
        wu = wu.mean(axis=0)
        wv = ncdata['V_GRD_L103']
        wv = wv.mean(axis=0)

        spd = np.sqrt(wu**2 + wv**2)

        # normalizar vetores
        wun = (wu/spd)
        wvn = (wv/spd)

        # now I have the u and v components from wind. What we must do next
        # is plot this data in a Basemap instance and save the final figure
        # in FIGU_DIR to visualize after
        fig, ax = plt.subplots(figsize=(15,13))

        m = oceano.make_map(ax)

        cb = m.contourf(lon,lat,spd,latlon=True)
        m.quiver(lon[::3,::3],lat[::3,::3],wun[::3,::3],wvn[::3,::3],latlon=True,scale=50)

        plt.title('Daily mean '+fname[-26:-18], fontsize=25)

        cbar = plt.colorbar(cb)
        cbar.set_clim(0,20.)

        fout = fname[-26:-18]+'.png'
        plt.savefig(NEW_DIR+fout)

        plt.close('all')
