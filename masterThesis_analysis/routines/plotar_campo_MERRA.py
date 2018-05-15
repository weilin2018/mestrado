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

BASE_DIR = '/media/danilo/Danilo/mestrado/ventopcse/data/MERRA/'
SAVE_DIR = BASE_DIR+'output/'

##### SUMMER 2014
nfiles = glob.glob(BASE_DIR+'data/cut/*.nc')
nfiles.sort()
SAVE_DIR = BASE_DIR+'data/cut//outputs/'

for f in nfiles:
    ncdata = xr.open_dataset(f)

    x = ncdata['lon'].values
    y = ncdata['lat'].values
    lon,lat = np.meshgrid(x,y)

    wu = ncdata['U10M'].values
    wu = wu.mean(axis=0) #media diaria
    wv = ncdata['V10M'].values
    wv = wv.mean(axis=0) #media diaria

    spd = np.sqrt(wu**2 + wv**2) # compute intensity

    # normalizar vetores
    wun = wu/spd
    wvn = wv/spd

    fig, ax = plt.subplots(figsize=(15,13))
    m = oceano.make_map(ax)

    cb = m.contourf(lon,lat,spd,latlon=True)
    m.quiver(lon,lat,wun,wvn,latlon=True)

    plt.title('Daily mean '+f[-15:-7])
    cbar = plt.colorbar(cb)
    cbar.set_clim(0,20.)

    fout = f[-15:-7]+'.png'
    plt.savefig(SAVE_DIR+fout)

    plt.close('all')
