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
def load_multivariate_eof(DATA_DIR):
    return False

def reconstructFieldWind(eofs):

    mode1st = np.sqrt(eofs[0][0,:,:]**2 + eofs[1][0,:,:]**2)


def getLatLon(fname='/media/danilo/Danilo/mestrado/ventopcse/data/CFSR/1992_2011/wnd10m.gdas.201012.grb2.nc'):

    ncin = xr.open_dataset(fname)

    lon = ncin['lon'].values - 360
    lat = ncin['lat'].values

    lon[lon == 0.] = np.nan
    lat[lat == 0.] = np.nan

    return lat,lon

def normalizeVectors(spd,wu,wv):

    wun = wu/spd
    wvn = wv/spd

    return wun*100,wvn*100

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()

DATA_DIR = '/home/danilo/Dropbox/mestrado/eof_analysis/'

# import latitude and longitude
lat,lon = getLatLon()
lon,lat = np.meshgrid(lon,lat)

# load multivariate eofs
nfiles = glob.glob(DATA_DIR+'combinada/jan/*.pickle')
nfiles.sort()

fname = nfiles[0]

plt.ion()
fig,ax = plt.subplots(ncols=2)

for i in range(ax.shape[0]):

    if i == 0:
        mode = 'First'
    elif i == 1:
        mode = 'Second'

    ax[i].set_title(mode+' mode')

m1 = oceano.make_map(ax[0],resolution='i')
m2 = oceano.make_map(ax[1],resolution='i')
# m3 = oceano.make_map(ax[2],resolution='i')

for fname in nfiles[-1:]:
    ax[0].clear()
    ax[1].clear()
    # ax[2].clear()

    m1 = oceano.make_map(ax[0],resolution='i')
    m2 = oceano.make_map(ax[1],resolution='i')
    # m3 = oceano.make_map(ax[2],resolution='i')

    data = pickle.load(open(fname,'r'))
    eofs = data['eofs']
    pcs  = data['pcs']

    firstMode  = np.sqrt(eofs[0][0,:,:]**2 + eofs[1][0,:,:]**2)
    wu_first,wv_first = eofs[0][0,:,:],eofs[1][0,:,:]
    # normalizing
    wu_first,wv_first = normalizeVectors(firstMode,wu_first,wv_first)

    secondMode = np.sqrt(eofs[0][1,:,:]**2 + eofs[1][1,:,:]**2)
    wu_second,wv_second = eofs[0][1,:,:],eofs[1][1,:,:]
    # normalizing
    wu_second,wv_second = normalizeVectors(secondMode,wu_second,wv_second)
    #
    # thirdMode  = np.sqrt(eofs[0][2,:,:]**2 + eofs[1][2,:,:]**2)
    # wu_third,wv_third = eofs[0][2,:,:],eofs[1][2,:,:]
    # # normalizing
    # wu_third,wv_third = normalizeVectors(thirdMode,wu_third,wv_third)


    m1.contourf(lon,lat,firstMode,latlon=True)
    m1.quiver(lon,lat,wu_first,wv_first,latlon=True)

    m2.contourf(lon,lat,secondMode,latlon=True)
    m2.quiver(lon,lat,wu_second,wv_second,latlon=True)

    # m3.contourf(lon,lat,thirdMode,latlon=True)
    # m3.quiver(lon,lat,wu_third,wv_third,latlon=True)

    plt.pause(0.5)
