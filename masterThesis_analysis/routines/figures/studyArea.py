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

plt.ion()
##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def import_data(fname):
    ncin = xr.open_dataset(fname)
    lon = ncin.lon.values
    lat = ncin.lat.values
    lon[lon==0] = np.nan
    lat[lat==0] = np.nan

    dep = ncin.depth.values

    ncin.close()

    return lon,lat,dep


def make_map(ax,llat=-30,ulat=-20,llon=-50,ulon=-39,resolution='l',nmeridians=3,nparallels=2,labels=[True,False,False,True]):

    if resolution == 'f':
        m = pickle.load(open('pickles/basemap.p','r'))
    else:
        m = Basemap(projection='merc', llcrnrlat=llat, urcrnrlat=ulat, llcrnrlon=llon, urcrnrlon=ulon, resolution=resolution)

    m.ax = ax

    m.drawcoastlines(linewidth=.1,color='k')
    m.drawmapboundary()
    m.fillcontinents(color='#c0c0c0')
    m.drawstates(linewidth=.01,color='gray')
	# definir meridianos e paralelos para plotar no mapa
    meridians=np.arange(llon,ulon,nmeridians)
    parallels=np.arange(llat,ulat,nparallels)
	# desenhar meridianos e paralelos conforme definido acima
    m.drawparallels(parallels,labels=labels,fontsize=8,color='gray',linewidth=.2)
    m.drawmeridians(meridians,labels=labels,fontsize=8,color='gray',linewidth=.2)

    return m
##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/output')

### FIGURE CONFIGURATION
figsize      = (8.4/2.54, 10./2.54)
cax_position = [0.205,0.12,0.76,0.02]

# import depth
lon,lat,depth   = import_data(DATA_DIR+"EC1.cdf")

# plotando mapa geral
fig,ax = plt.subplots(figsize=figsize)
mGeneral = make_map(ax,resolution='f')
mGeneral.contourf(lon,lat,depth,cmap=cmo.cm.deep,latlon=True)

# criando minimapas
