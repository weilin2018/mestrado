"""
Script to plot SBB's model_grid and CFSv2 grid and detect some point between
those two grids in the same coordinate or closest.

we need to do this, so we can compare wind data in the output file
with the wind data from cfsv2, to be sure we're forcing with the same
wind.
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
BASE_DIR = oceano.make_dir()
SBB_GRID = BASE_DIR.replace('github/', 'ventopcse/output/')
CFS_GRID = "/home/danilo/Dropbox/mestrado/data/data2model/JF2014/tuv/"

# load sbb grid
sbb = xr.open_dataset(glob.glob(SBB_GRID+"*.cdf")[0])
# load cfsv2 grid
cfs = xr.open_dataset(glob.glob(CFS_GRID+"*.nc")[0])

# extracting grids from ncfiles: lon_m,lat_m (sbb model_grid) and
# lon_r,lat_r (cfsv2 grid)
lon_m = sbb['lon'].values
lat_m = sbb['lat'].values

lon_r = cfs['lon'].values - 360
lat_r = cfs['lat'].values

# clearing nan values
lon_m[lon_m == 0.] = np.nan
lat_m[lat_m == 0.] = np.nan
lon_r[lon_r == 0.] = np.nan
lat_r[lat_r == 0.] = np.nan

# regridding reanalysis grid
lon_r,lat_r = np.meshgrid(lon_r,lat_r)

# plot grid

plt.ion()

fig,ax = plt.subplots(figsize=(15,15))
m = oceano.make_map(ax)

# m.plot(lon_m,lat_m,'k',alpha=.5,latlon=True)
# m.plot(lon_m.T,lat_m.T,'k',alpha=.5,latlon=True)

# m.plot(lon_r,lat_r,'r',alpha=.5,latlon=True)
# m.plot(lon_r.T,lat_r.T,'r',alpha=.5,latlon=True)

"""
Pensando ...

preciso fazer dois loops para poder percorrer cada par i,j na matriz CFSv2 e
buscar pelos pontos mais próximos dessa coordenada na matriz SBB.

armazenar os resultados, tanto os indices da matriz SBB, mas também a
distância euclidiana de cada ponto.

usar a menor distância euclidiana para definir o ponto mais próximo entre as
grades. A ideia é pegar um ponto que bata a intersecção das grades
(quase) perfeitamente.


Ponto de partida: analisar a rotina find_nearest e observar como ela é
executada.

"""

for i in np.arange(0,10):
    for j in np.arange(0,5):
        m.scatter(x[i,j]. y[i,j], marker='o', c='k', latlon=True)
        plt.pause(1)

for i in np.arange(0,lon_r.shape[0]):
    for j in np.arange(0,lon_r.shape[1]):
        m.scatter(lon_r[j,i],lat_r[j,i],c='k',marker='o',latlon=True)

        plt.pause(1)
