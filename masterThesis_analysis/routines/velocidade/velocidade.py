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

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def formatGrid_plot(grid,fname):
    import numpy as np
    ij=np.load(fname)
    # for a 2D array (lon,lat)
    if len(grid.shape)==2:
        grid=grid[ij[1], :]
        grid=grid[:, ij[0]]
    # if grid is a 3D array (temp,salt,speed)
    if len(grid.shape)==3:
        grid=grid[:,ij[1], ij[0]]
    return grid



##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/output')

experimento = 'EA1.cdf'
fname = DATA_DIR + experimento

ncin = xr.open_dataset(fname)
# extraindo grid
lon = ncin.lon.values
lat = ncin.lat.values
lon[lon == 0.] = np.nan
lat[lat == 0.] = np.nan

# extraindo velocidades na superficie
u = ncin.u[:,0,:,:]
v = ncin.v[:,0,:,:]
s = np.sqrt(u**2 + v**2)

wu = ncin.wu[:,:,:]
wv = ncin.wv[:,:,:]

# normalizando os vetores
un = u/s
vn = v/s

# timestep
t = -1

# criando variaveis para plot mais organizado
xplot = formatGrid_plot(lon,'/media/danilo/Danilo/mestrado/github/masterThesis_analysis/routines/index_list.npy')
yplot = formatGrid_plot(lat,'/media/danilo/Danilo/mestrado/github/masterThesis_analysis/routines/index_list.npy')
uplot = formatGrid_plot(un[t,:,:],'/media/danilo/Danilo/mestrado/github/masterThesis_analysis/routines/index_list.npy')
vplot = formatGrid_plot(vn[t,:,:],'/media/danilo/Danilo/mestrado/github/masterThesis_analysis/routines/index_list.npy')

wuplot = formatGrid_plot(wu[t,:,:],'/media/danilo/Danilo/mestrado/github/masterThesis_analysis/routines/index_list.npy')
wvplot = formatGrid_plot(wv[t,:,:],'/media/danilo/Danilo/mestrado/github/masterThesis_analysis/routines/index_list.npy')

fig,ax = plt.subplots()
m = oceano.make_map(ax)

m.contourf(lon,lat,s[t,:,:],latlon=True,cmap=cmo.cm.speed)
m.quiver(xplot,yplot,uplot,vplot,scale=70,width=0.001,headwidth=8,headlength=8,alpha=0.6,latlon=True)
m.quiver(xplot,yplot,wuplot,wvplot,latlon=True,color='k', alpha=0.3, scale=300, pivot='middle', headlength=3, headaxislength=2.8)
