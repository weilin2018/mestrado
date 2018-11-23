# BULK STRATIFICATION

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
import seawater as sw


import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def bulk_stratification(ncin,ilat=[82,86,90,94,98,101],ilon=70,timestep=None):

    # importing variables
    Tsurf = ncin.temp[timestep,0,ilat,:ilon]
    Ssurf = ncin.salt[timestep,0,ilat,:ilon]
    Tbott =ncin.temp[timestep,-1,ilat,:ilon]
    Sbott =ncin.temp[timestep,-1,ilat,:ilon]

    # calculating density based on seawater.dens0
    Dsurf = sw.dens0(Ssurf.values,Tsurf.values)
    Dbott = sw.dens0(Sbott.values,Tbott.values)

    # calculating bulk stratification by subtracting the surface water density from near-bottom
    B = Dbott - Dsurf

    return B

def calc_distance(ncin):
    # importing coordinates (already 2D)
    lon = ncin.lon.values
    lat = ncin.lat.values
    # cleaning non-values
    lon[lon == 0.] = np.nan
    lat[lat == 0.] = np.nan

    # calculating distance
    dist,angle = sw.dist(lon,lat)

    return dist





##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
SAVE_FIG = BASE_DIR + 'masterThesis_analysis/figures/experiments_outputs/temperature/crossSection_EA5/'
fname = glob.glob(DATA_DIR+"*.cdf")

# select which experiment you want to plot:
exp = 'EC1.cdf'

for f in fname:
    if exp in f:
        experiment = f

ncin = xr.open_dataset(f)

ilat=[82,86,90,94,98,101]
ilon=70

# visualize locations
fig,ax = plt.subplots()
m = oceano.make_map(ax,ulat=-23.122243,llat=-25.129797,ulon=-44.089583,llon=-45.6101478,resolution='f')
m.scatter(ncin.lon[ilats,:ilon].values,ncin.lat[ilats,:ilon].values,latlon=True)
m.contour(ncin.lon.values,ncin.lat.values,ncin.depth.values,levels=[100,150,200],colors=('k'),linestyles=('--'),latlon=True)

# plotting for the first day of the blocking event
B = bulk_stratification(ncin,timestep=46)
dist = calc_distance(ncin)
dist = np.tile(dist[ilats[-1],:ilon],(6,1))

fig,ax = plt.subplots()
ax.scatter(dist,B)
