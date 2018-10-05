# using my package for model's visualization

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

from modelVisualization.interface import Experiment

exp = 'hotstarts/run_coldProg_summer'
fname = '/media/danilo/Danilo/mestrado/ventopcse/output/%s.cdf'%(exp)

# instanciating Experiment
warm = Experiment(fname,timeStart='2013-12-21',timeEnd='2013-12-30',region='pcse')

# plot the last timestep
t = -1
plt.ion()

# SST
contours = np.arange(17,28,0.1)
fig,ax = plt.subplots()
m = oceano.make_map(ax)

cf = m.contourf(warm.lon,warm.lat,warm.ncin.temp[t,0,:,:].values,contours,latlon=True,cmap=cmo.cm.thermal)
cbar = plt.colorbar(cf)

# SSS
contours = np.arange(33,38,.2)
fig,ax = plt.subplots()
m = oceano.make_map(ax)

cf = m.contourf(warm.lon,warm.lat,warm.ncin.salt[t,0,:,:].values,contours,latlon=True,cmap=cmo.cm.haline)
cbar = plt.colorbar(cf)

# SSV - sea surface velocity
u = warm.ncin.u[t,0,:,:].values
v = warm.ncin.v[t,0,:,:].values
s = np.sqrt(u**2 + v**2)

un = u/s
vn = v/s

fig,ax = plt.subplots()
m = oceano.make_map(ax)

# cf = m.contourf(warm.lon,warm.lat,s,latlon=True,cmap=cmo.cm.speed)
cf = m.contourf(warm.lon,warm.lat,warm.ncin.elev[t,:,:].values,latlon=True,cmap='RdBu')
qv = m.quiver(warm.lon,warm.lat,u,v,latlon=True)
