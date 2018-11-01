#-*-coding;utf-8-*-
"""

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
import cmocean as cmo
import gsw

# pacotes para minimap
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

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
# if BASE_DIR.split("/")[2] == 'tparente':
#     DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output_modelo/exp03_variables/')
#     fname = glob.glob(DATA_DIR+"*.cdf")
# else:
#     DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
#     fname = glob.glob(DATA_DIR+"*.cdf")
plt.ion()

DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
fname = glob.glob(DATA_DIR+"*.cdf")

# select which experiment you want to plot:
exp = 'EC1.cdf'
SAVE_FIG = BASE_DIR + 'masterThesis_analysis/figures/experiments_outputs/temperature/crossSection_EA1/'

for f in fname:
    if exp in f:
        experiment = f

fname = experiment
ncin = xr.open_dataset(fname)
# extract variables
lon,lat = ncin['lon'].values, ncin['lat'].values
lon[lon == 0.] = np.nan
lat[lat == 0.] = np.nan
depth = ncin['depth'].values
sigma = ncin['sigma'].values

os.system('clear')
loc = input('Selecione a regiao para plotar a secao vertical: [0] - Ubatuba, [1] - Santos e [2] Cananeia')

if loc == 0:
    ind = 99# inor = 99
    loc = 'Ubatuba'
if loc == 1:
    ind = 28# icen = 28
    loc = 'Santos'
if loc == 2:
    ind = 19 # isul = 19
    loc = 'Cananeia'

# definindo algumas variaveis
eixoX_m = gsw.distance(lon[ind,:],lat[ind,:])
# eixoX_m = np.squeeze(eixoX_m)
# eixoX_m = np.append(eixoX_m,np.nan)
x,prof,sig = oceano.create_newDepth(lon,depth,sigma,ind)      # create new depth
liminf,limsup = 5,83               # limits with non-nan values

fig,axes = plt.subplots(ncols=2,figsize=(13,5))
axes[0].plot(lon[ind,liminf:limsup],-depth[ind,liminf:limsup],'k')
axes[0].fill_between(lon[ind,liminf:limsup], -200, -depth[ind,liminf:limsup],color='#c0c0c0')
axes[0].set_ylim([-200,1])
if ind == 99:
    # ubatuba
    axes[0].set_xlim([-45.,-44.4])

axes[0].margins(0)
axes[0].set_ylabel(u'Depth [m]',fontsize=18)
axes[0].set_title('Experimento Controle',fontsize=18)

axes[1].plot(lon[ind,liminf:limsup],-depth[ind,liminf:limsup],'k')
axes[1].fill_between(lon[ind,liminf:limsup], -200, -depth[ind,liminf:limsup],color='#c0c0c0')
axes[1].set_ylim([-200,1])
if ind == 99:
    # ubatuba
    axes[1].set_xlim([-45.,-44.4])

axes[1].margins(0)
axes[1].set_ylabel(u'Depth [m]',fontsize=18)
axes[1].set_title(u'Experimento Anômalo',fontsize=18)

# defining the begin and the end to plot
tBegin = 0 # climatologic position
tFinal = 304 # final do evento em estudo


# plotting Experimento Controle (left)

# plotting the evolution of the isotherm line
# for t in np.arange(46,303,1):
#     temp = ncin.temp[t,:,ind,:]
#     cs  = axes[0].contour(x[:,liminf:limsup],sig[:,liminf:limsup],temp[:,liminf:limsup],levels=[18.],colors=('#c0c0c0'),linestyles=('--'))
#     plt.pause(0.1)

temp = ncin.temp[tBegin,:,ind,:]
cs  = axes[0].contour(x[:,liminf:limsup],sig[:,liminf:limsup],temp[:,liminf:limsup],levels=[18.],colors=('r'),linestyles=('--'))

temp = ncin.temp[tFinal,:,ind,:]
cs  = axes[0].contour(x[:,liminf:limsup],sig[:,liminf:limsup],temp[:,liminf:limsup],levels=[18.],colors=('g'),linestyles=('--'))

# plotting Experimento Anomalo
ncin = xr.open_dataset(fname.replace('EC1','EA1'))

# plotting the evolution of the isotherm line
# for t in np.arange(46,303,1):
#     temp = ncin.temp[t,:,ind,:]
#     cs  = axes[1].contour(x[:,liminf:limsup],sig[:,liminf:limsup],temp[:,liminf:limsup],levels=[18.],colors=('#c0c0c0'),linestyles=('--'))
#     plt.pause(0.1)

temp = ncin.temp[tBegin,:,ind,:]
cs  = axes[1].contour(x[:,liminf:limsup],sig[:,liminf:limsup],temp[:,liminf:limsup],levels=[18.],colors=('r'),linestyles=('--'))

temp = ncin.temp[tFinal,:,ind,:]
cs  = axes[1].contour(x[:,liminf:limsup],sig[:,liminf:limsup],temp[:,liminf:limsup],levels=[18.],colors=('g'),linestyles=('--'))

plt.suptitle(loc,fontsize=24)
