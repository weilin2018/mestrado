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
def create_Structure(ncin,indexes):
    lon,lat = ncin['lon'].values, ncin['lat'].values
    lon[lon == 0.] = np.nan
    lat[lat == 0.] = np.nan
    depth = ncin['depth'].values
    sigma = ncin['sigma'].values

    liminf,limsup = 5,100               # limits with non-nan values

    fig,axes = plt.subplots(nrows=3,ncols=2,figsize=(17.4/2.54, 18/2.54))

    for ind in indexes:
        if ind == 99:
            axesInd = 0
            axes[axesInd,0].set_title('Experimento Controle',fontsize=8)
            axes[axesInd,1].set_title(u'Experimento Anômalo',fontsize=8)
        if ind == 28:
            axesInd = 1
        if ind == 19:
            axesInd = 2
            axes[axesInd,0].set_xlabel(u'Longitude',fontsize=8)
            axes[axesInd,1].set_xlabel(u'Longitude',fontsize=8)

        # definindo algumas variaveis
        eixoX_m = gsw.distance(lon[ind,:],lat[ind,:])
        # eixoX_m = np.squeeze(eixoX_m)
        # eixoX_m = np.append(eixoX_m,np.nan)
        x,prof,sig = oceano.create_newDepth(lon,depth,sigma,ind)      # create new depth

        axes[axesInd,0].plot(lon[ind,liminf:limsup],-depth[ind,liminf:limsup],'k')
        axes[axesInd,0].fill_between(lon[ind,liminf:limsup], -200, -depth[ind,liminf:limsup],color='#c0c0c0')
        axes[axesInd,0].set_ylim([-200,1])
        # if ind == 99:
        #     # ubatuba
        #     axes[axesInd,0].set_xlim([-45.,-44.4])

        axes[axesInd,0].margins(0)
        axes[axesInd,0].set_ylabel(u'Profundidade [m]',fontsize=8)

        axes[axesInd,1].plot(lon[ind,liminf:limsup],-depth[ind,liminf:limsup],'k')
        axes[axesInd,1].fill_between(lon[ind,liminf:limsup], -200, -depth[ind,liminf:limsup],color='#c0c0c0')
        axes[axesInd,1].set_ylim([-200,1])

        axes[axesInd,1].margins(0)
        # axes[axesInd,1].set_ylabel(u'Profundidade [m]',fontsize=18)



    return fig,axes

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
plt.ion()

# configurações do plot
figsize = (17.4/2.54, 10/2.54)

DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
fname = glob.glob(DATA_DIR+"*.cdf")

# select which experiment you want to plot:
exp = 'EC1_10DTI.cdf'
SAVE_FIG = BASE_DIR + 'masterThesis_analysis/figures/experiments_outputs/temperature/crossSection_EA1/'

for f in fname:
    if exp in f:
        experiment = f

fname = experiment
ncin = xr.open_dataset(fname)

lon,lat = ncin['lon'].values, ncin['lat'].values
lon[lon == 0.] = np.nan
lat[lat == 0.] = np.nan
depth = ncin['depth'].values
sigma = ncin['sigma'].values

# index para as latitudes das seções
indexes = [99,28,19]

fig,axes = create_Structure(ncin,indexes)

title = u'Posição inicial (vermelho) e final (verde) da isohalina de 36'
plt.suptitle(title,fontsize=10)

# defining the begin and the end to plot
tBegin = 0 # climatologic position
tFinal = 303 # final do evento em estudo

######## Experimento Controle!!!!
for ind in indexes:
    if ind == 99:
        axesInd = 0
    if ind == 28:
        axesInd = 1
    if ind == 19:
        axesInd = 2

    # definindo algumas variaveis
    eixoX_m = gsw.distance(lon[ind,:],lat[ind,:])
    # eixoX_m = np.squeeze(eixoX_m)
    # eixoX_m = np.append(eixoX_m,np.nan)
    x,prof,sig = oceano.create_newDepth(lon,depth,sigma,ind)      # create new depth
    liminf,limsup = 5,100               # limits with non-nan values

    salt = ncin.salt[tBegin,:,ind,:]
    cs  = axes[axesInd,0].contour(x[:,liminf:limsup],sig[:,liminf:limsup],salt[:,liminf:limsup],levels=[36.],colors=('r'),linestyles=('--'))
    salt = ncin.salt[tFinal,:,ind,:]
    cs  = axes[axesInd,0].contour(x[:,liminf:limsup],sig[:,liminf:limsup],salt[:,liminf:limsup],levels=[36.],colors=('g'),linestyles=('--'))

######## Experimento Anomalo!!!!
ncin = xr.open_dataset(fname.replace('EC1','EA1'))

for ind in indexes:
    if ind == 99:
        axesInd = 0
    if ind == 28:
        axesInd = 1
    if ind == 19:
        axesInd = 2

    # definindo algumas variaveis
    eixoX_m = gsw.distance(lon[ind,:],lat[ind,:])
    # eixoX_m = np.squeeze(eixoX_m)
    # eixoX_m = np.append(eixoX_m,np.nan)
    x,prof,sig = oceano.create_newDepth(lon,depth,sigma,ind)      # create new depth
    liminf,limsup = 5,100               # limits with non-nan values


    salt = ncin.salt[tBegin,:,ind,:]
    cs  = axes[axesInd,1].contour(x[:,liminf:limsup],sig[:,liminf:limsup],salt[:,liminf:limsup],levels=[36.],colors=('r'),linestyles=('--'))
    salt = ncin.salt[tFinal,:,ind,:]
    cs  = axes[axesInd,1].contour(x[:,liminf:limsup],sig[:,liminf:limsup],salt[:,liminf:limsup],levels=[36.],colors=('g'),linestyles=('--'))

plt.subplots_adjust(top=0.925,bottom=0.06,left=0.115,right=0.95,hspace=0.2,wspace=0.28)
plt.savefig('/media/danilo/Danilo/mestrado/github/masterThesis_analysis/figures/experiments_outputs/salinity/isohalinaPosition/Secao_All.pdf')
