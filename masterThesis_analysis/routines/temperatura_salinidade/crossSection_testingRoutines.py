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
from scipy import interpolate

import matplotlib

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################

def plot1(fname,ind=99):
    grid_x,grid_z,Tinterp,stdl,dist,dist2,sig,h1 = oceano.crossSection(fname,ind=ind)

    ind = 99 # section
    contours = np.arange(14,35,0.5)
    Tplot = np.nanmean(Tinterp[:,:,ind,:],axis=0)
    ndepth = np.linspace(0,100,100)

    newDist,grid_x2,grid_z2,Tinterp = oceano.interpDistance(Tplot,h1,80,ndepth)

    fig,ax = plt.subplots(ncols=2,figsize=(16./2.54,8./2.54),sharey=True)
    ax[0].contourf(grid_x,-grid_z,Tplot,contours,cmap=cmo.cm.thermal)
    ax[0].plot(dist2[-1,:],sig[-1,:],'k')
    ax[0].set_xlim([0,100])
    ax[0].set_ylim([-100,0])
    ax[0].contour(grid_x,-grid_z,Tplot,levels=[18.],colors=('k'))
    ax[0].set_title('Coordenada Z, Distancia irregular',fontsize=10)
    # plotting grid
    # ax[0].plot(grid_x,-grid_z,'k',alpha=.3);
    # ax[0].plot(grid_x.T,-grid_z.T,'k',alpha=.3);

    ax[1].contourf(grid_x2,-grid_z2,Tinterp,contours,cmap=cmo.cm.thermal)
    ax[1].plot(dist2[-1,:],sig[-1,:],'k')
    ax[1].set_xlim([0,100])
    ax[1].set_ylim([-100,0])
    ax[1].contour(grid_x,-grid_z,Tplot,levels=[18.],colors=('k'))
    ax[1].set_title('Coordenada Z, Distancia Regular',fontsize=10)

    plt.tight_layout()
    plt.subplots_adjust(top=0.911,bottom=0.103,left=0.08,right=0.98,hspace=0.18,wspace=0.084)

def plot2(ax,fname,horizResolution=10000,vertResolution=100,depRef=100,ind=99):
    ncin = xr.open_dataset(fname)

    # extract grid and other general variables
    # important: lat and lon already gridded
    lon   = ncin.lon.values
    lon[lon == 0.] = np.nan
    depth = ncin.depth.values
    sigma = ncin.sigma.values
    h1    = ncin['h1'].values
    temp = ncin.temp[296:306,:,:,:]
    T = np.nanmean(temp[:,:,ind,:],axis=0)

    # Tplot,ndist,ndepth,dist2,sig,depth = crossSection_optimized(lon,depth,sigma,h1,T,cutLon,ind=ind)
    Tplot,ndist,ndepth,dist2,sig,depth = oceano.crossSection_optimized(lon,depth,sigma,h1,T,horizResolution=horizResolution,vertResolution=vertResolution,depRef=depRef,ind=ind)

    xgrid,zgrid = np.meshgrid(ndist,ndepth)

    contours = np.arange(14,35,0.5)

    ax.set_title('ind = %s'%(str(ind)))
    ax.contourf(xgrid,-zgrid,Tplot,contours,cmap=cmo.cm.thermal,extend='both')
    ax.fill_between(dist2[-1,:], -depRef, sig[-1,:],color='#c0c0c0')
    ax.plot(dist2[-1,:],sig[-1,:],'k')
    ax.contour(xgrid,-zgrid,Tplot,levels=[18.],colors=('k'))
    ax.set_xlim([0,150000])
    ax.set_ylim([-depRef,0])

    plt.tight_layout()
    plt.subplots_adjust(top=0.911,bottom=0.103,left=0.08,right=0.98,hspace=0.18,wspace=0.084)

    return ax

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
SAVE_DIR = BASE_DIR + 'masterThesis_analysis/figures/experiments_outputs/temperature/'
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
plt.ion()

# select which experiment you want to plot:
fname = DATA_DIR + 'EA2.cdf'

# mantendo uma consistencia entre a profundidade maxima utilizada e a resolucao vertical
# o resultado sera melhor.
fig,ax = plt.subplots()
ax = plot2(ax,fname,ind=99,depRef=800,vertResolution=800)



labels = [item.get_text() for item in ax.get_xticklabels()]

ax.set_xticklabels(labels)
# plot1(fname,ind=18)
