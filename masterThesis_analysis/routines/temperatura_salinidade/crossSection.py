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

import scipy.ndimage

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def crossSection_optimized(lon,depth,sigma,h1,variable,cutLon,ind=99):
    # setting the new grid
    ndepth = np.linspace(0,100,100) # new vertical resolution
    ndist  = np.linspace(0,100000,10000) # new horizontal resolution

    # interpolating horizontal/distance
    newDist,grid_x,grid_z,Tinterp = oceano.interpDistance(variable,h1,cutLon,ndepth)

    # interpolatin depth, to have same len of ndist
    ndep = oceano.interpDepth(depth,h1,ind,ndist,cutLon)

    # interpolating vertically
    stdl,Tplot = oceano.interpSigma(Tinterp,sigma,ndep,100)

    # create grid to plot transect
    xgrid,zgrid = np.meshgrid(ndist/1000,ndepth) # with distance in km

    dist = np.cumsum(h1[ind,:])/1000
    x,prof,sig = oceano.create_newDepth(lon,depth,sigma,ind)
    dist2 = np.tile(dist,(21,1))

    return Tplot,ndist/1000,ndepth,dist2,sig,depth

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

def plot2(fname,ind=99,cutLon=80):
    ncin = xr.open_dataset(fname)

    # extract grid and other general variables
    # important: lat and lon already gridded
    lon   = ncin.lon.values
    lon[lon == 0.] = np.nan
    depth = ncin.depth.values
    sigma = ncin.sigma.values
    h1    = ncin['h1'].values
    temp = ncin.temp[300:303,:,:,:]
    T = np.nanmean(temp[:,:,ind,:cutLon],axis=0)

    Tplot,ndist,ndepth,dist2,sig,depth = crossSection_optimized(lon,depth,sigma,h1,T,cutLon,ind=ind)

    xgrid,zgrid = np.meshgrid(ndist,ndepth)
    contours = np.arange(14,35,0.5)

    fig,ax = plt.subplots()
    ax.contourf(xgrid,-zgrid,Tplot,contours,cmap=cmo.cm.thermal)
    ax.fill_between(dist2[-1,:], -100, -depth[ind,:],color='#c0c0c0')
    ax.plot(dist2[-1,:],sig[-1,:],'k')
    ax.set_xlim([0,100])
    ax.set_ylim([-100,0])
    ax.contour(xgrid,-zgrid,Tplot,levels=[18.],colors=('k'))

    plt.tight_layout()
    plt.subplots_adjust(top=0.911,bottom=0.103,left=0.08,right=0.98,hspace=0.18,wspace=0.084)

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
SAVE_DIR = BASE_DIR + 'masterThesis_analysis/figures/experiments_outputs/temperature/'
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
plt.ion()

# select which experiment you want to plot:
fname = DATA_DIR + 'EA1.cdf'



# plot1(fname,ind=18)
