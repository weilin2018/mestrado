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
import seawater as sw

import matplotlib
matplotlib.use('PS')

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
            axes[axesInd,0].set_xlabel(u'Distância [km]',fontsize=8)
            axes[axesInd,1].set_xlabel(u'Distância [km]',fontsize=8)

        x,prof,sig = oceano.create_newDepth(lon,depth,sigma,ind)      # create new depth
        dist,inds = calcDistance(x,sigma)
        liminf,limsup = 5,83               # limits with non-nan values

        d = np.delete(depth,inds,axis=1)

        axes[axesInd,0].plot(dist[0,:],-d[ind,:],'k')
        axes[axesInd,0].fill_between(dist[0,:], -200, -d[ind,:],color='#c0c0c0')
        axes[axesInd,0].set_ylim([-200,1])
        axes[axesInd,0].set_xlim([0,160])

        # if ind == 99:
        #     # ubatuba
        #     axes[axesInd,0].set_xlim([-45.,-44.4])

        axes[axesInd,0].margins(0)
        axes[axesInd,0].set_ylabel(u'Profundidade [m]',fontsize=8)

        axes[axesInd,1].plot(dist[0,:],-d[ind,:],'k')
        axes[axesInd,1].fill_between(dist[0,:], -200, -d[ind,:],color='#c0c0c0')
        axes[axesInd,1].set_ylim([-200,1])
        axes[axesInd,1].set_xlim([0,160])

        # if ind == 99:
        #     # ubatuba
        #     axes[axesInd,1].set_xlim([-45.,-44.4])

        axes[axesInd,1].margins(0)
        # axes[axesInd,1].set_ylabel(u'Profundidade [m]',fontsize=18)

    return fig,axes

def create_Structure_3(ncin,indexes):
    lon,lat = ncin['lon'].values, ncin['lat'].values
    lon[lon == 0.] = np.nan
    lat[lat == 0.] = np.nan
    depth = ncin['depth'].values
    sigma = ncin['sigma'].values

    fig,axes = plt.subplots(nrows=3,ncols=3,figsize=(19.4/2.54, 20/2.54))

    for ind in indexes:
        if ind == 99:
            axesInd = 0
            axes[axesInd,0].set_title('Climatologia',fontsize=8)
            axes[axesInd,1].set_title(u'EA1',fontsize=8)
            axes[axesInd,2].set_title(u'EA2',fontsize=8)
        if ind == 28:
            axesInd = 1
        if ind == 19:
            axesInd = 2
            axes[axesInd,0].set_xlabel(u'Distância [km]',fontsize=8)
            axes[axesInd,1].set_xlabel(u'Distância [km]',fontsize=8)
            axes[axesInd,2].set_xlabel(u'Distância [km]',fontsize=8)

        x,prof,sig = oceano.create_newDepth(lon,depth,sigma,ind)      # create new depth
        dist,inds = calcDistance(x,sigma)
        liminf,limsup = 5,83               # limits with non-nan values
        d = np.delete(depth,inds,axis=1)

        # configurating axes (labels, sizes, limits)
        axes[axesInd,0].plot(dist[0,:],-d[ind,:],'k')
        axes[axesInd,0].fill_between(dist[0,:], -100, -d[ind,:],color='#c0c0c0')
        axes[axesInd,0].set_ylim([-100,1])
        axes[axesInd,0].set_xlim([0,100])
        axes[axesInd,0].margins(0)
        axes[axesInd,0].set_ylabel(u'Profundidade [m]',fontsize=8)

        axes[axesInd,1].plot(dist[0,:],-d[ind,:],'k')
        axes[axesInd,1].fill_between(dist[0,:], -100, -d[ind,:],color='#c0c0c0')
        axes[axesInd,1].set_ylim([-100,1])
        axes[axesInd,1].set_xlim([0,100])

        axes[axesInd,2].plot(dist[0,:],-d[ind,:],'k')
        axes[axesInd,2].fill_between(dist[0,:], -100, -d[ind,:],color='#c0c0c0')
        axes[axesInd,2].set_ylim([-100,1])
        axes[axesInd,2].set_xlim([0,100])

        # hiding ticks labels
        axes[0,0].xaxis.set_major_formatter(plt.NullFormatter())

        axes[0,1].xaxis.set_major_formatter(plt.NullFormatter())
        axes[0,1].yaxis.set_major_formatter(plt.NullFormatter())

        axes[0,2].xaxis.set_major_formatter(plt.NullFormatter())
        axes[0,2].yaxis.set_major_formatter(plt.NullFormatter())

        axes[1,0].xaxis.set_major_formatter(plt.NullFormatter())

        axes[1,1].xaxis.set_major_formatter(plt.NullFormatter())
        axes[1,1].yaxis.set_major_formatter(plt.NullFormatter())

        axes[1,2].xaxis.set_major_formatter(plt.NullFormatter())
        axes[1,2].yaxis.set_major_formatter(plt.NullFormatter())

        axes[2,1].yaxis.set_major_formatter(plt.NullFormatter())
        axes[2,2].yaxis.set_major_formatter(plt.NullFormatter())

    return fig,axes

def calcDistance(x,sig,xx=110):

    inds = np.where(np.isnan(x[0,:])) # aonde e nan
    lats = np.ones([xx])*11

    # removendo nan's
    x =  np.delete(x[0,:],inds)
    lats=np.delete(lats,inds)

    dist2 = np.cumsum(np.append(0,sw.dist(lats,x)[0]))
    dist = np.tile(dist2,(len(sig),1))

    return dist,inds

def plotSection_experiment(fname):
    ncin = xr.open_dataset(fname)
    # extracting general variables
    lon,lat = ncin['lon'].values, ncin['lat'].values
    lon[lon == 0.] = np.nan
    lat[lat == 0.] = np.nan
    depth = ncin['depth'].values
    sigma = ncin['sigma'].values

    # setting latitude index for each section
    indexes = [99,28,19]
    # setting time indexes to create mean
    tBegin = 46
    tFinal = 303
    # create figure's structure
    fig,axes = create_Structure(ncin,indexes)
    # create contourf resolution
    contours= np.arange(10,30,0.1)
    # axis adjustments
    for i in range(2):
        for j in range(3):
            axes[j,i].set_xlim([0,100])
            axes[j,i].set_ylim([-100,0])
    # inserting suptitle
    title = u'Posição inicial (vermelho) e final (verde) da isoterma de 18' + r'$^o$C'
    plt.suptitle(title,fontsize=10)

    # plotting control experiment
    for ind in indexes:
        if ind == 99:
            axesInd = 0
        if ind == 28:
            axesInd = 1
        if ind == 19:
            axesInd = 2

        # definindo algumas variaveis
        x,prof,sig = oceano.create_newDepth(lon,depth,sigma,ind)      # create new depth
        liminf,limsup = 5,83               # limits with non-nan values
        dist,inds = calcDistance(x,sigma)

        d       = np.delete(depth[ind,:],inds)
        newSig  = np.delete(sig,inds,axis=1)

        temp = ncin.temp[tBegin:tFinal,:,ind,:]
        temp = np.nanmean(temp,axis=0)
        t    = np.delete(temp,inds,axis=1)
        cf   = axes[axesInd,0].contourf(dist,newSig,t,contours,cmap=cmo.cm.thermal,extend='max')
        cs   = axes[axesInd,0].contour(dist,newSig,t,levels=[18.],colors=('k'),linestyles=('--'))

    # reading anomalous experiment
    ncin = xr.open_dataset(fname.replace('EC','EA'))

    # plotting anomalous experiment
    for ind in indexes:
        if ind == 99:
            axesInd = 0
        if ind == 28:
            axesInd = 1
        if ind == 19:
            axesInd = 2

        # definindo algumas variaveis
        x,prof,sig = oceano.create_newDepth(lon,depth,sigma,ind)      # create new depth
        liminf,limsup = 5,83               # limits with non-nan values
        dist,inds = calcDistance(x,sigma)

        d       = np.delete(depth[ind,:],inds)
        newSig  = np.delete(sig,inds,axis=1)

        temp = ncin.temp[tBegin:tFinal,:,ind,:]
        temp = np.nanmean(temp,axis=0)
        t    = np.delete(temp,inds,axis=1)
        cf   = axes[axesInd,1].contourf(dist,newSig,t,contours,cmap=cmo.cm.thermal,extend='max')
        cs   = axes[axesInd,1].contour(dist,newSig,t,levels=[18.],colors=('k'),linestyles=('--'))

def plotsection_climatology(fname,timestep=0):
    ncin = xr.open_dataset(fname)
    # extracting general variables
    lon,lat = ncin['lon'].values, ncin['lat'].values
    lon[lon == 0.] = np.nan
    lat[lat == 0.] = np.nan
    depth = ncin['depth'].values
    sigma = ncin['sigma'].values

    # setting latitude index for each section
    indexes = [99,27,19]

    # extracting temperature
    temp = ncin.temp[timestep,:,:,:]

    # create figure's structure
    fig,axes = create_Structure_3(ncin,indexes)
    # create contourf resolution
    contours= np.arange(10,30,0.5)
    # axis adjustments
    for i in range(3):
        for j in range(3):
            axes[j,i].set_xlim([0,100])
            axes[j,i].set_ylim([-100,0])

    # subplot configuration
    # axes[0,0].set_title('Climatologia',fontsize=8)
    # axes[0,1].set_title('Experimento %s'%(fname.split('/')[-1][:3]),fontsize=8)
    # inserting suptitle
    # title = u'Posição inicial (vermelho) e final (verde) da isoterma de 18' + r'$^o$C'
    title = u'Seção vertical de temperatura climatológica e nos experimentos anômalos,'\
          + u'\n com destaque para a isoterma de 18' + r'$^o$C'
    plt.suptitle(title,fontsize=10)

    # plotting control experiment
    for ind in indexes:
        if ind == 99:
            axesInd = 0
        if ind == 27:
            axesInd = 1
        if ind == 19:
            axesInd = 2

        # definindo algumas variaveis
        x,prof,sig = oceano.create_newDepth(lon,depth,sigma,ind)      # create new depth
        liminf,limsup = 5,83               # limits with non-nan values
        dist,inds = calcDistance(x,sigma)

        d       = np.delete(depth[ind,:],inds)
        newSig  = np.delete(sig,inds,axis=1)

        temp = ncin.temp[timestep,:,ind,:]
        t    = np.delete(temp,inds,axis=1)
        cf1  = axes[axesInd,0].contourf(dist,newSig,t,contours,cmap=cmo.cm.thermal,extend='max')
        cs   = axes[axesInd,0].contour(dist,newSig,t,levels=[18.],colors=('k'),linestyles=('--'))

        for c in cf1.collections:
            c.set_edgecolor('face')
            c.set_linewidth(0.00000000001)

    # plotting anomalous period - EA1
    for ind in indexes:
        if ind == 99:
            axesInd = 0
        if ind == 27:
            axesInd = 1
        if ind == 19:
            axesInd = 2

        # definindo algumas variaveis
        x,prof,sig = oceano.create_newDepth(lon,depth,sigma,ind)      # create new depth
        liminf,limsup = 5,83               # limits with non-nan values
        dist,inds = calcDistance(x,sigma)

        d       = np.delete(depth[ind,:],inds)
        newSig  = np.delete(sig,inds,axis=1)

        temp = ncin.temp[46:303,:,ind,:]
        temp = np.nanmean(temp,axis=0)
        t    = np.delete(temp,inds,axis=1)
        cf2  = axes[axesInd,1].contourf(dist,newSig,t,contours,cmap=cmo.cm.thermal,rasterized=True,extend='max')
        cs   = axes[axesInd,1].contour(dist,newSig,t,levels=[18.],colors=('k'),linestyles=('--'))

        for c in cf2.collections:
            c.set_edgecolor('face')
            c.set_linewidth(0.00000000001)

    # plotting anomalous period - EA2
    if fname.split('/')[-1][:3] == 'EA1':
        ncin = xr.open_dataset(fname.replace('1','2'))

    for ind in indexes:
        if ind == 99:
            axesInd = 0
        if ind == 27:
            axesInd = 1
        if ind == 19:
            axesInd = 2

        # definindo algumas variaveis
        x,prof,sig = oceano.create_newDepth(lon,depth,sigma,ind)      # create new depth
        liminf,limsup = 5,83               # limits with non-nan values
        dist,inds = calcDistance(x,sigma)

        d       = np.delete(depth[ind,:],inds)
        newSig  = np.delete(sig,inds,axis=1)

        temp = ncin.temp[46:303,:,ind,:]
        temp = np.nanmean(temp,axis=0)
        t    = np.delete(temp,inds,axis=1)
        cf3  = axes[axesInd,2].contourf(dist,newSig,t,contours,cmap=cmo.cm.thermal,rasterized=True,extend='max')
        cs   = axes[axesInd,2].contour(dist,newSig,t,levels=[18.],colors=('k'),linestyles=('--'))

        for c in cf3.collections:
            c.set_edgecolor('face')
            c.set_linewidth(0.00000000001)

    # matplotib trick to remove white thin lines when saving contourf in pdf
    for c in cf1.collections:
        c.set_edgecolor('face')
        c.set_linewidth(0.00000000001)

    for c in cf2.collections:
        c.set_edgecolor('face')
        c.set_linewidth(0.00000000001)

    for c in cf3.collections:
        c.set_edgecolor('face')
        c.set_linewidth(0.00000000001)

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
SAVE_DIR = BASE_DIR + 'masterThesis_analysis/figures/experiments_outputs/temperature/'
#plt.ion()

# configurações do plot
figsize = (17.4/2.54, 10/2.54)

DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
fname = glob.glob(DATA_DIR+"*.cdf")

# select which experiment you want to plot:
exp = 'EA1.cdf'

for f in fname:
    if exp in f:
        experiment = f

fname = experiment

plotsection_climatology(fname)

plt.subplots_adjust(top=0.915,bottom=0.05,left=0.09,right=0.98,hspace=0.09,wspace=0.13)
plt.savefig(SAVE_DIR + 'temperatura_secao.eps')
