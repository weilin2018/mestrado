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
import seawater as sw

# pacotes para minimap
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

import matplotlib
# matplotlib.style.use('ggplot')
matplotlib.use('PS')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def struc(indexes):

    fig,axes = plt.subplots(nrows=3,ncols=3,figsize=(25.4/2.54,20./2.54))

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


def create_Structure_3(ncin,indexes):
    lon,lat = ncin['lon'].values, ncin['lat'].values
    lon[lon == 0.] = np.nan
    lat[lat == 0.] = np.nan
    depth = ncin['depth'].values
    sigma = ncin['sigma'].values
    dx    = ncin['h1'].values

    fig,axes = plt.subplots(nrows=3,ncols=3,figsize=(25.4/2.54, 20/2.54))

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
        dist,inds = calcDistance(x,sigma,dx[ind,:])
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

def calcDistance(x,sig,dx,xx=110):

    inds = np.where(np.isnan(x[0,:])) # aonde e nan
    lats = np.ones([xx])*11

    # removendo nan's
    x =  np.delete(x[0,:],inds)
    lats=np.delete(lats,inds)

    # dist2 = np.cumsum(np.append(0,sw.dist(lats,x)[0]))
    # dist = np.tile(dist2,(len(sig),1))

    dist      = np.cumsum(dx)/1000
    dist      = np.tile(dist,(len(sig),1))
    dist      = np.delete(dist,inds,axis=1)

    return dist,inds

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


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
#plt.ion()

# configurações do plot
figsize = (17.4/2.54, 10/2.54)

DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
fname = glob.glob(DATA_DIR+"*.cdf")

# select which experiment you want to plot:
exp = 'EA1.cdf'
# SAVE_FIG = BASE_DIR + 'masterThesis_analysis/figures/experiments_outputs/temperature/crossSection_EC2/'

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
h1    = ncin['h1'].values
temp  = ncin.temp.values
# index para as latitudes das seções
indexes = [99,28,19]

fig,axes = struc(indexes)
for i in range(3):
    for j in range(3):
        axes[j,i].set_xlim([0,150])
        axes[j,i].set_ylim([-100,0])

axes[0,0].set_title(u'Climatologia',fontsize=8)
axes[0,1].set_title(u'Experimento EA1',fontsize=8)
axes[0,2].set_title(u'Experimento EA2',fontsize=8)

title = u'Seção vertical de temperatura climatológica e nos experimentos anômalos,'\
      + u'\n com destaque para a isoterma de 18' + r'$^o$C'
plt.suptitle(title,fontsize=10)

# defining the begin and the end to plot
tBegin = 46 # climatologic position
tFinal = 410 # final do evento em estudo

os.system('clear')
print('# ----- PLOTTING CLIMATOLOGY ----- #')

for ind in indexes:
    if ind == 99:
        axesInd = 0
        cutLon = 80
    if ind == 28:
        axesInd = 1
        cutLon = 80
    if ind == 19:
        axesInd = 2
        cutLon = 90

    T = np.nanmean(temp[:3,:,ind,:cutLon],axis=0)

    Tplot,ndist,ndepth,dist2,sig,depth = crossSection_optimized(lon,depth,sigma,h1,T,cutLon,ind=ind)

    xgrid,zgrid = np.meshgrid(ndist,ndepth)
    contours = np.arange(14,35,0.5)

    cf1  = axes[axesInd,0].contourf(xgrid,-zgrid,Tplot,contours,cmap=cmo.cm.thermal,extend='max')
    cs   = axes[axesInd,0].contour(xgrid,-zgrid,Tplot,levels=[18.],colors=('k'),linestyles=('--'))
    axes[axesInd,0].fill_between(dist2[-1,:], -100, -depth[ind,:],color='#c0c0c0')
    axes[axesInd,0].plot(dist2[-1,:],sig[-1,:],'k')

    for c in cf1.collections:
        c.set_edgecolor('face')
        c.set_linewidth(0.00000000001)

print('# ----- PLOTTING ANOMALY 1 ----- #')

for ind in indexes:
    if ind == 99:
        axesInd = 0
        cutLon = 80
    if ind == 28:
        axesInd = 1
        cutLon = 80
    if ind == 19:
        axesInd = 2
        cutLon = 80

    T = np.nanmean(temp[tBegin-3:tBegin+3,:,ind,:cutLon],axis=0)

    Tplot,ndist,ndepth,dist2,sig,depth = crossSection_optimized(lon,depth,sigma,h1,T,cutLon,ind=ind)
    # begin: 18 isotherm position
    cs   = axes[axesInd,1].contour(xgrid,-zgrid,Tplot,levels=[18.],colors=('k'),linestyles=('--'))
    # final position and vertical structure
    T = np.nanmean(temp[tFinal-3:tFinal+3,:,ind,:cutLon],axis=0)
    Tplot,ndist,ndepth,dist2,sig,depth = crossSection_optimized(lon,depth,sigma,h1,T,cutLon,ind=ind)
    cs   = axes[axesInd,1].contour(xgrid,-zgrid,Tplot,levels=[18.],colors=('w'),linestyles=('--'))
    cf2  = axes[axesInd,1].contourf(xgrid,-zgrid,Tplot,contours,cmap=cmo.cm.thermal,extend='max')

    axes[axesInd,1].fill_between(dist2[-1,:], -100, -depth[ind,:],color='#c0c0c0')
    axes[axesInd,1].plot(dist2[-1,:],sig[-1,:],'k')

    for c in cf2.collections:
        c.set_edgecolor('face')
        c.set_linewidth(0.00000000001)


print('# ----- PLOTTING ANOMALY 2 ----- #')
ncin = xr.open_dataset(fname.replace('1','2'))
temp  = ncin.temp.values

for ind in indexes:
    if ind == 99:
        axesInd = 0
        cutLon = 80
    if ind == 28:
        axesInd = 1
        cutLon = 80
    if ind == 19:
        axesInd = 2
        cutLon = 80

    T = np.nanmean(temp[tBegin-3:tBegin+3,:,ind,:cutLon],axis=0)

    Tplot,ndist,ndepth,dist2,sig,depth = crossSection_optimized(lon,depth,sigma,h1,T,cutLon,ind=ind)
    # begin: 18 isotherm position
    cs   = axes[axesInd,2].contour(xgrid,-zgrid,Tplot,levels=[18.],colors=('k'),linestyles=('--'))
    # final position and vertical structure
    T = np.nanmean(temp[tFinal-3:tFinal+3,:,ind,:cutLon],axis=0)
    Tplot,ndist,ndepth,dist2,sig,depth = crossSection_optimized(lon,depth,sigma,h1,T,cutLon,ind=ind)
    cs   = axes[axesInd,2].contour(xgrid,-zgrid,Tplot,levels=[18.],colors=('w'),linestyles=('--'))
    cf3  = axes[axesInd,2].contourf(xgrid,-zgrid,Tplot,contours,cmap=cmo.cm.thermal,extend='max')

    axes[axesInd,2].fill_between(dist2[-1,:], -100, -depth[ind,:],color='#c0c0c0')
    axes[axesInd,2].plot(dist2[-1,:],sig[-1,:],'k')

    for c in cf3.collections:
        c.set_edgecolor('face')
        c.set_linewidth(0.00000000001)

plt.tight_layout()
plt.subplots_adjust(top=0.905,bottom=0.059,left=0.068,right=0.987,hspace=0.11,wspace=0.068)

plt.savefig(BASE_DIR+ 'masterThesis_analysis/figures/experiments_outputs/temperature/secao3x3.eps')
