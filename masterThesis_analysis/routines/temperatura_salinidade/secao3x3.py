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

    axes[0,0].set_ylabel('Profundidade [m]')
    axes[1,0].set_ylabel('Profundidade [m]')
    axes[2,0].set_ylabel('Profundidade [m]')

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

def search_information(ncin,ind,nstepBegin,nstepFinal,loc,var):
    # based on an ind value, return a dictionary with informations for each
    # cross section, such as location, variable position, etc.

    # set some variables
    if var == 'temp':
        sigma = -1
        value = 18.
    if var == 'salt':
        sigma = 0
        value = 36.

    # iniatilize dictionary
    info = {
        'location': loc,
        'beginPos': oceano.find_distance_of_a_value(ncin,ind,nstepBegin,sigma,var,value)[0],
        'finalPos': oceano.find_distance_of_a_value(ncin,ind,nstepFinal,sigma,var,value)[0]
    }

    return info


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
plt.ion()

# configurações do plot
figsize = (17.4/2.54, 10/2.54)
# index para as latitudes das seções
indexes = [99,28,19]
# configuracoes para as secoes verticais
horizResolution = 10000
vertResolution  = 100 # isso significa que teremos uma resolucao vertical de 1m
depRef          = 200 # profundidade de referencia para interpolacao

limiteEixoX = 300000 # limite, em metros, do eixo X para a secao vertical
contours = np.arange(5,35,0.5)

fname = DATA_DIR + "EA1.cdf"
ncin = xr.open_dataset(fname)

lon,lat = ncin['lon'].values, ncin['lat'].values
lon[lon == 0.] = np.nan
lat[lat == 0.] = np.nan
depth = ncin['depth'].values
sigma = ncin['sigma'].values
h1    = ncin['h1'].values
temp  = ncin.temp.values

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
nstepBegin = np.arange(40,52,1) # begin of anomalous period
nstepFinal = np.arange(297,303,1) # final of anomalous period

os.system('clear')
print('# ----- PLOTTING CLIMATOLOGY ----- #')

for ind in indexes:
    if ind == 99:
        axesInd = 0
        location = 'Ubatuba'
    if ind == 28:
        axesInd = 1
        location = 'Santos'
    if ind == 19:
        axesInd = 2
        location = u'Cananéia'

    T = np.nanmean(temp[:3,:,ind,:],axis=0)

    Tplot,ndist,ndepth,dist2,sig,depth = oceano.crossSection_optimized(lon,depth,sigma,h1,T,horizResolution=horizResolution,vertResolution=vertResolution,depRef=depRef,ind=ind)

    xgrid,zgrid = np.meshgrid(ndist,ndepth)

    cf1  = axes[axesInd,0].contourf(xgrid,-zgrid,Tplot,contours,cmap=cmo.cm.thermal,extend='both')
    cs   = axes[axesInd,0].contour(xgrid,-zgrid,Tplot,levels=[18.],colors=('k'),linestyles=('--'))
    axes[axesInd,0].fill_between(dist2[-1,:], -depRef, sig[-1,:],color='#c0c0c0')
    axes[axesInd,0].plot(dist2[-1,:],sig[-1,:],'k')
    axes[axesInd,0].set_xlim([0,limiteEixoX])
    axes[axesInd,0].set_ylim([-depRef,0])

    for c in cf1.collections:
        c.set_edgecolor('face')
        c.set_linewidth(0.00000000001)

    # plot text box
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    textstr = u'%s' % (location)
    axes[axesInd,0].text(0.18, 0.28, textstr, transform=axes[axesInd,0].transAxes, fontsize=8,
        va='top', ha='center',bbox=props)


print('# ----- PLOTTING ANOMALY 1 ----- #')

for ind in indexes:
    if ind == 99:
        axesInd = 0
        infos = search_information(ncin,ind,nstepBegin,nstepFinal,'Ubatuba','temp')
    if ind == 28:
        axesInd = 1
        infos = search_information(ncin,ind,nstepBegin,nstepFinal,'Santos','temp')
    if ind == 19:
        axesInd = 2
        infos = search_information(ncin,ind,nstepBegin,nstepFinal,u'Cananéia','temp')

    T = np.nanmean(temp[nstepBegin,:,ind,:],axis=0)

    Tplot,ndist,ndepth,dist2,sig,depth = oceano.crossSection_optimized(lon,depth,sigma,h1,T,horizResolution=horizResolution,vertResolution=vertResolution,depRef=depRef,ind=ind)

    xgrid,zgrid = np.meshgrid(ndist,ndepth)

    # begin: 18 isotherm position
    cs   = axes[axesInd,1].contour(xgrid,-zgrid,Tplot,levels=[18.],colors=('k'),linestyles=('--'))
    # final position and vertical structure
    T = np.nanmean(temp[nstepFinal,:,ind,:],axis=0)
    Tplot,ndist,ndepth,dist2,sig,depth = oceano.crossSection_optimized(lon,depth,sigma,h1,T,horizResolution=horizResolution,vertResolution=vertResolution,depRef=depRef,ind=ind)
    cs   = axes[axesInd,1].contour(xgrid,-zgrid,Tplot,levels=[18.],colors=('w'),linestyles=('--'))
    cf2  = axes[axesInd,1].contourf(xgrid,-zgrid,Tplot,contours,cmap=cmo.cm.thermal,extend='both')
    axes[axesInd,1].fill_between(dist2[-1,:], -depRef, sig[-1,:],color='#c0c0c0')
    axes[axesInd,1].plot(dist2[-1,:],sig[-1,:],'k')
    axes[axesInd,1].set_xlim([0,limiteEixoX])
    axes[axesInd,1].set_ylim([-depRef,0])

    # plot text box
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    textstr = u'Inicio\n %.1f km \n Final\n %.1f km' % (infos['beginPos'],infos['finalPos'])
    axes[axesInd,1].text(0.17, 0.32, textstr, transform=axes[axesInd,1].transAxes, fontsize=8,
        va='top', ha='center',bbox=props)

    for c in cf2.collections:
        c.set_edgecolor('face')
        c.set_linewidth(0.00000000001)


print('# ----- PLOTTING ANOMALY 2 ----- #')
ncin = xr.open_dataset(fname.replace('1','2'))
temp  = ncin.temp.values

for ind in indexes:
    if ind == 99:
        axesInd = 0
        infos = search_information(ncin,ind,nstepBegin,nstepFinal,'Ubatuba','temp')
    if ind == 28:
        axesInd = 1
        infos = search_information(ncin,ind,nstepBegin,nstepFinal,'Santos','temp')
    if ind == 19:
        axesInd = 2
        infos = search_information(ncin,ind,nstepBegin,nstepFinal,u'Cananéia','temp')

    T = np.nanmean(temp[nstepBegin,:,ind,:],axis=0)

    Tplot,ndist,ndepth,dist2,sig,depth = oceano.crossSection_optimized(lon,depth,sigma,h1,T,horizResolution=horizResolution,vertResolution=vertResolution,depRef=depRef,ind=ind)

    xgrid,zgrid = np.meshgrid(ndist,ndepth)

    # begin: 18 isotherm position
    cs   = axes[axesInd,2].contour(xgrid,-zgrid,Tplot,levels=[18.],colors=('k'),linestyles=('--'))
    # final position and vertical structure
    T = np.nanmean(temp[nstepFinal,:,ind,:],axis=0)
    Tplot,ndist,ndepth,dist2,sig,depth = oceano.crossSection_optimized(lon,depth,sigma,h1,T,horizResolution=horizResolution,vertResolution=vertResolution,depRef=depRef,ind=ind)
    cs   = axes[axesInd,2].contour(xgrid,-zgrid,Tplot,levels=[18.],colors=('w'),linestyles=('--'))
    cf3  = axes[axesInd,2].contourf(xgrid,-zgrid,Tplot,contours,cmap=cmo.cm.thermal,extend='both')
    axes[axesInd,2].fill_between(dist2[-1,:], -depRef, sig[-1,:],color='#c0c0c0')
    axes[axesInd,2].plot(dist2[-1,:],sig[-1,:],'k')
    axes[axesInd,2].set_xlim([0,limiteEixoX])
    axes[axesInd,2].set_ylim([-depRef,0])

    # plot text box
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    textstr = u'Inicio\n %.1f km \n Final\n %.1f km' % (infos['beginPos'],infos['finalPos'])
    axes[axesInd,2].text(0.17, 0.32, textstr, transform=axes[axesInd,2].transAxes, fontsize=8,
        va='top', ha='center',bbox=props)

    for c in cf3.collections:
        c.set_edgecolor('face')
        c.set_linewidth(0.00000000001)

# ajusta a figura antes de se ajustar os labels do eixo x, pois aumenta-se a quantidade de
# ticks no eixo quando dá tight_layout
plt.tight_layout()
plt.subplots_adjust(top=0.905,bottom=0.059,left=0.073,right=0.987,hspace=0.11,wspace=0.068)

# updating x tick labels
labels = [item.get_text() for item in axes[2,0].get_xticklabels()]
newlabels = []
for lab in labels:
    l = float(lab)/1000
    newlabels.append(int(l))

axes[2,0].set_xticklabels(newlabels)
axes[2,1].set_xticklabels(newlabels)
axes[2,2].set_xticklabels(newlabels)

plt.savefig('/home/danilo/Dropbox/mestrado/figuras/secao3x3.eps')
