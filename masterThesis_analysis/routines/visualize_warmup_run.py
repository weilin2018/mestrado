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

def create_figuresThesis(ncin,timestep,var='temp',title=u'Condição Inicial de Temperatura'):

    """
        >>> contours = np.arange(33,38,.2)
        >>> create_figuresThesis(ncin,final,var='salt',title=u'Condição Inicial de Salinidade')

    """

    # create dictionary with colorbars for each variable
    contours = {
        'temp': {
            'surf' : np.arange(21,28,.1),
            'mean' : np.arange(2,28,.1),
            'bott' : np.arange(2,28,.1),
            'colorbar': cmo.cm.thermal,
            'tickInterval': [2,7,7],
            'cb_title': r'Temperatura ($^o$C)',
        },
        'salt': {
            'surf' : np.arange(33,37,.01),
            'mean' : np.arange(33,37,.01),
            'bott' : np.arange(33,37,.01),
            'colorbar': cmo.cm.haline,
            'tickInterval': [1,1,1],
            'cb_title': 'Salinidade',
        }
    }

    lon = ncin.lon.values
    lon[lon == 0.] = np.nan
    lat = ncin.lat.values
    lat[lat == 0.] = np.nan

    zpos = ncin.dims['sigma'] # qtd de niveis sigma

    varSurface = ncin[var][timestep,0,:,:].values
    varInterme = ncin[var][timestep,zpos/2,:,:].values
    varBottom  = ncin[var][timestep,-1,:,:].values

    fig,ax = plt.subplots(ncols=3,figsize=(20,6))

    # plotting
    caxSurface = fig.add_axes([0.14, 0.82, 0.05, 0.015])
    mSurface = oceano.make_map(ax[0],resolution='f')
    cb_ticks = np.arange(contours[var]['surf'].min(),round(contours[var]['surf'].max()),contours[var]['tickInterval'][0])
    csurface = mSurface.contourf(lon,lat,varSurface,contours[var]['surf'],latlon=True,cmap=contours[var]['colorbar'])
    cbarSurf = plt.colorbar(csurface,cax=caxSurface,orientation='horizontal',ticks=cb_ticks)
    plt.title(contours[var]['cb_title'],fontsize=10)
    ax[0].set_title(u'Superfície',fontsize=15)

    caxInterme = fig.add_axes([0.413,0.82,0.05,0.015])
    mInterme = oceano.make_map(ax[1],resolution='f')
    cb_ticks = np.arange(contours[var]['surf'].min(),round(contours[var]['surf'].max()),contours[var]['tickInterval'][1])
    cinterme = mInterme.contourf(lon,lat,varInterme,contours[var]['surf'],latlon=True,cmap=contours[var]['colorbar'])
    cbarInte = plt.colorbar(cinterme,cax=caxInterme,orientation='horizontal',ticks=cb_ticks)
    plt.title(contours[var]['cb_title'],fontsize=10)
    ax[1].set_title(u'Meio',fontsize=15)

    caxBottom = fig.add_axes([0.685,0.82,0.05,0.015])
    mBottom = oceano.make_map(ax[2],resolution='f')
    cb_ticks = np.arange(contours[var]['bott'].min(),round(contours[var]['bott'].max()),contours[var]['tickInterval'][2])
    cbottom = mBottom.contourf(lon,lat,varBottom,contours[var]['bott'],latlon=True,cmap=contours[var]['colorbar'])
    cbarBott = plt.colorbar(cbottom,cax=caxBottom,orientation='horizontal',ticks=cb_ticks)
    plt.title(contours[var]['cb_title'],fontsize=10)
    ax[2].set_title('Fundo',fontsize=15)

    plt.suptitle(title,fontsize=24)


def create(ncin,timestep,title=u'Condição Inicial de Temperatura, Salinidade e Corrente na Superfície'):

    """
        >>> contours = np.arange(33,38,.2)
        >>> create_figuresThesis(ncin,final,var='salt',title=u'Condição Inicial de Salinidade')

    """

    # create dictionary with colorbars for each variable
    contours = {
        'temp': {
            'surf' : np.arange(21,28,.1),
            'mean' : np.arange(2,28,.1),
            'bott' : np.arange(2,28,.1),
            'colorbar': cmo.cm.thermal,
            'tickInterval': [2,7,7],
            'cb_title': r'Temperatura ($^o$C)',
        },
        'salt': {
            'surf' : np.arange(33,37,.01),
            'mean' : np.arange(33,37,.01),
            'bott' : np.arange(33,37,.01),
            'colorbar': cmo.cm.haline,
            'tickInterval': [1,1,1],
            'cb_title': 'Salinidade',
        },
        'current': {
            'surf' : np.arange(0,0.16,.01),
            'mean' : np.arange(33,37,.01),
            'bott' : np.arange(33,37,.01),
            'colorbar': cmo.cm.speed,
            'tickInterval': [0.05,0.05,0.05],
            'cb_title': r'Corrente (m s$^{-1}$)',
        }
    }

    lon = ncin.lon.values
    lon[lon == 0.] = np.nan
    lat = ncin.lat.values
    lat[lat == 0.] = np.nan

    zpos = ncin.dims['sigma'] # qtd de niveis sigma

    varSST = ncin['temp'][timestep,0,:,:].values
    varSSS = ncin['salt'][timestep,0,:,:].values
    u,v    = ncin['u'][timestep,0,:,:].values,ncin['v'][timestep,0,:,:].values
    varSSC      = np.sqrt(u**2 + v**2)

    # normalizing u,v vectors
    un,vn = u/varSSC, v/varSSC

    if timestep == 0:

        fig,ax = plt.subplots(ncols=2,figsize=(20,6))

        # plotting
        caxTemperature = fig.add_axes([0.14, 0.78, 0.05, 0.015])
        mTemp = oceano.make_map(ax[0],resolution='c')
        cb_ticks = np.arange(contours['temp']['surf'].min(),round(contours['temp']['surf'].max()),contours['temp']['tickInterval'][0])
        ctemperature = mTemp.contourf(lon,lat,varSST,contours['temp']['surf'],latlon=True,cmap=contours['temp']['colorbar'])
        cbarSurf = plt.colorbar(ctemperature,cax=caxTemperature,orientation='horizontal',ticks=cb_ticks)
        plt.title(contours['temp']['cb_title'],fontsize=10)
        ax[0].set_title(u'Temperatura',fontsize=15)

        caxSalt = fig.add_axes([0.413,0.78,0.05,0.015])
        mSalt = oceano.make_map(ax[1],resolution='c')
        cb_ticks = np.arange(contours['salt']['mean'].min(),round(contours['salt']['mean'].max()),contours['salt']['tickInterval'][1])
        csalt = mSalt.contourf(lon,lat,varSSS,contours['salt']['mean'],latlon=True,cmap=contours['salt']['colorbar'])
        cbarInte = plt.colorbar(csalt,cax=caxSalt,orientation='horizontal',ticks=cb_ticks)
        plt.title(contours['salt']['cb_title'],fontsize=10)
        ax[1].set_title(u'Salinidade',fontsize=15)

        plt.subplots_adjust(top=0.86,bottom=0.045,left=0.125,right=0.9,hspace=0.2,wspace=0.2)

    else:

        fig,ax = plt.subplots(ncols=3,figsize=(20,6))

        # plotting
        caxTemperature = fig.add_axes([0.14, 0.78, 0.05, 0.015])
        mTemp = oceano.make_map(ax[0],resolution='c')
        cb_ticks = np.arange(contours['temp']['surf'].min(),round(contours['temp']['surf'].max()),contours['temp']['tickInterval'][0])
        ctemperature = mTemp.contourf(lon,lat,varSST,contours['temp']['surf'],latlon=True,cmap=contours['temp']['colorbar'])
        cbarSurf = plt.colorbar(ctemperature,cax=caxTemperature,orientation='horizontal',ticks=cb_ticks)
        plt.title(contours['temp']['cb_title'],fontsize=10)
        ax[0].set_title(u'Temperatura',fontsize=15)

        caxSalt = fig.add_axes([0.413,0.78,0.05,0.015])
        mSalt = oceano.make_map(ax[1],resolution='c')
        cb_ticks = np.arange(contours['salt']['mean'].min(),round(contours['salt']['mean'].max()),contours['salt']['tickInterval'][1])
        csalt = mSalt.contourf(lon,lat,varSSS,contours['salt']['mean'],latlon=True,cmap=contours['salt']['colorbar'])
        cbarInte = plt.colorbar(csalt,cax=caxSalt,orientation='horizontal',ticks=cb_ticks)
        plt.title(contours['salt']['cb_title'],fontsize=10)
        ax[1].set_title(u'Salinidade',fontsize=15)

        var = 'current'
        caxCurrent = fig.add_axes([0.685,0.78,0.05,0.015])
        mCurr = oceano.make_map(ax[2],resolution='c')
        # cb_ticks = np.arange(contours[var]['surf'].min(),round(contours[var]['surf'].max()),contours[var]['tickInterval'][2])
        cc = mCurr.contour(lon,lat,ncin.depth.values,levels=[100,200,300],colors=('k'),latlon=True,ls='--',alpha=.4)
        cb_ticks = [0.0,0.07,0.14]
        ccurrent = mCurr.contourf(lon,lat,varSSC,contours[var]['surf'],latlon=True,cmap=contours[var]['colorbar'])
        qv = mCurr.quiver(lon[::2,::5],lat[::2,::5],un[::2,::5],vn[::2,::5],latlon=True,scale=60,minshaft=2)
        cbarBott = plt.colorbar(ccurrent,cax=caxCurrent,orientation='horizontal',ticks=cb_ticks)
        plt.title(contours[var]['cb_title'],fontsize=10)
        ax[2].set_title('Corrente',fontsize=15)

        plt.subplots_adjust(top=0.86,bottom=0.045,left=0.125,right=0.9,hspace=0.2,wspace=0.2)


    plt.suptitle(title,fontsize=24)

warm = 'warmupControle.cdf'
fname = '/media/danilo/Danilo/mestrado/ventopcse/output/%s'%(warm)

# instanciating Experiment
ncin = xr.open_dataset(fname)

plt.ion()
# # plotar a climatologia elaborada e lida pelo modelo no instante inicial:
# create_figuresThesis(ncin,0,var='temp',title='Climatologia de Temperatura')
# create_figuresThesis(ncin,0,var='salt',title='Climatologia de Salinidade')
#
# # plotar condição inicial com campos baroclínicos estabilizados
# create_figuresThesis(ncin,-1,var='temp',title=u'Condição Inicial de Temperatura')
# create_figuresThesis(ncin,-1,var='salt',title=u'Condição Inicial de Salinidade')

# def create_figuresThesis_velocidade(ncin):

lon = ncin.lon.values
lon[lon == 0.] = np.nan
lat = ncin.lat.values
lat[lat == 0.] = np.nan

# SSV - sea surface velocity
u = ncin.u[-1,0,:,:].values
v = ncin.v[-1,0,:,:].values
s = np.sqrt(u**2 + v**2)

un = u/s
vn = v/s

fig,ax = plt.subplots()
m = oceano.make_map(ax)
cf = m.contourf(lon,lat,s,latlon=True,cmap=cmo.cm.speed)
cc = m.contour(lon,lat,ncin.depth.values,levels=[100,200,300],colors=('k'),latlon=True,ls='--',alpha=.4)

#qv = m.quiver(lon[::2,::5],lat[::2,::5],un[::2,::5],vn[::2,::5],latlon=True,scale=60,minshaft=2)

lonext = lon.copy()
latext = lat.copy()
uext = un.copy()
vext = vn.copy()
lonext[40:120,:] = np.nan
latext[40:120,:] = np.nan
uext[40:120,:] = np.nan
vext[40:120,:] = np.nan

qv = m.quiver(lonext[::2,::5],latext[::2,::5],uext[::2,::5],vext[::2,::5],latlon=True,scale=60,minshaft=2,alpha=.4)


# ajustar somente para o espaco do canal. lon[40:120,:],lat[40:120,:],ucsb[40:120,:],vcsb[40:120,:]
loncsb = lon[40:120,:]
latcsb = lat[40:120,:]
ucsb   = un[40:120,:]
vcsb   = vn[40:120,:]

qv = m.quiver(loncsb[::5,::10],latcsb[::5,::10],ucsb[::5,::10],vcsb[::5,::10],latlon=True,scale=60,minshaft=2,alpha=.4)

# plot some vectors with full colors and larger
