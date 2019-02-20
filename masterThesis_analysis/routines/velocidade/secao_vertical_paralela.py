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

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano
import masterThesisPack.plots as ocplot

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def create_Structure(figsize):

    fig,axes = plt.subplots(nrows=3,ncols=2,figsize=figsize)

    axes[0,0].set_title('Paralela',fontsize=8)
    axes[0,1].set_title('Perpendicular',fontsize=8)

    axes[0,0].xaxis.set_major_formatter(plt.NullFormatter())
    axes[0,1].xaxis.set_major_formatter(plt.NullFormatter())
    axes[1,0].xaxis.set_major_formatter(plt.NullFormatter())
    axes[1,1].xaxis.set_major_formatter(plt.NullFormatter())
    axes[0,1].yaxis.set_major_formatter(plt.NullFormatter())
    axes[1,1].yaxis.set_major_formatter(plt.NullFormatter())
    axes[2,1].yaxis.set_major_formatter(plt.NullFormatter())

    axes[0,0].set_ylabel(u'Profundidade [m]',fontsize=8)
    axes[1,0].set_ylabel(u'Profundidade [m]',fontsize=8)
    axes[2,0].set_ylabel(u'Profundidade [m]',fontsize=8)
    axes[2,0].set_xlabel(u'Distância [km]', fontsize=8)
    axes[2,1].set_xlabel(u'Distância [km]', fontsize=8)

    cax = fig.add_axes([.35,.05,.35,.02])

    return fig,axes,cax

def rotate_velocityField(u,v,ang):

    import decomp
    ur = np.zeros(u.shape)*np.nan
    vr = np.zeros(v.shape)*np.nan

    for j in range(u.shape[0]):
        U,V = u[j,:],v[j,:]
        angle = ang[j,:]

        INT,DIR = decomp.uv2intdir(U,V,0,angle)
        uro,vro = decomp.intdir2uv(INT,DIR,0,angle)
        ur[j,:] = uro
        vr[j,:] = vro

    return ur,vr

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
plt.ion()

exp = raw_input('Digite o experimento para plotar secao vertical: ')

DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
fname = DATA_DIR + exp + '.cdf'

ncin = xr.open_dataset(fname)

lon,lat = ncin['lon'].values, ncin['lat'].values
lon[lon == 0.] = np.nan
lat[lat == 0.] = np.nan
depth = ncin['depth'].values
sigma = ncin['sigma'].values
h1    = ncin['h1'].values
angle = ncin['ang'].values
# configurações do plot
figsize = (15.4/2.54, 15/2.54)
# index para as latitudes das seções
indexes = [99,28,19]
# configuracoes para as secoes verticais
horizResolution = 10000
vertResolution  = 100 # isso significa que teremos uma resolucao vertical de 1m
depRef          = 200 # profundidade de referencia para interpolacao

limiteEixoX = 300000 # limite, em metros, do eixo X para a secao vertical
contours = np.arange(34.1,37.1,0.1)

fig,axes,cax = create_Structure(figsize)
title = u'Seção vertical de salinidade em Ubatuba (esquerda), Santos (meio) e Cananéia (direita),\n'\
      + u'no dia 13 de Fevereiro em %s.\n'% (exp)
plt.suptitle(title,fontsize=10)

nstepBegin = np.arange(0,9,1) # climatolofic data
nstepFinal = np.arange(280,289,1) # final of anomalous period

# rotating vectors
u = np.nanmean(ncin.u[nstepFinal,:,indexes,:],axis=0)
v = np.nanmean(ncin.v[nstepFinal,:,indexes,:],axis=0)
ur,vr = rotate_velocityField(u,v,angle[indexes,:])

for t in range(umean.shape[0]):
    urot,vrot,spdrot = tratando_corrente(umean[t,:,:],vmean[t,:,:],depth,(-1)*angle)
    ur[t,:,:] = urot
    vr[t,:,:] = vrot
    spd[t,:,:]= spdrot

for ind in indexes:
    if ind == 99:
        axesInd = 0
    if ind == 28:
        axesInd = 1
    if ind == 19:
        axesInd = 2

    print('# ----- PLOTTING [secao: %i] 14 JAN, 2014 ----- #'%(ind))
    U = np.nanmean(ncin.u[nstepFinal,:,ind,:],axis=0)
    Tplot,ndist,ndepth,dist2,sig,depth = oceano.crossSection_optimized(lon,depth,sigma,h1,U,horizResolution=horizResolution,vertResolution=vertResolution,depRef=depRef,ind=ind)
    # gridding vertical section
    xgrid,zgrid = np.meshgrid(ndist,ndepth)

    # begin: 18 isotherm position
    cf1  = axes[axesInd].contourf(xgrid,-zgrid,Tplot,contours,cmap=cmo.cm.haline,extend='max')
    cs   = axes[axesInd].contour(xgrid,-zgrid,Tplot,levels=[36.],colors=('k'),linestyles=('--'))
    axes[axesInd].fill_between(dist2[-1,:], -depRef, sig[-1,:],color='#c0c0c0')
    axes[axesInd].plot(dist2[-1,:],sig[-1,:],'k')
    axes[axesInd].set_xlim([0,limiteEixoX
    ])
    axes[axesInd].set_ylim([-depRef,0])
    # plot text box
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    deltax = infos['finalPos_X']-infos['beginPos_X']
    textstr = r'$\Delta$x = %.1f km'%(deltax)
    axes[axesInd].text(0.18, 0.22, textstr, transform=axes[axesInd].transAxes, fontsize=8,
        va='top', ha='center',bbox=props)
