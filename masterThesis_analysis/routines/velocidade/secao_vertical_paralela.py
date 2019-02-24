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

    axes[0,0].set_title('Perpendicular',fontsize=8)
    axes[0,1].set_title('Paralela',fontsize=8)

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

    cax = fig.add_axes([.3,.05,.45,.02])

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

def tratando_corrente(u,v,depth,angle):

    ur,vr = rotate_velocityField(u,v,angle)
    spd = np.sqrt(ur**2+vr**2)
    spd = np.where(depth < 200, spd,np.nan)

    return ur,vr,spd

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
figsize = (15.4/2.54, 18/2.54)
# index para as latitudes das seções
indexes = [99,28,19]
# configuracoes para as secoes verticais
horizResolution = 10000
vertResolution  = 100 # isso significa que teremos uma resolucao vertical de 1m
depRef          = 200 # profundidade de referencia para interpolacao

limiteEixoX = 300000 # limite, em metros, do eixo X para a secao vertical
contours = np.arange(-.4,.4,0.01)

fig,axes,cax = create_Structure(figsize)
title = u'Seção vertical de velocidade em Ubatuba (superior), Santos (meio) e \nCananéia (inferior),'\
      + u'no dia 15 de Janeiro em %s.\n'% (exp)
plt.suptitle(title,fontsize=10)

# figure adjustments
plt.tight_layout()
plt.subplots_adjust(top=0.892,bottom=0.158,left=0.113,right=0.979,hspace=0.075,wspace=0.085)

nstepBegin = np.arange(48,57,1)   # first day
nstepFinal = np.arange(280,289,1) # final of anomalous period

# rotating vectors
u = np.nanmean(ncin.u[nstepBegin,:,:,:],axis=0)
v = np.nanmean(ncin.v[nstepBegin,:,:,:],axis=0)
ur,vr = np.zeros(u.shape)*np.nan,np.zeros(u.shape)*np.nan

for t in range(u.shape[0]):
    urot,vrot,spdrot = tratando_corrente(u[t,:,:],v[t,:,:],depth,(-1)*angle)
    ur[t,:,:] = urot
    vr[t,:,:] = vrot


for ind in indexes:
    if ind == 99:
        axesInd = 0
    if ind == 28:
        axesInd = 1
    if ind == 19:
        axesInd = 2

    print('# ----- PLOTTING [secao: %i] 14 JAN, 2014 ----- #'%(ind))
    U = ur[:,ind,:]
    Tplot,ndist,ndepth,dist2,sig,depth = oceano.crossSection_optimized(lon,depth,sigma,h1,U,horizResolution=horizResolution,vertResolution=vertResolution,depRef=depRef,ind=ind)
    # gridding vertical section
    xgrid,zgrid = np.meshgrid(ndist,ndepth)

    # begin: 18 isotherm position
    cf1  = axes[axesInd,0].contourf(xgrid,-zgrid,Tplot,contours,cmap=cmo.cm.delta,extend='both')
    axes[axesInd,0].fill_between(dist2[-1,:], -depRef, sig[-1,:],color='#c0c0c0')
    axes[axesInd,0].plot(dist2[-1,:],sig[-1,:],'k')
    axes[axesInd,0].set_xlim([0,limiteEixoX
    ])
    axes[axesInd,0].set_ylim([-depRef,0])

for ind in indexes:
    if ind == 99:
        axesInd = 0
    if ind == 28:
        axesInd = 1
    if ind == 19:
        axesInd = 2

    print('# ----- PLOTTING [secao: %i] 14 JAN, 2014 ----- #'%(ind))
    U = vr[:,ind,:]
    Tplot,ndist,ndepth,dist2,sig,depth = oceano.crossSection_optimized(lon,depth,sigma,h1,U,horizResolution=horizResolution,vertResolution=vertResolution,depRef=depRef,ind=ind)
    # gridding vertical section
    xgrid,zgrid = np.meshgrid(ndist,ndepth)

    # begin: 18 isotherm position
    cf1  = axes[axesInd,1].contourf(xgrid,-zgrid,Tplot,contours,cmap=cmo.cm.delta,extend='both')
    axes[axesInd,1].fill_between(dist2[-1,:], -depRef, sig[-1,:],color='#c0c0c0')
    axes[axesInd,1].plot(dist2[-1,:],sig[-1,:],'k')
    axes[axesInd,1].set_xlim([0,limiteEixoX
    ])
    axes[axesInd,1].set_ylim([-depRef,0])


# plotting colorbar
cbar = plt.colorbar(cf1,orientation='horizontal',cax=cax,format='%.2f')
# setting colorbar tick labels
from matplotlib import ticker
tick_locator = ticker.MaxNLocator(nbins=6)
cbar.locator = tick_locator
cbar.update_ticks()

cbar.ax.axes.tick_params(axis='both',which='both',labelsize=8)
cbar.ax.set_title(r'Velocidade [m.s$^{-1}$]',fontsize=8)

for c in cbar.ax.collections:
    c.set_edgecolor('face')
    c.set_linewidth(0.00000000001)


# updating x tick labels
labels = [item.get_text() for item in axes[2,0].get_xticklabels()]
newlabels = []
for lab in labels:
    l = float(lab)/1000
    newlabels.append(int(l))

axes[2,0].set_xticklabels(newlabels)
axes[2,1].set_xticklabels(newlabels)

plt.savefig('/home/danilo/Dropbox/mestrado/figuras/secoes_verticais/evolucao_speed_%s.eps'%(exp))
plt.savefig('/home/danilo/Dropbox/mestrado/figuras/lowResolution/secoes_verticais/evolucao_speed_%s.png'%(exp))
plt.close()
%reset -f
