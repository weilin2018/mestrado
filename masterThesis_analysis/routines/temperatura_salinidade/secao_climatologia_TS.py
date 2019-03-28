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
# insert functions here
def create_structure_2x3(indexes,exp):

    fig,axes = plt.subplots(nrows=2,ncols=3,figsize=(20.4/2.54, 15/2.54))

    axes[0,0].set_title(u'Seção Norte',fontsize=8)
    axes[0,1].set_title(u'Seção Centro',fontsize=8)
    axes[0,2].set_title(u'Seção Sul',fontsize=8)

    axes[1,0].set_xlabel(u'Distância [km]', fontsize=8)
    axes[1,1].set_xlabel(u'Distância [km]', fontsize=8)
    axes[1,2].set_xlabel(u'Distância [km]', fontsize=8)

    axes[0,0].set_ylabel('Profundidade [m]',fontsize=8)
    axes[1,0].set_ylabel('Profundidade [m]',fontsize=8)

    # hiding ticks labels
    axes[0,0].xaxis.set_major_formatter(plt.NullFormatter())

    axes[0,1].xaxis.set_major_formatter(plt.NullFormatter())
    axes[0,1].yaxis.set_major_formatter(plt.NullFormatter())
    axes[0,2].xaxis.set_major_formatter(plt.NullFormatter())
    axes[0,2].yaxis.set_major_formatter(plt.NullFormatter())

    axes[1,1].yaxis.set_major_formatter(plt.NullFormatter())
    axes[1,2].yaxis.set_major_formatter(plt.NullFormatter())

    # adding colorbar axes (cbaxes)
    caxTemp = fig.add_axes([.115,.04,.405,.02])
    caxSalt = fig.add_axes([.565,.04,.405,.02])

    return fig,axes,caxTemp,caxSalt


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
plt.ion()

exp = 'EA1'

# configurações do plot
figsize = (17.4/2.54, 10/2.54)
# index para as latitudes das seções
indexes = [99,28,19]
# configuracoes para as secoes verticais
horizResolution = 10000
vertResolution  = 100 # isso significa que teremos uma resolucao vertical de 1m
depRef          = 200 # profundidade de referencia para interpolacao

limiteEixoX = 300000 # limite, em metros, do eixo X para a secao vertical

# criando contours
contours_temp = np.arange(11,35,0.5)
contours_salt = np.arange(34.5,36.01,0.05)
contours_velo = np.arange(0,1.5,0.01)

DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
fname = DATA_DIR + exp + '.cdf'

ncin = xr.open_dataset(fname)

lon,lat = ncin['lon'].values, ncin['lat'].values
lon[lon == 0.] = np.nan
lat[lat == 0.] = np.nan
depth = ncin['depth'].values
sigma = ncin['sigma'].values
h1    = ncin['h1'].values

nstepBegin = np.arange(0,9,1) # climatolofic data

# fig,axes,caxTemp,caxSalt,caxVelo = create_Structure_3(ncin,indexes)
fig,axes,caxTemp,caxSalt = create_structure_2x3(ncin,indexes)
# title = u'Seção vertical de temperatura (superior) e salinidade (inferior) \nem Ubatuba (esquerda), Santos (meio) e Cananéia (direita)'
# plt.suptitle(title,fontsize=10)

# plotando primeiro a temperatura
os.system('clear')
for ind in indexes:
    if ind == 99:
        axesInd = 0
    if ind == 28:
        axesInd = 1
    if ind == 19:
        axesInd = 2
    print('Plotting temperature in [%i]' % (ind))
    data = np.nanmean(ncin.temp[nstepBegin,:,ind,:],axis=0)
    Dplot,ndist,ndepth,dist2,sig,depth = oceano.crossSection_optimized(lon,depth,sigma,h1,data,horizResolution=horizResolution,vertResolution=vertResolution,depRef=depRef,ind=ind)
    # gridding vertical section
    xgrid,zgrid = np.meshgrid(ndist,ndepth)

    # begin: 18 isotherm position
    cf1  = axes[0,axesInd].contourf(xgrid,-zgrid,Dplot,contours_temp,cmap=cmo.cm.thermal,extend='max')
    cs   = axes[0,axesInd].contour(xgrid,-zgrid,Dplot,levels=[18.],colors=('k'),linestyles=('--'))
    axes[0,axesInd].fill_between(dist2[-1,:], -depRef, sig[-1,:],color='#c0c0c0')
    axes[0,axesInd].plot(dist2[-1,:],sig[-1,:],'k')
    axes[0,axesInd].set_xlim([0,limiteEixoX])
    axes[0,axesInd].set_ylim([-depRef,0])

    # for c in axes[0,axesInd].collections:
    #     c.set_edgecolor('face')
    #     c.set_linewidth(0.00000000001)

# plotting colorbar
cbar = plt.colorbar(cf1,orientation='horizontal',cax=caxTemp,format='%i')
# setting colorbar tick labels
from matplotlib import ticker
tick_locator = ticker.MaxNLocator(nbins=6)
cbar.locator = tick_locator
cbar.update_ticks()

cbar.ax.axes.tick_params(axis='both',which='both',labelsize=8)
cbar.ax.set_title(r'Temperatura ($^o$C)',fontsize=8)


# plotando primeiro a salinidade
for ind in indexes:
    if ind == 99:
        axesInd = 0
    if ind == 28:
        axesInd = 1
    if ind == 19:
        axesInd = 2
    print('Plotting salinity in [%i]' % (ind))
    data = np.nanmean(ncin.salt[nstepBegin,:,ind,:],axis=0)
    Dplot,ndist,ndepth,dist2,sig,depth = oceano.crossSection_optimized(lon,depth,sigma,h1,data,horizResolution=horizResolution,vertResolution=vertResolution,depRef=depRef,ind=ind)
    # gridding vertical section
    xgrid,zgrid = np.meshgrid(ndist,ndepth)

    # begin: 18 isotherm position
    cf1  = axes[1,axesInd].contourf(xgrid,-zgrid,Dplot,contours_salt,cmap=cmo.cm.haline,extend='both')
    cs   = axes[1,axesInd].contour(xgrid,-zgrid,Dplot,levels=[36.],colors=('k'),linestyles=('--'))
    axes[1,axesInd].fill_between(dist2[-1,:], -depRef, sig[-1,:],color='#c0c0c0')
    axes[1,axesInd].plot(dist2[-1,:],sig[-1,:],'k')
    axes[1,axesInd].set_xlim([0,limiteEixoX])
    axes[1,axesInd].set_ylim([-depRef,0])

    # for c in axes[1,axesInd].collections:
    #     c.set_edgecolor('face')
    #     c.set_linewidth(0.00000000001)

# plotting colorbar
cbar = plt.colorbar(cf1,ticks=[34.5,34.9,35.25,35.6,36.],orientation='horizontal',cax=caxSalt,format='%.1f')
# setting colorbar tick labels
from matplotlib import ticker
tick_locator = ticker.MaxNLocator(nbins=6)
cbar.locator = tick_locator
cbar.update_ticks()

cbar.ax.axes.tick_params(axis='both',which='both',labelsize=8)
cbar.ax.set_title('Salinidade',fontsize=8)

plt.tight_layout()
plt.subplots_adjust(top=0.946,bottom=0.165,left=0.085,right=0.985,hspace=0.072,wspace=0.072)

# updating x tick labels
labels = [item.get_text() for item in axes[1,2].get_xticklabels()]
newlabels = []
for lab in labels:
    l = float(lab)/1000
    newlabels.append(int(l))

axes[1,0].set_xticklabels(newlabels)
axes[1,1].set_xticklabels(newlabels)
axes[1,2].set_xticklabels(newlabels)

plt.savefig('/home/danilo/Dropbox/mestrado/figuras/secoes_verticais/secao_climatologia_TS.eps')
plt.savefig('/home/danilo/Dropbox/mestrado/figuras/lowResolution/secoes_verticais/secao_climatologia_TS.png',dpi=150)
plt.close()
%reset -f
