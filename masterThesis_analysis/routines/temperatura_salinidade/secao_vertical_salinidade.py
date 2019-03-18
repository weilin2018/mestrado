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
def create_structure_1row(indexes,exp):

    fig,axes = plt.subplots(nrows=1,ncols=3,figsize=(20.4/2.54, 10/2.54))

    axes[0].set_title(u'Seção Norte',fontsize=8)
    axes[1].set_title(u'Seção Centro',fontsize=8)
    axes[2].set_title(u'Seção Sul',fontsize=8)

    axes[0].set_xlabel(u'Distância [km]', fontsize=8)
    axes[1].set_xlabel(u'Distância [km]', fontsize=8)
    axes[2].set_xlabel(u'Distância [km]', fontsize=8)

    axes[0].set_ylabel(u'Profundidade [m]',fontsize=8)

    # hiding ticks labels
    axes[1].yaxis.set_major_formatter(plt.NullFormatter())
    axes[2].yaxis.set_major_formatter(plt.NullFormatter())

    # adding colorbar axes (cbaxes)
    cax = fig.add_axes([.35,.05,.35,.02])

    return fig,axes,cax

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

# configurações do plot
figsize = (17.4/2.54, 10/2.54)
# index para as latitudes das seções
indexes = [99,28,19]
# configuracoes para as secoes verticais
horizResolution = 10000
vertResolution  = 100 # isso significa que teremos uma resolucao vertical de 1m
depRef          = 200 # profundidade de referencia para interpolacao

limiteEixoX = 300000 # limite, em metros, do eixo X para a secao vertical
contours = np.arange(34.5,36.01,0.05)

fig,axes,cax = create_structure_1row(indexes,exp)
title = u'Seção vertical de salinidade em Ubatuba (esquerda), Santos (meio) e Cananéia (direita),\n'\
      + u'no dia 13 de Fevereiro em %s.\n'% (exp)
plt.suptitle(title,fontsize=10)

nstepBegin = np.arange(0,9,1) # climatolofic data
nstepFinal = np.arange(280,289,1) # final of anomalous period

for ind in indexes:
    if ind == 99:
        axesInd = 0
        infos = ocplot.search_information(ncin,ind,nstepBegin,nstepFinal,'Ubatuba','salt',dz=False)
    if ind == 28:
        axesInd = 1
        infos = ocplot.search_information(ncin,ind,nstepBegin,nstepFinal,'Santos','salt',dz=False)
    if ind == 19:
        axesInd = 2
        infos = ocplot.search_information(ncin,ind,nstepBegin,nstepFinal,u'Cananéia','salt',dz=False)

    print('# ----- PLOTTING [secao: %i] 14 JAN, 2014 ----- #'%(ind))
    T = np.nanmean(ncin.salt[nstepFinal,:,ind,:],axis=0)
    Tplot,ndist,ndepth,dist2,sig,depth = oceano.crossSection_optimized(lon,depth,sigma,h1,T,horizResolution=horizResolution,vertResolution=vertResolution,depRef=depRef,ind=ind)
    # gridding vertical section
    xgrid,zgrid = np.meshgrid(ndist,ndepth)

    # begin: 18 isotherm position
    cf1  = axes[axesInd].contourf(xgrid,-zgrid,Tplot,contours,cmap=cmo.cm.haline,extend='both')
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

# plotting colorbar
cbar = plt.colorbar(cf1,ticks=[34.5,34.9,35.25,35.6,36.],orientation='horizontal',cax=cax,format='%.1f')
# setting colorbar tick labels
# from matplotlib import ticker
# tick_locator = ticker.MaxNLocator(nbins=4)
# cbar.locator = tick_locator
# cbar.update_ticks()

cbar.ax.axes.tick_params(axis='both',which='both',labelsize=8)
cbar.ax.set_title('Salinidade',fontsize=8)

# ajusta a figura antes de se ajustar os labels do eixo x, pois aumenta-se a quantidade de
# ticks no eixo quando dá tight_layout
plt.tight_layout()
# plt.subplots_adjust(top=0.905,bottom=0.059,left=0.073,right=0.987,hspace=0.11,wspace=0.068)
# plt.subplots_adjust(top=0.925,bottom=0.06,left=0.115,right=0.95,hspace=0.2,wspace=0.28)
plt.subplots_adjust(top=0.841,bottom=0.237,left=0.087,right=0.982,hspace=0.099,wspace=0.118)

labels = [item.get_text() for item in axes[0].get_xticklabels()]
newlabels = []
for lab in labels:
    l = float(lab)/1000
    newlabels.append(int(l))

axes[0].set_xticklabels(newlabels)
axes[1].set_xticklabels(newlabels)
axes[2].set_xticklabels(newlabels)

plt.savefig('/home/danilo/Dropbox/mestrado/figuras/secoes_verticais/secao_salinidade_%s.eps'%(exp))
plt.savefig('/home/danilo/Dropbox/mestrado/figuras/lowResolution/secoes_verticais/secao_salinidade_%s.png'%(exp))
plt.close()
%reset -f
