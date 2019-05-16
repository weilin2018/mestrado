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

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano
import masterThesisPack.plots as ocplot
##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def create_structure():

    fig,axes = plt.subplots(nrows=2,ncols=3,figsize=(20.4/2.54, 13/2.54))
    cax = fig.add_axes([0.3,0.05,0.51,0.02])

    # setting titles
    axes[0,0].set_title('Norte/Ubatuba',fontsize=8)
    axes[0,1].set_title('Centro/Santos',fontsize=8)
    axes[0,2].set_title(u'Sul/Cananéia',fontsize=8)

    for i in range(axes.shape[0]):
        axes[i,0].set_ylabel('Profundidade [m]',fontsize=8)

    for j in range(axes.shape[1]):
        axes[1,j].set_xlabel(u'Distância [km]',fontsize=8)

    # hiding ticklabels
    for j in np.arange(0,axes.shape[1]):
        if j > 0:
            # hide y axis
            axes[0,j].yaxis.set_major_formatter(plt.NullFormatter())
            axes[1,j].yaxis.set_major_formatter(plt.NullFormatter())

        # hide x axis
        axes[0,j].xaxis.set_major_formatter(plt.NullFormatter())

    plt.tight_layout()
    plt.subplots_adjust(top=0.96,bottom=0.191,left=0.117,right=0.969,hspace=0.049,wspace=0.11)

    return fig,axes,cax

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
plt.ion()


os.system('clear')
exp = 'EA1' #raw_input('Digite o experimento a ser plotado: ')
save= int(raw_input('Digite 0 para salvar em PNG e 1 para salvar em PDF: '))

if save == 0:
    extension = '.png'
else:
    extension = '.pdf'

# configurações do plot
# index para as latitudes das seções
indexes = [99,28,19]
# configuracoes para as secoes verticais
horizResolution = 10000
vertResolution  = 100 # isso significa que teremos uma resolucao vertical de 1m
depRef          = 200 # profundidade de referencia para interpolacao

limiteEixoX = 300000 # limite, em metros, do eixo X para a secao vertical

DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
fname = DATA_DIR + exp + '.cdf'

ncin = xr.open_dataset(fname)

lon,lat = ncin['lon'].values, ncin['lat'].values
lon[lon == 0.] = np.nan
lat[lat == 0.] = np.nan
depth = ncin['depth'].values
sigma = ncin['sigma'].values
h1    = ncin['h1'].values

# property's configuration
var = 'salt'                          # property to be ploted
cmap= cmo.cm.haline                   # type of colorbar to be used in contourf
contours = np.arange(34.5,36.01,0.05) # contours to be used as resolution of the colorbar
levels = [36.]                        # to be used in contour (plot isoline)
outName= var + 'Section' + extension

fig,axes,cax = create_structure()
# title = u'Seção vertical de salinidade e, Ubatuba (superior), Santos (meio) e Cananéia (inferior),\n'\
#       + u'com a estrutura climatológica à esquerda e a estrutura no dia 13 de Fevereiro em %s\n'% (exp)\
#       + u'à direita.'
# plt.suptitle(title,fontsize=10)

# defining the begin and the end to plot
nstepBegin = np.arange(0,9,1) # climatolofic data
nstepFinal = np.arange(280,289,1) # final of anomalous period

for ind in indexes:
    if ind == 99:
        axesInd = 0
        infos = ocplot.search_information(ncin,ind,nstepBegin,nstepFinal,'Ubatuba',var,dz=True)
    if ind == 28:
        axesInd = 1
        infos = ocplot.search_information(ncin,ind,nstepBegin,nstepFinal,'Santos',var,dz=True)
    if ind == 19:
        axesInd = 2
        infos = ocplot.search_information(ncin,ind,nstepBegin,nstepFinal,u'Cananéia',var,dz=True)


    print('# ----- PLOTTING [secao: %i] CLIMATOLOGY ----- #'%(ind))
    V = np.nanmean(ncin[var][nstepBegin,:,ind,:].values,axis=0)
    Vplot,ndist,ndepth,dist2,sig,depth = oceano.crossSection_optimized(lon,depth,sigma,h1,V,horizResolution=horizResolution,vertResolution=vertResolution,depRef=depRef,ind=ind)
    # gridding vertical section
    xgrid,zgrid = np.meshgrid(ndist,ndepth)

    # begin: 18 isotherm position
    cf1  = axes[0,axesInd].contourf(xgrid,-zgrid,Vplot,contours,cmap=cmap,extend='max')
    axes[0,axesInd].fill_between(dist2[-1,:], -depRef, sig[-1,:],color='#c0c0c0')
    axes[0,axesInd].plot(dist2[-1,:],sig[-1,:],'k')
    axes[0,axesInd].set_xlim([0,limiteEixoX])
    axes[0,axesInd].set_ylim([-depRef,0])

    if save:
        for c in cf1.collections:
            c.set_edgecolor('face')
            c.set_linewidth(0.00000000001)

    cs   = axes[0,axesInd].contour(xgrid,-zgrid,Vplot,levels=levels,colors=('k'),linestyles=('--'))

    print('# ----- PLOTTING [secao: %i] ANOMALOUS ----- #'%(ind))
    V = np.nanmean(ncin[var][nstepFinal,:,ind,:].values,axis=0)
    Vplot,ndist,ndepth,dist2,sig,depth = oceano.crossSection_optimized(lon,depth,sigma,h1,V,horizResolution=horizResolution,vertResolution=vertResolution,depRef=depRef,ind=ind)
    # gridding vertical section
    xgrid,zgrid = np.meshgrid(ndist,ndepth)

    # begin: 18 isotherm position
    cf1  = axes[1,axesInd].contourf(xgrid,-zgrid,Vplot,contours,cmap=cmap,extend='max')
    axes[1,axesInd].fill_between(dist2[-1,:], -depRef, sig[-1,:],color='#c0c0c0')
    axes[1,axesInd].plot(dist2[-1,:],sig[-1,:],'k')
    axes[1,axesInd].set_xlim([0,limiteEixoX])
    axes[1,axesInd].set_ylim([-depRef,0])

    if save:
        for c in cf1.collections:
            c.set_edgecolor('face')
            c.set_linewidth(0.00000000001)

    cs   = axes[1,axesInd].contour(xgrid,-zgrid,Vplot,levels=levels,colors=('k'),linestyles=('--'))

# plotting colorbar
cbar = plt.colorbar(cf1,ticks=[34.5,34.9,35.25,35.6,36.],orientation='horizontal',cax=cax,format='%.1f')
# setting colorbar tick labels
# from matplotlib import ticker
# tick_locator = ticker.MaxNLocator(nbins=4)
# cbar.locator = tick_locator
# cbar.update_ticks()

cbar.ax.axes.tick_params(axis='both',which='both',labelsize=8)
cbar.ax.set_title('Salinidade',fontsize=8)

# updating x tick labels
labels = [item.get_text() for item in axes[1,0].get_xticklabels()]
newlabels = []
for lab in labels:
    l = float(lab)/1000
    newlabels.append(int(l))

axes[1,0].set_xticklabels(newlabels)
axes[1,1].set_xticklabels(newlabels)
axes[1,2].set_xticklabels(newlabels)

if save == 0:
    plt.savefig('/media/danilo/Danilo/mestrado/gitlab/mestrado/dissertacao/presentation/figures/results/'+outName,dpi=300)
else:
    plt.savefig('/media/danilo/Danilo/mestrado/gitlab/mestrado/dissertacao/presentation/figures/results/'+outName)
