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
def create_Structure_3(ncin,indexes):

    fig,axes = plt.subplots(nrows=3,ncols=3,figsize=(25.4/2.54, 20/2.54))

    for ind in indexes:
        if ind == 99:
            axesInd = 0
            axes[axesInd,0].set_title('Climatologia',fontsize=8)
            axes[axesInd,1].set_title(u'EA1 - 15 de Fevereiro',fontsize=8)
            axes[axesInd,2].set_title(u'EA2 - 15 de Fevereiro',fontsize=8)
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

    # adding colorbar axes (cbaxes)
    cax = fig.add_axes([.3,.04,.45,.02])

    return fig,axes,cax

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
    idx_begin = oceano.find_distance_of_a_value(ncin,ind,nstepBegin,sigma,var,value)[1][0][0]
    idx_final = oceano.find_distance_of_a_value(ncin,ind,nstepFinal,sigma,var,value)[1][0][0]

    info = {
        'location': loc,
        'beginPos_X': oceano.find_distance_of_a_value(ncin,ind,nstepBegin,sigma,var,value)[0],
        'finalPos_X': oceano.find_distance_of_a_value(ncin,ind,nstepFinal,sigma,var,value)[0],
        'beginPos_Z': oceano.find_depth_of_a_value(ncin,ind,nstepBegin,idx_begin,sigma,var,value),
        'finalPos_Z': oceano.find_depth_of_a_value(ncin,ind,nstepFinal,idx_final,sigma,var,value)
    }

    return info

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
contours = np.arange(33,36,0.1)


DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
fname = DATA_DIR + exp + '.cdf'

ncin = xr.open_dataset(fname)

lon,lat = ncin['lon'].values, ncin['lat'].values
lon[lon == 0.] = np.nan
lat[lat == 0.] = np.nan
depth = ncin['depth'].values
sigma = ncin['sigma'].values
h1    = ncin['h1'].values
salt  = ncin.salt

fig,axes,caxes = create_Structure_3(ncin,indexes)

# axes[0,0].set_title(u'Climatologia',fontsize=8)
# axes[0,1].set_title(u'Experimento EA1',fontsize=8)
# axes[0,2].set_title(u'Experimento EA2',fontsize=8)

title = u'Seção vertical de salinidade em Ubatuba (superior), Santos (meio) e Cananéia (inferior),\n'\
      + u'com a estrutura climatológica à esquerda e a estrutura no dia 15 de Fevereiro em EA1 (meio) e EA2 (direita).\n'
plt.suptitle(title,fontsize=10)

# defining the begin and the end to plot
nstepBegin = np.arange(0,4,1) # begin of anomalous period
axes[0,0].set_title('Climatologia',fontsize=8)
nstepFinal = np.arange(297,303,1) # final of anomalous period

os.system('clear')
for ind in indexes:
    if ind == 99:
        axesInd = 0
        infos = search_information(ncin,ind,nstepBegin,nstepFinal,'Ubatuba','salt')
    if ind == 28:
        axesInd = 1
        infos = search_information(ncin,ind,nstepBegin,nstepFinal,'Santos','salt')
    if ind == 19:
        axesInd = 2
        infos = search_information(ncin,ind,nstepBegin,nstepFinal,u'Cananéia','salt')

    S = np.nanmean(salt[nstepBegin,:,ind,:],axis=0)
    Splot,ndist,ndepth,dist2,sig,depth = oceano.crossSection_optimized(lon,depth,sigma,h1,S,horizResolution=horizResolution,vertResolution=vertResolution,depRef=depRef,ind=ind)
    # gridding vertical section
    xgrid,zgrid = np.meshgrid(ndist,ndepth)

    # begin: 18 isotherm position
    cf1  = axes[axesInd,0].contourf(xgrid,-zgrid,Splot,contours,cmap=cmo.cm.haline,extend='max')
    cs   = axes[axesInd,0].contour(xgrid,-zgrid,Splot,levels=[36.],colors=('k'),linestyles=('--'))
    axes[axesInd,0].fill_between(dist2[-1,:], -depRef, sig[-1,:],color='#c0c0c0')
    axes[axesInd,0].plot(dist2[-1,:],sig[-1,:],'k')
    axes[axesInd,0].set_xlim([0,limiteEixoX])
    axes[axesInd,0].set_ylim([-depRef,0])

    # plot text box
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    textstr = u'%s' % (infos['location'])
    axes[axesInd,0].text(0.17, 0.32, textstr, transform=axes[axesInd,0].transAxes, fontsize=8,
        va='top', ha='center',bbox=props)

    S = np.nanmean(salt[nstepFinal,:,ind,:],axis=0)
    Splot,ndist,ndepth,dist2,sig,depth = oceano.crossSection_optimized(lon,depth,sigma,h1,S,horizResolution=horizResolution,vertResolution=vertResolution,depRef=depRef,ind=ind)

    # final 18 isotherm position
    cf2  = axes[axesInd,1].contourf(xgrid,-zgrid,Splot,contours,cmap=cmo.cm.haline,extend='max')
    cs   = axes[axesInd,1].contour(xgrid,-zgrid,Splot,levels=[36.],colors=('k'),linestyles=('--'))
    axes[axesInd,1].fill_between(dist2[-1,:], -depRef, sig[-1,:],color='#c0c0c0')
    axes[axesInd,1].plot(dist2[-1,:],sig[-1,:],'k')
    axes[axesInd,1].set_xlim([0,limiteEixoX])
    axes[axesInd,1].set_ylim([-depRef,1])

    # plot text box
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    deltax = infos['finalPos_X']-infos['beginPos_X']
    deltaz = np.abs(infos['finalPos_Z']-infos['beginPos_Z'])
    textstr = r'$\Delta$x = %.1f km'%(deltax)+ '\n'+ r'$\Delta$z = %.1f m' % (deltaz)
    axes[axesInd,1].text(0.17, 0.32, textstr, transform=axes[axesInd,1].transAxes, fontsize=8,
        va='top', ha='center',bbox=props)

    for c in cf1.collections:
        c.set_edgecolor('face')
        c.set_linewidth(0.00000000001)

    for c in cf2.collections:
        c.set_edgecolor('face')
        c.set_linewidth(0.00000000001)

######## Experimento Anomalo 2!!!!
ncin = xr.open_dataset(fname.replace('1','2'))
salt = ncin.salt

for ind in indexes:
    if ind == 99:
        axesInd = 0
        infos = search_information(ncin,ind,nstepBegin,nstepFinal,'Ubatuba','salt')
    if ind == 28:
        axesInd = 1
        infos = search_information(ncin,ind,nstepBegin,nstepFinal,'Santos','salt')
    if ind == 19:
        axesInd = 2
        infos = search_information(ncin,ind,nstepBegin,nstepFinal,u'Cananéia','salt')

    S = np.nanmean(salt[nstepFinal,:,ind,:],axis=0)
    Splot,ndist,ndepth,_,sig,depth = oceano.crossSection_optimized(lon,depth,sigma,h1,S,horizResolution=horizResolution,vertResolution=vertResolution,depRef=depRef,ind=ind)

    # final 18 isotherm position
    cf3  = axes[axesInd,2].contourf(xgrid,-zgrid,Splot,contours,cmap=cmo.cm.haline,extend='max')
    cs   = axes[axesInd,2].contour(xgrid,-zgrid,Splot,levels=[36.],colors=('k'),linestyles=('--'))
    axes[axesInd,2].fill_between(dist2[-1,:], -depRef, sig[-1,:],color='#c0c0c0')
    axes[axesInd,2].plot(dist2[-1,:],sig[-1,:],'k')
    axes[axesInd,2].set_xlim([0,limiteEixoX])
    axes[axesInd,2].set_ylim([-depRef,1])

    # plot text box
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    deltax = infos['finalPos_X']-infos['beginPos_X']
    deltaz = np.abs(infos['finalPos_Z']-infos['beginPos_Z'])
    textstr = r'$\Delta$x = %.1f km'%(deltax)+ '\n'+ r'$\Delta$z = %.1f m' % (deltaz)
    axes[axesInd,2].text(0.17, 0.32, textstr, transform=axes[axesInd,2].transAxes, fontsize=8,
        va='top', ha='center',bbox=props)

    for c in cf3.collections:
        c.set_edgecolor('face')
        c.set_linewidth(0.00000000001)


# plotting colorbar
cbar = plt.colorbar(cf1,orientation='horizontal',cax=caxes,format='%i')
# setting colorbar tick labels
from matplotlib import ticker
tick_locator = ticker.MaxNLocator(nbins=6)
cbar.locator = tick_locator
cbar.update_ticks()

cbar.ax.axes.tick_params(axis='both',which='both',labelsize=8)
cbar.ax.set_title(r'Temperatura ($^o$C)',fontsize=8)

for c in cbar.ax.collections:
    c.set_edgecolor('face')
    c.set_linewidth(0.00000000001)


# ajusta a figura antes de se ajustar os labels do eixo x, pois aumenta-se a quantidade de
# ticks no eixo quando dá tight_layout
plt.tight_layout()
# plt.subplots_adjust(top=0.905,bottom=0.059,left=0.073,right=0.987,hspace=0.11,wspace=0.068)
# plt.subplots_adjust(top=0.925,bottom=0.06,left=0.115,right=0.95,hspace=0.2,wspace=0.28)
plt.subplots_adjust(top=0.900,bottom=0.144,left=0.073,right=0.987,hspace=0.075,wspace=0.058)

# updating x tick labels
labels = [item.get_text() for item in axes[2,0].get_xticklabels()]
newlabels = []
for lab in labels:
    l = float(lab)/1000
    newlabels.append(int(l))

axes[2,0].set_xticklabels(newlabels)
axes[2,1].set_xticklabels(newlabels)
axes[2,2].set_xticklabels(newlabels)

# plt.savefig(BASE_DIR.replace('mestrado/github','Dropbox/mestrado/figuras')+ 'secao_3x3plots.eps')
# plt.savefig('/home/danilo/Dropbox/mestrado/figuras/secoes_verticais/evolucao_salt_3x3.pdf')
