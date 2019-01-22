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
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

from matplotlib import dates
import datetime
import cmocean as cmo

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

plt.ion()
##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def import_data(fname):
    ncin = xr.open_dataset(fname)
    lon = ncin.lon.values
    lat = ncin.lat.values
    lon[lon==0] = np.nan
    lat[lat==0] = np.nan

    dep = ncin.depth.values

    ncin.close()

    return lon,lat,dep

def make_map(ax,llat=-30,ulat=-20,llon=-50,ulon=-39,resolution='l',nmeridians=3,nparallels=2,labels=[True,False,False,True]):

    if resolution == 'f':
        m = pickle.load(open('pickles/basemap.p','r'))
    else:
        m = Basemap(projection='merc', llcrnrlat=llat, urcrnrlat=ulat, llcrnrlon=llon, urcrnrlon=ulon, resolution=resolution)

    m.ax = ax

    m.drawcoastlines(linewidth=.1,color='k')
    m.drawmapboundary(fill_color='#e6ffff')
    m.fillcontinents(color='#c0c0c0')
    m.drawstates(linewidth=.01,color='gray')
	# definir meridianos e paralelos para plotar no mapa
    meridians=np.arange(llon,ulon,nmeridians)
    parallels=np.arange(llat,ulat,nparallels)
	# desenhar meridianos e paralelos conforme definido acima
    m.drawparallels(parallels,labels=labels,fontsize=8,color='gray',linewidth=.2)
    m.drawmeridians(meridians,labels=labels,fontsize=8,color='gray',linewidth=.2)

    return m

def make_minimap(ax,llat=-30,ulat=-20,llon=-50,ulon=-39,resolution='l',nmeridians=3,nparallels=2,labels=[True,False,False,True]):
    m = Basemap(projection='merc', llcrnrlat=llat, urcrnrlat=ulat, llcrnrlon=llon, urcrnrlon=ulon, resolution=resolution)

    m.ax = ax

    m.drawcoastlines(linewidth=.1)
    m.drawmapboundary()
    m.fillcontinents(color='#c0c0c0')

    return m

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/output')

# localizacoes [lat,lon]
# dados de vento in situ
loc_Eslaje      = [-24.32000, -46.180, 'Laje de Santos']
loc_pnboia      = [-25.28000, -44.925, 'PNBOIA']
# dados de vento reanalise/modelo global
loc_ncep_laje   = [-24.22487, -46.228, 'NCEP-LS']
loc_ncep_pnboia = [-25.24701, -45.000, 'NCEP-ST']
# dados de corrente
loc_ES08        = [-25.19,-45.82,'ECOSAN']
loc_css         = [-23.832,-45.405,'TRADICASS']

### FIGURE CONFIGURATION
figsize      = (8.4/2.54, 10./2.54)
# figsize      = (15.4/2.54, 14./2.54)

# import depth
lon,lat,depth   = import_data(DATA_DIR+"EC1.cdf")

# plotando mapa geral
bbox = [1.,-.08] # bbox_to_anchor - legend


fig,ax = plt.subplots(figsize=figsize) #ja em polegadas para artigo
# cax = fig.add_axes(cax_position)

m = make_map(ax,ulat=-23.,llat=-26.2,llon=-47.5,ulon=-44.,resolution='h',nmeridians=1,nparallels=1)
contours = np.arange(5,700.,0.1)
# cf = m.contourf(lon,lat,depth,contours,latlon=True,cmap=cmo.cm.deep,extend='max')
cs = m.contour(lon,lat,depth,levels=[50.,100.,300.],colors=('k'),linestyles=('--'),linewidths=.1,latlon=True)
plt.clabel(cs, fontsize=9, inline=1,fmt='%i')

# insitu meoteorlogico
m.scatter(loc_Eslaje[1],loc_Eslaje[0],s=15,c='g',marker='s',latlon=True,label=loc_Eslaje[2])
m.scatter(loc_pnboia[1],loc_pnboia[0],s=15,c='b',marker='s',latlon=True,label=loc_pnboia[2])

# insitu oceanografico
m.scatter(loc_ES08[1],loc_ES08[0],s=15,c='k',marker='o',latlon=True,label=loc_ES08[2])
m.scatter(loc_css[1],loc_css[0],s=15,c='#e6ffff',marker='o',latlon=True,label=loc_css[2])

m.scatter(loc_ncep_laje[1],loc_ncep_laje[0],s=15,c='r',marker='^',latlon=True,label=loc_ncep_laje[2])
m.scatter(loc_ncep_pnboia[1],loc_ncep_pnboia[0],s=15,c='r',marker='^',latlon=True,label=loc_ncep_pnboia[2])

# criando minimapas
axins = zoomed_inset_axes(ax, 5, loc='upper right')
axins.set_xlim(-47,-44)
axins.set_ylim(-26,-23)

plt.xticks(visible=False)
plt.yticks(visible=False)

map2 = make_minimap(axins,ulat=-23.7,llat=-23.9,llon=-45.51,ulon=-45.3,resolution='f',nmeridians=1,nparallels=1)

map2.contourf(lon,lat,depth,np.arange(5,50,0.1),latlon=True,cmap=cmo.cm.deep)
map2.scatter(loc_css[1],loc_css[0],s=15,c='b',marker='o',latlon=True,label=loc_css[2])

# ax.legend(loc='lower center',fontsize=8,ncol=2)
ax.legend(bbox_to_anchor=bbox,fontsize=7,ncol=3)
# axins.legend(loc='best',fontsize=8)

plt.tight_layout()
plt.subplots_adjust(top=0.977,bottom=0.094,left=0.121,right=0.969,hspace=0.0,wspace=0.2)

"""
Apos salvar em masterThesis_analysis/figuras/, o arquivo .eps salvo foi aberto no
inkscape e então tratado para se obter o resultado final que foi para o relatório
anual de 2019.

arquivo final: masterThesis_analysis/figuras/mapas/pontos_dados_relatorioAnual2019.eps

"""
