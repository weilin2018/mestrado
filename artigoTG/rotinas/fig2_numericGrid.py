#!/usr/bin/env python2
#-*-coding:utf-8-*-

import numpy as np # atribuir as arrays do .cdf para numpy's arrays
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset, inset_axes
import xray as xr

import cmocean as cmo

import pickle

from artigoTGpack import artigoTGpack as oceano


# # diretorio base
BASE_DIR = oceano.make_dir()

SIMS_DIR = BASE_DIR.replace('github/', 'artigo_data/simulacoes/')
SIMS_DIR = '/media/danilo/Danilo/mestrado/artigo_data/simulacoes/2ndPhase/'

# importar batimetria direto de uma saida do modelo
#data = xr.open_dataset('/media/danilo/Danilo/mestrado/artigo_data/simulacoes/val3/gcmplt.cdf')
data = xr.open_dataset(SIMS_DIR+'reRUN00.cdf')
ncin = xr.open_dataset(SIMS_DIR+'depth_rerun00.cdf')
depth = ncin['depth'].values
lon = data['lon'].values
lat = data['lat'].values
lon[lon == 0.] = np.nan
lat[lon == 0.] = np.nan


# definindo coordenadas para plotar a região desejada
llon = -44.9;
ulon = -43.45;

llat = -23.6;
ulat = -22.75;

piraquaraf    = [-44.440706, -23.010962]
picinguaba    = [-44.830344, -23.378115]
guaratiba     = [-43.561320, -23.064399]
MambucabaArrow= [-44.525,-23.03, -2000,2000,-5000,5000]

## scale position
scaleLon      = [-44.65]
scaleLat      = [-23.55]

############### FIGURE PARAMETERS
figWdith  = 17.4/2.54
figHeight = 12./2.54

fig = plt.figure(figsize=(figWdith,figHeight))
ax = fig.add_subplot(111)

# -----------------------------------------------------------------------------
# ---------------------    PLOTANDO A BASE DO MAPA   --------------------------
# -----------------------------------------------------------------------------
m = Basemap(projection='merc', llcrnrlat=llat, urcrnrlat=ulat, llcrnrlon=llon, urcrnrlon=ulon, resolution='f')
m.ax = ax
# plotar outras coisas do mapa
m.drawcoastlines(linewidth=0.2) #linha de costa em alta resolução
m.drawmapboundary(fill_color='#f2f2f2') # fill_color colore o oceano
# m.fillcontinents(color='#ffd480') # colorir o continente
m.fillcontinents(color='0.85')
m.drawrivers()
# definir meridianos e paralelos para plotar no mapa
meridians=np.arange(-44.9,-33.4,0.25)
parallels=np.arange(-23.6,-22.75,0.1)
# desenhar meridianos e paralelos conforme definido acima
m.drawparallels(parallels,labels=[True,False,False,True],fontsize=13,
                                    fontweight='bold',color='gray',linewidth=0.5)
m.drawmeridians(meridians,labels=[True,False,False,True],fontsize=13,
                                    fontweight='bold',color='gray',linewidth=0.5)

m.drawmapscale(scaleLon[0], scaleLat[0], scaleLon[0], scaleLat[0], 20, barstyle='fancy', yoffset=1000)

m.plot(lon,lat,'k',alpha=.3,latlon=True)
m.plot(lon.T,lat.T,'k',alpha=.3,latlon=True)

x,y = m(piraquaraf[0],piraquaraf[1])
plt.scatter(x,y,marker='*',color='black',s=80,label='Piraquara de Fora Inlet - NCAAA Outtake Point',zorder=2)

x,y=m(picinguaba[0]+0.05, picinguaba[1]+0.015)
ax.text(x,y,u'Picinguaba',color='k',fontsize=13, ha='center',va='center')

x,y=m(guaratiba[0]+0.035, guaratiba[1]+0.025)
ax.text(x,y,u'Barra de\nGuaratiba',color='k',fontsize=13, ha='center',va='center')

x,y=m(MambucabaArrow[0],MambucabaArrow[1])
ax.arrow(x,y,MambucabaArrow[2],MambucabaArrow[3],head_width=0.5,head_length=0.1,fc='k',ec='k')
ax.text(x+MambucabaArrow[4],y+MambucabaArrow[5],u'Mambucaba\nRiver Mouth',fontsize=13,ha='center',va='center')

lg = plt.legend(loc='upper right',fontsize=12,numpoints=1,scatterpoints=1)
lg.get_frame().set_alpha(.4)


plt.show()
