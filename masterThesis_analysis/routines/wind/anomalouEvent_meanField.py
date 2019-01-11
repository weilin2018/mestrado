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
plt.ion()

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def make_map(ax,llat=-30,ulat=-20,llon=-50,ulon=-39,resolution='l',nmeridians=3,nparallels=2,labels=[True,False,False,True]):

    m = Basemap(projection='merc', llcrnrlat=llat, urcrnrlat=ulat, llcrnrlon=llon, urcrnrlon=ulon, resolution=resolution)

    m.ax = ax

    m.drawcoastlines(linewidth=.1)
    m.drawmapboundary()
    m.fillcontinents(color='#c0c0c0')
	# definir meridianos e paralelos para plotar no mapa
    meridians=np.arange(llon,ulon,nmeridians)
    parallels=np.arange(llat,ulat,nparallels)
	# desenhar meridianos e paralelos conforme definido acima
    m.drawparallels(parallels,labels=labels,fontsize=8,color='gray',linewidth=.2)
    m.drawmeridians(meridians,labels=labels,fontsize=8,color='gray',linewidth=.2)

    return m

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/data')

ncin = xr.open_dataset(DATA_DIR+'wind_20140109_20140301.nc')

# calculating speed based on U and V components
spd = np.sqrt(ncin['U_GRD_L103']**2 + ncin['V_GRD_L103']**2)
lon = ncin.lon.values - 360
lat = ncin.lat.values
# mean in time
sme = spd.mean(axis=0)
ume = np.mean(ncin.U_GRD_L103.values,axis=0)/sme
vme = np.mean(ncin.V_GRD_L103.values,axis=0)/sme

# gridding coordinates to plot
lon,lat = np.meshgrid(lon,lat)

# plotting speed just to check
fig,ax = plt.subplots(figsize=(8.4/2.54,7./2.54))
plt.title(u'Velocidade Média de 15 de Janeiro\na 13 de Fevereiro de 2014',fontsize=8)
# left bottom width height
cax = fig.add_axes([0.92,0.16,0.03,0.66])

m = make_map(ax)
cf = m.contourf(lon,lat,sme,cmap=cmo.cm.speed,latlon=True)
qv = m.quiver(lon[::3,::2],lat[::3,::2],ume[::3,::2],vme[::3,::2],latlon=True)
cb = plt.colorbar(cf,cax=cax,orientation='vertical')
cb.set_label(u'Velocidade Média',fontsize=8)
#
# rect = (0,0.08,1.,0.95)
# plt.tight_layout(rect=rect) # box for tight_subplot_layout
# plt.subplots_adjust(top=0.87,bottom=0.12,left=0.15,right=0.865,hspace=0.13,wspace=0.0)
#
