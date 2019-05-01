# add some description here

import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import pandas as pd
import os
import pickle
from mpl_toolkits.basemap import Basemap
import cmocean as cmo
import datetime

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def make_map(ax,llat=-30,ulat=-20,llon=-50,ulon=-39,resolution='l',nmeridians=3,nparallels=2):

    m = Basemap(projection='merc', llcrnrlat=llat, urcrnrlat=ulat, llcrnrlon=llon, urcrnrlon=ulon, resolution=resolution)

    m.ax = ax

    m.fillcontinents(color='#c0c0c0')

    m.drawcoastlines(linewidth=.1)
    m.drawmapboundary()

	# definir meridianos e paralelos para plotar no mapa
    meridians=np.arange(llon,ulon,nmeridians)
    parallels=np.arange(llat,ulat,nparallels)
	# desenhar meridianos e paralelos conforme definido acima
    m.drawparallels(parallels,labels=[True,False,False,True],fontsize=13,fontweight='bold',color='gray',linewidth=0.0)
    m.drawmeridians(meridians,labels=[True,False,False,True],fontsize=13,fontweight='bold',color='gray',linewidth=0.0)

    return m

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
plt.ion()

# configurações do plot
figsize = (17.4/2.54, 10/2.54)

DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
fname = DATA_DIR + "warmupControle.cdf"

timestep = 0

ncin = xr.open_dataset(fname)

temp = ncin.temp[timestep,0,:,:].values # temperatura
salt = ncin.salt[timestep,0,:,:].values # salinidade
lon = ncin.lon.values
lat = ncin.lat.values
lon[lon == 0] = np.nan
lat[lat == 0] = np.nan

fig,ax = plt.subplots()
mTemp = make_map(ax)

mTemp.contourf(lon,lat,temp,np.arange(21,29,0.001),cmap=cmo.cm.thermal,latlon=True)
plt.savefig('/media/danilo/Danilo/mestrado/gitlab/mestrado/dissertacao/presentation/figures/temp_climatology.png',dpi=300,bbox_inches='tight',transparent=True)
plt.close()

fig,ax = plt.subplots()
mSalt = make_map(ax)

mSalt.contourf(lon,lat,salt,np.arange(33,37,0.01),cmap=cmo.cm.haline,latlon=True)
plt.savefig('/media/danilo/Danilo/mestrado/gitlab/mestrado/dissertacao/presentation/figures/salt_climatology.png',dpi=300,bbox_inches='tight',transparent=True)
plt.close()
