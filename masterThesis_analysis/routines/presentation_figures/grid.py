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

# opening netcdf file
ncin = xr.open_dataset(fname)

# extracting variables
lon = ncin.lon.values
lat = ncin.lat.values
dep = ncin.depth.values
lon[lon == 0.0] = np.nan
lat[lat == 0.0] = np.nan

# create figure instance
fig,ax = plt.subplots()
m = make_map(ax,resolution='f')
m.drawstates()
# plot depth as contourf and contour of 100, 200, 1000, and 2000 meters depth
cf = m.contourf(lon,lat,dep*(-1),np.arange(-3000,0,10),latlon=True,cmap=cmo.cm.deep_r)
cr = m.contour(lon,lat,dep*(-1),colors=('k'),levels=[-2000,-1000,-200,-100],latlon=True)
# plot grid
m.plot(lon,lat,'k',alpha=.1,latlon=True)
m.plot(lon.T,lat.T,'k',alpha=.1,latlon=True)
# saving figure in presentation's folder, then close
plt.savefig('/media/danilo/Danilo/mestrado/gitlab/mestrado/dissertacao/presentation/figures/numericalGrid.png',dpi=300,bbox_inches='tight',transparent=True)
plt.close()
