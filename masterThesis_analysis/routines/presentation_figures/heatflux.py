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

    # m.fillcontinents(color='#c0c0c0')
    #
    # m.drawcoastlines(linewidth=.1)
    # m.drawmapboundary()

	# definir meridianos e paralelos para plotar no mapa
    meridians=np.arange(llon,ulon,nmeridians)
    parallels=np.arange(llat,ulat,nparallels)
	# desenhar meridianos e paralelos conforme definido acima
    m.drawparallels(parallels,labels=[True,False,False,True],fontsize=13,fontweight='bold',color='gray',linewidth=0.0)
    m.drawmeridians(meridians,labels=[True,False,False,True],fontsize=13,fontweight='bold',color='gray',linewidth=0.0)

    return m

def calcQnet(ncin):

    DSWRF = ncin['DSWRF_L1'].values
    DLWRF = ncin['DLWRF_L1'].values
    USWRF = ncin['USWRF_L1'].values
    ULWRF = ncin['ULWRF_L1'].values

    qsw = DSWRF - USWRF
    qlw = DLWRF - ULWRF

    lat = ncin['LHTFL_L1'].values
    sen = ncin['SHTFL_L1'].values

    # net heat flux is a sum of the 4 terms from the heat balance equations
    hflx  = qsw + qlw + lat + sen

    return hflx

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
plt.ion()

# configurações do plot
figsize = (17.4/2.54, 10/2.54)

# DATA_DIR = "/media/danilo/Danilo/mestrado/ventopcse/data/CFSv2/thflx_presentation/"
DATA_DIR = '/home/danilo/Dropbox/mestrado/data/data2model/2014/run/hflx/'
fname = DATA_DIR + "cdas1.20140210.sfluxgrbf.grb2.nc"
ncin  = xr.open_dataset(fname)

lon = ncin.lon.values - 360
lat = ncin.lat.values
lon,lat = np.meshgrid(lon,lat)

hflx= calcQnet(ncin)

# select only one time
hflx = hflx[0,:,:]

fig,ax = plt.subplots(figsize=(10,10))
mhflx = make_map(ax,resolution='f')

mhflx.contourf(lon,lat,hflx,cmap=cmo.cm.thermal,latlon=True)
mhflx.drawcoastlines(linewidth=.8)
mhflx.drawmapboundary()


ncin = xr.open_dataset('/home/danilo/Dropbox/mestrado/data/data2model/2014/run/wind/arquivos_semuso/cdas1.20140210.sfluxgrbf.grb2.nc')
lon = ncin.lon.values - 360
lat = ncin.lat.values
lon,lat = np.meshgrid(lon,lat)
wu  = ncin.U_GRD_L103[0,:,:].values
wv  = ncin.V_GRD_L103[0,:,:].values

mhflx.quiver(lon,lat,wu,wv,latlon=True)

plt.savefig('/media/danilo/Danilo/mestrado/gitlab/mestrado/dissertacao/presentation/figures/heatflux.png',dpi=300)
plt.close()
