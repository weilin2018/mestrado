#-*-coding:utf-8-*-
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
import matplotlib.gridspec as gridspec

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
def make_map(ax,llat=-30,ulat=-20,llon=-50,ulon=-39,resolution='l',nmeridians=3,nparallels=2,labels=[True,False,False,True]):
    # criar mapa sem continente colorido, somente a linha de costa
    # nao colocar meridianos e paralelos
    # nao colocar coordenadas no eixo (sumir com os eixos)

    m = Basemap(projection='merc',llcrnrlat=llat, urcrnrlat=ulat, llcrnrlon=llon, urcrnrlon=ulon, resolution='h')
    # m = pickle.load(open('pickles/basemap.p','r'))
    m.ax = ax
    m.drawcoastlines(linewidth=.2)
    m.fillcontinents(color='white',alpha=0)

    meridians=np.arange(llon,ulon,nmeridians)
    parallels=np.arange(llat,ulat,nparallels)

    return m,meridians,parallels

def export_data(fname,timestep=0):
    # plotting climatologic data: t = 0, k = 0
    ncin = xr.open_dataset(fname)

    lon,lat = ncin.lon.values, ncin.lat.values
    depth = ncin.depth.values
    sigma = ncin.sigma.values
    lon[lon == 0.] = np.nan
    lat[lat == 0.] = np.nan
    depth = ncin.depth.values

    # extracting salinity data, in a specific timestep
    elev = np.nanmean(ncin.elev[timestep,:,:],axis=0)
    elev = np.where(depth < 200, elev,np.nan)

    return lon,lat,elev,depth

def configurePlot(axes):

    axes[0,0].set_title(u'EC1 - 14/Janeiro',fontsize=8)
    axes[1,0].set_title(u'EC2 - 15/Fevereiro',fontsize=8)
    axes[0,1].set_title(u'EA1 - 14/Janeiro',fontsize=8)
    axes[1,1].set_title(u'EA2 - 15/Fevereiro',fontsize=8)

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
fname = DATA_DIR + 'EA1.cdf'

plt.ion()

fig,axes = plt.subplots(ncols=2,nrows=2)
configurePlot(axes)
# timestep = [46,303]
timestep = [np.arange(65,73,1),np.arange(280,289,1)]
contours = np.arange(-.5,.5,0.01)
cont = 0 # axis to plot

# plot EA1
### begin
lon,lat,elev,depth = export_data(fname,timestep=timestep[0])
elev[:5,:] = np.nan
elev[-5:,:] = np.nan

m,meridians,parallels = make_map(axes[0,0],ulon=np.nanmax(lon)-.5,llon=np.nanmin(lon)-.2,ulat=np.nanmax(lat)+.2,llat=np.nanmin(lat))

cf = m.contourf(lon,lat,elev,contours,cmap='RdBu',latlon=True,rasterized=True)

for c in cf.collections:
    c.set_edgecolor('face')
    c.set_linewidth(0.00000000001)

### final
lon,lat,elev,depth = export_data(fname,timestep=timestep[1])
elev[:5,:] = np.nan
elev[-5:,:] = np.nan

m,meridians,parallels = make_map(axes[1,0],ulon=np.nanmax(lon)-.5,llon=np.nanmin(lon)-.2,ulat=np.nanmax(lat)+.2,llat=np.nanmin(lat))

cf = m.contourf(lon,lat,elev,contours,cmap='RdBu',latlon=True,rasterized=True)

for c in cf.collections:
    c.set_edgecolor('face')
    c.set_linewidth(0.00000000001)

# plot EA2
### begin
lon,lat,elev,depth = export_data(fname.replace('1','2'),timestep=timestep[0])
elev[:5,:] = np.nan
elev[-5:,:] = np.nan

m,meridians,parallels = make_map(axes[0,1],ulon=np.nanmax(lon)-.5,llon=np.nanmin(lon)-.2,ulat=np.nanmax(lat)+.2,llat=np.nanmin(lat))

cf = m.contourf(lon,lat,elev,contours,cmap='RdBu',latlon=True,rasterized=True)

for c in cf.collections:
    c.set_edgecolor('face')
    c.set_linewidth(0.00000000001)

### final
lon,lat,elev,depth = export_data(fname.replace('1','2'),timestep=timestep[1])
elev[:5,:] = np.nan
elev[-5:,:] = np.nan

m,meridians,parallels = make_map(axes[1,1],ulon=np.nanmax(lon)-.5,llon=np.nanmin(lon)-.2,ulat=np.nanmax(lat)+.2,llat=np.nanmin(lat))

cf = m.contourf(lon,lat,elev,contours,cmap='RdBu',latlon=True,rasterized=True)

for c in cf.collections:
    c.set_edgecolor('face')
    c.set_linewidth(0.00000000001)
