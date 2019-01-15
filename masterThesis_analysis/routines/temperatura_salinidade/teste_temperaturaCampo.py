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
matplotlib.style.use('ggplot')

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

    m = Basemap(projection='merc', llcrnrlat=llat, urcrnrlat=ulat, llcrnrlon=llon, urcrnrlon=ulon, resolution=resolution)
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

    # extracting temperature data, in a specific timestep
    temp = ncin.temp[timestep,:,:,:]
    temp = np.where(depth < 100, temp,np.nan)

    lon,lat = ncin.lon.values, ncin.lat.values
    depth = ncin.depth.values
    sigma = ncin.sigma.values
    lon[lon == 0.] = np.nan
    lat[lat == 0.] = np.nan
    depth = ncin.depth.values

    # extracting temperature data, in a specific timestep
    temp = ncin.temp[timestep,:,:,:]
    temp = np.where(depth < 100, temp,np.nan)

    return lon,lat,temp


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
fname = DATA_DIR + 'EA1.cdf'

lon,lat,temp = export_data(fname)

plt.figure(figsize=(12/2.54,12/2.54))
gs = gridspec.GridSpec(3,3)

# creating axis
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[1,1])
ax3 = plt.subplot(gs[2,2])

# creating basemap instance
m1,meridians,parallels = make_map(ax1,ulon=np.nanmax(lon)-.5,llon=np.nanmin(lon)-.2,ulat=np.nanmax(lat)+.2,llat=np.nanmin(lat))
m2,_,_ = make_map(ax2,ulon=np.nanmax(lon)-.5,llon=np.nanmin(lon)-.2,ulat=np.nanmax(lat)+.2,llat=np.nanmin(lat))
m3,_,_ = make_map(ax3,ulon=np.nanmax(lon)-.5,llon=np.nanmin(lon)-.2,ulat=np.nanmax(lat)+.2,llat=np.nanmin(lat))

# positiong each axes
ax1.set_position([.03,.48,.6,.5])
ax2.set_position([.22,.28,.6,.5])
ax3.set_position([.41,.08,.6,.5])

contours = np.arange(13,35,0.1)

m1.contourf(lon,lat,temp[0,:,:],contours,cmap=cmo.cm.thermal,latlon=True)
m2.contourf(lon,lat,temp[10,:,:],contours,cmap=cmo.cm.thermal,latlon=True)
m3.contourf(lon,lat,temp[20,:,:],contours,cmap=cmo.cm.thermal,latlon=True)

m1.drawmeridians(meridians,labels=[True,False,False,False],fontsize=8,color='gray',linewidth=.2)
m1.drawparallels(parallels,labels=[True,False,False,False],fontsize=8,color='gray',linewidth=.2)

m2.drawmeridians(meridians,labels=[False,False,False,False],fontsize=8,color='gray',linewidth=.2)
m2.drawparallels(parallels,labels=[False,False,False,False],fontsize=8,color='gray',linewidth=.2)

m3.drawmeridians(meridians,labels=[False,False,False,False],fontsize=8,color='gray',linewidth=.2)
m3.drawparallels(parallels,labels=[False,False,False,False],fontsize=8,color='gray',linewidth=.2)
