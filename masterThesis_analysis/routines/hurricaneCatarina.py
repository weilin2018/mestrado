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

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

os.system('clear')
BASE_DIR = oceano.make_dir()
SAVE_DIR = BASE_DIR + 'masterThesis_analysis/routines/pickles/'
FIGU_DIR = BASE_DIR + 'masterThesis_analysis/figures/'
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/data/Furacao_Catarina/')

#################################
#           PART I              #
#################################

# reading and opening ncfile
nfiles = glob.glob(DATA_DIR+"*.nc")[0]
ncdata = xr.open_dataset(nfiles)

# extracting data from ncfile
wu = ncdata['U_GRD_L103'].data
wv = ncdata['V_GRD_L103'].data
time = ncdata['time'].data
lon = ncdata['lon'].data - 360
lat = ncdata['lat'].data

def plot_evolution(lon,lat,wu,wv,time):

    lon,lat = np.meshgrid(lon,lat)

    # create basemap instance with new limits
    fig,ax = plt.subplots(nrows=2,ncols=2)

    m = Basemap(projection='merc', llcrnrlat=lat.min(), urcrnrlat=lat.max(), llcrnrlon=lon.min(), urcrnrlon=lon.max(), resolution='l')
    m.ax = ax[0,0]
    m.drawcoastlines(linewidth=0.8)
    u=wu[6,:,:]
    v=wv[6,:,:]
    s=np.sqrt(u*u + v*v)

    m.contourf(lon,lat,s,latlon=True)
    m.quiver(lon[::2,::2],lat[::2,::2],u[::2,::2],v[::2,::2],latlon=True)

    ax[0,0].set_title(str(time[6]))

    m = Basemap(projection='merc', llcrnrlat=lat.min(), urcrnrlat=lat.max(), llcrnrlon=lon.min(), urcrnrlon=lon.max(), resolution='l')
    m.ax = ax[0,1]
    m.drawcoastlines(linewidth=0.8)
    u=wu[12,:,:]
    v=wv[12,:,:]
    s=np.sqrt(u*u + v*v)

    m.contourf(lon,lat,s,latlon=True)
    m.quiver(lon[::2,::2],lat[::2,::2],u[::2,::2],v[::2,::2],latlon=True)

    ax[0,1].set_title(str(time[12]))

    m = Basemap(projection='merc', llcrnrlat=lat.min(), urcrnrlat=lat.max(), llcrnrlon=lon.min(), urcrnrlon=lon.max(), resolution='l')
    m.ax = ax[1,0]
    m.drawcoastlines(linewidth=0.8)
    u=wu[18,:,:]
    v=wv[18,:,:]
    s=np.sqrt(u*u + v*v)

    m.contourf(lon,lat,s,latlon=True)
    m.quiver(lon[::2,::2],lat[::2,::2],u[::2,::2],v[::2,::2],latlon=True)

    ax[1,0].set_title(str(time[18]))

    m = Basemap(projection='merc', llcrnrlat=lat.min(), urcrnrlat=lat.max(), llcrnrlon=lon.min(), urcrnrlon=lon.max(), resolution='l')
    m.ax = ax[1,1]
    m.drawcoastlines(linewidth=0.8)
    u=wu[-1,:,:]
    v=wv[-1,:,:]
    s=np.sqrt(u*u + v*v)

    m.contourf(lon,lat,s,latlon=True)
    m.quiver(lon[::2,::2],lat[::2,::2],u[::2,::2],v[::2,::2],latlon=True)

    ax[1,1].set_title(str(time[-1]))

    plt.suptitle('Hurricane Catarina',fontsize=40)

def plot_animation(lon,lat,wu,wv,time):

    fig,ax = plt.subplots()

    for t in np.arange(0,len(time)):
        m = Basemap(projection='merc', llcrnrlat=lat.min(), urcrnrlat=lat.max(), llcrnrlon=lon.min(), urcrnrlon=lon.max(), resolution='l')

        m.ax = ax

        m.drawcoastlines(linewidth=0.8)


        u=wu[t,:,:]
        v=wv[t,:,:]
        s=np.sqrt(u*u + v*v)

        m.contourf(lon,lat,s,latlon=True)
        m.quiver(lon[::2,::2],lat[::2,::2],u[::2,::2],v[::2,::2],latlon=True)

        plt.pause(2)
        plt.gca().clear()


def evolution(lon,lat,u,v,t,ax):
    lon,lat = np.meshgrid(lon,lat)

    # create basemap instance with new limits
    m = Basemap(projection='merc', llcrnrlat=lat.min(), urcrnrlat=lat.max(), llcrnrlon=lon.min(), urcrnrlon=lon.max(), resolution='l')
    m.ax = ax
    m.drawcoastlines(linewidth=0.8)
    s=np.sqrt(u*u + v*v)

    m.contourf(lon,lat,s,latlon=True)
    m.quiver(lon[::2,::2],lat[::2,::2],u[::2,::2],v[::2,::2],latlon=True)

    m.scatter(lon[13,18],lat[13,18],marker='D',latlon=True)

    ax.set_title(str(t))


def plot_everything(lon,lat,wind,rotated,time):

    import matplotlib.gridspec as gridspec

    gs = gridspec.GridSpec(2,4)

    ax = plt.subplot(gs[0,0])
    evolution(lon,lat,wind['wu'][6,:,:],wind['wv'][6,:,:],time[6],ax)

    ax = plt.subplot(gs[0,1])
    evolution(lon,lat,wind['wu'][12,:,:],wind['wv'][12,:,:],time[12],ax)

    ax = plt.subplot(gs[0,2])
    evolution(lon,lat,wind['wu'][18,:,:],wind['wv'][18,:,:],time[18],ax)

    ax = plt.subplot(gs[0,3])
    evolution(lon,lat,wind['wu'][-1,:,:],wind['wv'][-1,:,:],time[-1],ax)

    ax = plt.subplot(gs[1,:])

    ax.plot(rotated.index,rotated.along.values,label='long shore')
    ax.plot(rotated.index,rotated.cross.values,label='cross shore')

    # vertical lines
    times = [6,12,18,23]
    for i in times:
        ax.axvline(x=time[i],color='k',alpha=.4)

    # horizontal lines
    ax.axhline(color='k',alpha=.4)

    ax.margins(0)
    plt.legend(loc='best')

    plt.suptitle('Hurricane Catarina',fontsize=40)


def visualize_timeseries(rotated):

    # visualizing data
    fig, ax = plt.subplots(nrows=3,ncols=1,sharex=True,figsize=(15,10))

    ax[0].plot(rotated.cross,label='CFSR')
    ax[0].margins(0)
    ax[0].set_ylim(-15,15)
    ax[0].set_title('Cross Shore')

    ax[0].legend(loc='lower left')

    ax[1].plot(rotated.along,label='CFSR')
    ax[1].margins(0)
    ax[1].set_ylim(-15,15)
    ax[1].set_title('Along Shore')

    ax[1].legend(loc='lower left')
    ax[1].axhline(color='k',alpha=.5)

    ax[2] = oceano.stickplot(serie,ax[2])
    ax[2].set_ylim(-0.7,1.)
    ax[2].set_title('CFSR')

    plt.suptitle('CFSR',fontsize=26)


#################################
#           PART II             #
#################################
# extract data
u = np.squeeze(ncdata['U_GRD_L103'].data[:,13,18])
v = np.squeeze(ncdata['V_GRD_L103'].data[:,13,18])
time = ncdata['time'].data

serie       = pd.DataFrame({'wu':u,'wv':v},index=pd.DatetimeIndex(time))
cross,along = oceano.rotateVectors(serie.wu.values,serie.wv.values,40.)
rotated     = pd.DataFrame({'along':along,'cross':cross},index=serie.index)



wind = {'wu':wu,'wv':wv}

plot_everything(lon,lat,wind,rotated,time)



visualize_timeseries(rotated)
