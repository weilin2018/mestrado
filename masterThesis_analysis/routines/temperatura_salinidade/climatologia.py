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
import gsw

# pacotes para minimap
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

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

    m.drawcoastlines(linewidth=.1)
    m.drawmapboundary()
    m.fillcontinents(color='#c0c0c0')
	# definir meridianos e paralelos para plotar no mapa
    meridians=np.arange(llon,ulon,nmeridians)
    parallels=np.arange(llat,ulat,nparallels)
	# desenhar meridianos e paralelos conforme definido acima
    m.drawparallels(parallels,labels=[True,False,False,True],fontsize=8,color='gray',linewidth=.4)
    m.drawmeridians(meridians,labels=[True,False,False,True],fontsize=8,color='gray',linewidth=.4)

    return m

def create_Structure(fname,timestep=0,savefig=False):

    sigmaLevels = [0,10,20] # which sigma levels to plot

    fig,axes = plt.subplots(nrows=3,ncols=2,figsize=(13/2.54, 16/2.54))
    # cax = fig.add_axes([0.2,0.05,0.61,0.02])
    axes_pos = [
        [0.18,0.05,0.30,0.02],
        [0.57,0.05,0.28,0.02],
    ]

    cax_temp = fig.add_axes(axes_pos[0])
    cax_salt = fig.add_axes(axes_pos[1])

    m = {}

    for i in range(3):
        for j in range(2):
            key = "%s%s"%(i,j)
            m[key] = make_map(axes[i,j])

    # plotting climatologic data: t = 0, k = 0
    ncin = xr.open_dataset(fname)

    lon,lat = ncin.lon.values, ncin.lat.values
    lon[lon == 0.] = np.nan
    lat[lat == 0.] = np.nan

    # extracting temperature data, in a specific timestep
    temp = ncin.temp[timestep,:,:,:]
    contours = np.arange(3,30,1)

    # key for axes in m
    col1 = ['00','10','20']

    for key,k in zip(col1,sigmaLevels):
        a = m[key]
        cf_temp = a.contourf(lon,lat,temp[k,:,:],contours,latlon=True,cmap=cmo.cm.thermal)

    # plotting anomalous experiment at the final
    salt = ncin.salt[timestep,:,:,:]
    contours = np.arange(33,37,0.01)

    col1 = ['01','11','21']
    for key,k in zip(col1,sigmaLevels):
        a = m[key]
        cf_salt = a.contourf(lon,lat,salt[k,:,:],contours,latlon=True,cmap=cmo.cm.haline)

    axes[0,0].set_title('Temperatura',fontsize=8)
    axes[0,1].set_title(u'Salinidade',fontsize=8)

    # setting colorbar configuration
    cb_temp = plt.colorbar(cf_temp,orientation='horizontal',cax=cax_temp,format='%i',ticks=np.arange(3,30,6.5))
    cb_salt = plt.colorbar(cf_salt,orientation='horizontal',cax=cax_salt,format='%0.1f',ticks=np.arange(33,37,0.9))
    fig.text(0.25,0.075,r'Temperatura ($^o$C)',fontsize=8)
    fig.text(0.65,0.075,r'Salinidade',fontsize=8)

    # title and some figure adjusts
    d = pd.to_datetime(ncin.time[timestep].values)
    plt.suptitle(u'Climatologia de Temperatura (à esquerda) e Salinidade \n(à direita), para a Superfície, Meio e Fundo',fontsize=10)
    rect = (0,0.08,1.,0.95)
    plt.tight_layout(rect=rect) # box for tight_subplot_layout
    plt.subplots_adjust(top=0.9,bottom=0.12,left=0.15,right=0.895,hspace=0.13,wspace=0.0)

    if savefig:
        plt.savefig('/media/danilo/Danilo/mestrado/github/masterThesis_analysis/figures/climatologia_temp_salt.png',dpi=300)

    return fig,axes

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
plt.ion()

# configurações do plot
figsize = (17.4/2.54, 10/2.54)

DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
fname = glob.glob(DATA_DIR+"*.cdf")

# select which experiment you want to plot:
exp = 'warmupControle.cdf'
SAVE_FIG = BASE_DIR + 'masterThesis_analysis/figures/experiments_outputs/temperature/crossSection_EA1/'

for f in fname:
    if exp in f:
        experiment = f

fname = experiment

timestep = 0#input('Type which timestep to plot: ')

fig,axes = create_Structure(fname,timestep=int(timestep),savefig=True)
