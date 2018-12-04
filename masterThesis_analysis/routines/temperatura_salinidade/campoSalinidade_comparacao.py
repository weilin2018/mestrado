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

import plots as plotDanilo

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

def create_Structure(fname,timestep=0,savefig=False):

    sigmaLevels = [0,10,20] # which sigma levels to plot

    fig,axes = plt.subplots(nrows=3,ncols=2,figsize=(13/2.54, 16/2.54))
    cax = fig.add_axes([0.2,0.05,0.61,0.02])

    m = {}

    for i in range(3):
        for j in range(2):
            key = "%s%s"%(i,j)
            m[key] = make_map(axes[i,j])

    contours = np.arange(34,36.1,0.01)

    # plotting climatologic data: t = 0, k = 0
    ncin = xr.open_dataset(fname)

    lon,lat = ncin.lon.values, ncin.lat.values
    depth = ncin.depth.values
    sigma = ncin.sigma.values
    lon[lon == 0.] = np.nan
    lat[lat == 0.] = np.nan
    depth = ncin.depth.values

    # extracting temperature data, in a specific timestep
    salt = ncin.salt[timestep,:,:,:]
    salt = np.where(depth < 100, salt,np.nan)

    # key for axes in m
    col1 = ['00','10','20']

    for key,k in zip(col1,sigmaLevels):
        a = m[key]
        cf = a.contourf(lon,lat,salt[k,:,:],contours,latlon=True,cmap=cmo.cm.haline)
        if k == 0:
            cs = a.contour(lon,lat,salt[k,:,:],levels=[36.],latlon=True,colors=('black'),linewidths=(0.5))

    # plotting anomalous experiment at the final
    ncin = xr.open_dataset(fname.replace('EC1','EA1'))
    salt = ncin.salt[timestep,:,:,:]
    salt = np.where(depth < 100, salt,np.nan)

    col1 = ['01','11','21']
    for key,k in zip(col1,sigmaLevels):
        a = m[key]
        cf = a.contourf(lon,lat,salt[k,:,:],contours,latlon=True,cmap=cmo.cm.haline)
        if k == 0:
            cs = a.contour(lon,lat,salt[k,:,:],levels=[36.],latlon=True,colors=('black'),linewidths=(0.5))

    axes[0,0].set_title('Experimento Controle',fontsize=8)
    axes[0,1].set_title(u'Experimento Anômalo',fontsize=8)

    # setting colorbar configuration
    cb = plt.colorbar(cf,orientation='horizontal',ticks=[34,34.3,34.6,34.9,35.2,35.5,35.8,36.1],cax=cax,format='%2.1f')
    fig.text(0.55,0.075,r'Salinidade',fontsize=8)

    # title and some figure adjusts
    d = pd.to_datetime(ncin.time[timestep].values)
    plt.suptitle(u'Salinidade nas camadas de superfície, meio e fundo, \n' \
                  u'no Experimento Controle (esquerda) e Anômalo (Direita). \n' \
                  '%s de %s'%(d.strftime('%d'),d.strftime('%B')),fontsize=10)
    rect = (0,0.08,1.,0.95)
    plt.tight_layout(rect=rect) # box for tight_subplot_layout
    plt.subplots_adjust(top=0.87,bottom=0.12,left=0.15,right=0.865,hspace=0.13,wspace=0.0)

    if savefig:
        plt.savefig('/media/danilo/Danilo/mestrado/github/masterThesis_analysis/figures/experiments_outputs/salinity/salinidade_superf_meio_fundo_timestep_%s.png'%(str(timestep)),dpi=300)

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
exp = 'EC1.cdf'

for f in fname:
    if exp in f:
        experiment = f

fname = experiment

# salinidade
contours = np.arange(34,36.1,0.01)

timestep = input('Type which timestep to plot: ')

if timestep == 999.:
    timestep = [0,46,303]

    for nstep in timestep:
        fig,axes = plotDanilo.create_Structure_horizontal(fname,contours,property='salt',timestep=int(nstep),savefig=True)
    plt.close()
else:
    fig,axes = plotDanilo.create_Structure_horizontal(fname,contours,property='salt',timestep=int(timestep),savefig=True)
