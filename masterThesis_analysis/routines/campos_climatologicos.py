#-*-coding;utf-8-*-
"""
    Plotar Figura 2.4 da dissertação
"""

import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import pandas as pd
import os
import pickle
from mpl_toolkits.basemap import Basemap
from matplotlib import dates
import datetime
import cmocean as cmo
import gsw

import matplotlib
matplotlib.use('PS')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here


def make_map(ax,llat=-30,ulat=-20,llon=-50,ulon=-39,resolution='l',nmeridians=3,nparallels=2,labels=[True,False,False,True]):

    if resolution == 'f':
        m = pickle.load(open('pickles/basemap.p','r'))
    else:
        m = Basemap(projection='merc', llcrnrlat=llat, urcrnrlat=ulat, llcrnrlon=llon, urcrnrlon=ulon, resolution=resolution)

    m.ax = ax

    m.drawcoastlines(linewidth=.01,color='#c0c0c0')
    m.drawmapboundary()
    m.fillcontinents(color='#c0c0c0')
    m.drawstates(linewidth=.01,color='gray')
	# definir meridianos e paralelos para plotar no mapa
    meridians=np.arange(llon,ulon,nmeridians)
    parallels=np.arange(llat,ulat,nparallels)
	# desenhar meridianos e paralelos conforme definido acima
    m.drawparallels(parallels,labels=labels,fontsize=8,color='gray',linewidth=.02)
    m.drawmeridians(meridians,labels=labels,fontsize=8,color='gray',linewidth=.02)

    return m

def create_Structure(fname,timestep=0,savefig=False):

    #fig,axes = plt.subplots(nrows=2,ncols=3,figsize=(16/2.54, 13/2.54))
    fig,axes = plt.subplots(nrows=2,ncols=3,figsize=(16/2.54, 13/2.54))
    # cax = fig.add_axes([0.2,0.05,0.61,0.02])

    # dictionary containing labels for subplots
    labels_dict = {
        '00': [True,False,False,False],
        '01': [False,False,False,False],
        '02': [False,False,False,False],
        '10': [True,False,False,True],
        '11': [False,False,False,True],
        '12': [False,False,False,True],
    }

    m = {}

    for j in range(3):
        for i in range(2):
            key = "%s%s"%(i,j)
            m[key] = make_map(axes[i,j],labels=labels_dict[key],ulat=-21,llat=-29,ulon=-40,resolution='f')
            axes[i,j].spines['left'].set_linewidth(0.2)
            axes[i,j].spines['right'].set_linewidth(0.2)
            axes[i,j].spines['bottom'].set_linewidth(0.2)
            axes[i,j].spines['top'].set_linewidth(0.2)

    sigmaLevels = [0,10,20] # which sigma levels to plot

    # fig,axes = plt.subplots(nrows=2,ncols=3,figsize=(16/2.54, 13/2.54))
    # cax = fig.add_axes([0.2,0.05,0.61,0.02])
    axes_pos = [
        [0.18,0.05,0.30,0.02],
        [0.57,0.05,0.28,0.02],
    ]

    cax_temp = fig.add_axes(axes_pos[0])
    cax_salt = fig.add_axes(axes_pos[1])
    #
    # m = {}
    #
    # for i in range(2):
    #     for j in range(3):
    #         key = "%s%s"%(i,j)
    #         m[key] = make_map(axes[i,j])

    # plotting climatologic data: t = 0, k = 0
    ncin = xr.open_dataset(fname)

    lon,lat = ncin.lon.values, ncin.lat.values
    lon[lon == 0.] = np.nan
    lat[lat == 0.] = np.nan

    # extracting temperature data, in a specific timestep
    temp = ncin.temp[timestep,:,:,:]
    contours = np.arange(3,30,1)

    # key for axes in m
    col1 = ['00','01','02']

    for key,k in zip(col1,sigmaLevels):
        a = m[key]
        cf_temp = a.contourf(lon,lat,temp[k,:,:],contours,latlon=True,cmap=cmo.cm.thermal)

    # plotting anomalous experiment at the final
    salt = ncin.salt[timestep,:,:,:]
    contours = np.arange(33,37,0.01)

    col1 = ['10','11','12']
    for key,k in zip(col1,sigmaLevels):
        a = m[key]
        cf_salt = a.contourf(lon,lat,salt[k,:,:],contours,latlon=True,cmap=cmo.cm.haline)

    # axes[0,0].set_title('Temperatura',fontsize=8)
    # axes[0,1].set_title(u'Salinidade',fontsize=8)

    # setting colorbar configuration
    cb_temp = plt.colorbar(cf_temp,orientation='horizontal',cax=cax_temp,format='%i',ticks=np.arange(3,30,6.5))
    from matplotlib import ticker
    tick_locator = ticker.MaxNLocator(nbins=6)
    cb_temp.locator = tick_locator
    cb_temp.update_ticks()
    cb_temp.ax.axes.tick_params(axis='both',which='both',labelsize=8)
    cb_temp.ax.set_title(r'Temperatura ($^o$C)',fontsize=8)

    cb_salt = plt.colorbar(cf_salt,orientation='horizontal',cax=cax_salt,format='%0.1f',ticks=np.arange(33,37,0.9))
    # fig.text(0.25,0.075,r'Temperatura ($^o$C)',fontsize=8)
    fig.text(0.65,0.075,r'Salinidade',fontsize=8)

    # title and some figure adjusts
    d = pd.to_datetime(ncin.time[timestep].values)
    plt.suptitle(u'Climatologia de Temperatura (superior) e Salinidade \n(inferior), para a Superfície, Meio e Fundo',fontsize=10)
    rect = (0,0.08,1.,0.95)
    plt.tight_layout(rect=rect) # box for tight_subplot_layout
    plt.subplots_adjust(top=0.94,bottom=0.09,left=0.06,right=0.99,hspace=0.0,wspace=0.18)

    if savefig:
        plt.savefig('/media/danilo/Danilo/mestrado/github/masterThesis_analysis/figures/climatologia_temp_salt_EDITARINKSCAPE.eps')

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
exp = 'warmup_plot.cdf'
SAVE_FIG = BASE_DIR + 'masterThesis_analysis/figures/experiments_outputs/temperature/crossSection_EA1/'

for f in fname:
    if exp in f:
        experiment = f

fname = experiment

timestep = 0#input('Type which timestep to plot: ')

fig,axes = create_Structure(fname,timestep=int(timestep),savefig=False)
plt.savefig(BASE_DIR + "/masterThesis_analysis/figures/climatologia_temp_salt_EDITARINKSCAPE.eps")
