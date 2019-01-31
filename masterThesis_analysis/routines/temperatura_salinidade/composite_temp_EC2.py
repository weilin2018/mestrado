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

    # extracting temperature data, in a specific timestep
    temp = ncin.temp[timestep,:,:,:]
    temp = np.where(depth < 100, temp,np.nan)

    return lon,lat,temp,depth

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
fname = DATA_DIR + 'EC2.cdf'

timestep = [46,303]

for nstep in timestep:
    plt.close()
    lon,lat,temp,depth = export_data(fname,timestep=nstep)
    # mask first and last 5 rows
    temp[:,:5,:] = np.nan
    temp[:,-5:,:] = np.nan

    fig = plt.figure(figsize=(12/2.54,12/2.54))
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
    ax1.set_position([.03,.46,.6,.5])
    ax2.set_position([.1,.28,.6,.5])
    ax3.set_position([.17,.1,.6,.5])

    contours = np.arange(13,35,0.1)

    cf1 = m1.contourf(lon,lat,temp[0,:,:],contours,cmap=cmo.cm.thermal,latlon=True,rasterized=True,extend='max')
    cf2 = m2.contourf(lon,lat,temp[10,:,:],contours,cmap=cmo.cm.thermal,latlon=True,rasterized=True,extend='max')
    cf3 = m3.contourf(lon,lat,temp[20,:,:],contours,cmap=cmo.cm.thermal,latlon=True,rasterized=True,extend='max')

    # matplotib trick to remove white thin lines when saving contourf in pdf
    for c in cf1.collections:
        c.set_edgecolor('face')
        c.set_linewidth(0.00000000001)

    for c in cf2.collections:
        c.set_edgecolor('face')
        c.set_linewidth(0.00000000001)

    for c in cf3.collections:
        c.set_edgecolor('face')
        c.set_linewidth(0.00000000001)

    # colorbar
    cax = fig.add_axes([.1,.08,.35,.02])
    cbar = plt.colorbar(cf1,orientation='horizontal',cax=cax,format='%i')

    # # figure's title
    # plt.suptitle(u'Temperatura nas camadas de superfície, meio e fundo no Experimento Controle (esquerda) e Anômalo (direita) em 14 de Janeiro')

    # setting colorbar tick labels
    from matplotlib import ticker
    tick_locator = ticker.MaxNLocator(nbins=6)
    cbar.locator = tick_locator
    cbar.update_ticks()

    cbar.ax.axes.tick_params(axis='both',which='both',labelsize=8)
    cbar.ax.set_title(r'Temperatura ($^o$C)',fontsize=8)

    output_fname = fname.split('/')[-1].replace('.cdf','_'+str(nstep))
    plt.savefig('/home/danilo/Pictures/composicao/temp/EC2/%s.eps'%(output_fname))

"""
Nota:

    Apos salvar a figura em .pdf, e' necessario ir ao inkscape,
    realizar o ungroup dos objetos de imagem e deletar algumas coisas como:

    . cor cinza
    . elementos dentro do continente
    . mover colorbar pro local adequado
    . nomear cada contourf
"""
