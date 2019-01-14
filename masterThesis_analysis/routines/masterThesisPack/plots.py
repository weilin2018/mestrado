#-*-coding;utf-8-*-

# arquivo contendo funcoes utilizadas no artigo
import numpy as np
import xray as xr
import pandas as pd
from scipy.spatial import cKDTree
from scipy import signal, fftpack
import scipy
import socket
import matplotlib.pyplot as plt
import glob
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import pickle
import math
import cmocean as cmo

# importing root functions
import masterThesisPack

def make_map(ax,llat=-30,ulat=-20,llon=-50,ulon=-39,resolution='l',nmeridians=3,nparallels=2,labels=[True,False,False,True]):

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

def create_Structure_horizontal(fname,contours,property='temp',timestep=0,savefig=False):
    """Funcao para criar uma estrutura horizontal de subplots, com os campos termohalinos da saida do modelo.

    Parameters
    ----------
    fname : str
        Caminho para o arquivo netcdf a ser lido.
    contours : np.ndarray
        Contours adequado para a propriedade a ser plotada.
    property : str
        Propriedade a ser plotada. Deve ser o mesmo padrao do arquivo netcdf
    timestep : int
        Timestep a ser plotado
    savefig : boolean
        Se deseja salvar ou nao a figura.

    Returns
    -------
    type
        Description of returned object.

    """

    sigmaLevels = [0,10,20]

    P = {
      'temp': 'Temperatura',
      'salt': 'Salinidade'
    }

    colorbarTitle = {
      'temp': r'Temperatura ($^o$C)',
      'salt': 'Salinidade'
    }

    colormap = {
        'temp':cmo.cm.thermal,
        'salt':cmo.cm.haline
    }

    #fig,axes = plt.subplots(nrows=2,ncols=3,figsize=(16/2.54, 13/2.54))
    fig,axes = plt.subplots(nrows=2,ncols=3,figsize=(11.69,8.27))
    cax = fig.add_axes([0.2,0.05,0.61,0.02])

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
            m[key] = make_map(axes[i,j],labels=labels_dict[key],ulat=-21,llat=-29,ulon=-40)
            axes[i,j].spines['left'].set_linewidth(0.2)
            axes[i,j].spines['right'].set_linewidth(0.2)
            axes[i,j].spines['bottom'].set_linewidth(0.2)
            axes[i,j].spines['top'].set_linewidth(0.2)

    # plotting climatologic data: t = 0, k = 0
    ncin = xr.open_dataset(fname)

    lon,lat = ncin.lon.values, ncin.lat.values
    depth = ncin.depth.values
    sigma = ncin.sigma.values
    lon[lon == 0.] = np.nan
    lat[lat == 0.] = np.nan
    depth = ncin.depth.values

    # extracting temperature data, in a specific timestep
    salt = ncin[property][timestep,:,:,:]
    salt = np.where(depth < 100, salt,np.nan)

    # key for axes in m
    col1 = ['00','01','02']

    for key,k in zip(col1,sigmaLevels):
        a = m[key]
        cf = a.contourf(lon,lat,salt[k,:,:],contours,latlon=True,cmap=colormap[property])
        if k == 0:
            cs = a.contour(lon,lat,salt[k,:,:],levels=[36.],latlon=True,colors=('black'),linewidths=(0.5))

    # plotting anomalous experiment at the final
    ncin = xr.open_dataset(fname.replace('EC1','EA1'))
    salt = ncin[property][timestep,:,:,:]
    salt = np.where(depth < 100, salt,np.nan)

    col1 = ['10','11','12']
    for key,k in zip(col1,sigmaLevels):
        a = m[key]
        cf = a.contourf(lon,lat,salt[k,:,:],contours,latlon=True,cmap=colormap[property])
        if k == 0:
            cs = a.contour(lon,lat,salt[k,:,:],levels=[36.],latlon=True,colors=('black'),linewidths=(0.5))

    axes[0,1].set_title('Experimento Controle',fontsize=8)
    axes[1,1].set_title(u'Experimento Anomalo',fontsize=8)

    # setting colorbar configuration
    cb = plt.colorbar(cf,orientation='horizontal',cax=cax,format='%i')
    fig.text(0.45,0.075,colorbarTitle[property],fontsize=8)

    # title and some figure adjusts
    d = pd.to_datetime(ncin.time[timestep].values)
    plt.suptitle(u'%s nas camadas de superficie, meio e fundo, no Experimento\n' \
                  u'Controle (superior) e Anomalo (inferior) em ' \
                  '%s de %s'%(P[property],d.strftime('%d'),d.strftime('%B')),fontsize=10)
    rect = (0,0.08,1.,0.95)
    plt.tight_layout(rect=rect) # box for tight_subplot_layout
    # plt.subplots_adjust(top=0.886,bottom=0.109,left=0.054,right=0.995,hspace=0.0,wspace=0.045)
    plt.subplots_adjust(top=0.915,bottom=0.11,left=0.036,right=0.999,hspace=0.082,wspace=0.061)

    if savefig:
        savefig_dir = masterThesisPack.make_dir()
        # plt.savefig('/media/danilo/Danilo/mestrado/github/masterThesis_analysis/figures/experiments_outputs/temperature/temperatura_superf_meio_fundo_timestep_%s.png'%(str(timestep)),dpi=300)
        if property == 'temp':
          plt.savefig(savefig_dir+'masterThesis_analysis/figures/experiments_outputs/temperature/temperatura_superf_meio_fundo_timestep_%s.eps'%(str(timestep)),orientation='landscape')
        if property == 'salt':
          plt.savefig(savefig_dir+'masterThesis_analysis/figures/experiments_outputs/salinity/salinidade_superf_meio_fundo_timestep_%s.eps'%(str(timestep)),orientation='landscape')

    return fig,axes
