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

from dateutil import parser

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
def find_timestep(ncin,timestep):

	refDate = parser.parse(timestep)

	time = ncin.time.values
	df   = pd.DataFrame(np.arange(0,len(time)),index=time)

	d = df.iloc[df.index.get_loc(refDate,method='nearest')]

	return d

def extrair_dados(ncin,i,j,t):

	temp = ncin.temp[t,:,i,j]
	salt = ncin.salt[t,:,i,j]

	return temp,salt

def struct2x6(nrows=2,ncols=6,figsize=None):
    fig,axes = plt.subplots(nrows=nrows,ncols=ncols,figsize=figsize,sharey=True)
    # customizacoes para um perfil vertical em subplots
    axes[0,0].set_ylabel('Profundidade [m]',fontsize=8)
    axes[1,0].set_ylabel('Profundidade [m]',fontsize=8)
    # set ylim do primeiro subplot
    axes[0,0].set_ylim([-15,0])
    axes[1,0].set_ylim([-15,0])
    # removendo labelleft dos subplots
    for i in np.arange(1,6,1):
        axes[0,i].tick_params(axis='y',labelleft='off')
        axes[1,i].tick_params(axis='y',labelleft='off')
    # demais configuracoes dos graficos
    for i in range(nrows):
        for j in range(ncols):
            axes[i,j].invert_yaxis()
            axes[i,j].xaxis.tick_top()
            axes[0,j].set_xlim([23,30]) # temperatura
            axes[1,j].set_xlim([34,36]) # salinidade
            axes[i,j].tick_params(axis='both',which='major',labelsize=8)
            axes[i,j].tick_params(axis='x',which='major',pad=.3)
            axes[i,j].margins(y=0)

    return fig,axes

def plot_2x6Figure(dct,keys,BASE_DIR,figsize=(20./2.54,12./2.54)):

    fig,axes = struct2x6(figsize=figsize)

    # definindo titulos
    for i in range(6):
        # extraindo dados do dicionario
        k = keys[i]
        data = dct[k]
        k = parser.parse(keys[i])
        # plotando temperatura
        axes[0,i].plot(data.Temperature.values,data.index.values,'grey')
        axes[0,i].text(25,-1.5,k.strftime('%H:%M').replace(' ','\n'),horizontalalignment='center',verticalalignment='center',fontsize=8)
        axes[1,i].plot(data.Salinity.values,data.index.values,'k')
        # axes[1,i].text(35.5,-1.5,k.strftime('%d/%m/%Y %H:%M').replace(' ','\n'),horizontalalignment='center',verticalalignment='center',fontsize=8)

    for i in np.arange(0,6,2):
        k = parser.parse(keys[i])
        axes[0,i].text(31,1.9,k.strftime('%d/%m/%Y').replace(' ','\n'),horizontalalignment='center',verticalalignment='center',fontsize=8)

    plt.tight_layout()
    plt.subplots_adjust(top=0.868,bottom=0.012,left=0.094,right=0.972,hspace=0.123,wspace=0.319)

    #plt.suptitle(u'Perfis Verticais no Canal de São Sebastião (CEBIMAR) - [%1.2f,%1.2f]'%(meanLocation[0],meanLocation[1]),fontsize=10)
    plt.suptitle(u'Perfis Verticais no Canal de São Sebastião (Modelo)',fontsize=10)

    plt.savefig(BASE_DIR+'masterThesis_analysis/figures/experiments_outputs/araca/EA1_cebimar_all.pdf',orientation='landscape')

def verticalProfile_structure(nrows,ncols,figsize,xlim):

    # criando figura
    fig,axes = plt.subplots(nrows=nrows,ncols=ncols,figsize=figsize)

    # customizacoes para um perfil vertical em subplots
    axes[0,0].set_ylabel('Profundidade [m]',fontsize=8)
    axes[1,0].set_ylabel('Profundidade [m]',fontsize=8)

    # removendo ticklabels da esquerda
    axes[0,1].tick_params(axis='y',labelleft='off') # remove ticklabels
    axes[0,2].tick_params(axis='y',labelleft='off') # remove ticklabels
    axes[1,1].tick_params(axis='y',labelleft='off') # remove ticklabels
    axes[1,2].tick_params(axis='y',labelleft='off') # remove ticklabels

    # x-axis on the top
    for i in range(nrows):
        for j in range(ncols):
            axes[i,j].invert_yaxis()
            axes[i,j].xaxis.tick_top()
            axes[i,j].set_xlim(xlim)
            axes[i,j].tick_params(axis='both',which='major',labelsize=8)
            axes[i,j].tick_params(axis='x',which='major',pad=.3)
            axes[i,j].margins(y=0)

            # axes[i,j].set_ylim([0,-15])

    axes[1,0].tick_params(axis='x',labeltop='off')
    axes[1,1].tick_params(axis='x',labeltop='off')
    axes[1,2].tick_params(axis='x',labeltop='off')

    return fig,axes

def plotMapa(i,j):
    # plotando mapa para visualizar a localizacao e comparacao da profundidade do modelo com a
    # profundidade do ponto observado
    fig,ax = plt.subplots()
    m = oceano.make_map(ax,resolution='f')
    m.plot(lon,lat,'k',alpha=.3,latlon=True);
    m.plot(lon.T,lat.T,'k',alpha=.3,latlon=True);
    cf = m.contourf(lon,lat,ncin.depth.values,np.arange(1,30,0.1),latlon=True,cmap=cmo.cm.deep);
    plt.colorbar(cf)

    m.scatter(lon[i,j],lat[i,j],s=30,c='r',latlon=True)

    m.scatter(-45.40085130,-23.81778822,s=30,c='g',latlon=True)

    plt.show()

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/output')

fname = DATA_DIR + 'EA1.cdf'
ncin  = xr.open_dataset(fname)

# localizando ponto de grade mais proximo do Araca
lon = ncin.lon.values
lat = ncin.lat.values
sig = ncin.sigma.values
dep = ncin.depth.values
lon[lon == 0.] = np.nan
lat[lat == 0.] = np.nan

######################## PERFIL VERTICAL ########################################

ponto = oceano.procurar_pontos_grade(-45.40085130,-23.81778822,lon,lat,n=1)
i     = int(ponto[0][0])
j     = int(ponto[0][1])

j     = 7

profundidade = dep[i,j]*sig

# localizando instantes de tempos correspondentes aos dados obtidos pelo CEBIMAR
instantesCEBIMAR = ['2014-01-14 13:18','2014-01-14 17:15','2014-01-15 12:39','2014-01-15 16:46','2014-01-16 11:09','2014-01-16 16:13']

nstep = []

for instante in instantesCEBIMAR:
    d = find_timestep(ncin,instante)
    nstep.append(d[0])

nstep = np.asarray(nstep)

# extraindo dados termohalinos e armazenando em um dicionario
dct = {}
ins = []

for d in nstep:
    tempo = pd.to_datetime(ncin.time[d].values)
    tempo = tempo.strftime('%Y-%m-%d %H:%M')
    ins.append(tempo)
    t,s = extrair_dados(ncin,i,j,d)
    df = pd.DataFrame({'Temperature':t.values,'Salinity':s.values},index=profundidade)
    dct[tempo] = df

# plot_2x6Figure(dct,ins,BASE_DIR)

###################### SERIE TEMPORAL ##########################################

ponto = oceano.procurar_pontos_grade(-45.403225,-23.816398,lon,lat,n=1)
i     = int(ponto[0][0])
j     = int(ponto[0][1])

# extraindo dados de temperatura e salinidade no ponto observado
localDepth = ncin.depth[i,j].values
time       = ncin.time.values

tSurf = ncin.temp[:,0,i,j]
tBott = ncin.temp[:,-1,i,j]

sSurf = ncin.salt[:,0,i,j]
sBott = ncin.salt[:,-1,i,j]

# passando tudo pra dataframe
temp = pd.DataFrame({'Surface':tSurf,'Bottom':tBott},index=pd.DatetimeIndex(time))
salt = pd.DataFrame({'Surface':sSurf,'Bottom':sBott},index=pd.DatetimeIndex(time))

fig,ax = plt.subplots(nrows=2,figsize=(20./2.54,10./2.54))
#ax[0].plot(salt.Surface,'#c0c0c0',label='Surface Salinity')
ax[0].plot(salt.Surface,'#c0c0c0',label=u'Superfície')
ax[0].plot(salt.Bottom,'k',label='Fundo')

ax[0].legend(loc='best')
ax[0].set_ylabel('Salinidade',fontsize=8)
ax[0].set_ylim([33.35, 35.55])
ax[0].tick_params(axis='x',which='major',labelbottom='off')

ax[1].plot(temp.Surface,'#c0c0c0',label=u'Superfície')
ax[1].plot(temp.Bottom,'#000000',label=u'Fundo')

ax[1].legend(loc='best')
ax[1].set_ylabel('Temperatura',fontsize=8)
ax[1].set_ylim([15.4, 31.1])
ax[1].set_xlabel('Time',fontsize=8)
ax[1].tick_params(axis='x',which='major',rotation=15.)

# colocando duas linhas verticais indicando o comeco e final do evento
ax[0].axvline('2014-01-15',color='black',alpha=.3,linestyle='dashed')
ax[1].axvline('2014-01-15',color='black',alpha=.3,linestyle='dashed')
ax[0].axvline('2014-02-14',color='black',alpha=.3,linestyle='dashed')
ax[1].axvline('2014-02-14',color='black',alpha=.3,linestyle='dashed')

plt.suptitle(u'Série Temporal de Janeiro à Março de 2014\nProduto EA1',fontsize=10)
plt.tight_layout()
plt.subplots_adjust(top=0.885,bottom=0.162,left=0.069,right=0.989,hspace=0.127,wspace=0.2)

plt.savefig(BASE_DIR+'masterThesis_analysis/figures/experiments_outputs/araca/EA1_araca_2015.pdf')
