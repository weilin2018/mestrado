# leitura dos dados enviados pelo CEBIMAR e geracao da figura cebimar_all.pdf

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

from dateutil import parser

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

def createDates(dates,times):

    date = []

    for d,h in zip(dates,times):
        t = d + ' ' + h
        date.append(parser.parse(t))

    return date

def struct2x6(nrows=2,ncols=6,figsize=None):
    fig,axes = plt.subplots(nrows=nrows,ncols=ncols,figsize=figsize,sharey=True)
    # customizacoes para um perfil vertical em subplots
    axes[0,0].set_ylabel('Profundidade [m]',fontsize=8)
    axes[1,0].set_ylabel('Profundidade [m]',fontsize=8)
    # set ylim do primeiro subplot
    # axes[0,0].set_ylim([-15,0])
    # axes[1,0].set_ylim([-15,0])
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
        # plotando temperatura
        axes[0,i].plot(data.Temperature.values,-data.index.values,'grey')
        axes[0,i].text(25,-1.5,k.strftime('%H:%M').replace(' ','\n'),horizontalalignment='center',verticalalignment='center',fontsize=8)
        axes[1,i].plot(data.Salinity.values,-data.index.values,'k')
        # axes[1,i].text(35.5,-1.5,k.strftime('%d/%m/%Y %H:%M').replace(' ','\n'),horizontalalignment='center',verticalalignment='center',fontsize=8)

    for i in np.arange(0,6,2):
        k = keys[i]
        axes[0,i].text(31,1.9,k.strftime('%d/%m/%Y').replace(' ','\n'),horizontalalignment='center',verticalalignment='center',fontsize=8)

    plt.tight_layout()
    plt.subplots_adjust(top=0.868,bottom=0.012,left=0.094,right=0.972,hspace=0.123,wspace=0.319)

    plt.suptitle(u'Perfis Verticais no Canal de São Sebastião (CEBIMAR) - [%1.2f,%1.2f]'%(meanLocation[0],meanLocation[1]),fontsize=10)

    plt.savefig(BASE_DIR+'masterThesis_analysis/figures/dados_observados/cebimar_all.pdf',orientation='landscape')


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

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/data/Araca/Aurea_CEBIMAR')

fname = DATA_DIR + 'casts_bined_4Dott.txt'

df = pd.read_csv(fname,sep='	')
df.drop(['Type','Station'],inplace=True,axis=1) # removendo algumas colunas
df.set_index('Cruise',inplace=True) # definindo index

# trocando o nome das colunas para algo mais simples
df.columns = ['Date','Time','Lon','Lat','Depth','Temperature','Salinity']

# converter date e time em datetime para podermos usar como condicao de separacao
dates = df.Date.values
times = df.Time.values

# removendo as antigas colunas de tempo
df.drop(['Date','Time'],inplace=True,axis=1)

date = createDates(dates,times)

# inserindo as novas datas no dataframe
df['datetime'] = date

# selecionando somente os perfis de Jan/2014 (designados pelo Cruise/index 201401)
jan2014 = df.loc[201401]

# definindo a localizacao media dos perfis
meanLocation = [jan2014.Lat.mean(),jan2014.Lon.mean()]

# Selecionando os dados pelo cruzeiro realizado em cada dia e hora
# indexes = jan2014.datetime
# indexes = indexes.drop_duplicates()
# nPerfis = len(indexes) # quantidade de perfis obtidos

# agrupando dados
d = jan2014.drop(['Lon','Lat'],axis=1)
d.set_index('Depth',inplace=True)

# armazenando em caso de necessidade
oldDepths = d.index.values

grouped = d.groupby('datetime')

dct = {}

for name,group in grouped:
    dct[name] = group

keys = dct.keys()
keys.sort()

plot_2x6Figure(dct,keys,BASE_DIR)
