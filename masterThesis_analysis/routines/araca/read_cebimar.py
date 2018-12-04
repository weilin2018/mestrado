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

# plotando
fig,ax = verticalProfile_structure(nrows=2,ncols=3,figsize=(15./2.54,12./2.54),xlim=[23,30])
plt.suptitle(u'Perfis Verticais no Canal de São Sebastião (CEBIMAR) - [%1.2f,%1.2f]'%(meanLocation[0],meanLocation[1]),fontsize=10)

ax[0,0].plot(dct[keys[0]].Temperature.values,dct[keys[0]].index.values,label=keys[0])
ax[0,0].text(25.,1.5,str(keys[0]).replace(' ','\n'),horizontalalignment='center',verticalalignment='center',fontsize=8)

ax[0,1].plot(dct[keys[1]].Temperature.values,dct[keys[1]].index.values,label=keys[1])
ax[0,1].text(25.,1.5,str(keys[1]).replace(' ','\n'),horizontalalignment='center',verticalalignment='center',fontsize=8)

ax[0,2].plot(dct[keys[2]].Temperature.values,dct[keys[2]].index.values,label=keys[2])
ax[0,2].text(25.,1.5,str(keys[2]).replace(' ','\n'),horizontalalignment='center',verticalalignment='center',fontsize=8)

ax[1,0].plot(dct[keys[3]].Temperature.values,dct[keys[3]].index.values,label=keys[3])
ax[1,0].text(25.,1.5,str(keys[3]).replace(' ','\n'),horizontalalignment='center',verticalalignment='center',fontsize=8)

ax[1,1].plot(dct[keys[4]].Temperature.values,dct[keys[4]].index.values,label=keys[4])
ax[1,1].text(25.,1.5,str(keys[4]).replace(' ','\n'),horizontalalignment='center',verticalalignment='center',fontsize=8)

ax[1,2].plot(dct[keys[5]].Temperature.values,dct[keys[5]].index.values,label=keys[5])
ax[1,2].text(25.,1.5,str(keys[5]).replace(' ','\n'),horizontalalignment='center',verticalalignment='center',fontsize=8)

#plt.savefig(BASE_DIR+'masterThesis_analysis/figures/dados_observados/cebimar_temperatura.eps',orientation='landscape')

# plotando
fig,ax = verticalProfile_structure(nrows=2,ncols=3,figsize=(15./2.54,12./2.54),xlim=[34,36])
plt.suptitle(u'Perfis Verticais no Canal de São Sebastião (CEBIMAR) - [%1.2f,%1.2f]'%(meanLocation[0],meanLocation[1]),fontsize=10)

ax[0,0].plot(dct[keys[0]].Salinity.values,dct[keys[0]].index.values,label=keys[0])
ax[0,0].text(35.5,1.5,str(keys[0]).replace(' ','\n'),horizontalalignment='center',verticalalignment='center',fontsize=8)

ax[0,1].plot(dct[keys[1]].Salinity.values,dct[keys[1]].index.values,label=keys[1])
ax[0,1].text(35.5,1.5,str(keys[1]).replace(' ','\n'),horizontalalignment='center',verticalalignment='center',fontsize=8)

ax[0,2].plot(dct[keys[2]].Salinity.values,dct[keys[2]].index.values,label=keys[2])
ax[0,2].text(35.5,1.5,str(keys[2]).replace(' ','\n'),horizontalalignment='center',verticalalignment='center',fontsize=8)

ax[1,0].plot(dct[keys[3]].Salinity.values,dct[keys[3]].index.values,label=keys[3])
ax[1,0].text(35.5,1.5,str(keys[3]).replace(' ','\n'),horizontalalignment='center',verticalalignment='center',fontsize=8)

ax[1,1].plot(dct[keys[4]].Salinity.values,dct[keys[4]].index.values,label=keys[4])
ax[1,1].text(35.5,1.5,str(keys[4]).replace(' ','\n'),horizontalalignment='center',verticalalignment='center',fontsize=8)

ax[1,2].plot(dct[keys[5]].Salinity.values,dct[keys[5]].index.values,label=keys[5])
ax[1,2].text(35.5,1.5,str(keys[5]).replace(' ','\n'),horizontalalignment='center',verticalalignment='center',fontsize=8)

#plt.savefig(BASE_DIR+'masterThesis_analysis/figures/dados_observados/cebimar_salinidade.eps',orientation='landscape')
