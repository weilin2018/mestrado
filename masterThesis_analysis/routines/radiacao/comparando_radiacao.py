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

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def fill_spaces(df,ref,key1,key2):

    index = df.index.values
    data1 = df[key1].values
    data2 =ref[key2].values

    nans = np.where(np.isnan(data1))
    data1[nans] = data2[nans]

    return pd.DataFrame({'full':data1},index=index)

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/data/radiacao_ondacurta')

ncin = xr.open_dataset(BASE_DIR+'masterThesis_analysis/routines/radiacao/radiacao.nc')
df   = ncin.to_dataframe()

# ts=[]
# for t in df.index.values:
#     ts.append(t[0])
# df.index = ts

# importar dados de radiacao do CFSv2
ncin = xr.open_dataset(DATA_DIR+"ubatuba.nc")
cfsv2 = pd.DataFrame({'UBATUBA_reanalise':np.squeeze(ncin.DSWRF_L1.values)},index=pd.DatetimeIndex(np.squeeze(ncin.time.values)))

ncin = xr.open_dataset(DATA_DIR+"paraty.nc")
cfsv2['PARATY_reanalise'] = np.squeeze(ncin.DSWRF_L1.values)

# selecionando os dados
reanalise = cfsv2['2013':'2017'].copy()
# reamostrando dados observados para bater com os dados de reanalise
observado = df.resample('6H').mean()

# investigando os outliers visualmente
# import seaborn as sns
# sns.boxplot(observado.PARATY_rad)

# removendo dados 3 desvios padrões acima, classificados como espúrios
stdv_ubatuba = observado.UBATUBA_rad.std()
stdv_paraty  = observado.PARATY_rad.std()
# observado.UBATUBA_rad[observado.UBATUBA_rad > 2*stdv_ubatuba] = np.nan
observado.PARATY_rad[observado.PARATY_rad > 2*stdv_paraty] = np.nan

# plotando dados para cada localizacao, mas passando um filtro quinzenal pra facilitar
fig,ax = plt.subplots(nrows=2,figsize=(15/2.54,10/2.54))

ax[0].set_title(u'Radiação (W m' +r'${^2}$'+ u') em Ubatuba (azul) e CFSv2 (vermelho)',fontsize=8)

reanalise.UBATUBA_reanalise.resample('D').mean().plot(ax=ax[0],label='CFSv2')
observado.UBATUBA_rad.resample('D').mean().plot(ax=ax[0],label='Observado')
observado.UBATUBA_rad.resample('1W').mean().plot(ax=ax[0],color='k',label='Weekly Mean')
reanalise.UBATUBA_reanalise.resample('1W').mean().plot(ax=ax[0],label='CFSv2 Weekly')

ax[0].legend(fontsize=8)
ax[0].set_ylim([0,1000])

reanalise.PARATY_reanalise.resample('D').mean().plot(ax=ax[1],label='CFSv2')
observado.PARATY_rad.resample('D').mean().plot(ax=ax[1],label='Observado')
observado.PARATY_rad.resample('1W').mean().plot(ax=ax[1],color='k',label='Weekly Mean')
ax[1].legend(fontsize=8,loc='upper left')
ax[1].set_ylim([0,1000])

# preenchendo lacunas nos dados observado em Ubatuba
uba = fill_spaces(observado,reanalise,'UBATUBA_rad','UBATUBA_reanalise')
uba.columns = ['Ubatuba']

# plotando dados para cada localizacao, mas passando um filtro quinzenal pra facilitar
fig,ax = plt.subplots(figsize=(16/2.54,8/2.54))

ax.set_title(u'Radiação (W m' +r'${^2}$'+ u') em Ubatuba (azul) e CFSv2 (vermelho)',fontsize=8)

# reanalise.UBATUBA_reanalise.resample('D').mean().plot(ax=ax,label='CFSv2')
uba.Ubatuba.resample('D').mean().plot(ax=ax,label='Observado')
uba.Ubatuba.resample('1W').mean().plot(ax=ax,color='k',label='Weekly Mean')
ax.legend(fontsize=8)
ax.set_ylim([0,500])
ax.axes.tick_params(axis='both',which='both',labelsize=8)


# calculando a media no trimestre DJF de todos os anos
def media_trimestre(df):

    medias = pd.DataFrame()

    for y in np.arange(df.index.year[0]+1,2018,1):
        start = '%i-12-01'%(y-1)
        final = '%i-03-01'%(y)

        medias[y] = df[start:final].mean()

    medias = medias.T

    return medias

#
m = media_trimestre(uba)
m.plot(kind='bar')
plt.title('media trimestral DJF')
