# reading data from Araca Bay, sended by Marcelo Dottori

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
import seawater as sw

import matplotlib
# matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/data/Araca/Dottori_etal_2015')

# importando primeira parte dos dados
ncin = xr.open_dataset(DATA_DIR + 'dados_tratados.nc')

df_1 = ncin.to_dataframe()

# criando um dataframe intermediario, para o periodo que nao houve registro de dados
dt_space = pd.date_range(start='2014-02-05 11:30:00',end='2014-02-11 11:30:00',freq='30min')
n = np.zeros(len(dt_space))*np.nan
dct = {'saltSurf': n,'tempSurf': n,'saltBott': n,'tempBott': n}

df_2 = pd.DataFrame(dct,index=pd.DatetimeIndex(dt_space))

# importando segunda partes dos dados
ncin = xr.open_dataset(DATA_DIR+'dados_nao_tratados.nc')
df_3 = ncin.to_dataframe()

# combinando tudo

df = pd.concat([df_1,df_2,df_3],axis=0)


# visualizando os dados
fig,ax = plt.subplots(nrows=2,figsize=(20./2.54,10./2.54))

# primeiro salinidade
ax[0].set_title('Salinidade',fontsize=8)
df.saltSurf.plot(ax=ax[0],color='gray',label=u'Superfície',linewidth=.6)
df.saltBott.plot(ax=ax[0],color='black',label='Fundo',linewidth=.6)

# segundo temperatura
ax[1].set_title('Temperatura',fontsize=8)
df.tempSurf.plot(ax=ax[1],color='gray',label=u'Superfície',linewidth=.6)
df.tempBott.plot(ax=ax[1],color='black',label='Fundo',linewidth=.6)

# legendas
ax[0].legend(loc='best')
ax[1].legend(loc='best')


# colocando duas linhas verticais indicando o comeco e final do evento
ax[0].axvline('2014-01-15',color='black',alpha=.3,linestyle='dashed')
ax[1].axvline('2014-01-15',color='black',alpha=.3,linestyle='dashed')
ax[0].axvline('2014-02-14',color='black',alpha=.3,linestyle='dashed')
ax[1].axvline('2014-02-14',color='black',alpha=.3,linestyle='dashed')

# y-labels
ax[0].set_ylabel(u'Salinidade',fontsize=8)
ax[1].set_ylabel(u'Temperatura',fontsize=8)
# ax[1].set_xlabel(u'Tempo (dias)',fontsize=8)

ax[0].tick_params(axis='both',which='minor',labelsize=8)
ax[1].tick_params(axis='both',which='minor',labelsize=8)

ax[0].tick_params(axis='y',pad=.2)
ax[1].tick_params(axis='y',pad=.2)

# removendo as margens
ax[0].margins(0)
ax[1].margins(0)

# removendo algumas informaçoes do grafico
ax[0].tick_params(axis='x',labelbottom='off') # remove ticklabels
plt.xticks(fontsize=8)

# ajustes finais
plt.tight_layout()
plt.subplots_adjust(top=0.855,bottom=0.142,left=0.069,right=0.989,hspace=0.357,wspace=0.2)

plt.suptitle(u'Série temporal de Jan a Mar de 2014, na Baáa do Araçá',fontsize=10)

plt.savefig(BASE_DIR+'masterThesis_analysis/figures/dados_observados/araca_2015.eps',orientation='landscape')
