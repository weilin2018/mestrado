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

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/data/INMET/netcdf')

fname = DATA_DIR + 'inmet_paraty_2006-2017.nc'
ncin  = xr.open_dataset(fname)
# converting Dataset into Dataframe
df = ncin.to_dataframe()

# converting from kJ to W:
df['radiacao'] = df.radiation / 3.6

newdf = df['2013':].copy()
newdf.iloc[np.where(df.radiacao == 0.)] = np.nan
media = newdf.mean()
std   = newdf.std()

# visualizando
fig,ax = plt.subplots(figsize=(15,5))

newdf.radiacao.plot(ax=ax,title=u'Radiação Solar Registrada (W m-2), em vermelho, e Média Semanal, em azul,  entre 2013 e 2018 em Paraty.')

ax.axvline('2014-01-14 00:00',color='k')
ax.axvline('2014-02-15 00:00',color='k')

ax.axvline('2014-01-01 00:00',color='k',linestyle='--')
ax.axvline('2014-03-01 00:00',color='k',linestyle='--')

# ax.margins(0)

# plotando um resample a cada 15 dias
newdf.radiacao.resample('1W').mean().plot(ax=ax)

plt.margins(0)
plt.tight_layout()

plt.savefig('/media/danilo/Danilo/mestrado/github/masterThesis_analysis/figures/dados_observados/INMET/radiacao_paraty.pdf')
