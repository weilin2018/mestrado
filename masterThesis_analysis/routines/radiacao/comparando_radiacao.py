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


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/data/radiacao_ondacurta')

ncin = xr.open_dataset(BASE_DIR+'masterThesis_analysis/routines/radiacao/radiacao.nc')
df   = ncin.to_dataframe()

ts=[]
for t in df.index.values:
    ts.append(t[0])
df.index = ts

# importar dados de radiacao do CFSv2
ncin = xr.open_dataset(DATA_DIR+"ubatuba.nc")
cfsv2 = pd.DataFrame({'UBATUBA_reanalise':np.squeeze(ncin.DSWRF_L1.values)},index=pd.DatetimeIndex(np.squeeze(ncin.time.values)))

ncin = xr.open_dataset(DATA_DIR+"paraty.nc")
cfsv2['PARATY_reanalise'] = np.squeeze(ncin.DSWRF_L1.values)

# selecionando os dados
reanalise = cfsv2['2013':'2017'].copy()
# reamostrando dados observados para bater com os dados de reanalise
observado = df.resample('6H').mean()

# plotando dados para cada localizacao, mas passando um filtro quinzenal pra facilitar
fig,ax = plt.subplots(nrows=2)

reanalise.UBATUBA_reanalise.resample('D').mean().plot(ax=ax[0])
observado.UBATUBA_rad.resample('D').mean().plot(ax=ax[0])

reanalise.PARATY_reanalise.resample('D').mean().plot(ax=ax[1])
observado.PARATY_rad.resample('D').mean().plot(ax=ax[1])
