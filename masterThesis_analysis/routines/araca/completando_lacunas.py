# extracao de dados do ghrsst proximo ao ponto da baia do araca

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


def load_data_observado():

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

	# controle de qualidade, principalmente para salinidade que esta com valores
	# absurdos
	df.saltBott['2014-01-20 22:30':'2014-02-05'] = np.nan
	df.saltSurf['2014-01-20 22:30':'2014-02-05'] = np.nan

	df.saltBott = remove_2std(df.saltBott)
	df.saltSurf = remove_2std(df.saltSurf)

	return df

################# importar dados observados ###############################
data = load_data_observado()
obs  = data.tempSurf
obs = obs.resample('1D').mean()

################## importar dados do satélite #############################
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/data/GHRSST')

fname = DATA_DIR + 'ghrsst_summer2014.nc'
ncin  = xr.open_dataset(fname)

# localizando ponto de grade mais proximo do Araca
lon = ncin.lon.values
lat = ncin.lat.values

ilon = -45.403225
ilat = -23.816398

i = np.where(lon == oceano.find_nearest_1D(lon,ilon))[0][0]
j = np.where(lat == oceano.find_nearest_1D(lat,ilat))[0][0]

sst = pd.DataFrame({'sst':ncin.analysed_sst[:,j,i].values - 273.15},index=pd.DatetimeIndex(ncin.time.values))

###################### merge nos dados, criando um novo dataframe com tudo #
extracted = obs.values
t         = sst['2014-01-29':'2014-02-14'].values
extracted[np.where(np.isnan(extracted))] = np.squeeze(t[:12])

data = load_data_observado()
obs  = data.tempSurf
obs = obs.resample('1D').mean()

completed = pd.DataFrame({'completed':extracted,'observed':obs.values},index=obs.index)

completed.plot()
plt.show()
