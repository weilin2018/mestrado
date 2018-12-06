# fonte: https://machinelearningmastery.com/time-series-trends-in-python/

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


from sklearn.linear_model import LinearRegression

def remove_2std(data):

    limsuperior = data.mean() + 2*data.std()
    liminferior = data.mean() - 2*data.std()
    data[data < liminferior] = np.nan
    data[data > limsuperior] = np.nan

    return data

def load_data():


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

data = load_data()

# selecionando primeiro trecho dos dados até 28/1
series = data[:'2014-01-28'].tempSurf

# fit linear model
X = [i for i in range(0, len(series))]
X = numpy.reshape(X, (len(X), 1))
y = series.values
model = LinearRegression()
model.fit(X, y)
# calculate trend
trend = model.predict(X)

# visualizar
plt.plot(y)
plt.plot(trend)
plt.show()

# selecionando a segunda parte dos dados, a partir de 11/02
series = data['2014-02-12':].tempSurf
series = series.dropna()

# fit linear model
X = [i for i in range(0, len(series))]
X = numpy.reshape(X, (len(X), 1))
y = series.values
model = LinearRegression()
model.fit(X, y)
# calculate trend
trend = model.predict(X)

# visualizar
plt.plot(y)
plt.plot(trend)
plt.show()

# continuar o codigo pensando numa forma de colocar as duas tendencia em cima da mesma série temporal
