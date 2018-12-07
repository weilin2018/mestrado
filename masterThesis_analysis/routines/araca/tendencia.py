# analise de tendencia dos dados observados (dottori2015) e modelados (EA1)
# gerando as figuras tendencia_araca_2015.eps e tendencia_EA1_araca_2015.eps
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

def load_data_modelado():
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

    return temp,salt

def viewTrend(y,trend):
    plt.plot(y)
    plt.plot(trend)
    plt.show()

data = load_data_observado()

# selecionando primeiro trecho dos dados até 28/1
series = data[:'2014-01-28'].tempSurf

# fit linear model
X = [i for i in range(0, len(series))]
X = np.reshape(X, (len(X), 1))
y = series.values
model = LinearRegression()
model.fit(X, y)
# calculate trend
trend_first = model.predict(X)
time_first  = series.index

# selecionando a segunda parte dos dados, a partir de 11/02
series = data['2014-02-12':].tempSurf
series = series.dropna()

# fit linear model
X = [i for i in range(0, len(series))]
X = np.reshape(X, (len(X), 1))
y = series.values
model = LinearRegression()
model.fit(X, y)
# calculate trend
trend_second = model.predict(X)
time_second  = series.index

# continuar o codigo pensando numa forma de colocar as duas tendencia em cima da mesma série temporal

fig,ax = plt.subplots(figsize=(15,4))
# plotting all time series
ax.plot(data.index,data.tempSurf.values,'#c0c0c0')

# plotting first trend
ax.plot(time_first,trend_first,'k')

# plotting 2nd trend
ax.plot(time_second,trend_second,'k')

ax.axvline('2014-01-15',color='black',alpha=.3,linestyle='dashed')
ax.axvline('2014-02-14',color='black',alpha=.3,linestyle='dashed')

ax.set_ylabel(r'Temperatura ($^o$C)',fontsize=8)
ax.set_xlabel(u'Tempo (dias)',fontsize=8)
plt.suptitle(u'Tendência de aquecimento e resfriamento da Temperatura da Superfície do Mar - Dados observados na Baía do Araçá',fontsize=10)

plt.tight_layout()
plt.subplots_adjust(top=0.924,bottom=0.116,left=0.042,right=0.988,hspace=0.2,wspace=0.2)

plt.savefig(BASE_DIR+'masterThesis_analysis/figures/dados_observados/tendencia_araca_2015.eps')
plt.show()


##### dados modelados
temp,__ = load_data_modelado()


# selecionando primeiro trecho dos dados até 28/1
series = temp[:'2014-01-28'].Surface

# fit linear model
X = [i for i in range(0, len(series))]
X = np.reshape(X, (len(X), 1))
y = series.values
model = LinearRegression()
model.fit(X, y)
# calculate trend
trend_first = model.predict(X)
time_first  = series.index

# selecionando a segunda parte dos dados, a partir de 11/02
series = temp['2014-02-12':].Surface
series = series.dropna()

# fit linear model
X = [i for i in range(0, len(series))]
X = np.reshape(X, (len(X), 1))
y = series.values
model = LinearRegression()
model.fit(X, y)
# calculate trend
trend_second = model.predict(X)
time_second  = series.index

# continuar o codigo pensando numa forma de colocar as duas tendencia em cima da mesma série temporal

fig,ax = plt.subplots(figsize=(15,4))
# plotting all time series
ax.plot(temp.index,temp.Surface.values,'#c0c0c0')

# plotting first trend
ax.plot(time_first,trend_first,'k')

# plotting 2nd trend
ax.plot(time_second,trend_second,'k')

ax.axvline('2014-01-15',color='black',alpha=.3,linestyle='dashed')
ax.axvline('2014-02-14',color='black',alpha=.3,linestyle='dashed')

ax.set_ylabel(r'Temperatura ($^o$C)',fontsize=8)
ax.set_xlabel(u'Tempo (dias)',fontsize=8)
plt.suptitle(u'Tendência de aquecimento e resfriamento da Temperatura da Superfície do Mar - Produto EA1',fontsize=10)

plt.tight_layout()
plt.subplots_adjust(top=0.924,bottom=0.116,left=0.042,right=0.988,hspace=0.2,wspace=0.2)

plt.savefig(BASE_DIR+'masterThesis_analysis/figures/experiments_outputs/araca/tendencia_EA1_araca_2015.eps')
