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

from dateutil import parser

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def converting_datetime(dates,hours):
	date = []

	for d,h in zip(dates,hours):
		if h == 0:
			date.append(str(d))
		else:
			date.append(str(d)+' '+str(h))

	return pd.to_datetime(date)

def plot_allSerie(df,newdf):
	# visualizando
	fig,ax = plt.subplots(figsize=(15,5))

	df.plot(ax=ax,title=u'Radiação Solar Registrada (W m-2), em vermelho, e Média Semanal, em azul,  entre 2013 e 2018 em Ubatuba.')

	ax.axvline('2014-01-14 00:00',color='k')
	ax.axvline('2014-02-15 00:00',color='k')

	ax.axvline('2014-01-01 00:00',color='k',linestyle='--')
	ax.axvline('2014-03-01 00:00',color='k',linestyle='--')

	# ax.margins(0)

	# plotando um resample a cada 15 dias
	newdf.radiacao.resample('1W').mean().plot(ax=ax,color='k')

	plt.margins(0)
	plt.tight_layout()


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/data/IOUSP/labdados')

data = pd.Series()
time = pd.Series()

for year in ['2013','2014','2015','2016','2017','2018']:
	fList = glob.glob(DATA_DIR+'u/um%s*.csv'%(year))
	fList.sort()

	# variaveis auxiliares
	aux1 = pd.DataFrame()

	for fname in fList:
		ncin = pd.read_csv(fname,sep=';',header=0,usecols=[0,2])
		ncin.columns = ['datetime','radiacao']
		# ncin.index = pd.DatetimeIndex(ncin.datetime)
		# ncin.drop('datetime',axis=1,inplace=True)

		# replacing comma by dots
		ncin.radiacao = ncin.radiacao.apply(lambda x : str(x).replace(',','.'))

		# converting radiacao data from object to float
		ncin.radiacao = ncin.radiacao.apply(float)

		aux1 = pd.concat([aux1,ncin])

	data = pd.concat([data,aux1.radiacao])
	time = pd.concat([time,aux1.datetime])

# criando dataframe com toda a série
df = pd.DataFrame({'radiacao':data.values},index=pd.DatetimeIndex(time.values))

# removendo qualquer valor negativo
df.loc[np.where(df.radiacao < 0)] = np.nan

# removendo dados faltantes em 2013
df['2013-06-03 11:45':'2013-09-11 00:15'] = np.nan

# corrigindo os dados aparentemente errados de 2015: parece que registrou Temperatura no lugar
df['2015-04-01 00:00':'2015-06-30 23:00'] = np.nan

# calculando a media sem os pontos zeros
newdf = df.copy()
newdf.iloc[np.where(df.radiacao == 0.)] = np.nan
media = newdf.mean()
std   = newdf.std()


# visualizando
fig,ax = plt.subplots(figsize=(15,5))

titulo = u'Radiação Solar Registrada (' + r'$W m^{-2}$' + u'), em vermelho, e Média Semanal em azul,  entre 2013 e 2015 em Ubatuba.'

df['2013-09-11':'2015-03-31'].plot(ax=ax,title=titulo)

ax.axvline('2014-01-14 00:00',color='k')
ax.axvline('2014-02-15 00:00',color='k')

ax.axvline('2014-01-01 00:00',color='k',linestyle='--')
ax.axvline('2014-03-01 00:00',color='k',linestyle='--')

# ax.margins(0)

# plotando um resample a cada 15 dias
newdf['2013-09-11':'2015-03-31'].radiacao.resample('1W').mean().plot(ax=ax)

plt.margins(0)
plt.tight_layout()

#plt.savefig('/media/danilo/Danilo/mestrado/github/masterThesis_analysis/figures/dados_observados/IOUSP/radiacao_ubatuba.pdf')
"""
converter as unidades para uma comparação mais precisa com os dados do INMET (kJ m-2)
"""
