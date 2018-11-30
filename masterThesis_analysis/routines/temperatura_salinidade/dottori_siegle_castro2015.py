#-*-coding;utf-8-*-
"""

"""
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
import seawater as sw

# pacotes para minimap
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

def calcTrend(df):

	coefficients, residuals, _, _, _ = np.polyfit(range(len(df.index)),df.values,1,full=True)
	mse = residuals[0]/(len(df.index))
	nrmse = np.sqrt(mse)/(df.max() - df.min())
	print('Slope ' + str(coefficients[0]))
	print('NRMSE: ' + str(nrmse))

	return coefficients


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/output')

experiment = 'EA2.cdf'
ncin = xr.open_dataset(DATA_DIR + experiment)


# localizando ponto de grade mais proximo do Araca
lon = ncin.lon.values
lat = ncin.lat.values
lon[lon == 0.] = np.nan
lat[lat == 0.] = np.nan

ponto = oceano.procurar_pontos_grade(-45.403225,-23.816398,lon,lat,n=1)
i     = int(ponto[0][0])
j     = int(ponto[0][1])

teste = raw_input('Deseja plotar para visualizar a localizacao do ponto? [y/N]: ')

if teste == 'y':
	fig,ax = plt.subplots()
	m = oceano.make_map(ax,resolution='f')
	x,y = m(lon,lat)

	m.plot(x,y,'k',alpha=.3)
	m.plot(x.T,y.T,'k',alpha=.3)
	m.scatter(x[i,j],y[i,j],s=30,c='r')

	plt.show()

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

# subamostrando para dados diários
tempDaily = temp.resample('D').mean()
saltDaily = salt.resample('D').mean()

fig,ax = plt.subplots(nrows=2)
#ax[0].plot(salt.Surface,'#c0c0c0',label='Surface Salinity')
ax[0].plot(saltDaily.Surface,'#c0c0c0',label='Surface Salinity')
ax[0].plot(saltDaily.Bottom,'k',label='Bottom Salinity')

ax[0].legend(loc='best')
ax[0].set_ylabel('Salinity',fontsize=8)
ax[0].set_ylim([30.5, 35.5])
ax[0].set_xlabel('Time',fontsize=8)

ax[1].plot(tempDaily.Surface,'#c0c0c0',label='Surface Temperature')
ax[1].plot(tempDaily.Bottom,'#000000',label='Bottom Temperature')

ax[1].legend(loc='best')
ax[1].set_ylabel('Temperature',fontsize=8)

ax[1].set_xlabel('Time',fontsize=8)

plt.show()

# removendo a tendencia das series temporais

# calculating trend for each timeseries
tSurf_coefs = calcTrend(tempDaily.Surface)
tBott_coefs = calcTrend(tempDaily.Bottom)

sSurf_coefs = calcTrend(saltDaily.Surface)
sBott_coefs = calcTrend(saltDaily.Bottom)

# como a tendencia é linear, podemos remove-la da serie usando uma diferenciacao:
# based on: https://machinelearningmastery.com/time-series-trends-in-python/
diff = list()
X = tempDaily.Surface.values
for i in range(1,len(X)):
	value = X[i] - X[i-1]
	diff.append(value)

