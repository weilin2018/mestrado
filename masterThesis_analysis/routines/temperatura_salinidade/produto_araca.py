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

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
def find_timestep(ncin,timestep):

	refDate = parser.parse(timestep)

	time = ncin.time.values
	df   = pd.DataFrame(np.arange(0,len(time)),index=time)

	d = df.iloc[df.index.get_loc(refDate,method='nearest')]

	return d

def extrair_dados(ncin,i,j,t):

	temp = ncin.temp[t.values,:,i,j]
	salt = ncin.salt[t.values,:,i,j]

	return temp,salt


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

######################## PERFIL VERTICAL ########################################

ponto = oceano.procurar_pontos_grade(-45.40085130,-23.81778822,lon,lat,n=1)
i     = int(ponto[0][0])
j     = int(ponto[0][1])

# localizando instantes de tempos correspondentes aos dados obtidos pelo CEBIMAR
instantesCEBIMAR = ['2014-01-14 13:18','2014-01-14 17:15']

d = find_timestep(ncin,'2014-01-16 16:13')

# extraindo dados termohalinos
temp,salt    = extrair_dados(ncin,i,j,d)
profundidade = dep[i,j]*sig

fig,ax = plt.subplots(ncols=2,figsize=(6,8),sharey=True)

ax[0].plot(np.squeeze(temp),profundidade)
ax[0].set_xlim([23,30])
ax[0].xaxis.tick_top()

ax[1].plot(np.squeeze(salt),profundidade)
ax[1].set_xlim([34,36])
ax[1].xaxis.tick_top()

plt.show()

###################### SERIE TEMPORAL ##########################################

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

fig,ax = plt.subplots(nrows=2)
#ax[0].plot(salt.Surface,'#c0c0c0',label='Surface Salinity')
ax[0].plot(salt.Surface,'#c0c0c0',label='Surface Salinity')
ax[0].plot(salt.Bottom,'k',label='Bottom Salinity')

ax[0].legend(loc='best')
ax[0].set_ylabel('Salinity',fontsize=8)
ax[0].set_ylim([30.5, 35.5])
ax[0].set_xlabel('Time',fontsize=8)

ax[1].plot(temp.Surface,'#c0c0c0',label='Surface Temperature')
ax[1].plot(temp.Bottom,'#000000',label='Bottom Temperature')

ax[1].legend(loc='best')
ax[1].set_ylabel('Temperature',fontsize=8)

ax[1].set_xlabel('Time',fontsize=8)

plt.show()
