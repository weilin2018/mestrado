"""
IMPORTANTE:

Para se rodar essa rotina, deve-se antes rodar a rotina produto_Tratamento.py,
para que os arquivos .nc com os dados do modelo tratado seja criado.

"""


import glob
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
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
import scipy.io as sio

import matplotlib
matplotlib.style.use('ggplot')
matplotlib.use('PS')

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
DATA_DIR = BASE_DIR.replace('github','ventopcse/data/IOUSP')
SAVE_DIR = BASE_DIR.replace('github','ventopcse')
FIG_DIR  = BASE_DIR + 'masterThesis_analysis/figures/validacao/'

##################### IMPORTING OBSERVATIONAL DATA
# reading file with CSS filtered data
dfFilt = xr.open_dataset(SAVE_DIR + "CSS_filtered.nc")
dfFilt = dfFilt.to_dataframe()
ts = []
for t in dfFilt.index.values:
    ts.append(t[0])
dfFilt.index = ts
# new datetimeIndex
dtRange = pd.date_range(start='2013-12-02 21:00',end='2014-03-06 19:00',freq='30Min')
dfFilt.index = dtRange

# cutting and resampling insitu data
dfTmp = dfFilt['2014-01-09 01:30:00':'2014-03-01 22:30:00'].copy()
dfObse = dfTmp.resample('3h').mean()
del dfTmp

##################### IMPORTING MODELED DATA
# importing product data, alread rotated and saved in a netcdf file
df5m = xr.open_dataset(SAVE_DIR + 'df5m_EA2.nc')
df5m = df5m.to_dataframe()
df15m = xr.open_dataset(SAVE_DIR + 'df15m_EA2.nc')
df15m = df15m.to_dataframe()
dfTemp = xr.open_dataset(SAVE_DIR + 'dfTemp_EA2.nc')
dfTemp = dfTemp.to_dataframe()

# putting all data in the same dataframe
dct = {
    'u5m': df5m['along 5m'].values,
    'v5m': df5m['cross 5m'].values,
    'u15m': df15m['along 15m'].values,
    'v15m': df15m['cross 15m'].values,
    'T5m': dfTemp['T5'].values,
    'T15m': dfTemp['T15'].values
}

dfModel = pd.DataFrame(dct,index=pd.date_range(start='2014-01-09 01:30:00',end='2014-03-01 22:30:00',freq='3h'))

# calculating skill coefficient for current
skillu5m  = oceano.skill_willmott(dfObse.u5m.values,dfModel.u5m.values)
skillv5m  = oceano.skill_willmott(dfObse.v5m.values,dfModel.v5m.values)
skillu15m = oceano.skill_willmott(dfObse.u15m.values,dfModel.u15m.values)
skillv15m = oceano.skill_willmott(dfObse.v15m.values,dfModel.v15m.values)
# correlation for temperature because is a better predictor of the behaviour
T5m_Corr = np.corrcoef(dfObse.T5m.values,dfModel.T5m.values)[0,1]
T15m_Corr = np.corrcoef(dfObse.T15m.values,dfModel.T15m.values)[0,1]


fig,axes = plt.subplots(nrows=6,sharex=True,figsize=(12./2.54,15/2.54))

axes[0].set_title('Componente Perpendicular - Prof. de 5m (Skill %0.2f)'%(skillu5m),fontsize=8)
axes[0].plot(dfObse.index.values,dfObse.u5m.values)
axes[0].plot(dfModel.index.values,dfModel.u5m.values)
axes[0].set_ylabel('u ['+r'$m^2$'+' s]',fontsize=8,labelpad=-2)
axes[0].legend(['Observado','Modelado'],fontsize=6,bbox_to_anchor=(0.,-0.1,.4,.5),mode='expand',ncol=2)

axes[1].set_title('Componente Paralela - Prof. de 5m (Skill %0.2f)'%(skillv5m),fontsize=8)
axes[1].plot(dfObse.index.values,dfObse.v5m.values)
axes[1].plot(dfModel.index.values,dfModel.v5m.values)
axes[1].set_ylabel('v ['+r'$m^2$'+' s]',fontsize=8,labelpad=-2)
axes[1].legend(['Observado','Modelado'],fontsize=6,bbox_to_anchor=(0.,-0.1,.4,.5),mode='expand',ncol=2)

axes[2].set_title('Componente Perpendicular - Prof. de 15m (Skill %0.2f)'%(skillu15m),fontsize=8)
axes[2].plot(dfObse.index.values,dfObse.u15m.values)
axes[2].plot(dfModel.index.values,dfModel.u15m.values)
axes[2].set_ylabel('u ['+r'$m^2$'+' s]',fontsize=8,labelpad=-2)
axes[2].legend(['Observado','Modelado'],fontsize=6,bbox_to_anchor=(0.,-0.1,.4,.5),mode='expand',ncol=2)

axes[3].set_title('Componente Paralela - Prof. de 15m (Skill %0.2f)'%(skillv15m),fontsize=8)
axes[3].plot(dfObse.index.values,dfObse.v15m.values)
axes[3].plot(dfModel.index.values,dfModel.v15m.values)
axes[3].set_ylabel('v ['+r'$m^2$'+' s]',fontsize=8,labelpad=-2)
axes[3].legend(['Observado','Modelado'],fontsize=6,bbox_to_anchor=(0.,-0.1,.4,.5),mode='expand',ncol=2)

# configuring figure
axes[0].set_ylim([-1.0,1.0])
axes[1].set_ylim([-1.0,1.0])
axes[2].set_ylim([-.5,.5])
axes[3].set_ylim([-.5,.5])

axes[4].set_title('Temperatura - Prof. de 5m (Corr. %0.2f)'%(T5m_Corr),fontsize=8)
axes[4].plot(dfObse.index.values,dfObse.T5m.values)
axes[4].plot(dfModel.index.values,dfModel.T5m.values)
axes[4].set_ylabel('T ['+r'$^o$'+'C]',fontsize=8,labelpad=-1)
axes[4].legend(['Observado','Modelado'],fontsize=6,bbox_to_anchor=(0.,-0.1,.4,.5),mode='expand',ncol=2)

axes[5].set_title('Temperatura - Prof. de 15m (Corr. %0.2f)'%(T15m_Corr),fontsize=8)
axes[5].plot(dfObse.index.values,dfObse.T15m.values)
axes[5].plot(dfModel.index.values,dfModel.T15m.values)
axes[5].set_ylabel('T ['+r'$^o$'+'C]',fontsize=8,labelpad=-1)
axes[5].legend(['Observado','Modelado'],fontsize=6,bbox_to_anchor=(0.,-0.1,.4,.5),mode='expand',ncol=2)

for ax in axes:
    ax.axes.tick_params(axis='both',which='both',labelsize=8)
    ax.margins(xmargin=0)

# automatic formatting, for xdate
fig.autofmt_xdate()

axes[-1].xaxis.set_major_locator(mdates.WeekdayLocator())
axes[-1].xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))

plt.suptitle(u'Parâmetros Estatísticos para EA2',fontsize=10)

plt.tight_layout()
# plt.subplots_adjust(top=0.945,bottom=0.11,left=0.126,right=0.976,hspace=0.587,wspace=0.2)
plt.subplots_adjust(top=0.91,bottom=0.06,left=0.126,right=0.976,hspace=0.587,wspace=0.2)

plt.savefig(FIG_DIR + 'EA2xCSS_validacao.pdf')
