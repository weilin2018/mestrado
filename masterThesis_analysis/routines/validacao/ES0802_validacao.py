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
import decomp
import cmocean as cmo

# from pandas import rolling_mean

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def find_verticalSigma(depth,sigma,ref):

    newDep = depth*sigma

    dep = oceano.find_nearest_1D(newDep,ref)

    return np.where(newDep == dep)

def quality_control(df):

    # removing values outside 3 standard deviations
    df[df.along > (df.along.mean() + 3*df.along.std())] = np.nan
    df[df.along < (df.along.mean() - 3*df.along.std())] = np.nan

    df[df.cross > (df.cross.mean() + 3*df.cross.std())] = np.nan
    df[df.cross < (df.cross.mean() - 3*df.cross.std())] = np.nan

    return df

def set_informations(mooring):

    if mooring == 0:
        observacao = 'ecosan_ES0802_56m.nc'
        i,j = 29,74
        k = 10
    elif mooring == 1:
        observacao = 'ecosan_ES0902_12m.nc'
        i,j = 23,27
        k = 13
    elif mooring == 2:
        observacao = 'ecosan_ES1002_12m.nc'
        i,j = 38,11
        k = 9

    return observacao,i,j,k

def rotaciona(u,v,angulo):
    """
    ATUALIZADA!
    Precisa de um angulo alpha definido
    Based on book Data Analysis in Physical Oceanography, page 425
    funcao traduzida por: Paula Birocchi
    Parameters:
    u = componente leste-oeste da corrente
    v = componente norte-sul da corrente
    angulo = angulo em radianos da costa em relacao a horizontal (direcao leste-oeste)
    Output:
    urotacionado = velocidade cross-shore (componente perpendicular a costa)
    vrotacionado = velocidade along-shore (componente paralela a costa)
    Observacao: o angulo eh sempre o mesmo neste caso. Para o Canal de
    Sao Sebastiao,por exemplo, nao eh o ideal, pois o canal eh mega curvilinear.
    Nestes casos, eh melhor utilizar rotacionando_vento_corrente_con_ang_saida.py
    """
    # Usado se a referencia for no LESTE com o ANGULO MEDIDO no sentido ANTI-HORARIO. Como está no Data Analysis:
    #   urotacionado = u * np.cos(angulo) + v * np.sin(angulo) este era o meu original# de acordo com Dottori tbm.
    #   vrotacionado = -u * np.sin(angulo) + v*np.cos(angulo) # este era o original! vamos ver agora. # de acordo com Dottori tbm.
    # Usado se a referencia for no NORTE com o ANGULO MEDIDO no sentido HORARIO. Canal de Sao Sebastiao: +51 graus.
    vrotacionado = u * np.sin(angulo) + v*np.cos(angulo) # de acordo com Ze Roberto
    urotacionado = u*np.cos(angulo) - v*np.sin(angulo) # de acordo com Ze Roberto
    return urotacionado, vrotacionado

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

#### DEFINING SOME VARIABLES
BASE_DIR   = oceano.make_dir()
ECOSAN_DIR = BASE_DIR.replace('github','ventopcse/data/ECOSAN')
MODELO_DIR = BASE_DIR.replace('github','ventopcse/output')
# ECOSAN_DIR = '/media/danilo/Danilo/mestrado/ventopcse/data/ECOSAN/'
# MODELO_DIR = '/media/danilo/Danilo/mestrado/ventopcse/output/'

magnetic = 24.
angleRot = 55.

# selecting which mooring to plot
os.system('clear')
mooring = input('Digite qual estacao plotar [0 - Santos, 1 - Peruibe e 2 - Monte Trigo]: ')
observacao,i,j,k = set_informations(mooring)

# nome dos arquivos para analisar
simulacao  = 'EC1.cdf'
# observacao = 'ecosan_ES0802_56m.nc' # Monte Trigo, prof 3m (-23.845, -45.66)

# fundeios localicazao ordenados por:
moorings = {
    'Santos':      [-25.19, -45.82,29,74],
    u'Peruíbe':    [-24.406,-46.89,23,27],
    'Monte Trigo': [-23.845,-45.66,38,11]
}

titulos = {
    0: 'Santos [56m]',
    1: 'Peruibe [11.4m]',
    2: 'Monte Trigo [12m]'
}

# ponto de grade do modelo mais proximo ao local do fundeio
# i=29
# j=74

#### LOADING OBSERVATIONAL DATA
observado = xr.open_dataset(ECOSAN_DIR+observacao)
u = observado.AVE.values
v = observado.AVN.values
alpha = np.deg2rad(magnetic+angleRot)
ur,vr = rotaciona(u,v,alpha) # rotating vectors
# creating dataframe
df_obs = pd.DataFrame({'along':vr,'cross':ur,'u':u,'v':v},index=pd.DatetimeIndex(observado.index.values))
# df_obs = df_obs['2006']
# removing some spurious values based on a basic criteria
df_obs = quality_control(df_obs)
df_obs = df_obs.interpolate(method='pchip')

#### LOADING MODELLED DATA
modelado  = xr.open_dataset(MODELO_DIR+simulacao)
lonMod = modelado.lon.values
latMod = modelado.lat.values
lonMod[lonMod == 0.] = np.nan
latMod[latMod == 0.] = np.nan
depth = modelado.depth[i,j].values
sigma = modelado.sigma.values

# k = 3 #10 # find_verticalSigma(depth,sigma,-23.)

u = modelado.u[:,k,i,j].values
v = modelado.v[:,k,i,j].values
alpha = np.deg2rad(0+angleRot)
ur,vr = rotaciona(u,v,alpha) # rotating vectors
df_mod = pd.DataFrame({'along':vr,'cross':ur,'u':u,'v':v},index=pd.date_range(start='2006-01-09 01:30',end='2006-03-01 22:30',freq='3H'))

# filtering using rolling_mean, from pandas
# df_filt_obs = df_obs.copy()
# df_filt_obs = rolling_mean(df_filt_obs, window=40, center=True,freq='3H')
df_filt_obs = df_obs.resample('3H').mean()

# df_filt_mod = df_mod.copy()
# df_filt_mod = rolling_mean(df_filt_mod, window=40, center=True,freq='3H')
df_filt_mod = df_mod.resample('3H').mean()


# cutting, by time, to match both dataframes
df_obs = df_filt_obs[df_filt_mod.index[0]:]
df_mod = df_filt_mod[:df_filt_obs.index[-1]]


# converting observed data to m/s. Apparently, the original is cm/s
df_obs /= 100


# filtering using rolling_mean, from pandas
df_filt_obs = df_obs.copy()
# df_filt_obs = rolling_mean(df_filt_obs, window=40, center=True,freq='3H')

df_filt_mod = df_mod.copy()
# df_filt_mod = rolling_mean(df_filt_mod, window=40, center=True,freq='3H')

# df_filt = df_obs.rolling(20,win_type='hamming').mean()

# df_obs = df_obs - df_filt
# df_filt_obs.cross *= -1

skill_along = oceano.skill_willmott(df_filt_obs.along.values,df_filt_mod.along.values)
skill_cross = oceano.skill_willmott(df_filt_obs.cross.values,df_filt_mod.cross.values)



fig,ax = plt.subplots(nrows=2,figsize=(15,8))

df_obs.along.plot(ax=ax[0],label='ECOSAN')
df_mod.along.plot(ax=ax[0],label='Modelo %s'%(simulacao.split('.')[0]))
# ax[0].legend(['','Skill: %1.3f'%(skill_along)])
ax[0].set_title('Componente Paralela - Skill: %1.3f'%(skill_along))
ax[0].legend()

df_obs.cross.plot(ax=ax[1],label='ECOSAN')
df_mod.cross.plot(ax=ax[1],label='Modelo %s'%(simulacao.split('.')[0]))
# ax[1].legend(['','Skill: %1.3f'%(skill_cross)])
ax[1].set_title('Componente Perpendicular - Skill: %1.3f'%(skill_cross))
ax[1].legend()

plt.suptitle(titulos[mooring], fontsize=25)
