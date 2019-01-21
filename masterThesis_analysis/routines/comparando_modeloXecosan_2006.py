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


def make_map_locations(grid,moorings):

    lon = grid[0]
    lat = grid[1]

    fig,ax = plt.subplots()
    m = oceano.make_map(ax,ulon=-44.,llon=-47.,ulat=-23.,llat=-26.,resolution='f')

    # plotting moorings
    for key in moorings.keys():
        ilon = moorings[key][1]
        ilat = moorings[key][0]
        i    = moorings[key][2]
        j    = moorings[key][3]
        m.scatter(ilon,ilat,s=30,marker='o',c='b',latlon=True,label=key)
        m.scatter(lon[i,j],lat[i,j],c='r',marker='^',latlon=True,label='SSB Grid - %s'%(key))

    plt.legend()

def qualityControl_insitu(df):
    limsup = df.UVE.mean() + 3*df.UVE.std()
    liminf = df.UVE.mean() - 3*df.UVE.std()

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR   = oceano.make_dir()
ECOSAN_DIR = BASE_DIR.replace('github','ventopcse/data/ECOSAN') #'/media/danilo/Danilo/mestrado/ventopcse/data/ECOSAN/'
MODELO_DIR = BASE_DIR.replace('github','ventopcse/output') #'/media/danilo/Danilo/mestrado/ventopcse/output/'

magnetic = 24.
angleRot = 55.


# nome dos arquivos para analisar
simulacao  = 'EC2.cdf'
observacao = 'ecosan_ES0802_56m.nc' # Monte Trigo, prof 3m (-23.845, -45.66)

# fundeios localicazao ordenados por:
moorings = {
    'Santos':      [-25.19, -45.82,29,74],
    u'Peruíbe':    [-24.406,-46.89,23,27],
    'Monte Trigo': [-23.845,-45.66,38,11]
}

# ponto de grade do modelo mais proximo ao local do fundeio
i=29
j=74

observado = xr.open_dataset(ECOSAN_DIR+observacao)
# u,v = decomp.intdir2uv(observado.ASPD.values,observado.AVDIR.values,24,55)
u = observado.AVE.values
v = observado.AVN.values

inte,dire = decomp.uv2intdir(u,v,magnetic,angleRot)
ur,vr     = decomp.intdir2uv(inte,dire,magnetic,angleRot)

df_obs = pd.DataFrame({'along':v,'cross':u},index=pd.DatetimeIndex(observado.index.values))
df_obs = df_obs['2006']
# extrair a posição do dado coletado dos atributos gerais do netcdf
# latObs = float(observado.attrs['Latitude'])
# lonObs = float(observado.attrs['Longitude'])

# extrair grade do modelo
timeMod = pd.date_range(start='2006-01-09 01:30:00',end='2006-03-01 22:30:00',freq='3H')

modelado  = xr.open_dataset(MODELO_DIR+simulacao)
lonMod = modelado.lon.values
latMod = modelado.lat.values
lonMod[lonMod == 0.] = np.nan
latMod[latMod == 0.] = np.nan

depth = modelado.depth[i,j].values
sigma = modelado.sigma.values

k = 10 #3 # find_verticalSigma(depth,sigma,-23.)

u = modelado.u[:,k,i,j].values
v = modelado.v[:,k,i,j].values
df_mod = pd.DataFrame({'u':u,'v':v},index=timeMod)
######################################################
# rotacionando os vetores de saída do modelo
#####################################################
inte,dire = decomp.uv2intdir(u,v,magnetic,angleRot)
ur,vr     = decomp.intdir2uv(inte,dire,magnetic,angleRot)

df_mod_rot = pd.DataFrame({'along':vr,'cross':ur},index=timeMod)
#
# # rotacionando e corrigindo declinação magnética dos dados ECOSAN
# ur,vr = decomp.intdir2uv(observado.ASPD.values,observado.AVDIR.values,magnetic,angleRot.)
#
# df_obs_rot = pd.DataFrame({'along':vr,'cross':ur},index=pd.DatetimeIndex(observado.index.values))
df_obs_rot = df_obs
# organizando as series temporais
df_obs_rot = df_obs_rot['2006-01-09 01:30:00':]
df_mod_rot = df_mod_rot[:'2006-02-09 11:00:00']

df_mod_rot *= 100
# df_obs_rot /= 10

re = df_obs_rot.along.resample('3H').mean()
mo = df_mod_rot.along.resample('3H').mean()

mo = mo['2006-01-09 00:00:00':'2006-02-01 00:00:00']
re = re[:'2006-02-01 00:00:00']

skill_along = oceano.skill_willmott(re,mo)

re = df_obs_rot.cross.resample('3H').mean()
mo = df_mod_rot.cross.resample('3H').mean()

mo = mo['2006-01-09 00:00:00':'2006-02-01 00:00:00']
re = re[:'2006-02-01 00:00:00']

skill_cross = oceano.skill_willmott(re,mo)

fig,ax = plt.subplots(nrows=2)

df_obs_rot.along.resample('3H').mean().plot(ax=ax[0])
df_mod_rot.along.resample('3H').mean().plot(ax=ax[0])
ax[0].legend(['','Skill: %1.3f'%(skill_along)])


df_obs_rot.cross.resample('3H').mean().plot(ax=ax[1])
df_mod_rot.cross.resample('3H').mean().plot(ax=ax[1])
ax[1].legend(['','Skill: %1.3f'%(skill_cross)])
