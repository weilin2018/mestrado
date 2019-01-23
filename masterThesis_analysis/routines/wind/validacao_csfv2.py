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
import decomp as dp
##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def convert_xarray2dataframe(df):

    df = df.to_dataframe()

    ts = []
    for t in df.index:
        ts.append(t[0])

    df.index = ts

    return df

def cfsv2_data(ncin):
    wu   = np.squeeze(ncin.U_GRD_L103.values)
    wv   = np.squeeze(ncin.V_GRD_L103.values)
    time = ncin.time.values

    df = pd.DataFrame({'wu':wu,'wv':wv},index=pd.DatetimeIndex(time))

    return df

def tratamento_cfsv2(df):

    # rotacao de 36 graus (maxima variancia)
    int,dir = dp.uv2intdir(df.wu.values,df.wv.values,0,0)
    wur,wvr = dp.intdir2uv(int,dir,0,36)

    df = pd.DataFrame({'wur':wur,'wvr':wvr},index=df.index.values)

    # controle de qualidade
    df[df > 3*df.std()] = np.nan
    df[df <-3*df.std()] = np.nan

    return df

# def getStatisticalAnalysis(df1,df2):
#
#     skill_cross = oceano.skill_willmott(df1.wind_cross.values,df2.wur.values)
#     skill_along = oceano.skill_willmott(df1.wind_along.values,df2.wvr.values)
#
#     return False



##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/data/')

# reading netcdf files already treated
lajeObs = xr.open_dataset(DATA_DIR+'Est_lajeSantos/lajesantos.nc')
lajeObs = convert_xarray2dataframe(lajeObs)

boiaObs = xr.open_dataset(DATA_DIR+'pnboia/pnboiaSantos.nc')
boiaObs = convert_xarray2dataframe(boiaObs)

lajeMod = xr.open_dataset(DATA_DIR+'serie_cfsv2/LajeDeSantos/2015/cfsv2_lajedesantos_2015.nc')
lajeMod = cfsv2_data(lajeMod)

boiaMod = xr.open_dataset(DATA_DIR+'serie_cfsv2/pnboia_santos/cfsv2_pnboiasantos_20112017.nc')
boiaMod = cfsv2_data(boiaMod)

# tratamento cfsv2
lajemod = tratamento_cfsv2(lajeMod)
boiamod = tratamento_cfsv2(boiaMod)

# tratamento insitu
lajeObs[lajeObs > 3*lajeObs.std()] = np.nan
lajeObs[lajeObs <-3*lajeObs.std()] = np.nan

# redimensionar para frequencia de 6 horas
lajeobs = lajeObs.resample('6H').mean()
boiaobs = boiaObs.resample('6H').mean()

# selecionando periodo
lajeobs = lajeobs['2015-03-20':]
lajemod = lajemod['2015-03-20':]
boiaobs = boiaobs['2014-11-17':'2015-03']
boiamod = boiamod['2014-11-17':'2015-03']

# plotando LAJE DE SANTOS
fig,ax = plt.subplots(nrows=2)

lajeobs.wind_along.plot(ax=ax[0],label='EMLS/along')
lajemod.wvr.plot(ax=ax[0],label='CFSv2/along')
ax[0].legend()

lajeobs.wind_cross.plot(ax=ax[1],label='EMLS/cross')
lajemod.wur.plot(ax=ax[1],label='CFSv2/cross')

# calculando parametros estatistico
skill_along,corr_along = oceano.getStatisticalAnalysis(lajeobs.wind_along,lajemod.wvr)
skill_cross,corr_cross = oceano.getStatisticalAnalysis(lajeobs.wind_cross,lajemod.wur)

os.system('clear')
print('# ------ ALONG SHORE ------ #')
print('Skill: %0.2f     Corr %0.2f' % (skill_along,corr_along))
print('# ------ CROSS SHORE ------ #')
print('Skill: %0.2f     Corr %0.2f' % (skill_cross,corr_cross))


### PLOTANDO PNBOIA
fig,ax = plt.subplots(nrows=2)

boiaobs.wind_along.plot(ax=ax[0],label='PNBOIA/along')
boiamod.wvr.plot(ax=ax[0],label='CFSv2/along')
ax[0].legend()

boiaobs.wind_cross.plot(ax=ax[1],label='PNBOIA/cross')
boiamod.wur.plot(ax=ax[1],label='CFSv2/cross')

# calculando parametros estatistico
skill_along,corr_along = oceano.getStatisticalAnalysis(boiaobs.wind_along,boiamod.wvr)
skill_cross,corr_cross = oceano.getStatisticalAnalysis(boiaobs.wind_cross,boiamod.wur)

# os.system('clear')
print('# ------ ALONG SHORE ------ #')
print('Skill: %0.2f     Corr %0.2f' % (skill_along,corr_along))
print('# ------ CROSS SHORE ------ #')
print('Skill: %0.2f     Corr %0.2f' % (skill_cross,corr_cross))
