# calculando momentos estatisticos para comparar com as tabelas 2 e 3 de Fontes e Castro (2017)

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


from oceans import lanc
import decomp as dp
##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def rotate(df,angrot,key,mag=-24.):
    INT,DIR = dp.uv2intdir(df[key].values,df[key.replace('u','v')].values,0,0)
    ur,vr    = dp.intdir2uv(INT,DIR,mag,angrot)

    return ur,vr

def print_estatistica(df,loc,prof,alpha):

    meansCross,maxsCross,minsCross,SDsCross = [],[],[],[]
    meansAlong,maxsAlong,minsAlong,SDsAlong = [],[],[],[]

    for p in prof:
        ur,vr = rotate(df,alpha,'u '+p)
        meansCross.append(np.nanmean(ur))
        maxsCross.append(np.nanmax(ur))
        minsCross.append(np.nanmin(ur))
        SDsCross.append(np.nanstd(ur))
        meansAlong.append(np.nanmean(vr))
        maxsAlong.append(np.nanmax(vr))
        minsAlong.append(np.nanmin(vr))
        SDsAlong.append(np.nanstd(vr))

    dct_cross = {
        'Min':minsCross,
        'Max':maxsCross,
        'Avg':meansCross,
        'SD':SDsCross,
    }

    cross = pd.DataFrame(dct_cross,index=prof)

    dct_along = {
        'Min':minsAlong,
        'Max':maxsAlong,
        'Avg':meansAlong,
        'SD':SDsAlong,
    }

    along = pd.DataFrame(dct_along,index=prof)

    return cross,along
##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/output')

# definicoes gerais
LS1       = [-24.3154,-46.163] # coordenada ponto LS1
LS1_alpha = np.deg2rad(-45.)   # angulo para rotacao dos vetores em LS1
LS2       = [-24.2865,-46.195] # coordenada ponto LS2
LS2_alpha = np.deg2rad(90.)    # angulo para rotacao dos vetores em LS2

ncin = xr.open_dataset(DATA_DIR + 'EA1.cdf')

# localizando os indices para cada estacao
lon,lat = ncin.lon.values,ncin.lat.values
indLS1  = oceano.procurar_pontos_grade(LS1[1],LS1[0],lon,lat)[0]
indLS2  = oceano.procurar_pontos_grade(LS2[1],LS2[0],lon,lat)[0]

# localizando a coordenada sigma
sigLS1 = [4,12,20]  # -5m, -14m, -26m
sigLS2 = [1,10,19] # -3m,-19m,-35m

# importando os dados

# LS1
dct  = {
    'u 8m': ncin.u[:,sigLS1[0],int(indLS1[0]),int(indLS1[1])].values,
    'u 23m': ncin.u[:,sigLS1[1],int(indLS1[0]),int(indLS1[1])].values,
    'u 35m': ncin.u[:,sigLS1[2],int(indLS1[0]),int(indLS1[1])].values,
    'v 8m': ncin.v[:,sigLS1[0],int(indLS1[0]),int(indLS1[1])].values,
    'v 23m': ncin.v[:,sigLS1[1],int(indLS1[0]),int(indLS1[1])].values,
    'v 35m': ncin.v[:,sigLS1[2],int(indLS1[0]),int(indLS1[1])].values,
}

# convertendo para dataframe
LS1_data = pd.DataFrame(dct,index=ncin.time.values)

# LS2
dct  = {
    'u 3m': ncin.u[:,sigLS2[0],int(indLS2[0]),int(indLS2[1])].values,
    'u 19m': ncin.u[:,sigLS2[1],int(indLS2[0]),int(indLS2[1])].values,
    'u 35m': ncin.u[:,sigLS2[2],int(indLS2[0]),int(indLS2[1])].values,
    'v 3m': ncin.v[:,sigLS2[0],int(indLS2[0]),int(indLS2[1])].values,
    'v 19m': ncin.v[:,sigLS2[1],int(indLS2[0]),int(indLS2[1])].values,
    'v 35m': ncin.v[:,sigLS2[2],int(indLS2[0]),int(indLS2[1])].values,
}

# convertendo para dataframe
LS2_data = pd.DataFrame(dct,index=ncin.time.values)


### tratamento das series temporais

# obtendo sinal sub inercial, com algo semelhante ao lanczos (33h)
dfLS1 = pd.rolling_mean(LS1_data,window=33,center=True,freq='3H')
dfLS2 = pd.rolling_mean(LS2_data,window=33,center=True,freq='3H')
# rotacionando o sistema de coordenadas para ficar alinhado as isobatas

stats_LS1_perpendicular,stats_LS1_paralela = print_estatistica(dfLS1,'LS1',['8m','23m','35m'],LS1_alpha)
stats_LS2_perpendicular,stats_LS2_paralela = print_estatistica(dfLS2,'LS2',['3m','19m','35m'],LS2_alpha)
