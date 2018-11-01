"""

Realizando analise estatistica feita por Fontes e Castro (2017), com dados
em um correntógrafo instalado nas mediações da Laje de Santos.

LS1 [-24.315, -46.162], local depth 39.3m
    Periodo:
        28/08/2013 10:24h
        04/07/2014 12:34h

Tratamento dos dados:

Correção da declinação magnética
Rotação dos vetores na direção da isóbata local: -45º
Filtragem passa baixa (Lanczos, 40h)


"""

%reset -f

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

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
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

def filtering(d):

    from oceans import lanc

    freq = 1./40  # Hours
    window_size = 96+1+96
    pad = np.zeros(window_size) * np.NaN

    wt = lanc(window_size, freq)

    lowSignal = np.convolve(wt,d,mode='same')

    return lowSignal


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')

# define mooring location
ilat = -24.315
ilon = -46.162

# select which experiment you want to plot:
exp = 'EA1'
if exp == 'EA1':
    typeRun = u'Anômalo - Qnet'
else:
    typeRun = 'Controle - Qnet'

# define file
fname = glob.glob(DATA_DIR+exp+".cdf")[0]
# load output data
ncin = xr.open_dataset(fname)

# select some variables
time=ncin.time
lon = ncin.lon
lat = ncin.lat
# find location in the model_grid
i,j = oceano.find_nearest(lon.values,lat.values,ilon,ilat)
# k = oceano.find_nearest1d

# extract data
dep = np.squeeze(ncin.depth[i,j].values) # local depth
sig = ncin.sigma # sigma level
prof= dep*sig # real depth

# extract current information
u = np.squeeze(ncin.u[:,:,i,j])
v = np.squeeze(ncin.v[:,:,i,j])

# corrigindo declinação magnética de -21.09º para 2014
degMag = np.deg2rad(-21.09)
u,v = rotaciona(u,v,degMag)

# rotacionando os vetores para a orientação da isóbata local (-45.)
alpha = np.deg2rad(-45.)
ur,vr = rotaciona(u,v,alpha)

# creating dataframe for each level
columns = 'level-8m,level-23m,level-35m'.split(',')
levels  = [4, 13, 19]
dct_u = {}
dct_v = {}

for key,k in zip(columns,levels):
    dct_u[key] = ur[:,k]
    dct_v[key] = vr[:,k]

ucomp = pd.DataFrame(dct_u,index=pd.DatetimeIndex(time.values))
vcomp = pd.DataFrame(dct_v,index=pd.DatetimeIndex(time.values))

# filtering tidal signal and high frequencies
# creating dataframe for each level
columns = 'level-8m,level-23m,level-35m'.split(',')
levels  = [4, 13, 19]
dct_u = {}
dct_v = {}

for key,k in zip(columns,levels):
    dct_u[key] = filtering(ur[:,k].values)
    dct_v[key] = filtering(vr[:,k].values)

ufilt = pd.DataFrame(dct_u,index=pd.DatetimeIndex(time.values))
vfilt = pd.DataFrame(dct_v,index=pd.DatetimeIndex(time.values))

# criando dataframe com dados estatisticos, dos dados filtrados (avg,max,min,std)
L = '-8.0,-23.0,-35.0'.split(',')
C = 'Min Max Avg SD'.split(' ')

# creating dataframe for average values
def media(x):
    return x.describe()['mean']
def Maximo(x):
    return x.describe()['max']
def Minimo(x):
    return x.describe()['min']
def STD(x):
    return x.describe()['std']
avg = pd.DataFrame(ufilt.apply(media))
Max = pd.DataFrame(ufilt.apply(Maximo))
Min = pd.DataFrame(ufilt.apply(Minimo))
Std = pd.DataFrame(ufilt.apply(STD))

stats_cross = pd.DataFrame({
    'avg':avg[0],
    'Max':Max[0],
    'Min':Min[0],
    'SD':Std[0]})

avg = pd.DataFrame(vfilt.apply(media))
Max = pd.DataFrame(vfilt.apply(Maximo))
Min = pd.DataFrame(vfilt.apply(Minimo))
Std = pd.DataFrame(vfilt.apply(STD))

stats_along = pd.DataFrame({
    'avg':avg[0],
    'Max':Max[0],
    'Min':Min[0],
    'SD':Std[0]})

# visualizando tudo
fig,axes = plt.subplots(nrows=2,ncols=3,sharex=True,figsize=(20,5))

## 8m depth
ucomp['level-8m'].plot(ax=axes[0,0],color='#c0c0c0')
ufilt['level-8m'].plot(ax=axes[0,0],color='k')

vcomp['level-8m'].plot(ax=axes[1,0],color='#c0c0c0')
vfilt['level-8m'].plot(ax=axes[1,0],color='k')

axes[0,0].set_title('8m depth')
axes[0,0].set_ylabel('cross')
axes[1,0].set_ylabel('along')

## 23m depth
ucomp['level-23m'].plot(ax=axes[0,1],color='#c0c0c0')
ufilt['level-23m'].plot(ax=axes[0,1],color='k')

vcomp['level-23m'].plot(ax=axes[1,1],color='#c0c0c0')
vfilt['level-23m'].plot(ax=axes[1,1],color='k')

axes[0,1].set_title('23m depth')
## 35m depth
ucomp['level-35m'].plot(ax=axes[0,2],color='#c0c0c0')
ufilt['level-35m'].plot(ax=axes[0,2],color='k')

vcomp['level-35m'].plot(ax=axes[1,2],color='#c0c0c0')
vfilt['level-35m'].plot(ax=axes[1,2],color='k')
axes[0,2].set_title('35m depth')

plt.suptitle(u'Experimento %s'%typeRun,fontsize=24)
plt.subplots_adjust(top=0.86,bottom=0.18,left=0.125,right=0.9,hspace=0.2,wspace=0.2)
