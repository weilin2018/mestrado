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

def load_model(modelado,i,j,k):

    #### LOADING MODELLED DATA
    # modelado  = xr.open_dataset(MODELO_DIR+simulacao)
    lonMod = modelado.lon.values
    latMod = modelado.lat.values
    lonMod[lonMod == 0.] = np.nan
    latMod[latMod == 0.] = np.nan
    depth = modelado.depth[i,j].values
    sigma = modelado.sigma.values

    u = modelado.u[:,k,i,j].values
    v = modelado.v[:,k,i,j].values
    # alpha = np.deg2rad(55)
    # ur,vr = rotaciona(u,v,alpha) # rotating vectors
    # df_mod = pd.DataFrame({'vr':vr,'ur':ur,'u':u,'v':v},index=pd.DatetimeIndex(modelado.time.values))
    df_mod = pd.DataFrame({'u':u,'v':v},index=pd.DatetimeIndex(modelado.time.values))

    return df_mod

def check_rotation(df):

    import mat4py

    ped = pd.DataFrame(mat4py.loadmat('/media/danilo/Danilo/mestrado/ventopcse/data/ECOSAN/morais2016/exp_sECOM/u_v_ecosan_100.mat'))

    # set index as datetime
    ped.index = pd.date_range(start='2006-01-06',end='2006-02-05',freq='6H')[:-1]

    # converting
    ped /= 100

    # cut and plot, based on df index
    df  = df[ped.index[0]:ped.index[-1]]

    fig,ax = plt.subplots(nrows=2,sharex=True)

    ped.u11.plot(ax=ax[0],label='matlab')
    df.cross.plot(ax=ax[0],label='python')

    ax[0].set_title('Skill: %1.2f'%(oceano.skill_willmott(ped.u11.values,df.cross.values)))

    ped.v11.plot(ax=ax[1],label='matlab')
    df.along.plot(ax=ax[1],label='python')

    ax[1].set_title('Skill: %1.2f'%(oceano.skill_willmott(ped.v11.values,df.along.values)))

    ax[1].legend()

    plt.suptitle('Verificacao da rotacao entre matlab e python')

    plt.show()


def filtering(d):

    from oceans import lanc

    freq = 1./40  # Hours
    window_size = 96+1+96
    pad = np.zeros(window_size) * np.NaN

    wt = lanc(window_size, freq)

    lowSignal = np.convolve(wt,d,mode='same')

    return lowSignal

def interpolating_model(df):

    # frequency of 1.5hour, starting from the first time given as output
    dt = pd.date_range(start='2006-01-09 01:30:00',end='2006-03-01 22:30:00',freq='1.5H')

    dct = {
        'u':      np.zeros(dt.shape[0])*np.nan,
        'v':      np.zeros(dt.shape[0])*np.nan,
        'along':  np.zeros(dt.shape[0])*np.nan,
        'cross':  np.zeros(dt.shape[0])*np.nan
    }

    newdf = pd.DataFrame(dct,index=pd.DatetimeIndex(dt))

    # merging dataframes
    final = pd.concat([df,newdf],axis=1)
    # removing columns with nan values
    final.dropna(axis=1,how='all',inplace=True)
    # interpolating missing values
    final.interpolate(inplace=True)

    return final

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

DATA_DIR = '/media/danilo/Danilo/mestrado/ventopcse/data/ECOSAN/Dados_mexidos/ES0802.Dat'
MODELO_DIR = '/media/danilo/Danilo/mestrado/ventopcse/output/'
simulacao  = 'EC1.cdf'

# set observed alpha, for rotating: alpha = declinacao mag + rotacao
alpha = np.deg2rad(24. + 55.)

##############################################
#       LEITURA DOS DADOS ES0802.Dat         #
##############################################
data = pd.read_csv(DATA_DIR,delimiter=',',usecols=[0,1,2,3,5,12],header=None)
data.columns = ['AVN','AVE','ASPD','AVDIR','datetime','TEMP']

data.index = pd.DatetimeIndex(data.datetime.values) # assigning datetimeindex
data.drop('datetime',axis=1,inplace=True) # removing column datetime

data.drop(['TEMP','ASPD','AVDIR'],axis=1,inplace=True) # removing temperature to validate only

# renaming columns
data.rename(columns={"AVN":'v',"AVE":'u'},inplace=True)

##############################################
#       REAMOSTRAGEM PARA FREQ DE 1H         #
##############################################
data = data.resample('1H').mean()

##############################################
#       CONVERSÃO DE cm/s PARA m/s           #
##############################################
# convertendo de cm/s para m/s
data /= 100

##############################################
#       ROTAÇÃO DOS VETORES u E v            #
##############################################
# rotating vectors, considering also magnetic declination
data['cross'],data['along'] = rotaciona(data.u.values,data.v.values,alpha)

# checking rotation with matlab's output of Morais (2016)
check_rotation(data.resample('6H').mean())

##############################################
#   FILTRANDO A MARÉ DO DADOS OBSERVADOS     #
##############################################
"""
    Como o modelo foi forçado SEM a maré, precisamos
    remover o sinal de maré astronômica dos dados observados,
    para poder fazer uma comparação justa entre as informações.
"""
# we need to remove tidal signal, by filtering
along_filt = filtering(data.along.values)
cross_filt = filtering(data.cross.values)
df = pd.DataFrame({'along_filt':along_filt,'cross_filt':cross_filt},index=data.index)

# plotting original vs filtered
fig,ax = plt.subplots(nrows=2)
ax[0].plot(data['2006'].index,data['2006'].cross,'#c0c0c0',label='original')
ax[0].plot(df['2006'].index,df['2006'].cross_filt,'k',label='filtered')
ax[0].set_title('Along shore component')
ax[0].legend()
ax[1].plot(data['2006'].index,data['2006'].along,'#c0c0c0',label='original')
ax[1].plot(df['2006'].index,df['2006'].along_filt,'k',label='filtered')
ax[1].set_title('Cross Shore Component')
ax[1].legend()

##############################################
#       REAMOSTRAGEM PARA FREQ DE 3H         #
##############################################
# same frequency from the model's output
data = df['2006-01-11':].resample('3H').mean()

##############################################
#   LEITURA DOS DADOS DE SAÍDA DO MODELO     #
##############################################
alpha = np.deg2rad(55.)

ncin = xr.open_dataset(MODELO_DIR+simulacao)
prod = load_model(ncin,29,74,10) # i=29,j=74,k=10
prod *= 10 # for a reason I cannot understand

# rotating vectors with alpha = 55. only
prod['cross'],prod['along'] = rotaciona(prod.u.values,prod.v.values,alpha)
newprod =  interpolating_model(prod) # with 1.5H frequency
newprod = newprod.resample('3H').mean()

# filtering
along_filt = filtering(newprod.along.values)
cross_filt = filtering(newprod.cross.values)
df = pd.DataFrame({'along_filt':along_filt,'cross_filt':cross_filt},index=newprod.index)

# ignoring first hours and cutting until last observed time
df = df['2006-01-11':'2006-02-07 09:00']
data = data[:-1]

fig,ax = plt.subplots(nrows=2)
ax[0].plot(data.index,data.along_filt,label='data')
ax[0].plot(df.index,df.along_filt,label='prod')
# prod.along.plot(ax=ax[0],label='prod')
ax[0].legend()
ax[1].plot(data.index,data.cross_filt,label='data')
ax[1].plot(df.index,df.cross_filt,label='prod')

plt.show()
