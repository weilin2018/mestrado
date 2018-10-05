'''
Arquivo para carregar, tratar e plotar os dados de vento de:
    . CFSv2
    . ERA-Interim
    . Laje de Santos


    Filtrar com Hamming, 30h não melhora muito a comparação dos dados.

    # filtrando usando hamming filter, com janela de 30h
    #
    # ncep_filtered = ncep.rolling(30,win_type='hamming').mean()
    # laje_filtered = laje.rolling(30,win_type='hamming').mean()
    #
    # laje_statisticalAnaysis(laje_filtered,ncep_filtered,'Filtrado, Hamming/30h')

'''

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

plt.ion()

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def laje_statisticalAnaysis(laje,cfsv2,whichSerie,labels=['Laje','CFSv2']):
    ''' '''
    # statistical analysis for Laje de Santos
    laje_cut = laje.copy()
    cfsv_cut = cfsv2.copy()

    laje_cut.interpolate(inplace=True,method='linear')
    cfsv_cut.interpolate(inplace=True,method='linear')

    skillWu = oceano.skill_willmott(laje_cut.wu.values, cfsv_cut.wu.values)
    skillWv = oceano.skill_willmott(laje_cut.wv.values, cfsv_cut.wv.values)

    corrWu = calculateCorr(laje_cut.wu.values, cfsv_cut.wu.values)[0]
    corrWv = calculateCorr(laje_cut.wv.values, cfsv_cut.wv.values)[0]

    mseWu  = calculateMSE(laje_cut.wu.values, cfsv_cut.wu.values)
    mseWv  = calculateMSE(laje_cut.wv.values, cfsv_cut.wv.values)

    # start and final date
    tstart = pd.to_datetime(laje_cut.index.values[0])
    tstart = tstart.strftime('%Y.%m.%d')
    tfinal = pd.to_datetime(laje_cut.index.values[-1])
    tfinal = tfinal.strftime('%Y.%m.%d')

    # plot data and skill
    fig, ax = plt.subplots(nrows=2,ncols=1,sharex=True)

    ax[0].plot(laje.wu,label='Laje')
    ax[0].plot(cfsv2.wu,label='CFSv2')
    ax[0].margins(0)
    ax[0].set_ylim(-15,15)
    ax[0].set_title('Along shore')

    wuText = r'Skill: %0.2f | Corr.: %0.2f | MSE: %0.2f' % (skillWu,corrWu,mseWu)
    ax[0].text('2015-06-07', 11.5, wuText, ha='center',va='center',bbox=dict(boxstyle='round', ec=(1.,0.5,0.5), fc=(1.,0.8,0.8)))
    ax[0].legend(labels,loc='lower left')

    ax[1].plot(laje.wv,label='Laje')
    ax[1].plot(cfsv2.wv,label='CFSv2')
    ax[1].margins(0)
    ax[1].set_ylim(-15,15)
    ax[1].set_title('Cross shore')

    wvText = r'Skill: %0.2f | Corr.: %0.2f | MSE: %0.2f' % (skillWv,corrWv,mseWv)
    ax[1].text('2015-06-07', 11.5, wvText, ha='center',va='center',bbox=dict(boxstyle='round', ec=(1.,0.5,0.5), fc=(1.,0.8,0.8)))
    ax[1].legend(labels,loc='lower left')

    # ax[2] = stickplot(laje*(-1),ax[2])
    # ax[2].set_title('Laje')
    # ax[3] = stickplot(cfsv2*(-1),ax[3])
    # ax[3].set_title('CFSv2')

    plt.suptitle('%s v %s [%s to %s]\n[%s]' % (labels[0],labels[1],tstart,tfinal,whichSerie),fontsize=26)
    plt.show()


##############################################################################
#                                STATS                                       #
##############################################################################
def calculateMSE(x,y):
    """Calculate mean squared error.

     Calculate mean squared error using the formula:

            MSE = 1/n * sum((x-y)^2)

    Parameters
    ----------
    x : array
        real data as an array
    y : array
        modeled data as an array

    Returns
    -------
    mse : float
        the mean squared error calculated

    """

    return np.mean((x - y)**2)

def calculateCorr(x,y):
    """Calculate correlation between two series with same size.

    Using the Pearson's method, from scipy.stats.pearsonr, this function
    calculate the correlation between x and y, two array with same size and
    without missing values.

    Parameters
    ----------
    x : array
        real data.
    y : array
        modeled data.

    Returns
    -------
    float
        correlation coefficient.

    """
    import scipy.stats as stats

    return stats.pearsonr(x,y)

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

# carregando dados da Laje de Santos já tratados e salvos em um netcdf pela
# funçao load_lajeSantosData.py
LAJE_DIR = '/media/danilo/Danilo/mestrado/ventopcse/data/lajeSantos_data.nc'
laje = xr.open_dataset(LAJE_DIR) # carregando netcdf
laje = laje.to_dataframe()       # convertendo para pd.DataFrame

# rename columns just for a moment
laje.columns = ['wu','wv']

# carregando dados do CFSv2, do ponto de grade mais próximo da Laje de Santos
NCEP_DIR = '/media/danilo/Danilo/mestrado/ventopcse/data/timeseries_cfsv2_lajedesantos_2015.nc'
ncep = xr.open_dataset(NCEP_DIR) # carregando netcdf
dct = {
    'wu': np.squeeze(ncep['U_GRD_L103'].values),
    'wv': np.squeeze(ncep['V_GRD_L103'].values),
}
ncep = pd.DataFrame(dct,index=ncep.time.values)       # convertendo para pd.DataFrame

"""
    Laje de Santos contempla o período de 2015-04-30 21:00 até 2015-12-31 20:00
        com frequência horária

    CFSv2 contempla o período de 2015-01-01 01:00 até 2015-12-31 19:00
        com frequência de 6h

    Abaixo pega-se as linhas de dados da Laje que batem com o mesmo instante
    dos dados do CFSv2
"""
# resample para alterar frequencia dos dados da Laje
laje = laje.resample('6H').mean()

# recortando para os dados começarem e terminarem no mesmo instante de tempo
ncep = ncep['2015-05-01 07:00':]
laje = laje['2015-05-01 01:00':'2015-12-31 19:00']

"""
    convenção trocada [ncep em convenção oceanográfica, in situ em meteorológica]
"""
# convenção trocada [ncep em convenção oceanográfica, in situ em meteorológica]
laje.wu *= -1

# interpolando a série observada, para preencher lacunas
laje = laje.interpolate(method='cubic')

"""
Rotacionando com -18º
"""
angRot = (54.*np.pi)/180
along,across = oceano.rotaciona(ncep.wu.values,ncep.wv.values,angRot)
ncep_rotacionado = pd.DataFrame({'wv':across,'wu':along},index=ncep.index)

laje_statisticalAnaysis(laje,ncep_rotacionado,'Rotacionado 54deg',['Laje de Santos','CFSv2'])

##############################################################################
#                               PNBOIA/SANTOS                                #
##############################################################################
PNBOIA_DIR = '/media/danilo/Danilo/mestrado/ventopcse/data/pnboiaSantos_data.nc'

boia = xr.open_dataset(PNBOIA_DIR)
boia = boia.to_dataframe()

# carregando dados do CFSv2, do ponto de grade mais próximo da Boia Santos
NCEP_DIR = '/media/danilo/Danilo/mestrado/ventopcse/data/timeseries_cfsv2_pnboia.nc'
ncep = xr.open_dataset(NCEP_DIR) # carregando netcdf
dct = {
    'wu': np.squeeze(ncep['U_GRD_L103'].values),
    'wv': np.squeeze(ncep['V_GRD_L103'].values),
}
ncep = pd.DataFrame(dct,index=ncep.time.values)       # convertendo para pd.DataFrame

# resample dos dados in situ
boia = boia.resample('6H').mean()

# recortando os dados para o mesmo período de informação
boia = boia['2015-05-01 06:00':'2015-12-31 19:00']
ncep = ncep['2015-05-01 07:00':'2015-12-31 19:00']

# convenção trocada [ncep em convenção oceanográfica, in situ em meteorológica]
boia *= -1

# laje_statisticalAnaysis(boia,ncep_rotacionado,'Sem Rotacionar')

angRot = (54.*np.pi)/180
along,across = oceano.rotaciona(ncep.wu.values,ncep.wv.values,angRot)
ncep_rotacionado = pd.DataFrame({'wv':across,'wu':along},index=ncep.index)


laje_statisticalAnaysis(boia,ncep_rotacionado,'Rotacionado 54deg',['PNBOIA/Santos','CFSv2'])


############# FILTRAGEM USANDO HAMMING, JANELA DE 30H

obs = boia.rolling(5,win_type='hamming').mean()
rea = ncep.rolling(5,win_type='hamming').mean()

# laje_statisticalAnaysis(obs['2015-05-08 07:00:00':],rea['2015-05-08 08:00:00':],'Hamming, 30h',['PNBOIA/Santos', 'CFSv2'])
