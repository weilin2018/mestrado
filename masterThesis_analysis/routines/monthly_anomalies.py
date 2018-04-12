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

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/', '/ventopcse/data/pickles/')

def anomaly(climatology,month):
    '''
        climatology: campo de vento médio tomado como estado básico
        monmth: campo de vento médio para o mês em questão
    '''
    return climatology - month
    #return month - climatology

##############################################################################
#                               MAIN CODE                                    #
##############################################################################

# listas para controlar os plots de todos os meses
meses = ['nov', 'dec','jan', 'feb', 'mar']
locs  = [(0,0), (0,2), (0,4), (1,1), (1,3)] # vetor com localizacao do subplots
pickles = ['climatology.pickle', 'summer2014.pickle', 'summer2015.pickle']

# importar as grades que serão utilizadas na interpolação
coarsedFile = glob.glob(BASE_DIR.replace('github/', 'ventopcse/data/CFSR/1992_2011/*.nc'))[0]
refinedFile = glob.glob(BASE_DIR.replace('github/', 'ventopcse/data/CFSv2/verao2014/*.nc'))[0]

coarsedGrid, refinedGrid = oceano.extrair_grades(coarsedFile, refinedFile)

##############################################################################
#                       CLIMATOLOGY - SUMMER/2014                            #
##############################################################################

# definindo títulos dos subplots
dctTitles = {
    'nov': 'Climatology - November/2013',
    'dec': 'Climatology - December/2013',
    'jan': 'Climatology - January/2014',
    'feb': 'Climatology - February/2014',
    'mar': 'Climatology - March/2014'
}

# importar pickles com os campos de vento médio
climatology = pickle.load(open(DATA_DIR+pickles[0],'r'))
summer      = pickle.load(open(DATA_DIR+pickles[1],'r'))

plt.figure(figsize=(16,10))

for mes,loc in zip(meses,locs):

    # importar os dados do mes em analise
    mes_climato = climatology[mes]
    mes_summer  = summer[mes]

    # interpolar os dados para a mesma grade
    wui = oceano.interpolar_grade(coarsedGrid,refinedGrid,mes_summer['wu'])
    wvi = oceano.interpolar_grade(coarsedGrid,refinedGrid,mes_summer['wv'])

    # calculando as anomalias
    wu_anomaly = anomaly(mes_climato['wu'], wui)
    wv_anomaly = anomaly(mes_climato['wv'], wvi)

    # calcular a velocidade anomala
    spd_anomaly = np.sqrt(wu_anomaly**2 + wv_anomaly**2)

    # realizar os plots
    ax = plt.subplot2grid(shape=(2,6), loc=loc, colspan=2)

    m = oceano.make_map(ax)

    contour_levels = np.arange(0,5.001,0.001)

    x,y = np.meshgrid(coarsedGrid['lon'],coarsedGrid['lat'])

    c = m.contourf(x,y,spd_anomaly,contour_levels,latlon=True,extend='max')
    q = m.quiver(x,y,wu_anomaly, wv_anomaly, latlon=True,
                                alpha=.7,scale=120,width=0.005,minshaft=2)
    m.ax.set_title(dctTitles[mes])

    # fazendo a colorbar ser do tamanho de cada subplot
    divider = make_axes_locatable(m.ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.2)
    cb = plt.colorbar(c,orientation='horizontal',ticks=[0,1,2,3,4,5],format='%d',cax=cax)
    # cb = plt.colorbar(c,orientation='horizontal',ticks=[0,2,4,6,8,10],format='%d',
    #     fraction=.057,pad=.06)

    cb.set_label(r'Wind Anomaly [$m.s^{-1}$]',fontsize=8, labelpad=-1)
    cb.ax.tick_params(labelsize=8)

plt.suptitle('Summer 2014',fontsize=24)
plt.show()

##############################################################################
#                       CLIMATOLOGY - SUMMER/2015                            #
##############################################################################

# definindo títulos dos subplots
dctTitles = {
    'nov': 'Climatology - November/2014',
    'dec': 'Climatology - December/2014',
    'jan': 'Climatology - January/2015',
    'feb': 'Climatology - February/2015',
    'mar': 'Climatology - March/2015'
}

# importar pickles com os campos de vento médio
climatology = pickle.load(open(DATA_DIR+pickles[0],'r'))
summer      = pickle.load(open(DATA_DIR+pickles[2],'r'))

plt.figure(figsize=(16,10))

for mes,loc in zip(meses,locs):

    # importar os dados do mes em analise
    mes_climato = climatology[mes]
    mes_summer  = summer[mes]

    # interpolar os dados para a mesma grade
    wui = oceano.interpolar_grade(coarsedGrid,refinedGrid,mes_summer['wu'])
    wvi = oceano.interpolar_grade(coarsedGrid,refinedGrid,mes_summer['wv'])

    # calculando as anomalias
    wu_anomaly = anomaly(mes_climato['wu'], wui)
    wv_anomaly = anomaly(mes_climato['wv'], wvi)

    # calcular a velocidade anomala
    spd_anomaly = np.sqrt(wu_anomaly**2 + wv_anomaly**2)

    # realizar os plots
    ax = plt.subplot2grid(shape=(2,6), loc=loc, colspan=2)

    m = oceano.make_map(ax)

    contour_levels = np.arange(0,5.001,0.001)

    x,y = np.meshgrid(coarsedGrid['lon'],coarsedGrid['lat'])

    c = m.contourf(x,y,spd_anomaly,contour_levels,latlon=True,extend='max')
    q = m.quiver(x,y,wu_anomaly, wv_anomaly, latlon=True,
                                alpha=.7,scale=150,width=0.005,minshaft=2)
    m.ax.set_title(dctTitles[mes])

    # fazendo a colorbar ser do tamanho de cada subplot
    divider = make_axes_locatable(m.ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.2)
    cb = plt.colorbar(c,orientation='horizontal',ticks=[0,1,2,3,4,5],format='%d',cax=cax)

    cb.set_label(r'Wind Anomaly [$m.s^{-1}$]',fontsize=8, labelpad=-1)
    cb.ax.tick_params(labelsize=8)

plt.suptitle('Summer 2015',fontsize=24)
plt.show()
