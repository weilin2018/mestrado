#!/usr/bin/env python2
#-*-coding:utf-8-*-
'''
Calcular produtos mensais do CFSv2 para Novembor, Dezembro, Janeiro, Fevereiro e
Março, do verão de 2015.

Este código:

	. lê os arquivos do diretório de dados para o verão
	. calcula a média mensal
	. plota  os valores para os meses NDJVM
	. salve um arquivo .pickle com essas médias

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

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano


def make_map(ax):

	m = Basemap(projection='merc', llcrnrlat=-30, urcrnrlat=-20, llcrnrlon=-50, urcrnrlon=-40, resolution='l')

	# m = pickle.load(open("/media/danilo/Danilo/mestrado/ventopcse/rotinas/sudesteBR.pkl", "r"))
	m.ax = ax

	m.drawcoastlines(linewidth=.8)
	m.drawmapboundary()

	# definir meridianos e paralelos para plotar no mapa
	meridians=np.arange(-50,-40,3)
	parallels=np.arange(-30,-20,2)
	# desenhar meridianos e paralelos conforme definido acima
	m.drawparallels(parallels,labels=[True,False,False,True],fontsize=13,fontweight='bold',color='gray')
	m.drawmeridians(meridians,labels=[True,False,False,True],fontsize=13,fontweight='bold',color='gray')


	return m

BASE_DIR = oceano.make_dir()

DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/data/CFSv2/verao2015/')

# extrair longitude e latitude
nfiles = glob.glob(DATA_DIR+"*.nc")[0]
ncdata = xr.open_dataset(nfiles)
lon    = ncdata['lon'].values - 360
lat    = ncdata['lat'].values

lon,lat = np.meshgrid(lon,lat)

# listas para controlar os plots de todos os meses
meses = ['nov', 'dec','jan', 'feb', 'mar']
locs  = [(0,0), (0,2), (0,4), (1,1), (1,3)] # vetor com localizacao do subplots
dates = ['201411','201412','201501','201502','201503']

dctTitles = {
    'nov': 'November/2014',
    'dec': 'December/2014',
    'jan': 'January/2015',
    'feb': 'February/2015',
    'mar': 'March/2015'
}

# dicionario para armazenar os dados
saveData = {}

fig = plt.figure(figsize=(16,8))

for mes,loc,date in zip(meses,locs,dates):

    # ler os dados do mes em loop
    wu,wv = oceano.read_month(date,DATA_DIR)

    # calcular velocidade
    spd = np.sqrt(wu**2 + wv**2)

    data = { 'wu': wu, 'wv': wv }  	# organizando dados para armazenamento
    saveData[mes] = data   			# dicionario para salvar os dados em pickle

    # realizar os plots
    ax = plt.subplot2grid(shape=(2,6), loc=loc, colspan=2)

    m = make_map(ax)

    contour_levels = np.arange(0,10.001,0.001)

    c = m.contourf(lon,lat,spd,contour_levels,latlon=True,extend='max')
    q = m.quiver(lon[::3,::3],lat[::3,::3],wu[::3,::3], wv[::3,::3], latlon=True,
                                alpha=.7,scale=150,width=0.005,minshaft=2)
    m.ax.set_title(dctTitles[mes])

    cb = plt.colorbar(c,orientation='horizontal',ticks=[0,2,4,6,8,10],format='%d',
        fraction=.057,pad=.06)
    cb.set_label(r'Wind [$m.s^{-1}$]',fontsize=8, labelpad=-1)
    cb.ax.tick_params(labelsize=8)

plt.show()


# salvar os dados em um pickle finalmente
fname = BASE_DIR.replace('github/', '/ventopcse/data/pickles/summer2015.pickle')

pickle.dump(saveData, open(fname,'w'))
