#!/usr/bin/env python2
#-*-coding:utf-8-*-
'''
Calcular produtos mensais do CFSR para os meses de Janeiro, Fevereiro e Março,
contemplando todo o período do conjunto de dados de reanálise (1979-2010).

Importante:

produto baixado em: https://rda.ucar.edu/datasets/ds093.1/index.html,
para uma área delimitada por:

		-19
-57				-37
		-29


Ainda preciso calcular o mesmo valor para Novembro e Dezembro, de 1979 a 2010


'''

import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import os
import pickle
from scipy.interpolate import griddata
from mpl_toolkits.basemap import Basemap

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

################################################################################
#  								                 FUNCOES               									     #
################################################################################

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


################################################################################
#                  								MAIN CODE 									                 #
################################################################################

BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/data/FTP/recortados/')

# importar latitude e longitude de um arquivo pickle
coarseGrid = pickle.load(open(DATA_DIR.replace('recortados/', 'coordenadas.pickle'), 'r'))
lon = coarseGrid['lon']
lat = coarseGrid['lat']
x,y = np.meshgrid(lon,lat)

meses = ['nov', 'dec','jan', 'feb', 'mar']
locs  = [(0,0), (0,2), (0,4), (1,1), (1,3)] # vetor com localizacao do subplots

dctTitles = {
    'nov': 'Nov/1992 - Nov/2010',
    'dec': 'Dez/1992 - Dez/2010',
    'jan': 'Jan/1993 - jan/2010',
    'feb': 'Fev/1993 - Fev/2010',
    'mar': 'Mar/1993 - Mar/2010'
}

fig = plt.figure(figsize=(16,8))

for mes,loc in zip(meses, locs):

    # para o mes de novembro
    nfiles = glob.glob(DATA_DIR+mes+'/*.nc')
    nfiles.sort()

    # determinando dimensoes das matrizes
    fmax = len(nfiles)
    imax, jmax = 32, 32

    mesWu, mesWv = np.zeros([fmax,imax,jmax]), np.zeros([fmax,imax,jmax])
    cont = 0

    for f in nfiles:
        ncdata = xr.open_dataset(f)
		# extrair dados
        wu     = np.squeeze(ncdata['10u'].values[::7,:,:])
        wv     = np.squeeze(ncdata['10v'].values[::7,:,:])
		# tomar média diária
        wu     = wu.mean(axis=0)
        wv     = wv.mean(axis=0)

        mesWu[cont,:,:] = wu[:,:]
        mesWv[cont,:,:] = wv[:,:]

        cont += 1
	# tomar média mensal
    mes_wu_mean = mesWu.mean(axis=0)
    mes_wv_mean = mesWv.mean(axis=0)

    ax = plt.subplot2grid(shape=(2,6), loc=loc, colspan=2)

    m = make_map(ax)

    contour_levels = np.arange(0,6.001,0.001)

    c = m.contourf(x,y,np.sqrt(mes_wu_mean**2 + mes_wv_mean**2), contour_levels,
        latlon=True,extend='max')
    q = m.quiver(x[::3,::3],y[::3,::3],wu[::3,::3], wv[::3,::3], latlon=True,
                                    alpha=.7,scale=150,width=0.005,minshaft=2)
    m.ax.set_title(dctTitles[mes])

    cb = plt.colorbar(c,orientation='horizontal',ticks=[0,2,4,6],format='%d',
        fraction=.057,pad=.06)
    cb.set_label(r'Wind [$m.s^{-1}$]',fontsize=8, labelpad=-1)
    cb.ax.tick_params(labelsize=8)


plt.show()
