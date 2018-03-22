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
import pandas as pd
import os
import pickle

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano


BASE_DIR = oceano.make_dir()

DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/data/CSFv2_wind/verao2014_completo/')


nfiles = glob.glob(DATA_DIR + '*.nc')			# reading all netcdf files
nfiles.sort() 									# sort by name

def media(sufixo, DATA_DIR):

	nfiles = glob.glob(DATA_DIR + sufixo)
	nfiles.sort()

	wu = np.zeros([len(nfiles), 21, 41])
	wv = np.zeros([len(nfiles), 21, 41])

	cont = 0

	for f in nfiles:
		ncdata = xr.open_dataset(f)
		u = ncdata['U_GRD_L103'].values[:9,:,:]
		v = ncdata['V_GRD_L103'].values[:9,:,:]

		wu[cont,:,:] = u.mean(axis=0)
		wv[cont,:,:] = v.mean(axis=0)

		cont += 1

	return wu.mean(axis=0), wv.mean(axis=0)


# calculando media para Novembro
novwu, novwv = media(sufixo='cdas1.201311*', DATA_DIR=DATA_DIR)

# calculando media para Dezembro
dezwu, dezwv = media(sufixo='cdas1.201312*', DATA_DIR=DATA_DIR)

# calculando media para Janeiro
janwu, janwv = media(sufixo='cdas1.201401*', DATA_DIR=DATA_DIR)

# calculando media para Fevereiro
febwu, febwv = media(sufixo='cdas1.201402*', DATA_DIR=DATA_DIR)

# calculando media para Marco
marwu, marwv = media(sufixo='cdas1.201403*', DATA_DIR=DATA_DIR)

ncdata = xr.open_dataset(nfiles[0])
lon    = ncdata['lon'].values
lat    = ncdata['lat'].values

del ncdata

lon, lat = np.meshgrid(lon,lat)

# plotando os produtos mensais

os.system('clear')
print('Plotando dados de Novembro')
plt.quiver(lon,lat, novwu, novwv)
plt.title('Novembro [2013]')
plt.show()

print('Plotando dados de Dezembro')
plt.quiver(lon,lat, dezwu, dezwv)
plt.title('Dezembro [2012]')
plt.show()

print('Plotando dados de Janeiro')
plt.quiver(lon,lat, janwu, janwv)
plt.title('Janeiro [2014]')
plt.show()

print('Plotando dados de Fevereiro')
plt.quiver(lon,lat, febwu, febwv)
plt.title('Fevereiro [2014]')
plt.show()

print('Plotando dados de Marco')
plt.quiver(lon,lat, marwu, marwv)
plt.title(u'Março [2014')
plt.show()
