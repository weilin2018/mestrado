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


def calculate_monthlyMean(ncdata,axis=0):
	''' 
		Calcula a media no eixo axis (0 é o default)
	'''
	wumean 		= ncdata['U_GRD_L103'].values
	wvmean 		= ncdata['V_GRD_L103'].values
	wumean 		= wumean.mean(axis=axis)
	wvmean 		= wvmean.mean(axis=axis)

	return wumean,wvmean


BASE_DIR = oceano.make_dir()

DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/data/data4climatology/CFSR_climatologia/')


nfiles = glob.glob(DATA_DIR + '*.nc')			# reading all netcdf files
nfiles.sort() 									# sort by name

#dec = np.zeros([31,32,33])						# 31 years, 32,33 coordinates
janwu = np.zeros([31,32,33])
janwv = np.zeros([31,32,33])
febwu = np.zeros([31,32,33])
febwv = np.zeros([31,32,33])
marwu = np.zeros([31,32,33])
marwv = np.zeros([31,32,33])

cont = 0

for i,k in zip(np.arange(0,95,3), np.arange(0,31,1)):
	# ler netcdf adequado
	ncJan = xr.open_dataset(nfiles[i])
	ncFeb = xr.open_dataset(nfiles[i+1])
	ncMar = xr.open_dataset(nfiles[i+2])

	wumean,wvmean = calculate_monthlyMean(ncJan)

	janwu[k,:,:] = wumean[:,:]
	janwv[k,:,:] = wvmean[:,:]

	wumean,wvmean = calculate_monthlyMean(ncFeb)

	febwu[k,:,:] = wumean[:,:]
	febwv[k,:,:] = wvmean[:,:]

	wumean,wvmean = calculate_monthlyMean(ncMar)

	marwu[k,:,:] = wumean[:,:]
	marwv[k,:,:] = wvmean[:,:]

# calculando os produtos mensais
janwu_mean = janwu.mean(axis=0)
janwv_mean = janwv.mean(axis=0)

febwu_mean = febwu.mean(axis=0)
febwv_mean = febwv.mean(axis=0)

marwu_mean = marwu.mean(axis=0)
marwv_mean = marwv.mean(axis=0)

# plotando os produtos mensais

os.system('clear')
print('Plotando dados de Janeiro')
plt.quiver(lon,lat, janwu_mean, janwv_mean)
plt.title('Janeiro [1979 - 2010]')
plt.show()

print('Plotando dados de Fevereiro')
plt.quiver(lon,lat, febwu_mean, febwv_mean)
plt.title('Fevereiro [1979 - 2010]')
plt.show()

print('Plotando dados de Marco')
plt.quiver(lon,lat, marwu_mean, marwv_mean)
plt.title(u'Março [1979 - 2010]')
plt.show()
