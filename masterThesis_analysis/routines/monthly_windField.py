#!/usr/bin/env python2
#-*-coding:utf-8-*-
''' read and plot CSFR grib file '''

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

DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/data/data4climatology/CFSR_climatologia/')
SAVE_DIR = BASE_DIR+'masterThesis_analysis/figures/monthly_wind_field/1979_2010/'

nfiles = glob.glob(DATA_DIR + '*.nc')			# reading all netcdf files
nfiles.sort() 									# sort by name

# importing basemap files
FILE_MAP = BASE_DIR + 'masterThesis_analysis/routines/sudesteBR.pkl'

# verificando arquivo a arquivo
i = 0 			# contador de imagens

for fname in nfiles[:2]:
	print('Reading file: %s'%fname)
	ncdata 		= xr.open_dataset(fname) 		# read netcdf file
	wumean 		= ncdata['U_GRD_L103'].values
	wvmean 		= ncdata['V_GRD_L103'].values
	wumean 		= wumean.mean(axis=0)
	wvmean 		= wvmean.mean(axis=0)

	lon,lat 	= np.meshgrid(ncdata['lon'].values, ncdata['lat'].values)

	fig,ax = plt.subplots(figsize=(12,8))
	m = pickle.load(open(FILE_MAP, "r"))

	m.ax = ax 

	m.drawcoastlines(linewidth=0.1)
	m.drawmapboundary(fill_color='#e5f2ff')
	m.fillcontinents(color='#ffd480')

	q = m.quiver(lon[::2,::2],lat[::2,::2],wumean[::2,::2],wvmean[::2,::2],latlon=True,pivot='middle')

	# criando o titulo YEAR/MONTH baseado no time do netcdf
	time  = ncdata['time'].values[0]
	year  = pd.to_datetime(time).year
	month = pd.to_datetime(time).month

	plt.title('%s-%s'%(year,month))
	outFig = SAVE_DIR+'%s.png' % str(i).zfill(6)
	#plt.savefig(outFig)
	plt.show()

	i += 1



'''
Algoritmo:

. ler todos os arquivos .nc do diretorio devido
. gridar dados
. tirar media mensal
. salvar imagem
'''