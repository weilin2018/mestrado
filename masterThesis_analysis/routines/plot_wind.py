#!/usr/bin/env python2
#-*-coding:utf-8-*-
''' read and plot CSFR grib file '''

import glob
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib
matplotlib.style.use('ggplot')

# directory with data
dataDir = "/home/tparente/danilo/mestrado/ventopcse/data/csfr2/"
# define coords location to extract data
# pnboia Santos location is: Latitude: -25.27 Longitude: -44.93
pnboiaSan = [-25.27, -44.93]

# import a grib file to extract latlon's index
grib = xr.open_dataset(dataDir+"cdas1.20170208.pgrbh.grb2",engine='pynio')

# extract latlon
lat = grib['lat_0'].data
lon = grib['lon_0'].data - 360

# find index
ilat = np.where(lat == -25.50)[0][0] # 13
ilon = np.where(lon == -45.00)[0][0] # 12

# select date interval
init_date = '20140101'
fina_date = '20161231'

### plotando 2014
print('plotando dados de 2014')

arquivo = glob.glob(dataDir+'cdas1.2014*')
arquivo.sort() # to organize files 

# grib = xr.open_data(arq,engine='pynio')
# uw = grib['UGRD_P0_L103_GLL0'].data[:,ilat,ilon]

# create numpy.ndarray to store all values
uwnd = []
vwnd = []

for arq in arquivo:
	grib = xr.open_dataset(arq,engine='pynio')
	uw = grib['UGRD_P0_L103_GLL0'].data[0,ilat,ilon]
	for i in range(1,4)
	if grib['VGRD_P0_L103_GLL1'].data.shape[0] == 2:
		vw = grib['VGRD_P0_L103_GLL1'].data[0,ilat,ilon]
	else:
		vw = grib['VGRD_P0_L103_GLL1'].data[ilat,ilon]

	uwnd.append(uw)
	vwnd.append(vw)

dias = np.arange(1,366,1)

'''
lat0,lon0 => somente para componente U 

lat1,lon1 até lat4,lon4 => componente V com variação de quantidade de acordo ao mês dos dados

Jan14 = 4v,1u


'''

# separando os meses de 2014
jan14 = arquivo[:31]
fev14 = arquivo[32:59] 
mar14 = arquivo[60:90]
abr14 = arquivo[91:120]
mai14 = arquivo[121:151]
jun14 = arquivo[152:181]
jul14 = arquivo[182:212]
ago14 = arquivo[213:243]
set14 = arquivo[244:273]
out14 = arquivo[274:304]
nov14 = arquivo[305:334]
dez14 = arquivo[335:]

'''
Rodar for para cade dia em um mes e verificar o tamanho len(grib.dims):
	se 10 => uma das variáveis V possui dois dados, então chamar função para extrair
	se 11 => possui 4 variáveis de V, chamar função para extrair

dims = gr['v'].data.shape 
if dims.shape == 3:
    # loop para horários 
else:
    # armazenar info⁠⁠⁠⁠

'''

def extractVW(grib, ilat, ilon, vars=4):
	''' 
		extract V wind component from grib file 
		ilat,ilon = indices para extrair dados
	'''
	vw_tmp = []
	for i in range(1,vars+1):
		parametro = 'VGRD_P0_L103_GLL'+str(i)
		if vars == 4:
			v = grib[parametro].data[ilat,ilon]
		elif vars == 3:
			# realizar testes para saber qual a dimensão das variáveis
			v_dim = grib[parametro].data.shape
			if len(v_dim) == 3:
				# loop para pegar os dois valores de v
			else:
				# pega somente um valor
				v = grib[parametro].data[ilat,ilon]

		vw_tmp.append(v)

		return vw_tmp # retorna uma lista com os 4 valores de V para cada instante [00, 06, 12, 18]





for f in arquivo:
	# lendo arquivo
	grib = xr.open_dataset(f, engine='pynio')

	# testar dimensoes
	if len(grib.dims) == 11:
		# chamar funcao para tratar dados de V com 4 variáveis
	elif len(grib.dims) == 10:
		# chamar funcao para tratar dados de V com 3 variáveis

