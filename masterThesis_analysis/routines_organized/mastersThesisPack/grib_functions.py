#!/usr/bin/env python2
#-*-coding:utf-8-*-

'''
only works if you run this in grib2 environment

'''

# arquivo com funcao para ler e recortar arquivos grib2, salvanod em seguinda em netcdf
import xarray as xr
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy import signal, fftpack
import scipy
import socket
import glob

# gerar diretorio base
def make_dir():
    '''
        Funcao para gerar o BASE_DIR baseado no computador em que esta rodando
    '''

    hostname = socket.gethostname()

    if hostname == 'oceano': # estou no meu pc pessoal
        BASE_DIR = '/media/danilo/Danilo/mestrado/github/'

        return BASE_DIR
    if hostname == 'tripoli':
        BASE_DIR = '/home/tparente/danilo/mestrado/github/'

        return BASE_DIR

# encontrar indices dos pontos mais proximo a uma coordenada
def find_nearest(lon,lat,ilon,ilat):
    '''
        lon,lat = lat e lon da grade
        ilon,ilat = ponto a ser encontrado
    '''

    # localizacao do terminal da ilha guaiba

    lo = lon.ravel()
    la = lat.ravel()

    coords = []

    for i,j in zip(la,lo):
        coords.append([i,j])

    coords = np.array(coords)

    locations_name = ['Terminal Ilha Guaiba']
    locations_posi = [[ilat,ilon]]

    locs = np.asarray(locations_posi)

    tree = cKDTree(coords)
    # procura em tree os pontos mais próximos dos pontos definidos acima
    dists,indexes = tree.query(locs,k=1)

    pontos = []

    for index in indexes:
        pontos.append(coords[index])

    # converter de lista para array
    pontos = np.asarray(pontos)

    # findind indexes from lat and lon
    ind = []

    for p in pontos:
        ind.append(np.where(lon == p[1]))

    ind = np.asarray(ind)

    # vetores para separar i e j para facilitar resgatar os dados de concentração
    iss=[]
    jss=[]

    for i,j in ind:
        iss.append(int(i))
        jss.append(int(j))

    return iss,jss


def read_grib2(nfiles,lats,lons,output):
	''' 

		@args
			nfiles: lista gerada pelo glob, contendo arquivos grib2
			lats: lista contendo a latitude superior e inferior para recortar
			lons: lista contendo a longitude oeste e leste para recortar
			output: diretorio para se salvar o arquivo netcdf gerado


	'''

	# read each grib file in a loop

	for f in nfiles:

		fname = f.split('/')[-1] 			# obter somente nome do arquivo
		fname = fname.replace('grb2','nc')	# trocar a extensão do arquivo

		# ponto superior da area selecionada
		ulat, ulon = lats[0], lons[1]

		# ponto inferior da area selecionada
		llat, llon = lats[1], lons[0]

		# importar o arquivo grib
		ncdata = xr.open_dataset(f, engine='pynio')

		# extrair variaveis
		lon, lat = ncdata['lon_0'].values, ncdata['lat_0'].values

		# converter lon de 0/360 para -180/180
		lon = lon - 180

		# determinando os indices correspondem a area selecionada, para facilitar
		# na hora de extrair os dados de vento e consumir menos espaço
		ilon = lon[lon < ulon]
		ilon = ilon[ilon > llon]
		ilat = lat[lat < ulat]
		ilat = ilat[ilat > llat]



















