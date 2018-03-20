#!/usr/bin/env python2
#-*-coding:utf-8-*-
''' 

	PROGRAMA PARA LER E PLOTAR O CAMPO DE VENTOS EXTRAÍDOS DO MODELO GLOBAL INTERIM

	arquivo baixado em NetCDF, localizado no diretório:
		/home/tparente/danilo/mestrado/ventopcse/data/2014_2016.nc
 '''

import glob
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import xray as xr
import pandas as pd
import matplotlib
import pickle 
import datetime

matplotlib.style.use('ggplot')

fname = xr.open_dataset('/home/tparente/danilo/mestrado/ventopcse/data/2012_2014.nc')

lati = fname['latitude'].data
long = fname['longitude'].data - 360
time = fname['time'].data 
wndU = fname['u10'].data 
wndV = fname['v10'].data 

# gridar coordenadas
lon,lat = np.meshgrid(long,lati)


# gerar instância do Basemap para gerar e salvar em pickle

# definir região
ur = [-17,-35] # upper right
ll = [-31,-53] # lower left

# m = Basemap(projection='merc', llcrnrlat=ll[0], urcrnrlat=ur[0], llcrnrlon=ll[1], urcrnrlon=ur[1], resolution='f')
# pickle.dump(m, open("sudesteBR.pkl","wb"))

# definir meridianos e paralelos para plotar no mapa
# meridians=np.arange(-53,-53,5)
# parallels=np.arange(-31,-17,5)
# desenhar meridianos e paralelos conforme definido acima
# m.drawparallels(parallels,labels=[True,False,False,True],fontsize=13,fontweight='bold',color='gray')
# m.drawmeridians(meridians,labels=[True,False,False,True],fontsize=13,fontweight='bold',color='gray')

# for i in range(len(time)):
# 	plt.quiver(wndU[i,:,:], wndV[i,:,:])
# 	outFig = '/home/tparente/danilo/mestrado/ventopcse/rotinas/output/%s.png' % str(i).zfill(6)
# 	plt.savefig(outFig)
# 	plt.close()


for i in range(len(time)):

	fig,ax = plt.subplots(figsize=(12,8))
	m = pickle.load(open("/home/tparente/danilo/mestrado/ventopcse/rotinas/sudesteBR.pkl", "r"))
	m.ax = ax 

	x,y = m(lon,lat)

	m.drawcoastlines()
	m.drawmapboundary(fill_color='#e5f2ff')
	m.fillcontinents(color='#ffd480')

	wu = wndU[i,:,:]
	wv = wndV[i,:,:]

	skip = 3

	q=m.quiver(x[::skip,::skip],y[::skip,::skip],wu[::skip,::skip],wv[::skip,::skip], pivot='middle')

	#plt.tight_layout()
	plt.title(time[i])

	outFig = '/home/tparente/danilo/mestrado/ventopcse/rotinas/output/%s.png' % str(i).zfill(6)
	plt.savefig(outFig)
	plt.close()

################################################################


def spdir2uv(spd, ang, deg=False):
	if deg:
		ang = np.deg2rad(ang)

	u = spd * np.sin(ang)
	v = spd * np.cos(ang)

	return u,v

def flip_grid(var, lons):
    fltr = lons >= 180
    # fltr =  [False False False ... True  True  True]
    newlons = np.concatenate(((lons - 360)[fltr], lons[~fltr]), axis=-1)
    # newlons = [-180 -177.5 -175 ... -5 -2.5 ] concatenated with [0 2.5 5 ... 175 177.5]
    # newlons = [-180 -177.5 -175 ... 175 177.5 180]
    if var.ndim == 2:
        newvar = np.concatenate((var[:, fltr], var[:, ~fltr]), axis=-1)
    elif var.ndim == 3:
        newvar = np.concatenate((var[:, :, fltr], var[:, :, ~fltr]), axis=-1)
    elif var.ndim == 4:
        newvar = np.concatenate((var[:, :, :, fltr], var[:, :, :, ~fltr]), axis=-1)        
        
    return newvar, newlons

