'''

do jeito que está, o código só calcula a EOF para um mês.

eu preciso elaborar uma rotina para ler os dados e ir tomando médias
mensais e alocando como numa matriz multidimensional, onde o eixo 0 é
o time.

aí eu posso falar de fazer uma EOF para obter uma climatologia.

mas para começar a brincar tá legal!!!

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


from eofs.xarray import Eof
from eofs.examples import example_data_path


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

DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/data/CFSR/1992_2011/')

# extrair longitude e latitude
nfiles = glob.glob(DATA_DIR+"*.nc")[0]
ncdata = xr.open_dataset(nfiles)
lon    = ncdata['lon'].values - 360
lat    = ncdata['lat'].values

lon,lat = np.meshgrid(lon,lat)

locs  = [(0,0), (0,2), (0,4)] # vetor com localizacao do subplots

wu = ncdata['U_GRD_L103'].values

# Create an EOF solver to do the EOF analysis. Square-root of cosine of
# latitude weights are applied before the computation of EOFs.
coslat = np.cos(np.deg2rad(ncdata['lat'].values))
wgts = np.sqrt(coslat)[..., np.newaxis]
solver = Eof(wu, weights=wgts)

# Retrieve the leading EOF, expressed as the correlation between the leading
# PC time series and the input SST anomalies at each grid point, and the
# leading PC time series itself.
eof1 = solver.eofsAsCorrelation(neofs=3)
pc1 = solver.pcs(npcs=2, pcscaling=1)

comp = 0

for loc in locs:

    # realizar os plots
    ax = plt.subplot2grid(shape=(2,6), loc=loc, colspan=2)

    m = make_map(ax)
    x,y = m(lon,lat)

    # Plot the leading EOF expressed as correlation in the Pacific domain.
    clevs = np.linspace(-1, 1, 11)
    fill = m.contourf(x,y,eof1[comp],levels=clevs,cmap=plt.cm.RdBu_r)

    cb = plt.colorbar(fill, orientation='horizontal')
    cb.set_label('correlation coefficient', fontsize=12)
    ax.set_title('EOF1 expressed as correlation', fontsize=13)

    comp += 1

ax = plt.subplot2grid(shape=(2,6), loc=(1,1), colspan=2)
ax.plot(pc1[:,0],color='k',linewidth=2)
ax.axhline(0, color='k')
ax.set_ylim(-3, 3)
ax.set_xlabel('Year')
ax.set_ylabel('Normalized Units')
ax.set_title('PC1 Time Series', fontsize=13)

ax = plt.subplot2grid(shape=(2,6), loc=(1,3), colspan=2)
ax.plot(pc1[:,1],color='k',linewidth=2)
ax.axhline(0, color='k')
ax.set_ylim(-3, 3)
ax.set_xlabel('Year')
ax.set_title('PC2 Time Series', fontsize=13)

plt.suptitle('Zonal component of 10m height wind from CFSv2',fontsize=25)
plt.show()

################################################################################
################################################################################
################################################################################
###### MULTIVARIADA
################################################################################
################################################################################
################################################################################

locs  = [(0,0), (0,1), (1,0),(1,1)] # vetor com localizacao do subplots


from eofs.multivariate.standard import MultivariateEof

wu = ncdata['U_GRD_L103'].values
wv = ncdata['V_GRD_L103'].values

solver = MultivariateEof((wu,wv),weights=None)

# Retrieve the leading EOF, expressed as the correlation between the leading
# PC time series and the input SST anomalies at each grid point, and the
# leading PC time series itself.
eof1 = solver.eofsAsCorrelation(neofs=2)
pc1 = solver.pcs(npcs=2, pcscaling=1)

fig,axes = plt.subplots(nrows=2,ncols=2,figsize=(15,15))

# 1st mode of WU
m = make_map(axes[0,0])
x,y = m(lon,lat)

# Plot the leading EOF expressed as correlation in the Pacific domain.
clevs = np.linspace(-1, 1, 11)
m.contourf(x,y,eof1[0][0,:,:],levels=clevs,cmap=plt.cm.RdBu_r)

# cb = plt.colorbar(fill, orientation='horizontal')
# cb.set_label('correlation coefficient', fontsize=12)
axes[0,0].set_title('1st mode for U-comp', fontsize=13)

#1st mode of WV
m = make_map(axes[0,1])

# Plot the leading EOF expressed as correlation in the Pacific domain.
clevs = np.linspace(-1, 1, 11)
m.contourf(x,y,eof1[1][0,:,:],levels=clevs,cmap=plt.cm.RdBu_r)

# cb = plt.colorbar(fill, orientation='horizontal')
# cb.set_label('correlation coefficient', fontsize=12)
axes[0,1].set_title('1st mode for V-comp', fontsize=13)

# 2nd mode of WU
m = make_map(axes[1,0])

# Plot the leading EOF expressed as correlation in the Pacific domain.
clevs = np.linspace(-1, 1, 11)
m.contourf(x,y,eof1[0][1,:,:],levels=clevs,cmap=plt.cm.RdBu_r)

# cb = plt.colorbar(fill, orientation='horizontal')
# cb.set_label('correlation coefficient', fontsize=12)
axes[1,0].set_title('2nd for U-comp', fontsize=13)

#2nd mode of WV
m = make_map(axes[1,1])

# Plot the leading EOF expressed as correlation in the Pacific domain.
clevs = np.linspace(-1, 1, 11)
m.contourf(x,y,eof1[1][1,:,:],levels=clevs,cmap=plt.cm.RdBu_r)

# cb = plt.colorbar(fill, orientation='horizontal')
# cb.set_label('correlation coefficient', fontsize=12)
axes[1,1].set_title('2nd for V-comp', fontsize=13)


plt.suptitle('Multivariate EOF for Zonal and Meridional \n Components of 10m height Wind',fontsize=25)
plt.show()
