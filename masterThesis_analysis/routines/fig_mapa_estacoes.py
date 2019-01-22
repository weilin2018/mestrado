#

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
import cmocean as cmo

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

plt.ion()
##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/output')

### FIGURE CONFIGURATION
figsize      = (8.4/2.54, 10./2.54)
cax_position = [0.205,0.12,0.76,0.02]

# localizacoes [lat,lon]
loc_Eslaje      = [-24.32000, -46.180, 'Laje de Santos']
loc_pnboia      = [-25.28000, -44.925, 'PNBOIA/Santos']
loc_ncep_laje   = [-24.22487, -46.228, 'NCEP/CFSv2']
loc_ncep_pnboia = [-25.24701, -45.000, 'CFSv2 - PNBOIA']

# importar batimetria de produto de modelo qualquer
ncin = xr.open_dataset(DATA_DIR + 'EC1.cdf')
depth = ncin.depth.values
lon   = ncin.lon.values
lat   = ncin.lat.values
# tratar as coordenadas
lon[lon == 0] = np.nan
lat[lat == 0] = np.nan

ncin.close()

fig,ax = plt.subplots(figsize=figsize) #ja em polegadas para artigo
# cax = fig.add_axes(cax_position)

m = oceano.make_map(ax,ulat=-23.,llat=-26.,llon=-47.,ulon=-44.,resolution='f',nmeridians=1,nparallels=1)
contours = np.arange(5,700.,0.1)
cf = m.contourf(lon,lat,depth,contours,latlon=True,cmap=cmo.cm.deep,extend='max')
cs = m.contour(lon,lat,depth,levels=[50.,100.,300.],colors=('k'),linestyle=('--'),latlon=True)
plt.clabel(cs, fontsize=9, inline=1,fmt='%i')
# cbar = plt.colorbar(cf,cax=cax,orientation='horizontal')
# cbar.set_label('Profundidade (m)')

m.scatter(loc_Eslaje[1],loc_Eslaje[0],s=30,c='b',marker='*',latlon=True,label=loc_Eslaje[2])
m.scatter(loc_pnboia[1],loc_pnboia[0],s=30,c='b',marker='s',latlon=True,label=loc_pnboia[2])

m.scatter(loc_ncep_laje[1],loc_ncep_laje[0],s=30,c='r',marker='^',latlon=True,label=loc_ncep_laje[2])
m.scatter(loc_ncep_pnboia[1],loc_ncep_pnboia[0],s=30,c='r',marker='^',latlon=True)

plt.legend(loc='best')
plt.subplots_adjust(top=0.962,bottom=0.038,left=0.156,right=0.995,hspace=0.2,wspace=0.2)

# plt.savefig('/media/danilo/Danilo/mestrado/github/masterThesis_analysis/figures/mapas/localizacao_pontos.png',dpi=300)
