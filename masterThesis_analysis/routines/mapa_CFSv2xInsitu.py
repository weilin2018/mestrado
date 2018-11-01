# add some description here

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

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

SAVE_DIR = '/media/danilo/Danilo/mestrado/github/masterThesis_analysis/figures/mapas/'
FIGNAME  = 'localizacao_pontos.png'

# load bathymetry
ncin = xr.open_dataset('/media/danilo/Danilo/mestrado/ventopcse/output/exp12.cdf')
depth= ncin.depth.values
lon  = ncin.lon.values
lat  = ncin.lat.values
lon[lon==0.]=np.nan
lat[lat==0.]=np.nan

# define locations and labels
locations_insitu = {
    'Laje de Santos':      [-24.32,-46.18,'*','b'],
    'PNBOIA/Santos':       [-25.28,-44.925,'s','b'],
}

locations_ncep = {
    'ncep_laje':           [-24.225,-46.277],
    'ncep_boia':           [-25.247,-45.]
}

# contour
contours = np.arange(5.,1000.,1.)

fig,ax = plt.subplots()

m = oceano.make_map(ax,ulat=-23.,llat=-26.,ulon=-44.,llon=-47.,nmeridians=1,nparallels=1,resolution='f')
cf = m.contourf(lon,lat,depth,contours,latlon=True,cmap=cmo.cm.deep,extend='max')
cs = m.contour(lon,lat,depth,latlon=True,levels=[50.,100.,300.],colors=('k'))
plt.clabel(cs,fmt='%i')

for key in locations_insitu.keys():
    m.scatter(locations_insitu[key][1],locations_insitu[key][0],s=30,c=locations_insitu[key][3],marker=locations_insitu[key][2],latlon=True,label=key)

m.scatter(-46.277,-24.225,s=30,c='r',marker='^',latlon=True)
m.scatter(-45.000,-25.247,s=30,c='r',marker='^',latlon=True,label='NCEP/CFSv2')

plt.legend()

plt.savefig(SAVE_DIR+FIGNAME,dpi=300)
os.system('convert -trim %s %s'%(SAVE_DIR+FIGNAME,SAVE_DIR+FIGNAME))
