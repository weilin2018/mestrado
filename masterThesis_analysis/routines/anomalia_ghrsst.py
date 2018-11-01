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
from scipy import interpolate
import cmocean as cmo

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def load_climato(fname):

    modeled = xr.open_dataset(fname)

    lonMod,latMod = modeled.lon.values,modeled.lat.values
    lonMod[lonMod == 0.] = np.nan
    latMod[latMod == 0.] = np.nan
    clim = modeled.temp[-1,0,:,:]

    return lonMod,latMod,clim,modeled.depth.values

def interp_ghrsst_to_ecom(x,y,data,xi,yi,cut=True,method='linear'):
    """
        x,y = grade original
        data= dado a ser interpolado
        xi,yi = grade para interpolar
    """

    # selecting some values inside a box to speed up interpolation
    if cut:
        ind = (x > -48.8) & (y < -21.663) & (y > -29.669)
        x = x[ind]
        y = y[ind]
        data = data[ind]

    x = x.ravel()
    y = y.ravel()

    interp = interpolate.griddata((x,y),data.ravel(),(xi,yi),method=method)

    return interp

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
GHRSST_DIR = BASE_DIR.replace('github/','ventopcse/data/GHRSST/')

FIGU_DIR = BASE_DIR + 'masterThesis_analysis/figures/experiments_outputs/ghrsst_v_ecom/'

### loading climatology used in the numerical experiments
lonClim,latClim,sstClim,depth = load_climato(DATA_DIR+'warmupControle.cdf')

### loading satellite sst
fname = GHRSST_DIR + 'ghrsst_JF2014.nc'
ncin = xr.open_dataset(fname)

lon = ncin.lon
lat = ncin.lat
lonSat,latSat = np.meshgrid(lon.values,lat.values)
# convertendo de Kelvin para Celsius
sst  = ncin['analysed_sst'] - 273.15

# media no periodo de interesse
sst_mean = sst.mean(dim='time')

# interpolando dados de satelite para a grade do modelo
import mat4py
data = mat4py.loadmat('/media/danilo/Danilo/mestrado/ventopcse/data/GHRSST/ghrsst_interpoladoSBB.mat')
sstSat = data['sstSat']

#sstSat = interp_ghrsst_to_ecom(lonSat,latSat,sst_mean.values,lonClim,latClim)

# calculando anomalia
sstAnom = sstSat - sstClim

##############################################
#     plotando anomalia meia     #
##############################################
fig,axes,m_axes,cbaxes = oceano.createPlot_structure_horizontal(ncols=2,figsize=(15,15))

cf1 = m_axes[0].contourf(lonSat,latSat,sst_mean.values,np.arange(18.,32.,1.),latlon=True,cmap=cmo.cm.thermal)
cb1 = plt.colorbar(cf1,cax=cbaxes[0],orientation='horizontal')
cb1.set_label(r'Temperatura ($^o$C)',fontsize=18)
cr1 = m_axes[0].contour(lonClim,latClim,depth,latlon=True,colors=('k'),linestyles=('--'),linewidths=[.4,.4,.4],levels=[100,200,1000])
plt.clabel(cr1,fontsize=9,inline=1,fmt='%i')

cf2 = m_axes[1].contourf(lonClim,latClim,sstAnom,np.arange(-1.,6.,.1),latlon=True,cmap='seismic')
cb2 = plt.colorbar(cf2,cax=cbaxes[1],orientation='horizontal')
cb2.set_label(r'Anomalia de Temperatura ($^o$C)',fontsize=18)
cr2 = m_axes[1].contour(lonClim,latClim,depth,latlon=True,colors=('k'),linestyles=('--'),linewidths=[.4,.4,.4],levels=[100,200,1000])
plt.clabel(cr2,fontsize=9,inline=1,fmt='%i')

m_axes[0].ax.set_title('GHRSST Summer Mean',fontsize=24)
m_axes[1].ax.set_title('GHRSST Summer Mean Anomaly',fontsize=24)

##############################################
#     plotando anomalia em alguns dias       #
##############################################
