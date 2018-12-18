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

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
sim1 = 'EA1.cdf'
sim2 = 'EA1_dti10.cdf'

MODELO_DIR = '/media/danilo/Danilo/mestrado/ventopcse/output/'

k,i,j = 0,38,11

#### LOADING MODELLED DATA
modelado  = xr.open_dataset(MODELO_DIR+sim1)

u = modelado.u[:,k,i,j].values
v = modelado.v[:,k,i,j].values
e = modelado.elev[:,i,j].values
t = modelado.temp[:,k,i,j].values
s = modelado.salt[:,k,i,j].values
EA1 = pd.DataFrame({'u':u,'v':v,'elev':e,'temp':t,'salt':s},index=pd.DatetimeIndex(modelado.time.values))

modelado  = xr.open_dataset(MODELO_DIR+sim2)

u = modelado.u[:,k,i,j].values
v = modelado.v[:,k,i,j].values
e = modelado.elev[:,i,j].values
t = modelado.temp[:,k,i,j].values
s = modelado.salt[:,k,i,j].values
EA1_dti10 = pd.DataFrame({'u':u,'v':v,'elev':e,'temp':t,'salt':s},index=pd.DatetimeIndex(modelado.time.values))

# plotando
fig,ax = plt.subplots(nrows=5,sharex=True)

ax[0].set_title('u-comp')
EA1.u.plot(ax=ax[0],label='EA1')
EA1_dti10.u.plot(ax=ax[0],label='EA1')

ax[1].set_title('v-comp')
EA1.v.plot(ax=ax[1],label='EA1')
EA1_dti10.v.plot(ax=ax[1],label='EA1')

ax[2].set_title('eta')
EA1.elev.plot(ax=ax[2],label='EA1')
EA1_dti10.elev.plot(ax=ax[2],label='EA1')

ax[3].set_title('temp')
EA1.temp.plot(ax=ax[3],label='EA1')
EA1_dti10.temp.plot(ax=ax[3],label='EA1')

ax[4].set_title('temp')
EA1.salt.plot(ax=ax[4],label='EA1')
EA1_dti10.salt.plot(ax=ax[4],label='EA1')

plt.legend()

plt.show()
