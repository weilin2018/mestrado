# localização: -23.8184; -45.40 -> grade PCSE [57,10]

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
import decomp

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
DATA_DIR = '/media/danilo/Danilo/mestrado/ventopcse/data/ECOSAN/EcoSanFundeios_Dados_Out2005Fev2006/ES12CanalSS/'
FILENAME1 = DATA_DIR + 'ES1201.Dat' # 5m depth
FILENAME2 = DATA_DIR + 'ES1202.Dat' # Xm depth

ncin = pd.read_csv(FILENAME1,skiprows=11,delimiter=',',usecols=[2,3,5,12],header=None)
ncin.columns = ['ASPD','AVDIR','datetime','STEMP']

ncin.index = pd.DatetimeIndex(ncin.datetime.values) # assigning datetimeindex
ncin.drop('datetime',axis=1,inplace=True) # removing column datetime

# we have a lot of values for each hour of each day, but it's not a pattern
# we can resample this timeseries is a hourly frequency
css5m = ncin.resample('H').mean()

# just a simple teste to visualize how things are going
fig,ax = plt.subplots()
css5m.STEMP.plot(ax=ax,c='#c0c0c0') # hourly frequency
css5m.resample('D').mean().STEMP.plot(ax=ax,c='k') # daily frequency
plt.title('5m depth temperature')

# converting intensity and direction, with magnetic declination correction and rotation
wu,wv = decomp.intdir2uv(css5m.ASPD.values,css5m.AVDIR.values,-21.09,0.)
d = pd.DataFrame({'wu':wu,'wv':wv},index=css5m.index)
d = d.resample('D').mean()

fig,ax = plt.subplots()
oceano.stickplot(d['2006'],ax)
