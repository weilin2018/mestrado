# reading data from Araca Bay, sended by Marcelo Dottori

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
import seawater as sw

import matplotlib
# matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

def create_dates(fname):

    ncin = np.loadtxt(fname,usecols=[0,1,2,3,4])

    dates = []
    years = ncin[:,0]
    month = ncin[:,1]
    days  = ncin[:,2]
    hours = ncin[:,3]
    minut = ncin[:,4]

    for y,m,d,h,mi in zip(years,month,days,hours,minut):
        dates.append(datetime.datetime(int(y),int(m),int(d),int(h),int(mi)))

    return np.asarray(dates)

##
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/data/Araca/Dottori_etal_2015')

names = ['surfSalt','surfSalt_flag',
		 'surfTemp','surfTemp_flag',
		 'bottSalt','bottSalt_flag',
		 'bottTemp','bottTemp_flag']

fname = DATA_DIR + 'ct_araca01.dat'

# creating datetimeindex
dt = create_dates(fname)

# reading data with appropriate index
df = pd.read_csv(fname, sep=' ',header=None,usecols=[7,8,9,10,11,12,13,14],names=names)
df.index = pd.DatetimeIndex(dt)

##
# Para facilitar o tratamento de cada variável, irei destrinchar em 4 dataframes [variavel, flag]
##
saltSurf = pd.read_csv(fname,sep=' ',header=None,usecols=[7,8],names=['saltSurf','flag'])
tempSurf = pd.read_csv(fname,sep=' ',header=None,usecols=[9,10],names=['tempSurf','flag'])
saltBott = pd.read_csv(fname,sep=' ',header=None,usecols=[11,12],names=['saltBott','flag'])
tempBott = pd.read_csv(fname,sep=' ',header=None,usecols=[13,14],names=['tempBott','flag'])

# quality control based on flags
#df.loc[df['column_name'] == some_value]
#saltSurf.loc[saltSurf['flag'] == 3] = np.nan
#tempSurf.loc[tempSurf['flag'] == 3] = np.nan
#saltBott.loc[saltBott['flag'] == 3] = np.nan
#tempBott.loc[tempBott['flag'] == 3] = np.nan

# creating only one dataframe with all values
df = pd.DataFrame({
    'tempSurf': tempSurf.tempSurf.values,
    'saltSurf': saltSurf.saltSurf.values,
    'tempBott': tempBott.tempBott.values,
    'saltBott': saltBott.saltBott.values
    },index=pd.DatetimeIndex(dt))

# plotting only 2014 data
df_2014 = df['2014'].copy()

# save as netcdf
d = df_2014.to_xarray()
d.to_netcdf(DATA_DIR+'dados_tratados.nc')

# resampling for daily measurements
df_2014_daily = df_2014.resample('D').mean()


# plotar
fig,ax = plt.subplots(nrows=2)

# salinity
df_2014['saltSurf'].plot(ax=ax[0])
df_2014['saltBott'].plot(ax=ax[0])

ax[0].legend(['Surface','Bottom'],loc='best')

df_2014['tempSurf'].plot(ax=ax[1])
df_2014['tempBott'].plot(ax=ax[1])

ax[1].legend(['Surface','Bottom'],loc='best')

savefig_dir = oceano.make_dir()
#plt.savefig(savefig_dir + 'masterThesis_analysis/figures/dados_observados/araca_TempSalt.png',dpi=300)

plt.show()
