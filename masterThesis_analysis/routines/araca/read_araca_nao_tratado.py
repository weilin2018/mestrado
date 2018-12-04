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

# superficie e fundo estao em arquivos diferentes

#############################################################################
#                               SUPERFICIE                                  #
#############################################################################
fname = DATA_DIR + 'TS_fev2014_marco2014_superficie.dat'

# criando datas
dates = create_dates(fname)

# criando dataframe
df_superf = pd.read_csv(fname, sep=' ',header=None,usecols=[6,7],names=['Temperatura','Salinidade'])
df_superf.index = pd.DatetimeIndex(dates)


#############################################################################
#                               FUNDO                                  #
#############################################################################
fname = DATA_DIR + 'TS_fev2014_marco2014_fundo.dat'

# criando datas
dates = create_dates(fname)

# criando dataframe
df_fundo = pd.read_csv(fname, sep=' ',header=None,usecols=[6,7],names=['Temperatura','Salinidade'])
df_fundo.index = pd.DatetimeIndex(dates)

# merging
dct = {
    'saltSurf': df_superf.Salinidade.values,
    'tempSurf': df_superf.Temperatura.values,
    'saltBott': df_fundo.Salinidade.values,
    'tempBott': df_fundo.Temperatura.values
}

df = pd.DataFrame(dct,index=pd.DatetimeIndex(dates))

# save as netcdf
d = df.to_xarray()
d.to_netcdf(DATA_DIR+'dados_nao_tratados.nc')
