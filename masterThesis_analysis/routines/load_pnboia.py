# Funcao para ler os dados da Estacao da Laje de Santos e salvar em um netcdf

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

# read files from laje de santos
def replace(x):
    return x.replace(',','.')

# convert columns to datetime
def convert2datetime(year,month,day,hour):
    '''
        receive some year,months,days and hours
        convert to datetime
        return a valid datetime
    '''

    import datetime
    t = []

    for y,m,d,h in zip(year,month,day,hour):
        t.append(datetime.datetime(*map(int,[y,m,d,h])))

    return np.asarray(t)

def loadData(DATA_DIR,pickle=False):

    import decomp

    if pickle:
        df = pickle.load(open(DATA_DIR,'r'))
    else:

        data = pd.read_csv(DATA_DIR,header=0,usecols=[1,2,3,4,5,6,10,13],decimal=".")

        dates = convert2datetime(data['Year'].values,data['Month'].values,data['Day'],data['Hour'])

        # df.columns = ['Lat', 'Lon', 'Wspd', 'Wdir']

        # replacing comma by dots
        data['Lat'] = data['Lat'].apply(replace)
        data['Lon'] = data['Lon'].apply(replace)
        data['Wspd'] = data['Wspd'].apply(replace)
        # convert intensity and direction to components (u,v)
        decliMag = -21.09
        angRot   = 54.

        # keeping arguments as arrays
        direction = data['Wdir'].apply(float).values
        intensity = data['Wspd'].apply(float).values

        # usando rotina da Carine, em decomp.py
        along,across = decomp.intdir2uv(intensity,direction,decliMag,angRot)

        # create pandas.dataframe
        i = pd.DatetimeIndex(dates)
        df = pd.DataFrame({'wu':along, 'wv':across},
                            index=i)

        # quality control
        df[df['wu'] > df['wu'].std()*3] = np.nan
        df[df['wu'] < df['wu'].std()*-3] = np.nan

        df[df['wv'] > df['wv'].std()*3] = np.nan
        df[df['wv'] < df['wv'].std()*-3] = np.nan

    return df

def saveNetCDF(df,path,attrs=False):

    if not attrs:
        attrs = {
            'Quality Control': 'removed data out of the limits of 3*standrd deviation',
            'Rotation': 'Rotation of 54 degrees, based on the coastline',
            'Magnetic Declination': '-21.09º',
            'Longitude': '46.1803°W',
            'Latitude' : '24.3194°S',
            'Convention': 'Already converted from Meteorological do Oceanographic'
        }

        # atributos para as variáveis!!!!
    attrs_variables = {
        'wu': {
            'units': u'm s-1',
            'standard_name' : u'alongshore_wind',
        },
        'wv': {
            'units': u'm s-1',
            'standard_name': u'cross_shore_wind',
        }
    }

    d = df.to_xarray() # convertendo para xarray instance
    d.attrs = attrs # inserindo atributos gerais do netcdf

    # inserindo atributos por variáveis
    for var in attrs_variables.keys():
        d[var].attrs = attrs_variables[var]

    # criando netcdf
    d.to_netcdf(path+'pnboiaSantos_data.nc')

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
os.system('clear')

# define some constants
BASE_DIR = oceano.make_dir()

# se não houver um netcdf ainda, usar a linha abaixo
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/data/pnboia/dados_validadosMarinha/santos.csv')

df = loadData(DATA_DIR)

saveNetCDF(df,'/media/danilo/Danilo/mestrado/ventopcse/data/')
