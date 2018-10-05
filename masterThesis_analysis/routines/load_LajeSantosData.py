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
def readFiles_Laje(nfiles):
    '''
        read file by file, extracting each column and appending in
        lists.

        args
            nfiles  (list): list with all files

        returns
            a bunch of lists ...
    '''

    # columns:
        # 1 - ano
        # 2 - mes
        # 3 - dia
        # 4 - hora
        # 5 - lat (errado)
        # 6 - lon (errado)
        # 7 - pressao
        # 8 - intensidade do vento
        # 9 - direcao do vento
        #10 - temperatura
        #11 - umidade relativa
        #12 - precipitacao
        #13 - radiacao solar (nao confiar, pois os passaros fazem sombra sobre o sensor)

    datas = []
    inten = []    # intensity
    direc = []    # directions

    for f in nfiles:
        # read file
        data = np.loadtxt(f)

        for i in np.arange(0,data.shape[0],1):
            dt = convert2datetime(data[i,0], data[i,1], data[i,2], data[i,3])

            datas.append(dt)
            inten.append(data[i,7])
            direc.append(data[i,8])

    return np.asarray(datas),np.asarray(inten),np.asarray(direc)


# convert columns to datetime
def convert2datetime(year,month,day,hour):
    '''
        receive some year,months,days and hours
        convert to datetime
        return a valid datetime
    '''

    import datetime

    return datetime.datetime(*map(int, [year,month,day,hour]))

def loadData(DATA_DIR,pickle=False):

    import decomp

    if pickle:
        df = pickle.load(open(DATA_DIR,'r'))
    else:
        # select files
        nfiles = glob.glob(DATA_DIR + '*.txt')
        nfiles.sort()

        # extract information from files
        dates,intensity,direction = readFiles_Laje(nfiles)
        # correction of magnetic declination
        # ndirection = correctionMagDec(direction,-21.03)
        # here we need tranform to oceanographic convention, by adding 180 degrees
        # ndirection = [ang+180. for ang in ndirection]
        # ndirection = np.asarray(ndirection)
        # ndirection *= np.pi/180

        # convert intensity and direction to components (u,v)
        decliMag = -21.09
        angRot   = 54.

        # keeping arguments as arrays
        direction = np.asarray(direction)
        intensity = np.asarray(intensity)

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

    d = df.to_xarray()
    d.attrs = attrs

    d.to_netcdf(path+'lajeSantos_data.nc')

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
os.system('clear')

# define some constants
BASE_DIR = oceano.make_dir()

# se não houver um netcdf ainda, usar a linha abaixo
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/data/Est_lajeSantos/2015/atualizado/')


df = loadData(DATA_DIR)

saveNetCDF(df,'/media/danilo/Danilo/mestrado/ventopcse/data/')
