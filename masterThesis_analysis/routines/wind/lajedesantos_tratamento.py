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


import decomp as dp
##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here

# convert columns to datetime
def convert2datetime(year,month,day,hour):
    '''
        receive some year,months,days and hours
        convert to datetime
        return a valid datetime
    '''

    import datetime

    return datetime.datetime(*map(int, [year,month,day,hour]))

# reading Laje de Santos files
def readFiles_Laje(nfiles):
    '''
        read file by file, extracting each column and appending in
        lists.

        args
            nfiles  (list): list with all files

        returns
            a bunch of lists ...
    '''

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

# function to facilitate the creation of netcdf files
def saveNetCDF(df,dir,attrs=None):
    df = df.to_xarray()
    if attrs:
        df.attrs = attrs

    df.to_netcdf(dir)

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/data/Est_lajeSantos/2015/atualizado/')

# lendo dados da Laje de Santos
nfiles = glob.glob(DATA_DIR + "*.txt")
nfiles.sort()
datas,ints,dirs = readFiles_Laje(nfiles)

# definicao dos angulos para correcao da declinacao magnetica e rotacao dos eixos
magDecli = -21.08
rotAngle = -36

# corrigindo declinacao e rotacionando o eixo cartesiano
# para alinhamento com a maxima variavencia (Mazzini, 2009)
wu,wv   = dp.intdir2uv(ints,dirs,0,0)
wv     *= -1
int,dir = dp.uv2intdir(wu,wv,0,0)
wur,wvr = dp.intdir2uv(int,dir,magDecli,-rotAngle)

# criando dataframe
dfLaje = pd.DataFrame({'wind_along':wvr,'wind_cross':wur},index=pd.DatetimeIndex(datas))

dirSave = BASE_DIR.replace('github','ventopcse/data/Est_lajeSantos') + 'lajesantos.nc'

saveNetCDF(dfLaje,dirSave)
