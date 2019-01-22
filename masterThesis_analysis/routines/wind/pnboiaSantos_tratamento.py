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
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here,
def extractColumnsIndexByName(columns_names,filepath):
    '''
        Function to locate the number of columns by name, returning those indexes
        to read_csv extract only that we want.

        args
            columns_names    (array): array with the names
            filepath         (string): full path to the file we want to load

        returns
            indexes          (list): list with indexes related to columns_names
    '''

    # load all file
    data = pd.read_csv(filepath)
    # extract columns names to a array
    all_columns_names = data.columns.values

    # locate the position of each column passed as argument
    # [i for i,x in enumerate(testlist) if x == 1]
    indexes = []
    for column in columns_names:
        values = [i for i,x in enumerate(all_columns_names) if x == column]
        indexes.append(values[0])

    return indexes

def readPNBOIA(ARGO_DIR,region='Santos'):

    if region == 'Santos':
        nfiles = ARGO_DIR + 'Bsantos_argos.csv'

    # columns positions
    cols_name = ['Wspd', 'Wdir','Wspdflag']
    columns_data = extractColumnsIndexByName(cols_name,nfiles)
    cols_date = ['Year', 'Month', 'Day', 'Hour', 'Minute']
    columns_date = extractColumnsIndexByName(cols_date,nfiles)

    data = pd.read_csv(nfiles,usecols=columns_date)

    dates = []

    for index,row in data.iterrows():
        # dt = convert2datetime(data[i,0], data[i,1], data[i,2], data[i,3])
        values = [row[0], row[1], row[2], row[3], row[4]]
        dt = datetime.datetime(*map(int, values))
        dates.append(dt)

    # convert dates to pd.DateTimeIndex
    i = pd.DatetimeIndex(dates)

    # read values, using dates as index
    data = pd.read_csv(nfiles,usecols=columns_data)
    data['datetime'] = i
    data.set_index('datetime',inplace=True)

    # basic quality control
    # using the flag to insert NaN values
    data[data['Wspdflag'] == 4] = np.nan

    # now we can remove the flag column
    data.drop(columns=['Wspdflag'], inplace=True)

    # replace comma by dots
    Wspd = replaceComma4dots(data['Wspd'])
    Wdir = replaceComma4dots(data['Wdir'])

    return data.index.values,Wspd,Wdir


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
ARGO_DIR = BASE_DIR.replace('github', 'ventopcse/data/pnboia')

# lendo o arquivo
datas,ints,dirs = readPNBOIA(ARGO_DIR,region='Santos')

# definicao dos angulos para correcao da declinacao magnetica e rotacao dos eixos
magDecli = -21.08
rotAngle = -36

# corrigindo declinacao e rotacionando o eixo cartesiano
# para alinhamento com a maxima variavencia (Mazzini, 2009)
wur,wvr = dp.intdir2uv(ints,dirs,magDecli,rotAngle)

# criando dataframe
dfPNBOIA = pd.DataFrame({'wind_along':wvr,'wind_cross':wur},index=pd.DatetimeIndex(datas))

dirSave = ARGO_DIR + 'pnboiaSantos.nc'

saveNetCDF(dfPNBOIA,dirSave)
