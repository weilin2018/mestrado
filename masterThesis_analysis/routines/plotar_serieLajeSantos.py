'''

df = pa.read_csv(path,
                 names = ['date', 'time', 'open', 'high', 'low', 'close', 'vol'],
                 parse_dates={'datetime':['date','time']},
                 keep_date_col = True,
                 index_col='datetime'
             )

criar datetimeindex usando as proprias colunas
df['datetime'] = df.apply(lambda row: datetime.datetime.strptime(row['date']+ ':' + row['time'], '%Y.%m.%d:%H:%M'), axis=1)

===

df = pd.read_csv(f,
                 names=['year', 'month', 'day','hour','x','y','p', 'intensity', 'direction','temp','i','-','-'],
                 parse_dates={'datetime':['year','month','day','hour']},
                 keep_date_col=True,
                 index_col='datetime'
        )

'''

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

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

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

def convert2datetime(year,month,day,hour):
    '''
        receive some year,months,days and hours
        convert to datetime
        return a valid datetime
    '''

    import datetime

    return datetime.datetime(*map(int, [year,month,day,hour]))

BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/', '/ventopcse/data/Est_lajeSantos/2015/')

# select files
nfiles = glob.glob(DATA_DIR + '*.txt')
nfiles.sort()

# extract information from files
dates,intensity,direction = readFiles_Laje(nfiles)

# convert intensity and direction to components (u,v)
wu,wv = oceano.spdir2uv(intensity, direction, deg=True)

# create pandas.dataframe
i = pd.DatetimeIndex(dates)
df = pd.DataFrame({'wu':wu, 'wv':wv},
                    index=i)

# quality control
df[df['wu'] > 10] = np.nan
df[df['wu'] <-10] = np.nan

df[df['wv'] > 10] = np.nan
df[df['wv'] <-10] = np.nan

df.plot(subplots=True, title='Laje de Santos Meteorological Station')
plt.show()

'''
Como os dados do NCEP são de 6 em 6 horas, a partir da 1h inicial, precisamos
redimensionar os dados para essa frequencia
'''
dfResampled = df.resample('6H')

dfResampled.plot(subplots=True)
plt.show()
