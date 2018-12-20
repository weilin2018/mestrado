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
import xlrd
from dateutil import parser

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
class Inmet(object):

    def __init__(self,fname,skiprows):
        self.fname = fname
        self.skiprows = skiprows

    def read_file(self):
        # le e armazena planilha em DataFrame
        self.sheet = pd.read_excel(self.fname,skiprows=self.skiprows,header=None)

    def extract_data(self,intCols,dirCols,radCols):
        self.intensity = self.create_timeseries(self.sheet.iloc[:,intCols].copy())
        self.direction = self.create_timeseries(self.sheet.iloc[:,dirCols].copy())
        self.radiation = self.create_timeseries(self.sheet.iloc[:,radCols].copy())

    def create_wind_dataframe(self,begin,final,freq):
        self.create_datetime(begin,final,freq)

        self.wind = pd.DataFrame({'intensity':self.intensity,'direction':self.direction},index=self.dtRange)

    def create_radiation_dataframe(self,begin,final,freq):
        self.create_datetime(begin,final,freq)

        self.solarrad = pd.DataFrame({'radiation':self.radiation},index=self.dtRange)

    def create_datetime(self,begin,final,freq):
        self.dtRange = pd.date_range(start=begin,end=final,freq=freq)

    def create_timeseries(self,newdata):
        # ler velocidade do vento
        data = []

        for i,row in newdata.iterrows():
            for j in row.values:
                if float(j):
                    data.append(float(j))
                else:
                    data.append(np.nan)

        return np.asarray(data)

    def convert_units(self,unit='kJ'):

        # convert from kJ to W
        if unit == 'kJ':
            self.solarrad['radiation_W'] = self.solarrad.radiation / 3.6

        # convert from W to kJ
        if unit == 'W':
            self.solarrad['radiation_kJ'] = self.solarrad.radiation * 3.6


    def save2netcdf(self,outFile,attrs=None):
        final_df = pd.DataFrame({'wInt':self.intensity,'wDir':self.direction,'rad': self.radiation},index=self.dtRange)

        final_df = final_df.to_xarray()
        final_df.attrs = attrs

        final_df.to_netcdf(outFile)


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/data/INMET')

nfiles = glob.glob(DATA_DIR+'*.xls')


#############################################################################
################### leitura do 1o arquivo: de 2006 a 2014
#############################################################################
fname = nfiles[0]

intCols = np.arange(1,25)
dirCols = np.arange(25,49)
radCols = np.arange(49,73)

paraty_file1 = Inmet(nfiles[0],11) # create object with two property: filename and skiprows
paraty_file1.read_file()          # read spreadsheet, creating sheet method as DataFrame
paraty_file1.extract_data(intCols,dirCols,radCols)

paraty_file1.create_wind_dataframe(begin='2006-11-19 00:00',final='2014-12-31 23',freq='1H')
paraty_file1.create_radiation_dataframe(begin='2006-11-19 00:00',final='2014-12-31 23',freq='1H')

# save data in netcdf format
attrs = {
    'Period': '2006-11-19 to 2014-12-31',
    'Database': 'INMET/CPTEC',
    'Location': 'Paraty/RJ',
    'lat': '23°13\'S',
    'lon': '44°43\'W',
    'Convention': 'Meteorological Coordinate System',

}

paraty_file1.save2netcdf('/media/danilo/Danilo/mestrado/ventopcse/data/INMET/netcdf/Paraty_2006_2014.nc',attrs)

#############################################################################
######################## read file for 2015 data
#############################################################################
"""
Neste arquivo, a ordem das variaveis esta trocada, logo:
"""
fname = nfiles[1]
intCols = np.arange(1,25)
dirCols = np.arange(49,73)
radCols = np.arange(25,49)


paraty_file2 = Inmet(fname,11)
paraty_file2.read_file()
paraty_file2.extract_data(intCols,dirCols,radCols)

paraty_file2.create_wind_dataframe(begin='2015-01-01 00:00',final='2017-12-31 23',freq='1H')
paraty_file2.create_radiation_dataframe(begin='2015-01-01 00:00',final='2017-12-31 23',freq='1H')

# save data in netcdf format
attrs = {
    'Period': '2015-01-01 to 2017-12-31',
    'Database': 'INMET/CPTEC',
    'Location': 'Paraty/RJ',
    'lat': '23°13\'S',
    'lon': '44°43\'W',
    'Convention': 'Meteorological Coordinate System',

}

paraty_file2.save2netcdf('/media/danilo/Danilo/mestrado/ventopcse/data/INMET/netcdf/Paraty_2015_2017.nc',attrs)

#############################################################################
######################## MERGING DATAFRAMES
#############################################################################
df_tmp = paraty_file1.wind.copy()
df_tmp['radiation'] = paraty_file1.solarrad.values

df_tmp2 = paraty_file2.wind.copy()
df_tmp2['radiation'] = paraty_file2.solarrad.values

data = pd.concat([df_tmp,df_tmp2])

df = data.to_xarray()
df.attrs = attrs

df.direction.attrs = {
    'long_name': 'wind_direction',
    'unit':      'degrees',
    'convention': 'meteorological',
}

df.intensity.attrs = {
    'long_name': 'wind_intensity',
    'unit':      'm s-2',
}

df.radiation.attrs = {
    'long_name': 'incident_solar_radiation',
    'unit':      'kJ m-2',
}

df.to_netcdf('/media/danilo/Danilo/mestrado/ventopcse/data/INMET/netcdf/inmet_paraty_2006-2017.nc')

del df,df_tmp,df_tmp2
################################################################################

# visualizando dados para o verao de 2014
# calculando a media sem os pontos zeros
newdf = newp.solarrad.copy()

# newdf.iloc[np.where(df.radiacao == 0.)] = np.nan
media = newdf.mean()
std   = newdf.std()


# visualizando
fig,ax = plt.subplots(figsize=(15,5))

newp.solarrad.plot(ax=ax,title=u'Radiação Solar Registrada (kJ m-2), em vermelho, e Média Semanal, em azul,  entre 2013 e 2018 em Ubatuba.')

ax.axvline('2014-01-14 00:00',color='k')
ax.axvline('2014-02-15 00:00',color='k')

ax.axvline('2014-01-01 00:00',color='k',linestyle='--')
ax.axvline('2014-03-01 00:00',color='k',linestyle='--')

# ax.margins(0)

# plotando um resample a cada 15 dias
newdf.radiation.resample('1W').mean().plot(ax=ax)

plt.margins(0)
plt.tight_layout()
