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
class Data(object):

    def __init__(self,fname,skiprows,usecols,dateIndex):
        self.fname = fname
        self.skiprows = skiprows
        self.usecols = usecols
        self.dateIndex = dateIndex

    def read_file(self):
        # le e armazena planilha em DataFrame
        self.sheet = pd.read_excel(self.fname,skiprows=self.skiprows,usecols=self.usecols,header=None)

    # def create_doubleHeader(self):
        # header = pd.read_excel(self.fname,skiprows=self.)

    def create_datetime(self,begin,final,freq):
        # dateBegin = self.sheet[0].values[0]
        # hourBegin = self.sheet.columns[1]
        # dateFinal = self.sheet[0].values[-1]
        # hourFinal = self.sheet.columns[-1]
        #
        # self.dates = [dateBegin,dateFinal]
        # self.hours = [hourBegin,hourFinal]

        # remover a linha e coluna
        # self.sheet = self.sheet.drop(self.sheet.columns[0],axis=1) # coluna
        self.dtRange = pd.date_range(start=begin,end=final,freq=freq)

    # def create_datetime4radiation(self,begin,final,hours):


    def create_timeseries(self):
        # ler velocidade do vento
        data = []

        for i,row in self.sheet.iterrows():
            for j in row.values:
                if float(j):
                    data.append(float(j))
                else:
                    data.append(np.nan)

        self.data = np.asarray(data)


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/data/INMET')

nfiles = glob.glob(DATA_DIR+'*.xls')

# reading first wind velocity [m/s] - colunas de 1 a 24
paraty_1 = Data(nfiles[0],11,np.arange(1,25),None)
paraty_1.read_file()
paraty_1.create_datetime(begin='2006-11-19 00:00',final='2014-12-31 23',freq='1H')
paraty_1.create_timeseries() # leitura dos dados de intensidade e organizacao do vetor

# armazenamento de dados de intensidade em dataframe
df_paraty = pd.DataFrame({'intensity':paraty_1.data},index=paraty_1.dtRange)

# reading now wind's direction, columns 25 to 48
paraty_1 = Data(nfiles[0],11,np.arange(25,49),None)
paraty_1.read_file()
paraty_1.create_timeseries() # leitura dos dados de intensidade e organizacao do vetor

# armazenando os dados de direcao
df_paraty['direction'] = paraty_1.data

# reading radiation data
paraty_1 = Data(nfiles[0],11,np.arange(49,64),None)
paraty_1.read_file()
paraty_1.create_timeseries() # leitura dos dados de intensidade e organizacao do vetor


# sheet = pd.read_excel(nfiles[0],skiprows=10,header=0,usecols=np.arange(0,25))

# cada coluna -> 1 horário
# cada linha  -> 1 dia

"""
Concatenar linha a linha uma na outra, criando um vetor de time (0 a 23) e,
ao mesmo tempo, criando um vetor de datas.

dica:

usar sheet.iterrows()

"""
