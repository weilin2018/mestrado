# extracao de dados do ghrsst proximo ao ponto da baia do araca

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

from dateutil import parser

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def converting_datetime(dates,hours):
	date = []

	for d,h in zip(dates,hours):
		if h == 0:
			date.append(str(d))
		else:
			date.append(str(d)+' '+str(h))

	return pd.to_datetime(date)



##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/data/IOUSP/labdados')

fname = DATA_DIR+'u/um2014t1.csv'

# importa dados, pulado o cabecalho do arquivo e usando a linha 0 (17 no arquifo) como header
# ncin = pd.read_csv(fname,sep=';',header=0,usecols=[0,2],parse_dates=True,infer_datetime_format=True)
ncin = pd.read_csv(fname,sep=';',header=0,usecols=[0,2])
ncin.columns = ['datetime','radiacao']
ncin.index = pd.DatetimeIndex(ncin.datetime)

# replacing comma by dots
ncin.radiacao = ncin.radiacao.apply(lambda x : str(x).replace(',','.'))

# converting radiacao data from object to float
ncin.radiacao = ncin.radiacao.apply(float)

ncin.plot(title='Radiacao em Ubatuba - Verao 2014 [W/m2]')
