'''

do jeito que está, o código só calcula a EOF para um mês.

eu preciso elaborar uma rotina para ler os dados e ir tomando médias
mensais e alocando como numa matriz multidimensional, onde o eixo 0 é
o time.

aí eu posso falar de fazer uma EOF para obter uma climatologia.

mas para começar a brincar tá legal!!!

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
from matplotlib import dates
import datetime

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano


from eofs.standard import Eof
from eofs.examples import example_data_path


BASE_DIR = oceano.make_dir()

DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/data/CFSR/1992_2011/')

# extrair longitude e latitude
nfiles = glob.glob(DATA_DIR+"*.nc")
nfiles.sort()
fname = nfiles[0]
ncdata = xr.open_dataset(fname)
lon    = ncdata['lon'].values# - 360
lat    = ncdata['lat'].values

lon,lat = np.meshgrid(lon,lat)

wu = ncdata['U_GRD_L103']


# Create an EOF solver to do the EOF analysis. Square-root of cosine of
# latitude weights are applied before the computation of EOFs.
coslat = np.cos(np.deg2rad(lat))
wgts = np.sqrt(coslat)[..., np.newaxis]
solver = Eof(wu, weights=wgts)

# Retrieve the leading EOF, expressed as the correlation between the leading
# PC time series and the input SST anomalies at each grid point, and the
# leading PC time series itself.
eof1 = solver.eofsAsCorrelation(neofs=1)
pc1 = solver.pcs(npcs=1, pcscaling=1)
