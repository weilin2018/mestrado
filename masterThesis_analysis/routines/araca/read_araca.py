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

##
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/Araca/Dottori_etal_2015')

names = ['year','month','day','hour','minute','seconds',
		 'surfSalt','surfSalt_flag',
		 'surfTemp','surfTemp_flag',
		 'bottSalt','bottSalt_flag',
		 'bottTemp','bottTemp_flag']

fname = '/home/tparente/danilo/mestrado/ventopcse/data/Araca/Dottori_etal_2015/ct_araca01.dat'
df = pd.read_csv(fname, sep=' ',header=None,names=names)

