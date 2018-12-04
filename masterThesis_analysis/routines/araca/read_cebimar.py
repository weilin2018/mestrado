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

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/data/Araca/Aurea_CEBIMAR')

fname = DATA_DIR + 'casts_bined_4Dott.txt'

ncin = pd.read_csv(fname,sep='	')
ncin.drop(['Type','Station'],inplace=True,axis=1) # removendo algumas colunas
ncin.set_index('Cruise',inplace=True) # definindo index

# pensar numa forma de agrupar pelo Cruise
