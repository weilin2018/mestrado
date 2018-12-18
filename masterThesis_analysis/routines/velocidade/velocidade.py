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
from modelVisualization.interface import Experiment


##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/output')

experimento = 'EA1.cdf'
fname = DATA_DIR + experimento
exp = Experiment(fname,timeStart='2014-01-15',timeEnd='2014-02-15',region='pcse')

exp.pcse()

# plotando animacao de velocidade
exp.plotAnim()
exp.velocidade(index_file='/media/danilo/Danilo/mestrado/github/masterThesis_analysis/routines/index_list.npy',sigma=0)
