

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
import gsw

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

import plots as plotDanilo
##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
plt.ion()

# configurações do plot
figsize = (17.4/2.54, 10/2.54)

DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
# select which experiment you want to plot:
exp     = 'EA1.cdf'
fname   = DATA_DIR+exp
savefig = False

# velocidade
contours = np.arange(0,1.5,0.01)

timestep = input('Type which timestep to plot (type 999 to plot three images): ')

if timestep == 999.:
    timestep = [0,46,303]

    for nstep in timestep:
        fig,axes = plotDanilo.create_Structure_horizontal(fname,contours,property='salt',timestep=int(nstep),savefig=savefig)
    # plt.close()
else:
    fig,axes = plotDanilo.create_Structure_horizontal(fname,contours,property='speed',timestep=int(timestep),savefig=savefig)
