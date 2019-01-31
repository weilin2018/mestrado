#-*-coding;utf-8-*-
"""

"""
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

# pacotes para minimap
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

import plots as plotDanilo

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
plt.ion()

# configurações do plot
figsize = (17.4/2.54, 10/2.54)

DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
fname = glob.glob(DATA_DIR+"*.cdf")

# select which experiment you want to plot:
exp = 'EC1.cdf'
savefig = True

for f in fname:
    if exp in f:
        experiment = f

fname = experiment

# temperatura
contours = np.arange(14,36,0.1)


timestep = input('Type which timestep to plot: ')

if timestep == 999.:
    timestep = [0,46,303]

    for nstep in timestep:
        fig,axes = plotDanilo.create_Structure_horizontal(fname,contours,property='temp',timestep=int(nstep),savefig=savefig)
    # plt.close()
else:
    fig,axes = plotDanilo.create_Structure_horizontal(fname,contours,property='temp',timestep=int(timestep),savefig=savefig)
