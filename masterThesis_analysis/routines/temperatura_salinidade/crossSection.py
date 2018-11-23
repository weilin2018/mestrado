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

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def structurePlots():

    fig,ax = plt.subplots(nrows=3,ncols=2)

    

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

ec = xr.open_dataset('/media/danilo/Danilo/mestrado/ventopcse/output/EC1.cdf')
ea = xr.open_dataset('/media/danilo/Danilo/mestrado/ventopcse/output/EA1.cdf')
