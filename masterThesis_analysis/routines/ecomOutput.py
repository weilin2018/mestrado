# using my package of model's visualization

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
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
# Type the name of the gcmplt you want to analyze, without the ext (e.g., exp06)
exp = 'control_2010'
fname = '/media/danilo/Danilo/mestrado/ventopcse/output/%s.cdf'%(exp)

control = Experiment(fname,timeStart='2010-01-15',timeEnd='2010-02-14',region='pcse')

# run graph
control.SBC()
control.plotGraph(var='temp')
control.graph()

# run field animation
control.plotAnim(var='temp',sigma=0)
d = {'cmap':cmo.cm.thermal,'latlon':True}
# control.field(**d)
