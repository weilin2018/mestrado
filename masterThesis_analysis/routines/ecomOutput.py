# using my package for model's visualization

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

def load_exp(fname,year='2010'):

    timeStart = '%s-01-15'%(year)
    timeEnd   = '%s-02-27'%(year)

    exp = Experiment(fname,timeStart=timeStart,timeEnd=timeEnd,region='pcse')

    return exp

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
# Type the name of the gcmplt you want to analyze, without the ext (e.g., exp06)
exp = 'exp11'
fname = '/media/danilo/Danilo/mestrado/ventopcse/output/%s.cdf'%(exp)

control = Experiment(fname,timeStart='2014-01-15',timeEnd='2014-02-27',region='pcse')

# run graph
plt.ion()

control.SBC() # define a location for timeseries
control.plotGraph(var='temp',sigma=-1)
control.graph(title='Sea Surface Temperature Timeseries in the Sao Sebastiao Channel')

# run field animation
control.plotAnim(var='temp',sigma=-1)
d = {'cmap':cmo.cm.thermal,'latlon':True}
control.field(title='Bottom Temperature',**d)

def exp06_v_control():
    """ testando a analise de anomalias """

    fname     = '/media/danilo/Danilo/mestrado/ventopcse/output/exp06.cdf'
    with_hflx = Experiment(fname,timeStart='2014-01-15',timeEnd='2014-02-14',region='pcse')

    fname = '/media/danilo/Danilo/mestrado/ventopcse/output/control_2010.cdf'
    control = Experiment(fname,timeStart='2010-01-15',timeEnd='2010-02-14',region='pcse')

    # creating boundaries
    with_hflx.Sp_coast()
    control.Sp_coast()

    # calculating thermal anomalies
    with_hflx.var = with_hflx.ncin.temp[control.timeStart[0]:control.timeEnd[0],0,:,:] - np.nanmean(with_hflx.ncin.temp[control.timeStart[0]:control.timeEnd[0],0,:,:])
    control.var   = control.ncin.temp[control.timeStart[0]:control.timeEnd[0],0,:,:] - np.nanmean(control.ncin.temp[control.timeStart[0]:control.timeEnd[0],0,:,:])

    with_hflx.var = with_hflx.var[-1,:,:]
    control.var   = control.var[-1,:,:]

    d = {'cmap':cmo.cm.thermal,'latlon':True}
    with_hflx.field(**d)
    control.field(**d)
