# colocar o código principal aqui: instanciar classe com Experiment(), gerar ncin, importar variaveis

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

# importing package
from modelVisualization.animation import Animation

##############################################################################
#                          [GEN] CLASSES                                     #
##############################################################################

class Experiment(Animation):

    def __init__(self,fname,timeStart,timeEnd):
        """
        Parameters
        ----------
        fname : string
            Full path to the file.
        timeStart : integer
            Index to begin the cut in the time axis of the data.
        timeEnd : type
            Index to end the cut in the time axis of the data.
        """
        self.fname = fname
        self.ncin  = xr.open_dataset(fname)
        self.timeStart = timeStart
        self.timeEnd   = timeEnd
        self.extractVariables()

    def extractVariables(self):
        # importing latitude and longitude, masking zero values with NaN
        lon = self.ncin.lon.values
        lat = self.ncin.lat.values
        lon[lon == 0.] = np.nan
        lat[lat == 0.] = np.nan

        self.lon = lon
        self.lat = lat

    def view_grid(self,figsize):
        fig,ax = plt.subplots(figsize=figsize)
        m = oceano.make_map(ax)
        x,y = m(self.lon,self.lat)
        m.plot(x,y,'k',alpha=.3)
        m.plot(x.T,y.T,'k',alpha=.3)
        plt.show()

    def anim(self,var='elev'):
        Animation.__init__(self,var=var)

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
# Type the name of the gcmplt you want to analyze, without the ext (e.g., exp06)
exp = 'control_2010'
fname = '/media/danilo/Danilo/mestrado/ventopcse/output/%s.cdf'%(exp)

control = Experiment(fname,timeStart=112,timeEnd=352)
control.anim()

d = {'cmap':'RdBu_r','latlon':True}
control.field(**d)
