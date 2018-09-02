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
from modelVisualization.places import mapa

##############################################################################
#                          [GEN] CLASSES                                     #
##############################################################################

class Experiment(Animation,mapa):

    def __init__(self,fname,timeStart,timeEnd,region='pcse'):
        """
        Parameters
        ----------
        fname : string
            Full path to the file.
        timeStart : string
            Index to begin the cut in the time axis of the data.
        timeEnd : string
            Index to end the cut in the time axis of the data.
        """
        self.fname = fname
        self.ncin  = xr.open_dataset(fname)
        self.region = region
        self.timeStart = self.findIndex_datetime(timeStart).values
        self.timeEnd   = self.findIndex_datetime(timeEnd).values
        self.extractVariables()
        self.definingRegionParameters()

    def extractVariables(self):
        # importing latitude and longitude, masking zero values with NaN
        lon = self.ncin.lon.values
        lat = self.ncin.lat.values
        lon[lon == 0.] = np.nan
        lat[lat == 0.] = np.nan

        self.lon = lon
        self.lat = lat

    def definingRegionParameters(self):
        # baseado na regiao instanciada, chama a funcao adequada em places.py
        # para setar os parametros para o basemap, em caso de plotagem de mapas
        placesAvailable = {
            'pcse': self.pcse,
            'sbc': self.Canal_ssb
        }

        placesAvailable[self.region]()


    def view_grid(self,figsize):
        fig,ax = plt.subplots(figsize=figsize)
        m = oceano.make_map(ax)
        x,y = m(self.lon,self.lat)
        m.plot(x,y,'k',alpha=.3)
        m.plot(x.T,y.T,'k',alpha=.3)
        plt.show()

    def findIndex_datetime(self,date):
        """Find index of a reference date in an array of datetimes from the model.

        Parameters
        ----------
        date : string
            Reference datetime to find, such as 2010-02-15.

        Returns
        -------
        d    : pandas.core.series.Series
            Closest date found, where d.name represent a timestamp instance
        and d.values are the index to be used.
        """
        try:
            from dateutil import parser
        except:
            print("Please install python-dateutils package to use this function")
            pass

        refDate = parser.parse(str(date)) # converting string into datetime instance
        # converting array of datetimes, from the model, into dataframe
        times = self.ncin.time.values
        df = pd.DataFrame(np.arange(0,len(times)),index=times)

        # performing search
        d = df.iloc[df.index.get_loc(refDate,method='nearest')]

        return d

    def anim(self,var='elev',sigma=0):
        Animation.__init__(self,var=var,sigma=sigma)

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
# Type the name of the gcmplt you want to analyze, without the ext (e.g., exp06)
exp = 'control_2010'
fname = '/media/danilo/Danilo/mestrado/ventopcse/output/%s.cdf'%(exp)

control = Experiment(fname,timeStart='2010-01-15',timeEnd='2010-02-14',region='pcse')
control.anim(var='temp',sigma=0)

d = {'cmap':cmo.cm.thermal,'latlon':True}

# d = {'cmap':'RdBu_r','latlon':True}
control.field(**d)
