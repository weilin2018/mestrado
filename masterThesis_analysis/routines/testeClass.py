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
class Experiment(object):
    def __init__(self,fname,timeStart=0,timeEnd=-1):
        self.ncin      = xr.open_dataset(fname)
        self.timeStart = timeStart
        self.timeEnd   = timeEnd


    def treatCoordinates(self,lon,lat):
        """ turn all 0.0 points in nan values """
        lon[lon == 0.] = np.nan
        lat[lat == 0.] = np.nan

        return lon,lat

    def importBasicInformations(self):
        """
            import from ncin all basic informations we need from
            the experiment, like longitude, latitude and time
        """
        # load coordinates
        lon = self.ncin.lon.values
        lat = self.ncin.lat.values
        self.time= self.ncin.time.values[self.timeStart:self.timeEnd]

        # treat coordinates
        self.lon,self.lat = self.treatCoordinates(lon,lat)


    def elevationPlot(self,i,j):
        """Simple test to plot an elevatino timeseries, based on i,j given location.

        Parameters
        ----------
        i : integer
            Description of parameter `i`.
        j : integer
            Description of parameter `j`.
        """
        # extract elevation data
        elev = self.ncin.elev.values[:,j,i]# add some description here

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
class Experiment(object):
    def __init__(self,fname,timeStart=0,timeEnd=-1):
        self.ncin      = xr.open_dataset(fname)
        self.timeStart = timeStart
        self.timeEnd   = timeEnd


    def treatCoordinates(self,lon,lat):
        """ turn all 0.0 points in nan values """
        lon[lon == 0.] = np.nan
        lat[lat == 0.] = np.nan

        return lon,lat

    def importBasicInformations(self):
        """
            import from ncin all basic informations we need from
            the experiment, like longitude, latitude and time
        """
        # load coordinates
        lon = self.ncin.lon.values
        lat = self.ncin.lat.values
        self.time= self.ncin.time.values[self.timeStart:self.timeEnd]

        # treat coordinates
        self.lon,self.lat = self.treatCoordinates(lon,lat)


    def elevationPlot(self,i,j):
        """Simple test to plot an elevatino timeseries, based on i,j given location.

        Parameters
        ----------
        i : integer
            Description of parameter `i`.
        j : integer
            Description of parameter `j`.
        """
        # extract elevation data
        elev = self.ncin.elev.values[:,j,i]
        # plot
        plt.plot(elev)
        plt.show()

    def windField(self):

        fig,ax = plt.subplots()

        spd = np.sqrt(self.ncin.wu.values[self.timeStart:self.timeEnd,:,:]**2 + self.ncin.wv.values[self.timeStart:self.timeEnd,:,:]**2)

        for i in range(self.time.shape[0]):
            ax.clear()
            m = oceano.make_map(ax,resolution='i')

            wu = self.ncin.wu.values[i,:,:]/spd
            wv = self.ncin.wv.values[i,:,:]/spd

            m.quiver(self.lon,self.lat,wu,wv,latlon=True)

            plt.pause(0.1)


        # criar

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

        # plot
        plt.plot(elev)
        plt.show()

    def windField(self,xStep=4,yStep=4):

        fig,ax = plt.subplots()

        spd = np.sqrt(self.ncin.wu.values[self.timeStart:self.timeEnd,:,:]**2 + self.ncin.wv.values[self.timeStart:self.timeEnd,:,:]**2)

        contour_levels = np.arange(-1.,1.,0.3/500)

        for i in range(self.time.shape[0]):
            ax.clear()
            m = oceano.make_map(ax,resolution='i')

            wu = self.ncin.wu.values[i,:,:]/spd[i,:,:]
            wv = self.ncin.wv.values[i,:,:]/spd[i,:,:]

            m.contourf(self.lon,self.lat,self.ncin.elev[i,:,:],contour_levels,latlon=True,cmap='coolwarm')
            m.quiver(self.lon[::xStep,::yStep],self.lat[::xStep,::yStep],wu[::xStep,::yStep],wv[::xStep,::yStep],latlon=True,scale=25)

            plt.pause(0.1)


        # criar

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
fname = '/media/danilo/Danilo/mestrado/ventopcse/output/exp06.cdf'

# instanciate an object, passing as argument, the full name
exp = Experiment(fname)
# import basic informations
exp.importBasicInformations()

# plot wind field animation
exp.windField()
