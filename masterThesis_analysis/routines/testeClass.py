# add some description here

import glob
import matplotlib.pyplot as plt
import numpy as np
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
        self.ncin      = xr.open_dataset(fname) # define xarray instance
        self.timeStart = timeStart              # define timeStart as attribute
        self.timeEnd   = timeEnd                # define timeEnd as attribute
        self.fname     = fname                  # define fname as attribute


    def treatCoordinates(self,lon,lat):
        """Remove 0.0 values of longitude and latitude, replacing by NaN.

        Parameters
        ----------
        lon : np.ndarray
            Longitude array.
        lat : np.ndarray
            Latitude array.

        Returns
        -------
        lon,lat : np.ndarray
            Longitude and Latitude arrays with zeros replaced by nan.

        """
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
        elev = self.ncin.elev.values[self.timeStart:self.timeEnd,j,i]# add some description here
        time = self.time
        # plot
        plt.plot(time,elev,label=self.fname[-9:-4])

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

        # remove variables from the memory
        del spd,contour_levels,wu,wv


    def wind_elev(self,xStep=4,yStep=4,beginPlot=0,endPlot=-1):

        # calculate speed wind to normalize arrows, based in all data for wind
        spd_norm = np.sqrt(self.ncin.wu[self.timeStart:self.timeEnd,:,:].values**2 + self.ncin.wv[self.timeStart:self.timeEnd,:,:].values**2)

        # define contour_levels for each contourf
        clevs_elev = np.arange(-0.3,0.3,0.01)
        clevs_wind = np.arange(-15,15,0.01)

        # start animation
        plt.ion()
        fig, ax = plt.subplots(ncols=2)
        # divider = make_axes_locatable(ax[0])
        # c_wind= divider.append_axes("right", size="5%",pad=0.05)
        # c_wind.set_ylim([-0.4,0.4])
        #
        # divider = make_axes_locatable(ax[1])
        # c_elev = divider.append_axes("right", size="5%", pad=0.05)
        # c_elev.set_ylim([-0.4,0.4])


        # checking if the object had the attributes needed
        if hasattr(self, 'time') != True:
            # if the object doesn't:
            self.importBasicInformations()
        else:
            for i in range(self.time[beginPlot:endPlot].shape[0]):
            # for i in range(10):

                # cleaning screen
                ax[0].clear()
                ax[1].clear()

                # printing time
                plt.suptitle(self.time[i])

                # defining basemap instance for each subplot
                m_wind = oceano.make_map(ax[0],resolution='i')
                m_elev = oceano.make_map(ax[1],resolution='i')

                # importing data in each time i
                wu  = self.ncin.wu.values[i,::xStep,::yStep]
                wv  = self.ncin.wv.values[i,::xStep,::yStep]
                spd = np.sqrt(wu**2 + wv**2)
                elev= self.ncin.elev.values[i,:,:]

                # normalizing vectors with the total speed
                wun = wu/spd_norm[i,::xStep,::yStep]
                wvn = wv/spd_norm[i,::xStep,::yStep]

                # plot wind velocity information in axis 0
                cw = m_wind.contourf(self.lon[::xStep,::yStep],self.lat[::xStep,::yStep],spd,clevs_wind,latlon=True)
                qw = m_wind.quiver(self.lon[::xStep,::yStep],self.lat[::xStep,::yStep],wun,wvn,latlon=True)

                # plot elevation information in axis 1
                ce = m_elev.contourf(self.lon,self.lat,elev,clevs_elev,latlon=True,cmap="RdBu_r")

                # if i == 0:
                #     cbar_elev = plt.colorbar(ce,orientation='vertical')
                #     cbar_wind = plt.colorbar(cw,orientation='vertical')

                plt.pause(0.1)




        # criar

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
fname = '/media/danilo/Danilo/mestrado/ventopcse/output/exp06.cdf'

# instanciate an object, passing as argument, the full name
exp = Experiment(fname,timeStart=112,timeEnd=352)
# import basic informations
exp.importBasicInformations()
