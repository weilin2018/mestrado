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


os.system('clear')
##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
def nearest_date(items,pivot):
    nearest=min(items, key=lambda x: abs(x - pivot))
    timedelta = abs(nearest - pivot)
    return nearest, timedelta

def locate_closest_date(time,day,month,year):

    # convert into datetimeindex
    dates = pd.DatetimeIndex(time)

    pivot = dates[0] + pd.DateOffset(days=day,month=month,year=year)

    # locate the closest datetime index
    nearestDate = nearest_date(dates,pivot)
    index = np.where(dates == nearestDate[0])

    return index

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
        self.fname     = fname                  # define fname as attribute

    def set_indexes4time(self,timeStart,timeEnd):
        self.timeStart = timeStart              # define timeStart as attribute
        self.timeEnd   = timeEnd                # define timeEnd as attribute

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

    def elevationPlot(self,i,j,points=None):
        """Simple test to plot an elevatino timeseries, based on i,j given location.

        Parameters
        ----------
        i : integer
            Description of parameter `i`.
        j : integer
            Description of parameter `j`.
        points : list
            Title locations for multiples plots.

        Example
        -------
        >> exp = Experiment('/media/danilo/Danilo/mestrado/ventopcse/output/gcmplt.cdf')
        >> exp.importBasicInformations()
        >> iss = [50,7,22]
        >> jss = [29,55,94]
        >> exp.elevationPlot(iss,jss,points=['Laje de Santos','CSB','Ubatuba'])
        """

        if type(i) == int: # plot just one location
            # extract elevation data
            elev = self.ncin.elev.values[self.timeStart:self.timeEnd,j,i]# add some description here
            time = self.time
            # plot
            plt.plot(time,elev,label=points)

            plt.show()
        else:
            # plot multiples location
            fig,axes = plt.subplots(nrows=len(i))
            for loc in range(len(i)):
                elev = exp.ncin.elev.values[exp.timeStart:exp.timeEnd,j[loc],i[loc]]
                axes[loc].plot(exp.time,elev,label=points[loc])
                plt.legend()

            plt.show()

    def plot_windField(self,xStep=4,yStep=4):

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

    def plot_windVectors_with_elevation(self,t=None,xStep=4,yStep=4,beginPlot=0,endPlot=-1,offset=0):
        """plot wind vectors over elevation field.

        Parameters
        ----------
        t : integer
            Index for some specific timestep to plot. If None, then will be plotted
            an animation.
        xStep : integer
            How much index to skip when plotting.
        yStep : integer
            How much index to skip when plotting.
        beginPlot : integer
            Index to begin the plot. Default is 0.
        endPlot : integer
            Index to end the plot. Default is the last index -1.
        offset : float
            Some value to add ou sub in contour levels of elevation.Could be
            positive, to add, or negative, to subtract.
        """

        # contour levels for elevation
        minValue = np.nanmin(self.ncin.elev.values[beginPlot:endPlot,:,:])+offset
        maxValue = np.nanmax(self.ncin.elev.values[beginPlot:endPlot,:,:])+offset
        clevs_elev = np.arange(-.3,0.3,0.001)

        # based on xStep and yStep, create a scheme to skip vectors
        # skipwind = (slice(None,None,xStep),slice(None,None,xStep))

        plt.ion()
        fig,ax = plt.subplots()

        if hasattr(self,'time') != True:
            self.importBasicInformations()
        else:
            if t:
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%",pad=0.05)
                cax.set_ylim([minValue,maxValue])

                m = oceano.make_map(ax,resolution='i')

                wu   = self.ncin.wu.values[t,::xStep,::yStep]
                wv   = self.ncin.wv.values[t,::xStep,::yStep]
                elev = self.ncin.elev.values[t,:,:]

                # plot elevation as contourf
                ce = m.contourf(self.lon,self.lat,elev,clevs_elev,latlon=True,cmap="RdBu_r")
                cbar = plt.colorbar(ce,cax=cax,orientation='vertical')
                cbar.set_label('Elevation [m]')

                # plot wind data as quiver
                qw = m.quiver(self.lon[::xStep,::yStep],self.lat[::xStep,::yStep],wu,wv,latlon=True,alpha=.3,scale=150,width=0.005,pivot='middle')
                ax.set_title(str(self.time[t]),fontsize=24)
            else:
                for i in range(self.time[beginPlot:endPlot].shape[0]):
                    ax.clear()
                    m = oceano.make_map(ax,resolution='i')

                    wu   = self.ncin.wu.values[i,::xStep,::yStep]
                    wv   = self.ncin.wv.values[i,::xStep,::yStep]
                    elev = self.ncin.elev.values[i,:,:]

                    # plot elevation as contourf
                    ce = m.contourf(self.lon,self.lat,elev,clevs_elev,latlon=True,cmap="RdBu_r")

                    if i == 0:
                        divider = make_axes_locatable(ax)
                        cax = divider.append_axes("right", size="5%",pad=0.05)
                        cax.set_ylim([minValue,maxValue])

                        cbar = plt.colorbar(ce,cax=cax,orientation='vertical')
                        cbar.set_label('Elevation [m]')

                    # plot wind data as quiver
                    qw = m.quiver(self.lon[::xStep,::yStep],self.lat[::xStep,::yStep],wu,wv,latlon=True,alpha=.3,scale=150,width=0.005,pivot='middle')
                    ax.set_title(str(self.time[i]),fontsize=24)
                    plt.pause(0.1)


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
# Type the name of the gcmplt you want to analyze, without the ext (e.g., exp06)
exp = 'control_2010'
fname = '/media/danilo/Danilo/mestrado/ventopcse/output/%s.cdf'%(exp)

# instanciate an object, passing as argument, the full name:,timeStart=112,timeEnd=352
exp = Experiment(fname)
# find indexes for start and final date to extract data
timeStart = locate_closest_date(exp.ncin.time.values,year=2010,month=01,day=14)[0][0]
timeEnd   = locate_closest_date(exp.ncin.time.values,year=2010,month=02,day=13)[0][0]
exp.set_indexes4time(timeStart=timeStart,timeEnd=timeEnd)

# import basic informations
exp.importBasicInformations()

# plotting multiple elevations locations
iss = [50,7,22]
jss = [29,55,94]
exp.elevationPlot(iss,jss,points=['Laje de Santos','CSB','Ubatuba'])
