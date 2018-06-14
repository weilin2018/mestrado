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

# pacotes para minimap
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################

def load_data(fname,vars='elev',startTime=112,endTime=352):
    """Function to load variables in vars from a netCDF file, given by fname.

    Parameters
    ----------
    fname : string
        String with full path to netCDF file to be read.
    vars : string
        Which variable to extract

    Returns
    -------
    lon,lat,variable : xarray.Dataset
    """
    # load file
    ncdata = xr.open_dataset(fname)

	# extract variables
    lon,lat = ncdata['lon'].values, ncdata['lat'].values
    lon[lon == 0.] = np.nan
    lat[lat == 0.] = np.nan
    time = ncdata['time'].values[startTime:endTime]

    data = ncdata[vars].values[startTime:endTime,:,:]

    return lon,lat,time,data

def plotData(m,lon,lat,data,contour_levels):

    cs = m.contourf(lon,lat,data,contour_levels,latlon=True)

def locatePCIandPCMlimits(depth):
    """Function to locate all indexes related to PCI and PCM limits, based on
    local depth.

    Parameters
    ----------
    depth : numpy.ndarray
        Variable with gridded depth.

    Returns
    -------
    ilons :
        Indexes from longitudes
    ilats :
        Indexes from latitudes
    """

    # find indexes for depth < 200

    return False

def elevationField(fname,savefig=None):
    """Plotting elevation field from experiment given by fname.

    Parameters
    ----------
    fname : string
        Full path to netCDF file

    """

    # check if data is already loaded
    if not 'elev' in locals():
        lon,lat,time,elev = load_data(fname)

    contour_levels = np.arange(-0.3,0.3,0.3/500)

    if debug:
        fig, ax = plt.subplots()

        for i in np.arange(0,len(time)):
            plt.gca().clear()
            m = oceano.make_map(ax, resolution='i')
            m.drawstates()
            cs = m.contourf(lon,lat,elev[i,:,:],latlon=True)

            plt.title(str(time[i]),fontsize=24)

            plt.pause(.5)
    else:
        for i in np.arange(0,len(time)):
            fig, ax = plt.subplots()

            m = oceano.make_map(ax,resolution='i')
            m.drawstates()
            cs = m.contourf(lon,lat,elev[i,:,:],contour_levels,latlon=True)

            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)

            plt.colorbar(cs,cax=cax)
            plt.title(str(time[i]),fontsize=24)

            if savefig:
                outname = str(i).zfill(4)+'.png'
                plt.savefig(savefig+outname)
            else:
                plt.show()

def elevationPlot(fname,ilat=55,jlon=7):

    # check if data is already loaded
    if not 'elev' in locals():
        lon,lat,time,elev = load_data(fname)

    eta = elev[:,ilat,jlon]

    df = pd.DataFrame({'eta':eta},index=pd.DatetimeIndex(time))

    df.plot(figsize=(15,20),title=u'Série Temporal de Elevação do Nível do Mar para o período de 15/01 à 13/02 de 2014')
    plt.margins(0)
    plt.show()

def tempMeanField(fname,savefig=None):

    # check if data is already loaded
    if not 'temp' in locals():
        lon,lat,time,temp = load_data(fname,vars='temp')

    # testing for shape
    if len(temp.shape) == 4:
        temp = temp[:,-1,:,:] # select last sigma level

    tmean = temp.mean(axis=0)

    fig,ax = plt.subplots()
    m = oceano.make_map(ax, resolution='i')
    csf = m.contourf(lon,lat,tmean,latlon=True,cmap=cmo.cm.thermal)
    cs  = m.contour(lon,lat,tmean,latlon=True,levels=[18.],colors=('k'),linestyles=('--'))

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cfs,cax=cax)

    plt.title('Mean Position for Isotherm of 18, \n During 15/Jan to 13/Feb of 2014',fontsize=24)

    # minimap
    axins = zoomed_inset_axes(m.ax,5,loc=4)
    axins.set_xlim(-10,0)
    axins.set_ylim(3,10)



    # m2 = oceano.make_minimap()
    pc = m2.contourf(lon,lat,tmean,latlon=True,cmap=cmo.cm.thermal)

    plt.show()



def temperatureField(fname,isotherm=18.,savefig=None):
    """Function to plot ONLY the location of isotherm line of 18ºC, associated with
    ACAS.

    Parameters
    ----------
    fname : string
        Full path of netCDF file.
    isotherm : float
        Which isotherm to plot.
    savefig : string
        Full path and name to save figure.
    """

    plt.ion()

    # check if data is already loaded
    if not 'temp' in locals():
        lon,lat,time,temp = load_data(fname,vars='temp')

    # testing for shape
    if len(temp.shape) == 4:
        temp = temp[:,-1,:,:] # select last sigma level

    # creating a mask for values greater than 19. and lower than 18.
    temp_masked = np.ma.masked_less(temp, 18.)
    temp_masked = np.ma.masked_greater(temp_masked, 19.)

    fig, ax = plt.subplots()

    for i in np.arange(0,temp.shape[0]):
        m = oceano.make_map(ax, resolution='i')
        csf = m.contourf(lon,lat,temp[i,:,:],latlon=True,cmap=cmo.cm.thermal)
        cs  = m.contour(lon,lat,temp[i,:,:],latlon=True,levels=[18.],colors=('k'),linestyles=('--'))
        # plt.clabel(cs, fmt='%2.1f',colors='k',fontsize=14)
        ts = pd.to_datetime(str(time[i]))
        plt.title(ts.strftime('%Y.%m.%d %H:%M'))
        plt.pause(1)
        plt.gca().clear()


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
FIGU_DIR = BASE_DIR + 'masterThesis_analysis/figures/experiments_outputs/elevation/'
#OUT_FILE = DATA_DIR+INP_FILE.replace('cdf','pickle')

fname = glob.glob(DATA_DIR+"*.cdf")[-1]
#
# elevationPlot(fname)
# elevationField(fname,savefig=FIGU_DIR)
