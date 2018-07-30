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
def crossSection(fname,DATA_DIR,savefig=None):
    """Main function to plot cross section of temperature.

    Parameters
    ----------
    DATA_DIR : string
        Full path to the local with netCDF's files (output model).
    savefig : string
        Full path to the directory.
    """
    # clear screen
    os.system('clear')

    ncdata = xr.open_dataset(fname)

    startTime=112
    endTime=352

    # extract variables
    lon,lat = ncdata['lon'].values, ncdata['lat'].values
    lon[lon == 0.] = np.nan
    lat[lat == 0.] = np.nan
    time = ncdata['time'].values[startTime:endTime]
    depth = ncdata['depth'].values
    sigma = ncdata['sigma'].values

    # settings configurations for the figure
    plt.figure(figsize=[16/2.54,19/2.54])

    grid = plt.GridSpec(3,3,wspace=0.5,hspace=0.3)

    northsec_axis = plt.subplot(grid[0,:2])
    northMap_axis = plt.subplot(grid[0,2])
    mnorth = oceano.make_map(northMap_axis,resolution='i')

    centralsec_axis = plt.subplot(grid[1,:2])
    centralMap_axis = plt.subplot(grid[1,2])
    mcentral = oceano.make_map(centralMap_axis,resolution='i')

    southsec_axis = plt.subplot(grid[2,:2])
    southMap_axis = plt.subplot(grid[2,2])
    msouth = oceano.make_map(southMap_axis,resolution='i')

    # select latitude index for cross section
    isul = 19
    icen = 28
    inor = 99

    for t in range(len(time)):
        print('Timestep: %i'%(t))
        temp = ncdata.temp[t]           # import data
        # plot data

        # cleaning axis
        southsec_axis.clear()
        centralsec_axis.clear()
        northsec_axis.clear()
        southMap_axis.clear()
        centralMap_axis.clear()
        northMap_axis.clear()
        # updating time
        plt.suptitle(time[t],fontsize=24)

        # plotting data for each location
        msouth = oceano.make_map(southMap_axis,resolution='i')
        mcentral = oceano.make_map(centralMap_axis,resolution='i')
        mnorth = oceano.make_map(northMap_axis,resolution='i')

        southsec_axis = plotCrossSection(southsec_axis,lon,lat,depth,sigma,isul,temp[:,:,:],limits=[5,83])
        southsec_axis.text(-47.4,-100,u'Cananéia',horizontalalignment='center')
        msouth.plot(lon[isul,5:83],lat[isul,5:83],'r',latlon=True)

        ################### santos
        centralsec_axis = plotCrossSection(centralsec_axis,lon,lat,depth,sigma,icen,temp[:,:,:],limits=[5,82])
        centralsec_axis.text(-46.3,-100,u'Santos',horizontalalignment='center')
        mcentral.plot(lon[icen,5:82],lat[icen,5:82],'r',latlon=True)

        ################### ubatuba
        northsec_axis = plotCrossSection(northsec_axis,lon,lat,depth,sigma,inor,temp[:,:,:],limits=[5,83])
        northsec_axis.text(-44.8,-100,u'Ubatuba',horizontalalignment='center')
        mnorth.plot(lon[inor,5:83],lat[inor,5:83],'r',latlon=True)

        # control time to the next plot
        if savefig:
            plt.savefig(savefig+str(t).zfill(4)+'.png')

        plt.pause(0.1)

def plotCrossSection(ax,lon,lat,depth,sigma,ind,temp,limits):
    """Apenas plota os dados em uma secao vertical, fazendo a conversão de nível
    sigma para profundidade.

    Parameters
    ----------
    ax : matplotlib.axes._subplots.AxesSubplot
        Axis to plot data.
    lon : numpy.ndarray
        Longitude vector 2D.
    depth : v
        Depth field (bathymetry).
    sigma : numpy.ndarray
        Array with sigma levels.
    ind : integer
        Index of each latitude is to be plotted.
    temp : numpy.ndarray
        Array with temperature data.
    limits : list
        Limits of non-nan values

    Returns
    -------
    ax : matplotlib.axes._subplots.AxesSubplot
    """

    # create contour levels
    contour_levels = np.arange(10,27,20/200.)
    ################### cananeia
    x,prof,sig = oceano.create_newDepth(lon,depth,sigma,ind)      # create new depth
    conc = temp[:,ind,:]            # extract cross section data
    liminf,limsup = limits[0],limits[1]               # limits with non-nan values

    cfs = ax.contourf(x[:,liminf:limsup],sig[:,liminf:limsup],conc[:,liminf:limsup],contour_levels,cmap=cmo.cm.thermal)
    cs  = ax.contour(x[:,liminf:limsup],sig[:,liminf:limsup],conc[:,liminf:limsup],levels=[18.],colors=('k'),linestyles=('--'))
    ax.plot(lon[ind,liminf:limsup],-depth[ind,liminf:limsup],'k')
    ax.fill_between(lon[ind,liminf:limsup], -200, -depth[ind,liminf:limsup],color='#c0c0c0')
    ax.set_ylim([-200,1])
    ax.margins(0)
    ax.set_ylabel(u'Depth [m]',fontsize=18)

    return ax

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

BASE_DIR = oceano.make_dir()
# if BASE_DIR.split("/")[2] == 'tparente':
#     DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output_modelo/exp03_variables/')
#     fname = glob.glob(DATA_DIR+"*.cdf")
# else:
#     DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
#     fname = glob.glob(DATA_DIR+"*.cdf")

DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
fname = glob.glob(DATA_DIR+"*.cdf")


FIGU_DIR = BASE_DIR + 'masterThesis_analysis/figures/experiments_outputs/elevation/'


# select which experiment you want to plot:
exp = 'exp06'
SAVE_FIG = BASE_DIR + 'masterThesis_analysis/figures/experiments_outputs/temperature/crossSection_%s/'%(exp)

for f in fname:
    if exp in f:
        experiment = f

crossSection(experiment,DATA_DIR,savefig=SAVE_FIG)
