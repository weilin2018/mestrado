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
def crossSection(fname,DATA_DIR,limits,region='pcse',savefig=None):
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

    # refined variables for CSS
    depth = depth[40:90,:40]
    lon   = lon[40:90,:40]
    lat   = lat[40:90,:40]

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

    # cutting coordinates only for CSS
    # lon = lon[40:90,:40]
    # lat = lat[40:90,:40]
    # selecting latitudes for each section plot, based on Birocchi (2018)
    isul = 8  #46
    icen = 13 #60
    inor = 32 #80
    llat =-24
    ulat =-23.6
    llon =-45.62
    ulon =-45.18

    liminf_sul = 8  #limits[0][0]
    limsup_sul = 22 #limits[0][1]
    liminf_cen = 4  #limits[1][0]
    limsup_cen = 22 #limits[1][1]
    liminf_nor = 5  #limits[2][0]
    limsup_nor = 18 #limits[2][1]

    # for t in range(len(time)):
    for t in range(5):
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
        msouth = oceano.make_map(southMap_axis,resolution='i',ulat=ulat,ulon=ulon,llon=llon,llat=llat)
        mcentral = oceano.make_map(centralMap_axis,resolution='i',ulat=ulat,ulon=ulon,llon=llon,llat=llat)
        mnorth = oceano.make_map(northMap_axis,resolution='i',ulat=ulat,ulon=ulon,llon=llon,llat=llat)

        southsec_axis = plotCrossSection(southsec_axis,lon,lat,depth,sigma,isul,temp[:,:,:],limits=[liminf_sul,limsup_sul],depthReference=-50)
        southsec_axis.text(-47.4,-100,u'Cananéia',horizontalalignment='center')
        msouth.plot(lon[isul,liminf_sul:limsup_sul],lat[isul,liminf_sul:limsup_sul],'r',latlon=True)

        ################### santos
        centralsec_axis = plotCrossSection(centralsec_axis,lon,lat,depth,sigma,icen,temp[:,:,:],limits=[liminf_cen,limsup_cen],depthReference=-35)
        centralsec_axis.text(-46.3,-100,u'Santos',horizontalalignment='center')
        mcentral.plot(lon[icen,liminf_cen:limsup_cen],lat[icen,liminf_cen:limsup_cen],'r',latlon=True)

        ################### ubatuba
        northsec_axis = plotCrossSection(northsec_axis,lon,lat,depth,sigma,inor,temp[:,:,:],limits=[liminf_nor,limsup_nor],depthReference=-20)
        northsec_axis.text(-44.8,-100,u'Ubatuba',horizontalalignment='center')
        mnorth.plot(lon[inor,liminf_nor:limsup_nor],lat[inor,liminf_nor:limsup_nor],'r',latlon=True)

        # control time to the next plot
        if savefig:
            outname = savefig+str(t).zfill(4)+'.png'
            plt.savefig(outname)
            os.system('convert -trim %s %s'%(outname,outname))

        plt.pause(0.1)

def plotCrossSection(ax,lon,lat,depth,sigma,ind,temp,limits,depthReference):
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
    cs  = ax.contour(x[:,liminf:limsup],sig[:,liminf:limsup],conc[:,liminf:limsup],levels=np.arange(18.,28.,2.),colors=('k'),linestyles=('--'))
    ax.clabel(cs,fmt='%i',inline=True)
    ax.plot(lon[ind,liminf:limsup],-depth[ind,liminf:limsup],'k')
    ax.fill_between(lon[ind,liminf:limsup], -50, -depth[ind,liminf:limsup],color='#c0c0c0')
    ax.set_ylim([-32,1])
    # ax.margins(0)
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
SAVE_FIG = BASE_DIR + 'masterThesis_analysis/figures/experiments_outputs/temperature/crossSection_%s/css/'%(exp)

for f in fname:
    if exp in f:
        experiment = f

# to plot cross sections in Cananeia, Santos and Ubatuba:
# crossSection(experiment,DATA_DIR,region='pcse',limits=[[5,83],[5,82],[5,83]],savefig=SAVE_FIG)

# only for Sao Sebastiao Channel, we have to change all limits, so:
crossSection(experiment,DATA_DIR,region='css',limits=[[2,32],[5,23],[3,28]])






# - --
depth = ncdata['depth'][40:90,0:40].data.copy()
sigma = ncdata['sigma'].data.copy()
ss =np.tile(sigma,(37,1))
ss=np.transpose(ss)
prof =depth[latt,:]* ss
ll = np.tile(lon[latt,:],(37,1))
