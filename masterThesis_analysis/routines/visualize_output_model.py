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

def load_data(fname,vars='elev',startTime=112,endTime=352,sigma=None):
    """Function to load variables in vars from a netCDF file, given by fname.

    Parameters
    ----------
    fname : string
        String with full path to netCDF file to be read.
    vars : string
        Which variable to extract
    startTime : integer
        index of the beginning of the cut in time.
    endTime : integer
        index of the end of the cut in time.
    sigma : integer or None
        sigma level, from 0 to 36.

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

    if sigma==None:
        data = ncdata[vars].values[startTime:endTime,:,:]
    else:
        data = ncdata[vars].values[startTime:endTime,sigma,:,:]

    return lon,lat,time,data

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
        lon,lat,time,temp = load_data(fname,vars='temp',sigma=-1) # bottom
        lon,lat,time,salt = load_data(fname,vars='salt',sigma=0)  # surface

    tmean = temp.mean(axis=0)
    smean = salt.mean(axis=0)

    fig = plt.figure(figsize=(20,15))
    ax = fig.add_subplot(111)

    m = oceano.make_map(ax, resolution='i')
    # csf = m.contourf(lon,lat,tmean,latlon=True,cmap=cmo.cm.thermal)
    cst  = m.contour(lon,lat,tmean,latlon=True,levels=[18.],colors=('k'),linestyles=('--'))
    css  = m.contour(lon,lat,smean,latlon=True,levels=[36.5],colors=('r'),linestyles=('--'))

    plt.title('Mean Position for Isotherm of 18 (black) and Isohaline \n of 36.5 (red), During 15/Jan to 13/Feb of 2014',fontsize=24)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="20%",pad=0.05)
    cax.plot(temp[:,55,7])
    cax.set_xlim([0,240])
    cax.set_ylim([22,25])

    plt.clabel(cst,cst.levels,inline=True,fontsize=10,fmt='%2.1f')
    plt.clabel(css,css.levels,inline=True,fontsize=10,fmt='%2.1f')



    # plt.colorbar(csf,cax=cax)

    # minimap
    axins = zoomed_inset_axes(ax, 1.5, loc=4)
    # axins.plot(np.arange(0,len(time)), temp[:,55,7])

    axins.set_xlim(-10,0)
    axins.set_ylim(3,10)

    axins.xaxis.set_visible('False')
    axins.yaxis.set_visible('False')
    m2 = oceano.make_map(axins,llat=-24.23,ulat=-22.62,llon=-46.75,ulon=-43.29,resolution='i')
    m2.contourf(lon,lat,tmean,latlon=True,cmap=cmo.cm.thermal)

    # m2 = oceano.make_minimap()
    # pc = m2.contourf(lon,lat,tmean,latlon=True,cmap=cmo.cm.thermal)

    if savefig:
        outname = str(i).zfill(4)+'.png'
        plt.savefig(savefig+outname)
    else:
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

    fig = plt.figure(figsize=(20,15))
    ax = fig.add_subplot(111)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="20%",pad=0.05)
    cax.set_xlim([0,240])
    cax.set_ylim([22,25])

    for i in np.arange(0,temp.shape[0]-230):
        ax.clear()
        cax.clear()
        cax.set_xlim([0,240])
        cax.set_ylim([22,25])
        m = oceano.make_map(ax, resolution='i')
        csf = m.contourf(lon,lat,temp[i,:,:],latlon=True,cmap=cmo.cm.thermal)
        cs  = m.contour(lon,lat,temp[i,:,:],latlon=True,levels=[18.],colors=('k'),linestyles=('--'))
        # plt.clabel(cs, fmt='%2.1f',colors='k',fontsize=14)
        ts = pd.to_datetime(str(time[i]))
        plt.title(ts.strftime('%Y.%m.%d %H:%M'))

        # plt.gca().clear()

        tmp = temp[:i,55,7]
        # cax.title('Temporal evolution of bottom temperature in SBC')
        cax.plot(tmp,'k')
        cax.fill_between(np.arange(0,len(tmp)), 0, tmp,color='k',alpha=0.4)

        plt.pause(0.5)



        # plt.clabel(cst,cst.levels,inline=True,fontsize=10,fmt='%2.1f')
        # plt.clabel(css,css.levels,inline=True,fontsize=10,fmt='%2.1f')



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
# lon,lat,time,temp = tempMeanField(fname)


fig, ax = plt.subplots()
ax.set_xlim([0,240])
ax.set_ylim([22,25])

for i in np.arange(0,240):
    tmp = temp[:i,55,7]

    ax.plot(tmp,'k')
    ax.fill_between(np.arange(0,len(tmp)),0,tmp,color='k')

    plt.pause(1)
