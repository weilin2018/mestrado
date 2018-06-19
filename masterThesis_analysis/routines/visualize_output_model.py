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

def load_data(fname,vars='elev',startTime=None,endTime=None,sigma=None):
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
    savefig : string
        Directory to save figures. If is None, then an animation will be displayed.
    """

    # check if data is already loaded
    if not 'elev' in locals():
        lon,lat,time,elev = load_data(fname[-1])

    contour_levels = np.arange(-0.3,0.3,0.3/500)

    if not savefig:
        plt.ion()

        fig = plt.figure(figsize=(20,15))
        ax = fig.add_subplot(111)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("bottom", size="5%",pad=0.05)

        for i in np.arange(0,elev.shape[0]):
            ax.clear()
            cax.clear()
            m = oceano.make_map(ax, resolution='i')
            m.drawstates()
            cs = m.contourf(lon,lat,elev[i,:,:],contour_levels,latlon=True,cmap="RdBu_r")
            cbar = plt.colorbar(cs,orientation='horizontal',cax=cax)
            cbar.set_label('Elevation [m]')
            plt.title(str(time[i]),fontsize=24)

            plt.pause(.5)
    else:
        for i in np.arange(0,elev.shape[0]):
            fig, ax = plt.subplots()
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)

            m = oceano.make_map(ax,resolution='i')
            m.drawstates()
            cs = m.contourf(lon,lat,elev[i,:,:],contour_levels,latlon=True,cmap="RdBu_r")
            cbar = plt.colorbar(cs,orientation='horizontal',cax=cax)
            cbar.set_label('Elevation [m]')

            plt.colorbar(cs,cax=cax)
            plt.suptitle(str(time[i]),fontsize=24)

            outname = str(i).zfill(4)+'.png'
            plt.savefig(savefig+outname)
            plt.close("all")

def elevationPlot(fname,ilat=55,jlon=7):

    # check if data is already loaded
    if not 'elev' in locals():
        lon,lat,time,elev = load_data(fname)

    # extract from a specific point
    eta = elev[:,ilat,jlon]

    # create a dataframe, because it's beautiful
    df = pd.DataFrame({'eta':eta},index=pd.DatetimeIndex(time))

    # plot and show
    df.plot(figsize=(15,20),title=u'Série Temporal de Elevação do Nível do Mar para o período de 15/01 à 13/02 de 2014')
    plt.margins(0)
    plt.show()

def tempMeanField(fname,savefig=None):

    fig = plt.figure(figsize=(20,15))
    ax = fig.add_subplot(111)

    m = oceano.make_map(ax, resolution='i')

    # plot isotherm
    lon,lat,time,temp = load_data(fname[0],vars='temp',sigma=-1) # bottom
    tmean = temp.mean(axis=0)
    # extract temporal serie from Sao Sebastial Channel
    temp = temp[:,55,7]
    csft = m.contourf(lon,lat,tmean,latlon=True,cmap=cmo.cm.thermal)
    cst  = m.contour(lon,lat,tmean,latlon=True,levels=[18.],colors=('k'),linestyles=('--'))
    

    # plot isohaline
    lon,lat,time,salt = load_data(fname[2],vars='salt',sigma=0)  # surface
    smean = salt.mean(axis=0)
    del salt
    csfs = m.contourf(lon,lat,smean,latlon=True,cmap=cmo.cm.haline)
    css  = m.contour(lon,lat,smean,latlon=True,levels=[36.5],colors=('r'),linestyles=('--'))

    plt.title('Mean Position for Isotherm of 18 (black) and Isohaline \n of 36.5 (red), Between 15/Jan to 13/Feb of 2014',fontsize=24)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="20%",pad=0.05)
    cax.plot(temp)
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
    cs = m2.contourf(lon,lat,tmean,latlon=True,cmap=cmo.cm.thermal)

    #plt.colorbar(cs,location='bottom')

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

    # check if data is already loaded
    if not 'temp' in locals():
        lon,lat,time,temp = load_data(fname[0],vars='temp',sigma=-1)

    # testing for shape
    if len(temp.shape) == 4:
        temp = temp[:,-1,:,:] # select last sigma level

    # creating a mask for values greater than 19. and lower than 18.
    temp_masked = np.ma.masked_less(temp, 18.)
    temp_masked = np.ma.masked_greater(temp_masked, 19.)

    if not savefig:
        plt.ion()
        # animation
        fig = plt.figure(figsize=(20,15))
        ax = fig.add_subplot(111)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("bottom", size="20%",pad=0.05)
        cax.set_xlim([0,240])
        cax.set_ylim([22,25])

        for i in np.arange(0,temp.shape[0]-235):
            ax.clear()
            cax.clear()
            cax.set_xlim([0,240])
            cax.set_ylim([22,25])
            m = oceano.make_map(ax, resolution='i')
            csf = m.contourf(lon,lat,temp_masked[i,:,:],latlon=True,cmap=cmo.cm.thermal)
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
    else:
        # save figures

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

            outname = str(i).zfill(4)+'.png'
            plt.savefig(savefig+outname)



        # plt.clabel(cst,cst.levels,inline=True,fontsize=10,fmt='%2.1f')
        # plt.clabel(css,css.levels,inline=True,fontsize=10,fmt='%2.1f')

def meanSubplots(fname,savefig=None):

    # importing temperature data and plotting
    lon,lat,time,temp = load_data(fname[0],vars='temp',sigma=-1) # bottom sigma level
    tmean = temp.mean(axis=0)
    # extract temporal serie from Sao Sebastial Channel
    temp = temp[:,55,7]

    # importing salinity data
    lon,lat,time,salt = load_data(fname[2],vars='salt',sigma=0) # surface sigma level
    smean = salt.mean(axis=0)
    # extract temporal serie
    salt = salt[:,55,7]

    # creating structure and plots configuration
    fig,axes = plt.subplots(ncols=2,nrows=1,figsize=(20,15))
    # axes for colorbar
    divider = make_axes_locatable(axes[0])
    ctemp   = divider.append_axes("right", size="5%", pad=0.05)
    ctemp.set_ylim([22,25])
    # axes for plot
    divider = make_axes_locatable(axes[0])
    cax1 = divider.append_axes("bottom", size="20%",pad=0.6)
    cax1.set_xlim([0,240])
    cax1.set_ylim([22,25])
    #cax1.clabel(r'Temperature [$^oC$]')

    # axes for colorbar (salinity)
    divider = make_axes_locatable(axes[1])
    csalt   = divider.append_axes("right", size="5%", pad=0.05)
    # axes for plot
    divider = make_axes_locatable(axes[1])
    cax2 = divider.append_axes("bottom", size="20%",pad=0.6)
    cax2.set_xlim([0,240])
    cax2.set_ylim([34.7,35.2])
    #cax2.clabel(r'Salinity')

    # temperature field
    m = oceano.make_map(axes[0], resolution='i')
    # plot isotherm
    csft = m.contourf(lon,lat,tmean,latlon=True,cmap=cmo.cm.thermal)
    cst  = m.contour(lon,lat,tmean,latlon=True,levels=[18.],colors=('k'),linestyles=('--'))
    plt.clabel(cst,cst.levels,inline=True,fontsize=10,fmt='%2.1f')
    cbar = plt.colorbar(csft,cax=ctemp)
    cbar.set_label(r'Temperature [$^oC$]')
    # plot location
    m.scatter(lon[55,7],lat[55,7],s=70,marker='o',color='black',alpha=0.5,latlon=True,zorder=10,linewidth=2.)
    # temperature plot
    cax1.plot(temp,'k')
    cax1.set_ylabel(r'Temperature [$^oC$]')
    cax1.set_xlabel(r'Time [in timesteps]')

    # salinity field
    m = oceano.make_map(axes[1], resolution='i')
    # plot isotherm
    csfs = m.contourf(lon,lat,smean,latlon=True,cmap=cmo.cm.haline)
    css  = m.contour(lon,lat,smean,latlon=True,levels=[36.5],colors=('r'),linestyles=('--'))
    plt.clabel(css,css.levels,inline=True,fontsize=10,fmt='%2.1f')
    cbar = plt.colorbar(csfs,cax=csalt)
    cbar.set_label(r'Salinity')
    # plot location
    m.scatter(lon[55,7],lat[55,7],s=70,marker='o',color='black',alpha=0.5,latlon=True,zorder=10,linewidth=2.)
    # salinity plot
    cax2.plot(salt,'r')
    cax2.set_ylabel(r'Salinity')
    cax2.set_xlabel(r'Time [in timesteps]')

    # setting title
    title=u"Mean Field of Bottom Temperature (left) and Surface Salinity (right) for the\n Period Between 15.01 to 13.04, 2014. Below the \n Temporal Evolution of Each Variable is Shown in the Same Period."
    plt.suptitle(title,fontsize=22)

    if savefig:
        outname = 'meanField_temp_and_salt.png'
        plt.savefig(savefig+outname)
    else:
        plt.show()

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

BASE_DIR = oceano.make_dir()
if BASE_DIR.split("/")[2] == 'tparente':
    DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output_modelo/exp03_variables/')
else:
    DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')

FIGU_DIR = BASE_DIR + 'masterThesis_analysis/figures/experiments_outputs/elevation/'

#OUT_FILE = DATA_DIR+INP_FILE.replace('cdf','pickle')

fname = glob.glob(DATA_DIR+"*.cdf")
#
# elevationPlot(fname)
elevationField(fname,savefig=FIGU_DIR)
# lon,lat,time,temp = tempMeanField(fname)
