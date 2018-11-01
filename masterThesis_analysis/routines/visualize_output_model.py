#-*-coding;utf-8-*-
%reset -f
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
        lon,lat,time,elev = load_data(fname)

    contour_levels = np.arange(-0.3,0.3,0.3/500)

    if not savefig:
        plt.ion()
        # creating structure and plots configuration
        fig,ax = plt.subplots(ncols=1,nrows=1,figsize=(20,15))

        # axes for colorbar
        divider = make_axes_locatable(ax)
        # celev   = divider.append_axes("right", size="5%", pad=0.05)
        #celev.set_ylim([22,25])
        # axes for plot
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("bottom", size="20%",pad=0.6)
        cax.set_xlim([0,240])
        cax.set_ylim([-0.4,0.4])

        for i in np.arange(0,time.shape[0]):
            ax.clear()
            cax.clear()
            plt.title(str(time[i]),fontsize=24)

            m = oceano.make_map(ax, resolution='i')
            m.drawstates()
            cs = m.contourf(lon,lat,elev[i,:,:],contour_levels,latlon=True,cmap="RdBu_r")
            # cbar = plt.colorbar(cs,orientation='horizontal',cax=cax)
            # cbar.set_label('Elevation [m]')


            tmp = elev[:i,55,7]
            # cax.title('Temporal evolution of bottom temperature in SBC')
            cax.set_xlim([0,240])
            cax.set_ylim([-0.3,0.3])
            cax.plot(tmp,'k')
            cax.fill_between(np.arange(0,len(tmp)), -1., tmp,color='k',alpha=0.4)


            plt.pause(0.5)
    else:
        for i in np.arange(0,elev.shape[0]):
            os.system('clear')
            print("Plotting timestep: %i" %(i))
            # creating structure and plots configuration
            fig,ax = plt.subplots(ncols=1,nrows=1,figsize=(20,15))

            # axes for colorbar
            divider = make_axes_locatable(ax)
            celev   = divider.append_axes("right", size="5%", pad=0.05)
            celev.set_ylim([-0.3,0.3])
            # axes for plot
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("bottom", size="20%",pad=0.6)
            cax.set_xlim([0,240])
            cax.set_ylim([-0.3,0.3])
            plt.suptitle(str(time[i]),fontsize=24)

            m = oceano.make_map(ax, resolution='i')
            m.drawstates()
            cs = m.contourf(lon,lat,elev[i,:,:],contour_levels,latlon=True,cmap="RdBu_r")
            cbar = plt.colorbar(cs,orientation='vertical',cax=celev)
            cbar.set_label('Elevation [m]')


            tmp = elev[:i,55,7]
            # cax.title('Temporal evolution of bottom temperature in SBC')
            cax.set_xlim([0,240])
            cax.set_ylim([-0.2,0.1])
            cax.plot(tmp,'k')
            cax.fill_between(np.arange(0,len(tmp)), -1., tmp,color='k',alpha=0.4)

            outname = str(i).zfill(4)+'.png'
            outname = savefig+outname
            plt.savefig(outname)
            os.system('convert -trim %s %s' % (outname,outname))

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

def temperatureField(fname,isotherm=18.,startTime=46,endTime=304,savefig=None):
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
        lon,lat,time,temp = load_data(fname,vars='temp',sigma=-1,startTime=startTime,endTime=endTime)

    # testing for shape
    if len(temp.shape) == 4:
        temp = temp[:,-1,:,:] # select last sigma level

    # creating a mask for values greater than 19. and lower than 18.
    # temp_masked = np.ma.masked_less(temp, 18.)
    # temp_masked = np.ma.masked_greater(temp_masked, 19.)
    temp_masked = temp.copy()

    if not savefig:
        plt.ion()
        # animation
        fig = plt.figure(figsize=(20,15))
        ax = fig.add_subplot(111)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("bottom", size="20%",pad=0.05)
        cax.set_xlim([0,240])
        cax.set_ylim([15,32])

        for i in np.arange(0,temp.shape[0]):
            ax.clear()
            cax.clear()
            cax.set_xlim([0,260])
            cax.set_ylim([15,32])
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
        os.system('clear')
        for i in np.arange(0,temp.shape[0]):
            print('plotting temp for timestep %s\n' % (i))
            fig,ax = plt.subplots(figsize=(20,15))

            divider = make_axes_locatable(ax)
            cax = divider.append_axes("bottom", size="20%",pad=0.05)
            cax.set_xlim([0,260])
            cax.set_ylim([15,32])

            m = oceano.make_map(ax, resolution='i')
            # csf = m.contourf(lon,lat,temp[i,:,:],latlon=True,cmap=cmo.cm.thermal)
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
            plt.close()

def salinityField(fname,isohaline=36.5,savefig=None):

    # check if data is already loaded
    if not 'salt' in locals():
        lon,lat,time,salt = load_data(fname,vars='salt',sigma=0,startTime=112,endTime=352)

    # testing for shape
    if len(salt.shape) == 4:
        salt = salt[:,0,:,:] # select first sigma level

    # creating a mask for values greater than 19. and lower than 18.
    salt_masked = np.ma.masked_less(salt, 18.)
    salt_masked = np.ma.masked_greater(salt_masked, 19.)

    contour_levels = np.arange(32,37,5/100.)

    if not savefig:
        plt.ion()
        # animation
        fig = plt.figure(figsize=(20,15))
        ax = fig.add_subplot(111)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("bottom", size="20%",pad=0.05)
        cax.set_xlim([0,240])
        # cax.set_ylim([22,25])

        for i in np.arange(0,salt.shape[0]):
            ax.clear()
            cax.clear()
            cax.set_xlim([0,240])
            # cax.set_ylim([22,25])
            m = oceano.make_map(ax, resolution='i')
            #csf = m.contourf(lon,lat,salt_masked[i,:,:],latlon=True,cmap=cmo.cm.haline)
            cs  = m.contour(lon,lat,salt[i,:,:],contour_levels,latlon=True,levels=[36.5],colors=('k'),linestyles=('--'))
            # plt.clabel(cs, fmt='%2.1f',colors='k',fontsize=14)
            ts = pd.to_datetime(str(time[i]))
            plt.title(ts.strftime('%Y.%m.%d %H:%M'))

            # plt.gca().clear()

            slt = salt[:i,55,7]
            # cax.title('Temporal evolution of bottom temperature in SBC')
            cax.plot(slt,'k')
            cax.fill_between(np.arange(0,len(slt)), 0, slt,color='k',alpha=0.4)

            plt.pause(0.5)
    else:
        # save figures
        os.system('clear')
        for i in np.arange(0,salt.shape[0]):
            print('plotting salt for timestep %s\n' % (i))
            fig,ax = plt.subplots(figsize=(20,15))
            divider = make_axes_locatable(ax)
            csalt   = divider.append_axes("right", size="5%", pad=0.05)
            csalt.set_ylim([33.0,37.0])

            m = oceano.make_map(ax, resolution='i')
            csf = m.contourf(lon,lat,salt[i,:,:],contour_levels,latlon=True,cmap=cmo.cm.haline)
            # cs  = m.contour(lon,lat,salt[i,:,:],latlon=True,levels=[36.5],colors=('k'),linestyles=('--'))
            # plt.clabel(cs, fmt='%2.1f',colors='k',fontsize=14)
            ts = pd.to_datetime(str(time[i]))
            plt.title(ts.strftime('%Y.%m.%d %H:%M'))

            cbar = plt.colorbar(csf,orientation='vertical',cax=csalt)
            cbar.set_label('Surface Salinity')

            outname = str(i).zfill(4)+'.png'
            plt.savefig(savefig+outname)

            os.system('convert -trim %s %s'%(savefig+outname,savefig+outname))

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

def velocityField(fname,savefig=None):

    # check if data is already loaded
    if not 'spd' in locals():
        lon,lat,time,u = load_data(fname,vars='u',startTime=112,endTime=352)
        lon,lat,time,u = load_data(fname,vars='u',startTime=112,endTime=352)

    contour_levels = np.arange(0.,5.,5./1000)

    if not savefig:
        plt.ion()
        # creating structure and plots configuration
        fig,ax = plt.subplots(ncols=1,nrows=1,figsize=(20,15))

        # axes for colorbar
        divider = make_axes_locatable(ax)
        # celev   = divider.append_axes("right", size="5%", pad=0.05)
        #celev.set_ylim([22,25])
        # axes for plot
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("bottom", size="20%",pad=0.6)
        cax.set_xlim([0,240])
        cax.set_ylim([-0.4,0.4])

        for i in np.arange(0,time.shape[0]):
            ax.clear()
            cax.clear()
            plt.title(str(time[i]),fontsize=24)

            m = oceano.make_map(ax, resolution='i')
            m.drawstates()
            cs = m.contourf(lon,lat,spd[i,:,:],contour_levels,latlon=True,cmap=cmo.cm.speed)
            cq = m.quiver(lon[::2,::2],lat[::2,::2],u[i,::2,::2],v[i,::2,::2],latlon=True)
            # cbar = plt.colorbar(cs,orientation='horizontal',cax=cax)
            # cbar.set_label('Elevation [m]')


            tmpu = u[:i,55,7]
            tmpv = v[:i,55,7]

            # cax.title('Temporal evolution of bottom temperature in SBC')
            cax.set_xlim([0,240])
            cax.set_ylim([-0.4,0.4])
            cax.plot(tmpu,'b')
            cax.plot(tmpv,'r')
            # cax.fill_between(np.arange(0,len(tmp)), -1., tmp,color='k',alpha=0.4)


            plt.pause(0.1)

def compareExperiments(fname1,fname2,var='temp',sigma=-1,savefig=None,colorbar_limits=[22,25]):
    """Routine to plot data field between two experiments.

    Parameters
    ----------
    fname1 : string
        Full path to experiment 1
    fname2 : string
        Full path to experiment 2
    var : string
        Which variable to plot. Default is Temperature associated to ..
    sigma : integer
        Which sigma level to plot. Default is ... -1, associated to the bottom.
    savefig : string
        FIGU_DIR with the full path for the directory to save figures. Default
        is None.
    colorbar_limits : vector
        Vector with inferior and superior limit for the colorbar.

    """

    # load exp03
    lon,lat,time1,temp1 = load_data(fname1,vars=var,sigma=sigma,startTime=112,endTime=352)
    # load exp04
    lon,lat,time2,temp2 = load_data(fname2,vars=var,sigma=sigma,startTime=32,endTime=272)

    if not savefig:
        plt.ion()
        # animation
        fig,ax = plt.subplots(figsize=(15,10), nrows=1, ncols=2)

        divider = make_axes_locatable(ax[0])
        cax1 = divider.append_axes("bottom", size="20%",pad=0.05)
        cax1.set_xlim([0,240])
        cax1.set_ylim(colorbar_limits)

        cbar1 = divider.append_axes("right", size="5%",pad=.5)
        cbar1.set_ylim(colorbar_limits)

        divider = make_axes_locatable(ax[1])
        cax2 = divider.append_axes("bottom", size="20%",pad=0.05)
        cax2.set_xlim([0,240])
        cax2.set_ylim(colorbar_limits)
        cbar2 = divider.append_axes("right", size="5%",pad=.5)
        cbar2.set_ylim(colorbar_limits)

        ax[0].set_title('Exp03')
        ax[1].set_title('Exp04')

        for i in np.arange(0,240-235):
            ax[0].clear()
            ax[1].clear()
            cax1.clear()
            cax1.set_xlim([0,240])
            cax1.set_ylim(colorbar_limits)
            cax2.clear()
            cax2.set_xlim([0,240])
            cax2.set_ylim(colorbar_limits)

            # plot first data
            m = oceano.make_map(ax[0], resolution='i')
            csf = m.contourf(lon,lat,temp1[i,:,:],latlon=True,cmap=cmo.cm.thermal)
            cs  = m.contour(lon,lat,temp1[i,:,:],latlon=True,levels=[18.],colors=('k'),linestyles=('--'))
            # plt.clabel(cs, fmt='%2.1f',colors='k',fontsize=14)
            plt.colorbar(csf,cax=cax1,orientation='horizontal')

            # plot second data
            m = oceano.make_map(ax[1], resolution='i')
            csf = m.contourf(lon,lat,temp2[i,:,:],latlon=True,cmap=cmo.cm.thermal)
            cs  = m.contour(lon,lat,temp2[i,:,:],latlon=True,levels=[18.],colors=('k'),linestyles=('--'))
            # plt.clabel(cs, fmt='%2.1f',colors='k',fontsize=14)
            plt.colorbar(csf,cax=cax2,orientation='horizontal')

            ts = pd.to_datetime(str(time2[i]))
            plt.suptitle(ts.strftime('%Y.%m.%d %H:%M'),fontsize=24)

            plt.pause(0.5)
        else:
            os.system("clear")
            for i in np.arange(0,240):
                print('Plotting timestep: %i' % (i))
                fig,ax = plt.subplots(figsize=(15,10), nrows=1, ncols=2)

                divider = make_axes_locatable(ax[0])
                cax1 = divider.append_axes("bottom", size="20%",pad=0.05)
                cax1.set_xlim([0,240])
                cax1.set_ylim(colorbar_limits)

                cbar1 = divider.append_axes("right", size="5%",pad=.5)
                cbar1.set_ylim(colorbar_limits)

                divider = make_axes_locatable(ax[1])
                cax2 = divider.append_axes("bottom", size="20%",pad=0.05)
                cax2.set_xlim([0,240])
                cax2.set_ylim(colorbar_limits)
                cbar2 = divider.append_axes("right", size="5%",pad=.5)
                cbar2.set_ylim(colorbar_limits)

                ax[0].set_title('Exp03')
                ax[1].set_title('Exp04')

                # plot first data
                m = oceano.make_map(ax[0], resolution='i')
                csf = m.contourf(lon,lat,temp1[i,:,:],latlon=True,cmap=cmo.cm.thermal)
                cs  = m.contour(lon,lat,temp1[i,:,:],latlon=True,levels=[18.],colors=('k'),linestyles=('--'))
                # plt.clabel(cs, fmt='%2.1f',colors='k',fontsize=14)
                plt.colorbar(csf,cax=cbar1,orientation='vertical')
                tmp1 = temp1[:i,55,7]
                cax1.plot(tmp1)
                cax1.fill_between(np.arange(0,len(tmp1)), 0, tmp1,color='k',alpha=0.4)

                # plot second data
                m = oceano.make_map(ax[1], resolution='i')
                csf = m.contourf(lon,lat,temp2[i,:,:],latlon=True,cmap=cmo.cm.thermal)
                cs  = m.contour(lon,lat,temp2[i,:,:],latlon=True,levels=[18.],colors=('k'),linestyles=('--'))
                # plt.clabel(cs, fmt='%2.1f',colors='k',fontsize=14)
                plt.colorbar(csf,cax=cbar2,orientation='vertical')
                tmp2 = temp2[:i,55,7]
                cax2.plot(tmp2)
                cax2.fill_between(np.arange(0,len(tmp2)), 0, tmp2,color='k',alpha=0.4)

                ts = pd.to_datetime(str(time2[i]))
                plt.suptitle(ts.strftime('%Y.%m.%d %H:%M'),fontsize=24)

                outname = str(i).zfill(4)+'.png'
                plt.savefig(savefig+outname)
                plt.close("all")
                os.system('convert -trim %s %s'%(savefig+outname,savefig+outname))

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

BASE_DIR = oceano.make_dir()
if BASE_DIR.split("/")[2] == 'tparente':
    DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output_modelo/exp03_variables/')
    fname = glob.glob(DATA_DIR+"*.cdf")
else:
    DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
    fname = glob.glob(DATA_DIR+"*.cdf")

FIGU_DIR = BASE_DIR + 'masterThesis_analysis/figures/experiments_outputs/elevation/'


# select which experiment you want to plot:
exp = 'EA1'
SAVE_FIG = BASE_DIR + 'masterThesis_analysis/figures/experiments_outputs/temperature/crossSection_%s/'%(exp)

for f in fname:
    if exp in f:
        experiment = f


# temperatureField(experiment,savefig=FIGU_DIR.replace('elevation', 'temperature/isotherm18_control2010'))

# - ----------------------------------------------

os.system('clear')
var = input("type which variable you want to plot: 1 - elevation, 2 - isotherm, 3 - Salinity Field and 0 - to exit: ")

sav = input('You want to [0] visualize or [1] save figures? ')

if var == 1:
    if sav == 0:
        elevationField(experiment)
    else:
        elevationField(experiment,savefig=FIGU_DIR)
elif var == 2:
    if sav == 0:
        temperatureField(experiment,startTime=46,endTime=304)
    else:
        temperatureField(experiment,savefig=FIGU_DIR.replace('elevation', 'temperature/isotherm18_EA1'))
elif var == 3:
    if sav == 0:
        salinityField(experiment)
    else:
        salinityField(experiment,savefig=FIGU_DIR.replace('elevation', 'salinity'))
else:
    exit




#######################
