# using my package for model's visualization

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
from scipy import interpolate

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

from modelVisualization.interface import Experiment

def interp_ghrsst_to_ecom(x,y,data,xi,yi,cut=True,method='linear'):
    """
        x,y = grade original
        data= dado a ser interpolado
        xi,yi = grade para interpolar
    """

    # selecting some values inside a box to speed up interpolation
    if cut:
        ind = (x > -48.8) & (y < -21.663) & (y > -29.669)
        x = x[ind]
        y = y[ind]
        data = data[ind]

    x = x.ravel()
    y = y.ravel()

    interp = interpolate.griddata((x,y),data.ravel(),(xi,yi),method=method)

    return interp

def createPlot_structure(nrows=1,ncols=3,figsize=(None,None)):
    fig,axes = plt.subplots(ncols=ncols,nrows=nrows,figsize=figsize)

    # creating basemap instances
    m_axes = []
    cbaxes = []
    axes_pos = [
        [0.12,0.25,0.23,0.02],
        [0.40,0.25,0.23,0.02],
        [0.67,0.25,0.23,0.02]
    ]
    for i in np.arange(0,3):
        m = oceano.make_map(axes[i])
        cax = fig.add_axes(axes_pos[i])
        cbaxes.append(cax)
        m_axes.append(m)

    m_axes = np.asarray(m_axes)
    cbaxes = np.asarray(cbaxes)

    return fig,axes,m_axes,cbaxes

def plotData(sat,mod1,mod2,lon,lat,depth,date):
    fig,axes,m_axes,cbaxes = createPlot_structure(nrows=1,ncols=3,figsize=(15,10))

    cf1 = m_axes[0].contourf(lon,lat,sat,np.arange(19.,33.,1.),latlon=True,cmap=cmo.cm.thermal)
    cb1 = plt.colorbar(cf1,cax=cbaxes[0],orientation='horizontal')
    cr1 = m_axes[0].contour(lon,lat,depth,latlon=True,colors=('k'),linestyles=('--'),linewidths=[.4,.4,.4],levels=[100,200,1000])
    plt.clabel(cr1,fontsize=9,inline=1,fmt='%i')

    cf2 = m_axes[1].contourf(lon,lat,mod1,np.arange(15.,33.,1.),latlon=True,cmap=cmo.cm.thermal)
    cb2 = plt.colorbar(cf2,cax=cbaxes[1],orientation='horizontal')
    cr2 = m_axes[1].contour(lon,lat,depth,latlon=True,colors=('k'),linestyles=('--'),linewidths=[.4,.4,.4],levels=[100,200,1000])
    plt.clabel(cr2,fontsize=9,inline=1,fmt='%i')

    cf3 = m_axes[2].contourf(lon,lat,mod2,np.arange(15.,33.,1.),latlon=True,cmap=cmo.cm.thermal)
    cb3 = plt.colorbar(cf3,cax=cbaxes[2],orientation='horizontal')
    cr3 = m_axes[2].contour(lon,lat,depth,latlon=True,colors=('k'),linestyles=('--'),linewidths=[.4,.4,.4],levels=[100,200,1000])
    plt.clabel(cr3,fontsize=9,inline=1,fmt='%i')


    m_axes[0].ax.set_title('GHRSST')
    m_axes[1].ax.set_title('EA1')
    m_axes[2].ax.set_title('EA2')

    plt.suptitle('Sea Surface Temperature [%s]'%(date),y=0.78,fontsize=25)

def plotData_anomalia(sat,mod1,anom,lon,lat,depth,date):
    fig,axes,m_axes,cbaxes = createPlot_structure(nrows=1,ncols=3,figsize=(15,10))

    cf1 = m_axes[0].contourf(lon,lat,sat,np.arange(19.,33.,1.),latlon=True,cmap=cmo.cm.thermal)
    cb1 = plt.colorbar(cf1,cax=cbaxes[0],orientation='horizontal')
    cr1 = m_axes[0].contour(lon,lat,depth,latlon=True,colors=('k'),linestyles=('--'),linewidths=[.4,.4,.4],levels=[100,200,1000])
    plt.clabel(cr1,fontsize=9,inline=1,fmt='%i')

    cf2 = m_axes[1].contourf(lon,lat,mod1,np.arange(15.,33.,1.),latlon=True,cmap=cmo.cm.thermal)
    cb2 = plt.colorbar(cf2,cax=cbaxes[1],orientation='horizontal')
    cr2 = m_axes[1].contour(lon,lat,depth,latlon=True,colors=('k'),linestyles=('--'),linewidths=[.4,.4,.4],levels=[100,200,1000])
    plt.clabel(cr2,fontsize=9,inline=1,fmt='%i')

    cf3 = m_axes[2].contourf(lon,lat,anom,np.arange(-5.,5.,.1),latlon=True,cmap='seismic')
    cb3 = plt.colorbar(cf3,cax=cbaxes[2],orientation='horizontal')
    cr3 = m_axes[2].contour(lon,lat,depth,latlon=True,colors=('k'),linestyles=('--'),linewidths=[.4,.4,.4],levels=[100,200,1000])
    plt.clabel(cr3,fontsize=9,inline=1,fmt='%i')


    m_axes[0].ax.set_title('GHRSST')
    m_axes[1].ax.set_title('EA1')
    m_axes[2].ax.set_title('Anomally')

    plt.suptitle('Sea Surface Temperature [%s]'%(date),y=0.78,fontsize=25)

##############################################################################
#                               MAIN CODE                                    #
##############################################################################

# global variables

BASE_DIR = oceano.make_dir()
if BASE_DIR.split("/")[2] == 'tparente':
    DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output_modelo/exp03_variables/')

else:
    DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
    GHRSST_DIR = BASE_DIR.replace('github/','ventopcse/data/GHRSST/')

FIGU_DIR = BASE_DIR + 'masterThesis_analysis/figures/experiments_outputs/ghrsst_v_ecom/'

fname_ghrsst = GHRSST_DIR + 'ghrsst_JF2014.nc'
fname_exp05  = DATA_DIR + 'EA1.cdf'
fname_exp11  = DATA_DIR + 'EA2.cdf'

# extracting data from GHRSST
sst_sat,time_sat,lon_sat,lat_sat = oceano.load_ghrsst(fname_ghrsst)

# extracting data from model's product
exp05 = Experiment(fname_exp05,timeStart='2014-01-15',timeEnd='2014-02-15',region='pcse')
exp11 = Experiment(fname_exp11,timeStart='2014-01-15',timeEnd='2014-02-15',region='pcse')
# extract sst data
exp05.sst = exp05.ncin['temp'][exp05.timeStart.item():exp05.timeEnd.item(),0,:,:]
exp11.sst = exp11.ncin['temp'][exp11.timeStart.item():exp11.timeEnd.item(),0,:,:]
time_model= exp05.ncin.time[exp11.timeStart.item():exp11.timeEnd.item()].values
exp05.extractVariables()
lon_mod = exp05.lon
lat_mod = exp05.lat

# selecting time indexes to plot
time_sat_begin = 2  # 2014-01-14
time_sat_middle= 19 # 2014-02-01
time_sat_final = 32 # 2014-02-14
time_mod_begin = exp05.findIndex_datetime(pd.to_datetime(time_sat[time_sat_begin]).strftime('%Y-%m-%d %H:%M'),time_model)[0]
time_mod_middle= exp05.findIndex_datetime(pd.to_datetime(time_sat[time_sat_middle]).strftime('%Y-%m-%d %H:%M'),time_model)[0]
time_mod_final = exp05.findIndex_datetime(pd.to_datetime(time_sat[time_sat_final]).strftime('%Y-%m-%d %H:%M'),time_model)[0]

# interpolating satellite data to model grid
sat1 = interp_ghrsst_to_ecom(lon_sat,lat_sat,sst_sat[time_sat_begin],lon_mod,lat_mod)
sat2 = interp_ghrsst_to_ecom(lon_sat,lat_sat,sst_sat[time_sat_middle],lon_mod,lat_mod)
sat3 = interp_ghrsst_to_ecom(lon_sat,lat_sat,sst_sat[time_sat_final],lon_mod,lat_mod)

# import depth data to contour some isobathymetric lines
depth = exp05.ncin.depth.values

##############################################################################
#                            PLOTTING PCSE                                   #
##############################################################################
plt.ion()

plotData(sat1,exp05.sst[time_mod_begin,:,:],exp11.sst[time_mod_begin,:,:],lon_mod,lat_mod,depth,pd.to_datetime(time_sat[time_sat_begin]).strftime('%Y-%m-%d'))
plotData(sat2,exp05.sst[time_mod_middle,:,:],exp11.sst[time_mod_middle,:,:],lon_mod,lat_mod,depth,pd.to_datetime(time_sat[time_sat_middle]).strftime('%Y-%m-%d'))
plotData(sat3,exp05.sst[time_mod_final,:,:],exp11.sst[time_mod_final,:,:],lon_mod,lat_mod,depth,pd.to_datetime(time_sat[time_sat_final]).strftime('%Y-%m-%d'))


# plotData(sat1-sst_t0,exp05.sst[time_mod_begin,:,:]-sst_t0,exp11.sst[time_mod_begin,:,:]-sst_t0,lon_mod,lat_mod,depth,pd.to_datetime(time_sat[time_sat_begin]).strftime('%Y-%m-%d'))

##############################################################################
#                            PLOTTING PCI/PCM                                #
##############################################################################
def masking_array(depth,data):
    data[np.where(depth > 100)] = np.nan
    return data

plt.ion()

# masking all values with depth > 100.
sate = masking_array(depth,sat1)
mod1 = masking_array(depth,exp05.sst[time_mod_begin,:,:])
mod2 = masking_array(depth,exp11.sst[time_mod_begin,:,:])
plotData(sate,mod1,mod2,lon_mod,lat_mod,depth,pd.to_datetime(time_sat[time_sat_begin]).strftime('%Y-%m-%d'))

sate = masking_array(depth,sat2)
mod1 = masking_array(depth,exp05.sst[time_mod_middle,:,:])
mod2 = masking_array(depth,exp11.sst[time_mod_middle,:,:])
plotData(sate,mod1,mod2,lon_mod,lat_mod,depth,pd.to_datetime(time_sat[time_sat_middle]).strftime('%Y-%m-%d'))

sate = masking_array(depth,sat3)
mod1 = masking_array(depth,exp05.sst[time_mod_final,:,:])
mod2 = masking_array(depth,exp11.sst[time_mod_final,:,:])
plotData(sate,mod1,mod2,lon_mod,lat_mod,depth,pd.to_datetime(time_sat[time_sat_final]).strftime('%Y-%m-%d'))

##############################################################################
#                            PLOTTING ANOMALY                                #
##############################################################################
# load seasonal climatology
ncin = xr.open_dataset('/media/danilo/Danilo/mestrado/ventopcse/output/warmupControle.cdf')
clim = ncin.temp[-1,0,:,:].values

anomI = sat3 - clim
anomII= exp05.sst[time_mod_middle,:,:] - clim

lon,lat = lon_mod,lat_mod
fig,axes,m_axes,cbaxes = createPlot_structure(nrows=1,ncols=3,figsize=(15,10))

cf1 = m_axes[0].contourf(lon,lat,anomI,np.arange(-6.,6.,.1),latlon=True,cmap='seismic')
cb1 = plt.colorbar(cf1,cax=cbaxes[0],orientation='horizontal')
cr1 = m_axes[0].contour(lon,lat,depth,latlon=True,colors=('k'),linestyles=('--'),linewidths=[.4,.4,.4],levels=[100,200,1000])
plt.clabel(cr1,fontsize=9,inline=1,fmt='%i')

cf2 = m_axes[1].contourf(lon,lat,anomII,np.arange(-6.,6.,.1),latlon=True,cmap='seismic')
cb2 = plt.colorbar(cf2,cax=cbaxes[1],orientation='horizontal')
cr2 = m_axes[1].contour(lon,lat,depth,latlon=True,colors=('k'),linestyles=('--'),linewidths=[.4,.4,.4],levels=[100,200,1000])
plt.clabel(cr2,fontsize=9,inline=1,fmt='%i')

m_axes[0].ax.set_title('GHRSST - CLIM')
m_axes[1].ax.set_title('PRODUCT - CLIM')



plotData_anomalia(sat3,exp05.sst[time_mod_final,:,:],anom,lon_mod,lat_mod,depth,pd.to_datetime(time_sat[time_sat_final]).strftime('%Y-%m-%d'))
