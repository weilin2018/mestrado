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
import gsw

# pacotes para minimap
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

import matplotlib
matplotlib.use('PS')
# matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here

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

    interp = griddata((x,y),data.ravel(),(xi,yi),method=method)

    return interp

def make_map(ax,llat=-30,ulat=-20,llon=-50,ulon=-39,resolution='l',nmeridians=3,nparallels=2):

    m = Basemap(projection='merc', llcrnrlat=llat, urcrnrlat=ulat, llcrnrlon=llon, urcrnrlon=ulon, resolution='h')

    # m = pickle.load(open('pickles/basemap.p','r'))
    m.ax = ax
    m.drawcoastlines(linewidth=.2)
    m.fillcontinents(color='white',alpha=0)

    return m

def createPlot_structure(lon,lat,nrows=1,ncols=3,figsize=(None,None)):
    fig,axes = plt.subplots(ncols=ncols,nrows=nrows,figsize=figsize)

    # creating basemap instances
    m_axes = []
    cbaxes = []
    axes_pos = [
        [0.03,0.1,0.29,0.02],
        [0.36,0.1,0.29,0.02],
        [0.68,0.1,0.29,0.02]
    ]
    for i in np.arange(0,ncols):
        m = make_map(axes[i],ulon=np.nanmax(lon)-.5,llon=np.nanmin(lon)-.2,ulat=np.nanmax(lat)+.2,llat=np.nanmin(lat))
        cax = fig.add_axes(axes_pos[i])
        cbaxes.append(cax)
        m_axes.append(m)

    m_axes = np.asarray(m_axes)
    cbaxes = np.asarray(cbaxes)

    return fig,axes,m_axes,cbaxes

def plotData(sat,mod1,mod2,lon,lat,depth,date):


    fig,axes,m_axes,cbaxes = createPlot_structure(lon,lat,figsize=(15/2.54,6.5/2.54))

    cf1 = m_axes[0].contourf(lon,lat,sat,np.arange(19.,33.,.8),latlon=True,cmap=cmo.cm.thermal,rasterized=True,extend='max')
    cb1 = plt.colorbar(cf1,cax=cbaxes[0],orientation='horizontal',ticks=np.arange(19,33,2))
    cr1 = m_axes[0].contour(lon,lat,depth,latlon=True,colors=('k'),linestyles=('--'),linewidths=[.4,.4],levels=[40,90])
    # plt.clabel(cr1,fontsize=9,inline=1,fmt='%i')
    # fig.text(0.23,0.11,r'Temperatura ($^o$C)',fontsize=8)

    cf2 = m_axes[1].contourf(lon,lat,mod1,np.arange(15.,30.,.8),latlon=True,cmap=cmo.cm.thermal,rasterized=True,extend='max')
    cb2 = plt.colorbar(cf2,cax=cbaxes[1],orientation='horizontal',ticks=np.arange(15,34,3))
    cr2 = m_axes[1].contour(lon,lat,depth,latlon=True,colors=('k'),linestyles=('--'),linewidths=[.4,.4],levels=[40,90])
    # plt.clabel(cr2,fontsize=8,inline=1,fmt='%i')
    # fig.text(0.65,0.11,r'Temperatura ($^o$C)',fontsize=8)

    cf3 = m_axes[2].contourf(lon,lat,mod2,np.arange(15.,30.,.8),latlon=True,cmap=cmo.cm.thermal,rasterized=True,extend='max')
    cb3 = plt.colorbar(cf3,cax=cbaxes[2],orientation='horizontal',ticks=np.arange(15,34,3))
    cr3 = m_axes[2].contour(lon,lat,depth,latlon=True,colors=('k'),linestyles=('--'),linewidths=[.4,.4],levels=[40,90])
    # plt.clabel(cr2,fontsize=8,inline=1,fmt='%i')
    # fig.text(0.65,0.11,r'Temperatura ($^o$C)',fontsize=8)

    # matplotib trick to remove white thin lines when saving contourf in pdf
    for c in cf1.collections:
        c.set_edgecolor('face')
        c.set_linewidth(0.00000000001)

    for c in cf2.collections:
        c.set_edgecolor('face')
        c.set_linewidth(0.00000000001)

    for c in cf3.collections:
        c.set_edgecolor('face')
        c.set_linewidth(0.00000000001)


    # setting colorbar tick labels
    from matplotlib import ticker
    tick_locator = ticker.MaxNLocator(nbins=5)
    cb1.locator = tick_locator
    cb1.update_ticks()
    cb1.ax.axes.tick_params(axis='both',which='both',labelsize=8)
    cb1.ax.set_title(r'Temperatura ($^o$C)',fontsize=8)

    cb2.locator = tick_locator
    cb2.update_ticks()
    cb2.ax.axes.tick_params(axis='both',which='both',labelsize=8)
    cb2.ax.set_title(r'Temperatura ($^o$C)',fontsize=8)

    cb3.locator = tick_locator
    cb3.update_ticks()
    cb3.ax.axes.tick_params(axis='both',which='both',labelsize=8)
    cb3.ax.set_title(r'Temperatura ($^o$C)',fontsize=8)

    plt.tight_layout()
    plt.subplots_adjust(top=0.941,bottom=0.124,left=0.031,right=0.969,hspace=0.2,wspace=0.13)

##############################################################################
#                               MAIN CODE                                    #
##############################################################################

# global variables
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/output')
GHRSST_DIR = BASE_DIR.replace('github/','ventopcse/data/GHRSST/')

FIGU_DIR = BASE_DIR + 'masterThesis_analysis/figures/experiments_outputs/ghrsst_v_ecom/'

fname_ghrsst = GHRSST_DIR + 'ghrsst_summer2014.nc'
fname_EA1  = DATA_DIR + 'EA1.cdf'
fname_EA2  = DATA_DIR + 'EA2.cdf'

# extracting data from GHRSST
sst_sat,time_sat,lon_sat,lat_sat = oceano.load_ghrsst(fname_ghrsst)

# extracting data from model's product
expEA1 = xr.open_dataset(fname_EA1)
expEA2 = xr.open_dataset(fname_EA2)
# extract sst data
sst1= expEA1.temp[295,0,:,:]
time= expEA1.time

sst2= expEA2.temp[295,0,:,:]

lon,lat = expEA1.lon.values, expEA1.lat.values
depth = expEA1.depth.values
sigma = expEA1.sigma.values
lon[lon == 0.] = np.nan
lat[lat == 0.] = np.nan
depth = expEA1.depth.values
lon_mod = lon
lat_mod = lat

# interpolando dados de satelite para a grade do modelo, para facilitar comparacao
# como eh um processo custoso, checamos se ja existe a variavel no ambiente
if 'sat' not in locals():
    sat = interp_ghrsst_to_ecom(lon_sat,lat_sat,sst_sat[32],lon_mod,lat_mod)

# masking data, where depth is higher than 1000m
maskCondition = np.greater(depth,100)
sat = np.ma.masked_where(maskCondition,sat)
sst1= np.ma.masked_where(maskCondition,sst1)
sst2= np.ma.masked_where(maskCondition,sst2)

# plotando os conjuntos de dados
plt.ion()

plotData(sat,sst1,sst2,lon,lat,depth,'0')

##### TODO 
def plotDiferenca(sat,mod1,mod2,lon,lat,depth,date):

    # calcular diferenca (GHRSST - ECOM)
    diff1 = sat - mod1
    diff2 = sat - mod2



    fig,axes,m_axes,cbaxes = createPlot_structure(lon,lat,figsize=(15/2.54,6.5/2.54))

    cf1 = m_axes[0].contourf(lon,lat,sat,np.arange(19.,33.,.8),latlon=True,cmap=cmo.cm.thermal,rasterized=True,extend='max')
    cb1 = plt.colorbar(cf1,cax=cbaxes[0],orientation='horizontal',ticks=np.arange(19,33,2))
    cr1 = m_axes[0].contour(lon,lat,depth,latlon=True,colors=('k'),linestyles=('--'),linewidths=[.4,.4],levels=[40,90])
    # plt.clabel(cr1,fontsize=9,inline=1,fmt='%i')
    # fig.text(0.23,0.11,r'Temperatura ($^o$C)',fontsize=8)

    cf2 = m_axes[1].contourf(lon,lat,mod1,np.arange(15.,30.,.8),latlon=True,cmap=cmo.cm.thermal,rasterized=True,extend='max')
    cb2 = plt.colorbar(cf2,cax=cbaxes[1],orientation='horizontal',ticks=np.arange(15,34,3))
    cr2 = m_axes[1].contour(lon,lat,depth,latlon=True,colors=('k'),linestyles=('--'),linewidths=[.4,.4],levels=[40,90])
    # plt.clabel(cr2,fontsize=8,inline=1,fmt='%i')
    # fig.text(0.65,0.11,r'Temperatura ($^o$C)',fontsize=8)

    cf3 = m_axes[2].contourf(lon,lat,mod2,np.arange(15.,30.,.8),latlon=True,cmap=cmo.cm.thermal,rasterized=True,extend='max')
    cb3 = plt.colorbar(cf3,cax=cbaxes[2],orientation='horizontal',ticks=np.arange(15,34,3))
    cr3 = m_axes[2].contour(lon,lat,depth,latlon=True,colors=('k'),linestyles=('--'),linewidths=[.4,.4],levels=[40,90])
    # plt.clabel(cr2,fontsize=8,inline=1,fmt='%i')
    # fig.text(0.65,0.11,r'Temperatura ($^o$C)',fontsize=8)

    # matplotib trick to remove white thin lines when saving contourf in pdf
    for c in cf1.collections:
        c.set_edgecolor('face')
        c.set_linewidth(0.00000000001)

    for c in cf2.collections:
        c.set_edgecolor('face')
        c.set_linewidth(0.00000000001)

    for c in cf3.collections:
        c.set_edgecolor('face')
        c.set_linewidth(0.00000000001)


    # setting colorbar tick labels
    from matplotlib import ticker
    tick_locator = ticker.MaxNLocator(nbins=5)
    cb1.locator = tick_locator
    cb1.update_ticks()
    cb1.ax.axes.tick_params(axis='both',which='both',labelsize=8)
    cb1.ax.set_title(r'Temperatura ($^o$C)',fontsize=8)

    cb2.locator = tick_locator
    cb2.update_ticks()
    cb2.ax.axes.tick_params(axis='both',which='both',labelsize=8)
    cb2.ax.set_title(r'Temperatura ($^o$C)',fontsize=8)

    cb3.locator = tick_locator
    cb3.update_ticks()
    cb3.ax.axes.tick_params(axis='both',which='both',labelsize=8)
    cb3.ax.set_title(r'Temperatura ($^o$C)',fontsize=8)

    plt.tight_layout()
    plt.subplots_adjust(top=0.941,bottom=0.124,left=0.031,right=0.969,hspace=0.2,wspace=0.13)

# plotData_2(sat,sst[:,:],lon_mod,lat_mod,depth,pd.to_datetime(time_sat[32]).strftime('%d/%m/%Y'),colorbar='EA2')

# plt.savefig(FIGU_DIR + 'EA1_ghrsst.eps')
