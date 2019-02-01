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

    m.ax = ax

    m.drawcoastlines(linewidth=.1)
    m.drawmapboundary()
    m.fillcontinents(color='#c0c0c0')
	# definir meridianos e paralelos para plotar no mapa
    meridians=np.arange(llon,ulon,nmeridians)
    parallels=np.arange(llat,ulat,nparallels)
	# desenhar meridianos e paralelos conforme definido acima
    # m.drawparallels(parallels,labels=[True,False,False,True],fontsize=8,color='gray',linewidth=.4)
    # m.drawmeridians(meridians,labels=[True,False,False,True],fontsize=8,color='gray',linewidth=.4)

    return m

def createPlot_structure(nrows=1,ncols=3,figsize=(None,None)):
    fig,axes = plt.subplots(ncols=ncols,nrows=nrows,figsize=figsize)

    # creating basemap instances
    m_axes = []
    cbaxes = []
    axes_pos = [
        [0.12,0.08,0.36,0.02],
        [0.55,0.08,0.36,0.02],
    ]
    for i in np.arange(0,ncols):
        m = make_map(axes[i])
        cax = fig.add_axes(axes_pos[i])
        cbaxes.append(cax)
        m_axes.append(m)

    m_axes = np.asarray(m_axes)
    cbaxes = np.asarray(cbaxes)

    return fig,axes,m_axes,cbaxes

def plotData_2(sat,mod1,lon,lat,depth,date):
    fig,axes,m_axes,cbaxes = createPlot_structure(nrows=1,ncols=2,figsize=(15/2.54,9/2.54))

    cf1 = m_axes[0].contourf(lon,lat,sat,np.arange(19.,33.,.8),latlon=True,cmap=cmo.cm.thermal)
    cb1 = plt.colorbar(cf1,cax=cbaxes[0],orientation='horizontal',ticks=np.arange(19,33,2))
    cr1 = m_axes[0].contour(lon,lat,depth,latlon=True,colors=('k'),linestyles=('--'),linewidths=[.4,.4],levels=[100,200])
    # plt.clabel(cr1,fontsize=9,inline=1,fmt='%i')
    fig.text(0.23,0.11,r'Temperatura ($^o$C)',fontsize=8)

    cf2 = m_axes[1].contourf(lon,lat,mod1,np.arange(15.,34.,.8),latlon=True,cmap=cmo.cm.thermal)
    cb2 = plt.colorbar(cf2,cax=cbaxes[1],orientation='horizontal',ticks=np.arange(15,34,3))
    cr2 = m_axes[1].contour(lon,lat,depth,latlon=True,colors=('k'),linestyles=('--'),linewidths=[.4,.4],levels=[100,200])
    # plt.clabel(cr2,fontsize=8,inline=1,fmt='%i')
    fig.text(0.65,0.11,r'Temperatura ($^o$C)',fontsize=8)

    m_axes[0].ax.set_title('GHRSST',fontsize=8)
    m_axes[1].ax.set_title('EA1',fontsize=8)

    # plt.suptitle(u'Temperatura da Superfície do Mar em %s'%(date),fontsize=10)
    plt.suptitle(u'Comparação da Temperatura da Superfície do Mar observada \n(GHRSST, à esquerda) com modelado (EA1, à direita), para o dia %s'%(date),fontsize=10)

    rect = (0,0.08,1.,0.95)
    plt.tight_layout(rect=rect) # box for tight_subplot_layout
    plt.subplots_adjust(top=0.935,bottom=0.085,left=0.12,right=0.91,hspace=0.13,wspace=0.2)

    #plt.savefig('/media/danilo/Danilo/mestrado/github/masterThesis_analysis/figures/experiments_outputs/against_ghrsst/%s.png'%(date.replace('/','-')),dpi=300)

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
expEA = xr.open_dataset(fname_EA1)
# extract sst data
sst = expEA.temp[295,0,:,:]
time= expEA.time

lon,lat = expEA.lon.values, expEA.lat.values
depth = expEA.depth.values
sigma = expEA.sigma.values
lon[lon == 0.] = np.nan
lat[lat == 0.] = np.nan
depth = expEA.depth.values
lon_mod = lon
lat_mod = lat

# interpolando dados de satelite para a grade do modelo, para facilitar comparacao
# como eh um processo custoso, checamos se ja existe a variavel no ambiente
if 'sat' not in locals():
    sat = interp_ghrsst_to_ecom(lon_sat,lat_sat,sst_sat[32],lon_mod,lat_mod)

# masking data, where depth is higher than 1000m
maskCondition = np.greater(depth,1000)
sat = np.ma.masked_where(maskCondition,sat)
sst = np.ma.masked_where(maskCondition,sst)

# plotando os conjuntos de dados
plt.ion()
plotData_2(sat,sst[:,:],lon_mod,lat_mod,depth,pd.to_datetime(time_sat[32]).strftime('%d/%m/%Y'))

plt.savefig(FIGU_DIR + 'EA1_ghrsst.eps')
