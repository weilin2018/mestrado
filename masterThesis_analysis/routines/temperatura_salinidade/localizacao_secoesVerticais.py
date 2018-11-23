# add some description here
%reset -f
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

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def make_map(ax,llat=-30,ulat=-20,llon=-50,ulon=-39,resolution='l',nmeridians=3,nparallels=2):

    m = Basemap(projection='merc', llcrnrlat=llat, urcrnrlat=ulat, llcrnrlon=llon, urcrnrlon=ulon, resolution=resolution)

    m.ax = ax

    m.drawcoastlines(linewidth=.1)
    m.drawmapboundary()
    m.fillcontinents(color='#c0c0c0')
	# definir meridianos e paralelos para plotar no mapa
    meridians=np.arange(llon,ulon,nmeridians)
    parallels=np.arange(llat,ulat,nparallels)
	# desenhar meridianos e paralelos conforme definido acima
    m.drawparallels(parallels,labels=[True,False,False,True],fontsize=8,color='gray',linewidth=.4)
    m.drawmeridians(meridians,labels=[True,False,False,True],fontsize=8,color='gray',linewidth=.4)

    return m

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
plt.ion()

# configurações do plot
figsize = (8/2.54, 8/2.54)

DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
fname = glob.glob(DATA_DIR+"EC1.cdf")[0]

ncin = xr.open_dataset(fname)

lon,lat = ncin.lon.values, ncin.lat.values
depth = ncin.depth.values
sigma = ncin.sigma.values
lon[lon == 0.] = np.nan
lat[lat == 0.] = np.nan
depth = ncin.depth.values

fig,ax = plt.subplots(figsize=figsize)

m = make_map(ax,ulat=-22.5,llat=-26.5,llon=-48,ulon=-43.5,nmeridians=3,nparallels=2,resolution='f')
x,y = m(lon,lat)

m.plot(x,y,'k',alpha=.4,linewidth=.3)
m.plot(x.T,y.T,'k',alpha=.4,linewidth=.3)

inor = 99
icen = 28
isul = 19

m.plot(x[inor,:83],y[inor,:83],'r',label='Ubatuba') # 63
m.plot(x[icen,:83],y[icen,:83],'g',label='Santos') # 75
m.plot(x[isul,:83],y[isul,:83],'b',label=u'Cananéia') # 46

plt.legend()

c = m.contour(x,y,depth,levels=[100],linestyles=('dashed'))
plt.clabel(c,inline=1,fmt='%i')


plt.suptitle(u'Localização das Seções Verticais, \nna Plataforma Continental de São Paulo',fontsize=10)
plt.tight_layout()
plt.subplots_adjust(top=0.927,bottom=0.0,left=0.134,right=0.991,hspace=0.2,wspace=0.2)
plt.savefig('/media/danilo/Danilo/mestrado/github/masterThesis_analysis/figures/experiments_outputs/loc_secoesVerticais.png',dpi=300)
