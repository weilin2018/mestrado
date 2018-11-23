# add some description here

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

    m.drawcoastlines(linewidth=0.5)
    m.drawmapboundary()
    m.drawstates(linewidth=0.2)
    m.drawcountries(linewidth=0.2)
    # m.drawrivers(linewidth=.1,color='#c0c0c0')
    # m.fillcontinents(color='#c0c0c0')

    # m.drawcoastlines(linewidth=.1)
    m.drawmapboundary()

	# definir meridianos e paralelos para plotar no mapa
    meridians=np.arange(llon,ulon,nmeridians)
    parallels=np.arange(llat,ulat,nparallels)
	# desenhar meridianos e paralelos conforme definido acima
    m.drawparallels(parallels,labels=[True,False,False,True],fontsize=8,fontweight='bold',color='gray',fmt='%2.2f')
    m.drawmeridians(meridians,labels=[True,False,False,True],fontsize=8,fontweight='bold',color='gray',fmt='%2.2f')

    return m

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
DATA_DIR = '/media/danilo/Danilo/RegCM/2nd_Workshop/RCMDATA/GPCP/'

flist = glob.glob(DATA_DIR + '*.nc')
flist.sort()

ncin = xr.open_dataset(flist[-1])

# recorte: lon: [-57.157,-35.522]
# recorte: lat: [-32.108,-15.072]
