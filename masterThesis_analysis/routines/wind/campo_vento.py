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
import cmocean as cmo

import decomp
from windrose import WindroseAxes

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano
import masterThesisPack.plots as ocplt
##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def make_map(ax,llat=-30,ulat=-20,llon=-50,ulon=-39,resolution='l',nmeridians=3,nparallels=2,labels=[True,False,False,True]):
    # criar mapa sem continente colorido, somente a linha de costa
    # nao colocar meridianos e paralelos
    # nao colocar coordenadas no eixo (sumir com os eixos)

    m = Basemap(projection='merc', llcrnrlat=llat, urcrnrlat=ulat, llcrnrlon=llon, urcrnrlon=ulon, resolution='h')
    # m = pickle.load(open('pickles/basemap.p','r'))
    m.ax = ax
    m.drawcoastlines(linewidth=.2)
    m.fillcontinents(color='white',alpha=0)

    meridians=np.arange(llon,ulon,nmeridians)
    parallels=np.arange(llat,ulat,nparallels)

    return m,meridians,parallels

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
FILE_DIR = BASE_DIR+'masterThesis_analysis/routines/index_list.npy'
os.system('clear')
experiment = raw_input('Digite o experimento controle a ser plotado [EC1, EC2]: ')
fname = DATA_DIR + experiment + '.cdf'
plt.ion()

outputFile = DATA_DIR + 'dados_interpolados_stdl/'

# importing general variables to create figures
ncin = xr.open_dataset(fname)

lon,lat = ncin.lon.values.copy(), ncin.lat.values.copy()
angle = ncin.ang.values.copy()
depth = ncin.depth.values.copy()
sigma = ncin.sigma.values.copy()
lon[lon == 0.] = np.nan
lat[lat == 0.] = np.nan

# which period user wants to visualize
nstepFirst = int(raw_input('Timestep to begin the visualization: '))
howManyOut = int(raw_input('How many plots/days to visualize [maximum is 55 days]: '))

os.system('clear')
print("#####")
animate = raw_input('Animation? y/[N]')


nsteps = np.arange(nstepFirst,howManyOut*8+nstepFirst,8)

if animate == 'y':

    # visualization
    fig,ax = plt.subplots()

    for t in nsteps:

        # extract wind velocity components, calculating the average in time
        wu = np.nanmean(ncin.u[t:t+8,0,:,:],axis=0)
        wv = np.nanmean(ncin.v[t:t+8,0,:,:],axis=0)
        spd= np.sqrt(wu**2 + wv**2)

        # normalize vectors by speed
        wun = wu/spd
        wvn = wv/spd

        # estabilish vectors
        xplot,yplot,wuplot,wvplot = ocplt.formatting_vectors(wun,wvn,lon,lat,FILE_DIR)

        m = oceano.make_map(ax)

        m.contourf(lon,lat,spd,cmap=cmo.cm.speed,latlon=True)
        m.quiver(xplot,yplot,wuplot,wvplot,scale=60,width=0.0015,headwidth=4,headlength=4,latlon=True)
        # plot title
        day = ncin.time[t].values
        day = np.datetime_as_string(day)
        plt.suptitle('Daily average wind field for %s'%(day))

        plt.pause(3)

# calculating spatial average
dailyMean_wu = []
dailyMean_wv = []

for t in nsteps:

    # extract wind velocity components, calculating the average in time
    wu = np.nanmean(ncin.wu[t:t+8,:,:],axis=0)
    wv = np.nanmean(ncin.wv[t:t+8,:,:],axis=0)

    wumean = np.nanmean(wu)
    wvmean = np.nanmean(wv)

    dailyMean_wu.append(wumean)
    dailyMean_wv.append(wvmean)

wumean = np.asarray(dailyMean_wu)
wvmean = np.asarray(dailyMean_wv)

# convert in direction and intensity
ws,wd = decomp.uv2intdir(wumean,wvmean,0,0)

df = pd.DataFrame({'speed':ws,'direction':wd})

# plot directional histogram
ax = WindroseAxes.from_ax()
ax.bar(wd,ws,normed=True, opening=0.8, edgecolor='white')
ax.set_legend()

plt.suptitle('Surface Velocity Daily Spatial Average for %i days simulated in %s'%(howManyOut,experiment))

# windrose.plot_windrose_df(df,kind='')
