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
matplotlib.style.use('ggplot')

# import sys
# sys.path.append('/home/danilo/mestrado/github/masterThesis_analysis/routines/masterThesisPack/')
#
# import masterThesisPack as oceano

# functions
def make_map(ax,llat=-30,ulat=-20,llon=-50,ulon=-39,resolution='l',nmeridians=3,nparallels=2):

    m = Basemap(projection='merc', llcrnrlat=llat, urcrnrlat=ulat, llcrnrlon=llon, urcrnrlon=ulon, resolution=resolution)

    m.ax = ax

    m.drawcoastlines(linewidth=0.1)
    m.drawmapboundary()
    m.fillcontinents(color='#c0c0c0')

    m.drawcoastlines(linewidth=.1)
    m.drawmapboundary()

	# definir meridianos e paralelos para plotar no mapa
    meridians=np.arange(llon,ulon,nmeridians)
    parallels=np.arange(llat,ulat,nparallels)
	# desenhar meridianos e paralelos conforme definido acima
    m.drawparallels(parallels,labels=[True,False,False,True],fontsize=13,fontweight='bold',color='gray')
    m.drawmeridians(meridians,labels=[True,False,False,True],fontsize=13,fontweight='bold',color='gray')

    return m

def plotSBB():
    fig,ax = plt.subplot()
    m = make_map(ax,ulat=-23.67,llat=-23.87,ulon=-45.31,llon=-45.50,resolution='f',nmeridians=1,nparallels=1)

    pickle.dump(m,open('/home/danilo/mestrado/github/lapeco2019/birocchi2018/ssbmap.pkl','w'))

    plt.close()

def create_subplots():

    fig,ax = plt.subplots(ncols=2,figsize=(24/2.54,15/2.54))

    # importing already created basemap instace
    mcrude = make_map(ax[0],ulat=-23.67,llat=-23.87,ulon=-45.31,llon=-45.50,resolution='f',nmeridians=1,nparallels=1)#pickle.load(open('/home/danilo/mestrado/github/lapeco2019/birocchi2018/ssbmap.pkl','r'))
    # mcrude.ax = ax[0]
    mrefin = make_map(ax[1],ulat=-23.67,llat=-23.87,ulon=-45.31,llon=-45.50,resolution='f',nmeridians=1,nparallels=1)#pickle.load(open('/home/danilo/mestrado/github/lapeco2019/birocchi2018/ssbmap.pkl','r'))
    # mcrude.ax = ax[1]

    ax[0].set_title('Crude Grid')
    ax[1].set_title('Refined Grid')

    return mcrude,mrefin

def extract_coordinates(ncin):
    lon = ncin.lon.values
    lat = ncin.lat.values
    lon[lon==0.] = np.nan
    lat[lat==0.] = np.nan

    return lon,lat

def estimate_total_concentration(conc,axis=0):
    return np.nansum(conc,axis=axis)

# main program
BASE_DIR = os.getcwd()
SAVE_DIR = BASE_DIR + 'figures/'

# import netcdf files
crude = xr.open_dataset(BASE_DIR + 'crude.cdf')
refin = xr.open_dataset(BASE_DIR + 'refined.cdf')

lon_crude,lat_crude = extract_coordinates(crude)
lon_refin,lat_refin = extract_coordinates(refin)

for t in np.arange(0,crude.time.shape[0],4): # total of 52 figures
    print('Plotting timestep: %i'%t)
    # t = 40
    # extract concentration, calculating the total amount on water column
    conc_crude = estimate_total_concentration(crude.conc[t],axis=0)
    conc_refin = estimate_total_concentration(refin.conc[t],axis=0)

    # auxiliar variables
    dif = np.nanmax(conc_refin) - np.nanmin(conc_refin)
    contours = np.arange(np.nanmin(conc_refin), np.nanmax(conc_refin), dif/100)

    # plotting
    mcrude,mrefin = create_subplots()

    mcrude.contourf(lon_crude,lat_crude,conc_crude,contours,cmap=cmo.cm.matter,latlon=True)
    mrefin.contourf(lon_refin,lat_refin,conc_refin,contours,cmap=cmo.cm.matter,latlon=True)

    # saving
    output = SAVE_DIR + str(t).zfill(5) + '.png'
    plt.savefig(output,dpi=150)

    # cleaning env
    plt.close()
    del conc_crude,conc_refin,mcrude,mrefin
