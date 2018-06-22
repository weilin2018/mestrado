"""
    Funcao para visualizar os dados dos periodos selecionados como períodos para a  
    simulação controle:

    Anos 2007, 2008 e 2010.
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
def load_data(fname):

	# extracting data from January
	ncdata = xr.open_dataset(fname[0])

	lon,lat = ncdata['lon'].values,ncdata['lat'].values
	# correcting longitude data from 0/360 to -180/180
	lon -= 360
	# gridding lon and lat
	lon,lat = np.meshgrid(lon,lat)
	time_jan = ncdata['time'].values

	wu_jan = ncdata['U_GRD_L103'].values
	wv_jan = ncdata['V_GRD_L103'].values

	# extracting data from February
	ncdata = xr.open_dataset(fname[1])
	time_feb = ncdata['time'].values

	wu_feb = ncdata['U_GRD_L103'].values
	wv_feb = ncdata['V_GRD_L103'].values

	# concatenating
	wu = np.vstack((wu_jan,wu_feb))
	wv = np.vstack((wv_jan,wv_feb))
	time = np.concatenate((time_jan,time_feb))

	return lon,lat,time,wu,wv

def createDataFrame(wu,wv,time):

	dct = {'wu':wu,'wv':wv}

	return pd.DataFrame(dct,index=pd.DatetimeIndex(time))

def plot_windField(fname,debug=True,savefig=None):

	# extract all data
	lon,lat,time,wu,wv = load_data(fname)

	# calculate speed
	spd = np.sqrt(wu*wu + wv*wv)
	# normalizing vectors
	wun = wu/spd
	wvn = wv/spd

	# testing for debug
	if debug:
		endTime = len(time)-10
	else:
		endTime = 0

	if not savefig:
		# animation
		plt.ion()

		# creating figure to plot
		fig,ax = plt.subplots()
		# create axis for colorbar
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("bottom", size="10%",pad=0.05)

		for i in np.arange(0,time.shape[0]-endTime):
			ax.clear()
			cax.clear()
			cax.set_xlim([0,240])
			cax.set_ylim([-16,18])
			cax.axhline(xmin=0,xmax=240,color='k',alpha=0.5)
			m = oceano.make_map(ax, resolution='i')
			csf = m.contourf(lon,lat,spd[i,:,:],latlon=True,cmap=cmo.cm.speed)
			q = m.quiver(lon,lat,wun[i,:,:],wvn[i,:,:],latlon=True)
			# plt.clabel(cs, fmt='%2.1f',colors='k',fontsize=14)
			ts = pd.to_datetime(str(time[i]))
			plt.suptitle(ts.strftime('%Y.%m.%d %H:%M'))

			m.scatter(lon[15,15],lat[15,15],latlon=True,marker='o',s=100)

			if i > 3:
				cax.plot(wu[:i,15,15],color='#12BE0C',label='WU')
				cax.plot(wv[:i,15,15],color='#4768B0',label='WV')

				cax.legend(loc='upper right')
			elif i > 230:

				cax.plot(wu[:i,15,15],color='#12BE0C',label='WU')
				cax.plot(wv[:i,15,15],color='#4768B0',label='WV')

				cax.legend(loc='upper left')

			plt.pause(0.1)


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

BASE_DIR = oceano.make_dir()
if BASE_DIR.split("/")[2] == 'tparente':
	DATA_DIR = '/home/tparente/Dropbox/mestrado/data/control_experiment/selected_periods/'
else:
	DATA_DIR = '/home/danilo/Dropbox/mestrado/data/control_experiment/selected_periods/'

# clear screen
os.system("clear")

# asking which period to view
year = input("Select which period to view: 2007, 2008 or 2010: ")

# reading files in directory
fname = glob.glob(DATA_DIR+str(year)+"/*.nc")

# time to plot
