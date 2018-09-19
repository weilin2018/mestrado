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
import glob
import matplotlib
from datetime import date

matplotlib.style.use('ggplot')

import marineHeatWaves as mhw

import sys
sys.path.append('../masterThesisPack/')

import masterThesisPack as oceano

plt.ion()

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

# define which variable extract from netcdf files
var = 'POT_L160_Avg_1'

#################### data from CFSR
data_dir = '/media/danilo/Danilo/mestrado/ventopcse/HeatBudget/data/'
nfiles = glob.glob(data_dir+"*.nc")
nfiles.sort()

sst_cfsr = []

for file_i,fname in enumerate(nfiles):
    ncin = xr.open_dataset(fname)
    data = ncin[var].values

    # create array to store daily means
    # dailySST = np.zeros(len(ncin.time)/4)
    # cont = 0

    for i in np.arange(0,len(ncin.time),4):
        sst_cfsr.append(np.nanmean(data[i:i+3,0,0]) - 273.15)

sst_cfsr = np.asarray(sst_cfsr)

# creating time array
dt   = pd.date_range(start='1982-01-01',end='2010-12-31',freq='1D')

cfsr = pd.DataFrame({'sst':sst_cfsr},index=dt)

###########################################################################
# alguns arquivos CFSv2 diferenciados precisam ser tratados de outra forma
##########################################################################
data_dir = '/media/danilo/Danilo/mestrado/ventopcse/HeatBudget/begin/'

# loading files
nfiles = glob.glob(data_dir+"*.nc")
nfiles.sort()

# loading data
sst_begin = []

for file_i,fname in enumerate(nfiles):
    ncin = xr.open_dataset(fname)
    data = ncin[var].values

    for i in np.arange(0,len(ncin.time),4):
        sst_begin.append(np.nanmean(data[i:i+3,0,0]) - 273.15)

sst_begin = np.asarray(sst_begin)

#################### data from CFSv2
data_dir = '/media/danilo/Danilo/mestrado/ventopcse/HeatBudget/'

# loading files
nfiles = glob.glob(data_dir+"*.nc")
nfiles.sort()

# loading data
sst = np.zeros(len(nfiles))

for file_i,fname in enumerate(nfiles):
    ncin  = xr.open_dataset(fname)
    data  = ncin[var].values

    # calculating daily mean
    daily = np.nanmean(data)

    # converting from Kelvin to Celsius
    sst[file_i] = np.squeeze(daily) - 273.15

# creating time array
dt = pd.date_range(start='2011-04-01',end='2018-09-13',freq='1D')

cfsv2 = pd.DataFrame({'sst':sst},index=dt)

#################### joining two dataset for a long serie
sst_total = []

for i in sst_cfsr:
    sst_total.append(i)

for i in sst_begin:
    sst_total.append(i)

for i in sst:
    sst_total.append(i)

sst_total = np.asarray(sst_total)

dt = pd.date_range(start='1982-01-01',end='2018-09-13',freq='1D')

sst = pd.DataFrame({'sst': sst_total},index=dt)

###############################################################################
###############################################################################
#                           CALCULATING HEAT BUDGET                           #
###############################################################################
###############################################################################
t = np.arange(date(1982,1,1).toordinal(),date(2014,12,31).toordinal()+1)
dates = [date.fromordinal(tt.astype(int)) for tt in t]

mhws, clim = mhw.detect(t, sst,climatologyPeriod=[1982,2012])
