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
import decomp # pacote da Carine

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def create_dates(ncin):

    dates = []
    years = ncin[:,0]
    month = ncin[:,1]
    days  = ncin[:,2]
    hours = ncin[:,3]

    for y,m,d,h in zip(years,month,days,hours):
        dates.append(datetime.datetime(int(y),int(m),int(d),int(h)))

    return np.asarray(dates)


# class Ecosan():

    # def __init__(self):


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = '/media/danilo/Danilo/mestrado/ventopcse/data/ECOSAN/EcoSanFundeios_Dados_Out2005Fev2006/EstMetLaje/'
FILENAME1 = BASE_DIR+'lage_0106.txt'
FILENAME2 = BASE_DIR+'lage_0206.txt'

# loading Jan
ncin = np.loadtxt(FILENAME1,skiprows=4,usecols=[0,1,2,3,5,6,8,9,10])
dates = create_dates(ncin)
# create dataframe
df1 = pd.DataFrame(ncin,index=dates)

# loading Feb
ncin = np.loadtxt(FILENAME2,skiprows=4,usecols=[0,1,2,3,5,6,8,9,10])
dates = create_dates(ncin)
df2 = pd.DataFrame(ncin,index=dates)

# concate both dfs
dfs = [df1,df2]
df = pd.concat(dfs)

# naming columns
df.columns = ['years','months','days','hours','lat','lon','speed','direction','temperature']
# remove some columns
df.drop(['years','months','days','hours'],axis=1,inplace=True)

# replace rows with zeroes by nan
df[df.lat == 0.] = np.nan
# basic quality control
df[df.speed > 3*df.speed.std()] = np.nan
df[df.speed <-3*df.speed.std()] = np.nan
df[df.temperature > 40] = np.nan

# convert intensity,direction into along and cross shore components,
# by rotating vectors in 54º and correcting magnetic declination -21.09º
along,cross = decomp.intdir2uv(df.speed,df.direction,-21.09,54.)
# df['along'] = along
# df['cross'] = cross

# Laje de Santos Station build in meteorological convention and, for a stickplot,
# is better to visualize wind in an oceanographic convention. So ...
df['wu'] = along * -1
df['wv'] = cross * -1

# check stickplots
fig,ax = plt.subplots()
oceano.stickplot(df,ax)

# check temperature
df.temperature.interpolate(method='linear',inplace=True)

fig,ax = plt.subplots()
df.temperature.plot(ax=ax)
ax.axhline(y=df.temperature.mean(),c='k',alpha=.4)
ax.axhline(y=df.temperature.std()*3+df.temperature.mean(),c='k',linestyle='--',alpha=.4)
ax.axhline(y=-df.temperature.std()*3+df.temperature.mean(),c='k',linestyle='--',alpha=.4)


# converting dataframe into xarray and saving into a netcdf file
pathout= '/media/danilo/Danilo/mestrado/ventopcse/data/LajedeSantos_ECOSAN_2006.nc'
netcdf = df.to_xarray()

# inserting attributes for each variable
netcdf.speed.attrs     = {'units':u'm s-1','long_name':u'wind_intensity'}
netcdf.direction.attrs = {'units':u'degrees','long_name':u'wind_direction'}
netcdf.temperature.attrs={'units':u'degrees C','long_name':u'air_temperature'}
netcdf.wu.attrs         ={'units':u'm s-1','long_name':u'along_shore_wind_component'}
netcdf.wv.attrs         ={'units':u'm s-1','long_name':u'cross_shore_wind_component'}

netcdf.attrs = {
    'Quality Control': u'Removed all values outside 3 standard deviations interval',
    'Statistical Treatment': u'along and cross shore componentes corrected for magnetic declination and rotated in 54oW',
    'Informations': u'Data available in ECOSAN dataset, for Jan and Feb, 2006'
}

netcdf.to_netcdf(pathout)

# to read:
ncin = xr.open_dataset(pathout)
df = ncin.to_dataframe()
