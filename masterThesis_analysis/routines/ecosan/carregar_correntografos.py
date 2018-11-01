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
def load_data(fname,rotate=False):
    ncin = pd.read_csv(fname,skiprows=11,delimiter=',',usecols=[0,1,2,3,5,12],header=None)
    ncin.columns = ['AVN','AVE','ASPD','AVDIR','datetime','TEMP']

    ncin.index = pd.DatetimeIndex(ncin.datetime.values) # assigning datetimeindex
    ncin.drop('datetime',axis=1,inplace=True) # removing column datetime

    # we have a lot of values for each hour of each day, but it's not a pattern
    # we can resample this timeseries is a hourly frequency
    ncin = ncin.resample('H').mean()

    if rotate:
        # rotating vectors
        angleRot = 55
        magnetic = 24
        u = ncin.AVE.values
        v = ncin.AVN.values
        ncin['ur'],ncin['vr'] = oceano.rot(u,v,angleRot+magnetic)

    return ncin


def load_data_old(fname):

    ncin = pd.read_csv(fname,skiprows=11,delimiter=',',usecols=[2,3,5,12],header=None)
    ncin.columns = ['ASPD','AVDIR','datetime','TEMP']

    ncin.index = pd.DatetimeIndex(ncin.datetime.values) # assigning datetimeindex
    ncin.drop('datetime',axis=1,inplace=True) # removing column datetime

    # we have a lot of values for each hour of each day, but it's not a pattern
    # we can resample this timeseries is a hourly frequency
    ncin = ncin.resample('H').mean()

    return ncin

def quality_control(df):
    return False

def save_netcdf(pathout,df,depth,lat,lon):

    attrs_overall = {
        'Deploy_number': '1360',
        'Mooring': 'ES08',
        'Serial_number': '1360-2D',
        'Maximum_depth': '100m',
        'Latitude': '%s'%(str(lat)),
        'Longitude': '%s'%(str(lon))
    }

    attrs_varirables = {
        'ASPD': {
            'long_name':       'current speed',
            'unit':            'm s-1',
        },
        'AVDIR': {
            'long_name':       'current direction',
            'unit':            'deg',
        },
        'TEMP':     {
            'long_name':       '%sm depth temperature'%(str(depth)),
            'unit':            'degC',
        }
    }

    d = df.to_xarray()

    d.attrs = attrs_overall

    for k in d.keys():
        if k in attrs_varirables.keys():
            d[k].attrs = attrs_varirables[k]

    d.to_netcdf(pathout)

def load_santos100m(DATA_DIR,SAVE_DIR,kind='netcdf',lat=-25.19,lon=-45.82):

    if kind != 'netcdf':
        FUNDEIO  = 'ES08Santos100m/'
        # each file represent a depth, as described in the report 0601
        FILE23m  = DATA_DIR + FUNDEIO + 'ES0801.Dat'
        FILE56m  = DATA_DIR + FUNDEIO + 'ES0802.Dat'
        FILE90m  = DATA_DIR + FUNDEIO + 'ES0803.Dat'

        # reading files
        depth23m = load_data(FILE23m)
        depth56m = load_data(FILE56m)
        depth90m = load_data(FILE90m)

        # saving netcdf files with the data
        save_netcdf(SAVE_DIR+'ecosan_ES0801_23m.nc',depth23m,23,lat,lon)
        save_netcdf(SAVE_DIR+'ecosan_ES0802_56m.nc',depth56m,56,lat,lon)
        save_netcdf(SAVE_DIR+'ecosan_ES0803_90m.nc',depth90m,90,lat,lon)
    else:
        depth23m = xr.open_dataset(SAVE_DIR+'ecosan_ES0801_23m.nc')
        depth56m = xr.open_dataset(SAVE_DIR+'ecosan_ES0802_56m.nc')
        depth90m = xr.open_dataset(SAVE_DIR+'ecosan_ES0803_90m.nc')

        depth23m = depth23m.to_dataframe()
        depth56m = depth56m.to_dataframe()
        depth90m = depth90m.to_dataframe()

    # plotting temperature
    fig,ax = plt.subplots(nrows=3,figsize=(15,15))

    depth23m['2006'].TEMP.plot(c='#c0c0c0',label='_',ax=ax[0])
    depth23m['2006'].resample('D').mean().TEMP.plot(ax=ax[0],c='k',label='Temperature at 23m depth')
    ax[0].set_ylim([11,27])
    ax[0].legend()

    depth56m['2006'].TEMP.plot(c='#c0c0c0',label='_',ax=ax[1])
    depth56m['2006'].resample('D').mean().TEMP.plot(ax=ax[1],c='k',label='Temperature at 56m depth')
    ax[1].set_ylim([11,27])
    ax[1].legend()

    depth90m['2006'].TEMP.plot(c='#c0c0c0',label='_',ax=ax[2])
    depth90m['2006'].resample('D').mean().TEMP.plot(ax=ax[2],c='k',label='Temperature at 90m depth')
    ax[2].set_ylim([11,27])
    ax[2].legend()

    plt.suptitle('Daily Temperature in the 2006 Summer \nSantos Mooring located at [-25.084, -45.7]',fontsize=24)

    return depth23m,depth56m,depth90m

def load_peruibe20m(DATA_DIR,SAVE_DIR,kind='netcdf',lat=-24.406,lon=-46.89):

    if kind != 'netcdf':
        FUNDEIO  = 'ES09Peruibe/'
        # each file represent a depth, as described in the report 0601
        FILE3m  = DATA_DIR + FUNDEIO + 'ES0901.Dat'
        FILE12m  = DATA_DIR + FUNDEIO + 'ES0902.Dat'


        # reading files
        depth3m = load_data(FILE3m)
        depth12m = load_data(FILE12m)

        # if everything are ok, save:
        # saving netcdf files with the data
        save_netcdf(SAVE_DIR+'ecosan_ES0901_03m.nc',depth3m,3,lat,lon)
        save_netcdf(SAVE_DIR+'ecosan_ES0902_12m.nc',depth12m,12,lat,lon)
    else:
        depth3m = xr.open_dataset(SAVE_DIR+'ecosan_ES0901_03m.nc')
        depth12m = xr.open_dataset(SAVE_DIR+'ecosan_ES0902_12m.nc')

        depth3m  = depth3m.to_dataframe()
        depth12m = depth12m.to_dataframe()

    # plotting temperature
    fig,ax = plt.subplots(nrows=2,figsize=(15,15))

    depth3m['2006'].TEMP.plot(c='#c0c0c0',label='_',ax=ax[0])
    depth3m['2006'].resample('D').mean().TEMP.plot(ax=ax[0],c='k',label='Temperature at 3m depth')
    ax[0].set_ylim([18,30])
    ax[0].legend()

    depth12m['2006'].TEMP.plot(c='#c0c0c0',label='_',ax=ax[1])
    depth12m['2006'].resample('D').mean().TEMP.plot(ax=ax[1],c='k',label='Temperature at 12m depth')
    ax[1].set_ylim([18,30])
    ax[1].legend()

    plt.suptitle('Daily Temperature in the 2006 Summer \nPeruibe Mooring located at [-24.4, -46.9]',fontsize=24)

    return depth3m,depth12m

def load_montetrigo20m(DATA_DIR,SAVE_DIR,kind='netcdf',lat=-23.845,lon=-45.66):

    if kind != 'netcdf':
        FUNDEIO = 'ES10MonteTrigo/'
        # each file represent a depth, as described in the report 0601
        FILE3m  = DATA_DIR + FUNDEIO + 'ES1001.Dat'
        FILE12m  = DATA_DIR + FUNDEIO + 'ES1002.Dat'

        # reading files
        depth3m = load_data(FILE3m)
        depth12m = load_data(FILE12m)

        # if everything are ok, save:
        # saving netcdf files with the data
        save_netcdf(SAVE_DIR+'ecosan_ES1001_03m.nc',depth3m,3,lat,lon)
        save_netcdf(SAVE_DIR+'ecosan_ES1002_12m.nc',depth12m,12,lat,lon)
    else:
        depth3m = xr.open_dataset(SAVE_DIR+'ecosan_ES1001_03m.nc')
        depth12m = xr.open_dataset(SAVE_DIR+'ecosan_ES1002_12m.nc')

        depth3m  = depth3m.to_dataframe()
        depth12m = depth12m.to_dataframe()

    # plotting temperature
    fig,ax = plt.subplots(nrows=2,figsize=(15,15))

    depth3m['2006'].TEMP.plot(c='#c0c0c0',label='_',ax=ax[0])
    depth3m['2006'].resample('D').mean().TEMP.plot(ax=ax[0],c='k',label='Temperature at 3m depth')
    ax[0].set_ylim([16,30])
    ax[0].legend()

    depth12m['2006'].TEMP.plot(c='#c0c0c0',label='_',ax=ax[1])
    depth12m['2006'].resample('D').mean().TEMP.plot(ax=ax[1],c='k',label='Temperature at 12m depth')
    ax[1].set_ylim([16,30])
    ax[1].legend()

    plt.suptitle('Daily Temperature in the 2006 Summer \nMonte Trigo Mooring located at [-23.84, -45.67]',fontsize=24)

    return depth3m,depth12m

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
DATA_DIR = '/media/danilo/Danilo/mestrado/ventopcse/data/ECOSAN/EcoSanFundeios_Dados_Out2005Fev2006/'
SAVE_DIR = '/media/danilo/Danilo/mestrado/ventopcse/data/ECOSAN/'

depth23m,depth56m,depth90m      = load_santos100m(DATA_DIR,SAVE_DIR,kind='other')
depth3m_per,depth12m_per        = load_peruibe20m(DATA_DIR,SAVE_DIR,kind='other')
depth3m_mont,depth12m_mont      = load_montetrigo20m(DATA_DIR,SAVE_DIR,kind='other')
