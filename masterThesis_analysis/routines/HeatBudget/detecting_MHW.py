"""

TODO:: adaptar codigo para pegar um ponto de interesse
na grade do GHRSST e extrair os dados neste ponto e realizar
os calculos.

Objetivo:

    . usando dados do GHRSST em um único ponto na porcao Sul da South Brazil Bight,
        verificar se, em 2014 e 2015, houveram episodios de marine heatwaves,
        pelo metodo usado por Oliver et al. (2017).

    . testar o mesmo procedimento para medicoes in situ, se conseguir

    . periodo de analise: verao2014, verao2015, sendo verao definido como
    periodo de Dez a Mar.

    Dataset: NCEP/NCAR (CFSR e CFSv2)
    Localizacao: Ponto mais proximo a Cananeia

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
def deseason_harmonic(dat, K, L):
    '''
        deseasoned_data, season, beta = deseason_harmonic(dat, K, L)
        Subtracts the seasonal cycle (season) from the data (data). Season
        is calculated by fitting K harmonics of the annual cycle (as well as the
        mean) to the data. Assumes on year has L elements (i.e., 365 for daily data,
        73 for pentad data, 52 for weekly data, etc.).
        Outputs the deseasonalized data, the season, and the fitting parameters (beta)
        Handles missing values as np.nan's

        Written by Eric Oliver, Dalhousie University, 2007-2011
        Adapted from original MATLAB script on 28 November 2012
    '''

#   Valid (non-NaN) elements
    valid = ~np.isnan(dat)

#   ensure that mat is a matrix and a column vector
    dat = np.mat(dat)
    if dat.shape[1]!=0:
        dat = dat.T

#   length of time series and generate time vector
    n = len(dat)
    time = np.mat(np.arange(1,n+1)/(1.*L))

#   set up mean and harmonics to fit data
    P = np.mat(np.ones((n,1)))
    for k in range(1,K+1):
        P = np.concatenate((P, np.cos(k*2*np.pi*time.T)), 1)
        P = np.concatenate((P, np.sin(k*2*np.pi*time.T)), 1)

#   Remove seasonal cycle by harmonic regression
    beta = (np.linalg.inv(P[valid,:].T*P[valid,:])*P[valid,:].T)*dat[valid]
    season = P*beta
    dat_ds = dat - season

    return dat_ds, season, beta

# visualizing
def visualizing_MHW(mhws,clim,ev=ev):
    plt.figure(figsize=(14,10))
    plt.subplot(2,1,1)
    # Plot SST, seasonal cycle, and threshold
    plt.plot(dates, sst, 'k-')
    plt.plot(dates, clim['thresh'], 'g-')
    plt.plot(dates, clim['seas'], 'b-')
    plt.title('SST (black), seasonal climatology (blue), \
              threshold (green), detected MHW events (shading)')
    plt.xlim(t[0], t[-1])
    plt.ylim(sst.min()-0.5, sst.max()+0.5)
    plt.ylabel(r'SST [$^\circ$C]')
    plt.subplot(2,1,2)
    # Find indices for all ten MHWs before and after event of interest and shade accordingly
    for ev0 in np.arange(ev-2, ev+2, 1):
        t1 = np.where(t==mhws['time_start'][ev0])[0][0]
        t2 = np.where(t==mhws['time_end'][ev0])[0][0]
        plt.fill_between(dates[t1:t2+1], sst[t1:t2+1], clim['thresh'][t1:t2+1], \
                         color=(1,0.6,0.5))
    # Find indices for MHW of interest and shade accordingly
    t1 = np.where(t==mhws['time_start'][ev])[0][0]
    t2 = np.where(t==mhws['time_end'][ev])[0][0]
    plt.fill_between(dates[t1:t2+1], sst[t1:t2+1], clim['thresh'][t1:t2+1], \
                     color='r')
    # Plot SST, seasonal cycle, threshold, shade MHWs with main event in red
    plt.plot(dates, sst, 'k-', linewidth=2)
    plt.plot(dates, clim['thresh'], 'g-', linewidth=2)
    plt.plot(dates, clim['seas'], 'b-', linewidth=2)
    plt.title('SST (black), seasonal climatology (blue), \
              threshold (green), detected MHW events (shading)')
    plt.xlim(mhws['time_start'][ev]-150, mhws['time_end'][ev]+150)
    plt.ylim(clim['seas'].min() - 1, clim['seas'].max() + mhws['intensity_max'][ev] + 0.5)
    plt.ylabel(r'SST [$^\circ$C]')

def describe(mhws,clim):
    os.system('clear')

    ev = np.argmax(mhws['intensity_max']) # Find largest event
    print 'Maximum intensity:', mhws['intensity_max'][ev], 'deg. C'
    print 'Average intensity:', mhws['intensity_mean'][ev], 'deg. C'
    print 'Cumulative intensity:', mhws['intensity_cumulative'][ev], 'deg. C-days'
    print 'Duration:', mhws['duration'][ev], 'days'
    print 'Start date:', mhws['date_start'][ev].strftime("%d %B %Y")
    print 'End date:', mhws['date_end'][ev].strftime("%d %B %Y")

    return ev

def load_sst(dataset,path):

    if dataset == 'NCEP':
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

        df_sst = pd.DataFrame({'sst': sst_total},index=dt)


    elif dataset == 'GHRSST':
        fname = '/media/danilo/Danilo/mestrado/WesternAtlantic_MHWs/data/OISST/OISSTv2.1981.2018.nc'
        ncin  = xr.open_dataset(fname)
        # select coordinate for Cananeia
        ilon = 9
        ilat = 17
        sst  = ncin.sst[:,ilon,ilat].values
        time = pd.date_range(start='1981-09-01',end='2018-09-15',freq='1D')

        df_sst = pd.DataFrame({'sst':sst},index=time)

    return df_sst

def annual_average(mhws,clim):

    # Annual averages
    years = mhwBlock['years_centre']
    plt.figure(figsize=(13,7))
    plt.subplot(2,2,2)
    plt.plot(years, mhwBlock['count'], 'k-')
    plt.plot(years, mhwBlock['count'], 'ko')
    if np.abs(trend['count']) - dtrend['count'] > 0:
         plt.title('Frequency (trend = ' + '{:.2}'.format(10*trend['count']) + '* per decade)')
    else:
         plt.title('Frequency (trend = ' + '{:.2}'.format(10*trend['count']) + ' per decade)')
    plt.ylabel('[count per year]')
    plt.grid()
    plt.subplot(2,2,1)
    plt.plot(years, mhwBlock['duration'], 'k-')
    plt.plot(years, mhwBlock['duration'], 'ko')
    if np.abs(trend['duration']) - dtrend['duration'] > 0:
        plt.title('Duration (trend = ' + '{:.2}'.format(10*trend['duration']) + '* per decade)')
    else:
        plt.title('Duration (trend = ' + '{:.2}'.format(10*trend['duration']) + ' per decade)')
    plt.ylabel('[days]')
    plt.grid()
    plt.subplot(2,2,4)
    plt.plot(years, mhwBlock['intensity_max'], '-', color=col_evMax)
    plt.plot(years, mhwBlock['intensity_mean'], 'k-')
    plt.plot(years, mhwBlock['intensity_max'], 'o', color=col_evMax)
    plt.plot(years, mhwBlock['intensity_mean'], 'ko')
    plt.legend(['Max', 'mean'], loc=2)
    if (np.abs(trend['intensity_max']) - dtrend['intensity_max'] > 0) * (np.abs(trend['intensity_mean']) - dtrend['intensity_mean'] > 0):
        plt.title('Intensity (trend = ' + '{:.2}'.format(10*trend['intensity_max']) + '* (max), ' + '{:.2}'.format(10*trend['intensity_mean'])  + '* (mean) per decade)')
    elif (np.abs(trend['intensity_max']) - dtrend['intensity_max'] > 0):
        plt.title('Intensity (trend = ' + '{:.2}'.format(10*trend['intensity_max']) + '* (max), ' + '{:.2}'.format(10*trend['intensity_mean'])  + ' (mean) per decade)')
    elif (np.abs(trend['intensity_mean']) - dtrend['intensity_mean'] > 0):
        plt.title('Intensity (trend = ' + '{:.2}'.format(10*trend['intensity_max']) + ' (max), ' + '{:.2}'.format(10*trend['intensity_mean'])  + '* (mean) per decade)')
    else:
        plt.title('Intensity (trend = ' + '{:.2}'.format(10*trend['intensity_max']) + ' (max), ' + '{:.2}'.format(10*trend['intensity_mean'])  + ' (mean) per decade)')
    plt.ylabel(r'[$^\circ$C]')
    plt.grid()
    plt.subplot(2,2,3)
    plt.plot(years, mhwBlock['intensity_cumulative'], 'k-')
    plt.plot(years, mhwBlock['intensity_cumulative'], 'ko')
    if np.abs(trend['intensity_cumulative']) - dtrend['intensity_cumulative'] > 0:
        plt.title('Cumulative intensity (trend = ' + '{:.2}'.format(10*trend['intensity_cumulative']) + '* per decade)')
    else:
        plt.title('Cumulative intensity (trend = ' + '{:.2}'.format(10*trend['intensity_cumulative']) + ' per decade)')
    plt.ylabel(r'[$^\circ$C$\times$days]')
    plt.grid()

def duration_plot(mhws,clim):
    evs = np.argsort(mhws['duration'])[-10:]
    plt.clf()
    for i in range(10):
        ev = evs[-(i+1)]
        plt.subplot(5,2,i+1)
        # Find indices for all ten MHWs before and after event of interest and shade accordingly
        for ev0 in np.arange(max(ev-10,0), min(ev+11,mhws['n_events']-1), 1):
            t1 = np.where(t==mhws['time_start'][ev0])[0][0]
            t2 = np.where(t==mhws['time_end'][ev0])[0][0]
            plt.fill_between(dates[t1:t2+1], sst[t1:t2+1], clim['thresh'][t1:t2+1], color=col_ev)
        # Find indices for MHW of interest (2011 WA event) and shade accordingly
        t1 = np.where(t==mhws['time_start'][ev])[0][0]
        t2 = np.where(t==mhws['time_end'][ev])[0][0]
        plt.fill_between(dates[t1:t2+1], sst[t1:t2+1], clim['thresh'][t1:t2+1], color=col_evMax)
        # Plot SST, seasonal cycle, threshold, shade MHWs with main event in red
        plt.plot(dates, sst, 'k-', linewidth=2)
        plt.plot(dates, clim['thresh'], col_thresh, linewidth=2)
        plt.plot(dates, clim['seas'], col_clim, linewidth=2)
        plt.title('Number ' + str(i+1))
        plt.xlim(mhws['time_start'][ev]-150, mhws['time_end'][ev]+150)
        if coldSpells:
            plt.ylim(clim['seas'].min() + mhws['intensity_max'][ev] - 0.5, clim['seas'].max() + 1)
        else:
            plt.ylim(clim['seas'].min() - 1, clim['seas'].max() + mhws['intensity_max'][ev] + 0.5)
    plt.ylabel(r'SST [$^\circ$C]')

def cumulative_intensity(mhws,clim):
    evs = np.argsort(np.abs(mhws['intensity_cumulative']))[-10:]
    plt.clf()
    for i in range(10):
        ev = evs[-(i+1)]
        plt.subplot(5,2,i+1)
        # Find indices for all ten MHWs before and after event of interest and shade accordingly
        for ev0 in np.arange(max(ev-10,0), min(ev+11,mhws['n_events']-1), 1):
            t1 = np.where(t==mhws['time_start'][ev0])[0][0]
            t2 = np.where(t==mhws['time_end'][ev0])[0][0]
            plt.fill_between(dates[t1:t2+1], sst[t1:t2+1], clim['thresh'][t1:t2+1], color=col_ev)
        # Find indices for MHW of interest (2011 WA event) and shade accordingly
        t1 = np.where(t==mhws['time_start'][ev])[0][0]
        t2 = np.where(t==mhws['time_end'][ev])[0][0]
        plt.fill_between(dates[t1:t2+1], sst[t1:t2+1], clim['thresh'][t1:t2+1], color=col_evMax)
        # Plot SST, seasonal cycle, threshold, shade MHWs with main event in red
        plt.plot(dates, sst, 'k-', linewidth=2)
        plt.plot(dates, clim['thresh'], col_thresh, linewidth=2)
        plt.plot(dates, clim['seas'], col_clim, linewidth=2)
        plt.title('Number ' + str(i+1))
        plt.xlim(mhws['time_start'][ev]-150, mhws['time_end'][ev]+150)
        if coldSpells:
            plt.ylim(clim['seas'].min() + mhws['intensity_max'][ev] - 0.5, clim['seas'].max() + 1)
        else:
            plt.ylim(clim['seas'].min() - 1, clim['seas'].max() + mhws['intensity_max'][ev] + 0.5)
    plt.ylabel(r'SST [$^\circ$C]')

def mean_intensity(mhws,clim):
    evs = np.argsort(np.abs(mhws['intensity_mean']))[-10:]
    plt.clf()
    for i in range(10):
        ev = evs[-(i+1)]
        plt.subplot(5,2,i+1)
        # Find indices for all ten MHWs before and after event of interest and shade accordingly
        for ev0 in np.arange(max(ev-10,0), min(ev+11,mhws['n_events']-1), 1):
            t1 = np.where(t==mhws['time_start'][ev0])[0][0]
            t2 = np.where(t==mhws['time_end'][ev0])[0][0]
            plt.fill_between(dates[t1:t2+1], sst[t1:t2+1], clim['thresh'][t1:t2+1], color=col_ev)
        # Find indices for MHW of interest (2011 WA event) and shade accordingly
        t1 = np.where(t==mhws['time_start'][ev])[0][0]
        t2 = np.where(t==mhws['time_end'][ev])[0][0]
        plt.fill_between(dates[t1:t2+1], sst[t1:t2+1], clim['thresh'][t1:t2+1], color=col_evMax)
        # Plot SST, seasonal cycle, threshold, shade MHWs with main event in red
        plt.plot(dates, sst, 'k-', linewidth=2)
        plt.plot(dates, clim['thresh'], col_thresh, linewidth=2)
        plt.plot(dates, clim['seas'], col_clim, linewidth=2)
        plt.title('Number ' + str(i+1))
        plt.xlim(mhws['time_start'][ev]-150, mhws['time_end'][ev]+150)
        if coldSpells:
            plt.ylim(clim['seas'].min() + mhws['intensity_max'][ev] - 0.5, clim['seas'].max() + 1)
        else:
            plt.ylim(clim['seas'].min() - 1, clim['seas'].max() + mhws['intensity_max'][ev] + 0.5)
    plt.ylabel(r'SST [$^\circ$C]')

def topEvents(mhws,clim):
    evs = np.argsort(np.abs(mhws['intensity_max']))[-10:]
    plt.figure(figsize=(23,16))
    for i in range(10):
        ev = evs[-(i+1)]
        plt.subplot(5,2,i+1)
        # Find indices for all ten MHWs before and after event of interest and shade accordingly
        for ev0 in np.arange(max(ev-10,0), min(ev+11,mhws['n_events']-1), 1):
            t1 = np.where(t==mhws['time_start'][ev0])[0][0]
            t2 = np.where(t==mhws['time_end'][ev0])[0][0]
            plt.fill_between(dates[t1:t2+1], sst[t1:t2+1], clim['thresh'][t1:t2+1], color=col_ev)
        # Find indices for MHW of interest (2011 WA event) and shade accordingly
        t1 = np.where(t==mhws['time_start'][ev])[0][0]
        t2 = np.where(t==mhws['time_end'][ev])[0][0]
        plt.fill_between(dates[t1:t2+1], sst[t1:t2+1], clim['thresh'][t1:t2+1], color=col_evMax)
        # Plot SST, seasonal cycle, threshold, shade MHWs with main event in red
        plt.plot(dates, sst, 'k-', linewidth=2)
        plt.plot(dates, clim['thresh'], col_thresh, linewidth=2)
        plt.plot(dates, clim['seas'], col_clim, linewidth=2)
        plt.title('Number ' + str(i+1))
        plt.xlim(mhws['time_start'][ev]-150, mhws['time_end'][ev]+150)
        if coldSpells:
            plt.ylim(clim['seas'].min() + mhws['intensity_max'][ev] - 0.5, clim['seas'].max() + 1)
        else:
            plt.ylim(clim['seas'].min() - 1, clim['seas'].max() + mhws['intensity_max'][ev] + 0.5)
    plt.ylabel(r'SST [$^\circ$C]')
##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
# df_sst_ncep = load_sst(dataset='NCEP',path='')
df_sst_sate = load_sst(dataset='GHRSST',path='')

df_sst = df_sst_sate

###############################################################################
###############################################################################
#                  REMOVING INTER ANNUAL AND SEASONAL CYCLE                   #
###############################################################################
###############################################################################
df_sst['deseason'],df_sst['season'],beta = deseason_harmonic(df_sst['sst'].values,2,365.25)

df_sst['trend'] = df_sst['deseason'].rolling(365).mean()

df_sst.plot()
plt.show()

###############################################################################
###############################################################################
#                  DETECTING MARINE HEATWAVES EVENTS                          #
###############################################################################
###############################################################################
sst   = df_sst['sst'].values
t     = np.arange(df_sst.index[0].toordinal(), df_sst.index[-1].toordinal()+1)
dates = [date.fromordinal(tt.astype(int)) for tt in t]

mhws, clim = mhw.detect(t, sst,climatologyPeriod=[1982,2017])

# checking some informations about dataset, returning the index of the
# most intense event
ev = describe(mhws,clim)

# visualizing all events found
# visualizing_MHW(mhw,clim,ev=ev)

# Searching for the 2014 (DJFM)
df = pd.DataFrame(mhws) # converting to pandas.DataFrame
aux = []

for i,value in enumerate(mhws['date_start']):
    if (value.year == 2014) and ((value.month == 1) or (value.month == 2)):
        aux.append(i)

for i in aux:
    print('####')
    print('Durantion: '+str(mhws['duration'][i]))
    print('Start: '+str(mhws['date_start'][i].isoformat()))
    print('End:   '+str(mhws['date_end'][i].isoformat()))

# visualizing 2014 event
visualizing_MHW(mhw,clim,ev=51)

"""
# plotar evento de 2014 [51]

Usar groupby para classificar os eventos por ano

"""
