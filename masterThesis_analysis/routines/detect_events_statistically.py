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
from datetime import date
import scipy.ndimage as ndimage

import decomp

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def testDirection(x):
    if (x >= 180.) and (x <= 270.):
        return True
    else:
        return False

def detect_atmosphericBlocking(df,convert_uv2intdir=False,minDuration=10,maxGaps=2,mergeEvents=False):
    """Detect presence of winds from the 1st quarter, during 10+ days.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with data. Index must be in datetime format.
    convert_uv2intdir : boolean
        Boolean for convert u,v components into intensity and direction.
    minDuration : int
        Minimum duration of each event.
    maxGap : int
        Maximum gap between sequencial events, in days.
    mergeEvents : boolean
        If True, will merge all events with gaps < maxGaps.

    Returns
    -------
    detected_events : dictionary
        Dictionary with all information from each event found.

    Credits
    -------
    Created by Danilo A. Silva (nilodna@gmail.com), based on the code
    of marineHeatWaves package, by Eric J. Oliver.

    """

    # dictionary to store all information found
    detected_events ={}
    detected_events['time_start'] = []
    detected_events['time_end']   = []
    detected_events['date_start'] = []
    detected_events['date_end']   = []
    detected_events['index_start']= []
    detected_events['index_end']  = []
    detected_events['duration']   = []

    # some configuration and constants
    # minDuration = 10                # minimum duration of each event
    # maxGap      = 2                 # maximum Gap between two sequencial events
    # mergeEvents = False             # link events with gap shorter than maxGap

    # convert to intensity and direction
    if convert_uv2intdir:
        intensity,direction = oceano.uv2intdir(df.wu.values,df.wv.values)
        d = pd.DataFrame({'intensity':intensity,'direction':direction},index=df.index)

    # convert frequency: from 6-hourly to daily
    d = df.resample('1D').mean()

    # convert pd.DataFrame.index into toordinal() sequence
    t = np.arange(d.index[0].toordinal(),d.index[-1].toordinal()+1)

    # boolean time series of 'True' and 'False'
    exceed_bool = d['direction'].apply(testDirection)
    # Find contiguous regions of exceed_bool = True
    events, n_events = ndimage.label(exceed_bool)

    # Find all events of duration >= minDuration (10 days)
    for ev in range(1,n_events+1):
        event_duration = (events == ev).sum()
        if event_duration < minDuration:
            continue
        detected_events['time_start'].append(t[np.where(events == ev)[0][0]])
        detected_events['time_end'].append(t[np.where(events == ev)[0][-1]])

    # linking events that occur before and a short gap (maxGap)
    if mergeEvents:
        gaps = np.array(detected_events['time_start'][1:]) - np.array(detected_events['time_end'][0:-1]) - 1
        if len(gaps) > 0:
            ev = np.where(gaps <= maxGaps)[0][0]
            detected_events['time_end'][ev] = detected_events[time_end][ev+1]
            del detected_events['time_start'][ev+1]
            del detected_events['time_end'][ev+1]

            gaps = np.array(detected_events['time_start'][1:]) - np.array(detected_events['time_end'][0:-1]) -1
            # if len(gaps) == 0:
            #     break

    detected_events['n_events'] = len(detected_events['time_start'])

    # creating informations about all events found, such as duration, dates, etc
    for ev in range(detected_events['n_events']):
        detected_events['date_start'].append(date.fromordinal(detected_events['time_start'][ev]))
        detected_events['date_end'].append(date.fromordinal(detected_events['time_end'][ev]))

        t_start = np.where(t == detected_events['time_start'][ev])[0][0]
        t_end   = np.where(t == detected_events['time_end'][ev])[0][0]
        duration = t_end - t_start

        detected_events['index_start'].append(t_start)
        detected_events['index_end'].append(t_end)
        detected_events['duration'].append(duration)

    return detected_events


def plot_bars(df):
    fig,ax = plt.subplots()
    evMax = np.argmax(evs.duration)

    ax.bar(range(df.n_events[0]),df.duration,width=.6,color=(.7,.7,.7))
    ax.bar(evMax,evs.duration[evMax],width=.6,color=(1.,.5,.5))

    ax.set_xlabel('Eventos')
    ax.set_ylabel('Duracao')

    return fig,ax

def plot_data(df,data):

    from datetime import timedelta,datetime

    # plotting the most durable event
    # evMax = np.argmax(df.duration)

    startOriginal = df.date_start - timedelta(days=10)
    finalOriginal = df.date_end + timedelta(days=10)

    fig,ax = plt.subplots(nrows=2)
    data[startOriginal:finalOriginal].wu.plot(ax=ax[0],label='Alongshore')
    data[df.date_start:df.date_end].wu.plot(ax=ax[0],color='k',label='Detected Event')
    ax[0].legend()

    data[startOriginal:finalOriginal].wv.plot(ax=ax[1],label='Cross Shore')
    data[df.date_start:df.date_end].wv.plot(ax=ax[1],color='k',label='Detected Event')
    ax[1].legend()

    plt.suptitle('Event duration: %s'%(df.duration),fontsize=24)

    return fig,ax

def plot_frequencia(df):
    y = pd.Series([d.year for d in df.date_start])
    howMany = y.value_counts()
    d = pd.DataFrame({'freq':howMany},index=howMany.index)

    # sort by increasing year
    d.sort_index(inplace=True)

    d.plot.bar()


def analise(df):
    # transform column date_start in DatetimeIndex
    df.index = pd.DatetimeIndex(df.date_start)

    freq_acumulada = df.duration.resample('Y').sum()
    freq_media     = df.duration.resample('Y').mean()

    # create a new dataframe only to turn this easier
    freqs = pd.DataFrame({'Acumulada':freq_acumulada,'Media':freq_media})
    freqs.index = freqs.index.year

    indAcumulada = np.argmax(freqs.Acumulada)
    indMedia     = np.argmax(freqs.Media)

    fig,ax = plt.subplots(ncols=2)

    ax[0].bar(freqs.index,freqs.Acumulada,width=.6,color=(.7,.7,.7))
    ax[0].bar(indAcumulada,freqs.max().Acumulada,width=.6,color=(1.,.5,.5))
    ax[0].set_title(u'Frequência Acumulada [dias]',fontsize=20)
    ax[0].set_ylabel(u'Dias acumulados')
    ax[0].set_xlabel(u'Anos')

    ax[1].bar(freqs.index,freqs.Media,width=.6,color=(.7,.7,.7))
    ax[1].bar(indMedia,freqs.max().Media,width=.6,color=(1.,.5,.5))
    ax[1].set_title(u'Frequência Média [dias]',fontsize=20)
    ax[1].set_ylabel(u'Dias em média')
    ax[1].set_xlabel(u'Anos')

    return freqs

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
os.system('clear')
BASE_DIR = oceano.make_dir()

NCEP_DIR = '/media/danilo/Danilo/mestrado/ventopcse/data/timeseries_cfsv2_pnboia.nc'
cfsv2 = xr.open_dataset(NCEP_DIR) # carregando netcdf
dct = {
    'wu': np.squeeze(cfsv2['U_GRD_L103'].values),
    'wv': np.squeeze(cfsv2['V_GRD_L103'].values),
}
cfsv2 = pd.DataFrame(dct,index=cfsv2.time.values)       # convertendo para pd.DataFrame

intensity,direction = decomp.uv2intdir(cfsv2.wu.values,cfsv2.wv.values,0,0)

df = pd.DataFrame({'intensity':intensity,'direction':direction},index=cfsv2.index)

events = detect_atmosphericBlocking(df,minDuration=10)

# converting evets into a dataframe
evs = pd.DataFrame(events)

# plotting the most durable events
evMax = np.argmax(evs.duration)
plot_data(evs.iloc[evMax,:],cfsv2)

# plotting bars with all event's duration
# plot_bars(evs)

# plotting the 5 most durable events
selected = evs.nlargest(5,columns=['duration'])
selected.index = range(5)

fig,ax = plt.subplots(nrows=5)

for i,row in selected.iterrows():
    oceano.stickplot(cfsv2[row[1]:row[0]],ax[i])
    ax[i].legend([str(row[2])])

# checando a direção dos dados
import windrose

windrose.plot_windrose_df(df,var_name='intensity')
plt.title('Para onde o vento vai ...',fontsize=24)





############ teste zone


##############################################################################
#                                  CFSR                                      #
##############################################################################

NCEP_DIR = '/media/danilo/Danilo/mestrado/ventopcse/data/timeseries_cfsr_lajedesantos_1979_2011.nc'
cfsv2 = xr.open_dataset(NCEP_DIR) # carregando netcdf
dct = {
    'wu': np.squeeze(cfsv2['U_GRD_L103'].values),
    'wv': np.squeeze(cfsv2['V_GRD_L103'].values),
}
cfsv2 = pd.DataFrame(dct,index=cfsv2.time.values)       # convertendo para pd.DataFrame

intensity,direction = decomp.uv2intdir(cfsv2.wu.values,cfsv2.wv.values,0,0)

df = pd.DataFrame({'intensity':intensity,'direction':direction},index=cfsv2.index)

events = detect_atmosphericBlocking(df,minDuration=10)

# converting evets into a dataframe
evs = pd.DataFrame(events)
