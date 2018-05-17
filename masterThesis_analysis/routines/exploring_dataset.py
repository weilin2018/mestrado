"""

Function to explore timeseries at the South Brazil Bight location,
looking for episodes of anomalous cold front, persisting more than 5 days.

In this program we use two dataset:
    MODERN-ERA RETROSPECTIVE ANALYSIS FOR RESEARCH AND APPLICATIONS (MERRA)
    CLIMATE FORECAST SYSTEM VERSION 2 (CFSv2)


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

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here

# read files with wind's field and extract a single point timeserie
def extract_timeseries_from_MERRA(nfiles,loc=[6,11]):
    """Read files with wind's field and extract a single point timeserie.

    Data downloaded from Simple Subset Wizard (SSW) from EarthData,
    dataset MERRA-2 tavg1 2D, with 1-hourly frequency.

    Link to Download
    ----------------
    https://disc.gsfc.nasa.gov/SSW/#keywords=tavg1_2d_ocn_N

    Credits
    -------
    Created by Danilo A. Silva <nilodna@gmail.com>

    Parameters
    ----------
    nfiles : numpy.ndarray
        Array with all files that must be read by this functions. Use glob and
        glob.sort() to create this array and sort by name.
    loc : list
        List with indexes for latitude and longitude to extract data.

    Returns
    -------
    wu,wv,time : numpy.ndarray
        Eastward and Northward 10m height wind and time.

    Example
    -------
    >>> nfiles = glob.glob('path/to/files/*.nc')
    >>> nfiles.sort()
    >>> wu,wv,time = extractData_from_MERRA(nfiles=nfiles, loc=[-46,-23])

    """
    nfiles.sort()               # sort files in case the user don't do that

    ilon = loc[0]               # extract index for longitude
    ilat = loc[1]               # and latitude

    wu,wv,time = [],[],[]       # create list to store the data extracted

    # loop to read each file in nfiles
    for f in nfiles:
        ncdata = xr.open_dataset(f)             # import netcdf file

        u = ncdata['U10M'].values[:,ilat,ilon]  # extract eastward,
        v = ncdata['V10M'].values[:,ilat,ilon]  # northward 10m wind
        t = ncdata['time'].values               # and time

        for i in np.arange(0,len(t)):           # loop to read each hour
            wu.append(u[i])
            wv.append(v[i])
            time.append(t[i])

    # convert lists to np.ndarray
    wu = np.asarray(wu)
    wv = np.asarray(wv)
    time = np.asarray(time)

    return wu,wv,time

# read files with timeseries in a single point from CFSv2
def read_cfsv2(nfiles):
    """Read and extract data from CFSv2 files.

    Parameters
    ----------
    nfiles : list
        Sorted by name list.

    Returns
    -------
    wu, wv, time: arrays


    """
    # define some variable
    wu   = []
    wv   = []
    time = []

    # read files
    for f in nfiles:
        ncdata = xr.open_dataset(f)
        u     = ncdata['U_GRD_L103'].values
        v     = ncdata['V_GRD_L103'].values
        t     = ncdata['time'].values

        for k in [0,1,2,3]:
            wu.append(np.squeeze(u[k,:,:]))
            wv.append(np.squeeze(v[k,:,:]))
            time.append(np.squeeze(t[k]))

    wu   = np.asarray(np.squeeze(wu))
    wv   = np.asarray(np.squeeze(wv))
    time = np.asarray(time)

    return wu,wv,time

# plot a stickplot using quiver
def stickplot(df,ax):
    """Create a stickplot.

    With the u and v components given as argument in pd.DataFrame df,
    this function plot a stickplot, using MPL quiver.

    Parameters
    ----------
    df : pandas.Dataframe
        dataframe containing wu and wv components and datetimeindex.
    ax : matplotlib.axis
        axis to stickplot

    Example
    -------
    >>>
    >>>

    Credits
    -------
    Stephane Raynaud
    http://permalink.gmane.org/gmane.comp.python.matplotlib.general/24155
    """

    # creating the date axis
    # dates = pd.to_datetime(df.index)
    # extracting components from dataframe
    u = df['cross'].values
    v = df['along'].values

    # calculating speed
    spd = np.sqrt(u**2 + v**2)
    maxSpd = np.nanmax(spd)

    # plotting
    # fig, ax = plt.subplots()

    qiv = ax.quiver(df.index, [[0]*len(df)], u, v, headlength=0, headwidth=0, headaxislength=0 )
    key = ax.quiverkey(qiv, 0.25, 0.75, maxSpd, "%0.2f $m^{2}s^{-1}$"%(maxSpd), labelpos='N', coordinates='axes' )

    # plot a horizontal line in y=0.0
    ax.axhline(y=0.0,xmin=df.index[0],xmax=df.index[-1],linewidth=1.,color='black')

    plt.setp(ax.get_yticklabels(), visible=False)
    ax.xaxis_date()

    # ax.set_xticks(['2012-01-07', '2012-01-21', '2012-02-04', '2012-02-18', '2012-03-03', '2012-03-17', '2012-03-31'])

    return ax

# plot two databases
def statisticalAnaysis(data1,data2,ind):
    ''' '''

    # interpoate to disapear with missing values
    data1.interpolate(inplace=True)
    data2.interpolate(inplace=True)

    # computing some statistical parameters, such as skill and correlation
    skillCross = 0.#oceano.skill_willmott(data1.cross.values,data2.cross.values)
    skillAlong = 0.#oceano.skill_willmott(data1.along.values,data2.along.values)

    # corrCross =
    # corrAlong =

    # plotting data and statistical analysis
    fig, ax = plt.subplots(nrows=2,ncols=1,sharex=True)

    ax[0].plot(data1.cross,label='MERRA')
    ax[0].plot(data2.cross,label='CFSv2')
    ax[0].margins(0)
    ax[0].set_ylim(-15,15)
    ax[0].set_title('Cross shore (skill = %0.2f)'%(skillCross))

    ax[0].legend(loc='best')

    ax[1].plot(data1.along,label='MERRA')
    ax[1].plot(data2.along,label='CFSv2')
    ax[1].margins(0)
    ax[1].set_ylim(-15,15)
    ax[1].set_title('Along shore (skill = %0.2f)'%(skillAlong))

    ax[1].legend(loc='best')

    plt.suptitle('Laje de Santos Location [%s] \n MERRA v CFSv2'%(ind), fontsize=26)
    plt.show()

# find all indexes where along along shore component is positive (northwards)
def findNorthwardsWind(data,condition='v'):
    """Find all indiexes where along shore is positive (northwards wind)
    and cross shore is negative (westward wind).

    Parameters
    ----------
    data : pandas.DataFrame
        Dataframe with along and cross shore data, with a datetimeindex.
    condition : string
        Which condition you want to use.

    Returns
    -------
    indexes : numpy.ndarray
        Vector with all indexes of all data that passed the criteria.

    """
    if condition == 'uv':
        # select only those data where:
        # along shore is positive AND cross shore is negative
        indexes = np.where((data.along.values > 0.) & (data.cross.values < 0.))
    elif condition == 'v':
        # select only those data where:
        # along shore is positive AND cross shore is negative
        indexes = np.where(data.along.values > 0.)
    elif condition == 'u':
        # select only those data where:
        # along shore is positive AND cross shore is negative
        indexes = np.where(data.cross.values < 0.)

    return indexes

def split_list(n):
    """will return the list index"""
    return [(x+1) for x,y in zip(n, n[1:]) if y-x != 1]

def get_sub_list(my_list):
    """will split the list base on the index

    credits:
    https://stackoverflow.com/questions/2361945/detecting-consecutive-integers-in-a-list

    """
    my_index = split_list(my_list)
    output = list()
    prev = 0
    for index in my_index:
        new_list = [ x for x in my_list[prev:] if x < index]
        output.append(new_list)
        prev += len(new_list)
    output.append([ x for x in my_list[prev:]])

    return output

# select all sequential data
def findSequence(indexes,above=12):
    """Select all sequential data with some criteria.

    In this function we assume the follow criteria:

        All data where along shore componenet is positive for more than
        3 days (12 datas considering a 6-hourly frequency) or 5 days (20 datas).

    Parameters
    ----------
    indexes : list
        List with all indexes of data.
    above : integer
        Number of data that must be sequential. Default is 12 datas, representing
        3 dyas in a 6-hourly frequency of data (from CFSv2, for instance).

    Returns
    -------
    sub : list
        List with all index of sequential data.

    """

    # ensure this is a list
    indexes = list(indexes)

    output = get_sub_list(indexes)

    sub = []

    for i in output:
        if len(i) > above:
            sub.append(i)

    return sub

# to vertical shade graph
def shadeSpaces(ax,index,sub):

    for sequential in sub:
        ax.axvspan(index[sequential][0], index[sequential][-1], facecolor='k',
            alpha=.5)

# stickplot
def stickplot(df,ax):
    """Create a stickplot.

    With the u and v components given as argument in pd.DataFrame df,
    this function plot a stickplot, using MPL quiver.

    Parameters
    ----------
    df : pandas.Dataframe
        dataframe containing wu and wv components and datetimeindex.
    ax : matplotlib.axis
        axis to stickplot

    Example
    -------
    >>>
    >>>

    Credits
    -------
    Stephane Raynaud
    http://permalink.gmane.org/gmane.comp.python.matplotlib.general/24155
    """

    # creating the date axis
    # dates = pd.to_datetime(df.index)
    # extracting components from dataframe
    u = df['cross'].values
    v = df['along'].values

    # calculating speed
    spd = np.sqrt(u**2 + v**2)
    maxSpd = np.nanmax(spd)

    # plotting
    # fig, ax = plt.subplots()

    qiv = ax.quiver(df.index, [[0]*len(df)], u, v, headlength=0, headwidth=0, headaxislength=0 )
    key = ax.quiverkey(qiv, 0.25, 0.75, maxSpd, "%0.2f $m^{2}s^{-1}$"%(maxSpd), labelpos='N', coordinates='axes' )

    # plot a horizontal line in y=0.0
    ax.axhline(y=0.0,xmin=df.index[0],xmax=df.index[-1],linewidth=1.,color='black')

    plt.setp(ax.get_yticklabels(), visible=False)
    ax.xaxis_date()

    # ax.set_xticks(['2012-01-07', '2012-01-21', '2012-02-04', '2012-02-18', '2012-03-03', '2012-03-17', '2012-03-31'])

    return ax

# plot one database only
def plotTimeSeries(rotated,data,period,figdir=None):

    # interpolate
    data.interpolate(inplace=True)

    # plotting data
    fig, ax = plt.subplots(nrows=3,ncols=1,sharex=True,figsize=(15,10))

    ax[0].plot(rotated.cross,label='CFSv2')
    ax[0].margins(0)
    ax[0].set_ylim(-15,15)
    ax[0].set_title('Cross Shore')

    ax[0].legend(loc='lower left')

    ax[1].plot(rotated.along,label='CFSv2')
    ax[1].margins(0)
    ax[1].set_ylim(-15,15)
    ax[1].set_title('Along Shore')

    ax[1].legend(loc='lower left')

    ax[2] = oceano.stickplot(data,ax[2])
    ax[2].set_ylim(-0.7,1.)
    ax[2].set_title('CFSv2')

    plt.suptitle('Laje de Santos Location\n CFSv2 dataset %s'%(period.replace('-','/')),fontsize=24)

    # plotar 00H
    zeroHour = rotated.iloc[rotated.index.hour == 0]
    ax[0].scatter(zeroHour.index,zeroHour['cross'].values,s=30,c='k',marker='*')
    ax[1].scatter(zeroHour.index,zeroHour['along'].values,s=30,c='k',marker='*')

    # plot horizontal line at y=0
    ax[0].axhline(c='k')
    ax[1].axhline(c='k')

    if figdir:
        outFile = figdir + str(period) + '.png'
        plt.savefig(outFile)
        # plt.show()
    else:
        plt.show()

    plt.close('all')


    # ax[2] = stickplot(data,ax[2])



##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

os.system('clear')
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/','ventopcse/data/MERRA/data/cut/')
SAVE_DIR = BASE_DIR + 'masterThesis_analysis/routines/pickles/'
FIGU_DIR = BASE_DIR + 'masterThesis_analysis/figures/exploring_datasets/'
NCEP_DIR = BASE_DIR.replace('github/', 'ventopcse/data/CFSv2/data/')

# ------------------------- reading merra --------------------------------- #

# try:
#     # load pickle file
#     merra = oceano.load_data(SAVE_DIR+'merraData.pickle')
#
# except:
#     # select files, creating a list with all files to be read
#     nfiles = glob.glob(DATA_DIR+'*.nc')
#     nfiles.sort()
#
#     # load data and save pickle
#     wu,wv,time = extract_timeseries_from_MERRA(nfiles,loc=[6,11])
#
#     # put data into pandas.DataFrame to better visualization
#     merra = pd.DataFrame({'wu':wu,'wv':wv},index=pd.DatetimeIndex(time))
#
#     oceano.save_data(merra,SAVE_DIR+'merraData.pickle')

# ---------------------------- reading cfsv2 ------------------------------- #

try:
    # load pickle file
    cfsv2 = oceano.load_data(SAVE_DIR+'cfsv2Data.pickle')
except:
    # select files, creating a list with all files to be read
    nfiles = glob.glob(NCEP_DIR+'cdas1*')
    nfiles.sort()

    ### extract timeseries
    wu,wv,time = read_cfsv2(nfiles)

    # put data into pandas.DataFrame to better visualization
    cfsv2 = pd.DataFrame({'wu':wu,'wv':wv},index=pd.DatetimeIndex(time))

    oceano.save_data(cfsv2,SAVE_DIR+'cfsv2Data.pickle')

# ---------------------------- rotate components --------------------------- #
# rotate u and v to cross and along shore components, using a phi of 50 degrees
# as used by Dottori and Belmiro (2009)

cross,along = oceano.rotateVectors(cfsv2.wu.values,cfsv2.wv.values,40.)

rotated = pd.DataFrame({'along':along,'cross':cross},index=cfsv2.index)

# -------------------------- smoothing timeseries -------------------------- #
# apply a digital filter, with a window of 30 hours, to remove high
# frequency signals.

filtered_data  = rotated.resample('12H').mean()
filtered_cfsv2 = cfsv2.resample('12H').mean()
# ---------------------------- data visualization -------------------------- #

####### PLOTTING TIMESERIES

# group by month AND year, so we can plot every month of every year separately
years  = np.arange(2012,2018)
months = np.arange(1,13)

dictData = {}

monthPlot = {}

for year in years[3:4]:
    for month in months:
        period  = '%s-%s'%(str(year),str(month).zfill(2))

        data = filtered_data[period]                   # slicing dataframe
        cfsv = filtered_cfsv2[period]

        # looking for the anomalous episodes
        indexes = findNorthwardsWind(data)[0]
        sub     = findSequence(indexes,above=12)

        # if you want to save figure, uncomment next line and comment the next one
        # plotTimeSeries(data,cfsv,period,figdir=FIGU_DIR)
        plotTimeSeries(data,cfsv,period)

        # plotTimeSeries(data,period,sub)
        if month == 1:
            # monthPlot[period+'_cross'] = data.cross
            monthPlot[period+'_along'] = data.along

        dictData[period] = len(sub)

####### CREATING BAR PLOT
index = pd.DatetimeIndex(dictData.keys())
freq = pd.DataFrame({'episodes': dictData.values()}, index=index)

fig, ax = plt.subplots()

for key in freq.index.sort_values():
    ax.bar(key,freq[key].values,label=key.strftime('%Y-%m'))

plt.bar(freq.index,freq.episodes.values,label=freq.index.strftime('%Y-%m'))
plt.xticks(freq.index, freq.index.strftime('%Y-%m'),rotation=60)
plt.show()





###############################
# plotar stickplot dos periodos ja selecionados para confirmar direcao


####### PLOTTING STICKPLOT
# load pickle file
cfsv2 = oceano.load_data(SAVE_DIR+'cfsv2Data.pickle')

# rotate u and v
cross,along = oceano.rotateVectors(cfsv2.wu.values,cfsv2.wv.values,40.)
rotated = pd.DataFrame({'along':along,'cross':cross},index=cfsv2.index)

periods = '20120327_20120404 20121001_20121017 20121117_20121120 20121127_20121205 20130117_20130126 20130304_20130430 20130821_20130829 20130916_20131008 20131224_20140101 20140101_20140406 20140805_20140821 20141013_20141027 20141105_20141115 20141127_20141202 20141206_20150107 20150123_20150129 20150204_20150221 20150307_20150312 20140307_20140312 20150316_20150321 20160101_20160110 20160114_20160125 20160208_20160217 20160417_20160105'.split(" ")

for per in periods:
    p = per.split("_")
    start = p[0][:4]+'-'+p[0][4:6]+'-'+p[0][6:]
    end   = p[1][:4]+'-'+p[1][4:6]+'-'+p[1][6:]
    cut = cfsv2[start:end]
    cut2 = rotated[start:end]

    if len(cut.index) != 0:
        # plot data and skill
        fig, ax = plt.subplots(nrows=3,ncols=1,sharex=True,figsize=(15,10))

        ax[0].plot(cut2.cross,label='CFSv2')
        ax[0].margins(0)
        ax[0].set_ylim(-15,15)
        ax[0].set_title('Cross Shore')

        ax[0].legend(loc='lower left')

        ax[1].plot(cut2.along,label='CFSv2')
        ax[1].margins(0)
        ax[1].set_ylim(-15,15)
        ax[1].set_title('Along Shore')

        ax[1].legend(loc='lower left')

        ax[2] = oceano.stickplot(cut,ax[2])
        ax[2].set_ylim(-0.7,1.)
        ax[2].set_title('CFSv2')

        plt.suptitle('CFSv2\n%s' % (per),fontsize=26)

        outFile = FIGU_DIR + str(per) + '.png'
        outFile = outFile.replace('[','')
        outFile = outFile.replace(']','')
        plt.savefig(outFile)
        plt.close("all")
