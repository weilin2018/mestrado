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
def read_reanalysis(nfiles):
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

        for k in np.arange(0,t.shape[0]):
            wu.append(np.squeeze(u[k,:,:]))
            wv.append(np.squeeze(v[k,:,:]))
            time.append(np.squeeze(t[k]))

    wu   = np.asarray(np.squeeze(wu))
    wv   = np.asarray(np.squeeze(wv))
    time = np.asarray(time)

    return wu,wv,time

def read_cfsr( nfiles):

    wu,wv,time = [],[],[]

    for f in nfiles:
        ncdata = xr.open_dataset(f)
        u = ncdata['U_GRD_L103'].values
        v = ncdata['V_GRD_L103'].values
        t = ncdata['time'].values

        for k in np.arange(0,t.shape[0]):
            wu.append(np.squeeze(u[k,:,:]))
            wv.append(np.squeeze(v[k,:,:]))
            time.append(np.squeeze(t[k]))

    wu   = np.asarray(np.squeeze(wu))
    wv   = np.asarray(np.squeeze(wv))
    time = np.asarray(time)

    return wu,wv,time



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
#
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
CFSR_DIR = BASE_DIR.replace('github/', 'ventopcse/data/CFSR/Laje/')

################################
#           PARTE I            #
################################

# loading files from CFSR and plotting stickplots
try:
    # load pickle file
    cfsr = oceano.load_data(SAVE_DIR+'cfsrData.pickle')
except:
    # select files, creating a list with all files to be read
    nfiles = glob.glob(CFSR_DIR+'*.nc')
    nfiles.sort()

    ### extract timeseries
    wu,wv,time = read_reanalysis(nfiles)

    # put data into pandas.DataFrame to better visualization
    cfsr = pd.DataFrame({'wu':wu,'wv':wv},index=pd.DatetimeIndex(time))

    oceano.save_data(cfsr,SAVE_DIR+'cfsrData.pickle')

# ---------------------------- rotate components --------------------------- #
# rotate u and v to cross and along shore components, using a phi of 40 degrees
# as used by Dottori and Castro (2009)

cross,along = oceano.rotateVectors(cfsr.wu.values,cfsr.wv.values,40.)

rotated_cfsr = pd.DataFrame({'along':along,'cross':cross},index=cfsr.index)

# -------------------------- smoothing timeseries -------------------------- #
# apply a digital filter, with a window of 30 hours, to remove high
# frequency signals.

# filtered_data  = rotated.resample('12H').mean()
filtered_cfsr = cfsr.resample('12H').mean()

# ---------------------------- data visualization -------------------------- #

################################
#           PARTE I            #
#   from 1981 to 2011 - CFSR   #
################################

# plotar stickplots de todos os anos (2012 a 2017) do periodo 01-01 a 02-28
years = list(set(cfsr.index.year))
years.sort()
years = years[2:]

def plot_years(years):
    i = 0               # contador de graficos

    fig,ax = plt.subplots(nrows=len(years),ncols=1,figsize=(15,10))

    for y in years[0:5]:
        start = str(y)+'-01-01'
        final = str(y)+'-02-28'
        cut = filtered_cfsr[start:final]

        ax[i] = oceano.stickplot(cut,ax[i])
        ax[i].set_ylim(-0.7,1.)
        ax[i].set_ylabel(y)

        if i != 5:
            ax[i].get_xaxis().set_ticks([])
        # else:
        #     ax[i].get_xaxis().set_ticks(cut.index.day)

        i += 1

    plt.suptitle(u'Período 01/01 a 28/02 - Laje de Santos', fontsize=24)
    OUTFIGURE = '0_JF_'+str(years[0])+str(years[-1])+'.png'
    plt.savefig(FIGU_DIR+OUTFIGURE)


# 80's
plot_years(years=years[:5])
plot_years(years=years[5:10])
plot_years(years=years[10:15])
plot_years(years=years[15:20])
plot_years(years=years[20:25])
plot_years(years=years[25:])



################################
#           PARTE II           #
#   from 2012 to 2017 - CFSv2  #
################################
try:
    # load pickle file
    cfsv2 = oceano.load_data(SAVE_DIR+'cfsv2Data.pickle')
except:
    # select files, creating a list with all files to be read
    nfiles = glob.glob(NCEP_DIR+'cdas1*')
    nfiles.sort()

    ### extract timeseries
    wu,wv,time = read_reanalysis(nfiles)

    # put data into pandas.DataFrame to better visualization
    cfsv2 = pd.DataFrame({'wu':wu,'wv':wv},index=pd.DatetimeIndex(time))

    oceano.save_data(cfsv2,SAVE_DIR+'cfsv2Data.pickle')

# ---------------------------- rotate components --------------------------- #
# rotate u and v to cross and along shore components, using a phi of 40 degrees
# as used by Dottori and Castro (2009)

cross,along = oceano.rotateVectors(cfsv2.wu.values,cfsv2.wv.values,40.)

rotated_cfsv2 = pd.DataFrame({'along':along,'cross':cross},index=cfsv2.index)

# -------------------------- smoothing timeseries -------------------------- #
# apply a digital filter, with a window of 30 hours, to remove high
# frequency signals.

# filtered_data  = rotated.resample('12H').mean()
filtered_cfsv2 = cfsv2.resample('12H').mean()

# ---------------------------- data visualization -------------------------- #

################################
#           PARTE III          #
#     longshore component      #
################################

# plotar stickplots de todos os anos (2012 a 2017) do periodo 01-01 a 02-28
years = list(set(cfsr.index.year))
years.sort()
years = years[2:]

def plot_years_long(years):
    i = 0               # contador de graficos

    fig,ax = plt.subplots(nrows=len(years),ncols=1,figsize=(15,10))

    for y in years[0:5]:
        start = str(y)+'-01-01'
        final = str(y)+'-02-28'
        cut = rotated_cfsr[start:final]

        ax[i].plot(cut.index, cut.along.values)
        ax[i].axhline(color='k',alpha=.4)
        ax[i].margins(0)
        ax[i].set_ylim(-15,15)
        ax[i].set_ylabel(y)

        # if i != 5:
        #     ax[i].get_xaxis().set_ticks([])
        # else:
        #     ax[i].get_xaxis().set_ticks(cut.index.day)

        i += 1

    title = 'Long shore velocity component at Laje de Santos \n JF from %i to %i' % (years[0], years[-1])
    plt.suptitle(title, fontsize=24)
    OUTFIGURE = 'long_JF_'+str(years[0])+str(years[-1])+'.png'
    # plt.savefig(FIGU_DIR+OUTFIGURE)


# 80's
plot_years_long(years=years[:5])
plot_years_long(years=years[5:10])
plot_years_long(years=years[10:15])
plot_years_long(years=years[15:20])
plot_years_long(years=years[20:25])
plot_years_long(years=years[25:])


################################
#           PARTE IV           #
################################
"""
    Procurando por períodos típicos de Janeiro e Fevereiro,
    com passagem de 4 a 5 frentes frias, que serão utilizados
    como nosso experimento controle.

    - procurar por stickplots de JF
    - destacar períodos com 2 ou 3 dias com ventos do quadrante sul 
    - contabilizar os eventos e, o período que tiver 4 ou 5 passagens, deve ser 
        classificado.
"""

