'''
Próximo passo:
    . identificado: problema está na leitura dos dados da laje.
        - pensar em como resolver isso!!!!

'''
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
import datetime

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################

# convert columns to datetime
def convert2datetime(year,month,day,hour):
    '''
        receive some year,months,days and hours
        convert to datetime
        return a valid datetime
    '''

    import datetime

    return datetime.datetime(*map(int, [year,month,day,hour]))

# read files with wind's field and extract a single point timeserie
def extractData_from_MERRA(nfiles,loc):
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

# read files from laje de santos
def readFiles_Laje(nfiles):
    """Read file from Laje dos Santos Meteorological Station.


    Parameters
    ----------
    nfiles : numpy.ndarray
        List with path of all files to be read.

    Returns
    -------
    date,intensity,direction : numpy,ndarray
        Description of returned object.

    """
    datas = []
    inten = []    # intensity
    direc = []    # directions

    for f in nfiles:
        # read file
        data = np.loadtxt(f)

        for i in np.arange(0,data.shape[0],1):
            dt = convert2datetime(data[i,0], data[i,1], data[i,2], data[i,3])

            datas.append(dt)
            inten.append(data[i,7])
            direc.append(data[i,8])

    return np.asarray(datas),np.asarray(inten),np.asarray(direc)

# function only to calculate wu and wv using direction and intensity
def windComponent(intensity, direction, dates, decliMag, deg=False):
    """function only to calculate wu and wv using direction and intensity.

    Parameters
    ----------
    intensity : numpy.ndarray
        Description of parameter `intensity`.
    direction : numpy.ndarray
        Description of parameter `direction`.
    dates : numpy.ndarray
        Description of parameter `dates`.
    decliMag : float
        Description of parameter `decliMag`.
    deg : boolean
        Description of parameter `deg`.

    Returns
    -------
    type
        Description of returned object.

    """
    # keeping arguments as arrays
    direction = np.asarray(direction)
    intensity = np.asarray(intensity)

    # wu,wv = oceano.spdir2uv(intensity,direction,deg=False)
    # wu,wv = oceano.compass2uv(direction,intensity,kind='meteo2uv')
    # usando rotina de conversão da Paulinha
    wu,wv = oceano.dirmag2uv(direction,intensity,decliMag,'norte')

    # create pandas.dataframe
    i = pd.DatetimeIndex(dates)
    df = pd.DataFrame({'wu':wu, 'wv':wv},
                        index=i)

    # quality control
    df[df['wu'] > 20] = np.nan
    df[df['wu'] <-20] = np.nan

    df[df['wv'] > 20] = np.nan
    df[df['wv'] <-20] = np.nan

    return df

# extract data from Laje dos Santos Meteorological Station
def lajeSantos(DATA_DIR):

    # select files
    nfiles = glob.glob(DATA_DIR + '*.txt')
    nfiles.sort()

    # extract information from files using an auxiliar function
    dates,intensity,directions = readFiles_Laje(nfiles)

    # corrects magnetic declination
    decliMag = -21.03 #degrees from real north
    laje = windComponent(intensity,directions,dates,decliMag)

    return laje

##############################################################################
#                          [GEN]PLOTS                                        #
##############################################################################
def load_MERRA_grid(nfile):

    # import netCDF file
    ncdata = xr.open_dataset(nfile)

    # extract latitude and longitude information
    glon = ncdata['lon'].values
    glat = ncdata['lat'].values

    return glon,glat

# plot MERRA grid and some given location
def plotGrid(point,nfile):
    '''
        Function to plot southeast brazilian coastline, as a background,
        and MERRA grid over the map. Also, this function plot some
        given location, to better visualize.

        args
            points      [tuple]: with array containing lat and lon

        return
            plt.show()

        example
            >>> points = ([-24.3194, -46.1803, 'Laje','s','b'], [-25.274, -44.93, 'PNBOIA','d','g'])
            >>> plotGrid(points)
    '''
    fig, ax = plt.subplots(figsize=(15,10))

    m = oceano.make_map(ax)

    # load grid
    glon,glat = load_MERRA_grid(nfile)
    glon,glat = np.meshgrid(glon,glat)

    # plot grid and other features
    m.plot(glon,glat,'k--',alpha=.4,latlon=True)
    m.plot(glon.T,glat.T,'k--',alpha=.4,latlon=True)

    for p in point:
        m.scatter(p[1], p[0], c=p[4], marker=p[3], s=50, label=p[2], latlon=True)

    plt.legend(loc='best')
    plt.title('MERRA grid \n Observed data location and nearest grid points')
    plt.show()

##############################################################################
#                                STATS                                       #
##############################################################################
def calculateMSE(x,y):
    """Calculate mean squared error.

     Calculate mean squared error using the formula:

            MSE = 1/n * sum((x-y)^2)

    Parameters
    ----------
    x : array
        real data as an array
    y : array
        modeled data as an array

    Returns
    -------
    mse : float
        the mean squared error calculated

    """

    return np.mean((x - y)**2)

def calculateCorr(x,y):
    """Calculate correlation between two series with same size.

    Using the Pearson's method, from scipy.stats.pearsonr, this function
    calculate the correlation between x and y, two array with same size and
    without missing values.

    Parameters
    ----------
    x : array
        real data.
    y : array
        modeled data.

    Returns
    -------
    float
        correlation coefficient.

    """
    import scipy.stats as stats

    return stats.pearsonr(x,y)

# smoothing the data using hamming window
def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
            w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y


##############################################################################
#                            MAIN CODE                                       #
##############################################################################

os.system('clear')

# define some constants
BASE_DIR = oceano.make_dir()
MERR_DIR = BASE_DIR.replace('github/', 'ventopcse/data/MERRA/')
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/data/Est_lajeSantos/2015/atualizado/')
ARGO_DIR = BASE_DIR.replace('github/', 'ventopcse/data/')

##############################################################################
#                          PLOTAR PONTOS                                     #
##############################################################################
points = (
    [-24.3194, -46.1803, 'Laje de Santos','s','b'],
    [-25.274, -44.93, 'PNBOIA/Santos',',','g'],
    # definir os pontos de grade mais próximo aos pontos de dados observados
    [-24.5, -46.25, '', 'o','k'],
    [-25.5 ,-45, '','o','k']
)

# para consulta posterior: os indices sao: Laje[6,11] e pnboia[8,9]

os.system('clear')
# select any nc file to send to the function
nfile = glob.glob(MERR_DIR+'summer2014_tavg1_2d_ocn/*.nc')[0]


print('Plotting map with locations')
plotGrid(points,nfile)

##############################################################################
#                          PLOTAR DADOS                                      #
##############################################################################
# importing MERRA data
nfiles = glob.glob(MERR_DIR+'summer2015_tavg1_2d_ocn/*.nc')
nfiles.sort()

wu,wv,time = extractData_from_MERRA(nfiles=nfiles,loc=[6,11])

# put data into pandas.DataFrame to better visualization
merra = pd.DataFrame({'wu':wu,'wv':wv},index=pd.DatetimeIndex(time))

# importing Laje dos Santos data
laje = lajeSantos(DATA_DIR=DATA_DIR)

# tentativa de uma filtragem das altas frequẽncias. Devo usar Hamming!!!!
# laje = laje.resample('30H').mean()
# merr = merra.resample('30H').mean()

# laje data only start in 2015-03-19, so we cut all data
lajeCut = laje['2015-04-21':'2015-12-01']
lajeCut['wu'] *= -1
merrCut = merr['2015-04-21':]

laje.interpolate(inplace=True)
merra.interpolate(inplace=True)

# create hamming window
w = np.hamming(30)

# filtering data using fftconvolve
wulaje = scipy.signal.fftconvolve(laje.wu.values,w,'same')
wvlaje = scipy.signal.fftconvolve(laje.wv.values,w,'same')
lajeFilt = pd.DataFrame({'wu':wulaje,'wv':wvlaje},index=laje.index)

wuMerr = scipy.signal.fftconvolve(merra.wu.values,w,'same')
wvMerr = scipy.signal.fftconvolve(merra.wv.values,w,'same')
merrFilt = pd.DataFrame({'wu':wuMerr,'wv':wvMerr},index=merra.index)


plotar(lajeFilt,merrFilt)

# funcao para facilitar a vida na visualizacao dos dados
def plotar(lajeCut,merrCut):
    ## comparing two datasets
    fig,ax = plt.subplots(nrows=2,ncols=1,sharex=True)

    ax[0].plot(lajeCut.index,lajeCut.wu.values,label='Laje dos Santos')
    ax[0].plot(merrCut.index,merrCut.wu.values,label='MERRA')
    ax[0].set_title('Eastward Component')

    ax[1].plot(lajeCut.index,lajeCut.wv.values,label='Laje dos Santos')
    ax[1].plot(merrCut.index,merrCut.wv.values,label='MERRA')
    ax[1].set_title('Northward Component')

    ax[0].legend(loc='best')
    ax[1].legend(loc='best')

    plt.suptitle('Laje dos Santos and MERRA nearest location data for DJFM/2014',fontsize=24)

    plt.show()
