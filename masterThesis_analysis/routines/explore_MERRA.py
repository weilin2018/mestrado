"""
    Explore the data from 2009-12 to 2017-12, looking for
    periods with a positive along shore velocity component persisting
    for more than 5 days, much more above the normal conditions of
    frontal systems (cold front).

    MODERN-ERA RETROSPECTIVE ANALYSIS FOR RESEARCH AND APPLICATIONS (MERRA)
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


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

os.system('clear')
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/','ventopcse/data/MERRA/data/')
SAVE_DIR = BASE_DIR + 'masterThesis_analysis/routines/pickles/'

# select files, creating a list with all files to be read
nfiles = glob.glob(DATA_DIR+'*.nc')
nfiles.sort()

### extract timeseries
os.system('Extracting data from a singlepoint')
wu,wv,time = extract_timeseries_from_MERRA(nfiles,loc=[6,11])

# put data into pandas.DataFrame to better visualization
os.system('Converting data into DataFrame')
merra = pd.DataFrame({'wu':wu,'wv':wv},index=pd.DatetimeIndex(time))




#
os.system('clear')
nfiles = glob.glob(DATA_DIR+'*.nc')
nfiles.sort()

for f in nfiles:
    os.system('Cutting %s'%(f))
    inFile = f
    newFile = f.replace('MERRA/data/','MERRA/data/cut/')
    os.system('cdo sellonlatbox,-60,-30,-30,-20 %s %s'%(inFile,newFile))
