"""
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

def load_data(fname,vars='elev'):
    """Function to load variables in vars from a netCDF file, given by fname.

    Parameters
    ----------
    fname : string
        String with full path to netCDF file to be read.
    vars : string
        Which variable to extract

    Returns
    -------
    lon,lat,variable : xarray.Dataset
    """

    # load file
	ncdata = xr.open_dataset(fname)

	# extract variables
	lon,lat = ncdata['lon'].values, ncdata['lat'].values
	lon[lon == 0.] = np.nan
	lat[lat == 0.] = np.nan

	data = ncdata[vars].values

	return lon,lat,data

