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

BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/data/')

# importing laje de santos data
raw = xr.open_dataset(DATA_DIR+'Est_lajeSantos/lajesantos.nc')
raw = raw.to_dataframe()

# cut only a period
raw = raw['2015-04':]
treat = raw['2015-04':].copy()
treat[treat > 3*treat.std()] = np.nan

std = treat.wind_along.std()

fig,ax = plt.subplots()

raw.wind_along.plot(ax=ax)
ax.axhline(y=3*std,c='k',ls='dashed')
ax.axhline(y=-3*std,c='k',ls='dashed')

plt.savefig('/media/danilo/Danilo/mestrado/gitlab/mestrado/dissertacao/presentation/figures/qualityControl.png')
