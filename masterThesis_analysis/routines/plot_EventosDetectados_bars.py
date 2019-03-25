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
matplotlib.use('PS')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
FIGU_DIR = BASE_DIR + 'masterThesis_analysis/figures/'

qtd = [14,11,7,6,5,4,4,1,1,1,1]
duration = [9,10,11,12,13,14,15,16,17,19,24]

data = pd.DataFrame({'Eventos':qtd},index=duration)

fig,ax = plt.subplots()
data.plot(ax=ax,kind='bar')

props = dict(boxstyle='round', facecolor='white', alpha=0.5)
fig.text(0.90,0.38,u'Verão de 2014',rotation=70,bbox=props)

plt.xlabel(u'Duração em dias',fontsize=8)
plt.ylabel(u'Número de Eventos Detectados',fontsize=8)

ax.axes.tick_params(axis='both',which='both',labelsize=8)
ax.margins(xmargin=0)

plt.tight_layout()

plt.savefig(FIGU_DIR + 'grafico_barra_EventosDetectados.eps')
