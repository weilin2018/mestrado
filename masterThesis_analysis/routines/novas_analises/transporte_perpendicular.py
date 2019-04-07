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
# matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def rotate_velocityField(u,v,ang):

    import decomp
    ur = np.zeros(u.shape)*np.nan
    vr = np.zeros(v.shape)*np.nan

    for j in range(u.shape[0]):
        U,V = u[j,:],v[j,:]
        angle = ang[j,:]

        INT,DIR = decomp.uv2intdir(U,V,0,angle)
        uro,vro = decomp.intdir2uv(INT,DIR,0,angle)
        ur[j,:] = uro
        vr[j,:] = vro

    return ur,vr

def tratando_corrente(u,v,depth,angle):

    ur,vr = rotate_velocityField(u,v,angle)
    spd = np.sqrt(ur**2+vr**2)
    # spd = np.where(depth < 200, spd,np.nan)

    return ur,vr,spd

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
plt.ion()

ncin = xr.open_dataset(DATA_DIR + 'EA1.cdf')

lon = ncin.lon.values.copy()
lat = ncin.lat.values.copy()
dep = ncin.depth.values.copy()
ang = ncin.ang.values.copy()


lon[ lon == 0. ] = np.nan
lat[ lat == 0. ] = np.nan

# limitando a regiao para extrair dados
jBegin, jFinal, i = 5,125,65

os.system('clear')
view = input('visualizar o plot? [0 para não, 1 para sim]')

if view:
    fig,ax = plt.subplots()
    m = oceano.make_map(ax,resolution='c')

    m.plot(lon,lat,'k',alpha=.3,latlon=True);
    m.plot(lon.T,lat.T,'k',alpha=.3,latlon=True);

    m.plot(lon[jBegin:jFinal,i],lat[jBegin:jFinal,i],'r',latlon=True);


# extraindo dados de velocidade (perpendicular, portanto precisa rotacionar)
u,v = ncin.u[0,:,:,:], ncin.v[0,:,:,:]
ur,vr,sr = tratando_corrente(u,v,dep,(-1)*ang)
