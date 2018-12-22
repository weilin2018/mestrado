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
import cmocean as cmo

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano
import decomp

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def formatGrid_plot(grid,fname):
    import numpy as np
    ij=np.load(fname)
    # for a 2D array (lon,lat)
    if len(grid.shape)==2:
        grid=grid[ij[1], :]
        grid=grid[:, ij[0]]
    # if grid is a 3D array (temp,salt,speed)
    if len(grid.shape)==3:
        grid=grid[:,ij[1], ij[0]]
    return grid


def rotate_velocityField(u,v,ang):

    ur = np.zeros(u.shape)*np.nan
    vr = np.zeros(v.shape)*np.nan

    for j in range(u.shape[0]):
        U,V = u[j,:].values,v[j,:].values
        angle = ang[j,:]

        INT,DIR = decomp.uv2intdir(U,V,0,angle)
        uro,vro = decomp.intdir2uv(INT,DIR,0,angle)
        ur[j,:] = uro
        vr[j,:] = vro

    return ur,vr



##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github','ventopcse/output')

experimento = 'EA1.cdf'
fname = DATA_DIR + experimento

# timestep
t = -1

ncin = xr.open_dataset(fname)
# extraindo grid
lon = ncin.lon.values.copy()
lat = ncin.lat.values.copy()
ang = ncin.ang.values.copy()
lon[lon == 0.] = np.nan
lat[lat == 0.] = np.nan

# extraindo velocidades na superficie
u  = ncin.u[t,0,:,:].copy()
v  = ncin.v[t,0,:,:].copy()
wu = ncin.wu[t,:,:].copy()
wv = ncin.wv[t,:,:].copy()

# rotacionar vetores de acordo com o proprio angulo da celula
ur,vr = rotate_velocityField(u[:,:],v[:,:],ang)

s = np.sqrt(ur**2 + vr**2)

# normalizando os vetores
un = ur/s
vn = vr/s

# criando variaveis para plot mais organizado
xplot = formatGrid_plot(lon,'/media/danilo/Danilo/mestrado/github/masterThesis_analysis/routines/index_list.npy')
yplot = formatGrid_plot(lat,'/media/danilo/Danilo/mestrado/github/masterThesis_analysis/routines/index_list.npy')
uplot = formatGrid_plot(un[:,:],'/media/danilo/Danilo/mestrado/github/masterThesis_analysis/routines/index_list.npy')
vplot = formatGrid_plot(vn[:,:],'/media/danilo/Danilo/mestrado/github/masterThesis_analysis/routines/index_list.npy')

wuplot = formatGrid_plot(wu[:,:],'/media/danilo/Danilo/mestrado/github/masterThesis_analysis/routines/index_list.npy')
wvplot = formatGrid_plot(wv[:,:],'/media/danilo/Danilo/mestrado/github/masterThesis_analysis/routines/index_list.npy')

fig,ax = plt.subplots()
m = oceano.make_map(ax)

m.contourf(lon,lat,s[:,:],latlon=True,cmap=cmo.cm.speed)
m.quiver(xplot[:,::3],yplot[:,::3],uplot[:,::3],vplot[:,::3],scale=70,width=0.001,headwidth=4,headlength=4,latlon=True)
# m.quiver(xplot[::2,::4],yplot[::2,::4],wuplot[::2,::4],wvplot[::2,::4],latlon=True,color='k', alpha=0.3, scale=150,width=0.005,pivot='middle')
m.quiver(xplot[::2,::4],yplot[::2,::4],wuplot[::2,::4],wvplot[::2,::4],latlon=True,color='k',alpha=.4,pivot='middle',headwidth=4,headlength=4,minshaft=2)
