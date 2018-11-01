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
from scipy import interpolate

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here

def view_newGrid(nx,ny):

    fig,ax = plt.subplots()

    m = oceano.make_map(ax,ulat=-22.)

    lons,lats,x,y = m.makegrid(nx,ny,returnxy=True)

    m.plot(x,y,'k',alpha=.3)
    m.plot(x.T,y.T,'k',alpha=.3)

    return lons,lats,x,y

def regularGrid(modelado,k,timestep,nx,ny):
    lon = modelado.lon.values
    lat = modelado.lat.values
    dep = modelado.depth.values
    elev = modelado.elev[timestep,:,:].values
    u   = modelado.u[timestep,k,:,:].values
    v   = modelado.v[timestep,k,:,:].values
    s   = np.sqrt(u**2 + v**2)

    lon[lon == 0.] = np.nan
    lat[lat == 0.] = np.nan
    u[u == 0.]     = np.nan
    v[v == 0.]     = np.nan

    # realizando copias das variáveis
    lon200 = lon.copy()
    lat200 = lat.copy()
    u200   = u.copy()
    v200   = v.copy()
    elev200= elev.copy()

    # tornando nan os valores em profundidades maiores que 200m
    ind = dep > 200.
    u200[ind] = np.nan
    v200[ind] = np.nan
    elev200[ind] = np.nan
    lon200[ind] = np.nan
    lat200[ind] = np.nan

    # removendo, como a Carine removeu, as 6 primeiras e ultimas linhas
    u200[:6,:]  = np.nan
    u200[-6:,:] = np.nan
    v200[:6,:]  = np.nan
    v200[-6:,:] = np.nan
    elev200[:6,:]  = np.nan
    elev200[-6:,:] = np.nan

    lonReg,latReg,x,y = view_newGrid(nx,ny)
    plt.close('all')

    # preparando variaveis de interpolacao
    ut     = np.ravel(u200)
    vt     = np.ravel(v200)

    X1,Y1  = np.ravel(lon),np.ravel(lat)
    X1     = X1[~np.isnan(vt)]
    Y1     = Y1[~np.isnan(vt)]

    ut     = ut[~np.isnan(vt)]
    vt     = vt[~np.isnan(vt)]

    points = np.array([X1,Y1]).T                # grade irregular
    uI = interpolate.griddata(points,ut,(lonReg,latReg),method='linear')
    vI = interpolate.griddata(points,vt,(lonReg,latReg),method='linear')

    sI = np.sqrt(uI**2 + vI**2)

    un = uI/sI
    vn = vI/sI

    contours = np.linspace(0,2.6,100)

    fig,ax = plt.subplots(ncols=2,figsize=(15,15))
    m = oceano.make_map(ax[0])

    m.contourf(x,y,sI,contours,cmap=cmo.cm.speed)
    m.quiver(x[::2,::2],y[::2,::2],uI[::2,::2]*100,vI[::2,::2]*100)
    m.contour(lon,lat,dep,levels=[100,200],colors='k',alpha=.4,latlon=True)
    ax[0].set_title('Regular Grid (100 pts)',fontsize=24)

    m = oceano.make_map(ax[1])
    m.contourf(lon,lat,s,contours,cmap=cmo.cm.speed,latlon=True)
    m.quiver(lon,lat,u,v,latlon=True)
    m.contour(lon,lat,dep,levels=[100,200],colors='k',alpha=.4,latlon=True)
    ax[1].set_title('Irregular Grid', fontsize=24)
    #
    # m.plot(x,y,'k',alpha=.1);
    # m.plot(x.T,y.T,'k',alpha=.1);

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

MODELO_DIR = '/media/danilo/Danilo/mestrado/ventopcse/output/'
simulacao  = 'EC2.cdf' # nome dos arquivos para analisar

nx         = 150        # numero de pontos na direção x da grade regular
ny         = 150        # numero de pontos na direção y da grade regular
timestep   = -1        # selecionar timestep para plotar
k          = 0         # selecionar qual nível sigma plotar dados


# extrair grade do modelo
modelado  = xr.open_dataset(MODELO_DIR+simulacao)

regularGrid(modelado,k,timestep,nx,ny)

#
# lon = modelado.lon.values
# lat = modelado.lat.values
# dep = modelado.depth.values
# elev = modelado.elev[timestep,:,:].values
# u   = modelado.u[timestep,k,:,:].values
# v   = modelado.v[timestep,k,:,:].values
# s   = np.sqrt(u**2 + v**2)
#
# lon[lon == 0.] = np.nan
# lat[lat == 0.] = np.nan
# u[u == 0.]     = np.nan
# v[v == 0.]     = np.nan
#
# # realizando copias das variáveis
# lon200 = lon.copy()
# lat200 = lat.copy()
# u200   = u.copy()
# v200   = v.copy()
# elev200= elev.copy()
#
# # tornando nan os valores em profundidades maiores que 200m
# ind = dep > 200.
# u200[ind] = np.nan
# v200[ind] = np.nan
# elev200[ind] = np.nan
# lon200[ind] = np.nan
# lat200[ind] = np.nan
#
# # removendo, como a Carine removeu, as 6 primeiras e ultimas linhas
# u200[:6,:]  = np.nan
# u200[-6:,:] = np.nan
# v200[:6,:]  = np.nan
# v200[-6:,:] = np.nan
# elev200[:6,:]  = np.nan
# elev200[-6:,:] = np.nan
#
# lonReg,latReg,x,y = view_newGrid(nx,ny)
# plt.close('all')
#
# # preparando variaveis de interpolacao
# ut     = np.ravel(u200)
# vt     = np.ravel(v200)
#
# X1,Y1  = np.ravel(lon),np.ravel(lat)
# X1     = X1[~np.isnan(vt)]
# Y1     = Y1[~np.isnan(vt)]
#
# ut     = ut[~np.isnan(vt)]
# vt     = vt[~np.isnan(vt)]
#
# points = np.array([X1,Y1]).T                # grade irregular
# uI = interpolate.griddata(points,ut,(lonReg,latReg),method='linear')
# vI = interpolate.griddata(points,vt,(lonReg,latReg),method='linear')
#
# sI = np.sqrt(uI**2 + vI**2)
#
# un = uI/sI
# vn = vI/sI
#
# contours = np.linspace(0,2.6,100)
#
# fig,ax = plt.subplots(ncols=2,figsize=(15,15))
# m = oceano.make_map(ax[0])
#
# m.contourf(x,y,sI,contours,cmap=cmo.cm.speed)
# m.quiver(x[::2,::2],y[::2,::2],un[::2,::2]*100,vn[::2,::2]*100)
# m.contour(lonMod,latMod,dep,levels=[100,200],colors='k',alpha=.4,latlon=True)
# ax[0].set_title('Regular Grid (100 pts)',fontsize=24)
#
# m = oceano.make_map(ax[1])
# m.contourf(lon,lat,s,contours,cmap=cmo.cm.speed,latlon=True)
# m.quiver(lon,lat,u,v,latlon=True)
# m.contour(lonMod,latMod,dep,levels=[100,200],colors='k',alpha=.4,latlon=True)
# ax[1].set_title('Irregular Grid', fontsize=24)
# #
# # m.plot(x,y,'k',alpha=.1);
# # m.plot(x.T,y.T,'k',alpha=.1);
