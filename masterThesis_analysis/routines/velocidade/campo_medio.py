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
matplotlib.use('PS')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano
import masterThesisPack.plots as ocplt

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def make_map(ax,llat=-30,ulat=-20,llon=-50,ulon=-39,resolution='l',nmeridians=3,nparallels=2,labels=[True,False,False,True]):
    # criar mapa sem continente colorido, somente a linha de costa
    # nao colocar meridianos e paralelos
    # nao colocar coordenadas no eixo (sumir com os eixos)

    m = Basemap(projection='merc', llcrnrlat=llat, urcrnrlat=ulat, llcrnrlon=llon, urcrnrlon=ulon, resolution='h')
    # m = pickle.load(open('pickles/basemap.p','r'))
    m.ax = ax
    m.drawcoastlines(linewidth=.2)
    m.fillcontinents(color='white',alpha=0)

    meridians=np.arange(llon,ulon,nmeridians)
    parallels=np.arange(llat,ulat,nparallels)

    return m,meridians,parallels

def export_data(fname,timestep=0):
    # plotting climatologic data: t = 0, k = 0
    ncin = xr.open_dataset(fname)

    lon,lat = ncin.lon.values.copy(), ncin.lat.values.copy()
    angle = ncin.ang.values.copy()
    depth = ncin.depth.values.copy()
    sigma = ncin.sigma.values.copy()
    lon[lon == 0.] = np.nan
    lat[lat == 0.] = np.nan

    # extracting temperature data, in a specific timestep
    u,v = ncin.u[timestep,:,:,:],ncin.v[timestep,:,:,:]
    # spd = np.sqrt(u**2 + v**2)
    # spd = np.where(depth < 100, spd,np.nan)

    return lon,lat,u,v,depth,angle

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
    spd = np.where(depth < 100, spd,np.nan)

    return ur,vr,spd

def calculateMeanvelocity(u,v,depth,angles):
    umean = np.nanmean(u,axis=0)
    vmean = np.nanmean(v,axis=0)
    umean,vmean,smean = tratando_corrente(umean,vmean,depth,angles)

    smean = np.where(depth<100,smean,np.nan)

    return umean,vmean,smean




##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
FILE_DIR = BASE_DIR+'masterThesis_analysis/routines/index_list.npy'
experiment = raw_input('Digite o experimento: ')
fname = DATA_DIR + experiment + '.cdf'
plt.ion()

# importing general variables to create figures
ncin = xr.open_dataset(fname)
ncin2 = xr.open_dataset(fname.replace('C','A'))

lon,lat = ncin.lon.values.copy(), ncin.lat.values.copy()
angle = ncin.ang.values.copy()
depth = ncin.depth.values.copy()
sigma = ncin.sigma.values.copy()
lon[lon == 0.] = np.nan
lat[lat == 0.] = np.nan

fig,axes = plt.subplots(nrows=2,ncols=2,figsize=(25.4/2.54, 20/2.54))

axes[0,0].set_title('%s - 14/Jan'%(experiment),fontsize=8)
axes[1,0].set_title('%s - 15/Fev'%(experiment),fontsize=8)
axes[0,1].set_title('%s - 14/Jan'%(experiment.replace('C','A')),fontsize=8)
axes[1,1].set_title('%s - 15/Fev'%(experiment.replace('C','A')),fontsize=8)

m1,meridians,parallels = make_map(axes[0,0],ulon=np.nanmax(lon)-.5,llon=np.nanmin(lon)-.2,ulat=np.nanmax(lat)+.2,llat=np.nanmin(lat))
m2,_,_ = make_map(axes[1,0],ulon=np.nanmax(lon)-.5,llon=np.nanmin(lon)-.2,ulat=np.nanmax(lat)+.2,llat=np.nanmin(lat))
m3,_,_ = make_map(axes[0,1],ulon=np.nanmax(lon)-.5,llon=np.nanmin(lon)-.2,ulat=np.nanmax(lat)+.2,llat=np.nanmin(lat))
m4,_,_ = make_map(axes[1,1],ulon=np.nanmax(lon)-.5,llon=np.nanmin(lon)-.2,ulat=np.nanmax(lat)+.2,llat=np.nanmin(lat))

# plotting EC1 - 14/Jan
u = np.nanmean( ncin.u[np.arange(48,57,1),:,:,:].values,axis=0)
v = np.nanmean( ncin.v[np.arange(48,57,1),:,:,:].values,axis=0)
umean,vmean,smean = calculateMeanvelocity(u,v,depth,angle)
xplot,yplot,uplot,vplot = ocplt.formatting_vectors(umean/smean,vmean/smean,lon,lat,FILE_DIR)

cf1 = m1.contourf(lon,lat,smean,latlon=True,cmap=cmo.cm.speed,rasterized=True)
qv1 = m1.quiver(xplot,yplot,uplot,vplot,scale=60,width=0.0015,headwidth=4,headlength=4,latlon=True)

# plotting EC1 - 15/Fev
u = np.nanmean(ncin.u[np.arange(280,289,1),:,:,:].values,axis=0)
v = np.nanmean(ncin.v[np.arange(280,289,1),:,:,:].values,axis=0)
umean,vmean,smean = calculateMeanvelocity(u,v,depth,angle)
xplot,yplot,uplot,vplot = ocplt.formatting_vectors(umean/smean,vmean/smean,lon,lat,FILE_DIR)

cf2 = m2.contourf(lon,lat,smean,latlon=True,cmap=cmo.cm.speed,rasterized=True)
qv2 = m2.quiver(xplot,yplot,uplot,vplot,scale=60,width=0.0015,headwidth=4,headlength=4,latlon=True)

# plotting EA1 - 14/Jan
u = np.nanmean(ncin2.u[np.arange(48,57,1),:,:,:].values,axis=0)
v = np.nanmean(ncin2.v[np.arange(48,57,1),:,:,:].values,axis=0)
umean,vmean,smean = calculateMeanvelocity(u,v,depth,angle)
xplot,yplot,uplot,vplot = ocplt.formatting_vectors(umean/smean,vmean/smean,lon,lat,FILE_DIR)

cf3 = m3.contourf(lon,lat,smean,latlon=True,cmap=cmo.cm.speed,rasterized=True)
qv3 = m3.quiver(xplot,yplot,uplot,vplot,scale=60,width=0.0015,headwidth=4,headlength=4,latlon=True)

# plotting EA1 - 15/Fev
u = np.nanmean(ncin2.u[np.arange(280,289,1),:,:,:].values,axis=0)
v = np.nanmean(ncin2.v[np.arange(280,289,1),:,:,:].values,axis=0)
umean,vmean,smean = calculateMeanvelocity(u,v,depth,angle)
xplot,yplot,uplot,vplot = ocplt.formatting_vectors(umean/smean,vmean/smean,lon,lat,FILE_DIR)

cf4 = m4.contourf(lon,lat,smean,latlon=True,cmap=cmo.cm.speed,rasterized=True)
qv4 = m4.quiver(xplot,yplot,uplot,vplot,scale=60,width=0.0015,headwidth=4,headlength=4,latlon=True)
