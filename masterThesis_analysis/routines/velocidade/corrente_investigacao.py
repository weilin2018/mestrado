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
import matplotlib.gridspec as gridspec

import matplotlib
# matplotlib.style.use('ggplot')
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

def export_data(fname,timestep=0,convertData=False,outputFile=None):
    # plotting climatologic data: t = 0, k = 0
    ncin = xr.open_dataset(fname)

    lon,lat = ncin.lon.values.copy(), ncin.lat.values.copy()
    angle = ncin.ang.values.copy()
    depth = ncin.depth.values.copy()
    sigma = ncin.sigma.values.copy()
    lon[lon == 0.] = np.nan
    lat[lat == 0.] = np.nan

    # transform data from sigma vertical coordinates to standard level
    if convertData:
	os.system('converting u and v for %s in timestep %i'%(exp,np.nanmean(timestep)))
        # interpolating data from sigma to standard level
        u,v = ncin.u[timestep,:,:,:],ncin.v[timestep,:,:,:]
        stdl,newu = oceano.sigma2stdl(u,sigma,23,depth,ncin.h1,lon,lat,'u component')
        umean = np.nanmean(newu,axis=0) # media no tempo
        np.save(outputFile+'%s_umean_%i.npy'%(exp,np.nanmean(timestep)),umean)

        stdl,newv = oceano.sigma2stdl(v,sigma,23,depth,ncin.h1,lon,lat,'v component')
        vmean = np.nanmean(newv,axis=0) # media no tempo
        np.save(outputFile+'%s_vmean_%i.npy'%(exp,np.nanmean(timestep)),vmean)

        u,v = umean,vmean
    else:
        u = np.load(outputFile+'%s_umean_%i.npy'%(exp,np.nanmean(timestep)))
        v = np.load(outputFile+'%s_vmean_%i.npy'%(exp,np.nanmean(timestep)))

        stdl = [0, 10, 25, 40, 50, 60, 70, 80, 100, 150, 300, 350, 400, 450, 500, 600, 700, 800, 1000, 1200, 1500, 1800, 2000]

    return lon,lat,u,v,depth,angle,stdl

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

    # ur,vr = rotate_velocityField(u,v,angle)
    spd = np.sqrt(u**2+v**2)
    spd = np.where(depth < 200, spd,np.nan)

    return u,v,spd

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
"""
    Importantissimo:

    Na eventualidade de rodar as simulações novamente e na ocasião de ser com novas
    configurações, é necessário converter os dados do modelo de coordenadas
    sigmas para standard levels novamente, atualizando assim os arquivos utilizados
    nesta rotina para plotagem.

    Importante manter isso em mente, pois se não for feito isso, a rotina continuará
    plotando os dados antigos.

    Para que isso ocorra, ative a seguinte opção na linha 101:
        convertData = True

"""

# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
FILE_DIR = BASE_DIR+'masterThesis_analysis/routines/index_list.npy'
os.system('clear')

convertData = False

exp = raw_input('Digite o experimento a ser plotado: ')
fname = DATA_DIR + exp +'.cdf'

plt.ion()

timestep = [np.arange(65,73,1),np.arange(280,289,1)]

nstep = timestep[1]
nstepWind = np.arange(272,280,1)
#
ncin = xr.open_dataset(fname)

lon,lat = ncin.lon.values.copy(), ncin.lat.values.copy()
angle = ncin.ang.values.copy()
depth = ncin.depth.values.copy()
sigma = ncin.sigma.values.copy()
lon[lon == 0.] = np.nan
lat[lat == 0.] = np.nan

# extracting current (sigma layer 0 - surface)
usurf = np.nanmean(ncin.u[nstep,0,:,:],axis=0)
vsurf = np.nanmean(ncin.v[nstep,0,:,:],axis=0)
ssurf = np.sqrt(usurf**2 + vsurf**2)

xplot,yplot,uplot,vplot = ocplt.formatting_vectors(usurf/ssurf,vsurf/ssurf,lon,lat,FILE_DIR)
ssurf = np.sqrt(uplot**2 + vplot**2)

# extracting wind
uwind = np.nanmean(ncin.wu[nstepWind,:,:],axis=0)/1.6
vwind = np.nanmean(ncin.wv[nstepWind,:,:],axis=0)/1.6
swind = np.sqrt(uwind**2 + vwind**2)

xplot,yplot,wuplot,wvplot = ocplt.formatting_vectors(uwind/swind,vwind/vwind,lon,lat,FILE_DIR)
swind = np.sqrt(wuplot**2 + wvplot**2)

# visualization
fig,axes = plt.subplots(ncols=2)
m0,_,_ = make_map(axes[0],ulon=np.nanmax(lon)-.5,llon=np.nanmin(lon)-.2,ulat=np.nanmax(lat)+.2,llat=np.nanmin(lat))
m1,_,_ = make_map(axes[1],ulon=np.nanmax(lon)-.5,llon=np.nanmin(lon)-.2,ulat=np.nanmax(lat)+.2,llat=np.nanmin(lat))

# view current on axes[0]
m0.contourf(lon,lat,ssurf,np.arange(0,1.5,0.1),cmap=cmo.cm.speed,latlon=True)
m0.quiver(xplot,yplot,uplot,vplot,scale=60,width=0.0015,headwidth=4,headlength=4,latlon=True)

# view wind on axes[1]
m1.contourf(lon,lat,swind,cmap=cmo.cm.speed,latlon=True)
m1.quiver(xplot,yplot,wuplot,wvplot,scale=60,width=0.0015,headwidth=4,headlength=4,latlon=True)
m1.ax.set_title('1 dia antes')
