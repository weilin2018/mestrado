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

    m = Basemap(projection='merc', llcrnrlat=llat, urcrnrlat=ulat, llcrnrlon=llon, urcrnrlon=ulon, resolution='i')
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
        U,V = u[j,:].values,v[j,:].values
        angle = ang[j,:]

        INT,DIR = decomp.uv2intdir(U,V,0,0)
        uro,vro = decomp.intdir2uv(INT,DIR,0,angle)
        ur[j,:] = uro
        vr[j,:] = vro

    return ur,vr

def tratando_corrente(u,v,depth,angle):

    ur,vr = rotate_velocityField(u,v,angle)
    spd = np.sqrt(ur**2+vr**2)
    spd = np.where(depth < 100, spd,np.nan)

    return ur,vr,spd

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
FILE_DIR = BASE_DIR+'masterThesis_analysis/routines/index_list.npy'
SAVE_DIR = BASE_DIR + 'masterThesis_analysis/figures/experiments_outputs/velocity/'
fname = DATA_DIR + 'EA2.cdf'

nstep = 303

# working with the data
lon,lat,u,v,depth,angles = export_data(fname,timestep=nstep)

# rotating vectors
ur = np.zeros(u.shape) * np.nan
vr = np.zeros(u.shape) * np.nan
spd= np.zeros(u.shape) * np.nan

for sigma in range(u.shape[0]):
    urot,vrot,spdrot = tratando_corrente(u[sigma,:,:],v[sigma,:,:],depth,(-1)*angles)
    ur[sigma,:,:] = urot
    vr[sigma,:,:] = vrot
    spd[sigma,:,:]= spdrot


# calculando a variancia integrada na profundidade
ur_var = np.nanvar(ur,axis=0)
vr_var = np.nanvar(vr,axis=0)

# mascarando dados em profundidades maior que 100
maskCondition = np.greater(depth,100)
masked_ur     = np.ma.masked_where(maskCondition,ur_var)
masked_vr     = np.ma.masked_where(maskCondition,vr_var)

# formatando os dados para um formato de visualizacao melhor e passando os vetores
# normalizados pela velocidade
# xplot,yplot,uplot,vplot = ocplt.formatting_vectors(ur_var,vr_var,lon,lat,FILE_DIR)

fig,ax = plt.subplots(ncols=2,figsize=(15./2.54,8./2.54))
ax[0].set_title('cross shore variance',fontsize=8)
ax[1].set_title('along shore variance',fontsize=8)
contours_u = np.arange(0,0.09,0.001)
contours_v = np.arange(0,0.1,0.01)

m1,_,_ = make_map(ax[0])

cf1 = m1.contourf(lon,lat,masked_ur,contours_u,cmap='YlOrBr',latlon=True)
cs1 = m1.contour(lon,lat,depth,levels=[80.,100.],latlon=True,colors=('black'),linewidths=(0.5))

m2,_,_ = make_map(ax[1])

cf2 = m2.contourf(lon,lat,masked_vr,contours_v,cmap='YlOrBr',latlon=True)
cs2 = m2.contour(lon,lat,depth,levels=[80.,100.],latlon=True,colors=('black'),linewidths=(0.5))
plt.tight_layout()

plt.savefig(SAVE_DIR + 'variancia.eps')
