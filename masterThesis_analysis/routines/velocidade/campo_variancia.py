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

    m = Basemap(projection='merc', llcrnrlat=llat, urcrnrlat=ulat, llcrnrlon=llon, urcrnrlon=ulon, resolution=resolution)
    # m = pickle.load(open('pickles/basemap.p','r'))
    m.ax = ax
    m.drawcoastlines(linewidth=.1)
    m.fillcontinents(color='white',alpha=0)

    meridians=np.arange(llon,ulon,nmeridians)
    parallels=np.arange(llat,ulat,nparallels)

    return m,meridians,parallels

def export_data(fname,timestep):
    # plotting climatologic data: t = 0, k = 0
    ncin = xr.open_dataset(fname)

    lon,lat = ncin.lon.values.copy(), ncin.lat.values.copy()
    angle = ncin.ang.values.copy()
    depth = ncin.depth.values.copy()
    sigma = ncin.sigma.values.copy()
    lon[lon == 0.] = np.nan
    lat[lat == 0.] = np.nan

    # extracting temperature data, in a specific timestep
    u = np.nanmean(ncin.u[timestep,:,:,:],axis=1)
    v = np.nanmean(ncin.v[timestep,:,:,:],axis=1)

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

def load_and_treat(ncin):
    # working with the data
    lon,lat,u,v,depth,angles = export_data(fname,timestep=np.arange(46,303,1))

    # calculando a variancia integrada na profundidade
    u_var = np.nanvar(u,axis=0)
    v_var = np.nanvar(v,axis=0)

    # mascarando dados em profundidades maior que 100
    maskCondition = np.greater(depth,200)
    masked_u      = np.ma.masked_where(maskCondition,u_var)
    masked_v      = np.ma.masked_where(maskCondition,v_var)

    # mascarando dados das 5 primeiras e ultimas linhas da grade
    masked_u[:5,:] = np.nan
    masked_v[:5,:] = np.nan
    masked_u[-5:,:] = np.nan
    masked_v[-5:,:] = np.nan

    return lon,lat,depth,masked_u,masked_v
##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
FILE_DIR = BASE_DIR+'masterThesis_analysis/routines/index_list.npy'
SAVE_DIR = BASE_DIR + 'masterThesis_analysis/figures/experiments_outputs/velocity/'

experiment = raw_input('Digite o experimento: ')
fname = DATA_DIR + experiment + '.cdf'
plt.ion()

# load control experiment
lon,lat,depth,ucontrol_var,vcontrol_var = load_and_treat(xr.open_dataset(fname))
# load anomalous experiment
lon,lat,depth,uanomalo_var,vanomalo_var = load_and_treat(xr.open_dataset(fname.replace('C','A')))

# generating figure's structure, with a 2x2 matrix of subplots

plt.close('all')
fig,ax = plt.subplots(ncols=2,nrows=2,figsize=(15./2.54,15./2.54))
ax[0,0].set_title('Componente u',fontsize=8)
ax[0,1].set_title('Componente v',fontsize=8)
plt.suptitle(u'Variância das componentes u e v integradas verticalmente',fontsize=10)

contours = np.arange(0,.05,.001)

# plotting control experiment
m1,_,_ = make_map(ax[0,0],ulon=-42,llon=-49,ulat=-22.3,llat=-29,resolution='i')
cf1 = m1.contourf(lon,lat,ucontrol_var,contours,cmap='YlOrBr',latlon=True)
cs1 = m1.contour(lon,lat,depth,levels=[80.],latlon=True,colors=('black'),linewidths=(0.2))
plt.clabel(cs1,[80],fmt='%i',inline=1,fontsize=8,manual=True)

m2,_,_ = make_map(ax[0,1],ulon=-42,llon=-49,ulat=-22.3,llat=-29,resolution='i')
cf2 = m2.contourf(lon,lat,vcontrol_var,contours,cmap='YlOrBr',latlon=True)
cs2 = m2.contour(lon,lat,depth,levels=[80.],latlon=True,colors=('black'),linewidths=(0.2))
plt.clabel(cs2,[80],fmt='%i',inline=1,fontsize=8,manual=True)

# plotting anomalous experiment
m1,_,_ = make_map(ax[1,0],ulon=-42,llon=-49,ulat=-22.3,llat=-29,resolution='i')
cf1 = m1.contourf(lon,lat,uanomalo_var,contours,cmap='YlOrBr',latlon=True)
cs1 = m1.contour(lon,lat,depth,levels=[80.],latlon=True,colors=('black'),linewidths=(0.2))
plt.clabel(cs1,[80],fmt='%i',inline=1,fontsize=8,manual=True)

m2,_,_ = make_map(ax[1,1],ulon=-42,llon=-49,ulat=-22.3,llat=-29,resolution='i')
cf2 = m2.contourf(lon,lat,vanomalo_var,contours,cmap='YlOrBr',latlon=True)
cs2 = m2.contour(lon,lat,depth,levels=[80.],latlon=True,colors=('black'),linewidths=(0.2))
plt.clabel(cs2,[80],fmt='%i',inline=1,fontsize=8,manual=True)

plt.tight_layout()
plt.subplots_adjust(top=0.907,bottom=0.031,left=0.025,right=0.975,hspace=0.065,wspace=0.0)

# inserting colorbar
caxes1 = fig.add_axes([.12,.1,.3,.03])
caxes2 = fig.add_axes([.62,.1,.3,.03])

cbar1 = plt.colorbar(cf1,cax=caxes1,orientation='horizontal')
cbar2 = plt.colorbar(cf2,cax=caxes2,orientation='horizontal')

# setting colorbar tick labels
from matplotlib import ticker
tick_locator = ticker.MaxNLocator(nbins=6)
cbar1.locator = tick_locator
cbar1.update_ticks()
cbar2.locator = tick_locator
cbar2.update_ticks()

cbar1.ax.axes.tick_params(axis='both',which='both',labelsize=8)
cbar1.ax.set_title(u'Variância ['+r'm$^{2}$.s$^{-2}$'+' ]',fontsize=8)
cbar2.ax.axes.tick_params(axis='both',which='both',labelsize=8)
cbar2.ax.set_title(u'Variância ['+r'm$^{2}$.s$^{-2}$'+' ]',fontsize=8)

#plt.savefig(SAVE_DIR + 'variancia.eps')
