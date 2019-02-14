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

    return lon,lat,depth,angle

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
experiment = 'EC2'
fname = DATA_DIR + experiment + '.cdf'
plt.ion()

# import global variables
lon,lat,depth,angle = export_data(fname)
ncin = xr.open_dataset(fname)

# import velocity components
umean,vmean = np.nanmean(ncin.u[:,:,:,:],axis=1),np.nanmean(ncin.v[:,:,:,:],axis=1)

# rotating vectors
ur = np.zeros(umean.shape) * np.nan
vr = np.zeros(umean.shape) * np.nan
spd= np.zeros(umean.shape) * np.nan

for t in range(umean.shape[0]):
    urot,vrot,spdrot = tratando_corrente(umean[t,:,:],vmean[t,:,:],depth,(-1)*angle)
    ur[t,:,:] = urot
    vr[t,:,:] = vrot
    spd[t,:,:]= spdrot

# calculando a variancia no temporal, da velocidade media
ur_var = np.nanvar(ur,axis=0)
vr_var = np.nanvar(vr,axis=0)


# mascarando dados em profundidades maior que 100
maskCondition = np.greater(depth,100)
masked_ur     = np.ma.masked_where(maskCondition,ur_var)
masked_vr     = np.ma.masked_where(maskCondition,vr_var)

# mascarando dados das 5 primeiras e ultimas linhas da grade
masked_ur[:5,:] = np.nan
masked_vr[:5,:] = np.nan
masked_ur[-5:,:] = np.nan
masked_vr[-5:,:] = np.nan

# formatando os dados para um formato de visualizacao melhor e passando os vetores
# normalizados pela velocidade
# xplot,yplot,uplot,vplot = ocplt.formatting_vectors(ur_var,vr_var,lon,lat,FILE_DIR)
plt.close('all')
fig,ax = plt.subplots(ncols=2,nrows=2,figsize=(15./2.54,8./2.54))
ax[0].set_title('Perpendicular',fontsize=8)
ax[1].set_title('Paralela',fontsize=8)
plt.suptitle(u'Variância das componentes perpendicular e paralela integradas verticalmente',fontsize=10)

contours = np.arange(0,.1,.001)

contours_u = np.arange(0,0.09,0.001)
contours_v = np.arange(0,0.1,0.01)

m1,_,_ = make_map(ax[0],ulon=-42,llon=-49,ulat=-22.3,llat=-29,resolution='i')
cf1 = m1.contourf(lon,lat,masked_ur,contours,cmap='YlOrBr',latlon=True)
cs1 = m1.contour(lon,lat,depth,levels=[80.],latlon=True,colors=('black'),linewidths=(0.2))
plt.clabel(cs1,[80],fmt='%i',inline=1,fontsize=8,manual=True)

m2,_,_ = make_map(ax[1],ulon=-42,llon=-49,ulat=-22.3,llat=-29,resolution='i')
cf2 = m2.contourf(lon,lat,masked_vr,contours,cmap='YlOrBr',latlon=True)
cs2 = m2.contour(lon,lat,depth,levels=[80.],latlon=True,colors=('black'),linewidths=(0.2))
plt.clabel(cs2,[80],fmt='%i',inline=1,fontsize=8,manual=True)

plt.tight_layout()
plt.subplots_adjust(top=0.86,bottom=0.018,left=0.021,right=0.969,hspace=0.0,wspace=0.0)

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
