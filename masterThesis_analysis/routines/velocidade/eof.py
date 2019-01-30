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

from eofs.standard import Eof
from eofs.multivariate.standard import MultivariateEof

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

def mask_byDepth(data,depth,reference=100):
    maskCond = np.greater(depth,reference)
    maskedData = np.ma.masked_where(maskCond,data)

    return maskedData

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
FILE_DIR = BASE_DIR+'masterThesis_analysis/routines/index_list.npy'
SAVE_DIR = BASE_DIR + 'masterThesis_analysis/figures/experiments_outputs/velocity/'
fname = DATA_DIR + 'EA2.cdf'

ncin = xr.open_dataset(fname)

# import variables and mask them based on depth > 100
depth = ncin.depth
u     = ncin.u[:,0,:,:]
v     = ncin.v[:,0,:,:]
lon   = ncin.lon.values
lat   = ncin.lat.values
lon[lon == 0.] = np.nan
lat[lat == 0.] = np.nan

# mascarando dados em profundidades maior que 100
masked_u,masked_v = np.zeros(u.shape),np.zeros(v.shape)
maskCondition = np.greater(depth.values,100)

for t in range(u.shape[0]):
    masked_u[t,:,:] = np.ma.masked_where(maskCondition,u[t,:,:])
    masked_v[t,:,:] = np.ma.masked_where(maskCondition,v[t,:,:])

# weight based on latitude
weight_array = np.sqrt(np.cos(np.deg2rad(lat)))
# create solver instance, calculating the multivariate EOF, with both variables with weighting
solver = MultivariateEof([masked_u,masked_v], weights=[weight_array, weight_array])
# calculating statistical modes
# set how many eigenvalues
eigenvalues = 100
eofs        = solver.eofs(neofs=eigenvalues)
pcs         = solver.pcs(npcs=eigenvalues,pcscaling=0)
varfrac     = solver.varianceFraction(neigs=eigenvalues)
lambdas     = solver.eigenvalues(neigs=eigenvalues)

plt.figure()
eof_num = range(1, eigenvalues)
plt.plot(eof_num, varfrac[0:eigenvalues-1], linewidth=2)
plt.plot(eof_num, varfrac[0:eigenvalues-1], linestyle='None', marker="o", color='r', markersize=8)
plt.axhline(0, color='k')
plt.xticks(range(1, eigenvalues))
# plt.title('Fraction of the total variance represented by each EOF')
plt.title(u'Fração da variância total representada por cada EOF')
plt.xlabel('EOF #')
plt.ylabel('Variance Fraction')
plt.xlim(1, eigenvalues-1)
plt.ylim(np.min(varfrac), np.max(varfrac)+0.01)
plt.show()

# for u
masked_eof4u = np.zeros(eofs[0].shape)
for i in range(3):
    masked_eof4u[i,:,:] = mask_byDepth(eofs[0][i,:,:],depth)

# for v
masked_eof4v = np.zeros(eofs[1].shape)
for i in range(3):
    masked_eof4v[i,:,:] = mask_byDepth(eofs[1][i,:,:],depth)

contours = np.arange(-0.4,0.4,0.01)

# visualizando a EOF
fig,axes = plt.subplots(ncols=3,nrows=2,figsize=(15./2.54,10./2.54))

# eof para u
m0,_,_ = make_map(axes[0,0],ulon=-42,llon=-49,ulat=-22.3,llat=-29,resolution='i')
m1,_,_ = make_map(axes[0,1],ulon=-42,llon=-49,ulat=-22.3,llat=-29,resolution='i')
m2,_,_ = make_map(axes[0,2],ulon=-42,llon=-49,ulat=-22.3,llat=-29,resolution='i')

cf = m0.contourf(lon,lat,masked_eof4u[0,:,:],contours,latlon=True)
cf = m1.contourf(lon,lat,masked_eof4u[1,:,:],contours,latlon=True)
cf = m2.contourf(lon,lat,masked_eof4u[2,:,:],contours,latlon=True)
# eof para v
m0,_,_ = make_map(axes[1,0],ulon=-42,llon=-49,ulat=-22.3,llat=-29,resolution='i')
m1,_,_ = make_map(axes[1,1],ulon=-42,llon=-49,ulat=-22.3,llat=-29,resolution='i')
m2,_,_ = make_map(axes[1,2],ulon=-42,llon=-49,ulat=-22.3,llat=-29,resolution='i')

m0.contourf(lon,lat,masked_eof4v[0,:,:],contours,latlon=True)
m1.contourf(lon,lat,masked_eof4v[1,:,:],contours,latlon=True)
cf=m2.contourf(lon,lat,masked_eof4v[2,:,:],contours,latlon=True)

cbax = fig.add_axes([0.12,0.05,0.78,.05])
cbar = plt.colorbar(cf,cax=cbax,orientation='horizontal')

cbar.ax.axes.tick_params(axis='both',which='both',labelsize=8)
# cbar.ax.set_title(u'EOF',fontsize=8)

plt.suptitle('EOF multivariada para corrente na superficie \n Linha superior para U e inferior para V',fontsize=10)

plt.tight_layout()
plt.subplots_adjust(top=0.898,bottom=0.117,left=0.027,right=0.973,hspace=0.055,wspace=0.0)
