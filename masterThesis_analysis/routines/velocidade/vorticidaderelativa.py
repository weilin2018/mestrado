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
    # m.drawcoastlines(linewidth=.2)
    # m.fillcontinents(color='white',alpha=0)
    #
    # meridians=np.arange(llon,ulon,nmeridians)
    # parallels=np.arange(llat,ulat,nparallels)

    return m

def import_interpolated_data(fname,timestep,outputFile):

    u = np.load(outputFile+'%s_umean_%i.npy'%(exp,np.nanmean(timestep)))
    v = np.load(outputFile+'%s_vmean_%i.npy'%(exp,np.nanmean(timestep)))

    stdl = [0, 10, 25, 40, 50, 60, 70, 80, 100, 150, 300, 350, 400, 450, 500, 600, 700, 800, 1000, 1200, 1500, 1800, 2000]

    return u,v,stdl

#!-----------------------------------------------!#
def createNewgrid(jr,ir,lon,lat):
    plt.ioff()
    fig,ax = plt.subplots()
    m = make_map(ax)
    plt.close()
    # new grid
    lons,lats,x,y = m.makegrid(jr,ir,returnxy=True)

    plt.ion()

    return lons,lats,x,y

def modelgrid2regular(jr,ir,lon,lat,u,v,depth):

    # creating new regular grid
    lons,lats,x,y = createNewgrid(jr,ir,lon,lat)

    # preparing variables to interpolate
    ut,vt = np.ravel(u),np.ravel(v)
    X1,Y1 = np.ravel(lon),np.ravel(lat)
    inds = np.where(~np.isnan(X1))
    X1 = X1[inds]
    Y1 = Y1[inds]
    ut = ut[inds]
    vt = vt[inds]
    zt = depth.ravel()[inds]

    # create points array
    points = np.array([X1,Y1]).T
    # interpolate data
    uI = interpolate.griddata(points,ut,(lons,lats),method='linear')
    vI = interpolate.griddata(points,vt,(lons,lats),method='linear')
    zI = interpolate.griddata(points,zt,(lons,lats),method='linear')

    return uI,vI,zI,lons,lats,x,y

def relativeVorticity(uI,vI,lons,lats):
    # calculating size of grid cell
    dx = np.mean(np.diff(lons,axis=1)*111000)
    dy = np.mean(np.diff(lats,axis=0)*111000)

    # compute term from relative vorticity equation
    dvdx = np.gradient(vI,dx,axis=0)
    dudy = np.gradient(uI,dy,axis=0)
    # relative vorticity (zeta)
    zeta = dvdx - dudy

    return zeta


#!-----------------------------------------------!#
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
outputFile = DATA_DIR + 'dados_interpolados_stdl/'

os.system('clear')

convertData = False

exp = raw_input('Digite o experimento a ser plotado: ')
fname = DATA_DIR + exp +'.cdf'

plt.ion()

timestep = [np.arange(65,289,1),np.arange(280,289,1)]
ncin = xr.open_dataset(fname)
# extracting grid (lon,lat)
lon,lat = ncin.lon.values.copy(),ncin.lat.values.copy()
# extracting depth
depth= ncin.depth.values.copy()
# extracting sizes of each grid cell (h1,h2)
h1,h2 = ncin.h1.values.copy(),ncin.h2.values.copy()
# removing zeros coordinates from grid
lon[lon == 0.] = np.nan
lat[lat == 0.] = np.nan
# extracting current components already interpolated to standard levels
umean,vmean,stdl = import_interpolated_data(fname,timestep[0],outputFile)
u,v,stdl = import_interpolated_data(fname.replace('C','A'),timestep[1],outputFile)

u = umean - u
v = vmean - v

# # recreating grid, based on lon,lat, but remember:
# # model runs with arakawa c-grid, so the coordinates are given in the central
# # point of a cell and current components are given between those central points.
# dx = 0.5 * (h1[1:,1:] + h1[:-1,:-1])
# dy = 0.5 * (h2[1:,1:] + h2[:-1,:-1])
#
# dvdx = np.zeros(dx.shape)*np.nan
# dudy = np.zeros(dy.shape)*np.nan


######## traveling on each cell, computing dvdx

# select first z level, that represent sea surface
uplot,vplot = u[0,:,:],v[0,:,:]

# interpolating data into a regular grid
uI,vI,zI,lons,lats,x,y = modelgrid2regular(200,200,lon,lat,uplot,vplot,depth)

# computing relative vorticity
zeta = relativeVorticity(uI,vI,lons,lats)

# removing values above 200 meters depth and
# dividing by 10e-5 just to turn things easier to visualize#
nzeta = zeta/10e-05
zetaPlot = np.where(zI < 200, nzeta, np.nan)

# visualization
fig,ax = plt.subplots()
m = make_map(ax)
contours = np.arange(-0.4,0.6,0.001)
cf = m.contourf(x,y,zetaPlot,contours,cmap=cmo.cm.curl)

m.drawcoastlines(linewidth=.2)
m.fillcontinents()

cbar = plt.colorbar(cf)
cbar.set_label(r'Vorticidade Relativa x10$^{-5}$ [$s^-1{}$]')
#!-----------------------------------------------!#

plt.savefig('/media/danilo/Danilo/mestrado/github/masterThesis_analysis/figures/experiments_outputs/vorticity/vorticidadeRelativa_%s___novo.png'%(exp),dpi=300)
plt.close()
