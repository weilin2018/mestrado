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


def gradient_sphere(f, *varargs):
    """
    Return the gradient of a 2-dimensional array on a sphere given a latitude
    and longitude vector.
    The gradient is computed using central differences in the interior
    and first differences at the boundaries. The returned gradient hence has
    the same shape as the input array.
    Parameters
    ----------
    f : A 2-dimensional array containing samples of a scalar function.
    latvec: latitude vector
    lonvec: longitude vector
    Returns
    -------
    g : dfdx and dfdy arrays of the same shape as `f` giving the derivative of `f` with
        respect to each dimension.
    Examples
    --------
    temperature = temperature(pressure,latitude,longitude)
    levs = pressure vector
    lats = latitude vector
    lons = longitude vector
    >>> tempin = temperature[5,:,:]
    >>> dfdlat, dfdlon = gradient_sphere(tempin, lats, lons)
    >>> dfdp, dfdlat, dfdlon = gradient_sphere(temperature, levs, lats, lons)
    based on gradient function from /usr/lib64/python2.6/site-packages/numpy/lib/function_base.py
    """

    R_earth = 6371200.
    N = len(f.shape)  # number of dimensions
    n = len(varargs)
    argsin = list(varargs)

    if N != n:
       raise SyntaxError("dimensions of input array must match the number of remaining arguments")

    df = np.gradient(f)

    if n == 1:
        lats = argsin[0]
        dfdy = df[0]
    elif n == 2:
        lats = argsin[0]
        lons = argsin[1]
        dfdy = df[0]
        dfdx = df[1]
    elif n == 3:
        levs = argsin[0]
        lats = argsin[1]
        lons = argsin[2]
        dfdz = df[0]
        dfdy = df[1]
        dfdx = df[2]
    else:
        raise SyntaxError(
                "invalid number of arguments")

    otype = f.dtype.char
    if otype not in ['f', 'd', 'F', 'D']:
        otype = 'd'

    latarr = np.zeros_like(f).astype(otype)
    lonarr = np.zeros_like(f).astype(otype)
    if N == 1:
       nlat = np.shape(f)
       for jj in range(0,nlat):
          latarr[jj,ii] = lats[jj]
       lonarr = latarr
       lons = lats
    elif N == 2:
       nlat, nlon = np.shape(f)
       for jj in range(0,nlat):
          for ii in range(0,nlon):
             latarr[jj,ii] = lats[jj]
             lonarr[jj,ii] = lons[ii]
    else:
       nz, nlat, nlon = np.shape(f)
       for kk in range(0,nz):
          for jj in range(0,nlat):
             for ii in range(0,nlon):
                latarr[kk,jj,ii] = lats[jj]
                lonarr[kk,jj,ii] = lons[ii]

    latrad = latarr*(pi/180)

    # use central differences on interior and first differences on endpoints

    outvals = []

    dlats = np.zeros_like(lats).astype(otype)
    dlats[1:-1] = (lats[2:] - lats[:-2])
    dlats[0] = (lats[1] - lats[0])
    dlats[-1] = (dlats[-2] - dlats[-1])

    dlons = np.zeros_like(lons).astype(otype)
    dlons[1:-1] = (lons[2:] - lons[:-2])
    dlons[0] = (lons[1] - lons[0])
    dlons[-1] = (dlons[-2] - dlons[-1])

    # Since we differenced in the reverse direction, change the sign
    #dlats = -1*dlats

    dlatarr = np.tile(dlats,[nlon,1])
    dlatarr = np.reshape(dlatarr,[nlat,nlon])

    dlonarr = np.zeros_like(f).astype(otype)
    if N==2:
       for jj in range(0,nlat):
          for ii in range(0,nlon):
              dlonarr[jj,ii] = dlons[ii]
    elif N==3:
       for kk in range(0,nz):
          for jj in range(0,nlat):
             for ii in range(0,nlon):
                 dlonarr[kk,jj,ii] = dlons[ii]

    dlatsrad = dlatarr*(pi/180)
    dlonsrad = dlonarr*(pi/180)
    latrad = latarr*(pi/180)

    if n==1:
       dx1 = R_earth * dlatsrad
       dfdy = dfdy/dx1

       return dfdy
    elif n==2:
       dx1 = R_earth * dlatsrad
       dx2 = R_earth * np.cos(latrad) * dlonsrad

       dfdy = dfdy/dx1
       dfdx = dfdx/dx2

       return dfdy,dfdx
    elif n==3:
       dx1 = R_earth * dlatsrad
       dx2 = R_earth * np.cos(latrad) * dlonsrad

       dfdy = dfdy/dx1
       dfdx = dfdx/dx2

       nzz = np.shape(levs)
       if not nzz:
          nzz=0

       if nzz>1:
           zin = levs
           dz = np.zeros_like(zin).astype(otype)
           dz[1:-1] = (zin[2:] - zin[:-2])/2
           dz[0] = (zin[1] - zin[0])
           dz[-1] = (zin[-1] - zin[-2])
           if zin[0,1,1] > zin[1,1,1]:
               dz = dz*-1 # assume the model top is the first index and the lowest model is the last index

           dx3 = np.ones_like(f).astype(otype)
           for kk in range(0,nz):
               dx3[kk,:,:] = dz[kk]
       else:
           dx3 = np.ones_like(f).astype(otype)
           dx3[:] = dx[0]

       dfdz = dfdz/dx3

       return dfdz,dfdy,dfdx



def relativeVorticity(ncin,nstep):
    u = np.nanmean(ncin.u[nstep,:,:,:],axis=0)
    v = np.nanmean(ncin.v[nstep,:,:,:],axis=0)
    h1= ncin.h1.values.copy()
    h2= ncin.h2.values.copy()
    lon,lat = ncin.lon.values.copy(),ncin.lat.values.copy()
    lon[lon == 0.] = np.nan
    lat[lat == 0.] = np.nan

    dx = np.gradient(h1)
    dy = np.gradient(h2)
    du = np.gradient()

    return zeta

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
ncin = xr.open_dataset(fname)
