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
from scipy.spatial import cKDTree
from matplotlib import dates
import datetime

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here

# encontrar indices dos pontos mais proximo a uma coordenada
def find_nearest(lon,lat,ilon,ilat,k=2):
    '''
        lon,lat = lat e lon da grade
        ilon,ilat = ponto a ser encontrado
    '''

    # localizacao do terminal da ilha guaiba

    lo = lon.ravel()
    la = lat.ravel()

    coords = []

    for i,j in zip(la,lo):
        coords.append([i,j])

    coords = np.array(coords)

    locations_posi = [[ilat,ilon]]

    locs = np.asarray(locations_posi)

    tree = cKDTree(coords)
    # procura em tree os pontos mais próximos dos pontos definidos acima
    dists,indexes = tree.query(locs,k=k)

    pontos = []

    for index in indexes:
        pontos.append(coords[index])

    # converter de lista para array
    pontos = np.asarray(pontos)[0]

    # findind indexes from lat and lon
    ind = []

    for p in pontos:
        ind.append(np.where(lon == p[1]))

    ind = np.asarray(ind)

    # vetores para separar i e j para facilitar resgatar os dados de concentração
    iss=[]
    jss=[]

    # if k > 1
    if k > 1:
        for dim in range(k):
            iss.append(int(ind[dim][0]))
            jss.append(int(ind[dim][1]))


    # # se tiver retornado mais de um ponto, entao ind tera um shape > 1
    # if ind.ndim > 2:
    #     for dim in range(ind.shape[0]):
    #         for i,j in ind[dim]:
    #             iss.append(int(i))
    #             jss.append(int(j))
    # else:
    #     for i,j in ind:
    #         iss.append(int(i))
    #         jss.append(int(j))
    #
    return iss,jss

def extract_Coordinates(ncin):
    lat = ncin.lat.values.copy()
    lon = ncin.lon.values.copy()

    # remove zeros
    lat[lat == 0.0] = np.nan
    lon[lon == 0.0] = np.nan

    return lon,lat

def getBoundaries(lon,lat):

    # only for BIG model_grid
    westBC  = [lon[2,:],lat[2,:]]
    eastBC  = [lon[-2,:],lat[-2,:]]
    southBC = [lon[:,-2],lat[:,-2]]

    return westBC,eastBC,southBC

def visualize_grids(pcse,big):

    fig,ax = plt.subplots()

    m = oceano.make_map(ax,resolution='i')

    # import coordinates from coarser grid
    lon,lat = extract_Coordinates(pcse)
    m.plot(lon,lat,'k',latlon=True,alpha=.3)
    m.plot(lon.T,lat.T,'k',latlon=True,alpha=.3)

    # import coordinates from refined grid
    lon,lat = extract_Coordinates(big)
    # get boundaries coordinates
    westBC,eastBC,southBC = getBoundaries(lon,lat)

    m.plot(westBC[0],westBC[1],'r',latlon=True)
    m.plot(eastBC[0],eastBC[1],'r',latlon=True)
    m.plot(southBC[0],southBC[1],'r',latlon=True)

    plt.show()

    return westBC,eastBC,southBC

def find_nearest_coordinates(boundary,lon,lat):
    # boundary: list with [0]: longs, [1]:latis, e.g. westBC
    # lon,lat: pcse grid

    i,j = [],[]

    for ilon,ilat in zip(boundary[0],boundary[1]):
        # select only non-NaN coordinates
        if (~np.isnan(ilon))&(~np.isnan(ilat)):
            iss,jss = find_nearest(lon,lat,ilon,ilat)
            # print(iss[0],jss[0])
            i.append(iss[0])
            j.append(jss[0])

    return np.asarray(i),np.asarray(j)

def viewGrid(lon,lat):
    fig,ax = plt.subplots()

    m = oceano.make_map(ax,resolution='i')

    m.plot(lon,lat,'k',latlon=True)
    m.plot(lon.T,lat.T,'k',latlon=True)

    plt.show()

def findNearest(lonBnd,latBnd,k=1):

    iss = []
    jss = []

    for lonref,latref in zip(lonBnd,latBnd):
        lo = lon_coarsed.ravel()
        la = lat_coarsed.ravel()

        coords = []

        for i,j in zip(la,lo):
            coords.append([i,j])

        coords = np.array(coords)

        locations_posi = [[latref,lonref]]

        locs = np.asarray(locations_posi)

        tree = cKDTree(coords)
        # procura em tree os pontos mais próximos dos pontos definidos acima
        dists,indexes = tree.query(locs,k=k)

        # searching coordinates based on the indexes
        pontos = [coords[index] for index in indexes]
        pontos = np.squeeze(np.asarray(pontos))

        # findind indexes from lat and lon
        if pontos.ndim == 1: # if there is only 1 point:
            ind = np.where(lon_coarsed == pontos[1])
        else: # for multiples points
            ind = [np.where(lon_coarsed == p[1]) for p in pontos]
        # convert list to np.ndarray
        ind = np.squeeze(np.asarray(ind))

        # vetores para separar i e j para facilitar resgatar os dados de concentração
        if ind.ndim == 1:
            iss.append(int(ind[0]))
            jss.append(int(ind[1]))
        elif ind.ndim >= 2:
            # iss = []
            # jss = []

            dims = ind.ndim
            for dim in range(dims):
                iss.append(ind[dim][0])
                jss.append(ind[dim][1])

    return iss,jss


def find(lon_coarsed,lat_coarsed,lonBnd,latBnd):
    """Find nearest points between coarsed and refined grids.

    Parameters
    ----------
    lon_coarsed : numpy.ndarray
        Longitude of coarsed grid.
    lat_coarsed : numpy.ndarray
        Latitude of coarsed grid
    lonBnd : numpy.ndarray
        Boundaries longitude of refined grid.
    latBnd : numpy.ndarray
        Boundaries latitudes of refined grid.

    Returns
    -------
    iss,jss : list
        List with indexes of coarsed grid, containing all
        nearest points to boundaries locations from refined grid.

    """

    # removing nan points from coarsed grid
    indNan_coarsed = np.where(~np.isnan(lon_coarsed))
    x,y = lon_coarsed[indNan_coarsed], lat_coarsed[indNan_coarsed]
    # removing from refined grid
    indNan = np.where(~np.isnan(lonBnd))
    xi,yi = lonBnd[indNan],latBnd[indNan]

    # creating array with coordinates
    mdgrid = np.array(zip( x.ravel(),y.ravel() ))
    # instancing KDTree object
    kdt = kd.KDTree(mdgrid)

    # locating points (distances and indexes)
    mdi_dist, mdi = kdt.query(np.array(zip(xi,yi)))

    # locating coordinates based on indexes in mdi
    points = np.asarray([mdgrid[index] for index in mdi])

    # based on points, find Is and Js
    ind    = np.asarray([np.where(lon_coarsed == p[0]) for p in points])

    # classifying all point found, separating in iss and jss, with same size
    for index in ind:
        if index[-1].shape == 1:
            iss.append(int(index[0]))
            jss.append(int(index[1]))
        else:
            iss.append(index[0][-1])
            jss.append(index[1][-1])

    return iss,jss

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR   = '/media/danilo/Danilo/mestrado/artigo_data/simulacoes/2ndPhase/'
# define PCSE file and BIG file
fname_pcse = BASE_DIR + "pcse_elevation.cdf"
fname_big  = BASE_DIR + "reRUN00.cdf"

# read netcdf file
ncin_pcse  = xr.open_dataset(fname_pcse)
ncin_big   = xr.open_dataset(fname_big)

# importing grid pcse
lon_coarsed,lat_coarsed = extract_Coordinates(ncin_pcse)
# visualizing grid
# viewGrid(lon_coarsed,lat_coarsed)

# importing grid big
lon_refined,lat_refined = extract_Coordinates(ncin_big)
# visualizing grid big
# viewGrid(lon_refined,lat_refined)

# select boundaries for refined grid
westBC  = [lon_refined[2,:],lat_refined[2,:]]
eastBC  = [lon_refined[-2,:],lat_refined[-2,:]]
southBC = [lon_refined[:,-2],lat_refined[:,-2]]

boundaries = np.asarray([westBC,southBC,eastBC])

#### --------------- finding nearest points of all boundaries created

# finding nearest points to the west boundary
i = 0# for west boundary
lonBnd,latBnd = boundaries[i][0],boundaries[i][1]

iss_west,jss_west = find(lon_coarsed,lat_coarsed,lonBnd,latBnd)

# finding nearest points to the south boundary
i = 1# for south boundary
lonBnd,latBnd = boundaries[i][0],boundaries[i][1]

iss_south,jss_south = find(lon_coarsed,lat_coarsed,lonBnd,latBnd)

# finding nearest points to the east boundary
i = 2# for east boundary
lonBnd,latBnd = boundaries[i][0],boundaries[i][1]

iss_east,jss_east = find(lon_coarsed,lat_coarsed,lonBnd,latBnd)

#### --------------- visualize both grids together
offset = 0.25                # offset for coordinates

fig,ax = plt.subplots()
m = oceano.make_map(ax,llat=np.nanmin(lat_refined)-offset, ulat=np.nanmax(lat_refined), llon=np.nanmin(lon_refined)-offset, ulon=np.nanmax(lon_refined)+offset,resolution='i')

m.plot(lon_coarsed,lat_coarsed,'k',latlon=True,alpha=.3);
m.plot(lon_coarsed.T,lat_coarsed.T,'k',latlon=True,alpha=.3);

# m.plot(lon_refined,lat_refined,'r',latlon=True,alpha=.3);
# m.plot(lon_refined.T,lat_refined.T,'r',latlon=True,alpha=.3);

# plotting boundaries of refined grid
os.system('clear')
for i in range(3):
    print('####################### Boundary: %i'%(i))
    lonBnd = boundaries[i][0]
    latBnd = boundaries[i][1]
    # m.plot(lonBnd,latBnd,'r',latlon=True,alpha=.9)
    m.scatter(lonBnd,latBnd,s=30,c='r',latlon=True)

# plotting in m
# first west boundary
for i,j in zip(iss_west,jss_west):
    m.scatter(lon_coarsed[i,j],lat_coarsed[i,j],s=50,c='g',latlon=True)

# first south boundary
for i,j in zip(iss_south,jss_south):
    m.scatter(lon_coarsed[i,j],lat_coarsed[i,j],s=50,c='g',latlon=True)

# first east boundary
for i,j in zip(iss_east,jss_east):
    m.scatter(lon_coarsed[i,j],lat_coarsed[i,j],s=50,c='g',latlon=True)

#### --------------- loading data from pcse_gcmplt.cdf and creating file
