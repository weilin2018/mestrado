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

    ind = np.asarray(ind)[0]

    # vetores para separar i e j para facilitar resgatar os dados de concentração
    iss=[]
    jss=[]

    # se tiver retornado mais de um ponto, entao ind tera um shape > 1
    if ind.ndim > 2:
        for dim in range(ind.shape[0]):
            for i,j in ind[dim]:
                iss.append(int(i))
                jss.append(int(j))
    else:
        for i,j in ind:
            iss.append(int(i))
            jss.append(int(j))

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
##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

# define PCSE file and BIG file
fname_pcse = '/media/danilo/Danilo/mestrado/ventopcse/output/cut_exp06.cdf'
fname_big  = '/media/danilo/Danilo/mestrado/artigo_data/simulacoes/2ndPhase/run01_validation/run01_validation.cdf'

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

# visualize both grid and defining boundaries, in coordinates, for
# big model_grid
westBC,eastBC,southBC = visualize_grids(ncin_pcse,ncin_big)
