"""
Routine to plot sea surface temperature (SST) from ECOM (model) and
GHRSST (Group for High Resolution SST) satellite data.

PS: only surface!!!!!

"""

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
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

# define which experiment to load
# trabalhar somente com os experimentos em cold start, até
# eu conseguir avaliar qual é o mais adequado
run = 'exp06'
startTime=112 # 32
endTime=353 # 272

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def load_ghrsst(fname):

    # load netcdf file
    ncdata = xr.open_dataset(fname)
    # extract data, converting from K to oC
    sst = ncdata['analysed_sst'].values - 273.15
    time = ncdata['time'].values

    return sst,time

def get_LatLong(fname):
    """Extract latitude and longitude from netcdf file.

    Parameters
    ----------
    fname : string
        Full path and name for a file.

    Returns
    -------
    lat,lon : array
        np.ndarray with latitude and longitudes.
    """

    # load file
    ncdata = xr.open_dataset(fname)
    # extract data
    lon = ncdata['lon'].values
    lat = ncdata['lat'].values
    # return
    return lon,lat

def create_matriz_ghrsst(nfiles):
    """Routine to read all files from GHRSST and create a
    3D-array with [time,lon,lat]. Also returning lat,lon gridded.

    Parameters
    ----------
    nfiles : list
        List from glob.glob with all files to read.

    Returns
    -------
    lon,lat : ndarray
        An array object with lon and lat coordinates gridded.
    all_sst : ndarray
        An 3D array object with all SST data extracted from GHRSST.
    """

    # criar matriz 3D [time,lon,lat]
    nlon = 910
    nlat = 911
    all_sst = np.zeros([len(nfiles),nlat,nlon])
    time = []

    # ensure temporal order
    nfiles.sort()

    for i,fname in enumerate(nfiles):

        if i == 0:
            lon,lat = get_LatLong(fname)

        # load sst data
        sst,t = load_ghrsst(fname)
        time.append(t)

        # store data into a 3D matrix
        all_sst[i,:,:] = sst

    # regrid lon,lat
    lon,lat = np.meshgrid(lon,lat)
    time = np.asarray(time)

    # return
    return lon,lat,time,all_sst

def load_ecom(fname):

    ncdata = xr.open_dataset(fname)

    lon = ncdata['lon'].values
    lat = ncdata['lat'].values
    sst = ncdata['temp'][startTime:endTime,0,:,:]
    time = ncdata['time'][startTime:endTime].values

    # calculating daily mean
    # mean_sst = sst.groupby(('time.day')).mean('time')
    # time = time[::8]

    # instead of calculate mean of each day, introducing an error
    # because of the nocturnal data, we extract only the sst information
    # from the closest timestep of GHRSST time.
    ind_tmod = np.arange(3,len(time),8)     # vector for indexes of modTime

    sst = sst[ind_tmod,:,:]
    time = time[ind_tmod]

    return lon,lat,time,sst

def find_nearest(lon,lat,ilon,ilat):
    from scipy import spatial

    locations_posi = [[ilat,ilon]]
    locs = np.asarray(locations_posi)

    lo = lon.ravel()
    la = lat.ravel()

    coords = []

    for i,j in zip(la,lo):
        coords.append([i,j])

    coords = np.array(coords)
    tree = spatial.KDTree(coords)

    dists,indexes = tree.query(locs,k=1)
    pontos = []

    for index in indexes:
        pontos.append(coords[index])

    # converter de lista para array
    pontos = np.asarray(pontos)

    idx,idy = np.where((lon == pontos[0][1]) & (lat == pontos[0][0]))
    idx = np.squeeze(idx)
    idy = np.squeeze(idy)

    return idx,idy

##############################################################################
#                   VISUALIZATION FUNCTIONS                                  #
##############################################################################
def plot_animation(i,lon_sat,lat_sat,sst_sat,lon_mod,lat_mod,sst_mod,savefig=None):

    # selecionar os indices dos locais para geracao de serie temporal no modelo
    # indexes for model_grid
    mod_jcan,mod_ican = 19, 70
    mod_jsan,mod_isan = 28, 70
    mod_juba,mod_iuba = 99, 70
    mod_jcab,mod_icab = 126,70

    # buscar os pontos mais proximos na grade do satelite aos pontos do modelo
    # indexes for closest position in satellite grid
    sat_jcan,sat_ican = find_nearest(lon_sat,lat_sat,lon_mod[mod_jcan,mod_ican],lat_mod[mod_jcan,mod_ican])
    sat_jsan,sat_isan = find_nearest(lon_sat,lat_sat,lon_mod[mod_jsan,mod_isan],lat_mod[mod_jsan,mod_isan])
    sat_juba,sat_iuba = find_nearest(lon_sat,lat_sat,lon_mod[mod_juba,mod_iuba],lat_mod[mod_juba,mod_iuba])
    sat_jcab,sat_icab = find_nearest(lon_sat,lat_sat,lon_mod[mod_jcab,mod_icab],lat_mod[mod_jcab,mod_icab])

    contour_levels = np.arange(15.,32.,17./500)

    # plotar
    fig,ax = plt.subplots(ncols=2,figsize=[16/2.54,19/2.54])

    # plot sattelite
    m1 = oceano.make_map(ax[0],resolution='i')
    divider = make_axes_locatable(ax[0])
    cax1 = divider.append_axes("right", size="5%",pad=0.05)
    cax1.set_ylim([15,32])
    divider = make_axes_locatable(ax[0])
    sattemp = divider.append_axes("bottom", size="20%",pad=0.3)
    sattemp.set_xlim([0,31])
    sattemp.set_ylim([15,30])
    # ax[0].set_title('GHRSST - %s')

    # plot modelled
    m2 = oceano.make_map(ax[1],resolution='i')
    divider = make_axes_locatable(ax[1])
    cax2 = divider.append_axes("right", size="5%",pad=0.05)
    cax2.set_ylim([15,32])
    divider = make_axes_locatable(ax[1])
    modtemp = divider.append_axes("bottom", size="20%",pad=0.3)
    modtemp.set_xlim([0,31])
    modtemp.set_ylim([15,30])

    # plot satellite data
    # for i in np.arange(0,32,1):
        # cleaning screen
        # ax[0].clear()
        # ax[1].clear()
        # sattemp.clear()
        # modtemp.clear()
        # sattemp.set_xlim([0,30])
        # sattemp.set_ylim([15,30])
        # modtemp.set_xlim([0,30])
        # modtemp.set_ylim([15,30])
        # m1 = oceano.make_map(ax[0],resolution='i')
        # m2 = oceano.make_map(ax[1],resolution='i')

    # plot sattelite data
    cb1 = m1.contourf(lon_sat,lat_sat,sst_sat[i,:,:],contour_levels,latlon=True,cmap=cmo.cm.thermal)
    plt.colorbar(cb1,cax=cax1)

    # # cax.title('Temporal evolution of bottom temperature in SBC')
    sattemp.plot(sst_sat[:i,sat_jcan,sat_ican],'k',label='Cananeia')
    sattemp.plot(sst_sat[:i,sat_jsan,sat_isan],'r',label='Santos')
    sattemp.plot(sst_sat[:i,sat_juba,sat_iuba],'b',label='Ubatuba')
    sattemp.plot(sst_sat[:i,sat_jcab,sat_icab],'y',label='Cabo Frio')

    sattemp.legend(bbox_to_anchor=(0.265,6.2))
    # cax.fill_between(np.arange(0,len(tmp)), 0, tmp,color='k',alpha=0.4)

    # plot modelling data
    cb2 = m2.contourf(lon_mod,lat_mod,sst_mod[i,:,:],contour_levels,latlon=True,cmap=cmo.cm.thermal)
    plt.colorbar(cb2,cax=cax2)

    modtemp.plot(sst_mod[:i,mod_jcan,mod_ican],'k',label='Cananeia')
    modtemp.plot(sst_mod[:i,mod_jsan,mod_isan],'r',label='Santos')
    modtemp.plot(sst_mod[:i,mod_juba,mod_iuba],'b',label='Ubatuba')
    modtemp.plot(sst_mod[:i,mod_jcab,mod_icab],'y',label='Cabo Frio')

    plt.show()

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

BASE_DIR = oceano.make_dir()
if BASE_DIR.split("/")[2] == 'tparente':
    DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output_modelo/exp03_variables/')

else:
    DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
    GHRSST_DIR = BASE_DIR.replace('github/','ventopcse/data/GHRSST/')

FIGU_DIR = BASE_DIR + 'masterThesis_analysis/figures/experiments_outputs/ghrsst_v_ecom/'

#### READ FILES FROM GHRSST ####
nfiles = glob.glob(GHRSST_DIR+"*.nc")        # read all files from GHRSST
nfiles.sort()                                # sort by name

#### READ FILE FROM ECOM ####
fname = glob.glob(DATA_DIR+run+'.cdf')[0]    # always catching the RUN file
experiment = run                             # storing the experiment name

# extract data from files
lon_sat,lat_sat,time_sat,sst_sat = create_matriz_ghrsst(nfiles)

# here we exclude the last value of GHRSST
time_sat = time_sat[:-1]
sst_sat  = sst_sat[:-1,:,:]

# extract data from files
lon_mod,lat_mod,time_mod,sst_mod = load_ecom(fname)

# basic treatment to remove 0. values from coordinates
lon_sat[lon_sat == 0] = np.nan
lat_sat[lat_sat == 0] = np.nan
lon_mod[lon_mod == 0] = np.nan
lat_mod[lat_mod == 0] = np.nan

########### PLOTTING last timestep only to visualize how the figure will be

plot_animation(-1,lon_sat,lat_sat,sst_sat,lon_mod,lat_mod,sst_mod)

################### testing temperature calibration
plt.ion()
fig,ax = plt.subplots(ncols=2)

sst_mod_alterado = sst + 4.

ax[0].set_ylim([15,32])
ax[1].set_ylim([15,32])

# # cax.title('Temporal evolution of bottom temperature in SBC')
ax[0].plot(sst_sat[:i,sat_jcan,sat_ican],'k',label='Cananeia')
ax[0].plot(sst_sat[:i,sat_jsan,sat_isan],'r',label='Santos')
ax[0].plot(sst_sat[:i,sat_juba,sat_iuba],'b',label='Ubatuba')
ax[0].plot(sst_sat[:i,sat_jcab,sat_icab],'y',label='Cabo Frio')

ax[1].plot(sst_mod_alterado[:i,mod_jcan,mod_ican],'k',label='Cananeia')
ax[1].plot(sst_mod_alterado[:i,mod_jsan,mod_isan],'r',label='Santos')
ax[1].plot(sst_mod_alterado[:i,mod_juba,mod_iuba],'b',label='Ubatuba')
ax[1].plot(sst_mod_alterado[:i,mod_jcab,mod_icab],'y',label='Cabo Frio')
