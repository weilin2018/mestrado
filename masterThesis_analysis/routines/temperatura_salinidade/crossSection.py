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

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

import scipy.ndimage

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here

def calcDistance(x,sig,dx,xx=110):

    inds = np.where(np.isnan(x[0,:])) # aonde e nan
    lats = np.ones([xx])*11

    # removendo nan's
    x =  np.delete(x[0,:],inds)
    lats=np.delete(lats,inds)

    # dist2 = np.cumsum(np.append(0,sw.dist(lats,x)[0]))
    # dist = np.tile(dist2,(len(sig),1))

    dist      = np.cumsum(dx)/1000
    dist      = np.tile(dist,(len(sig),1))
    dist      = np.delete(dist,inds,axis=1)

    return dist,inds

def crossSection(fname):

    ind = 99 # latitude for transect
    ncin = xr.open_dataset(fname)

    # extract grid and other general variables
    # important: lat and lon already gridded
    lon   = ncin.lon.values
    lat   = ncin.lat.values
    lon[lon == 0.] = np.nan
    lat[lat == 0.] = np.nan
    depth = ncin.depth.values
    sigma = ncin.sigma.values
    h1    = ncin['h1'].values
    dx = h1[ind,:]
    temp = ncin.temp[300:303,:,:,:]

    # cutting variable data, to speed up interpolation
    cutLon = 80
    Tcutted = temp[:,:,:,:cutLon]

    # computing distance based on h1
    dist = np.cumsum(dx)/1000

    # bathymetry
    x,prof,sig = oceano.create_newDepth(lon,depth,sigma,ind)
    dist2 = np.tile(dist,(21,1))

    # interpolating sigma to standard level
    nstdl = 100
    stdl,Tinterp = oceano.sigma2stdl(Tcutted,sigma,nstdl,depth,h1,lon,lat,'Temperatura')

    # creating vertical matrix, for visualization
    grid_x,grid_z = np.meshgrid(dist[:cutLon],np.linspace(0,100,100))

    var2plot = np.nanmean(Tinterp[:,:,ind,:],axis=0)

    fig,ax = plt.subplots()
    ax.contourf(grid_x,-grid_z,var2plot)
    ax.plot(grid_x,-grid_z,'k',alpha=.3);
    ax.plot(grid_x.T,-grid_z.T,'k',alpha=.3);

    ax.plot(dist2[-1,:],sig[-1,:],'k')

    ax.set_xlim([0,100])
    ax.set_ylim([-100,0])

    return grid_x,grid_z,Tinterp,stdl,dist,dist2,sig

def sigma2stdl_teste(variable,sigma,nlines,depth,h1,lon,lat,name):

    os.system('clear')
    print('\ninterpolating sigma to standard levels: ' + name)

    # creating 4D matrix to store interpolated data
    nsteps,nstdl,lines,columns = variable.shape
    # nsteps, sigmalevels, lines, columns = variable.shape
    variableI = np.zeros((nsteps,nstdl,nlines,columns)) * np.nan
    # variableI = np.zeros((nsteps, nstdl, lines, columns))
    # variableI[:] = np.NAN

    # creating new resolution on x axis
    if nlines == 1000: # adapted for inner and middle shelf
        # vertical resolution: 1km
        newline = np.linspace(0,100000,1000)

    ### --- vectorize data to decrease number of for loops --- ###

    vecvar = np.zeros((lines,nsteps*nstdl*columns))
    vecvarI= np.zeros((nlines,nsteps*nstdl*columns))

    #############################################################
    # comecar a pensar na interpolacao horizontal, linha a linha#
    #############################################################
    locallon = np.zeros((nsteps*nstdl*columns))
    cont = 0

    for n in np.arange(0, nsteps):
    	# print('nstep ' + str(n) + ' of ' + str(nsteps-1))
    	for z in np.arange(0, nstdl):
    		for c in np.arange(0, columns):
    			vecvar[:, cont] = variable[n,z,:, c]
    			locallon[cont] = lon[l, c]
    			cont += 1

    for i i nnp.arange(0,nsteps*nstdl*columns):
        if ~np.isnan()




    # vecvar = np.zeros((sigmalevels, nsteps*lines*columns))
    # vecvarI = np.zeros((nstdl, nsteps*lines*columns))
    # localdep = np.zeros((nsteps*lines*columns))
    # cont = 0
    #
    # for n in np.arange(0, nsteps):
    # 	# print('nstep ' + str(n) + ' of ' + str(nsteps-1))
    # 	for l in np.arange(0, lines):
    # 		for c in np.arange(0, columns):
    # 			vecvar[:, cont] = variable[n, :, l, c]
    # 			localdep[cont] = depth[l, c]
    # 			cont += 1

    ### --- interpolate sigma to standard levels --- ###
    # plot counter
    kplot = int(nsteps*lines*columns/10)
    k = 0
    for i in np.arange(0, nsteps*lines*columns):
    	k += 1
    	if k == kplot:
    		print(str(np.round(i/(nsteps*lines*columns)*100)) + '%')
    		k = 0

    	if ~np.isnan(localdep[i]):
    		# print('yes 1/2')

    		# levels of data (m)
    		depsigma = -localdep[i]*sigma

    		# include surface with same value of first sigma level m to interpolate
    		D = list(depsigma)
    		D.insert(0, 0)

    		# select profile and include surface
    		profile = np.zeros(sigmalevels+1)
    		profile[1:] = vecvar[:, i]
    		profile[0] = profile[1]

    		# watercolumn positions only
    		watercolumn = stdl <= localdep[i]
    		stdl2interp = np.array(stdl)[watercolumn]

    		# interpolate to the same standard levels
    		fsigma2stdl = interpolate.interp1d(D, profile)
    		profileI = fsigma2stdl(stdl2interp)

    		# stores at vectorized variable
    		vecvarI[watercolumn, i] = profileI
    		vecvarI[~watercolumn, i] = np.NAN

    ### --- back to original shape --- ###

    cont = 0

    for n in np.arange(0, nsteps):
    	for l in np.arange(0, lines):
    		for c in np.arange(0, columns):
    			variableI[n, :, l, c] = vecvarI[:, cont]
    			cont += 1

    # corners to NaN
    variableI[:, :, 1, -2] = np.NAN
    variableI[:, :, -2, -2] = np.NAN

    return stdl,variableI


##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
SAVE_DIR = BASE_DIR + 'masterThesis_analysis/figures/experiments_outputs/temperature/'
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
plt.ion()

# select which experiment you want to plot:
fname = DATA_DIR + 'EA1.cdf'
grid_x,grid_z,Tinterp,stdl,dist,dist2,sig = crossSection(fname)
