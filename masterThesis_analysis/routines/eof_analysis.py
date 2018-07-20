
#add some description here

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

from eofs.standard import Eof
from eofs.multivariate.standard import MultivariateEof 

import matplotlib
#matplotlib.style.use('ggplot')

#import sys
#sys.path.append('masterThesisPack/')

#import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here

def createMatrix(nlon,nlat,nt):
    matrix = np.empty((nt,nlat,nlon))

    return matrix

def calculateEof(data):
    solver = Eof(data)
    eofs   = solver.eofs(neofs=3)
    pcs    = solver.pcs(npcs=3,pcscaling=0)

    return eofs,pcs

def calculateEofMulti(dt1,dt2):
    solver = MultivariateEof([dt1,dt2])
    eofs   = solver.eofs(neofs=3)
    pcs    = solver.pcs(npcs=3,pcscaling=0)

    return eofs,pcs

def savePickle(filename,dct):
    pickle.dump(dct,open(filename,'w'))

def processMonth(month,ntime,nfiles):
    nt = len(nfiles)
    for i in np.arange(0,ntime):
        print('Processing timestep %i'%(i))
        # create matrix with nlon,nlat,nt
        data = createMatrix(nlon,nlat,nt)
        # variable to control which timestep
        id = 0
        # read each file in nfiles
        for fname in nfiles:
            # load netcdf file
            ncin = xr.open_dataset(fname)
            # extract the variable for the time i
            d    = ncin.variables['U_GRD_L103'][i,:,:]
            # save extracted data into matrix
            data[id,:,:] = d
            id += 1
            ncin.close()

	# performing eof analysis
        eofs,pcs = calculateEof(data)
        # save eofs for each timestemp of input
        outFile = DATA_DIR.replace('data/','output/%s_%i.pickle'%(month,i))
        savePickle(outFile,{'eofs':eofs,'pcs':pcs})

def processMulti(month,ntime,nfiles):
    nt = len(nfiles)
    for i in np.arange(0,ntime):
        print('Processing timestep %i'%(i))
        # create matrix with nlon,nlat,nt
        data_u = createMatrix(nlon,nlat,nt)
        data_v = createMatrix(nlon,nlat,nt)
       # variable to control which timestep
        id = 0
        # read each file in nfiles
        for fname in nfiles:
            # load netcdf file
            ncin = xr.open_dataset(fname)
            # extract the variable for the time i
            wu    = ncin.variables['U_GRD_L103'][i,:,:]
            wv    = ncin.variables['V_GRD_L103'][i,:,:]
            # save extracted data into matrix
            data_u[id,:,:] = wu
            data_v[id,:,:] = wv
            id += 1
            ncin.close()

        # performing eof analysis
        eofs,pcs = calculateEofMulti(data_u,data_v)
        # save eofs for each timestemp of input
        outFile = DATA_DIR.replace('data/','output/%s_%i.pickle'%(month,i))
        savePickle(outFile,{'eofs':eofs,'pcs':pcs})

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

DATA_DIR = '/home/danilos/mestrado/processing_climatology/data/'
nlon     = 65
nlat     = 64


# leitura dos arquivos para Jan
nfiles          = glob.glob(DATA_DIR+"jan/*.nc")
nfiles.sort()

#processMonth(month='jan',ntime=124,nfiles=nfiles)
# perform a multivariate eof
processMulti(month='jan',ntime=124,nfiles=nfiles)

###### perfoming eof analysis for february data
nfiles         = glob.glob(DATA_DIR+"fev/*.nc")
nfiles.sort()

# process each file, for each day, for each 6-hourly input
# important: hence I'll not compare with a bissext year (2014), then I'm not
# interested in data for Feb, 29.
#processMonth(month='feb',ntime=112,nfiles=nfiles)
# perform a multivariate eof
processMulti(month='fev',ntime=112,nfiles=nfiles)

###### perfoming eof analysis for march data
nfiles         = glob.glob(DATA_DIR+"mar/*.nc")
nfiles.sort()

# process each file, for each day, for each 6-hourly input
#processMonth(month='mar',ntime=124,nfiles=nfiles)
# perform a multivariate eof
processMulti(month='mar',ntime=124,nfiles=nfiles)
