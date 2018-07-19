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

from eofs.standard import Eof
import eofs.multivariate.standard as eofmulti

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

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

def savePickle(filename,dct):
    pickle.dump(dct,open(filename,'w'))

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

nt              = len(nfiles)   # quantos anos

id              = 0

# processa cada arquivo
for i in np.arange(0,124):
    # criar matriz com nt,nlat,nlon
    data = createMatrix(nlon,nlat,nt)
    id = 0
    # processar cada fname, armazenando os dados e calculando a EOF
    for fname in nfiles:
        ncin = xr.open_dataset(fname)
        d    = ncin.variables['U_GRD_L103'][i,:,:]
        data[id,:,:] = d
        id += 1
        ncin.close()

    # com todos os arquivos processados e armazenados em uma matriz,
    # calcular EOF
    eofs,pcs = calculateEof(data)

    # armazenar os dados em um pickle com numero do tempo
    outFile = DATA_DIR.replace('data/','jan_%i.pickle'%(i))
    savePickle(outFile,{'eofs':eofs,'pcs':pcs})
