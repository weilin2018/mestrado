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

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def make_map(ax,llat=-30,ulat=-20,llon=-50,ulon=-39,resolution='l',nmeridians=3,nparallels=2):

    m = Basemap(projection='merc', llcrnrlat=llat, urcrnrlat=ulat, llcrnrlon=llon, urcrnrlon=ulon, resolution=resolution)

    m.ax = ax

    m.drawcoastlines(linewidth=0.5)
    m.drawmapboundary()
    m.drawstates(linewidth=0.2)
    m.drawcountries(linewidth=0.2)
    # m.drawrivers(linewidth=.1,color='#c0c0c0')
    # m.fillcontinents(color='#c0c0c0')

    # m.drawcoastlines(linewidth=.1)
    m.drawmapboundary()

	# definir meridianos e paralelos para plotar no mapa
    meridians=np.arange(llon,ulon,nmeridians)
    parallels=np.arange(llat,ulat,nparallels)
	# desenhar meridianos e paralelos conforme definido acima
    m.drawparallels(parallels,labels=[True,False,False,True],fontsize=8,fontweight='bold',color='gray',fmt='%2.2f')
    m.drawmeridians(meridians,labels=[True,False,False,True],fontsize=8,fontweight='bold',color='gray',fmt='%2.2f')

    return m

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
FILE_DIR = '/media/danilo/Danilo/RegCM/2nd_Workshop/Project_output/output_precipitation/output/'
SAVE_DIR = '/media/danilo/Danilo/RegCM/2nd_Workshop/Project_output/output_precipitation/figures/'
nfiles = glob.glob(FILE_DIR+"*.nc")
nfiles.sort()

variable = 'SRF'

# select only files with SHF
flist = []
for n in nfiles:
    if variable in n:
        flist.append(n)

flist = np.asarray(flist)

for fname,month in zip(flist,['Nov/1979','Dec/1979','Jan/1980','Feb/1980','Mar/1980']):

    # importing data
    ncin = xr.open_dataset(fname)
    pr = ncin.pr

    # importing grid
    xlon = ncin.xlon
    ylat = ncin.xlat

    # converting from flux [kg m-2 s-1] to accumulated precipitation in mm,
    # by multiplying the monthly mean by the seconds in the month related.
    monthly = pr.mean(dim='time') # mean
    # monthly *= pr.shape[0]*3600 # total of seconds of each month
    monthly *= 86400
    pr.attrs['units'] = u'mm day-1' # also changing the unit

    # plotting a monthly mean
    ulon = xlon.values.max()
    llon = xlon.values.min()
    ulat = ylat.values.max()
    llat = ylat.values.min()

    contours = np.arange(0., 15,15/100.)

    fig,ax = plt.subplots()
    cax = fig.add_axes([0.17,0.12,0.66,0.02])

    m = make_map(ax,ulon=ulon,ulat=ulat,llon=llon,llat=llat,nmeridians=5,nparallels=4)
    # m.plot(xlon.values,ylat.values,'k',latlon=True,alpha=.3);
    # m.plot(xlon.values.T,ylat.values.T,'k',latlon=True,alpha=.3);

    cf = m.contourf(xlon.values,ylat.values,monthly,contours,latlon=True)
    cbar = plt.colorbar(cf,cax=cax,orientation='horizontal')
    cbar.set_label('%s [%s]'%(pr.attrs['long_name'],pr.attrs['units']))

    figureName = fname.split('/')[-1][12:-3]

    plt.suptitle('Monthly Precipitation for %s'%(month),fontsize=10)
    plt.tight_layout()
    plt.subplots_adjust(top=0.941,bottom=0.184,left=0.023,right=0.977,hspace=0.2,wspace=0.2)


    outFile = SAVE_DIR + figureName + '.png'
    plt.savefig(outFile,dpi=200)

    os.system('convert -trim %s %s'%(outFile,outFile))
