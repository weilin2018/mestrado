#!  /usr/bin/env python2
#

"""
Compute and plot the leading EOF of sea surface temperature in the
central Pacific during winter time.

The spatial pattern of this EOF is the canonical El Nino pattern, and
the associated time series shows large peaks and troughs for well-known
El Nino and La Nina events.

This example uses the plain numpy interface.
"""

from matplotlib.backends.backend_pdf import PdfPages
# from netCDF4 import Dataset
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import os
import glob
from scipy import signal
import calendar
from datetime import date, timedelta
import time
import pickle
import math
import numpy.polynomial.polynomial as poly

from eofs.standard import Eof
from eofs.examples import example_data_path

#*** Declare parameters ****************************************************
nlon = 34
nlat = 17
ndy  = 36

startY = 1982
endY   = 2000

ny = endY - startY + 1

nm = 12
nd = 366
nt = ny*ndy

lats = [-19.875+x*2.5 for x in range(0, nlat)]
lons = [-154.875+x*2.5 for x in range(0, nlon)]

lats=np.asarray(lats)
lons=np.asarray(lons)

#*** Declare variables ****************************************************
sst_raw=np.empty((nt,nlat,nlon))
sst_all=np.empty((ny,ndy,nlat,nlon))
sst_diff=np.empty((ny,ndy,nlat,nlon))
sst_final=np.empty((nt,nlat,nlon))

sst_detrend=np.empty((nt,nlat,nlon))
sst_coeffs=np.empty((2,nlat,nlon))
sst_season=np.empty((37,nlat,nlon))

sst_detrend[:,:,:] = np.nan
sst_all[:,:,:,:] = np.nan

id = 0
for i in range(startY,startY+ny):
  for j in range(1,13):
    for filename in glob.glob('/home/danilo/Dropbox/mestrado/data/eof_data/'+str(i)+"{0:0>2}".format(j)+'*.nc'):
       ncin = xr.open_dataset(filename)
       sst = ncin.variables['analysed_sst'][:]
       sst_raw[id,:,:] = sst[:,:]
       #sst_raw[id,:,:] = ncin.variables['analysed_sst'][:]
       id = id + 1
       ncin.close()

#*** Detrend **************************************************************
x=np.empty((nt))
for i in range(0,nt):
  x[i] = i

for i in range(0,nlat):
  for j in range(0,nlon):
    ytemp = np.copy(sst_raw[:,i,j])
    y = sst_raw[:,i,j]
    b = ~np.isnan(y)
    coefs = poly.polyfit(x[b], y[b], 1)
    sst_coeffs[0,i,j] = coefs[0]
    sst_coeffs[1,i,j] = coefs[1]
    ffit = poly.polyval(x[b], coefs)
    sst_detrend[b,i,j] = y[b] - ffit

#*** Rearrange data for seasonal removal *******************************************************
id1 = 0
for i in range(startY,startY+ny):
   id2 = id1 + ndy
   sst_all[i-startY,0:id2-id1,:,:] = sst_detrend[id1:id2,:,:]
   id1 = id2

#*** Calculate seasonal cycle *******************************************************
sst_season = np.mean(sst_all, axis=0)

#*** Remove seasonal cycle *******************************************************
for i in range(0,ny):
   sst_diff[i,:,:,:] = sst_all[i,:,:,:] - sst_season[:,:,:]

#*** Rearrange array for EOF *******************************************************
id1 = 0
for i in range(startY,startY+ny):
   id2 = id1 + ndy
   sst_final[id1:id2,:,:] = sst_diff[i-startY,0:id2-id1,:,:]
   id1 = id2

#*** Create an EOF solver to do the EOF analysis. Square-root of cosine of **********
#*** latitude weights are applied before the computation of EOFs.          **********
coslat = np.cos(np.deg2rad(lats))
wgts = np.sqrt(coslat)[..., np.newaxis]
solver = Eof(sst_final, weights=wgts)

# Retrieve the leading EOF, expressed as the correlation between the leading
# PC time series and the input SST anomalies at each grid point, and the
# leading PC time series itself.

eof1 = solver.eofs(neofs=10)
pc1 = solver.pcs(npcs=10, pcscaling=0)
varfrac = solver.varianceFraction()
lambdas = solver.eigenvalues()

# Plot the leading EOF expressed as correlation in the Pacific domain.

parallels = np.arange(-90,90,10.)
meridians = np.arange(-180,180,20)

for i in range(0,5):
    fig=plt.figure()
    plt.subplot(211)
    #ax=fig.add_axes([0.1,0.1,0.8,0.8])
    m = Basemap(projection='cyl', llcrnrlon=min(lons), llcrnrlat=min(lats),
        urcrnrlon=max(lons), urcrnrlat=max(lats))
    x, y = m(*np.meshgrid(lons, lats))
    clevs = np.linspace(np.min(eof1[i,:,:].squeeze()), np.max(eof1[i,:,:].squeeze()), 11)
    cs=m.contourf(x, y, eof1[i,:,:].squeeze(), clevs, cmap=plt.cm.RdBu_r)
    m.drawcoastlines()
    m.fillcontinents(color='#000000',lake_color='#99ffff')
    m.drawparallels(parallels,labels=[1,0,0,0])
    m.drawmeridians(meridians,labels=[1,0,0,1])

    #cb = plt.colorbar(cs, orientation='horizontal')
    cb = m.colorbar(cs, 'right', size='5%', pad='2%')
    cb.set_label('EOF', fontsize=12)
    plt.title('EOF ' + str(i+1), fontsize=16)

    plt.subplot(212)
    days = [startY+(x*10+1)/365.0 for x in range(0, nt)]
    plt.plot(days, pc1[:,i], linewidth=2)
    plt.xticks(range(startY, endY), rotation='vertical')
    plt.axhline(0, color='k')
    plt.xlabel('Year')
    plt.ylabel('PC Amplitude')
    plt.xlim(startY, endY)
    plt.ylim(np.min(pc1.squeeze()), np.max(pc1.squeeze()))
    plt.xticks(range(startY, endY))
    plt.tight_layout()
    plt.show()

plt.figure()
eof_num = range(1, 16)
plt.plot(eof_num, varfrac[0:15], linewidth=2)
plt.plot(eof_num, varfrac[0:15], linestyle='None', marker="o", color='r', markersize=8)
plt.axhline(0, color='k')
plt.xticks(range(1, 16))
plt.title('Fraction of the total variance represented by each EOF')
plt.xlabel('EOF #')
plt.ylabel('Variance Fraction')
plt.xlim(1, 15)
plt.ylim(np.min(varfrac), np.max(varfrac)+0.01)
plt.show()
