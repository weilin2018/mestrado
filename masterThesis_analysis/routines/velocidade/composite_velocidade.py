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
import matplotlib.gridspec as gridspec

import matplotlib
# matplotlib.style.use('ggplot')
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

def export_data(fname,timestep=0):
    # plotting climatologic data: t = 0, k = 0
    ncin = xr.open_dataset(fname)

    lon,lat = ncin.lon.values.copy(), ncin.lat.values.copy()
    angle = ncin.ang.values.copy()
    depth = ncin.depth.values.copy()
    sigma = ncin.sigma.values.copy()
    lon[lon == 0.] = np.nan
    lat[lat == 0.] = np.nan

    # extracting temperature data, in a specific timestep
    u,v = ncin.u[timestep,:,:,:],ncin.v[timestep,:,:,:]
    # spd = np.sqrt(u**2 + v**2)
    # spd = np.where(depth < 100, spd,np.nan)

    return lon,lat,u,v,depth,angle

def tratando_corrente(u,v,depth,angle):

    ur,vr = ocplt.rotate_velocityField(u,v,angle)
    spd = np.sqrt(ur**2+vr**2)
    spd = np.where(depth < 100, spd,np.nan)

    return ur,vr,spd

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
FILE_DIR = BASE_DIR+'masterThesis_analysis/routines/index_list.npy'
experiment = 'EC2'
fname = DATA_DIR + experiment +'.cdf'

timestep = [46,303]

for nstep in timestep:
    plt.close()
    # working with the data
    lon,lat,u,v,depth,angles = export_data(fname,timestep=nstep)

    # rotating vectors
    ur_surf,vr_surf,spd_surf = tratando_corrente(u[0,:,:],v[0,:,:],depth,angles)
    ur_meio,vr_meio,spd_meio = tratando_corrente(u[10,:,:],v[10,:,:],depth,angles)
    ur_fund,vr_fund,spd_fund = tratando_corrente(u[19,:,:],v[19,:,:],depth,angles)

    # formatando os dados para um formato de visualizacao melhor e passando os vetores
    # normalizados pela velocidade
    xplot,yplot,uplot_surf,vplot_surf = ocplt.formatting_vectors(ur_surf/spd_surf,vr_surf/spd_surf,lon,lat,FILE_DIR)
    xplot,yplot,uplot_meio,vplot_meio = ocplt.formatting_vectors(ur_meio/spd_meio,vr_meio/spd_meio,lon,lat,FILE_DIR)
    xplot,yplot,uplot_fund,vplot_fund = ocplt.formatting_vectors(ur_fund/spd_fund,vr_fund/spd_fund,lon,lat,FILE_DIR)

    # plotando no grafico
    fig = plt.figure(figsize=(12/2.54,12/2.54))
    gs = gridspec.GridSpec(3,3)

    # creating axis
    ax1 = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[1,1])
    ax3 = plt.subplot(gs[2,2])

    # creating basemap instance
    m1,meridians,parallels = make_map(ax1,ulon=np.nanmax(lon)-.5,llon=np.nanmin(lon)-.2,ulat=np.nanmax(lat)+.2,llat=np.nanmin(lat))
    m2,_,_ = make_map(ax2,ulon=np.nanmax(lon)-.5,llon=np.nanmin(lon)-.2,ulat=np.nanmax(lat)+.2,llat=np.nanmin(lat))
    m3,_,_ = make_map(ax3,ulon=np.nanmax(lon)-.5,llon=np.nanmin(lon)-.2,ulat=np.nanmax(lat)+.2,llat=np.nanmin(lat))

    # positiong each axes
    ax1.set_position([.03,.46,.6,.5])
    ax2.set_position([.1,.28,.6,.5])
    ax3.set_position([.17,.1,.6,.5])

    contours = np.arange(0,1.5,0.1)

    # plotando velocidade
    cf1 = m1.contourf(lon,lat,spd_surf,contours,cmap=cmo.cm.speed,latlon=True,rasterized=True,extend='max')
    cf2 = m2.contourf(lon,lat,spd_meio,contours,cmap=cmo.cm.speed,latlon=True,rasterized=True,extend='max')
    cf3 = m3.contourf(lon,lat,spd_fund,contours,cmap=cmo.cm.speed,latlon=True,rasterized=True,extend='max')
    # plotando vetores
    qv1 = m1.quiver(xplot,yplot,uplot_surf,vplot_surf,scale=60,width=0.0015,headwidth=4,headlength=4,latlon=True)
    qv2 = m2.quiver(xplot,yplot,uplot_meio,vplot_meio,scale=60,width=0.0015,headwidth=4,headlength=4,latlon=True)
    qv3 = m3.quiver(xplot,yplot,uplot_fund,vplot_fund,scale=60,width=0.0015,headwidth=4,headlength=4,latlon=True)


    # matplotib trick to remove white thin lines when saving contourf in pdf
    for c in cf1.collections:
        c.set_edgecolor('face')
        c.set_linewidth(0.00000000001)

    for c in cf2.collections:
        c.set_edgecolor('face')
        c.set_linewidth(0.00000000001)

    for c in cf3.collections:
        c.set_edgecolor('face')
        c.set_linewidth(0.00000000001)

    # colorbar
    cax = fig.add_axes([.1,.08,.35,.02])
    cbar = plt.colorbar(cf1,orientation='horizontal',cax=cax,format='%0.2f')

    # # figure's title
    # plt.suptitle(u'Temperatura nas camadas de superfície, meio e fundo no Experimento Controle (esquerda) e Anômalo (direita) em 14 de Janeiro')

    # setting colorbar tick labels
    from matplotlib import ticker
    tick_locator = ticker.MaxNLocator(nbins=6)
    cbar.locator = tick_locator
    cbar.update_ticks()

    cbar.ax.axes.tick_params(axis='both',which='both',labelsize=8)
    cbar.ax.set_title(r'Velocidade (m s$^{-1}$)',fontsize=8)

    output_fname = fname.split('/')[-1].replace('.cdf','_'+str(nstep))
    plt.savefig('/home/danilo/Dropbox/mestrado/figuras/composicao/speed/%s/%s.eps'%(experiment,output_fname))
