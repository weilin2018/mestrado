# BULK STRATIFICATION

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
import seawater as sw


import matplotlib
# matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def bulk_stratification(ncin,ilat=[82,86,90,94,98,101],ilon=70,timestep=None):

    # importing variables
    Tsurf = ncin.temp[timestep,0,ilat,:ilon]
    Ssurf = ncin.salt[timestep,0,ilat,:ilon]
    Tbott =ncin.temp[timestep,-1,ilat,:ilon]
    Sbott =ncin.temp[timestep,-1,ilat,:ilon]

    # calculating density based on seawater.dens0
    Dsurf = sw.dens0(Ssurf.values,Tsurf.values)
    Dbott = sw.dens0(Sbott.values,Tbott.values)

    # calculating bulk stratification by subtracting the surface water density from near-bottom
    B = Dbott - Dsurf

    return B

def calcDistance(ncin,ind,lenght):

    lon,lat = ncin['lon'].values, ncin['lat'].values
    lon[lon == 0.] = np.nan
    lat[lat == 0.] = np.nan
    depth = ncin['depth'].values
    sigma = ncin['sigma'].values

    xx = lon.shape[1]

    x,prof,sig = oceano.create_newDepth(lon,depth,sigma,ind)

    inds = np.where(np.isnan(x[0,:])) # aonde e nan
    lats = np.ones([xx])*11

    # removendo nan's
    x =  np.delete(x[0,:],inds)
    lats=np.delete(lats,inds)

    dist2 = np.cumsum(np.append(0,sw.dist(lats,x)[0]))
    if lenght != 0:
        dist = np.tile(dist2,(lenght,1))
    else:
        dist = dist2

    return dist,inds

def plot_nearBottomTemp(ncin,exp,SAVE_FIG,location,ilat=[99]):

    ###### near-bottom temperature
    fig,ax = plt.subplots(nrows=2,figsize=(17./2.54,12./2.54))

    plt.suptitle(u'Temperatura próximo do fundo ao largo de Ubatuba',fontsize=10)

    ax[0].set_title('%s - 14/01/2014'%(exp.replace('.cdf','')),fontsize=8)
    ax[0].set_ylabel(r'Temperatura [$^o$C]',fontsize=8)
    # ax[0].set_xlabel('Distance [km]',fontsize=8)
    ax[0].set_ylim([0,30])
    ax[0].set_xlim([0,140])

    ax[0].tick_params(axis='x',labelbottom='off') # remove ticklabels

    ax[1].set_title('%s - 15/02/2014'%(exp.replace('.cdf','')),fontsize=8)
    ax[1].set_ylabel(r'Temperature [$^o$C]',fontsize=8)
    ax[1].set_xlabel(u'Distância da costa [km]',fontsize=8)
    ax[1].set_ylim([0,30])
    ax[1].set_xlim([0,140])

    # inserting lines inside the graphic
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax[0].axvline(x=location[0],color='k',ls='--')
    ax[0].text(location[0]+.2, 3., '%s'%(location[0]),fontsize=10,verticalalignment='top', bbox=props)
    ax[0].axhline(y=18,color='k',ls='--')
    ax[1].axvline(x=location[1],color='k',ls='--')
    ax[1].text(location[1]+.2, 3., '%s'%(location[1]),fontsize=10,verticalalignment='top', bbox=props)
    ax[1].axhline(y=18,color='k',ls='--')

    for i in ilat:
        t = ncin.temp[46,-1,i,:].values
        dist,inds = calcDistance(ncin,i,0)

        t = np.delete(t,inds,axis=0)
        ax[0].plot(dist,t)
        # ax[0].scatter(dist,t,s=2,alpha=.4,c='k',label='%s'%(i))

    for i in ilat:
        t = ncin.temp[303,-1,i,:].values
        dist,inds = calcDistance(ncin,i,0)

        t = np.delete(t,inds,axis=0)
        ax[1].plot(dist,t)
        # ax[1].scatter(dist,t,s=2,alpha=.4,c='k',label='%s'%(i))

    plt.tight_layout()
    plt.subplots_adjust(top=0.894,bottom=0.1,left=0.093,right=0.958,hspace=0.16,wspace=0.2)
    plt.savefig(SAVE_FIG+"posicao_%s_tempBottom.pdf"%(exp.replace('.cdf','')))

def plot_surfaceSalt(ncin,exp,SAVE_FIG,location,ilat=[99]):

    ###### near-bottom temperature
    fig,ax = plt.subplots(nrows=2,figsize=(17./2.54,12./2.54))

    plt.suptitle(u'Salinidade na superfície ao largo de Ubatuba',fontsize=10)

    # first axes configuration
    ax[0].set_title('%s - 14/01/2014'%(exp.replace('.cdf','')),fontsize=8)
    ax[0].set_ylabel('Salinidade',fontsize=8)
    # ax[0].set_xlabel('Distance [km]',fontsize=8)
    ax[0].set_ylim([33,38])
    ax[0].set_xlim([0,150])

    ax[0].tick_params(axis='both', which='major', labelsize=8) # fontsize ticklabels
    ax[0].tick_params(axis='x',labelbottom='off') # remove ticklabels

    # second axis configuration
    # fig,ax = plt.subplots()
    ax[1].set_title('%s - 15/02/2014'%(exp.replace('.cdf','')),fontsize=8)
    ax[1].set_ylabel('Salinidade',fontsize=8)
    ax[1].set_xlabel(u'Distância da costa [km]',fontsize=8)
    ax[1].set_ylim([33,38])
    ax[1].set_xlim([0,150])

    ax[1].tick_params(axis='both', which='major', labelsize=8) # fontsize ticklabels

    # inserting lines inside the graphic
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax[0].axvline(x=location[0],color='k',ls='--')
    ax[0].text(location[0], 33.8, '%s'%(location[0]),fontsize=10,verticalalignment='top', bbox=props)
    ax[0].axhline(y=36.,color='k',ls='--')
    ax[1].axvline(x=location[1],color='k',ls='--')
    ax[1].text(location[1], 33.8, '%s'%(location[1]),fontsize=10,verticalalignment='top', bbox=props)
    ax[1].axhline(y=36.,color='k',ls='--')

    # plotting data
    for i in ilat:
        t = ncin.salt[46,0,i,:].values
        dist,inds = calcDistance(ncin,i,0)

        t = np.delete(t,inds,axis=0)
        ax[0].plot(dist,t)
        # ax[0].scatter(dist,t,s=2,alpha=.4,c='k',label='%s'%(i))


    for i in ilat:
        t = ncin.salt[303,0,i,:].values
        dist,inds = calcDistance(ncin,i,0)

        t = np.delete(t,inds,axis=0)
        ax[1].plot(dist,t)
        # ax[1].scatter(dist,t,s=2,alpha=.4,c='k',label='%s'%(i))

    plt.tight_layout()
    plt.subplots_adjust(top=0.894,bottom=0.1,left=0.093,right=0.958,hspace=0.16,wspace=0.2)

    # ax[0].set_rasterized(True)
    # ax[1].set_rasterized(True)
    plt.savefig(SAVE_FIG+"posicao_%s_saltSurf.pdf"%(exp.replace('.cdf','')))

def isoterma(f1,f2,SAVE_FIG):

    exp1 = f1.split('/')[-1].replace('.cdf','')
    exp2 = f2.split('/')[-1].replace('.cdf','')

    fig,ax = plt.subplots(nrows=2,ncols=2,figsize=(17./2.54,12./2.54))
    fig.suptitle(u'Seção Vertical de Temperatura ao Largo de Ubatuba',fontsize=10)

    ax[0,0].set_title('%s - 14/01/2014'%(exp1),fontsize=8)
    ax[0,1].set_title('%s - 14/01/2014'%(exp2),fontsize=8)

    ax[1,0].set_title('%s - 15/02/2014'%(exp1),fontsize=8)
    ax[1,1].set_title('%s - 15/02/2014'%(exp2),fontsize=8)

    ax[0,0].set_ylabel('Profundidade [m]',fontsize=8)
    ax[1,0].set_ylabel('Profundidade [m]',fontsize=8)
    ax[1,0].set_xlabel(u'Distância da Costa [km]',fontsize=8)
    ax[1,1].set_xlabel(u'Distância da Costa [km]',fontsize=8)

    ax[0,0].tick_params(axis='both', which='major', labelsize=8) # fontsize ticklabels
    ax[1,0].tick_params(axis='both', which='major', labelsize=8) # fontsize ticklabels
    ax[0,1].tick_params(axis='both', which='major', labelsize=8) # fontsize ticklabels
    ax[1,1].tick_params(axis='both', which='major', labelsize=8) # fontsize ticklabels
    ax[0,0].tick_params(axis='x',labelbottom='off') # remove ticklabels
    ax[0,1].tick_params(axis='x',labelbottom='off') # remove ticklabels

    ncin = xr.open_dataset(f1)

    tempBegin = ncin.temp[46,:,99,:]
    tempFinal = ncin.temp[303,:,99,:]
    sigma = ncin.sigma.values
    depth = ncin.depth.values

    # create new arrays
    x,prof,sig = oceano.create_newDepth(ncin.lon.values,depth,sigma,99)
    dist,inds = calcDistance(ncin,99,21)

    # removendo nan values
    tBegin = np.delete(tempBegin,inds,axis=1)
    tFinal = np.delete(tempFinal,inds,axis=1)
    d      = np.delete(depth,inds,axis=1)
    s      = np.delete(sig,inds,axis=1)

    ax[0,0].fill_between(dist[0,:],-3000,-d[99,:],color='#c0c0c0')
    ax[0,0].set_xlim([0,150])
    ax[0,0].set_ylim([-400,0])
    ax[1,0].fill_between(dist[0,:],-3000,-d[99,:],color='#c0c0c0')
    ax[1,0].set_xlim([0,150])
    ax[1,0].set_ylim([-400,0])

    ax[0,1].fill_between(dist[0,:],-3000,-d[99,:],color='#c0c0c0')
    ax[0,1].set_xlim([0,150])
    ax[0,1].set_ylim([-400,0])
    ax[1,1].fill_between(dist[0,:],-3000,-d[99,:],color='#c0c0c0')
    ax[1,1].set_xlim([0,150])
    ax[1,1].set_ylim([-400,0])

    # plotar isotermas
    cs = ax[0,0].contour(dist,s,tBegin,levels=np.arange(12,30,2),colors=('k'),linestyles=('--'),linewidths=.3)
    cs = ax[0,0].contour(dist,s,tBegin,levels=[18],colors=('k'))
    ax[0,0].clabel(cs, inline=1, fontsize=8)


    # plotar isotermas
    cs = ax[1,0].contour(dist,s,tFinal,levels=np.arange(12,30,2),colors=('k'),linestyles=('--'),linewidths=.3)
    cs = ax[1,0].contour(dist,s,tFinal,levels=[18],colors=('k'))
    ax[1,0].clabel(cs, inline=1, fontsize=8)

    ###

    ncin = xr.open_dataset(f2)

    tempBegin = ncin.temp[46,:,99,:]
    tempFinal = ncin.temp[303,:,99,:]

    # removendo nan values
    tBegin = np.delete(tempBegin,inds,axis=1)
    tFinal = np.delete(tempFinal,inds,axis=1)
    d      = np.delete(depth,inds,axis=1)
    s      = np.delete(sig,inds,axis=1)

    # plotar isotermas
    cs = ax[0,1].contour(dist,s,tBegin,levels=np.arange(12,30,2),colors=('k'),linestyles=('--'),linewidths=.3)
    cs = ax[0,1].contour(dist,s,tBegin,levels=[18],colors=('k'))
    # ax[0,1].clabel(cs, inline=1, fontsize=8)


    # plotar isotermas
    cs = ax[1,1].contour(dist,s,tFinal,levels=np.arange(12,30,2),colors=('k'),linestyles=('--'),linewidths=.3)
    cs = ax[1,1].contour(dist,s,tFinal,levels=[18],colors=('k'))
    # ax[1,1].clabel(cs, inline=1, fontsize=8)

    plt.tight_layout()
    plt.subplots_adjust(top=0.894,bottom=0.082,left=0.106,right=0.962,hspace=0.167,wspace=0.169)

    comparacao = '%sv%s'%(f1.split('/')[-1].replace('.cdf',''),f2.split('/')[-1].replace('.cdf',''))
    plt.savefig(SAVE_FIG+'secao_%s_temp.pdf'%(comparacao))

def isohalina(f1,f2,SAVE_FIG):

    exp1 = f1.split('/')[-1].replace('.cdf','')
    exp2 = f2.split('/')[-1].replace('.cdf','')

    fig,ax = plt.subplots(nrows=2,ncols=2,figsize=(17./2.54,12./2.54))
    fig.suptitle(u'Seção Vertical de Salinidade ao Largo de Ubatuba',fontsize=10)

    ax[0,0].set_title('%s - 14/01/2014'%(exp1),fontsize=8)
    ax[0,1].set_title('%s - 14/01/2014'%(exp2),fontsize=8)

    ax[1,0].set_title('%s - 15/02/2014'%(exp1),fontsize=8)
    ax[1,1].set_title('%s - 15/02/2014'%(exp2),fontsize=8)

    ax[0,0].set_ylabel('Profundidade [m]',fontsize=8)
    ax[1,0].set_ylabel('Profundidade [m]',fontsize=8)
    ax[1,0].set_xlabel(u'Distância da Costa [km]',fontsize=8)
    ax[1,1].set_xlabel(u'Distância da Costa [km]',fontsize=8)

    ax[0,0].tick_params(axis='both', which='major', labelsize=8) # fontsize ticklabels
    ax[1,0].tick_params(axis='both', which='major', labelsize=8) # fontsize ticklabels
    ax[0,1].tick_params(axis='both', which='major', labelsize=8) # fontsize ticklabels
    ax[1,1].tick_params(axis='both', which='major', labelsize=8) # fontsize ticklabels
    ax[0,0].tick_params(axis='x',labelbottom='off') # remove ticklabels
    ax[0,1].tick_params(axis='x',labelbottom='off') # remove ticklabels

    ncin = xr.open_dataset(f1)

    saltBegin = ncin.salt[46,:,99,:]
    saltFinal = ncin.salt[303,:,99,:]
    sigma = ncin.sigma.values
    depth = ncin.depth.values

    # create new arrays
    x,prof,sig = oceano.create_newDepth(ncin.lon.values,depth,sigma,99)
    dist,inds = calcDistance(ncin,99,21)

    # removendo nan values
    tBegin = np.delete(saltBegin,inds,axis=1)
    tFinal = np.delete(saltFinal,inds,axis=1)
    d      = np.delete(depth,inds,axis=1)
    s      = np.delete(sig,inds,axis=1)

    ax[0,0].fill_between(dist[0,:],-3000,-d[99,:],color='#c0c0c0')
    ax[0,0].set_xlim([0,150])
    ax[0,0].set_ylim([-400,0])
    ax[1,0].fill_between(dist[0,:],-3000,-d[99,:],color='#c0c0c0')
    ax[1,0].set_xlim([0,150])
    ax[1,0].set_ylim([-400,0])

    ax[0,1].fill_between(dist[0,:],-3000,-d[99,:],color='#c0c0c0')
    ax[0,1].set_xlim([0,150])
    ax[0,1].set_ylim([-400,0])
    ax[1,1].fill_between(dist[0,:],-3000,-d[99,:],color='#c0c0c0')
    ax[1,1].set_xlim([0,150])
    ax[1,1].set_ylim([-400,0])

    # plotar isotermas
    cs = ax[0,0].contour(dist,s,tBegin,levels=np.arange(35.0,36.8,.2),colors=('k'),linestyles=('--'),linewidths=.3)
    cs = ax[0,0].contour(dist,s,tBegin,levels=[36],colors=('k'))
    ax[0,0].clabel(cs, inline=1, fontsize=8)


    # plotar isotermas
    cs = ax[1,0].contour(dist,s,tFinal,levels=np.arange(35.0,36.8,.2),colors=('k'),linestyles=('--'),linewidths=.3)
    cs = ax[1,0].contour(dist,s,tFinal,levels=[36],colors=('k'))
    ax[1,0].clabel(cs, inline=1, fontsize=8)

    ###

    ncin = xr.open_dataset(f2)

    saltBegin = ncin.salt[46,:,99,:]
    saltFinal = ncin.salt[303,:,99,:]

    # removendo nan values
    tBegin = np.delete(saltBegin,inds,axis=1)
    tFinal = np.delete(saltFinal,inds,axis=1)
    d      = np.delete(depth,inds,axis=1)
    s      = np.delete(sig,inds,axis=1)

    # plotar isotermas
    cs = ax[0,1].contour(dist,s,tBegin,levels=np.arange(35.0,36.8,.2),colors=('k'),linestyles=('--'),linewidths=.3)
    cs = ax[0,1].contour(dist,s,tBegin,levels=[36],colors=('k'))
    ax[0,1].clabel(cs, inline=1, fontsize=8)


    # plotar isotermas
    cs = ax[1,1].contour(dist,s,tFinal,levels=np.arange(35.0,36.8,.2),colors=('k'),linestyles=('--'),linewidths=.3)
    cs = ax[1,1].contour(dist,s,tFinal,levels=[36],colors=('k'))
    ax[1,1].clabel(cs, inline=1, fontsize=8)

    plt.tight_layout()
    plt.subplots_adjust(top=0.894,bottom=0.082,left=0.106,right=0.962,hspace=0.167,wspace=0.169)

    comparacao = '%sv%s'%(f1.split('/')[-1].replace('.cdf',''),f2.split('/')[-1].replace('.cdf',''))
    plt.savefig(SAVE_FIG+'secao_%s_salt.pdf'%(comparacao))
##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
SAVE_FIG = BASE_DIR + 'masterThesis_analysis/figures/experiments_outputs/castro2014/'
fname = glob.glob(DATA_DIR+"*.cdf")

# select which experiment you want to plot:
exp = 'EA1.cdf'

for f in fname:
    if exp in f:
        experiment = f

# plotando near bottom temperature

ncin = xr.open_dataset(experiment)
location = [7.6, 23.85]
plot_nearBottomTemp(ncin,exp,SAVE_FIG,location)

ncin = xr.open_dataset(experiment.replace('EA1','EA5'))
location = [7.6, 19.22]
plot_nearBottomTemp(ncin,exp.replace('EA1','EA5'),SAVE_FIG,location)


# plotando surface salinity
ncin = xr.open_dataset(experiment)
location = [68.12, 139.03]
plot_surfaceSalt(ncin,exp,SAVE_FIG,location)

ncin = xr.open_dataset(experiment.replace('EA1','EA5'))
location = [68.12, 107.73]
plot_surfaceSalt(ncin,exp.replace('EA1','EA5'),SAVE_FIG,location)

# plotando contour secao vertical para comparar com Figura 5 de Castro (2014)
isoterma(experiment,experiment.replace('EA1','EA5'),SAVE_FIG)
isohalina(experiment,experiment.replace('EA1','EA5'),SAVE_FIG)

# plotando contour secao vertical comparando controle e anomalo para o mestrado
isoterma(experiment.replace('EA1','EC1'),experiment,SAVE_FIG)
isohalina(experiment.replace('EA1','EC1'),experiment,SAVE_FIG)
