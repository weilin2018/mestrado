#-*-coding;utf-8-*-

# arquivo contendo funcoes utilizadas no artigo
import numpy as np
import xray as xr
import pandas as pd
from scipy.spatial import cKDTree
from scipy import signal, fftpack
import scipy
import socket
import matplotlib.pyplot as plt
import glob
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import pickle
import math
import cmocean as cmo
import decomp
# importing root functions
import masterThesisPack

def make_map(ax,llat=-30,ulat=-20,llon=-50,ulon=-39,resolution='l',nmeridians=3,nparallels=2,labels=[True,False,False,True]):

    if resolution == 'f':
        m = pickle.load(open('pickles/basemap.p','r'))
    else:
        m = Basemap(projection='merc', llcrnrlat=llat, urcrnrlat=ulat, llcrnrlon=llon, urcrnrlon=ulon, resolution=resolution)

    m.ax = ax

    m.drawcoastlines(linewidth=.01,color='#c0c0c0')
    m.drawmapboundary()
    m.fillcontinents(color='#c0c0c0')
    m.drawstates(linewidth=.01,color='gray')
	# definir meridianos e paralelos para plotar no mapa
    meridians=np.arange(llon,ulon,nmeridians)
    parallels=np.arange(llat,ulat,nparallels)
	# desenhar meridianos e paralelos conforme definido acima
    m.drawparallels(parallels,labels=labels,fontsize=8,color='gray',linewidth=.02)
    m.drawmeridians(meridians,labels=labels,fontsize=8,color='gray',linewidth=.02)

    return m

def create_Structure_horizontal(fname,contours,property='temp',timestep=0,savefig=False):
    """Funcao para criar uma estrutura horizontal de subplots, com os campos termohalinos da saida do modelo.

    Parameters
    ----------
    fname : str
        Caminho para o arquivo netcdf a ser lido.
    contours : np.ndarray
        Contours adequado para a propriedade a ser plotada.
    property : str
        Propriedade a ser plotada. Deve ser o mesmo padrao do arquivo netcdf
    timestep : int
        Timestep a ser plotado
    savefig : boolean
        Se deseja salvar ou nao a figura.

    Returns
    -------
    type
        Description of returned object.

    """

    sigmaLevels = [0,10,19]

    P = {
      'temp': 'Temperatura',
      'salt': 'Salinidade',
      'speed': 'Velocidade'
    }

    colorbarTitle = {
      'temp': r'Temperatura ($^o$C)',
      'salt': 'Salinidade',
      'speed': r'Velocidade (m s$^-1$)'
    }

    colormap = {
        'temp':cmo.cm.thermal,
        'salt':cmo.cm.haline,
        'speed':cmo.cm.speed
    }


    #fig,axes = plt.subplots(nrows=2,ncols=3,figsize=(16/2.54, 13/2.54))
    fig,axes = plt.subplots(nrows=2,ncols=3,figsize=(11.69,8.27))
    cax = fig.add_axes([0.2,0.05,0.61,0.02])

    # dictionary containing labels for subplots
    labels_dict = {
        '00': [True,False,False,False],
        '01': [False,False,False,False],
        '02': [False,False,False,False],
        '10': [True,False,False,True],
        '11': [False,False,False,True],
        '12': [False,False,False,True],
    }

    m = {}

    for j in range(3):
        for i in range(2):
            key = "%s%s"%(i,j)
            m[key] = make_map(axes[i,j],labels=labels_dict[key],ulat=-21,llat=-29,ulon=-40,resolution='f')
            axes[i,j].spines['left'].set_linewidth(0.2)
            axes[i,j].spines['right'].set_linewidth(0.2)
            axes[i,j].spines['bottom'].set_linewidth(0.2)
            axes[i,j].spines['top'].set_linewidth(0.2)

    # plotting climatologic data: t = 0, k = 0
    ncin = xr.open_dataset(fname)

    lon,lat = ncin.lon.values, ncin.lat.values
    depth = ncin.depth.values
    sigma = ncin.sigma.values
    lon[lon == 0.] = np.nan
    lat[lat == 0.] = np.nan
    depth = ncin.depth.values

    # extracting temperature data, in a specific timestep
    if property == 'speed':
        ncin['speed'] = np.sqrt(ncin.u**2 + ncin.v**2)

    data = ncin[property][timestep,:,:,:]
    data = np.where(depth < 100, data,np.nan)

    # key for axes in m
    col1 = ['00','01','02']

    for key,k in zip(col1,sigmaLevels):
        a = m[key]
        cf = a.contourf(lon,lat,data[k,:,:],contours,latlon=True,cmap=colormap[property])

        if (property == 'salt') and (k == 0):
            cs = a.contour(lon,lat,data[k,:,:],levels=[36.],latlon=True,colors=('black'),linewidths=(0.5))
        if (property == 'temp') and (k == 20):
            cs = a.contour(lon,lat,data[k,:,:],levels=[18.],latlon=True,colors=('black'),linewidths=(0.5))

    # plotting anomalous experiment at the final
    ncin = xr.open_dataset(fname.replace('EC','EA'))
    if property == 'speed':
        ncin['speed'] = np.sqrt(ncin.u**2 + ncin.v**2)

    data = ncin[property][timestep,:,:,:]
    data = np.where(depth < 100, data,np.nan)

    col1 = ['10','11','12']
    for key,k in zip(col1,sigmaLevels):
        a = m[key]
        cf = a.contourf(lon,lat,data[k,:,:],contours,latlon=True,cmap=colormap[property])
        if k == 0:
            cs = a.contour(lon,lat,data[k,:,:],levels=[36.],latlon=True,colors=('black'),linewidths=(0.5))

    axes[0,1].set_title('Experimento Controle',fontsize=8)
    axes[1,1].set_title(u'Experimento Anomalo',fontsize=8)

    # setting colorbar configuration
    cb = plt.colorbar(cf,orientation='horizontal',cax=cax,format='%i')
    fig.text(0.45,0.075,colorbarTitle[property],fontsize=8)

    # title and some figure adjusts
    d = pd.to_datetime(ncin.time[timestep].values)
    plt.suptitle(u'%s nas camadas de superficie, meio e fundo, no Experimento\n' \
                  u'Controle (superior) e Anomalo (inferior) em ' \
                  '%s de %s'%(P[property],d.strftime('%d'),d.strftime('%B')),fontsize=10)
    rect = (0,0.08,1.,0.95)
    plt.tight_layout(rect=rect) # box for tight_subplot_layout
    # plt.subplots_adjust(top=0.886,bottom=0.109,left=0.054,right=0.995,hspace=0.0,wspace=0.045)
    plt.subplots_adjust(top=0.915,bottom=0.11,left=0.036,right=0.999,hspace=0.082,wspace=0.061)

    if savefig:
        savefig_dir = masterThesisPack.make_dir()
        # plt.savefig('/media/danilo/Danilo/mestrado/github/masterThesis_analysis/figures/experiments_outputs/temperature/temperatura_superf_meio_fundo_timestep_%s.png'%(str(timestep)),dpi=300)
        if property == 'temp':
          plt.savefig(savefig_dir+'masterThesis_analysis/figures/experiments_outputs/temperature/temperatura_superf_meio_fundo_timestep_%s.eps'%(str(timestep)),orientation='landscape')
        if property == 'salt':
          plt.savefig(savefig_dir+'masterThesis_analysis/figures/experiments_outputs/salinity/salinidade_superf_meio_fundo_timestep_%s.eps'%(str(timestep)),orientation='landscape')
        if property == 'speed':
          plt.savefig(savefig_dir+'masterThesis_analysis/figures/experiments_outputs/velocity/valocity_superf_meio_fundo_timestep_%s.eps'%(str(timestep)),orientation='landscape')

    return fig,axes

def formatGrid_plot(grid,fname):
    import numpy as np
    ij=np.load(fname)
    # for a 2D array (lon,lat)
    if len(grid.shape)==2:
        grid=grid[ij[1], :]
        grid=grid[:, ij[0]]
    # if grid is a 3D array (temp,salt,speed)
    if len(grid.shape)==3:
        grid=grid[:,ij[1], ij[0]]
    return grid

def formatting_vectors(u,v,lon,lat,auxFile):
    xplot = formatGrid_plot(lon,auxFile)
    yplot = formatGrid_plot(lat,auxFile)
    uplot = formatGrid_plot(u,auxFile)
    vplot = formatGrid_plot(v,auxFile)

    return xplot,yplot,uplot,vplot

def rotate_velocityField(u,v,ang):

    import decomp
    ur = np.zeros(u.shape)*np.nan
    vr = np.zeros(v.shape)*np.nan

    for j in range(u.shape[0]):
        U,V = u[j,:].values,v[j,:].values
        angle = ang[j,:]

        INT,DIR = decomp.uv2intdir(U,V,0,angle)
        uro,vro = decomp.intdir2uv(INT,DIR,0,angle)
        ur[j,:] = uro
        vr[j,:] = vro

    return ur,vr

def create_Structure_horizontal_Quiver(fname,contours,FILE_DIR,property='speed',timestep=0,savefig=False):
    """Funcao para criar uma estrutura horizontal de subplots, com os campos de velocidade
    da saida do modelo.

    Parameters
    ----------
    fname : str
        Caminho para o arquivo netcdf a ser lido.
    contours : np.ndarray
        Contours adequado para a propriedade a ser plotada.
    property : str
        Propriedade a ser plotada. Deve ser o mesmo padrao do arquivo netcdf
    timestep : int
        Timestep a ser plotado
    savefig : boolean
        Se deseja salvar ou nao a figura.

    Returns
    -------
    type
        Description of returned object.

    """

    import decomp as dp

    sigmaLevels = [0,10,19]

    P = {
      'speed': 'Velocidade'
    }

    colorbarTitle = {
      'speed': r'Velocidade (m s$^{-1}$)'
    }

    colormap = {
        'speed':cmo.cm.speed
    }


    #fig,axes = plt.subplots(nrows=2,ncols=3,figsize=(16/2.54, 13/2.54))
    fig,axes = plt.subplots(nrows=2,ncols=3,figsize=(11.69,8.27))
    cax = fig.add_axes([0.2,0.05,0.61,0.02])

    # dictionary containing labels for subplots
    labels_dict = {
        '00': [True,False,False,False],
        '01': [False,False,False,False],
        '02': [False,False,False,False],
        '10': [True,False,False,True],
        '11': [False,False,False,True],
        '12': [False,False,False,True],
    }

    m = {}

    for j in range(3):
        for i in range(2):
            key = "%s%s"%(i,j)
            m[key] = make_map(axes[i,j],labels=labels_dict[key],ulat=-21,llat=-29,ulon=-40,resolution='f')
            axes[i,j].spines['left'].set_linewidth(0.2)
            axes[i,j].spines['right'].set_linewidth(0.2)
            axes[i,j].spines['bottom'].set_linewidth(0.2)
            axes[i,j].spines['top'].set_linewidth(0.2)

    # plotting climatologic data: t = 0, k = 0
    ncin = xr.open_dataset(fname)

    lon,lat = ncin.lon.values.copy(), ncin.lat.values.copy()
    depth = ncin.depth.values.copy()
    sigma = ncin.sigma.values.copy()
    angle = ncin.ang.values.copy()
    lon[lon == 0.] = np.nan
    lat[lat == 0.] = np.nan

    # extracting temperature data, in a specific timestep
    u = ncin['u'][timestep,:,:,:].copy()
    v = ncin['v'][timestep,:,:,:].copy()

    # key for axes in m
    col1 = ['00','01','02']

    for key,k in zip(col1,sigmaLevels):
        # rotate vectors and calculate speed
        ur,vr = rotate_velocityField(u[k,:,:],v[k,:,:],angle)
        speed = np.sqrt(ur**2 + vr**2)
        speed = np.where(depth < 100, speed,np.nan)
        un,vn = ur/speed,vr/speed

        xplot,yplot,uplot,vplot = formatting_vectors(un,vn,lon,lat,FILE_DIR)


        a = m[key]
        cf = a.contourf(lon,lat,speed,contours,latlon=True,cmap=colormap[property])
        qv = a.quiver(xplot,yplot,uplot,vplot,scale=70,width=0.001,headwidth=4,headlength=4,latlon=True)

    # clearing space
    del ncin,u,v,xplot,yplot,uplot,vplot

    # plotting anomalous experiment at the final
    ncin = xr.open_dataset(fname.replace('EC','EA'))
    speed = np.sqrt(ncin.u[timestep,:,:,:]**2 + ncin.v[timestep,:,:,:]**2)

    u = ncin['u'][timestep,:,:,:]
    v = ncin['v'][timestep,:,:,:]
    speed = np.where(depth < 100, speed,np.nan)

    col1 = ['10','11','12']
    for key,k in zip(col1,sigmaLevels):
        # rotate vectors and calculate speed
        ur,vr = rotate_velocityField(u[k,:,:],v[k,:,:],angle)
        speed = np.sqrt(ur**2 + vr**2)
        speed = np.where(depth < 100, speed,np.nan)
        un,vn = ur/speed,vr/speed

        xplot,yplot,uplot,vplot = formatting_vectors(un,vn,lon,lat,FILE_DIR)


        a = m[key]
        cf = a.contourf(lon,lat,speed,contours,latlon=True,cmap=colormap[property])
        qv = a.quiver(xplot,yplot,uplot,vplot,scale=70,width=0.001,headwidth=4,headlength=4,latlon=True)

    axes[0,1].set_title('Experimento Controle',fontsize=8)
    axes[1,1].set_title(u'Experimento Anomalo',fontsize=8)

    # setting colorbar configuration
    cb = plt.colorbar(cf,orientation='horizontal',cax=cax,format='%0.2f')
    fig.text(0.45,0.075,colorbarTitle[property],fontsize=8)

    # title and some figure adjusts
    d = pd.to_datetime(ncin.time[timestep].values)
    plt.suptitle(u'%s nas camadas de superficie, meio e fundo, no Experimento\n' \
                  u'Controle (superior) e Anomalo (inferior) em ' \
                  '%s de %s'%(P[property],d.strftime('%d'),d.strftime('%B')),fontsize=10)
    rect = (0,0.08,1.,0.95)
    plt.tight_layout(rect=rect) # box for tight_subplot_layout
    # plt.subplots_adjust(top=0.886,bottom=0.109,left=0.054,right=0.995,hspace=0.0,wspace=0.045)
    plt.subplots_adjust(top=0.915,bottom=0.11,left=0.036,right=0.999,hspace=0.082,wspace=0.061)

    if savefig:
        savefig_dir = masterThesisPack.make_dir()
        plt.savefig(savefig_dir+'masterThesis_analysis/figures/experiments_outputs/velocity/velocity_superf_meio_fundo_timestep_%s.eps'%(str(timestep)),orientation='landscape')

    return fig,axes

def search_information(ncin,ind,nstepBegin,nstepFinal,loc,var):
    # based on an ind value, return a dictionary with informations for each
    # cross section, such as location, variable position, etc.

    # set some variables
    if var == 'temp':
        sigma = -1
        value = 18.
    if var == 'salt':
        sigma = 0
        value = 36.

    # iniatilize dictionary
    idx_begin = masterThesisPack.find_distance_of_a_value(ncin,ind,nstepBegin,sigma,var,value)[1][0][0]
    idx_final = masterThesisPack.find_distance_of_a_value(ncin,ind,nstepFinal,sigma,var,value)[1][0][0]

    info = {
        'location': loc,
        'beginPos_X': masterThesisPack.find_distance_of_a_value(ncin,ind,nstepBegin,sigma,var,value)[0],
        'finalPos_X': masterThesisPack.find_distance_of_a_value(ncin,ind,nstepFinal,sigma,var,value)[0],
        'beginPos_Z': masterThesisPack.find_depth_of_a_value(ncin,ind,nstepBegin,idx_begin,sigma,var,value),
        'finalPos_Z': masterThesisPack.find_depth_of_a_value(ncin,ind,nstepFinal,idx_final,sigma,var,value)
    }

    return info
