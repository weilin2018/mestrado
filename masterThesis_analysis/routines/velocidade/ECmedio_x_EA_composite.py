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

def export_data(fname,timestep=0,convertData=False,outputFile=None):
    # plotting climatologic data: t = 0, k = 0
    ncin = xr.open_dataset(fname)

    lon,lat = ncin.lon.values.copy(), ncin.lat.values.copy()
    angle = ncin.ang.values.copy()
    depth = ncin.depth.values.copy()
    sigma = ncin.sigma.values.copy()
    lon[lon == 0.] = np.nan
    lat[lat == 0.] = np.nan

    # transform data from sigma vertical coordinates to standard level
    if convertData:
	os.system('converting u and v for %s in timestep %i'%(exp,np.nanmean(timestep)))
        # interpolating data from sigma to standard level
        u,v = ncin.u[timestep,:,:,:],ncin.v[timestep,:,:,:]
        stdl,newu = oceano.sigma2stdl(u,sigma,23,depth,ncin.h1,lon,lat,'u component')
        umean = np.nanmean(newu,axis=0) # media no tempo
        np.save(outputFile+'%s_umean_%i.npy'%(exp,np.nanmean(timestep)),umean)

        stdl,newv = oceano.sigma2stdl(v,sigma,23,depth,ncin.h1,lon,lat,'v component')
        vmean = np.nanmean(newv,axis=0) # media no tempo
        np.save(outputFile+'%s_vmean_%i.npy'%(exp,np.nanmean(timestep)),vmean)

        u,v = umean,vmean
    else:
        u = np.load(outputFile+'%s_umean_%i.npy'%(exp,np.nanmean(timestep)))
        v = np.load(outputFile+'%s_vmean_%i.npy'%(exp,np.nanmean(timestep)))

        stdl = [0, 10, 25, 40, 50, 60, 70, 80, 100, 150, 300, 350, 400, 450, 500, 600, 700, 800, 1000, 1200, 1500, 1800, 2000]

    return lon,lat,u,v,depth,angle,stdl

def rotate_velocityField(u,v,ang):

    import decomp
    ur = np.zeros(u.shape)*np.nan
    vr = np.zeros(v.shape)*np.nan

    for j in range(u.shape[0]):
        U,V = u[j,:],v[j,:]
        angle = ang[j,:]

        INT,DIR = decomp.uv2intdir(U,V,0,angle)
        uro,vro = decomp.intdir2uv(INT,DIR,0,angle)
        ur[j,:] = uro
        vr[j,:] = vro

    return ur,vr

def tratando_corrente(u,v,depth,angle):

    # ur,vr = rotate_velocityField(u,v,angle)
    spd = np.sqrt(u**2+v**2)
    spd = np.where(depth < 200, spd,np.nan)

    return u,v,spd

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
"""
    Importantissimo:

    Na eventualidade de rodar as simulações novamente e na ocasião de ser com novas
    configurações, é necessário converter os dados do modelo de coordenadas
    sigmas para standard levels novamente, atualizando assim os arquivos utilizados
    nesta rotina para plotagem.

    Importante manter isso em mente, pois se não for feito isso, a rotina continuará
    plotando os dados antigos.

    Para que isso ocorra, ative a seguinte opção na linha 101:
        convertData = True
"""
# beginnig of the main code
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/output/')
FILE_DIR = BASE_DIR+'masterThesis_analysis/routines/index_list.npy'
os.system('clear')

convertData = False

exp = raw_input('Digite o experimento a ser plotado: ')
fname = DATA_DIR + exp +'.cdf'

plt.ion()

# o primeiro timestep sera para o controle, o segundo para o anomalo
timestep = [np.arange(65,289,1),np.arange(280,289,1)]

for nstep in timestep:
    outputFile = DATA_DIR + 'dados_interpolados_stdl/'
    lon,lat,u,v,depth,angle,stdl = export_data(fname,timestep=nstep,convertData=convertData,outputFile=outputFile)

    # rotating vectors
    ur_surf,vr_surf,spd_surf = tratando_corrente(u[0,:,:],v[0,:,:],depth,angle)
    ur_meio,vr_meio,spd_meio = tratando_corrente(u[3,:,:],v[3,:,:],depth,angle)
    ur_fund,vr_fund,spd_fund = tratando_corrente(u[7,:,:],v[7,:,:],depth,angle)

    # formatando os dados para um formato de visualizacao melhor e passando os vetores
    # normalizados pela velocidade
    xplot,yplot,uplot_surf,vplot_surf = ocplt.formatting_vectors(ur_surf/spd_surf,vr_surf/spd_surf,lon,lat,FILE_DIR)
    xplot,yplot,uplot_meio,vplot_meio = ocplt.formatting_vectors(ur_meio/spd_meio,vr_meio/spd_meio,lon,lat,FILE_DIR)
    xplot,yplot,uplot_fund,vplot_fund = ocplt.formatting_vectors(ur_fund/spd_fund,vr_fund/spd_fund,lon,lat,FILE_DIR)

    # mask first and last 5 rows
    spd_surf[:5,:] = np.nan
    spd_surf[-5:,:] = np.nan
    spd_meio[:5,:] = np.nan
    spd_meio[-5:,:] = np.nan
    spd_fund[:5,:] = np.nan
    spd_fund[-5:,:] = np.nan
    uplot_surf[:3,:],vplot_surf[:3,:] = np.nan,np.nan
    uplot_surf[-3:,:],vplot_surf[-3:,:] = np.nan,np.nan
    uplot_meio[:3,:],vplot_meio[:3,:] = np.nan,np.nan
    uplot_meio[-3:,:],vplot_meio[-3:,:] = np.nan,np.nan
    uplot_fund[:3,:],vplot_fund[:3,:] = np.nan,np.nan
    uplot_fund[-3:,:],vplot_fund[-3:,:] = np.nan,np.nan

    # plotando no grafico
    fig = plt.figure(figsize=(12/2.54,12/2.54))
    fig.patch.set_visible(False)

    gs = gridspec.GridSpec(3,3)

    # creating axis
    ax1 = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[1,1])
    ax3 = plt.subplot(gs[2,2])
    # removing outside border
    ax1.axis('off')
    ax2.axis('off')
    ax3.axis('off')

    # creating basemap instance
    m1,meridians,parallels = make_map(ax1,ulon=np.nanmax(lon)-.5,llon=np.nanmin(lon)-.2,ulat=np.nanmax(lat)+.2,llat=np.nanmin(lat))
    m2,_,_ = make_map(ax2,ulon=np.nanmax(lon)-.5,llon=np.nanmin(lon)-.2,ulat=np.nanmax(lat)+.2,llat=np.nanmin(lat))
    m3,_,_ = make_map(ax3,ulon=np.nanmax(lon)-.5,llon=np.nanmin(lon)-.2,ulat=np.nanmax(lat)+.2,llat=np.nanmin(lat))

    # positiong each axes
    ax1.set_position([.03,.46,.6,.5])
    ax2.set_position([.1,.28,.6,.5])
    ax3.set_position([.17,.1,.6,.5])

    contours = np.arange(0,1.1,0.1)

    # plotando velocidade
    cf1 = m1.contourf(lon,lat,spd_surf,contours,cmap=cmo.cm.speed,latlon=True,rasterized=True,extend='max')
    cf2 = m2.contourf(lon,lat,spd_meio,contours,cmap=cmo.cm.speed,latlon=True,rasterized=True,extend='max')
    cf3 = m3.contourf(lon,lat,spd_fund,contours,cmap=cmo.cm.speed,latlon=True,rasterized=True,extend='max')

    # cr1 = m1.contour(lon,lat,depth,levels=[40,80],colors=('#c0c0c0'),linewidths=(0.1,0.1),latlon=True,alpha=.3)
    # cr2 = m2.contour(lon,lat,depth,levels=[40,80],colors=('#c0c0c0'),linewidths=(0.1,0.1),latlon=True,alpha=.3)
    # cr3 = m3.contour(lon,lat,depth,levels=[40,80],colors=('#c0c0c0'),linewidths=(0.1,0.1),latlon=True,alpha=.3)
    # plotando vetores
    qv1 = m1.quiver(xplot,yplot,uplot_surf,vplot_surf,scale=60,width=0.0015,headwidth=4,headlength=4,latlon=True)
    qv2 = m2.quiver(xplot,yplot,uplot_meio,vplot_meio,scale=60,width=0.0015,headwidth=4,headlength=4,latlon=True)
    qv3 = m3.quiver(xplot,yplot,uplot_fund,vplot_fund,scale=60,width=0.0015,headwidth=4,headlength=4,latlon=True)

    # inserting some texts
    ax1.text(648156,895737,'00 m',fontsize=8,va='center',ha='center')
    ax2.text(639807,895573,'40 m',fontsize=8,va='center',ha='center')
    ax3.text(647828,887225,'80 m',fontsize=8,va='center',ha='center')

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
    cax = fig.add_axes([.32,.22,.35,.02])
    cbar = plt.colorbar(cf1,orientation='horizontal',cax=cax,format='%0.1f')

    # # figure's title
    # plt.suptitle(u'Temperatura nas camadas de superfície, meio e fundo no Experimento Controle (esquerda) e Anômalo (direita) em 14 de Janeiro')

    # setting colorbar tick labels
    from matplotlib import ticker
    tick_locator = ticker.MaxNLocator(nbins=6)
    cbar.locator = tick_locator
    cbar.update_ticks()

    cbar.ax.axes.tick_params(axis='both',which='both',labelsize=8)
    cbar.ax.set_title(r'Velocidade (m s$^{-1}$)',fontsize=8)

    output_fname = fname.split('/')[-1].replace('.cdf','_'+str(int(np.mean(nstep))))
    plt.savefig('/home/danilo/Dropbox/mestrado/figuras/composicao/std_level/speed/%s/%s__17Jans.eps'%(exp,output_fname))
