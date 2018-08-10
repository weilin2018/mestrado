#!/usr/bin/env python2
#-*-coding:utf-8-*-

import numpy as np # atribuir as arrays do .cdf para numpy's arrays
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset, inset_axes
import xray as xr

import cmocean as cmo

import pickle

from artigoTGpack import artigoTGpack as oceano

# ----- functions
def make_map(ax,region='bigs',resolution='i',parallels_label=[True,False,False,True],rotation=None):
    # definindo coordenadas para plotar a região desejada

    if region == 'bigs':
        llon = -44.8;
        ulon = -43.45;
        llat = -23.55;
        ulat = -22.85;

        parallels = [-22.85,-23.0,-23.2,-23.4,-23.55]
        meridians = [-44.9,-44.4,-43.9,-43.4]

    if region == 'ribeira':
        # -22.920421, -44.275737
        # -23.067316, -44.539019
        llon = -44.5390919;
        ulon = -44.275737;
        llat = -23.067316;
        ulat = -22.920421;

        parallels = [-23.065316, -22.9938685, -22.925421]
        meridians = [-44.5390919, -44.404115, -44.275737]

    # -----------------------------------------------------------------------------
    # ---------------------    PLOTANDO A BASE DO MAPA   --------------------------
    # -----------------------------------------------------------------------------
    m = Basemap(projection='merc', llcrnrlat=llat, urcrnrlat=ulat, llcrnrlon=llon, urcrnrlon=ulon, resolution=resolution)
    m.ax = ax
    # plotar outras coisas do mapa
    m.drawcoastlines(linewidth=0.1,color='#c0c0c0') #linha de costa em alta resolução
    m.drawmapboundary(fill_color='#f2f2f2') # fill_color colore o oceano
    # m.fillcontinents(color='#ffd480') # colorir o continente
    m.fillcontinents(color='0.85')
    m.drawrivers()
    # desenhar meridianos e paralelos conforme definido acima
    m.drawparallels(parallels,labels=parallels_label,fontsize=8,
                                        color='gray',linewidth=0.5,fmt='%0.1f')
    m.drawmeridians(meridians,labels=parallels_label,fontsize=8,
                                        color='gray',linewidth=0.5,fmt='%0.1f',rotation=rotation)

    return m

def create_Figure_structure(figsize=(20,15),resolution='i'):

    fig,ax = plt.subplots(nrows=2,figsize=figsize)

    m_allArea  = make_map(ax[0],region='bigs',resolution=resolution)
    m_zoomBay  = make_map(ax[1],region='ribeira',resolution=resolution)

    # divider = make_axes_locatable(ax[1])
    # cax = divider.append_axes("bottom",size='5%',pad=0.15)
    cax = fig.add_axes([0.135,0.12,0.76,0.02])

    # for orient in ['top','bottom','left','right']:
    #     m_allArea.ax.spines[orient].set_visible(False)
    #     m_zoomBay.ax.spines[orient].set_visible(False)

    return m_allArea, m_zoomBay, cax

def create_Figure_structure_2ndVersion(figsize=(20,15),resolution='i',rotation=None):

    fig,ax = plt.subplots(nrows=2,figsize=figsize)

    m_allArea  = make_map(ax[0],region='ribeira',resolution=resolution,parallels_label=[True,False,False,False],rotation=rotation)
    m_zoomBay  = make_map(ax[1],region='ribeira',resolution=resolution)

    # define position of colorbar's axis
    cax = fig.add_axes([0.14,0.12,0.765,0.02])

    return m_allArea, m_zoomBay, cax

def create_Figure_structure_4plots(figsize=(20,15),resolution='i',region='bigs',rotation=None,cax_parameters=[0.12,0.10,0.78,0.02]):

    fig,ax = plt.subplots(ncols=2,nrows=2,figsize=figsize)

    m1  = make_map(ax[0,0],region=region,resolution=resolution,parallels_label=[True,False,False,False],rotation=rotation)
    m2  = make_map(ax[0,1],region=region,resolution=resolution,parallels_label=[False,False,False,False],rotation=rotation)
    m3  = make_map(ax[1,0],region=region,resolution=resolution,rotation=rotation)
    m4  = make_map(ax[1,1],region=region,resolution=resolution,parallels_label=[False,False,False,True],rotation=rotation)

    cax = fig.add_axes(cax_parameters)

    return m1,m2,m3,m4,cax

# ----- Defining class
class Experiment(object):

    def __init__(self,fname,figsize,region='bigs'):
        """
            fname : string with directory and netcdf filename
            figsize : tuple with width and height of figures (in cm)
        """
        self.ncin  = xr.open_dataset(fname)
        self.fname = fname
        # convert figsize from cm to inches dividing by 2.54
        self.figSize = (figsize[0]/2.54, figsize[1]/2.54)
        self.region  = region

    def importVariables_basic(self):
        """ importing basic variables, as latitude and longitude """
        self.lon = self.ncin.lon.values
        self.lat = self.ncin.lat.values

        self.lon[self.lon == 0.] = np.nan
        self.lat[self.lat == 0.] = np.nan

    def importvariables_Circulation(self,i=None,sigma=None):
        """ importing zonal and meridional velocity """
        if i == None:
            u = self.ncin.u.values[:,0,:,:]
            v = self.ncin.v.values[:,0,:,:]
        else:
            u = self.ncin.u.values[i,:,:,:]
            v = self.ncin.v.values[i,:,:,:]

        # removing nan values
        u[u == np.nanmin(u)] = np.nan
        v[v == np.nanmin(v)] = np.nan

        # inserting as attribute
        self.u = u
        self.v = v

    def calculateMean_verticalVelocity(self,axis=0):
        """ calculating mean vertical velocity """
        self.umean = np.nanmean(self.u,axis=axis)
        self.vmean = np.nanmean(self.v,axis=axis)
        self.imean = np.sqrt(self.umean ** 2 + self.vmean ** 2)
        # normalizing velocity vectors by the intensity
        self.un = self.umean/self.imean
        self.vn = self.vmean/self.imean

    def calculateSpeed(self):

        self.imean = np.sqrt(self.u**2 + self.v**2)
        self.un    = self.u/self.imean
        self.vn    = self.v/self.imean

    def calcMax(self,data):
        # condition to select some specific region
        ilonlat = ((self.lon > -44.4) & (self.lon < -43.9) & (self.lat > -23.2) & (self.lat < -23))

        d = data[ilonlat]

        return np.nanmax(d)

    def cutData4Ribeira(self):
        """ recortar dados somente para a região da Baia da Ribeira """

        self.lonRib = self.lon[80:130,:30]
        self.latRib = self.lat[80:130,:30]
        self.spdRib  = self.imean[:,80:130,:30]
        # cut also velocities vectors
        self.uRib = self.un[:,80:130,:30]
        self.vRib = self.vn[:,80:130,:30]

    def define_adjustSubplot_Parameters(self,top,bottom,left,right,hspace,wspace):
        self.subAdjust_top     = top
        self.subAdjust_bottom  = bottom
        self.subAdjust_left    = left
        self.subAdjust_right   = right
        self.subAdjust_hspace  = hspace
        self.subAdjust_wspace  = wspace

    def plot_Figure_WindDriven_Circulation(self):
        """ plot figure 3 related to experiment I output - circulation """
        if ~hasattr(self,'umean'):
            self.importVariables_basic()
            self.importvariables_Circulation(i=-1)
            self.calculateMean_verticalVelocity()

        # define slice parameter, to skip vectors
        skip_all = (slice(None,None,15),slice(None,None,15))
        skip_bay = (slice(None,None,4),slice(None,None,4))

        # plot normalized vectors
        # create subplot structure
        self.m_all, self.m_rib, self.cax = create_Figure_structure(figsize=self.figSize,resolution='f')
        contour_levels = np.arange(0,0.3,0.01)

        cf = self.m_all.contourf(self.lon,self.lat,self.imean,contour_levels,latlon=True,cmap=cmo.cm.speed)
        self.m_rib.contourf(self.lon,self.lat,self.imean,contour_levels,latlon=True,cmap=cmo.cm.speed)

        self.m_all.quiver(self.lon[skip_all],self.lat[skip_all],self.un[skip_all],self.vn[skip_all],latlon=True,scale=30,minshaft=2)
        self.m_rib.quiver(self.lon[skip_bay],self.lat[skip_bay],self.un[skip_bay],self.vn[skip_bay],latlon=True,scale=20,minshaft=4)

        # plot colorbar using the first contourf
        cb = plt.colorbar(cf,orientation='horizontal',ticks=[0.0,0.1,0.2],cax=self.cax,format='%.1f')
        cb.set_label('Average Circulation'+ ' (ms' +r'$^{-1}$)',fontsize=10)

        # insert panels identification (a) and (b)
        self.m_all.ax.text(0.03,0.85,'(a)',transform=self.m_all.ax.transAxes)
        self.m_rib.ax.text(0.03,0.85,'(b)',transform=self.m_rib.ax.transAxes)

        rect = (0.14,0.13,1.,1.)
        plt.tight_layout(rect=rect)
        plt.subplots_adjust(top=1.0,bottom=0.185,left=0.130,right=0.90,wspace=0.,hspace=0.105)

        plt.savefig('/media/danilo/Danilo/mestrado/github/artigoTG/figures/%s.eps'%(self.figname))

    def plot_Figure5(self,spring_flood,spring_ebb,neap_flood,neap_ebb):
        """ this function plot only the surface circulation """
        if ~hasattr(self,'imean'):
            self.importVariables_basic()
            self.importvariables_Circulation() # import all timestep data
            self.calculateSpeed()

        # define contours
        contour_levels = np.arange(0,0.7,0.001)

        # plot only normalized vector [un,vn]
        self.m1,self.m2,self.m3,self.m4, self.cax = create_Figure_structure_4plots(figsize=self.figSize,resolution='f',region=self.region)

        # plotting spring flood
        cf = self.m1.contourf(self.lon,self.lat,self.imean[spring_flood,:,:],contour_levels,latlon=True,cmap=cmo.cm.speed)
        qv = self.m1.quiver(self.lon[::15,::15],self.lat[::15,::15],self.un[spring_flood,::15,::15],self.vn[spring_flood,::15,::15],latlon=True)
        # self.m1.ax.set_title('Spring Flood, with max of %0.2f'%(self.calcMax(np.squeeze(self.imean[spring_flood,:,:]))))
        maxValue = self.calcMax(np.squeeze(self.imean[spring_flood,:,:]))
        textMax = 'Max Curr\n%0.2f'%(maxValue) + 'ms' +r'$^{-1}$'
        self.m1.ax.text(0.90,0.91,textMax,transform=self.m1.ax.transAxes,ha='center',va='center',fontsize=8)
        self.m1.ax.text(0.03,0.85,'(a)',transform=self.m1.ax.transAxes)
        self.m1.ax.set_title('Spring Flood',fontsize=10)

        cb = plt.colorbar(cf,orientation='horizontal',ticks=[0.0,0.2,0.4,0.6],cax=self.cax,format='%.1f')
        cb.set_label('Surface Circulation'+r' (ms$^{-1}$)',fontsize=10,labelpad=-1)

        # plotting spring ebb
        self.m2.contourf(self.lon,self.lat,self.imean[spring_ebb,:,:],contour_levels,latlon=True,cmap=cmo.cm.speed)
        self.m2.quiver(self.lon[::15,::15],self.lat[::15,::15],self.un[spring_ebb,::15,::15],self.vn[spring_ebb,::15,::15],latlon=True)
        # self.m2.ax.set_title('Spring Ebb, with max of %0.2f'%(self.calcMax(np.squeeze(self.imean[spring_ebb,:,:]))))
        maxValue = self.calcMax(np.squeeze(self.imean[spring_ebb,:,:]))
        textMax = 'Max Curr\n%0.2f'%(maxValue) + 'ms' +r'$^{-1}$'
        self.m2.ax.text(0.90,0.91,textMax,transform=self.m2.ax.transAxes,ha='center',va='center',fontsize=8)
        self.m2.ax.text(0.03,0.85,'(b)',transform=self.m2.ax.transAxes)
        self.m2.ax.set_title('Spring Ebb',fontsize=10)

        # plotting neap flood
        self.m3.contourf(self.lon,self.lat,self.imean[neap_flood,:,:],contour_levels,latlon=True,cmap=cmo.cm.speed)
        self.m3.quiver(self.lon[::15,::15],self.lat[::15,::15],self.un[neap_flood,::15,::15],self.vn[neap_flood,::15,::15],latlon=True)
        maxValue = self.calcMax(np.squeeze(self.imean[neap_flood,:,:]))
        textMax = 'Max Curr\n%0.2f'%(maxValue) + 'ms' +r'$^{-1}$'
        self.m3.ax.text(0.90,0.91,textMax,transform=self.m3.ax.transAxes,ha='center',va='center',fontsize=8)
        self.m3.ax.text(0.03,0.85,'(c)',transform=self.m3.ax.transAxes)
        self.m3.ax.set_title('Neap Flood',fontsize=10)

        # plotting neap ebb
        self.m4.contourf(self.lon,self.lat,self.imean[neap_ebb,:,:],contour_levels,latlon=True,cmap=cmo.cm.speed)
        self.m4.quiver(self.lon[::15,::15],self.lat[::15,::15],self.un[neap_ebb,::15,::15],self.vn[neap_ebb,::15,::15],latlon=True)
        maxValue = self.calcMax(np.squeeze(self.imean[neap_ebb,:,:]))
        textMax = 'Max Curr\n%0.2f'%(maxValue) + 'ms' +r'$^{-1}$'
        self.m4.ax.text(0.90,0.91,textMax,transform=self.m4.ax.transAxes,ha='center',va='center',fontsize=8)
        self.m4.ax.text(0.03,0.85,'(d)',transform=self.m4.ax.transAxes)
        self.m4.ax.set_title('Neap Ebb',fontsize=10)

        if hasattr(self,'subAdjust_left'):
            rect = (0.13,0.16,1.,1.)
            plt.tight_layout(rect=rect)
            plt.subplots_adjust(top=self.subAdjust_top,bottom=self.subAdjust_bottom,left=self.subAdjust_left,right=self.subAdjust_right,wspace=self.subAdjust_wspace,hspace=self.subAdjust_hspace)

        if hasattr(self,'figname'):
            plt.savefig('/media/danilo/Danilo/mestrado/github/artigoTG/figures/%s.png'%(self.figname),dpi=600)
        else:
            plt.show()

    def plot_Figure6(self,spring_flood,spring_ebb,neap_flood,neap_ebb):
        """ this function plot only the surface circulation """
        if ~hasattr(self,'imean'):
            self.importVariables_basic()
            self.importvariables_Circulation() # import all timestep data
            self.calculateSpeed()
            self.cutData4Ribeira()

        # define contours
        contour_levels = np.arange(0,0.7,0.001)

        self.m1,self.m2,self.m3,self.m4, self.cax = create_Figure_structure_4plots(figsize=self.figSize,resolution='f',region='ribeira',rotation=20,cax_parameters=[0.12,0.10,0.78,0.02])
        # plotting spring flood
        cf = self.m1.contourf(self.lonRib,self.latRib,self.spdRib[spring_flood,:,:],contour_levels,latlon=True,cmap=cmo.cm.speed)
        qv = self.m1.quiver(self.lonRib[::4,::4],self.latRib[::4,::4],self.uRib[spring_flood,::4,::4],self.vRib[spring_flood,::4,::4],latlon=True)
        # self.m1.ax.set_title('Spring Flood, with max of %0.2f'%(self.calcMax(np.squeeze(self.s[spring_flood,:,:]))))
        maxValue = np.nanmax(self.spdRib[spring_flood,:,:])
        textMax = 'Max Curr\n%0.2f'%(maxValue) + 'ms' +r'$^{-1}$'
        self.m1.ax.text(0.90,0.90,textMax,transform=self.m1.ax.transAxes,ha='center',va='center',fontsize=8)
        self.m1.ax.text(0.03,0.85,'(a)',transform=self.m1.ax.transAxes)
        self.m1.ax.set_title('Spring Flood',fontsize=10,y=.97)

        cb = plt.colorbar(cf,orientation='horizontal',ticks=[0.0,0.2,0.4,0.6],cax=self.cax,format='%.1f')
        cb.set_label('Surface Circulation'+r' (ms$^{-1}$)',fontsize=10,labelpad=-1)

        # plotting spring ebb
        self.m2.contourf(self.lonRib,self.latRib,self.spdRib[spring_ebb,:,:],contour_levels,latlon=True,cmap=cmo.cm.speed)
        self.m2.quiver(self.lonRib[::4,::4],self.latRib[::4,::4],self.uRib[spring_ebb,::4,::4],self.vRib[spring_ebb,::4,::4],latlon=True)
        maxValue = np.nanmax(self.spdRib[spring_ebb,:,:])
        textMax = 'Max Curr\n%0.2f'%(maxValue) + 'ms' +r'$^{-1}$'
        self.m2.ax.text(0.90,0.90,textMax,transform=self.m2.ax.transAxes,ha='center',va='center',fontsize=8)
        self.m2.ax.text(0.03,0.85,'(b)',transform=self.m2.ax.transAxes)
        self.m2.ax.set_title('Spring Ebb',fontsize=10,y=.97)

        # plotting neap flood
        self.m3.contourf(self.lonRib,self.latRib,self.spdRib[neap_flood,:,:],contour_levels,latlon=True,cmap=cmo.cm.speed)
        self.m3.quiver(self.lonRib[::4,::4],self.latRib[::4,::4],self.uRib[neap_flood,::4,::4],self.vRib[neap_flood,::4,::4],latlon=True)
        maxValue = np.nanmax(self.spdRib[neap_flood,:,:])
        textMax = 'Max Curr\n%0.2f'%(maxValue) + 'ms' +r'$^{-1}$'
        self.m3.ax.text(0.90,0.90,textMax,transform=self.m3.ax.transAxes,ha='center',va='center',fontsize=8)
        self.m3.ax.text(0.03,0.85,'(c)',transform=self.m3.ax.transAxes)
        self.m3.ax.set_title('Neap Flood',fontsize=10,y=.97)

        # plotting neap ebb
        self.m4.contourf(self.lonRib,self.latRib,self.spdRib[neap_ebb,:,:],contour_levels,latlon=True,cmap=cmo.cm.speed)
        self.m4.quiver(self.lonRib[::4,::4],self.latRib[::4,::4],self.uRib[neap_ebb,::4,::4],self.vRib[neap_ebb,::4,::4],latlon=True)
        maxValue = np.nanmax(self.spdRib[neap_ebb,:,:])
        textMax = 'Max Curr\n%0.2f'%(maxValue) + 'ms' +r'$^{-1}$'
        self.m4.ax.text(0.90,0.90,textMax,transform=self.m4.ax.transAxes,ha='center',va='center',fontsize=8)
        self.m4.ax.text(0.03,0.85,'(d)',transform=self.m4.ax.transAxes)
        self.m4.ax.set_title('Neap Ebb',fontsize=10,y=.97)

        if hasattr(self,'subAdjust_left'):
            rect = (0.13,0.11,1.,1.)
            plt.tight_layout(rect=rect)
            plt.subplots_adjust(top=self.subAdjust_top,bottom=self.subAdjust_bottom,left=self.subAdjust_left,right=self.subAdjust_right,wspace=self.subAdjust_wspace,hspace=self.subAdjust_hspace)

        if hasattr(self,'figname'):
            plt.savefig('/media/danilo/Danilo/mestrado/github/artigoTG/figures/%s.png'%(self.figname),dpi=600)
        else:
            plt.show()

    def plot_Figure6_2ndVersion(self,spring_flood,spring_ebb,neap_flood,neap_ebb):
        if ~hasattr(self,'imean'):
            self.importVariables_basic()
            self.importvariables_Circulation() # import all timestep data
            self.calculateSpeed()
            self.cutData4Ribeira()

        # define contours
        contour_levels = np.arange(0,0.3,0.001)

        self.m1,self.m2,self.cax = create_Figure_structure_2ndVersion(figsize=self.figSize,rotation=20,resolution='f')

        cf = self.m1.contourf(self.lonRib,self.latRib,self.spdRib[spring_ebb,:,:],contour_levels,latlon=True,cmap=cmo.cm.speed)
        qv = self.m1.quiver(self.lonRib[::4,::4],self.latRib[::4,::4],self.uRib[spring_ebb,::4,::4],self.vRib[spring_ebb,::4,::4],latlon=True)
        # self.m1.ax.set_title('Spring Flood, with max of %0.2f'%(self.calcMax(np.squeeze(self.s[spring_flood,:,:]))))
        maxValue = np.nanmax(self.spdRib[spring_ebb,:,:])
        textMax = 'Max Curr\n%0.2f'%(maxValue) + 'ms' +r'$^{-1}$'
        self.m1.ax.text(0.88,0.88,textMax,transform=self.m1.ax.transAxes,ha='center',va='center',fontsize=8)
        self.m1.ax.text(0.03,0.85,'(a)',transform=self.m1.ax.transAxes)

        cb = plt.colorbar(cf,orientation='horizontal',ticks=[0.0,0.1,0.2],cax=self.cax,format='%.1f')
        cb.set_label('Surface Circulation'+r' (ms$^{-1}$)',fontsize=10,labelpad=-1)

        self.m2.contourf(self.lonRib,self.latRib,self.spdRib[neap_ebb,:,:],contour_levels,latlon=True,cmap=cmo.cm.speed)
        self.m2.quiver(self.lonRib[::4,::4],self.latRib[::4,::4],self.uRib[neap_ebb,::4,::4],self.vRib[neap_ebb,::4,::4],latlon=True)
        maxValue = np.nanmax(self.spdRib[neap_ebb,:,:])
        textMax = 'Max Curr\n%0.2f'%(maxValue) + 'ms' +r'$^{-1}$'
        self.m2.ax.text(0.88,0.88,textMax,transform=self.m2.ax.transAxes,ha='center',va='center',fontsize=8)
        self.m2.ax.text(0.03,0.85,'(b)',transform=self.m2.ax.transAxes)

        if hasattr(self,'subAdjust_left'):
            rect = (0.13,0.11,1.,1.)
            plt.tight_layout(rect=rect)
            plt.subplots_adjust(top=self.subAdjust_top,bottom=self.subAdjust_bottom,left=self.subAdjust_left,right=self.subAdjust_right,wspace=self.subAdjust_wspace,hspace=self.subAdjust_hspace)

        if hasattr(self,'figname'):
            plt.savefig('/media/danilo/Danilo/mestrado/github/artigoTG/figures/%s.eps'%(self.figname))
        else:
            plt.show()

# ----- load file from Oceano Computer ----- #
DATA_DIR = '/media/danilo/Danilo/mestrado/artigo_data/simulacoes/sims_dispersao/'

# ----- showing information in a Basemap instance
def fig3():
    expI = Experiment(DATA_DIR+'expI.cdf',figsize=(8.4,10.))
    expI.figname = 'Fig3'
    expI.plot_Figure_WindDriven_Circulation()

def fig4():
    expII = Experiment(DATA_DIR+'expII.cdf',figsize=(8.4,10.))
    expII.figname = 'Fig4'
    expII.plot_Figure_WindDriven_Circulation()

def fig5():
    expIII = Experiment(DATA_DIR+'expIII.cdf',figsize=(17.4,12))
    expIII.figname = 'Fig5'
    # define some values for tight_subplot_layout
    expIII.subAdjust_top    = 0.973
    expIII.subAdjust_bottom = 0.132
    expIII.subAdjust_left   = 0.057
    expIII.subAdjust_right  = 0.995
    expIII.subAdjust_hspace = 0.003
    expIII.subAdjust_wspace = 0.047
    expIII.importVariables_basic() # import lat and lon
    expIII.importvariables_Circulation() # import velocities components
    expIII.calculateSpeed() # calculate speed and normalize vectors
    expIII.plot_Figure5(229,218,290,293)

def fig6():
    expIII = Experiment(DATA_DIR+'expIII.cdf',figsize=(17.4,12))
    expIII.importVariables_basic()
    expIII.importvariables_Circulation()
    expIII.calculateSpeed()
    expIII.cutData4Ribeira()
    # if you want to plot 4 subplots for each tidal condition and phase, uncomment the 7 lines below
    # expIII.figSize = (17.4/2.54,12./2.54)
    # expIII.figname = 'Fig6'
    # expIII.subAdjust_top    = 0.977
    # expIII.subAdjust_bottom = 0.172
    # expIII.subAdjust_left   = 0.072
    # expIII.subAdjust_right  = 0.973
    # expIII.subAdjust_hspace = 0.072
    # expIII.subAdjust_wspace = 0.082
    # expIII.plot_Figure6(229,218,290,293)

    # if you want to plot 4 subplots for each tidal condition and phase, uncomment the 7 lines below
    expIII.figSize = (8.4/2.54,10./2.54)
    expIII.figname = 'Fig6_2ndVersion'
    expIII.define_adjustSubplot_Parameters(top=0.977,bottom=0.172,left=0.072,right=0.973,hspace=0.072,wspace=0.082)
    expIII.plot_Figure6_2ndVersion(229,218,290,293)




# -------------- TESTING ZONE
#
#
# # testing plot experiment III / Figure 5
# expIII = Experiment(DATA_DIR+'expIII.cdf',figsize=(8.4,12))
# expIII.subAdjust_top    = 0.977
# expIII.subAdjust_bottom = 0.172
# expIII.subAdjust_left   = 0.072
# expIII.subAdjust_right  = 0.973
# expIII.subAdjust_hspace = 0.072
# expIII.subAdjust_wspace = 0.082
# expIII.importVariables_basic()
# expIII.importvariables_Circulation()
# expIII.calculateSpeed()
# expIII.cutData4Ribeira()
# # expIII.calculateMean_verticalVelocity(axis=1)
#
#
# expIII.m1,expIII.m3,expIII.cax = create_Figure_structure_2ndVersion(figsize=expIII.figSize)
#
#
#
#
# spring_flood = 229 #229 # 66.9 cm/s (210)
# spring_ebb   = 218 #218 #
# neap_flood   = 290 #290 # 0.28 m/s
# neap_ebb     = 293 #293 # 0.12 m/s
#
# # define contours
# contour_levels = np.arange(0,0.7,0.001)
#
# # expIII.m1,expIII.m2,expIII.m3,expIII.m4, expIII.cax = create_Figure_structure_4plots(figsize=expIII.figSize,resolution='i',region='ribeira',rotation=20,cax_parameters=[0.12,0.10,0.78,0.02])
# # plotting spring flood
# cf = expIII.m1.contourf(expIII.lonRib,expIII.latRib,expIII.spdRib[spring_flood,:,:],contour_levels,latlon=True,cmap=cmo.cm.speed)
# qv = expIII.m1.quiver(expIII.lonRib[skip_all],expIII.latRib[skip_all],expIII.uRib[skip_all2],expIII.vRib[skip_all2],latlon=True)
# # expIII.m1.ax.set_title('Spring Flood, with max of %0.2f'%(expIII.calcMax(np.squeeze(expIII.s[spring_flood,:,:]))))
# maxValue = np.nanmax(expIII.spdRib[spring_flood,:,:])
# textMax = 'Max Curr\n%0.2f'%(maxValue) + 'ms' +r'$^{-1}$'
# expIII.m1.ax.text(0.90,0.90,textMax,transform=expIII.m1.ax.transAxes,ha='center',va='center',fontsize=8)
# expIII.m1.ax.text(0.03,0.85,'(a)',transform=expIII.m1.ax.transAxes)
# expIII.m1.ax.set_title('Spring Flood',fontsize=10,y=.97)
#
# cb = plt.colorbar(cf,orientation='horizontal',ticks=[0.0,0.2,0.4,0.6],cax=expIII.cax,format='%.1f')
# cb.set_label('Surface Circulation'+r' [m s$^{-1}$]',fontsize=10,labelpad=-1)
#
# # plotting spring ebb
# expIII.m2.contourf(expIII.lonRib,expIII.latRib,expIII.spdRib[spring_ebb,:,:],contour_levels,latlon=True,cmap=cmo.cm.speed)
# expIII.m2.quiver(expIII.lonRib[::4,::4],expIII.latRib[::4,::4],expIII.uRib[spring_ebb,::4,::4],expIII.vRib[spring_ebb,::4,::4],latlon=True)
# maxValue = np.nanmax(expIII.spdRib[spring_ebb,:,:])
# textMax = 'Max Curr\n%0.2f'%(maxValue) + 'm s' +r'$^{-1}$'
# expIII.m2.ax.text(0.90,0.90,textMax,transform=expIII.m2.ax.transAxes,ha='center',va='center',fontsize=8)
# expIII.m2.ax.text(0.03,0.85,'(b)',transform=expIII.m2.ax.transAxes)
# expIII.m2.ax.set_title('Spring Ebb',fontsize=10,y=.97)
#
# # plotting neap flood
# expIII.m3.contourf(expIII.lonRib,expIII.latRib,expIII.spdRib[neap_flood,:,:],contour_levels,latlon=True,cmap=cmo.cm.speed)
# expIII.m3.quiver(expIII.lonRib[::4,::4],expIII.latRib[::4,::4],expIII.uRib[neap_flood,::4,::4],expIII.vRib[neap_flood,::4,::4],latlon=True)
# maxValue = np.nanmax(expIII.spdRib[neap_flood,:,:])
# textMax = 'Max Curr\n%0.2f'%(maxValue) + 'm s' +r'$^{-1}$'
# expIII.m3.ax.text(0.90,0.90,textMax,transform=expIII.m3.ax.transAxes,ha='center',va='center',fontsize=8)
# expIII.m3.ax.text(0.03,0.85,'(c)',transform=expIII.m3.ax.transAxes)
# expIII.m3.ax.set_title('Neap Flood',fontsize=10,y=.97)
#
# # plotting neap ebb
# expIII.m4.contourf(expIII.lonRib,expIII.latRib,expIII.spdRib[neap_ebb,:,:],contour_levels,latlon=True,cmap=cmo.cm.speed)
# expIII.m4.quiver(expIII.lonRib[::4,::4],expIII.latRib[::4,::4],expIII.uRib[neap_ebb,::4,::4],expIII.vRib[neap_ebb,::4,::4],latlon=True)
# maxValue = np.nanmax(expIII.spdRib[neap_ebb,:,:])
# textMax = 'Max Curr\n%0.2f'%(maxValue) + 'm s' +r'$^{-1}$'
# expIII.m4.ax.text(0.90,0.90,textMax,transform=expIII.m4.ax.transAxes,ha='center',va='center',fontsize=8)
# expIII.m4.ax.text(0.03,0.85,'(d)',transform=expIII.m4.ax.transAxes)
# expIII.m4.ax.set_title('Neap Ebb',fontsize=10,y=.97)
#
# rect = (0.13,0.11,1.,1.)
# plt.tight_layout(rect=rect)
# plt.subplots_adjust(top=0.975,bottom=0.197,left=0.057,right=0.995,wspace=0.047,hspace=0.00)
