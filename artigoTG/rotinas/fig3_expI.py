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
import pandas as pd
import cmocean as cmo
import os

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

def nearest_date(items,pivot):
    nearest=min(items, key=lambda x: abs(x - pivot))
    timedelta = abs(nearest - pivot)
    return nearest, timedelta

def locate_closest_date(time,day):

    # convert into datetimeindex
    dates = pd.DatetimeIndex(time)

    pivot = dates[0] + pd.DateOffset(days=day)

    # locate the closest datetime index
    nearestDate = nearest_date(dates,pivot)
    index = np.where(dates == nearestDate[0])

    return index

def convertData2Bq(c,unit='mg/L'):

    # constante = 0.3e+22
    c[c<0.0] = np.nan # valores negativos são removidos
    if unit == 'mg/L':
        # converter de mg para g
        c = c*.1e-3
    else:
        # converter de kg para g
        c = c*.1e3

    c = c*35334440479135.016 # converter de g para Bq

    return c

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

    def extract_concentration(self,time_i):
        """ extract concentration from time_i, calculate vertical integration
        and mask negative values [bad flag] """
        conc = np.nansum(self.ncin.conc[time_i,:,:,:], axis=0)
        conc[conc < 0] = np.nan

        # convert to Bq
        # k = convertData2Bq(conc)

        return conc

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

    def plot_Figure4(self,spring_flood,spring_ebb,neap_flood,neap_ebb):
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

    def animation_expIII(self):

        contour_levels = np.arange(0,0.7,0.001)

        fig,ax = plt.subplots(figsize=self.figSize)
        cax = fig.add_axes([0.135,0.12,0.76,0.02])

        for i in range(len(self.ncin.time.values)):
            ax.clear()
            cax.clear()
            m = make_map(ax,region='bigs',resolution='i')

            cf = m.contourf(self.lon,self.lat,self.imean[i,:,:],contour_levels, latlon=True,cmap=cmo.cm.speed)
            qv = m.quiver(self.lon[::15,::15],self.lat[::15,::15],self.un[i,::15,::15],self.vn[i,::15,::15],latlon=True)
            # self.m1.ax.set_title('Spring Flood, with max of %0.2f'%(self.calcMax(np.squeeze(self.imean[spring_flood,:,:]))))
            maxValue = self.calcMax(np.squeeze(self.imean[i,:,:]))
            textMax = 'Max Curr\n%0.2f'%(maxValue) + 'ms' +r'$^{-1}$'
            m.ax.text(0.90,0.91,textMax,transform=m.ax.transAxes,ha='center',va='center',fontsize=8)
            m.ax.text(0.03,0.85,'(a)',transform=m.ax.transAxes)
            m.ax.set_title('%s'%(str(self.ncin.time[i])),fontsize=10)

            cb = plt.colorbar(cf,orientation='horizontal',ticks=[0.0,0.2,0.4,0.6],cax=cax,format='%.1f')
            cb.set_label('Surface Circulation'+r' (ms$^{-1}$)',fontsize=10,labelpad=-1)

            if hasattr(self,'figname'):
                plt.savefig('/media/danilo/Danilo/mestrado/github/artigoTG/figures/animation/%s.png'%(str(i).zfill(3)),dpi=600)
                os.system('convert -trim %s %s' % ('/media/danilo/Danilo/mestrado/github/artigoTG/figures/animation/%s.png'%(str(i).zfill(3)),'/media/danilo/Danilo/mestrado/github/artigoTG/figures/animation/%s.png'%(str(i).zfill(3))))
            else:
                plt.pause(0.5)

    def plot_Figure5(self,spring_flood,spring_ebb,neap_flood,neap_ebb):
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

    def plot_Figure5_2ndVersion(self,spring_flood,spring_ebb,neap_flood,neap_ebb):
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

    def plot_Concentration(self,time1,time2,time3,time4):
        if ~hasattr(self,'lon'):
            self.importVariables_basic()

        # find index for each time passed as argument (in days)
        timestep_1 = locate_closest_date(self.ncin.time.values,time1)[0][0]
        timestep_2 = locate_closest_date(self.ncin.time.values,time2)[0][0]
        timestep_3 = locate_closest_date(self.ncin.time.values,time3)[0][0]
        timestep_4 = locate_closest_date(self.ncin.time.values,time4)[0][0]

        # load vertical integrated concentration
        self.conc1 = self.extract_concentration(timestep_1)*100
        self.conc2 = self.extract_concentration(timestep_2)*100
        self.conc3 = self.extract_concentration(timestep_3)*100
        self.conc4 = self.extract_concentration(timestep_4)*100

        # create contour_levels [0 to 100]
        contour_levels = np.arange(0,50.,0.01)

        # create structure
        self.m1,self.m2,self.m3,self.m4,self.cax = create_Figure_structure_4plots(figsize=self.figSize,resolution=self.resolution,cax_parameters=[0.125,0.12,0.75,0.02])

        cf = self.m1.contourf(self.lon,self.lat,self.conc1,contour_levels,latlon=True,cmap=cmo.cm.matter,extend='max')
        self.m1.ax.text(0.9,0.91,'%i days'%(time1),transform=self.m1.ax.transAxes,ha='center',va='center',fontsize=8)
        self.m1.ax.text(0.03,0.85,'(a)',transform=self.m1.ax.transAxes)

        cbar = plt.colorbar(cf,cax=self.cax,orientation='horizontal')
        cbar.set_label('Vertical Integrated Concentration (%)')

        self.m2.contourf(self.lon,self.lat,self.conc2,contour_levels,latlon=True,cmap=cmo.cm.matter,extend='max')
        self.m2.ax.text(0.9,0.91,'%i days'%(time2),transform=self.m2.ax.transAxes,ha='center',va='center',fontsize=8)
        self.m2.ax.text(0.03,0.85,'(b)',transform=self.m2.ax.transAxes)

        self.m3.contourf(self.lon,self.lat,self.conc3,contour_levels,latlon=True,cmap=cmo.cm.matter,extend='max')
        self.m3.ax.text(0.9,0.91,'%i days'%(time3),transform=self.m3.ax.transAxes,ha='center',va='center',fontsize=8)
        self.m3.ax.text(0.03,0.85,'(c)',transform=self.m3.ax.transAxes)

        self.m4.contourf(self.lon,self.lat,self.conc4,contour_levels,latlon=True,cmap=cmo.cm.matter,extend='max')
        self.m4.ax.text(0.9,0.91,'%i days'%(time4),transform=self.m4.ax.transAxes,ha='center',va='center',fontsize=8)
        self.m4.ax.text(0.03,0.85,'(d)',transform=self.m4.ax.transAxes)


        if hasattr(self,'subAdjust_left'):
            rect = (0.13,0.16,1.,1.)
            plt.tight_layout(rect=rect)
            plt.subplots_adjust(top=self.subAdjust_top,bottom=self.subAdjust_bottom,left=self.subAdjust_left,right=self.subAdjust_right,wspace=self.subAdjust_wspace,hspace=self.subAdjust_hspace)

        if hasattr(self,'figname'):
            plt.savefig('/media/danilo/Danilo/mestrado/github/artigoTG/figures/%s.png'%(self.figname),dpi=600)
        else:
            plt.show()


# ----- load file from Oceano Computer ----- #
DATA_DIR = '/media/danilo/Danilo/mestrado/artigo_data/simulacoes/sims_dispersao/'

# ----
def fig2():
    expI = Experiment(DATA_DIR+'expI.cdf',figsize=(8.4,10.))
    expI.figname = 'Fig2'
    expI.plot_Figure_WindDriven_Circulation()

def fig3():
    expII = Experiment(DATA_DIR+'expII.cdf',figsize=(8.4,10.))
    expII.figname = 'Fig3'
    expII.plot_Figure_WindDriven_Circulation()

def fig4():
    expIII = Experiment(DATA_DIR+'expIII.cdf',figsize=(17.4,12))
    # expIII.figname = 'Fig4'
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
    expIII.plot_Figure4(229,218,290,293)

def fig5():
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
    # expIII.plot_Figure5(229,218,290,293)

    # if you want to plot 4 subplots for each tidal condition and phase, uncomment the 7 lines below
    expIII.figSize = (8.4/2.54,10./2.54)
    expIII.figname = 'Fig5'
    expIII.define_adjustSubplot_Parameters(top=0.977,bottom=0.172,left=0.072,right=0.973,hspace=0.072,wspace=0.082)
    expIII.plot_Figure5_2ndVersion(229,218,290,293)

def fig6():
    instantes = [3, 10, 21, 60]
    expIV = Experiment(DATA_DIR+"expIV.cdf",figsize=(17.4,10.))
    expIV.resolution = 'f'
    expIV.importVariables_basic()
    expIV.figname = 'Fig6'
    expIV.subAdjust_top    = 0.995
    expIV.subAdjust_bottom = 0.160
    expIV.subAdjust_left   = 0.125
    expIV.subAdjust_right  = 0.875
    expIV.subAdjust_hspace = 0.000
    expIV.subAdjust_wspace = 0.026

    expIV.plot_Concentration(instantes[0],instantes[1],instantes[2],instantes[3])

def fig7():
    instantes = [3, 10, 27, 60]
    expV = Experiment(DATA_DIR+"expV.cdf",figsize=(17.4,10.))
    expV.resolution = 'f'
    expV.importVariables_basic()
    expV.figname = 'Fig7'
    expV.subAdjust_top    = 0.995
    expV.subAdjust_bottom = 0.160
    expV.subAdjust_left   = 0.125
    expV.subAdjust_right  = 0.875
    expV.subAdjust_hspace = 0.000
    expV.subAdjust_wspace = 0.026

    expV.plot_Concentration(instantes[0],instantes[1],instantes[2],instantes[3])

def fig8():
    instantes = [7, 12, 22, 30]
    expVI = Experiment(DATA_DIR+"expVI.cdf",figsize=(17.4,10.))
    expVI.resolution = 'f'
    expVI.importVariables_basic()
    expVI.figname = 'Fig8'
    expVI.subAdjust_top    = 0.995
    expVI.subAdjust_bottom = 0.160
    expVI.subAdjust_left   = 0.125
    expVI.subAdjust_right  = 0.875
    expVI.subAdjust_hspace = 0.000
    expVI.subAdjust_wspace = 0.026

    expVI.plot_Concentration(instantes[0],instantes[1],instantes[2],instantes[3])
