import matplotlib.pyplot as plt
import numpy as np
import cmocean as cmo

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano


class Animation:

    def __init__(self,var='elev',sigma=0):
        plt.ion()

        self.varname = var

        if len(self.ncin[var].shape) == 4:
            # is a 3D var (temp,salt,others) with sigma level in axis 1
            self.var = self.ncin[var][self.timeStart.item():self.timeEnd.item(),sigma,:,:]
        elif len(self.ncin[var].shape) == 3:
            # is a 2D var with time in axis 0
            self.var = self.ncin[var][self.timeStart.item():self.timeEnd.item(),:,:]

    def select_colormap(self):
        if self.varname == 'elev':
            cmap = 'RdBu_r'
        elif self.varname == 'temp':
            cmap = cmo.cm.thermal
        elif self.varname == 'salt':
            cmap = cmo.cm.haline

        return cmap

    def velocidade(self,index_file,intervalo=.2,sigma=0):

        # checar se existe u e v no objeto
        if ~hasattr(self,'s'):
            self.u = self.ncin['u'][self.timeStart.item():self.timeEnd.item(),sigma,:,:]
            self.v = self.ncin['v'][self.timeStart.item():self.timeEnd.item(),sigma,:,:]
            self.s = np.sqrt(self.u**2 + self.v**2)

        # normalizando as componentes
        un = self.u / self.s
        vn = self.v / self.s

        # preparar as variaveis para plot
        xplot = oceano.formatGrid_plot(self.lon,index_file)
        yplot = oceano.formatGrid_plot(self.lat,index_file)
        contour_levels = np.arange(np.nanmin(self.s),np.nanmax(self.s),0.2)

        # animacao
        fig,ax = plt.subplots()
        for t in np.arange(self.timeStart.item(),self.timeEnd.item(),1):
            ax.clear()
            m = oceano.make_map(ax)
            plt.suptitle('Velocity at: '+ str(self.ncin.time[t].values))

            uplot = oceano.formatGrid_plot(un[t,:,:],index_file)
            vplot = oceano.formatGrid_plot(vn[t,:,:],index_file)

            m.contourf(self.lon,self.lat,self.s[t,:,:],contour_levels,latlon=True,cmap=cmo.cm.speed)
            m.quiver(xplot,yplot,uplot,vplot,scale=80,width=0.001,headwidth=4,
                            headlength=4,alpha=0.6,latlon=True,pivot='middle')

            plt.pause(intervalo)

    def field(self,title='',**kwargs):

        # in case user don't send a cmap, we use the variable name
        # to select a proper colormap
        if 'cmap' not in kwargs.keys():
            kwargs['cmap'] = self.select_colormap()

        # we create contour_levels outside the for loop, so we can
        # use the extreme values in the entire time domain
        # dif = np.abs(np.nanmin(self.var.values)-np.nanmax(self.var.values))
        # contour_levels = np.arange(np.nanmin(self.var.values),np.nanmax(self.var.values),round(dif/10.))

        if (self.var.ndim > 2):
            # if self.var is a 3D array, plot an animation
            plt.ion()

            fig,ax = plt.subplots()
            for t in np.arange(self.timeStart.item(),self.timeEnd.item(),1):
                ax.clear()

                data = self.var[t]
                # m = oceano.make_map(ax)
                self.Mapa(ax)

                self.mapa.contourf(self.lon,self.lat,data,**kwargs)
                plt.title(title + '\n' +str(self.ncin.time[t].values))
                plt.pause(0.5)
        else:
            # if self.var is a 2D array, plot a static map
            plt.ion()

            fig,self.ax = plt.subplots()
            self.Mapa(self.ax)
            self.csf = self.mapa.contourf(self.lon,self.lat,self.var,**kwargs)
            plt.colorbar(self.csf)
            plt.title(title,fontsize=12)


    def cross_section(self,**kwargs):
        pass


class Visualization:

    def __init__(self,var='elev',sigma=0):
        self.varname = var

        if len(self.ncin[var].shape) == 4:
            # is a 3D var (temp,salt,others) with sigma level in axis 1
            self.var = self.ncin[var][self.timeStart.item():self.timeEnd.item(),sigma,self.ilat,self.ilon]
        elif len(self.ncin[var].shape) == 3:
            # is a 2D var with time in axis 0
            self.var = self.ncin[var][self.timeStart.item():self.timeEnd.item(),self.ilat,self.ilon]

        self.time = self.ncin.time[self.timeStart.item():self.timeEnd.item()]

    def graph(self,title='',**kwargs):
        fig,ax = plt.subplots()
        ax.plot(self.time,self.var,'k')
        ax.set_ylabel(self.var.attrs['long_name'] + " [" +self.var.attrs['units'] + "]",fontsize=8)
        ax.set_xlabel(self.time.attrs['long_name'],fontsize=8)
        ax.set_title(title,fontsize=12)
        ax.margins(0)
