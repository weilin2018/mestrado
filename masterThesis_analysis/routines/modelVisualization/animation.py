import matplotlib.pyplot as plt
import numpy as np
import cmocean as cmo

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano


class Animation:

    def __init__(self,var='elev',sigma=0):
        self.varname = var

        if len(self.ncin[var].shape) == 4:
            # is a 3D var (temp,salt,others) with sigma level in axis 1
            self.var = self.ncin[var][self.timeStart.item():self.timeEnd.item(),sigma,:,:].values
        elif len(self.ncin[var].shape) == 3:
            # is a 2D var with time in axis 0
            self.var = self.ncin[var][self.timeStart.item():self.timeEnd.item(),:,:].values

    def select_colormap(self):
        if self.varname == 'elev':
            cmap = 'RdBu_r'
        elif self.varname == 'temp':
            cmap = cmo.cm.thermal
        elif self.varname == 'salt':
            cmap = cmo.cm.haline

        return cmap

    def field(self,**kwargs):
        plt.ion()

        fig,ax = plt.subplots()
        for t in np.arange(self.timeStart.item(),self.timeEnd.item(),1):
            ax.clear()

            data = self.var[t]
            # m = oceano.make_map(ax)
            self.Mapa(ax)

            self.mapa.contourf(self.lon,self.lat,data,**kwargs)
            plt.title(str(self.ncin.time[t].values))
            plt.pause(0.5)

    def cross_section(self,**kwargs):
        pass


class Visualization:

    def __init__(self,ilon,ilat,var='elev',sigma=0):
        self.varname = var

        if len(self.ncin[var].shape) == 4:
            # is a 3D var (temp,salt,others) with sigma level in axis 1
            self.var = self.ncin[var][self.timeStart.item():self.timeEnd.item(),sigma,ilat,ilon].values
        elif len(self.ncin[var].shape) == 3:
            # is a 2D var with time in axis 0
            self.var = self.ncin[var][self.timeStart.item():self.timeEnd.item(),ilat,ilon].values
